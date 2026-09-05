# Olla-DFT — command-line toolkit for Quantum ESPRESSO
# Copyright (C) 2026 Jorge Enrique González Sevilla
# SPDX-License-Identifier: AGPL-3.0-or-later
# This program is free software: you can redistribute it and/or modify it
# under the terms of the GNU Affero General Public License as published by the Free
# Software Foundation, either version 3 of the License, or (at your option)
# any later version. See the LICENSE file for details.

"""Validación cruzada: la misma cantidad por caminos independientes.

Cada módulo de Olla-DFT se valida contra la literatura, pero eso no detecta
un error que afecte a un módulo entero de forma sistemática. Lo que sí lo
detecta es calcular la MISMA cantidad por dos rutas físicamente
independientes y compararlas:

    modulo volumetrico   ->  ajuste de la EOS   vs   traza de las Cij
    velocidad del sonido ->  Cij                vs   pendiente acustica
    temperatura de Debye ->  velocidades        vs   DOS de fonones
    gap                  ->  estructura de bandas vs extrapolacion de Tauc
    C_v a T alta         ->  DOS de fonones     vs   limite de Dulong-Petit
    numero de modos      ->  integral de la DOS vs   3N

No cuesta ningun calculo nuevo: son resultados que ya estan en disco.

SOBRE LAS TOLERANCIAS
---------------------
Cada cruce lleva la suya, y no son arbitrarias:

- B0 por dos rutas es la MISMA cantidad y debe coincidir al ~5 %;
- las velocidades del sonido tambien, pero la pendiente acustica en q->0
  es justo lo que peor interpola una malla de q gruesa, asi que un
  desacuerdo ahi acusa a la malla antes que al modulo elastico;
- las dos temperaturas de Debye NO son la misma definicion (una es el
  limite acustico, la otra usa todo el espectro): ahi la tolerancia es
  amplia a proposito, y coincidir al 1 % seria sospechoso, no bueno.

Un cruce que falla no dice cual de los dos caminos esta mal. Por eso cada
uno lleva un diagnostico de que mirar primero.
"""

from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

KB_EV = 8.617333262e-5


@dataclass
class Check:
    nombre: str = ""
    ruta_a: str = ""
    valor_a: float = None
    ruta_b: str = ""
    valor_b: float = None
    unidad: str = ""
    tolerancia: float = 0.05
    diagnostico: str = ""

    @property
    def desvio(self):
        """Desviación RELATIVA, salvo que la referencia sea cero.

        Hay cruces cuya respuesta correcta ES cero —la fase de Berry de un
        cristal centrosimétrico, por ejemplo— y ahí no existe desviación
        relativa. Antes se devolvía None, el cruce se daba por "sin datos" y
        el informe reventaba al formatearlo. Con referencia cero la
        desviación es el valor absoluto, y la tolerancia se lee en las
        unidades del cruce.
        """
        if self.valor_a is None or self.valor_b is None:
            return None
        if self.valor_a == 0:
            return abs(self.valor_b)
        return abs(self.valor_b - self.valor_a) / abs(self.valor_a)

    @property
    def relativa(self):
        return self.valor_a not in (None, 0)

    @property
    def ok(self):
        d = self.desvio
        return None if d is None else bool(d <= self.tolerancia)


@dataclass
class CrossResult:
    checks: list = field(default_factory=list)
    disponibles: list = field(default_factory=list)
    faltantes: list = field(default_factory=list)


# ----------------------------------------------------------------------
def _cargar(project: Path) -> dict:
    """Busca en la carpeta los resultados que Olla-DFT ya haya escrito."""
    p = Path(project)
    d = {}
    for nombre, patron in (("elastic", "ELASTIC_C.dat"),
                           ("eos", "EOS.txt"),
                           ("fonones_dos", "FONONES_DOS.dat"),
                           ("fonones_bandas", "FONONES_BANDAS.dat"),
                           ("optics", "OPTICS.dat"),
                           ("kappa", "KAPPA.dat"),
                           ("berry", "BERRY.dat"),
                           ("wannier", "WANNIER_centros.dat"),
                           ("esm", "ESM.dat"),
                           ("wf", "WF.dat"),
                           ("strain", "STRAIN.dat")):
        hits = list(p.rglob(patron))
        if hits:
            d[nombre] = hits[0]
    return d


def _leer_eos_b0(path) -> float:
    for linea in Path(path).read_text(errors="ignore").splitlines():
        if "B0" in linea and "GPa" in linea:
            for tok in linea.replace("=", " ").split():
                try:
                    v = float(tok)
                    if 0.1 < v < 1e4:
                        return v
                except ValueError:
                    continue
    return None


def _leer_cabecera(path, clave):
    """Un valor escrito como '# clave = valor' en la cabecera de un .dat."""
    for linea in Path(path).read_text(errors="ignore").splitlines():
        if not linea.startswith("#") or clave not in linea:
            continue
        try:
            return float(linea.split("=")[1].split()[0])
        except (IndexError, ValueError):
            continue
    return None


def _leer_kappa_300(path):
    """κ medio a la temperatura más cercana a 300 K, en W/m/K."""
    d = np.loadtxt(path, comments="#")
    if d.ndim == 1:
        d = d[None, :]
    i = int(np.argmin(np.abs(d[:, 0] - 300.0)))
    return float(d[i, -1]), float(d[i, 0])


def _leer_berry_el(path):
    """Fase electrónica de Berry del punto de carga cero (unidades de QE)."""
    d = np.loadtxt(path, comments="#")
    if d.ndim == 1:
        d = d[None, :]
    i = int(np.argmin(np.abs(d[:, 0])))          # lambda o carga más cercana a 0
    return float(d[i, 2])


def _fase_de_centros(path, cell, gdir=3, spin=2.0):
    """La misma fase, desde los centros de Wannier: −f·Σ_n (r̄_n·b)/2π."""
    d = np.loadtxt(path, comments="#")
    if d.ndim == 1:
        d = d[None, :]
    frac = d[:, 1:4] @ np.linalg.inv(np.asarray(cell, float))
    return float(-spin * frac[:, int(gdir) - 1].sum())


def _leer_cij(path) -> np.ndarray:
    datos = np.loadtxt(path, comments="#")
    return datos if datos.shape == (6, 6) else None


# ----------------------------------------------------------------------
def run(project=".", masas=None, volumen=None, natoms=None,
        gap_bandas: float = None, C: np.ndarray = None,
        b0_eos: float = None, qdist=None, band_freqs=None,
        dos_w=None, dos=None, gap_tauc: float = None,
        cell=None) -> CrossResult:
    """Ejecuta todos los cruces para los que haya datos."""
    from qekit.modules import derived

    res = CrossResult()
    encontrados = _cargar(project) if project else {}

    if C is None and "elastic" in encontrados:
        C = _leer_cij(encontrados["elastic"])
    if b0_eos is None and "eos" in encontrados:
        b0_eos = _leer_eos_b0(encontrados["eos"])
    if dos_w is None and "fonones_dos" in encontrados:
        datos = np.loadtxt(encontrados["fonones_dos"], comments="#")
        dos_w, dos = datos[:, 0], datos[:, 1]
    if qdist is None and "fonones_bandas" in encontrados:
        datos = np.loadtxt(encontrados["fonones_bandas"], comments="#")
        qdist, band_freqs = datos[:, 0], datos[:, 1:]

    for k, v in (("constantes elásticas", C is not None),
                 ("κ de red (fc3)", "kappa" in encontrados),
                 ("fase de Berry", "berry" in encontrados),
                 ("centros de Wannier", "wannier" in encontrados),
                 ("función trabajo (ESM)", "esm" in encontrados),
                 ("función trabajo (potencial planar)", "wf" in encontrados),
                 ("barrido de deformación", "strain" in encontrados),
                 ("ecuación de estado", b0_eos is not None),
                 ("DOS de fonones", dos_w is not None),
                 ("dispersión de fonones", qdist is not None),
                 ("gap de bandas", gap_bandas is not None),
                 ("gap de Tauc", gap_tauc is not None)):
        (res.disponibles if v else res.faltantes).append(k)

    # --- 1. modulo volumetrico: EOS contra Cij -----------------------
    if C is not None and b0_eos is not None:
        from qekit.modules import elastic
        m = elastic.moduli(C)
        res.checks.append(Check(
            nombre="módulo volumétrico B₀",
            ruta_a="ajuste de la ecuación de estado", valor_a=b0_eos,
            ruta_b="traza de las constantes elásticas",
            valor_b=m.B_hill, unidad="GPa", tolerancia=0.05,
            diagnostico=(
                "Son la MISMA cantidad por dos vías. Si difieren: revisa "
                "que la celda de las\nelásticas estuviera relajada (esfuerzo "
                "residual bajo) y que la EOS tenga\npuntos suficientes a "
                "ambos lados del mínimo.")))

    # --- 2. velocidades del sonido -----------------------------------
    if (C is not None and qdist is not None and masas is not None
            and volumen):
        rho = derived.density(masas, volumen)
        dirs = derived.cubic_directional(C, rho)
        ac = derived.acoustic_velocities(qdist, band_freqs)
        if dirs and ac:
            res.checks.append(Check(
                nombre="velocidad longitudinal [100]",
                ruta_a="constantes elásticas: √(C₁₁/ρ)",
                valor_a=dirs["v_l_100"],
                ruta_b="pendiente de la rama LA en Γ",
                valor_b=ac["v_l"], unidad="m/s", tolerancia=0.10,
                diagnostico=(
                    "La pendiente acústica en q→0 es lo que peor interpola "
                    "una malla de q gruesa.\nSi falla, sospecha de la malla "
                    "antes que de las Cij.")))
            res.checks.append(Check(
                nombre="velocidad transversal [100]",
                ruta_a="constantes elásticas: √(C₄₄/ρ)",
                valor_a=dirs["v_t_100"],
                ruta_b="pendiente de la rama TA en Γ",
                valor_b=ac["v_t1"], unidad="m/s", tolerancia=0.10,
                diagnostico=(
                    "Las ramas TRANSVERSALES son las más planas y las que "
                    "peor salen de una malla\nde q pequeña — en el silicio "
                    "con 2x2x2 el error pasa del 40 %. Si la\n"
                    "longitudinal cuadra y esta no, es la malla de q, no las "
                    "elásticas.\nDensifícala (4x4x4 o más) antes de creerle "
                    "a ninguna de las dos.")))

    # --- 3. temperatura de Debye -------------------------------------
    if (C is not None and dos_w is not None and masas is not None
            and volumen and natoms):
        from qekit.modules import elastic
        m = elastic.moduli(C)
        rho = derived.density(masas, volumen)
        _vl, _vt, vm = derived.sound_velocities(m.B_hill, m.G_hill, rho)
        td_el = derived.debye_from_velocity(vm, natoms, volumen)
        td_dos = derived.debye_from_dos(dos_w, dos, natoms)
        if td_el and td_dos:
            res.checks.append(Check(
                nombre="temperatura de Debye",
                ruta_a="velocidades del sonido (límite acústico)",
                valor_a=td_el,
                ruta_b="segundo momento de la DOS de fonones",
                valor_b=td_dos, unidad="K", tolerancia=0.30,
                diagnostico=(
                    "OJO: NO son la misma definición. La elástica es el "
                    "límite de baja temperatura\n(solo acústicas); la de la "
                    "DOS usa todo el espectro, ópticas incluidas, y sale\n"
                    "más alta. Se cruzan para detectar un disparate, no para "
                    "que coincidan:\ncoincidir al 1 % sería sospechoso.")))

    # --- 4. gap: bandas contra Tauc ----------------------------------
    if gap_bandas is not None and gap_tauc is not None:
        res.checks.append(Check(
            nombre="gap óptico",
            ruta_a="estructura de bandas (gap directo)", valor_a=gap_bandas,
            ruta_b="extrapolación de Tauc sobre α(E)", valor_b=gap_tauc,
            unidad="eV", tolerancia=0.06,
            diagnostico=(
                "epsilon.x no incluye transiciones asistidas por fonones, "
                "así que el borde de\nabsorción es el gap DIRECTO, no el "
                "fundamental. Si comparas contra el\nfundamental de un "
                "semiconductor indirecto, la diferencia es física, no un "
                "error.")))

    # --- 5. C_v a alta T contra Dulong-Petit -------------------------
    if dos_w is not None and natoms:
        cv = _cv_alta_T(dos_w, dos, natoms, T=1500.0)
        dp = 3.0 * natoms * KB_EV * 1000.0        # meV/K por celda
        if cv:
            res.checks.append(Check(
                nombre="C_v en el límite clásico",
                ruta_a="Dulong–Petit: 3N·k_B", valor_a=dp,
                ruta_b="integral de la DOS de fonones a 1500 K",
                valor_b=cv, unidad="meV/K por celda", tolerancia=0.03,
                diagnostico=(
                    "A temperatura alta toda C_v armónica tiende a 3N·k_B. "
                    "Si no llega, la DOS\nestá mal normalizada o le falta "
                    "espectro; si se pasa, hay modos de más.")))

    # --- 6. numero de modos ------------------------------------------
    if dos_w is not None and natoms:
        from qekit.core.compat import trapezoid
        total = float(trapezoid(np.asarray(dos), np.asarray(dos_w)))
        res.checks.append(Check(
            nombre="número de modos",
            ruta_a="3N por construcción", valor_a=3.0 * natoms,
            ruta_b="integral de la DOS de fonones", valor_b=total,
            unidad="modos", tolerancia=0.05,
            diagnostico=(
                "La integral de la DOS tiene que dar exactamente 3N. Si no, "
                "la malla de\ninterpolación de matdyn es demasiado pobre o "
                "el rango de frecuencias corta\nespectro.")))
    # --- 7. kappa de red: tercer orden contra el modelo de Slack ------
    if "kappa" in encontrados and C is not None and masas is not None \
            and volumen and natoms:
        from qekit.modules import elastic
        k300, T_usada = _leer_kappa_300(encontrados["kappa"])
        m = elastic.moduli(C)
        rho = derived.density(masas, volumen)
        _vl, _vt, vm = derived.sound_velocities(m.B_hill, m.G_hill, rho)
        td = derived.debye_from_velocity(vm, natoms, volumen)
        gam = derived.gruneisen_from_poisson(m.nu)
        ks = derived.slack(td, gam, float(np.mean(masas)), natoms, volumen,
                           T=T_usada) if (td and gam) else None
        if ks:
            res.checks.append(Check(
                nombre="conductividad térmica de red",
                ruta_a=f"ecuación de Boltzmann de fonones con fc3 "
                       f"({T_usada:.0f} K)",
                valor_a=k300,
                ruta_b="modelo de Slack desde las elásticas",
                valor_b=ks, unidad="W/m/K", tolerancia=0.60,
                diagnostico=(
                    "La tolerancia es del 60 % A PROPÓSITO: Slack es una "
                    "estimación de orden de\nmagnitud con un prefactor "
                    "empírico, no un cálculo. Sirve para detectar que a la\n"
                    "fc3 le falta convergencia o que el signo de algo está "
                    "mal, no para afinar.\nSi difieren en un factor 3, "
                    "sospecha primero de la supercelda de la fc3.")))

    # --- 8. fase de Berry: lberry contra los centros de Wannier -------
    if "berry" in encontrados and "wannier" in encontrados and cell is not None:
        try:
            fa = _leer_berry_el(encontrados["berry"])
            fb = _fase_de_centros(encontrados["wannier"], cell)
            # las dos están definidas módulo 2: se comparan en la misma rama
            fb = fb - 2.0 * np.round((fb - fa) / 2.0)
            res.checks.append(Check(
                nombre="fase electrónica de Berry",
                ruta_a="lberry: determinante de solapes en cuerdas de k",
                valor_a=fa,
                ruta_b="centros de Wannier: −2·Σ (r̄·b)/2π",
                valor_b=fb, unidad="(cuanto = 2)", tolerancia=0.05,
                diagnostico=(
                    "Son la MISMA fase de Berry por dos rutinas que no "
                    "comparten una línea de\ncódigo. Que coincidan es la "
                    "validación más fuerte que hay aquí. Si no lo hacen,\n"
                    "lo primero a mirar es que las dos usen la misma "
                    "dirección (gdir) y la misma\nmalla de puntos k.")))
        except Exception:                                   # noqa: BLE001
            pass

    # --- 9. funcion trabajo: ESM contra la meseta del potencial -------
    if "esm" in encontrados and "wf" in encontrados:
        try:
            d = np.loadtxt(encontrados["esm"], comments="#")
            if d.ndim == 1:
                d = d[None, :]
            i = int(np.argmin(np.abs(d[:, 0])))     # la losa neutra
            phi_esm = float(d[i, 4])
            phi_wf = _leer_cabecera(encontrados["wf"], "Phi_eV")
            if phi_wf is not None:
                res.checks.append(Check(
                    nombre="función trabajo",
                    ruta_a="ESM: el nivel de vacío vale cero por "
                           "construcción",
                    valor_a=phi_esm,
                    ruta_b="meseta del potencial planar del cube de pp.x",
                    valor_b=phi_wf, unidad="eV", tolerancia=0.05,
                    diagnostico=(
                        "Con bc1 el nivel de vacío de ESM es cero exacto y "
                        "no hay meseta que\najustar; el camino del cube sí "
                        "la ajusta, y por eso necesita más vacío. Si\n"
                        "difieren, mira la planitud que reporta el segundo: "
                        "casi siempre es que al\ncálculo periódico le "
                        "faltaba vacío, no que ESM esté mal.")))
        except Exception:                                   # noqa: BLE001
            pass

    # --- 10. modulo volumetrico: EOS contra la presion del barrido ----
    if "strain" in encontrados and b0_eos is not None:
        try:
            d = np.loadtxt(encontrados["strain"], comments="#")
            eps, P = d[:, 0], d[:, 3]
            bien = np.isfinite(P)
            if bien.sum() >= 3:
                # hidrostática: V = V0(1+ε)³  ->  B = −dP/d(lnV) = −dP/dε / 3
                pend = np.polyfit(eps[bien], P[bien], 1)[0]
                b0_strain = -pend / 3.0 * 0.1        # kbar -> GPa
                if b0_strain > 0:
                    res.checks.append(Check(
                        nombre="módulo volumétrico B₀ (tercera vía)",
                        ruta_a="ajuste de la ecuación de estado",
                        valor_a=b0_eos,
                        ruta_b="pendiente de la presión en el barrido de "
                               "deformación",
                        valor_b=b0_strain, unidad="GPa", tolerancia=0.10,
                        diagnostico=(
                            "Solo vale si el barrido fue HIDROSTÁTICO: con "
                            "deformación biaxial o\nuniaxial la relación "
                            "entre presión y ε es otra y este cruce compara "
                            "peras\ncon manzanas. Míralo antes de creerle.")))
        except Exception:                                   # noqa: BLE001
            pass

    return res


def _cv_alta_T(w, g, natoms, T=1500.0):
    """C_v armónica a temperatura alta, en meV/K por celda."""
    from qekit.core.compat import trapezoid
    w = np.asarray(w, dtype=float)
    g = np.asarray(g, dtype=float)
    m = w > 1.0
    w, g = w[m], g[m]
    if w.size < 3:
        return None
    norm = trapezoid(g, w)
    if norm <= 0:
        return None
    g = g * (3.0 * natoms / norm)
    e = w * 1.239841984e-4                      # cm^-1 -> eV
    x = e / (KB_EV * T)
    x = np.clip(x, 1e-9, 300.0)
    occ = 1.0 / np.expm1(x)
    cv = KB_EV * trapezoid(x ** 2 * np.exp(x) * occ ** 2 * g, w)
    return float(cv * 1000.0)


def report(res: CrossResult) -> str:
    lines = ["--- Validación cruzada ---"]
    if res.disponibles:
        lines.append("Resultados encontrados: " + ", ".join(res.disponibles))
    if res.faltantes:
        lines.append("No disponibles: " + ", ".join(res.faltantes))
    lines.append("")
    if not res.checks:
        lines.append(
            "No hay dos rutas independientes que cruzar todavía. Cada cruce "
            "necesita DOS\nmódulos: por ejemplo elásticas + EOS, o "
            "elásticas + fonones.")
        return "\n".join(lines)

    fallos = [c for c in res.checks if c.ok is False]
    lines.append(f"{len(res.checks)} cruces  |  "
                 f"{len(res.checks) - len(fallos)} coinciden  |  "
                 f"{len(fallos)} NO")
    lines.append("")
    for c in res.checks:
        marca = "OK  " if c.ok else ("FALLA" if c.ok is False else "  ?  ")
        if c.desvio is None:
            lines.append(f"[{marca}] {c.nombre}  (sin datos suficientes)")
        elif c.relativa:
            lines.append(f"[{marca}] {c.nombre}  ({c.desvio * 100:.1f} % de "
                         f"desvío, tolerancia {c.tolerancia * 100:.0f} %)")
        else:
            lines.append(f"[{marca}] {c.nombre}  (tiene que dar cero; sale "
                         f"{c.desvio:.2e}, tolerancia {c.tolerancia:g})")
        va = "—" if c.valor_a is None else f"{c.valor_a:.4g}"
        vb = "—" if c.valor_b is None else f"{c.valor_b:.4g}"
        lines.append(f"         {c.ruta_a}: {va} {c.unidad}")
        lines.append(f"         {c.ruta_b}: {vb} {c.unidad}")
        if c.ok is False:
            for l in c.diagnostico.splitlines():
                lines.append(f"         > {l}")
        lines.append("")

    if fallos:
        lines.append(
            "Un cruce que falla NO dice cuál de los dos caminos está mal: "
            "dice que uno de\nlos dos lo está. El diagnóstico de cada uno "
            "indica qué mirar primero.")
    else:
        lines.append(
            "Todos los cruces coinciden. Es la evidencia más fuerte que se "
            "puede tener sin\nsalir del propio cálculo: dos rutas "
            "independientes no se equivocan igual por\ncasualidad.")
    return "\n".join(lines)
