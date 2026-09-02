# Olla-DFT — command-line toolkit for Quantum ESPRESSO
# Copyright (C) 2026 Jorge Enrique González Sevilla
# SPDX-License-Identifier: GPL-3.0-or-later
# This program is free software: you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by the Free
# Software Foundation, either version 3 of the License, or (at your option)
# any later version. See the LICENSE file for details.

"""Energía de formación de defectos cargados.

La fórmula tiene cinco términos y cada uno esconde una trampa distinta:

    E_f[D^q](ε_F, μ) = E[D^q] − E[perfecto] − Σ n_i μ_i + q(ε_VBM + ε_F) + E_corr

  E[D^q] − E[perfecto]   las dos superceldas TIENEN que compartir celda,
                         cutoff, malla k y pseudos, o la resta no significa nada.
  Σ n_i μ_i              el potencial químico no es una constante del material:
                         depende de las condiciones de síntesis y acota E_f
                         entre dos límites en vez de fijarla.
  q(ε_VBM + ε_F)         el electrón que entra o sale viene de un reservorio
                         cuya energía es el nivel de Fermi, medido desde el VBM.
  E_corr                 una celda cargada no es neutra: QE mete un fondo
                         uniforme para poder resolver Poisson, y la carga
                         interacciona con sus propias imágenes. Sin corregirlo,
                         E_f puede errar por más de 1 eV en una supercelda
                         normal, y el error CRECE con q².

La constante de Madelung de la corrección se calcula aquí por suma de Ewald
sobre la celda real, no se toma de una tabla: así vale para cualquier forma
de supercelda y no hay que confiar en qué convenio usaba la tabla.
"""

from dataclasses import dataclass, field
from math import erfc, pi, sqrt
from pathlib import Path

import numpy as np

from qekit.core import provenance, qeout
from qekit.core.errors import ErrorDeUso, FaltanDatos
from qekit.core import style as qstyle
from qekit.modules import builder, sweep

# e²/(4πε₀) en eV·Å
KE = 14.399645

# Factor de forma de Lany-Zunger: recoge el término de orden superior del
# apantallamiento, que el monopolo de Makov-Payne no lleva. Varía poco con la
# forma de la celda (−0.369 cúbica simple, −0.343 fcc, −0.342 bcc), así que
# se usa un valor único y se dice que lo es.
C_SHAPE = -0.35

ESQUEMAS = ("ninguna", "makov-payne", "lany-zunger")


def etiqueta_q(q: int) -> str:
    """'−1', '0', '+1'.  Ni '+0' ni '-0', que se leen como erratas."""
    q = int(q)
    return "0" if q == 0 else f"{q:+d}"


# ----------------------------------------------------------------------
# Constante de Madelung por suma de Ewald
# ----------------------------------------------------------------------
def madelung_xi(cell, tol: float = 1e-10) -> float:
    """Energía electrostática de una carga unidad en un fondo neutralizante.

    Devuelve ξ (negativo) en Å⁻¹, tal que E = q²·ξ/2 en unidades e²/Å.
    Verificado contra el valor de libro de la red cúbica simple
    (α = 2.8372974) con siete cifras, y comprobado que no depende de la
    escala de la celda.
    """
    cell = np.asarray(cell, dtype=float)
    V = abs(np.linalg.det(cell))
    if V <= 0:
        raise ErrorDeUso("la celda tiene volumen cero; no puedo calcular la "
                         "corrección de imagen.")
    eta = sqrt(pi) / V ** (1.0 / 3.0)
    recip = 2 * pi * np.linalg.inv(cell).T

    rc = 1.0
    while erfc(eta * rc) / rc > tol and rc < 500:
        rc += 0.5
    gc = 2 * eta * sqrt(-np.log(tol))

    def _profundidad(vectores, corte):
        n = []
        for i in range(3):
            alt = abs(np.linalg.det(vectores)) / np.linalg.norm(
                np.cross(vectores[(i + 1) % 3], vectores[(i + 2) % 3]))
            n.append(int(np.ceil(corte / alt)) + 1)
        return n

    def _malla(n):
        ejes = [np.arange(-k, k + 1) for k in n]
        I, J, K = np.meshgrid(*ejes, indexing="ij")
        idx = np.stack([I.ravel(), J.ravel(), K.ravel()], axis=1)
        return idx[np.any(idx != 0, axis=1)]

    R = _malla(_profundidad(cell, rc)) @ cell
    r = np.linalg.norm(R, axis=1)
    r = r[r <= rc]
    s_real = float(np.sum(np.array([erfc(eta * x) for x in r]) / r))

    G = _malla(_profundidad(recip, gc)) @ recip
    g2 = np.sum(G * G, axis=1)
    g2 = g2[g2 <= gc ** 2]
    s_rec = float((4 * pi / V) * np.sum(np.exp(-g2 / (4 * eta ** 2)) / g2))

    return s_real + s_rec - 2 * eta / sqrt(pi) - pi / (eta ** 2 * V)


def constante_madelung(cell) -> float:
    """α_M de la celda, con la longitud de referencia L = V^(1/3)."""
    cell = np.asarray(cell, dtype=float)
    L = abs(np.linalg.det(cell)) ** (1.0 / 3.0)
    return float(-madelung_xi(cell) * L)


def correccion_imagen(q: int, cell, epsilon: float = None,
                      esquema: str = "lany-zunger") -> dict:
    """Corrección de tamaño finito para una celda con carga q.

    Devuelve un diccionario con el término de Makov-Payne, el de
    Lany-Zunger y el que corresponde al esquema pedido. La corrección es
    POSITIVA: la celda cargada sale artificialmente estabilizada por su
    propio fondo y hay que devolverle esa energía.
    """
    if esquema not in ESQUEMAS:
        raise ErrorDeUso(
            f"esquema de corrección desconocido '{esquema}'. "
            f"Opciones: {', '.join(ESQUEMAS)}.")
    if q == 0 or esquema == "ninguna":
        return {"E_mp": 0.0, "E_lz": 0.0, "E_corr": 0.0, "alpha": None,
                "esquema": "ninguna" if q == 0 else esquema}
    if not epsilon or epsilon <= 0:
        raise ErrorDeUso(
            "para corregir una celda cargada hace falta la constante "
            "dieléctrica del material (--epsilon). Sin apantallar, la "
            "corrección sale ε veces demasiado grande: en silicio (ε≈11.7) "
            "eso es un factor 12. Si de verdad quieres verla sin corregir, "
            "usa --correction ninguna y el reporte lo dirá.")
    cell = np.asarray(cell, dtype=float)
    L = abs(np.linalg.det(cell)) ** (1.0 / 3.0)
    alpha = constante_madelung(cell)
    e_mp = KE * q ** 2 * alpha / (2.0 * epsilon * L)
    e_lz = e_mp * (1.0 + C_SHAPE * (1.0 - 1.0 / epsilon))
    return {"E_mp": float(e_mp), "E_lz": float(e_lz),
            "E_corr": float(e_mp if esquema == "makov-payne" else e_lz),
            "alpha": float(alpha), "esquema": esquema}


# ----------------------------------------------------------------------
# Alineamiento de potencial (el término ΔV de Freysoldt)
# ----------------------------------------------------------------------
#: Unidades en que pp.x escribe el potencial electrostático (plot_num=11).
UNIDADES_POTENCIAL = {"Ry": qeout.RY_EV, "eV": 1.0}


def alineamiento(pot_defecto: str, pot_perfecto: str, eje: int = 2,
                 fraccion_lejos: float = 0.25,
                 unidades_cube: str = "Ry") -> dict:
    """ΔV entre la supercelda con defecto y la perfecta, lejos del defecto.

    Los dos cálculos no comparten origen de energías: al meter carga y un
    fondo uniforme, todo el potencial se desplaza. Ese desplazamiento entra
    en E_f multiplicado por q, así que hay que medirlo donde el defecto ya
    no se nota — en la zona más lejana de la celda — y no en promedio sobre
    toda ella, que es donde está el defecto.

    pp.x escribe el potencial en Ry (`unidades_cube="Ry"`); aquí se pasa a
    eV, que es la unidad de todo lo demás en E_f ('dV', 'sigma' y 'perfil'
    salen en eV). Si el cube ya viene en eV, pásalo con
    `unidades_cube="eV"`.

    La dispersión de ΔV en esa zona se devuelve también: si no es pequeña,
    la supercelda es demasiado chica y el número no vale.
    """
    from qekit.modules import fields

    if unidades_cube not in UNIDADES_POTENCIAL:
        raise ErrorDeUso(
            f"unidades del potencial '{unidades_cube}' desconocidas; "
            f"opciones: {', '.join(UNIDADES_POTENCIAL)}.")
    factor = UNIDADES_POTENCIAL[unidades_cube]

    cd = fields.read_cube(pot_defecto)
    cp = fields.read_cube(pot_perfecto)
    zd, vd = fields.planar_average(cd, axis=eje)
    zp, vp = fields.planar_average(cp, axis=eje)
    if len(vd) != len(vp):
        vp = np.interp(zd, zp, vp, period=zd[-1] - zd[0] + (zd[1] - zd[0]))
    dv = (np.asarray(vd) - np.asarray(vp)) * factor      # eV

    n = len(dv)
    # el defecto está donde |ΔV| es mayor; la zona buena es la opuesta
    centro = int(np.argmax(np.abs(dv - np.median(dv))))
    ancho = max(3, int(n * fraccion_lejos))
    lejos = [(centro + n // 2 + k - ancho // 2) % n for k in range(ancho)]
    muestra = dv[lejos]
    return {"dV": float(np.mean(muestra)),
            "sigma": float(np.std(muestra)),
            "perfil": dv, "z": np.asarray(zd),
            "indices_lejos": lejos, "unidades": "eV",
            "unidades_cube": unidades_cube}


def electrones(atoms, pseudos: dict):
    """Electrones de valencia de la celda, según los z_valence de los UPF."""
    total = 0.0
    for sym in atoms.get_chemical_symbols():
        z = (pseudos.get(sym) or {}).get("z_valence")
        if z is None:
            return None
        total += float(z)
    return total


# ----------------------------------------------------------------------
# Estructura del cálculo
# ----------------------------------------------------------------------
@dataclass
class DefectRun:
    kind: str = "vacancy"
    cargas: list = field(default_factory=list)
    jobs: list = field(default_factory=list)
    energies: dict = field(default_factory=dict)      # q -> eV
    converged: dict = field(default_factory=dict)
    E_perfecto: float = None
    vbm: float = None
    gap: float = None
    perf_ok: bool = None
    cell: np.ndarray = None
    n_especies: dict = field(default_factory=dict)    # símbolo -> n_i (añadidos)
    mu: dict = field(default_factory=dict)
    epsilon: float = None
    esquema: str = "lany-zunger"
    dV: float = 0.0
    dV_sigma: float = None
    natoms_perf: int = 0
    supercell: tuple = None
    aviso_mu: str = ""

    def correccion(self, q: int) -> dict:
        return correccion_imagen(q, self.cell, self.epsilon, self.esquema)

    def E_f(self, q: int, e_fermi: float = 0.0) -> float:
        """Energía de formación a un nivel de Fermi dado (medido desde el VBM)."""
        e = self.energies.get(q)
        if e is None or self.E_perfecto is None or self.vbm is None:
            return None
        termino_mu = sum(n * self.mu.get(s, 0.0)
                         for s, n in self.n_especies.items())
        return (e - self.E_perfecto - termino_mu
                + q * (self.vbm + e_fermi)
                + self.correccion(q)["E_corr"]
                + q * self.dV)


# ----------------------------------------------------------------------
# Preparación
# ----------------------------------------------------------------------
def prepare(atoms, kind: str = "vacancy", site: int = 0,
            new_element: str = None, supercell=(2, 2, 2), position=None,
            cargas=(0,), outdir: str = "defectos", pseudo_dir: str = None,
            insulator: bool = True, ecutwfc: float = None,
            ecutrho: float = None, kspacing: float = None,
            relax_ions: bool = True, epsilon: float = None,
            esquema: str = "lany-zunger", vdw: str = None,
            nspin: int = 1, magnetization: dict = None) -> tuple:
    """Escribe la supercelda perfecta y una por cada estado de carga."""
    from qekit.modules.inputgen import _estimate_nbnd

    cargas = sorted({int(q) for q in cargas})
    if esquema not in ESQUEMAS:
        raise ErrorDeUso(f"--correction es uno de: {', '.join(ESQUEMAS)}.")
    if any(q != 0 for q in cargas) and esquema != "ninguna" and not epsilon:
        raise ErrorDeUso(
            "hay estados de carga distintos de 0 y no diste --epsilon. La "
            "constante dieléctrica es lo que apantalla la interacción del "
            "defecto con sus imágenes; sin ella la corrección sale ε veces "
            "de más. Puedes obtenerla con 'olla-dft optics' (el límite "
            "ε₁(ω→0)), buscarla en la literatura, o pedir explícitamente "
            "--correction ninguna.")

    perfecto, info = builder.defect(atoms, kind=kind, site=site,
                                    new_element=new_element,
                                    supercell=supercell, position=position)
    defectuoso = info.atoms

    # n_i = átomos AÑADIDOS a la supercelda (negativo si se quitaron)
    n_especies = {}
    if info.especie_ida:
        n_especies[info.especie_ida] = n_especies.get(info.especie_ida, 0) - 1
    if info.especie_nueva:
        n_especies[info.especie_nueva] = n_especies.get(info.especie_nueva, 0) + 1

    # las dos celdas comparten TODO menos los átomos
    conjunto = perfecto.copy()
    for s in set(defectuoso.get_chemical_symbols()):
        if s not in conjunto.get_chemical_symbols():
            from ase import Atoms as _A
            conjunto += _A(s, positions=[[0.0, 0.0, 0.0]])
    common = sweep.prepare_common(conjunto, pseudo_dir, ecutwfc, ecutrho,
                                  insulator,
                                  prefix=perfecto.get_chemical_formula(
                                      mode="hill", empirical=True))
    out = Path(outdir); out.mkdir(parents=True, exist_ok=True)
    grid = sweep.default_grid(perfecto, kspacing)
    nbnd = _estimate_nbnd(perfecto, common["pseudos"])

    run = DefectRun(kind=kind, cargas=cargas,
                    cell=np.asarray(perfecto.cell.array, dtype=float),
                    n_especies=n_especies, epsilon=epsilon, esquema=esquema,
                    natoms_perf=len(perfecto),
                    supercell=tuple(int(x) for x in supercell))

    # ------------------------------------------------------------------
    # Electrones impares y occupations='fixed'.
    #
    # Quitarle un electrón a un aislante deja un número IMPAR de electrones,
    # y con occupations='fixed' pw.x no sabe repartirlos: aborta con "the
    # system is metallic, specify occupations", que suena a que el material
    # es un metal y no lo es. La salida correcta es hacer el cálculo con
    # espín y decirle la magnetización total (1 por cada electrón desparejado),
    # que además es la física correcta: un defecto con un electrón de más o de
    # menos TIENE momento. Se activa para TODOS los estados de carga, no solo
    # los impares, porque si no, las energías no serían comparables entre sí.
    # ------------------------------------------------------------------
    ne_def = electrones(defectuoso, common["pseudos"])
    mags = {}
    aviso_impar = None
    if insulator and ne_def is not None:
        impares = [q for q in cargas if int(round(ne_def - q)) % 2]
        if impares:
            nspin = 2
            for q in cargas:
                mags[q] = float(int(round(ne_def - q)) % 2)
            aviso_impar = (
                "Estados de carga con un número impar de electrones: "
                + ", ".join(etiqueta_q(q) for q in impares)
                + ".\n  Con occupations='fixed' pw.x no puede repartirlos y "
                  "aborta diciendo que el\n  sistema es metálico. Se activa "
                  "el cálculo con espín (nspin=2) en TODOS los\n  estados, "
                  "con tot_magnetization = 1 en los impares y 0 en los pares: "
                  "es la\n  física correcta y deja las energías comparables "
                  "entre sí.")

    extras = dict(vdw=vdw, nspin=nspin, magnetization=magnetization, nbnd=nbnd)
    ne_perf = electrones(perfecto, common["pseudos"])
    run.jobs.append(sweep.write_scf_job(
        perfecto, common, out / "_perfecto", "supercelda perfecta", grid,
        meta={"papel": "perf"}, calculation="scf",
        tot_magnetization=(float(int(round(ne_perf)) % 2)
                           if (mags and ne_perf is not None) else None),
        **extras))

    for q in cargas:
        etiqueta = f"q{'m' if q < 0 else 'p'}{abs(q)}"
        job = sweep.write_scf_job(
            defectuoso, common, out / etiqueta,
            f"q = {etiqueta_q(q)}", grid,
            meta={"papel": "def", "q": q},
            calculation="relax" if relax_ions else "scf",
            tot_charge=float(q), tot_magnetization=mags.get(q), **extras)
        run.jobs.append(job)

    sweep.write_run_script(run.jobs, out / "run.sh")

    detalle = {"vacancy": f"vacancia de {info.especie_ida}",
               "substitution": f"{info.especie_ida} sustituido por "
                               f"{info.especie_nueva}",
               "interstitial": f"{info.especie_nueva} intersticial"}.get(kind, kind)
    report = ["--- Defectos cargados ---",
              f"Defecto: {detalle}",
              f"Supercelda: {run.supercell[0]}x{run.supercell[1]}x"
              f"{run.supercell[2]}  ({len(perfecto)} átomos perfectos, "
              f"{len(defectuoso)} con el defecto)",
              f"Estados de carga: "
              f"{', '.join(etiqueta_q(q) for q in cargas)}",
              f"Malla k: {grid[0]}x{grid[1]}x{grid[2]}  |  "
              + ("posiciones relajadas en cada carga" if relax_ions
                 else "posiciones fijas")]
    lado = float(min(np.linalg.norm(run.cell, axis=1)))
    if any(q != 0 for q in cargas):
        alpha = constante_madelung(run.cell)
        L = abs(np.linalg.det(run.cell)) ** (1.0 / 3.0)
        report.append(f"Corrección de imagen: {esquema}, α_M = {alpha:.4f} "
                      f"(Ewald sobre esta celda), L = V^(1/3) = {L:.3f} Å"
                      + (f", ε = {epsilon:g}" if epsilon else ""))
        if epsilon:
            q_max = max(abs(q) for q in cargas)
            e1 = correccion_imagen(q_max, run.cell, epsilon, esquema)["E_corr"]
            report.append(f"  Para q = ±{q_max} la corrección vale "
                          f"{e1:.3f} eV. Si ese número es del tamaño de la "
                          f"E_f que\n  esperas, la supercelda es demasiado "
                          f"pequeña para fiarse.")
    if aviso_impar:
        report.append(aviso_impar)
    for w in info.warnings:
        report.append(f"AVISO: {w}")
    if lado < 10.0 and any(q != 0 for q in cargas):
        report.append(
            f"AVISO: {lado:.1f} Å de lado con defectos cargados es poco. El "
            "error de tamaño\n  finito va como q²/L y la corrección solo "
            "quita el término principal.")
    warn = sweep.missing_pseudo_warning(common)
    if warn:
        report.append(warn)
    report += ["", f"{len(run.jobs)} cálculos escritos en '{out.resolve()}'",
               "Córrelos con --run, o a mano con ./run.sh dentro de esa carpeta."]
    return run, "\n".join([r for r in report if r])


# ----------------------------------------------------------------------
# Recolección
# ----------------------------------------------------------------------
def collect(run: DefectRun, results: list = None, mu: dict = None) -> DefectRun:
    por_dir = {str(r.job.directory): r for r in (results or [])}

    def _leer(job):
        r = por_dir.get(str(job.directory))
        if r is not None and r.ok and r.result is not None:
            return r.result
        try:
            return qeout.read_xml(str(job.directory))
        except Exception:                                   # noqa: BLE001
            return None

    run.energies, run.converged = {}, {}
    for job in run.jobs:
        res = _leer(job)
        if job.meta.get("papel") == "perf":
            if res is not None:
                run.E_perfecto = res.total_energy
                run.perf_ok = res.converged
                run.vbm = res.homo
                if res.homo is not None and res.lumo is not None:
                    run.gap = res.lumo - res.homo
            continue
        q = job.meta["q"]
        run.energies[q] = res.total_energy if res else None
        run.converged[q] = res.converged if res else None

    # potencial químico
    run.mu = dict(mu or {})
    faltan = [s for s in run.n_especies if s not in run.mu]
    run.aviso_mu = ", ".join(faltan) if faltan else ""
    return run


def asignar_mu_elemental(run: DefectRun, simbolos_perfectos: list) -> bool:
    """μ = E(perfecto)/N para un cristal ELEMENTAL, si el usuario no lo dio.

    Solo vale si el cristal tiene una sola especie: ahí el reservorio es el
    propio material y μ está fijado sin ambigüedad. En un compuesto NO se
    hace, porque μ de cada especie se mueve entre dos límites y elegir uno
    en silencio daría una E_f con pinta de definitiva que no lo es.
    """
    unicos = set(simbolos_perfectos)
    if len(unicos) != 1 or run.E_perfecto is None or not run.natoms_perf:
        return False
    sym = unicos.pop()
    faltan = [s for s in run.n_especies if s not in run.mu]
    if faltan != [sym]:
        return False
    run.mu[sym] = run.E_perfecto / run.natoms_perf
    run.aviso_mu = ""
    return True


# ----------------------------------------------------------------------
# Niveles de transición
# ----------------------------------------------------------------------
def niveles_transicion(run: DefectRun) -> list:
    """ε(q/q') donde E_f(q) y E_f(q') se cruzan, en eV sobre el VBM.

    Se devuelve una entrada por cada par de estados de carga CONSECUTIVOS
    (ordenados por q), con la bandera `dentro` = el cruce cae dentro del
    gap. No se filtra por la envolvente inferior: un cruce entre dos
    estados que nunca son el más estable también sale en la lista, con
    `dentro` diciendo si cae en el gap. Quien necesite solo los niveles
    observables debe cruzar esta lista con `envolvente`.
    """
    qs = sorted(q for q in run.cargas if run.E_f(q) is not None)
    fuera = []
    for a, b in zip(qs, qs[1:]):
        fa, fb = run.E_f(a, 0.0), run.E_f(b, 0.0)
        if b == a:
            continue
        eps = (fa - fb) / (b - a)
        fuera.append({"q1": a, "q2": b, "eps": float(eps),
                      "dentro": (run.gap is not None and 0.0 <= eps <= run.gap)})
    return fuera


def envolvente(run: DefectRun, e_fermi):
    """E_f más baja en cada ε_F y qué carga la produce."""
    e_fermi = np.atleast_1d(np.asarray(e_fermi, dtype=float))
    qs = [q for q in run.cargas if run.E_f(q) is not None]
    if not qs:
        return None, None
    matriz = np.array([[run.E_f(q, ef) for ef in e_fermi] for q in qs])
    idx = np.argmin(matriz, axis=0)
    return matriz[idx, np.arange(len(e_fermi))], np.array([qs[i] for i in idx])


# ----------------------------------------------------------------------
# Reporte
# ----------------------------------------------------------------------
def _nombre_defecto(run: DefectRun) -> str:
    partes = []
    for s, n in sorted(run.n_especies.items()):
        partes.append(f"{'+' if n > 0 else '−'}{abs(n)} {s}")
    return f"{run.kind} ({', '.join(partes)})" if partes else run.kind


def report(run: DefectRun) -> str:
    if not run.energies:
        raise FaltanDatos(
            "no hay resultados todavía. Corre los cálculos (--run, o ./run.sh "
            "en la carpeta) y vuelve con --collect.")

    L = ["--- Energía de formación de defectos ---",
         f"Defecto: {_nombre_defecto(run)}   |   supercelda "
         f"{run.supercell[0]}x{run.supercell[1]}x{run.supercell[2]} "
         f"({run.natoms_perf} átomos)"]

    if run.E_perfecto is None:
        L += ["", "Falta la energía de la supercelda perfecta: sin ella no hay "
                  "resta que hacer.\n  Está en la carpeta _perfecto."]
        return "\n".join(L)
    L.append(f"E(perfecto) = {run.E_perfecto:.6f} eV")

    if run.vbm is None:
        L += ["", "No pude leer el VBM de la supercelda perfecta. Sin él, el "
                  "término q(ε_VBM + ε_F)\n  no se puede escribir y E_f de "
                  "los estados cargados no está definida.\n  Suele pasar en "
                  "un metal (no hay gap) o si el cálculo no tenía bandas "
                  "vacías."]
        return "\n".join(L)
    L.append(f"VBM = {run.vbm:.4f} eV"
             + (f"   |   gap = {run.gap:.4f} eV" if run.gap else
                "   |   gap: no se pudo leer"))

    if run.aviso_mu:
        L += ["",
              f"FALTA el potencial químico de: {run.aviso_mu}.",
              "  E_f depende linealmente de μ, así que sin él los números de "
              "abajo están\n  desplazados por una constante desconocida: las "
              "DIFERENCIAS entre cargas y\n  los niveles de transición sí "
              "valen, el valor absoluto de E_f no.",
              "  Dáselo con  --mu ELEMENTO=VALOR  (energía por átomo del "
              "reservorio, en eV)."]
    else:
        L.append("Potenciales químicos: "
                 + ", ".join(f"μ({s}) = {v:.4f} eV"
                             for s, v in sorted(run.mu.items())))

    if run.dV:
        L.append(f"Alineamiento de potencial ΔV = {run.dV:+.4f} eV"
                 + (f" (σ = {run.dV_sigma:.4f} eV en la zona lejana)"
                    if run.dV_sigma is not None else "")
                 + f"; entra en E_f como q·ΔV = {run.dV:+.4f} eV por unidad "
                   "de carga\n  (el potencial de pp.x viene en Ry y se pasó "
                   "a eV)")
        if run.dV_sigma is not None and run.dV_sigma > 0.3 * abs(run.dV):
            L.append("  La dispersión es grande comparada con el valor: el "
                     "defecto todavía se\n  nota en la zona 'lejana', o sea "
                     "que la supercelda es pequeña.")

    hay_carga = any(q != 0 for q in run.cargas)
    if hay_carga:
        alpha = constante_madelung(run.cell)
        Lc = abs(np.linalg.det(run.cell)) ** (1.0 / 3.0)
        L.append(f"Corrección de imagen: {run.esquema}, α_M = {alpha:.4f}, "
                 f"L = {Lc:.3f} Å"
                 + (f", ε = {run.epsilon:g}" if run.epsilon else ", SIN apantallar"))
        if run.esquema == "ninguna":
            L.append("  SIN CORREGIR: las E_f de los estados cargados están "
                     "sistemáticamente\n  bajas, y el error crece con q². "
                     "Solo sirven para comparar entre sí\n  cargas del mismo "
                     "valor absoluto.")

    L += ["", f"  {'q':>3s} {'E(defecto)':>14s} {'E_corr':>9s} "
              f"{'E_f(ε_F=0)':>12s}"
          + (f" {'E_f(ε_F=gap)':>13s}" if run.gap else "")]
    L.append("  " + "-" * (56 if run.gap else 42))
    for q in run.cargas:
        e = run.energies.get(q)
        if e is None:
            L.append(f"  {etiqueta_q(q):>3s} {'sin resultado':>14s}")
            continue
        c = run.correccion(q)["E_corr"]
        f0 = run.E_f(q, 0.0)
        fila = f"  {etiqueta_q(q):>3s} {e:>14.6f} {c:>9.4f} {f0:>12.4f}"
        if run.gap:
            fila += f" {run.E_f(q, run.gap):>13.4f}"
        if run.converged.get(q) is False:
            fila += "   << SIN CONVERGER"
        L.append(fila)

    trans = niveles_transicion(run)
    if trans:
        L += ["", "Niveles de transición de carga (eV sobre el VBM):"]
        for t in trans:
            marca = "" if t["dentro"] else "   << fuera del gap"
            L.append(f"  ε({etiqueta_q(t['q1'])}/{etiqueta_q(t['q2'])}) = "
                     f"{t['eps']:8.4f}{marca}")
        fuera = [t for t in trans if not t["dentro"]]
        if fuera and run.gap:
            L.append("  Un nivel fuera del gap no es un estado observable: "
                     "quiere decir que esa\n  carga nunca llega a ser la más "
                     "estable dentro del gap.")

    if run.gap:
        ef = np.linspace(0.0, run.gap, 201)
        env, qs = envolvente(run, ef)
        if env is not None:
            estables = []
            for q in dict.fromkeys(qs):
                estables.append(int(q))
            L += ["", "Cargas estables al recorrer el gap: "
                      + " → ".join(etiqueta_q(q) for q in estables)]
            L.append(f"  E_f mínima en el gap: {float(np.min(env)):.4f} eV")

    if run.perf_ok is False:
        L.append("\nAVISO: la supercelda perfecta no convergió; toda la "
                 "columna hereda ese error.")
    return "\n".join(L)


def export(run: DefectRun, outdir: str = ".") -> list:
    out = Path(outdir); out.mkdir(parents=True, exist_ok=True)
    f = out / "FORMACION.dat"
    lines = [provenance.header(
        f"energia de formacion, {_nombre_defecto(run)}",
        {"supercelda": "x".join(str(v) for v in run.supercell),
         "E_perfecto_eV": run.E_perfecto, "VBM_eV": run.vbm,
         "gap_eV": run.gap, "epsilon": run.epsilon, "esquema": run.esquema,
         "dV_eV": run.dV}),
        f"# {'q':>4s} {'E_def(eV)':>16s} {'E_corr(eV)':>12s} "
        f"{'E_f(eF=0)':>12s}"]
    for q in run.cargas:
        if run.energies.get(q) is None:
            continue
        lines.append(f"{q:6d} {run.energies[q]:16.8f} "
                     f"{run.correccion(q)['E_corr']:12.5f} "
                     f"{run.E_f(q, 0.0):12.5f}")
    if run.gap:
        lines += ["", "# E_f(eF) en el gap: eF  " +
                  "  ".join(f"q={q:+d}" for q in run.cargas)]
        for ef in np.linspace(0.0, run.gap, 51):
            vals = [run.E_f(q, ef) for q in run.cargas]
            lines.append(f"{ef:12.5f} " + " ".join(
                f"{v:12.5f}" if v is not None else f"{'nan':>12s}" for v in vals))
    f.write_text("\n".join(lines) + "\n")
    txt = out / "FORMACION.txt"
    txt.write_text(report(run) + "\n")
    return [str(f), str(txt)]


def plot(run: DefectRun, outfile: str = "formacion", formats="pdf,png",
         theme: str = None, size: str = None, family: str = None,
         background: str = None, palette=None, usetex: bool = None,
         width="single", journal: str = "generic", aspect: float = 0.78,
         mono: bool = False, dpi: int = None) -> list:
    """El diagrama clásico: E_f contra ε_F, una recta por carga."""
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError as exc:                              # pragma: no cover
        raise RuntimeError("matplotlib no está instalado.") from exc
    if run.gap is None or run.vbm is None:
        raise FaltanDatos(
            "para el diagrama hace falta el gap de la supercelda perfecta, y "
            "no se pudo leer. En un metal este diagrama no tiene sentido.")
    qs = [q for q in run.cargas if run.E_f(q) is not None]
    if not qs:
        raise FaltanDatos("no hay energías de formación que graficar.")

    st = qstyle.apply(theme, size=size, family=family, background=background,
                      palette=palette, usetex=usetex, mono=mono)
    fig, ax = qstyle.new_figure(width, journal, aspect)
    cols = qstyle.palette(max(len(qs), 3), mono=mono)
    ef = np.linspace(0.0, run.gap, 401)

    for i, q in enumerate(qs):
        y = [run.E_f(q, e) for e in ef]
        ax.plot(ef, y, lw=st["line"] * 0.8, color=cols[i % len(cols)],
                dashes=[3.0, 2.0], label=f"q = {etiqueta_q(q)}")

    env, qenv = envolvente(run, ef)
    if env is not None:
        ax.plot(ef, env, lw=st["line"] * 2.0, color=qstyle.INK_SOFT,
                solid_capstyle="round", zorder=1, alpha=0.85)

    for t in niveles_transicion(run):
        if not t["dentro"]:
            continue
        ax.axvline(t["eps"], color=qstyle.INK_FAINT, lw=st["axis_line"],
                   dashes=[1.5, 1.5], zorder=0)
        y_codo = run.E_f(t["q1"], t["eps"])
        ax.annotate(f"({etiqueta_q(t['q1'])}/{etiqueta_q(t['q2'])})",
                    xy=(t["eps"], y_codo), xytext=(4, -12),
                    textcoords="offset points", fontsize=st["legend"],
                    color=qstyle.INK_SOFT)

    ax.set_xlim(0.0, run.gap)
    ax.set_xlabel(r"$\varepsilon_F$ sobre el VBM (eV)")
    ax.set_ylabel(r"$E_f$ (eV)")
    ax.legend(frameon=False, fontsize=st["legend"], ncol=2)
    written = qstyle.save(fig, outfile, formats, dpi=dpi, modulo="formacion")
    plt.close(fig)
    return written
