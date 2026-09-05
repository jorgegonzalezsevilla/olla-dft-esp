# Olla-DFT — command-line toolkit for Quantum ESPRESSO
# Copyright (C) 2026 Jorge Enrique González Sevilla
# SPDX-License-Identifier: AGPL-3.0-or-later
# This program is free software: you can redistribute it and/or modify it
# under the terms of the GNU Affero General Public License as published by the Free
# Software Foundation, either version 3 of the License, or (at your option)
# any later version. See the LICENSE file for details.

"""Superficies cargadas: el medio de apantallamiento efectivo (ESM).

Una losa en una celda periódica está rodeada de copias de sí misma también
en la dirección del vacío. Mientras la losa sea neutra y simétrica eso solo
cuesta vacío de más. En cuanto deja de serlo —una superficie con dipolo, o
peor, una superficie CARGADA— el cálculo periódico deja de describir lo que
se quiere describir:

  - Con dipolo, el campo de las imágenes no se apaga: hace falta mucho vacío,
    y aun así queda un error que decae despacio.
  - Con carga neta, Quantum ESPRESSO añade un fondo uniforme compensador que
    llena TODO el volumen, vacío incluido. Ese fondo no existe en ningún
    experimento, y la energía que sale depende del tamaño de la celda de una
    forma que no converge a nada.

El ESM (Otani y Sugino, PRB 73, 115407) sustituye las imágenes en z por una
condición de contorno explícita: se resuelve la ecuación de Poisson en la
celda y se empalma con una solución analítica fuera. Tres opciones, y elegir
mal es el error más común:

  - **bc1: vacío / vacío.** Para losas NEUTRAS. El nivel de vacío queda
    fijado a cero por construcción, así que la función trabajo es
    directamente −E_F y no hay que ajustar ninguna meseta. Y deja de
    depender del vacío: en Al(111) la energía y E_F salen iguales con 8, 12
    y 16 Å (0.4 meV de diferencia). En un portátil eso es menos de la mitad
    de celda.
  - **bc2: metal / metal.** Un condensador. Admite `esm_efield`.
  - **bc3: vacío / metal.** El caso del electrodo: la contracarga se va al
    metal. Es el único que tiene sentido con carga neta junto con bc2.

Y la consecuencia que hay que entender antes de usarlo con carga: **con bc2
o bc3 el vacío deja de ser un parámetro de convergencia y pasa a ser física**.
Es la distancia al contraelectrodo, así que la energía y el potencial crecen
LINEALMENTE con él —medido: 0.1498 Ry y 20.4 eV por cada 4 Å en Al(111)—.
No hay que converger esa distancia: hay que elegirla, y decir cuál es.

Dos trampas más, las dos silenciosas:

  1. **ESM mide z desde el CENTRO de la celda.** Una losa centrada donde la
     deja ASE, en c/2, cae justo sobre la frontera de ESM, y lo que sale son
     energías de cientos de Ry sin ningún aviso. Aquí se centra en z = 0
     antes de escribir nada.
  2. **bc1 con carga neta no es un problema bien planteado.** Una losa
     cargada rodeada de vacío por los dos lados tiene un campo que llega al
     infinito y una energía divergente. QE calcula algo igualmente. Aquí se
     rechaza.
"""

from dataclasses import dataclass, field
from pathlib import Path
import re

import numpy as np

from qekit.core.errors import ErrorDeUso, FaltanDatos
from qekit.core import style as qstyle

BC = {
    "bc1": "vacío / vacío — losas neutras",
    "bc2": "metal / metal — condensador, admite campo aplicado",
    "bc3": "vacío / metal — electrodo, para superficies cargadas",
}
# margen (Å) desde el borde de la losa al empezar a promediar el vacío
MARGEN_VACIO = 2.0
# vacío mínimo aceptable: con ESM no hace falta más
VACIO_MINIMO = 6.0
# e/Å² a µF/cm²:  1.602176634e-19 C / (1e-8 cm)² = 1.602e-3 C/cm² = 1.602e3 µC/cm²
E_A2_A_UF_CM2 = 1.602176634e3


def centrar(atoms, eje=2):
    """Deja la losa centrada en z = 0, que es donde ESM la espera.

    Es la corrección que más cambia el resultado y la que nadie ve: con la
    losa en c/2 —lo que hace `center()` de ASE— ESM la parte por la mitad
    contra su propia frontera. La energía sale mal por cientos de Ry y el
    cálculo termina normalmente.
    """
    a = atoms.copy()
    z = a.get_positions()[:, int(eje)]
    a.positions[:, int(eje)] -= 0.5 * (z.min() + z.max())
    return a


def espesor_y_vacio(atoms, eje=2):
    """Grosor de la losa y vacío que le queda, en Å."""
    z = atoms.get_positions()[:, int(eje)]
    c = float(np.linalg.norm(atoms.cell.array[int(eje)]))
    esp = float(z.max() - z.min())
    return esp, c - esp


def comprobar(atoms, bc="bc1", cargas=(0.0,), eje=2):
    """Todo lo que tiene que cumplirse antes de escribir un input de ESM."""
    avisos = []
    bc = str(bc).lower()
    if bc not in BC:
        raise ErrorDeUso(
            f"condición de contorno '{bc}' desconocida. Las que hay:\n  "
            + "\n  ".join(f"{k}: {v}" for k, v in BC.items()))
    cargas = [float(q) for q in cargas]
    if bc == "bc1" and any(abs(q) > 1e-12 for q in cargas):
        raise ErrorDeUso(
            "bc1 es vacío por los dos lados, y una losa cargada rodeada de "
            "vacío tiene un campo\nque llega al infinito: la energía diverge "
            "y el problema no está bien planteado.\npw.x calcula algo de "
            "todas formas (salieron −379 y −677 Ry para la misma losa con "
            "dos\nvacíos distintos). Para una superficie cargada usa "
            "--bc bc3, que pone un electrodo\nmetálico al otro lado y le da "
            "la contracarga.")

    cell = np.array(atoms.cell.array, float)
    otros = [i for i in range(3) if i != int(eje)]
    for i in otros:
        if abs(cell[i, int(eje)]) > 1e-6:
            raise ErrorDeUso(
                "ESM exige que la celda sea ortogonal en z: los vectores a y "
                "b tienen que estar\nen el plano xy. El tuyo tiene componente "
                "z. Reorienta la celda antes.")
    esp, vac = espesor_y_vacio(atoms, eje)
    if vac < VACIO_MINIMO:
        avisos.append(
            f"Solo quedan {vac:.1f} Å de vacío. Con ESM no hace falta mucho, "
            f"pero sí lo bastante\n  para que la densidad electrónica llegue "
            f"a cero antes de la frontera: {VACIO_MINIMO:.0f} Å es el mínimo "
            f"razonable.")
    z = atoms.get_positions()[:, int(eje)]
    if abs(0.5 * (z.min() + z.max())) > 0.05:
        avisos.append(
            "La losa no estaba centrada en z = 0 y se ha centrado. ESM mide z "
            "desde el CENTRO\n  de la celda, no desde el origen: una losa "
            "dejada en c/2 (lo que hace ASE) cae\n  sobre la frontera de ESM "
            "y el resultado es basura sin ningún mensaje de error.")
    if bc in ("bc2", "bc3") and any(abs(q) > 1e-12 for q in cargas):
        avisos.append(
            f"Con {bc} y carga neta, el vacío NO es un parámetro de "
            f"convergencia: es la distancia\n  al contraelectrodo, y entra en "
            f"la respuesta. La energía y E_F crecen linealmente\n  con él "
            f"(medido en Al(111): 0.15 Ry y 20 eV por cada 4 Å). Elige la "
            f"distancia y\n  dila; no la 'converjas'.")
    return avisos


def leer_esm1(ruta):
    """El archivo prefix.esm1: z, carga y potenciales promediados en el plano."""
    p = Path(ruta)
    if p.is_dir():
        cand = sorted(p.rglob("*.esm1"))
        if not cand:
            raise FaltanDatos(
                f"en {p} no hay ningún archivo .esm1. Lo escribe pw.x cuando "
                f"corre con\n  assume_isolated = 'esm'; si no está, el "
                f"cálculo no usó ESM.")
        p = cand[0]
    d = np.loadtxt(p)
    return {"z": d[:, 0], "carga": d[:, 1], "v_hartree": d[:, 2],
            "v_local": d[:, 3], "v_total": d[:, 4], "archivo": str(p)}


def nivel_vacio(perfil, espesor=None, margen=MARGEN_VACIO, lado=None,
                tol=1e-3, margen_max=8.0):
    """Potencial medio en la región de vacío, y su dispersión.

    La región se elige por GEOMETRÍA (más allá del borde de la losa) y no por
    dónde hay carga: con bc3 la contracarga se acumula en la frontera
    metálica, así que un criterio basado en la carga se queda sin puntos justo
    en el caso que interesa.

    El margen crece hasta que el potencial sea de verdad PLANO (desviación
    típica por debajo de `tol`). La densidad electrónica se derrama unos
    angstroms más allá de los núcleos, y promediar demasiado cerca de la losa
    devuelve el potencial de la cola, no el del vacío: con 2 Å de margen sale
    −0.026 eV en Al(111), con 5 Å sale −5·10⁻⁵. Es un error de 26 meV en la
    función trabajo, silencioso, y del tamaño de las diferencias que se suelen
    discutir.

    `lado` = −1 o +1 restringe el promedio a un lado, que es lo que hay que
    hacer con bc3: el lado del metal no tiene nivel de vacío que promediar.
    """
    z, v = perfil["z"], perfil["v_total"]
    if espesor is None:
        chg = np.abs(perfil["carga"])
        espesor = (2 * np.abs(z[chg > 0.02 * chg.max()]).max()
                   if chg.max() > 0 else 0.0)
    lim = float(np.abs(z).max()) - espesor / 2.0
    tope = min(float(margen_max), max(float(margen), lim - 0.5))
    mejor = None
    usado = float(margen)
    while usado <= tope + 1e-9:
        m = np.abs(z) > espesor / 2.0 + usado
        if lado is not None:
            m &= (np.sign(z) == np.sign(lado))
        if m.sum() >= 5:
            mejor = m
            if float(v[m].std()) < tol:
                break
        usado += 0.5
    m = mejor
    if m is None:
        raise FaltanDatos(
            "no queda región de vacío suficiente para promediar el "
            "potencial: la losa llena\n  casi toda la celda. Con ESM basta "
            "poco vacío, pero hace falta ALGO: unos 6 Å.")
    return float(v[m].mean()), float(v[m].std()), int(m.sum())


def funcion_trabajo(perfil, fermi, espesor=None, lado=None):
    """Φ = V_vacío − E_F. Con bc1, V_vacío = 0 por construcción."""
    v, _s, _n = nivel_vacio(perfil, espesor=espesor, lado=lado)
    return float(v - float(fermi))


def linealidad(cargas, phi, tol=0.02):
    """¿Es Φ(q) una recta? Si no lo es, no hay capacitancia que reportar.

    En un condensador ideal el potencial es lineal en la carga. Cuando deja
    de serlo con ESM suele ser porque los electrones empiezan a derramarse
    sobre el electrodo de la frontera, y a partir de ahí ni Φ ni la energía
    describen ya la superficie que se quería describir. Medido en Al(111)
    con bc3: con ±0.15 e la pendiente entre puntos consecutivos cambia un
    factor 3, así que ahí ya no hay recta.
    """
    q = np.asarray(cargas, float)
    P = np.asarray(phi, float)
    if len(q) < 3:
        return None, None
    coef = np.polyfit(q, P, 1)
    aj = np.polyval(coef, q)
    rango = float(P.max() - P.min())
    if rango <= 0:
        return 1.0, 0.0
    desv = float(np.abs(P - aj).max()) / rango
    return (desv <= tol), desv


def capacitancia(cargas, potenciales, area):
    """C = dq/dΦ, en µF/cm². `area` en Å².

    Es la pendiente de la curva carga-potencial. Para una interfaz metal/vacío
    con un electrodo a unos pocos Å salen valores de decenas de µF/cm²; una
    doble capa electroquímica real está en 10-100. Si sale algo de otro orden,
    lo primero que hay que mirar es la distancia al contraelectrodo, porque
    con bc2/bc3 es ella la que la fija.
    """
    q = np.asarray(cargas, float)
    V = np.asarray(potenciales, float)
    if len(q) < 2:
        return None, None
    coef = np.polyfit(V, q, 1)
    r2 = None
    if len(q) > 2:
        aj = np.polyval(coef, V)
        ss = float(((q - q.mean()) ** 2).sum())
        r2 = 1.0 - float(((q - aj) ** 2).sum()) / ss if ss > 0 else 1.0
    return float(coef[0]) / float(area) * E_A2_A_UF_CM2, r2


def potencial_de_carga_cero(cargas, potenciales):
    """El Φ al que la superficie no tiene carga neta. Se interpola."""
    q = np.asarray(cargas, float)
    V = np.asarray(potenciales, float)
    if len(q) < 2:
        return None
    orden = np.argsort(q)
    return float(np.interp(0.0, q[orden], V[orden]))


def gran_canonico(energias, cargas, potencial):
    """Ω = E + q·Φ: la energía a potencial fijo en vez de a carga fija.

    Un electrodo real está a potencial constante, no a carga constante: es el
    potenciostato el que fija Φ y la carga se ajusta. Comparar energías de
    losas con cargas distintas SIN esta transformación es comparar sistemas
    con distinto número de electrones, que no es comparar nada.
    """
    E = np.asarray(energias, float)
    q = np.asarray(cargas, float)
    return E + q * float(potencial)


@dataclass
class EsmRun:
    """Un barrido de carga con ESM."""
    formula: str = ""
    bc: str = "bc1"
    cargas: list = field(default_factory=list)
    campo: float = 0.0
    espesor: float = 0.0
    vacio: float = 0.0
    area: float = 0.0
    jobs: list = field(default_factory=list)
    energias: list = field(default_factory=list)      # eV
    fermis: list = field(default_factory=list)        # eV
    perfiles: list = field(default_factory=list)
    vac: list = field(default_factory=list)           # nivel de vacío, eV
    phi: list = field(default_factory=list)           # función trabajo, eV
    avisos: list = field(default_factory=list)


def prepare(atoms, outdir="esm", bc="bc1", cargas=(0.0,), campo=0.0,
            esm_w=0.0, nfit=4, pseudo_dir=None, ecutwfc=None, ecutrho=None,
            kspacing=None, degauss=0.02, smearing="mv"):
    """Un scf por carga, todos con la losa centrada y ESM bien puesto."""
    from qekit.modules import inputgen, sweep
    from qekit.core.runner import Job

    cargas = [float(q) for q in cargas]
    avisos = comprobar(atoms, bc, cargas)
    at = centrar(atoms)
    esp, vac = espesor_y_vacio(at)
    cell = np.array(at.cell.array, float)
    area = float(abs(np.cross(cell[0], cell[1])[2]))

    run = EsmRun(formula=at.get_chemical_formula(), bc=str(bc).lower(),
                 cargas=cargas, campo=float(campo), espesor=esp, vacio=vac,
                 area=area, avisos=avisos)
    common = sweep.prepare_common(at, pseudo_dir, ecutwfc, ecutrho, False)
    out = Path(outdir); out.mkdir(parents=True, exist_ok=True)
    grid = sweep.default_grid(at, kspacing)
    grid = (grid[0], grid[1], 1)          # nunca k a lo largo del vacío

    for i, q in enumerate(cargas):
        d = out / f"q{i:02d}"
        d.mkdir(parents=True, exist_ok=True)
        txt = inputgen.build_pw_input(
            atoms=at, pseudos=common["pseudos"], calculation="scf",
            prefix=common["prefix"], pseudo_dir=common["pseudo_dir"],
            ecutwfc=common["ecutwfc"], ecutrho=common["ecutrho"],
            kcard=f"K_POINTS automatic\n  {grid[0]} {grid[1]} 1 0 0 0\n",
            insulator=False, degauss=degauss, smearing=smearing,
            conv_thr=1e-8, tot_charge=(q if abs(q) > 1e-12 else None))
        extra = ("  assume_isolated = 'esm'\n"
                 f"  esm_bc          = '{run.bc}'\n"
                 f"  esm_nfit        = {int(nfit)}\n")
        if abs(float(esm_w)) > 0:
            extra += f"  esm_w           = {float(esm_w):.6f}\n"
        if run.bc == "bc2" and abs(float(campo)) > 0:
            extra += f"  esm_efield      = {float(campo):.8f}\n"
        txt = re.sub(r"(&SYSTEM\n)", r"\1" + extra, txt, count=1)
        sweep.write_input(d / "pw.in", txt)
        run.jobs.append(Job(name=f"q = {q:+.3f} e", directory=d,
                            meta={"carga": q}))

    sweep.write_run_script(run.jobs, out / "run.sh")
    rep = [f"--- Superficie con ESM: {run.formula} ---",
           f"Condición de contorno: {run.bc} — {BC[run.bc]}",
           f"Losa de {esp:.2f} Å con {vac:.2f} Å de vacío   |   "
           f"área {area:.3f} Å²",
           f"Cargas: {', '.join(f'{q:+g}' for q in cargas)} e",
           ""]
    if run.bc == "bc1":
        rep += ["Con bc1 el nivel de vacío vale CERO por construcción: la "
                "función trabajo es",
                "directamente −E_F, sin ajustar ninguna meseta. Y deja de "
                "depender del vacío:",
                "en Al(111) la energía cambia 6·10⁻⁶ Ry entre 8 y 16 Å, así "
                "que se puede usar",
                "media celda y ahorrarse la mitad del cálculo.", ""]
    rep += [f"Archivos en '{out.resolve()}':",
            f"  q00..q{len(cargas) - 1:02d}/   un scf por carga",
            "  ./run.sh        los lanza todos",
            "",
            f"Luego: olla-dft esm <losa> --collect -o {out}"]
    warn = sweep.missing_pseudo_warning(common)
    if warn:
        rep.append(warn)
    return run, common, "\n".join(rep)


def collect(run, outdir="esm"):
    """Lee energías, E_F y el perfil .esm1 de cada carga."""
    from qekit.core import qeout
    out = Path(outdir)
    dirs = [Path(j.directory) for j in run.jobs] if run.jobs else \
        sorted(p for p in out.glob("q[0-9][0-9]") if p.is_dir())
    if not dirs:
        raise FaltanDatos(f"en {out} no hay carpetas q00, q01...")
    run.energias, run.fermis, run.perfiles, run.vac, run.phi = [], [], [], [], []
    for d in dirs:
        try:
            res = qeout.read_xml(str(d))
        except Exception as exc:                            # noqa: BLE001
            raise FaltanDatos(
                f"no puedo leer el resultado de {d.name}: {exc}.\n  ¿Corriste "
                f"los cálculos? `bash run.sh` dentro de la carpeta, o el "
                f"mismo comando con --run.") from None
        run.energias.append(float(res.total_energy)
                            if res.total_energy is not None else float("nan"))
        run.fermis.append(float(res.fermi) if res.fermi is not None
                          else float("nan"))
        perfil = leer_esm1(d)
        run.perfiles.append(perfil)
        # con bc3 el lado +z es metal: ahí no hay nivel de vacío que promediar
        lado = -1 if run.bc == "bc3" else None
        v, _s, _n = nivel_vacio(perfil, espesor=run.espesor, lado=lado)
        run.vac.append(v)
        run.phi.append(v - run.fermis[-1])
    if len(run.cargas) != len(dirs):
        run.cargas = [float(j.meta.get("carga", i))
                      for i, j in enumerate(run.jobs)] or \
            list(range(len(dirs)))
    return run


def report(run) -> str:
    L = [f"--- Superficie con ESM: {run.formula} ---",
         f"Condición de contorno: {run.bc} — {BC.get(run.bc, '')}",
         f"Losa de {run.espesor:.2f} Å con {run.vacio:.2f} Å de vacío   |   "
         f"área {run.area:.3f} Å²"]
    if run.bc == "bc2" and run.campo:
        L.append(f"Campo aplicado: {run.campo:g} Ry/u.a.")
    if not run.energias:
        L += ["", "Todavía no hay resultados."]
    else:
        L += ["", "   q (e)     E (eV)        E_F (eV)   V_vac (eV)   "
                  "Φ (eV)"]
        for q, e, f, v, p in zip(run.cargas, run.energias, run.fermis,
                                 run.vac, run.phi):
            L.append(f"  {q:+7.3f} {e:14.5f} {f:11.4f} {v:12.5f} "
                     f"{p:9.4f}")
        if run.bc == "bc1":
            peor = max(abs(v) for v in run.vac)
            L += ["", f"Nivel de vacío: |V| ≤ {peor:.1e} eV"
                      + ("  ✓ (bc1 lo fija a cero; que salga cero quiere "
                         "decir que la losa no toca la frontera)"
                         if peor < 1e-3 else
                         "  ← no es cero: la losa está tocando la frontera "
                         "de ESM")]
        if len(run.cargas) > 1:
            # La capacitancia sale del VOLTAJE de la celda (el nivel de vacío
            # respecto de la frontera de ESM), no de la función trabajo. Φ
            # mezcla el voltaje con el cambio del dipolo de superficie y no
            # tiene por qué ser lineal; V_vac sí: en Al(111) con bc3 sale
            # lineal a 4 cifras mientras Φ se desvía un 16 %.
            okv, desvv = linealidad(run.cargas, run.vac)
            if okv:
                Cv, r2v = capacitancia(run.cargas, run.vac, run.area)
                L += ["", "Voltaje de la celda (nivel de vacío respecto de "
                          "la frontera de ESM):",
                      f"  C = dq/dV = {abs(Cv):.3f} µF/cm²"
                      + (f"   (ajuste lineal, R² = {r2v:.6f})"
                         if r2v is not None else ""),
                      "  Es la capacitancia DE ESTE MONTAJE —losa más hueco "
                      "hasta el electrodo—, no una",
                      "  propiedad del material: cambia si cambias el vacío. "
                      "Compárala con ε₀/d antes",
                      "  de creértela."]
            ok, desv = linealidad(run.cargas, run.phi)
            if ok:
                C, r2 = capacitancia(run.cargas, run.phi, run.area)
                pzc = potencial_de_carga_cero(run.cargas, run.phi)
                L += ["", f"Capacitancia  C = dq/dΦ = {C:.2f} µF/cm²"
                          + (f"   (ajuste lineal, R² = {r2:.6f})"
                             if r2 is not None else ""),
                      f"Potencial de carga cero: Φ = {pzc:.4f} eV",
                      "  Con bc2 o bc3 la capacitancia depende de la "
                      "distancia al contraelectrodo:",
                      "  es una capacitancia DE ESTE MONTAJE, no una "
                      "propiedad del material."]
            elif ok is not None:
                run.avisos.append(
                    f"Φ(q) = V_vac − E_F NO es una recta: se desvía un "
                    f"{desv * 100:.0f} % del ajuste lineal, así que\n  no doy "
                    f"un potencial de carga cero sobre ella. Φ mezcla el "
                    f"voltaje de la celda con el\n  cambio del dipolo de "
                    f"superficie al cargarla, y no tiene por qué ser lineal; "
                    f"la que sí\n  lo es, y de la que sale la capacitancia, "
                    f"es V_vac.")
            run.avisos.append(
                "La capacitancia de arriba es la de ESTE montaje. Se "
                "comprobó que obedece la ley\n  del condensador plano: en "
                "Al(111), 1/C frente a la distancia al electrodo sale una "
                "recta\n  de pendiente 1/ε₀ con un 0.4 % de error y "
                "R² = 0.99998 entre 4 y 11 Å. Eso valida la\n  fórmula y la "
                "condición de contorno; NO valida que tu montaje represente "
                "la\n  interfaz que te interesa, que es otra cosa y "
                "depende de ti.")
            L += ["", "Sobre las ENERGÍAS de estas filas: con carga neta no "
                      "son comparables entre sí.",
                  "Cada una tiene distinto número de electrones, y además la "
                  "energía de ESM incluye",
                  "la interacción con la carga imagen del electrodo, que "
                  "crece como q². Este módulo",
                  "da el perfil de potencial y Φ(q), que sí están bien "
                  "definidos; no lo uses para",
                  "restar energías de superficies con cargas distintas."]
    for a in run.avisos:
        L += ["", f"AVISO: {a}"]
    return "\n".join(L)


def export(run, outdir="esm") -> list:
    out = Path(outdir); out.mkdir(parents=True, exist_ok=True)
    escritos = []
    if run.energias:
        f = out / "ESM.dat"
        np.savetxt(f, np.column_stack([run.cargas, run.energias, run.fermis,
                                       run.vac, run.phi]), fmt="%14.6f",
                   header="q(e)   E(eV)   E_F(eV)   V_vac(eV)   Phi(eV)")
        escritos.append(str(f))
    for i, p in enumerate(run.perfiles):
        f = out / f"ESM_perfil_q{i:02d}.dat"
        np.savetxt(f, np.column_stack([p["z"], p["carga"], p["v_total"]]),
                   fmt="%14.6f", header="z(A)   carga(e/A)   v_tot(eV)")
        escritos.append(str(f))
    f = out / "ESM.txt"
    f.write_text(report(run) + "\n", encoding="utf-8")
    escritos.append(str(f))
    return escritos


def plot(run, outfile="esm", formats="pdf,png", theme=None, size=None,
         family=None, background=None, palette=None, usetex=None,
         width="single", journal="generic", aspect=0.72, mono=False,
         dpi=None) -> list:
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError as exc:                              # pragma: no cover
        raise RuntimeError("matplotlib no está instalado.") from exc
    if not run.perfiles:
        raise FaltanDatos("no hay perfiles que dibujar.")
    st = qstyle.apply(theme, size=size, family=family, background=background,
                      palette=palette, usetex=usetex, mono=mono)
    cols = qstyle.palette(max(3, len(run.perfiles)), mono=mono)
    fig, ax = qstyle.new_figure(width, journal, aspect)
    for i, (p, q) in enumerate(zip(run.perfiles, run.cargas)):
        ax.plot(p["z"], p["v_total"], lw=st["line"], color=cols[i % len(cols)],
                label=f"q = {q:+g} e")
    zs = run.perfiles[0]["z"]
    ax.axvspan(-run.espesor / 2, run.espesor / 2, color="0.5", alpha=0.12,
               lw=0, zorder=0)
    ax.set_xlabel("z (Å)   —   la losa está en la banda gris")
    ax.set_ylabel("potencial promediado en el plano (eV)")
    ax.set_xlim(float(zs.min()), float(zs.max()))
    if len(run.perfiles) > 1:
        ax.legend(frameon=False, fontsize=st["legend"])
    escritos = qstyle.save(fig, outfile, formats, dpi=dpi, modulo="esm")
    plt.close(fig)

    if len(run.cargas) > 1:
        fig2, ax2 = qstyle.new_figure(width, journal, aspect)
        ax2.plot(run.phi, run.cargas, marker="o", ms=4, lw=st["line"],
                 color=cols[0])
        pzc = potencial_de_carga_cero(run.cargas, run.phi)
        if pzc is not None:
            ax2.axvline(pzc, color=cols[1], lw=0.8, dashes=[3.0, 2.0])
            ax2.axhline(0.0, color=cols[1], lw=0.8, dashes=[3.0, 2.0])
            ax2.annotate(f"PZC = {pzc:.3f} eV", (pzc, 0),
                         textcoords="offset points", xytext=(6, 6),
                         fontsize=st["legend"], color=cols[1])
        ax2.set_xlabel(r"$\Phi$ (eV)")
        ax2.set_ylabel("q (e por celda)")
        escritos += qstyle.save(fig2, str(outfile) + "_carga", formats,
                                dpi=dpi, modulo="esm")
        plt.close(fig2)
    return escritos
