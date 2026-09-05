# Olla-DFT — command-line toolkit for Quantum ESPRESSO
# Copyright (C) 2026 Jorge Enrique González Sevilla
# SPDX-License-Identifier: AGPL-3.0-or-later
# This program is free software: you can redistribute it and/or modify it
# under the terms of the GNU Affero General Public License as published by the Free
# Software Foundation, either version 3 of the License, or (at your option)
# any later version. See the LICENSE file for details.

"""Sitios de adsorción sobre una losa: enumerarlos, montarlos y compararlos.

La energía de adsorción es una resta de tres números:

    E_ads = E(losa + molécula) - E(losa) - n·E(molécula)

La resta ya vivía en `thermochem.adsorcion()`. Lo que faltaba, y es lo que
de verdad cuesta a mano, es todo lo de alrededor: encontrar los sitios de
la superficie, no repetir los que la simetría hace equivalentes, colocar la
molécula a una altura sensata en cada uno, y garantizar que los tres
cálculos son comparables. Ese último punto es el que más silenciosamente se
rompe: si la losa limpia y la losa con molécula no comparten celda, cutoff,
malla k y pseudos, la resta da un número perfectamente formado que no
significa nada. Aquí los tres cálculos se generan a la vez y de la misma
plantilla, precisamente para que no puedan divergir.

Tipos de sitio que se enumeran:
  top      encima de un átomo de la superficie
  bridge   sobre el punto medio de dos átomos vecinos
  hollow   sobre el centro de un triángulo de átomos (fcc, hcp, 4-fold...)
"""

from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

from qekit.core import provenance, qeout, structure
from qekit.core.errors import ErrorDeUso, FaltanDatos
from qekit.core import style as qstyle
from qekit.modules import sweep, thermochem

TIPOS = ("top", "bridge", "hollow")

# Distancia máxima (Å) para considerar vecinos dos átomos de la superficie al
# formar los puentes. Por encima de esto el "puente" no está entre nada.
R_VECINO = 3.6

# Tolerancia (Å) para decidir qué átomos forman la capa superior.
TOL_CAPA = 0.6


@dataclass
class Sitio:
    tipo: str
    xy: tuple                   # posición cartesiana en el plano (Å)
    z: float                    # altura de la superficie bajo el sitio (Å)
    huella: object = None       # firma de distancias, para deduplicar
    rotacion: float = 0.0       # grados alrededor de la normal
    etiqueta: str = ""


@dataclass
class AdsorbRun:
    jobs: list = field(default_factory=list)
    sitios: list = field(default_factory=list)
    energies: list = field(default_factory=list)     # eV, uno por sitio
    converged: list = field(default_factory=list)
    alturas: list = field(default_factory=list)      # Å tras relajar
    contactos: list = field(default_factory=list)    # Å, distancia mínima
    E_slab: float = None
    E_mol: float = None
    slab_ok: bool = None
    mol_ok: bool = None
    molecula: str = ""
    n_mol: int = 1
    natoms_slab: int = 0
    natoms_mol: int = 1
    altura_inicial: float = None
    relajado: bool = True

    @property
    def energias_ads(self) -> list:
        if self.E_slab is None or self.E_mol is None:
            return [None] * len(self.sitios)
        out = []
        for e in self.energies:
            out.append(None if e is None
                       else thermochem.adsorcion(e, self.E_slab, self.E_mol,
                                                 n=self.n_mol)["E_ads"])
        return out


# ----------------------------------------------------------------------
# Geometría de la superficie
# ----------------------------------------------------------------------
def atomos_superficie(slab, cara: str = "top", tol: float = TOL_CAPA) -> list:
    """Índices de los átomos de la capa expuesta."""
    z = slab.get_positions()[:, 2]
    if cara == "top":
        ref = z.max()
        return [i for i in range(len(slab)) if z[i] >= ref - tol]
    if cara == "bottom":
        ref = z.min()
        return [i for i in range(len(slab)) if z[i] <= ref + tol]
    raise ErrorDeUso(f"--face es 'top' o 'bottom'; recibí '{cara}'.")


def _replicas(slab, idx, n=1):
    """Posiciones de los átomos `idx` y sus réplicas periódicas en el plano."""
    cell = slab.cell.array
    pos = slab.get_positions()[idx]
    out, orig = [], []
    for i in range(-n, n + 1):
        for j in range(-n, n + 1):
            desp = i * cell[0] + j * cell[1]
            for k, p in enumerate(pos):
                out.append(p + desp)
                orig.append(idx[k])
    return np.array(out), np.array(orig)


# Cuántas distancias del entorno se comparan para decidir si dos sitios son
# el mismo, y con qué tolerancia en Å. Un número FIJO de vecinos y no un
# radio: con un radio, un átomo que cae justo en el borde entra en la firma
# de un sitio y no en la de su equivalente, y dos sitios idénticos salen
# distintos. Se comprobó en Al(111) 2x2, donde el radio daba 3 huecos
# distintos donde solo hay dos (fcc y hcp).
N_VECINOS_HUELLA = 24
TOL_HUELLA = 0.05


def _huella(slab, punto, k: int = N_VECINOS_HUELLA) -> np.ndarray:
    """Firma de un sitio: distancias a sus k vecinos más cercanos, ordenadas.

    Dos sitios equivalentes por simetría tienen el mismo entorno. Se incluyen
    TODAS las capas, no solo la superficial, porque es lo único que separa un
    hueco fcc de uno hcp: por arriba se ven idénticos y difieren en si hay
    un átomo debajo en la segunda capa.
    """
    cell = slab.cell.array
    anchos = [np.linalg.norm(cell[0]), np.linalg.norm(cell[1])]
    # réplicas suficientes para que la esfera de vecinos quepa entera
    n = max(2, int(np.ceil(12.0 / max(1e-6, min(anchos)))) + 1)
    todos, _ = _replicas(slab, list(range(len(slab))), n=n)
    d = np.sort(np.linalg.norm(todos - np.asarray(punto), axis=1))
    if len(d) < k:
        d = np.pad(d, (0, k - len(d)), constant_values=d[-1] if len(d) else 0.0)
    return d[:k]


def _misma_huella(a: np.ndarray, b: np.ndarray, tol: float = TOL_HUELLA) -> bool:
    return bool(np.max(np.abs(a - b)) < tol)


def sitios(slab, cara: str = "top", tipos=TIPOS, tol: float = TOL_CAPA,
           r_vecino: float = R_VECINO) -> list:
    """Enumera los sitios de adsorción no equivalentes de la cara expuesta."""
    from scipy.spatial import Delaunay

    idx = atomos_superficie(slab, cara, tol)
    if not idx:
        raise ErrorDeUso("no encontré átomos de superficie; ¿es una losa con "
                         "vacío? Córtala con 'olla-dft surface'.")
    pos = slab.get_positions()[idx]
    z_sup = pos[:, 2].max() if cara == "top" else pos[:, 2].min()
    rep, orig = _replicas(slab, idx, n=1)

    cand = []
    if "top" in tipos:
        for p in pos:
            cand.append(Sitio("top", (p[0], p[1]), p[2]))

    if "bridge" in tipos:
        for a in range(len(pos)):
            d = np.linalg.norm(rep - pos[a], axis=1)
            for b in np.where((d > 0.1) & (d <= r_vecino))[0]:
                m = (pos[a] + rep[b]) / 2.0
                cand.append(Sitio("bridge", (m[0], m[1]), float(m[2])))

    if "hollow" in tipos:
        plano = rep[:, :2]
        if len(plano) >= 4:
            try:
                tri = Delaunay(plano)
            except Exception:                               # noqa: BLE001
                tri = None
            if tri is not None:
                for simplex in tri.simplices:
                    v = rep[simplex]
                    lados = [np.linalg.norm(v[i] - v[j])
                             for i, j in ((0, 1), (1, 2), (0, 2))]
                    if max(lados) > r_vecino * 1.6:
                        continue    # triángulo estirado: no es un hueco real
                    c = v.mean(axis=0)
                    cand.append(Sitio("hollow", (c[0], c[1]), float(c[2])))

    # --- quedarse solo con los que están dentro de la celda y no se repiten ---
    cell2 = slab.cell.array[:2, :2]
    inv = np.linalg.inv(cell2.T)
    unicos = []
    for s in cand:
        f = inv @ np.array(s.xy)
        if not (-1e-6 <= f[0] < 1 - 1e-6 and -1e-6 <= f[1] < 1 - 1e-6):
            # traer a la celda de referencia antes de comparar
            f = f % 1.0
            xy = cell2.T @ f
            s = Sitio(s.tipo, (float(xy[0]), float(xy[1])), s.z)
        h = _huella(slab, (s.xy[0], s.xy[1], z_sup))
        if any(_misma_huella(h, u.huella) for u in unicos):
            continue
        s.huella = h
        unicos.append(s)

    orden = {"top": 0, "bridge": 1, "hollow": 2}
    unicos.sort(key=lambda s: (orden[s.tipo], s.xy))
    cuenta = {}
    for s in unicos:
        cuenta[s.tipo] = cuenta.get(s.tipo, 0) + 1
        s.etiqueta = f"{s.tipo}{cuenta[s.tipo]}"
    return unicos


def cargar_molecula(nombre: str):
    """Molécula por nombre (base de ASE) o desde un archivo."""
    p = Path(nombre)
    if p.exists():
        return structure.load(str(p))
    try:
        from ase.build import molecule
        return molecule(nombre)
    except Exception:                                       # noqa: BLE001
        try:
            from ase.collections import g2
            disponibles = ", ".join(sorted(g2.names)[:14])
        except Exception:                                   # noqa: BLE001
            disponibles = "CO, CO2, H2O, NH3, O2, CH4..."
        raise ErrorDeUso(
            f"no reconozco '{nombre}' ni como archivo ni como molécula de la "
            f"base de ASE. Algunas que sí: {disponibles}. También puedes "
            f"pasar un .xyz o .cif con la molécula.") from None


def colocar(slab, mol, sitio: Sitio, altura: float = 2.0,
            rotacion: float = 0.0, ancla: int = 0, cara: str = "top"):
    """Devuelve la losa con la molécula puesta en el sitio."""
    m = mol.copy()
    m.set_cell(slab.cell)
    m.set_pbc(slab.pbc)
    if rotacion:
        m.rotate(rotacion, "z", center=m.get_positions()[ancla])
    p = m.get_positions()
    ancla_pos = p[ancla]
    signo = 1.0 if cara == "top" else -1.0
    destino = np.array([sitio.xy[0], sitio.xy[1], sitio.z + signo * altura])
    m.set_positions(p + (destino - ancla_pos))
    out = slab.copy()
    out += m
    return out


# ----------------------------------------------------------------------
# Preparación
# ----------------------------------------------------------------------
def prepare(slab, molecula: str, outdir: str = "adsorb",
            altura: float = 2.0, tipos=TIPOS, cara: str = "top",
            rotaciones: int = 1, ancla: int = 0, pseudo_dir: str = None,
            insulator: bool = False, ecutwfc: float = None,
            ecutrho: float = None, kspacing: float = None,
            relax_ions: bool = True, vdw: str = None, dipolo: bool = False,
            nspin: int = 1, magnetization: dict = None) -> tuple:
    from qekit.core import kpoints as kp

    if 2 not in kp.direcciones_con_vacio(slab):
        raise ErrorDeUso(
            "esta estructura no tiene vacío en c: la energía de adsorción "
            "necesita una losa con vacío por encima. Córtala con "
            "'olla-dft surface -m \"1 1 1\" --vacuum 20'.")
    if rotaciones < 1:
        raise ErrorDeUso(f"--rotations debe ser al menos 1; recibí {rotaciones}.")
    tipos = tuple(tipos)
    malos = [t for t in tipos if t not in TIPOS]
    if malos:
        raise ErrorDeUso(
            f"tipo de sitio desconocido: {', '.join(malos)}. "
            f"Opciones: {', '.join(TIPOS)}.")

    mol = cargar_molecula(molecula)
    if ancla >= len(mol):
        raise ErrorDeUso(
            f"--anchor {ancla} no existe: la molécula tiene {len(mol)} átomos "
            f"(se numeran desde 0).")

    lista = sitios(slab, cara=cara, tipos=tipos)
    if rotaciones > 1 and len(mol) > 1:
        ampliada = []
        for s in lista:
            for k in range(rotaciones):
                s2 = Sitio(s.tipo, s.xy, s.z, s.huella,
                           rotacion=360.0 * k / rotaciones,
                           etiqueta=f"{s.etiqueta}_r{k}")
                ampliada.append(s2)
        lista = ampliada
    elif rotaciones > 1:
        # una molécula de un átomo no cambia al girar; girarla sería pagar
        # N veces el mismo cálculo para obtener N veces el mismo número.
        rotaciones = 1

    # Los pseudos y los cutoffs se resuelven sobre la UNIÓN de losa y
    # molécula: si se hicieran solo con la losa, el carbono y el oxígeno del
    # adsorbato no tendrían pseudo, y el cutoff sería el de los metales, que
    # es demasiado bajo para el oxígeno. Los tres cálculos comparten el
    # resultado, que es justo lo que hace la resta comparable.
    conjunto = slab + mol
    conjunto.set_cell(slab.cell)
    common = sweep.prepare_common(conjunto, pseudo_dir, ecutwfc, ecutrho,
                                  insulator,
                                  prefix=slab.get_chemical_formula(
                                      mode="hill", empirical=True))
    out = Path(outdir); out.mkdir(parents=True, exist_ok=True)
    grid = sweep.default_grid(slab, kspacing)
    calc = "relax" if relax_ions else "scf"

    run = AdsorbRun(sitios=lista, molecula=molecula, natoms_slab=len(slab),
                    natoms_mol=len(mol), altura_inicial=altura,
                    relajado=relax_ions)

    # Los tres tipos de cálculo comparten celda, cutoffs, malla k y pseudos:
    # es la única forma de que la resta signifique algo. La molécula va en la
    # MISMA celda que la losa, no en una caja aparte, por lo mismo.
    mol_sola = mol.copy()
    mol_sola.set_cell(slab.cell)
    mol_sola.set_pbc(slab.pbc)
    mol_sola.center()

    # --dipole: la sierra de la corrección dipolar (tefield/dipfield, edir=3)
    # va en los TRES cálculos, no solo en el de la losa con adsorbato: si la
    # referencia se calcula sin corregir, la resta arrastra el error.
    extras = dict(vdw=vdw, nspin=nspin, magnetization=magnetization,
                  dipole_correction=3 if dipolo else False)
    run.jobs.append(sweep.write_scf_job(
        slab, common, out / "_losa", "losa limpia", grid,
        meta={"papel": "slab"}, calculation=calc, **extras))
    run.jobs.append(sweep.write_scf_job(
        mol_sola, common, out / "_molecula", f"{molecula} aislada", grid,
        meta={"papel": "mol"}, calculation=calc, **extras))

    for s in lista:
        sistema = colocar(slab, mol, s, altura=altura, rotacion=s.rotacion,
                          ancla=ancla, cara=cara)
        job = sweep.write_scf_job(
            sistema, common, out / s.etiqueta, s.etiqueta, grid,
            meta={"papel": "ads", "sitio": s.etiqueta}, calculation=calc,
            **extras)
        run.jobs.append(job)

    sweep.write_run_script(run.jobs, out / "run.sh")

    cuenta = {}
    for s in lista:
        cuenta[s.tipo] = cuenta.get(s.tipo, 0) + 1
    report = ["--- Sitios de adsorción ---",
              f"Losa: {slab.get_chemical_formula()} ({len(slab)} átomos), "
              f"cara {cara}",
              f"Adsorbato: {molecula} ({len(mol)} átomos), "
              f"ancla = átomo {ancla} ({mol.get_chemical_symbols()[ancla]}), "
              f"altura inicial {altura:g} Å",
              "Sitios no equivalentes: "
              + ", ".join(f"{n} {t}" for t, n in sorted(cuenta.items()))
              + f"  ({len(lista)} cálculos"
              + (f", {rotaciones} rotaciones cada uno" if rotaciones > 1 else "")
              + ")",
              f"Malla k: {grid[0]}x{grid[1]}x{grid[2]}  |  "
              + ("posiciones relajadas" if relax_ions else "posiciones fijas"),
              "Referencias: losa limpia y molécula aislada, en la MISMA celda "
              "y con los\n  mismos cutoffs y malla, para que la resta sea "
              "válida."]
    if not vdw:
        report.append(
            "AVISO: sin corrección de van der Waals. En fisisorción (moléculas "
            "cerradas\n  sobre superficies) el enlace ES dispersión: sin --vdw "
            "la energía sale\n  cerca de cero y la geometría desligada.")
    if not dipolo and cara == "top":
        report.append(
            "Sugerencia: una molécula adsorbida en una sola cara deja la losa "
            "polar.\n  Con --dipole se cancela el dipolo artificial a través "
            "del vacío.")
    warn = sweep.missing_pseudo_warning(common)
    if warn:
        report.append(warn)
    report += ["", f"{len(run.jobs)} cálculos escritos en '{out.resolve()}'",
               "Córrelos con --run, o a mano con ./run.sh dentro de esa carpeta."]
    return run, "\n".join(report)


# ----------------------------------------------------------------------
# Recolección
# ----------------------------------------------------------------------
def _leer(job, por_dir):
    r = por_dir.get(str(job.directory))
    if r is not None and r.ok and r.result is not None:
        return r.result
    try:
        return qeout.read_xml(str(job.directory))
    except Exception:                                       # noqa: BLE001
        return None


def collect(run: AdsorbRun, results: list = None) -> AdsorbRun:
    por_dir = {str(r.job.directory): r for r in (results or [])}
    run.energies, run.converged, run.alturas, run.contactos = [], [], [], []
    n_slab = run.natoms_slab
    for job in run.jobs:
        papel = job.meta.get("papel")
        res = _leer(job, por_dir)
        if papel == "slab":
            run.E_slab = res.total_energy if res else None
            run.slab_ok = res.converged if res else None
            continue
        if papel == "mol":
            run.E_mol = res.total_energy if res else None
            run.mol_ok = res.converged if res else None
            continue
        run.energies.append(res.total_energy if res else None)
        run.converged.append(res.converged if res else None)
        if res is not None and res.positions is not None and len(res.positions) > n_slab:
            p = np.asarray(res.positions)
            losa, ads = p[:n_slab], p[n_slab:]
            run.alturas.append(float(ads[:, 2].min() - losa[:, 2].max()))
            d = np.linalg.norm(ads[:, None, :] - losa[None, :, :], axis=2)
            run.contactos.append(float(d.min()))
        else:
            run.alturas.append(None)
            run.contactos.append(None)
    return run


# ----------------------------------------------------------------------
# Reporte
# ----------------------------------------------------------------------
def report(run: AdsorbRun) -> str:
    if not run.energies:
        raise FaltanDatos(
            "no hay resultados todavía. Corre los cálculos (--run, o ./run.sh "
            "en la carpeta) y vuelve con --collect.")
    L = ["--- Energías de adsorción ---",
         f"Adsorbato: {run.molecula}   |   losa de {run.natoms_slab} átomos"]
    faltan = []
    if run.E_slab is None:
        faltan.append("la losa limpia")
    if run.E_mol is None:
        faltan.append("la molécula aislada")
    if faltan:
        L.append("")
        L.append("No se puede calcular E_ads: falta la energía de "
                 + " y ".join(faltan) + ".")
        L.append("  Sin las dos referencias la resta no existe; corre esos "
                 "dos cálculos\n  (están en _losa y _molecula) y vuelve con "
                 "--collect.")
        return "\n".join(L)

    L.append(f"E(losa) = {run.E_slab:.6f} eV    "
             f"E({run.molecula}) = {run.E_mol:.6f} eV")
    L.append("")
    L.append(f"  {'sitio':<12s} {'E_ads (eV)':>11s} {'altura (Å)':>11s} "
             f"{'contacto (Å)':>13s}")
    L.append("  " + "-" * 51)

    eads = run.energias_ads
    filas = sorted(
        [(i, s) for i, s in enumerate(run.sitios)],
        key=lambda t: (eads[t[0]] is None, eads[t[0]] if eads[t[0]] is not None else 0.0))
    for i, s in filas:
        e = eads[i]
        if e is None:
            L.append(f"  {s.etiqueta:<12s} {'sin resultado':>11s}")
            continue
        alt = run.alturas[i]
        con = run.contactos[i]
        fila = (f"  {s.etiqueta:<12s} {e:>11.4f} "
                f"{(alt if alt is not None else float('nan')):>11.3f} "
                f"{(con if con is not None else float('nan')):>13.3f}")
        if run.converged[i] is False:
            fila += "   << SIN CONVERGER"
        L.append(fila)

    validos = [(i, eads[i]) for i, _ in filas if eads[i] is not None]
    if validos:
        mejor_i, mejor_e = validos[0]
        s = run.sitios[mejor_i]
        L.append("")
        L.append(f"Sitio más favorable: {s.etiqueta} ({s.tipo}), "
                 f"E_ads = {mejor_e:.4f} eV")
        if mejor_e > 0:
            cerca = [run.contactos[i] for i, _ in filas
                     if run.contactos[i] is not None]
            if not run.relajado:
                L.append(
                    "  POSITIVA en todos los sitios, y el cálculo fue con "
                    "posiciones FIJAS:\n  lo más probable es que la altura "
                    f"inicial ({run.altura_inicial:g} Å) no sea la de "
                    "equilibrio\n  y estés midiendo la repulsión. Quita "
                    "--fixed-ions para que relaje.")
            elif cerca and min(cerca) < 1.2:
                L.append(
                    f"  POSITIVA, y el contacto más corto es {min(cerca):.2f} Å: "
                    "los átomos están\n  encima unos de otros. Revisa la "
                    "geometría antes de concluir nada.")
            else:
                L.append("  POSITIVA: en este nivel de teoría el adsorbato NO "
                         "se pega en ningún\n  sitio probado. Si esperabas "
                         "fisisorción, prueba con --vdw.")
        elif mejor_e > -0.30:
            L.append("  Es fisisorción débil (|E_ads| < 0.3 eV): a temperatura "
                     "ambiente la\n  molécula se desorbe. El número depende "
                     "mucho de la corrección de dispersión.")
        elif mejor_e < -2.0 and run.natoms_mol > 1:
            L.append("  Enlace muy fuerte (|E_ads| > 2 eV): normalmente hay "
                     "reacción, no\n  adsorción molecular. Mira la geometría "
                     "relajada: quizá la molécula\n  se disoció y estás "
                     "midiendo la energía de los fragmentos.")
        elif mejor_e < -2.0:
            L.append("  Es quimisorción fuerte. Ojo con la referencia: este "
                     "E_ads se mide contra\n  el ÁTOMO aislado, no contra la "
                     "molécula. Para compararlo con la\n  literatura de "
                     "moléculas diatómicas hay que restar media energía de\n"
                     "  disociación.")
        if len(validos) > 1:
            segundo = validos[1][1]
            L.append(f"  Diferencia con el segundo sitio: "
                     f"{abs(segundo - mejor_e):.4f} eV")
            if abs(segundo - mejor_e) < 0.05:
                L.append("  Los dos primeros están dentro de 50 meV: a esta "
                         "precisión no se puede\n  decir cuál gana; hacen falta "
                         "cutoffs y malla más finos para separarlos.")

    if run.slab_ok is False or run.mol_ok is False:
        L.append("")
        L.append("AVISO: alguna referencia no convergió; E_ads hereda ese error.")
    sin_conv = [run.sitios[i].etiqueta for i in range(len(run.sitios))
                if run.converged[i] is False]
    if sin_conv:
        L.append(f"SIN CONVERGER: {', '.join(sin_conv)}")
    return "\n".join(L)


def export(run: AdsorbRun, outdir: str = ".") -> list:
    out = Path(outdir); out.mkdir(parents=True, exist_ok=True)
    f = out / "ADSORCION.dat"
    eads = run.energias_ads
    lines = [provenance.header(
        f"energias de adsorcion de {run.molecula}",
        {"E_slab_eV": run.E_slab, "E_mol_eV": run.E_mol,
         "atomos_losa": run.natoms_slab}),
        f"# {'sitio':<12s} {'tipo':<8s} {'E_ads(eV)':>12s} "
        f"{'altura(A)':>11s} {'contacto(A)':>12s}"]
    nan = float("nan")
    for i, s in enumerate(run.sitios):
        if eads[i] is None:
            continue
        lines.append(
            f"  {s.etiqueta:<12s} {s.tipo:<8s} {eads[i]:>12.5f} "
            f"{(run.alturas[i] if run.alturas[i] is not None else nan):>11.3f} "
            f"{(run.contactos[i] if run.contactos[i] is not None else nan):>12.3f}")
    f.write_text("\n".join(lines) + "\n")
    txt = out / "ADSORCION.txt"
    txt.write_text(report(run) + "\n")
    return [str(f), str(txt)]


def plot(run: AdsorbRun, outfile: str = "adsorcion", formats="pdf,png",
         theme: str = None, size: str = None, family: str = None,
         background: str = None, palette=None, usetex: bool = None,
         width="single", journal: str = "generic", aspect: float = 0.70,
         mono: bool = False, dpi: int = None) -> list:
    """Barras de E_ads por sitio, ordenadas y coloreadas por tipo."""
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError as exc:                              # pragma: no cover
        raise RuntimeError("matplotlib no está instalado.") from exc

    eads = run.energias_ads
    datos = [(run.sitios[i].etiqueta, run.sitios[i].tipo, eads[i])
             for i in range(len(run.sitios)) if eads[i] is not None]
    if not datos:
        raise FaltanDatos("no hay energías de adsorción que graficar.")
    datos.sort(key=lambda t: t[2])

    st = qstyle.apply(theme, size=size, family=family, background=background,
                      palette=palette, usetex=usetex, mono=mono)
    fig, ax = qstyle.new_figure(width, journal, aspect)
    cols = qstyle.palette(3, mono=mono)
    color_de = {"top": cols[0], "bridge": cols[1], "hollow": cols[2]}
    x = np.arange(len(datos))
    ax.bar(x, [d[2] for d in datos], width=0.68,
           color=[color_de.get(d[1], cols[0]) for d in datos])
    ax.axhline(0.0, color=qstyle.INK_FAINT, lw=st["axis_line"])
    ax.set_xticks(x)
    ax.set_xticklabels([d[0] for d in datos], rotation=45, ha="right",
                       fontsize=st["legend"])
    ax.set_ylabel(r"$E_{\mathrm{ads}}$ (eV)")
    vistos = []
    for t in ("top", "bridge", "hollow"):
        if any(d[1] == t for d in datos):
            vistos.append(plt.Rectangle((0, 0), 1, 1, color=color_de[t], label=t))
    if len(vistos) > 1:
        ax.legend(handles=vistos, frameon=False, fontsize=st["legend"])
    written = qstyle.save(fig, outfile, formats, dpi=dpi, modulo="adsorcion")
    plt.close(fig)
    return written
