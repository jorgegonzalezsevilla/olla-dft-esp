# Olla-DFT — command-line toolkit for Quantum ESPRESSO
# Copyright (C) 2026 Jorge Enrique González Sevilla
# SPDX-License-Identifier: AGPL-3.0-or-later
# This program is free software: you can redistribute it and/or modify it
# under the terms of the GNU Affero General Public License as published by the Free
# Software Foundation, either version 3 of the License, or (at your option)
# any later version. See the LICENSE file for details.

"""Barrido de deformación: propiedades en función de la deformación aplicada.

Cuatro de los artículos que la auditoría de cobertura revisó hacen del barrido
de deformación su figura principal (gap contra deformación biaxial, banda de
conducción contra deformación, colapso del momento magnético). El cálculo en
sí es simple —deformar la celda, relajar y volver a mirar—, pero hacerlo a
mano son N carpetas, N inputs y una hoja de cálculo, y ahí es donde se cuelan
los errores: un punto sin converger que pasa desapercibido, una deformación
aplicada a la celda ya deformada del punto anterior, o comparar energías de
cálculos que no comparten cutoff.

Lo que sale de aquí:
  - E(ε), y el mínimo por parábola: si no cae en ε = 0, la estructura de
    partida no estaba relajada y todo lo demás hereda ese error.
  - gap(ε) y el POTENCIAL DE DEFORMACIÓN dEgap/dε, que es lo que se compara
    con el experimento.
  - momento magnético(ε), para las transiciones de espín.
  - para deformación biaxial en una lámina, el módulo biaxial 2D desde la
    curvatura de E(ε).
"""

from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

from qekit.core import provenance, qeout
from qekit.core.errors import ErrorDeUso, FaltanDatos
from qekit.core import style as qstyle
from qekit.modules import sweep


# Cada modo define qué componentes de la deformación se activan, en notación
# de Voigt (0..5 = xx, yy, zz, yz, xz, xy).
MODOS = {
    "biaxial":     ((0, 1),    "biaxial en el plano ab (εxx = εyy)"),
    "uniaxial-a":  ((0,),      "uniaxial a lo largo de a (εxx)"),
    "uniaxial-b":  ((1,),      "uniaxial a lo largo de b (εyy)"),
    "uniaxial-c":  ((2,),      "uniaxial a lo largo de c (εzz)"),
    "hidrostatica": ((0, 1, 2), "hidrostática (εxx = εyy = εzz)"),
    "cizalla":     ((5,),      "cizalla en el plano ab (εxy)"),
}

# Qué eje deja libre --relax-perp según el modo: el perpendicular al plano
# deformado. Para la hidrostática no hay perpendicular que valga.
_PERP = {"biaxial": "z", "uniaxial-a": "shape", "uniaxial-b": "shape",
         "uniaxial-c": "2Dxy", "cizalla": "shape"}

_VOIGT = [(0, 0), (1, 1), (2, 2), (1, 2), (0, 2), (0, 1)]


def matriz(modo: str, eps: float) -> np.ndarray:
    """Matriz de deformación 3x3 del modo dado, para una deformación `eps`."""
    if modo not in MODOS:
        raise ErrorDeUso(
            f"modo de deformación desconocido '{modo}'. "
            f"Opciones: {', '.join(sorted(MODOS))}")
    e = np.zeros((3, 3))
    for comp in MODOS[modo][0]:
        i, j = _VOIGT[comp]
        if i == j:
            e[i, j] = eps
        else:
            e[i, j] = e[j, i] = eps / 2.0
    return e


def rango(texto: str) -> list:
    """Convierte 'MIN:MAX:N' (en POR CIENTO) en una lista de fracciones.

    Se usa por ciento porque es como lo escriben los artículos ("±5 %"), y
    porque un rango en fracciones invita a teclear 5 queriendo decir 5 % y
    obtener una celda estirada al 500 %.
    """
    partes = str(texto).replace(",", ":").split(":")
    if len(partes) != 3:
        raise ErrorDeUso(
            "--range se escribe MIN:MAX:N en por ciento, por ejemplo "
            f"-5:5:11 (de -5 % a +5 % en 11 puntos); recibí '{texto}'.")
    try:
        lo, hi = float(partes[0]), float(partes[1])
        n = int(partes[2])
    except ValueError:
        raise ErrorDeUso(
            f"--range necesita dos números y un entero; recibí '{texto}'."
        ) from None
    if n < 3:
        raise ErrorDeUso(
            f"--range necesita al menos 3 puntos para poder ajustar una "
            f"curva; pediste {n}.")
    if hi <= lo:
        raise ErrorDeUso(
            f"el máximo de --range tiene que ser mayor que el mínimo; "
            f"recibí de {lo:g} % a {hi:g} %.")
    if max(abs(lo), abs(hi)) > 30.0:
        raise ErrorDeUso(
            f"±{max(abs(lo), abs(hi)):g} % es una deformación enorme: a esa "
            "escala la respuesta ya no es elástica y el cristal suele "
            "romperse o cambiar de fase. ¿Escribiste fracciones en vez de "
            "por ciento? --range va en POR CIENTO.")
    vals = [round(v / 100.0, 10) for v in np.linspace(lo, hi, n)]
    if not any(abs(v) < 1e-12 for v in vals):
        # el punto sin deformar es la referencia de todo lo demás
        vals = sorted(vals + [0.0])
    return vals


@dataclass
class StrainRun:
    modo: str = "biaxial"
    strains: list = field(default_factory=list)     # fracciones
    jobs: list = field(default_factory=list)
    energies: list = field(default_factory=list)    # eV/celda
    gaps: list = field(default_factory=list)        # eV o None
    homos: list = field(default_factory=list)
    lumos: list = field(default_factory=list)
    pressures: list = field(default_factory=list)   # GPa
    moments: list = field(default_factory=list)     # magnetones de Bohr
    converged: list = field(default_factory=list)
    natoms: int = 1
    area0: float = None        # Å², área ab de la celda sin deformar
    volume0: float = None      # Å³
    laminar: bool = False      # hay vacío en c
    relax_perp: bool = False

    @property
    def ok(self) -> list:
        """Índices de los puntos con energía leída."""
        return [i for i, e in enumerate(self.energies) if e is not None]


# ----------------------------------------------------------------------
# Preparación
# ----------------------------------------------------------------------
def prepare(atoms, modo: str = "biaxial", rangos: str = "-5:5:11",
            outdir: str = "strain", pseudo_dir: str = None,
            insulator: bool = False, ecutwfc: float = None,
            ecutrho: float = None, kspacing: float = None,
            relax_ions: bool = True, relax_perp: bool = False,
            nspin: int = 1, magnetization: dict = None,
            hubbard: dict = None, vdw: str = None) -> tuple:
    """Escribe un cálculo por cada deformación del rango."""
    from qekit.core import kpoints as kp

    if modo not in MODOS:
        raise ErrorDeUso(
            f"modo de deformación desconocido '{modo}'. "
            f"Opciones: {', '.join(sorted(MODOS))}")
    strains = rango(rangos)

    common = sweep.prepare_common(atoms, pseudo_dir, ecutwfc, ecutrho, insulator)
    out = Path(outdir); out.mkdir(parents=True, exist_ok=True)
    grid = sweep.default_grid(atoms, kspacing)

    # Bandas vacías: sin ellas no hay LUMO y la columna del gap sale vacía.
    # Con occupations='fixed' QE calcula SOLO las ocupadas y no avisa; el
    # usuario ve la tabla sin gap y no tiene forma de saber por qué.
    from qekit.modules.inputgen import _estimate_nbnd
    nbnd = _estimate_nbnd(atoms, common["pseudos"])
    if nbnd and nspin == 2:
        nbnd = int(nbnd * 1.2) + 2

    cell0 = np.asarray(atoms.cell.array, dtype=float)
    laminar = 2 in kp.direcciones_con_vacio(atoms)
    run = StrainRun(
        modo=modo, strains=strains, natoms=len(atoms),
        area0=float(np.linalg.norm(np.cross(cell0[0], cell0[1]))),
        volume0=float(abs(np.linalg.det(cell0))),
        laminar=laminar, relax_perp=relax_perp,
    )

    if relax_perp and modo == "hidrostatica":
        raise ErrorDeUso(
            "--relax-perp no tiene sentido con deformación hidrostática: se "
            "está deformando la celda en las tres direcciones, no queda eje "
            "perpendicular que relajar.")
    if modo == "biaxial" and not laminar:
        aviso_lam = ("AVISO: 'biaxial' deforma a y b y deja c fijo, que es lo "
                     "que se hace en una lámina.\n  Esta celda no tiene vacío "
                     "en c: si es material en bulto, quizá querías "
                     "'hidrostatica'.")
    else:
        aviso_lam = None

    if relax_ions and relax_perp:
        calc = "vc-relax"
        dofree = _PERP[modo]
    elif relax_ions:
        calc, dofree = "relax", None
    else:
        calc, dofree = "scf", None

    modo_txt = MODOS[modo][1]
    report = ["--- Barrido de deformación ---",
              f"Estructura: {atoms.get_chemical_formula()} ({len(atoms)} átomos)",
              f"Modo: {modo_txt}",
              f"Deformaciones: {len(strains)} puntos de "
              f"{min(strains) * 100:+.2f} % a {max(strains) * 100:+.2f} %",
              f"Malla k: {grid[0]}x{grid[1]}x{grid[2]}",
              f"Bandas: {nbnd} (con las vacías necesarias para leer el gap)"
              if nbnd else "Bandas: automáticas (sin UPF no se puede estimar; "
                           "puede que no haya gap en la tabla)",
              "Posiciones internas: " + ("relajadas" if relax_ions else "fijas")
              + (f"; celda con cell_dofree='{dofree}' (relajación de Poisson)"
                 if dofree else "; celda fija en la deformación impuesta")]
    if aviso_lam:
        report.append(aviso_lam)
    if nspin == 2 and not magnetization:
        report.append(
            "AVISO: --nspin 2 sin magnetización inicial suele converger a la\n"
            "  solución no magnética, y entonces el momento sale plano en todo\n"
            "  el barrido por un motivo numérico, no físico. Usa --mag.")
    warn = sweep.missing_pseudo_warning(common)
    if warn:
        report.append(warn)

    for eps in strains:
        e = matriz(modo, eps)
        deformed = atoms.copy()
        # SIEMPRE desde la celda original: aplicar la deformación sobre la
        # celda ya deformada del punto anterior acumula el error y el eje
        # deja de valer lo que dice la etiqueta.
        deformed.set_cell(cell0 @ (np.eye(3) + e), scale_atoms=True)
        etiqueta = f"{eps * 100:+.2f} %"
        nombre = f"e{eps * 100:+07.3f}".replace("+", "p").replace("-", "m").replace(".", "_")
        job = sweep.write_scf_job(
            deformed, common, out / nombre, etiqueta, grid,
            meta={"strain": eps}, calculation=calc, nbnd=nbnd,
            cell_dofree=dofree, nspin=nspin, magnetization=magnetization,
            hubbard=hubbard, vdw=vdw,
        )
        run.jobs.append(job)

    sweep.write_run_script(run.jobs, out / "run.sh")
    report += ["", f"{len(run.jobs)} cálculos escritos en '{out.resolve()}'",
               "Córrelos con --run, o a mano con ./run.sh dentro de esa carpeta."]
    return run, "\n".join(report)


# ----------------------------------------------------------------------
# Recolección
# ----------------------------------------------------------------------
def collect(run: StrainRun, results: list = None) -> StrainRun:
    for campo in ("energies", "gaps", "homos", "lumos", "pressures",
                  "moments", "converged"):
        setattr(run, campo, [])
    por_dir = {str(r.job.directory): r for r in (results or [])}
    for job in run.jobs:
        res = None
        r = por_dir.get(str(job.directory))
        if r is not None and r.ok:
            res = r.result
        if res is None:
            try:
                res = qeout.read_xml(str(job.directory))
            except Exception:                                # noqa: BLE001
                res = None
        if res is None:
            for campo in ("energies", "gaps", "homos", "lumos", "pressures",
                          "moments"):
                getattr(run, campo).append(None)
            run.converged.append(None)
            continue
        run.energies.append(res.total_energy)
        run.homos.append(res.homo)
        run.lumos.append(res.lumo)
        run.gaps.append(res.lumo - res.homo
                        if (res.homo is not None and res.lumo is not None)
                        else None)
        run.pressures.append(res.pressure)
        run.moments.append(res.total_magnetization)
        run.converged.append(res.converged)
    return run


# ----------------------------------------------------------------------
# Análisis
# ----------------------------------------------------------------------
def minimo(run: StrainRun) -> tuple:
    """Deformación de energía mínima por parábola. Devuelve (eps*, E*)."""
    idx = [i for i in run.ok if run.energies[i] is not None]
    if len(idx) < 3:
        return None, None
    x = np.array([run.strains[i] for i in idx])
    y = np.array([run.energies[i] for i in idx])
    # solo la zona cercana al mínimo: lejos, E(ε) deja de ser parabólica
    j = int(np.argmin(y))
    lo, hi = max(0, j - 3), min(len(x), j + 4)
    if hi - lo < 3:
        return None, None
    a, b, c = np.polyfit(x[lo:hi], y[lo:hi], 2)
    if a <= 0:
        return None, None
    eps = -b / (2 * a)
    return float(eps), float(a * eps ** 2 + b * eps + c)


def potencial_deformacion(run: StrainRun) -> tuple:
    """dEgap/dε por regresión lineal. Devuelve (pendiente_eV, r2)."""
    idx = [i for i in run.ok if run.gaps[i] is not None]
    if len(idx) < 3:
        return None, None
    x = np.array([run.strains[i] for i in idx])
    y = np.array([run.gaps[i] for i in idx])
    m, b = np.polyfit(x, y, 1)
    pred = m * x + b
    ss_res = float(np.sum((y - pred) ** 2))
    ss_tot = float(np.sum((y - y.mean()) ** 2))
    r2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else 0.0
    return float(m), float(r2)


def modulo_biaxial(run: StrainRun) -> float:
    """(1/A₀)·d²E/dε² para deformación biaxial en una lámina, en N/m.

    Es la combinación C11 + 2C12 + C22 (para una lámina isótropa, 2(C11+C12)),
    NO el módulo de Young. Se devuelve tal cual y el reporte lo dice, porque
    llamarlo Young sería cómodo y falso.
    """
    if run.modo != "biaxial" or not run.area0:
        return None
    idx = [i for i in run.ok if run.energies[i] is not None]
    if len(idx) < 3:
        return None
    x = np.array([run.strains[i] for i in idx])
    y = np.array([run.energies[i] for i in idx])
    cerca = np.abs(x) <= 0.03          # el ajuste armónico solo vale cerca de 0
    if cerca.sum() < 3:
        cerca = np.ones_like(x, dtype=bool)
    a = np.polyfit(x[cerca], y[cerca], 2)[0]
    # 2a = d²E/dε² en eV; a Å⁻² -> N/m con 1 eV/Å² = 16.0218 N/m
    return float(2.0 * a / run.area0 * 16.021766)


def cierre_de_gap(run: StrainRun) -> float:
    """Deformación a la que el gap se anula, si el barrido la cruza."""
    idx = [i for i in run.ok if run.gaps[i] is not None]
    if len(idx) < 2:
        return None
    for a, b in zip(idx, idx[1:]):
        ga, gb = run.gaps[a], run.gaps[b]
        if ga is None or gb is None:
            continue
        if (ga > 0.02) != (gb > 0.02):
            xa, xb = run.strains[a], run.strains[b]
            if abs(ga - gb) < 1e-12:
                continue
            return float(xa + (0.02 - ga) * (xb - xa) / (gb - ga))
    return None


# ----------------------------------------------------------------------
# Reporte
# ----------------------------------------------------------------------
def report(run: StrainRun) -> str:
    if not run.energies:
        raise FaltanDatos(
            "no hay resultados que leer todavía. Corre los cálculos "
            "(--run, o ./run.sh en la carpeta del barrido) y vuelve con "
            "--collect.")
    n = run.natoms
    L = ["--- Barrido de deformación: resultados ---",
         f"Modo: {MODOS[run.modo][1]}",
         f"Celda de referencia: {run.natoms} átomos, "
         f"V₀ = {run.volume0:.3f} Å³"
         + (f", A₀ = {run.area0:.3f} Å²" if run.laminar else ""),
         ""]

    e_ref = None
    for i, eps in enumerate(run.strains):
        if abs(eps) < 1e-12 and run.energies[i] is not None:
            e_ref = run.energies[i]
    if e_ref is None:
        vals = [e for e in run.energies if e is not None]
        e_ref = min(vals) if vals else 0.0

    hay_gap = any(g is not None for g in run.gaps)
    sin_gap = (not hay_gap) and any(h is not None for h in run.homos)
    hay_mag = any(m for m in run.moments)
    cab = f"  {'ε (%)':>8s} {'ΔE (meV/át)':>13s} {'P (GPa)':>10s}"
    if hay_gap:
        cab += f" {'gap (eV)':>10s}"
    if hay_mag:
        cab += f" {'M (μB)':>9s}"
    L.append(cab)
    L.append("  " + "-" * (len(cab) - 2))

    fallidos = []
    for i, eps in enumerate(run.strains):
        e = run.energies[i]
        if e is None:
            fallidos.append(f"{eps * 100:+.2f} %")
            continue
        fila = (f"  {eps * 100:>8.2f} {(e - e_ref) / n * 1000:>13.2f} "
                f"{(run.pressures[i] if run.pressures[i] is not None else float('nan')):>10.2f}")
        if hay_gap:
            g = run.gaps[i]
            fila += f" {g:>10.4f}" if g is not None else f" {'—':>10s}"
        if hay_mag:
            m = run.moments[i]
            fila += f" {m:>9.3f}" if m is not None else f" {'—':>9s}"
        if run.converged[i] is False:
            fila += "   << SIN CONVERGER"
        L.append(fila)

    L.append("")
    eps0, _ = minimo(run)
    if eps0 is not None:
        L.append(f"Mínimo de energía en ε = {eps0 * 100:+.3f} %")
        if abs(eps0) > 0.003:
            L.append(
                "  AVISO: el mínimo no cae en ε = 0. La estructura de partida\n"
                "  no estaba relajada, y todo el barrido está medido desde un\n"
                "  punto que no es el de equilibrio. Relájala primero:\n"
                "    olla-dft gen ESTRUCTURA -p vc-relax")

    if hay_gap:
        m, r2 = potencial_deformacion(run)
        if m is not None:
            L.append(f"Potencial de deformación dEgap/dε = {m:+.3f} eV "
                     f"(ajuste lineal, R² = {r2:.4f})")
            if r2 < 0.9:
                L.append("  El gap no responde de forma lineal en este rango "
                         "(R² bajo): la pendiente\n  es un promedio, no una "
                         "constante; mira la gráfica antes de citarla.")
        cierre = cierre_de_gap(run)
        if cierre is not None:
            L.append(f"El gap se cierra alrededor de ε = {cierre * 100:+.2f} % "
                     "(transición a metal)")

    if hay_mag:
        ms = [(run.strains[i], run.moments[i]) for i in run.ok
              if run.moments[i] is not None]
        if ms:
            lo = min(ms, key=lambda t: abs(t[1]))
            L.append(f"Momento magnético: de {max(abs(m) for _, m in ms):.3f} μB "
                     f"a un mínimo de {abs(lo[1]):.3f} μB en ε = {lo[0] * 100:+.2f} %")

    if sin_gap:
        L.append("No hay gap en la tabla: los cálculos no tienen bandas "
                 "vacías, así que\n  hay un HOMO pero no un LUMO. Vuelve a "
                 "preparar el barrido (sin --collect)\n  para que se incluyan, "
                 "o añade nbnd a mano en los inputs.")

    y2d = modulo_biaxial(run)
    if y2d is not None and run.laminar:
        L.append(f"Módulo biaxial 2D = (1/A₀)·d²E/dε² = {y2d:.1f} N/m")
        L.append("  Es la combinación C11 + 2C12 + C22, no el módulo de Young. "
                 "Para las Cij\n  por separado:  olla-dft elastic --2d")

    if fallidos:
        L.append("")
        L.append(f"Sin resultado en {len(fallidos)} punto(s): "
                 + ", ".join(fallidos))
    sin_conv = [f"{run.strains[i] * 100:+.2f} %" for i in range(len(run.strains))
                if run.converged[i] is False]
    if sin_conv:
        L.append(f"SIN CONVERGER en {len(sin_conv)} punto(s): "
                 + ", ".join(sin_conv)
                 + "\n  Esos puntos NO son comparables con el resto; "
                   "vuelve a correrlos antes de leer la curva.")
    return "\n".join(L)


# ----------------------------------------------------------------------
# Exportación y figura
# ----------------------------------------------------------------------
def export(run: StrainRun, outdir: str = ".") -> list:
    out = Path(outdir); out.mkdir(parents=True, exist_ok=True)
    f = out / "STRAIN.dat"
    lines = [provenance.header(
        f"barrido de deformación ({MODOS[run.modo][1]})",
        {"atomos": run.natoms, "V0_A3": f"{run.volume0:.4f}"}),
        f"# {'eps':>10s} {'E(eV)':>18s} {'gap(eV)':>12s} "
        f"{'P(GPa)':>12s} {'M(muB)':>10s}"]
    nan = float("nan")
    for i, eps in enumerate(run.strains):
        if run.energies[i] is None:
            continue
        g = run.gaps[i] if run.gaps[i] is not None else nan
        p = run.pressures[i] if run.pressures[i] is not None else nan
        m = run.moments[i] if run.moments[i] is not None else nan
        lines.append(f"{eps:12.6f} {run.energies[i]:18.8f} {g:12.5f} "
                     f"{p:12.4f} {m:10.4f}")
    f.write_text("\n".join(lines) + "\n")
    txt = out / "STRAIN.txt"
    txt.write_text(report(run) + "\n")
    return [str(f), str(txt)]


def plot(run: StrainRun, outfile: str = "strain", formats="pdf,png",
         theme: str = None, size: str = None, family: str = None,
         background: str = None, palette=None, usetex: bool = None,
         width="single", journal: str = "generic", aspect: float = 0.80,
         mono: bool = False, dpi: int = None) -> list:
    """E(ε) arriba y la propiedad electrónica abajo, compartiendo el eje x."""
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError as exc:                              # pragma: no cover
        raise RuntimeError("matplotlib no está instalado.") from exc

    idx = [i for i in run.ok if run.energies[i] is not None]
    if len(idx) < 2:
        raise FaltanDatos("hacen falta al menos dos puntos con energía "
                          "para graficar el barrido.")
    st = qstyle.apply(theme, size=size, family=family, background=background,
                      palette=palette, usetex=usetex, mono=mono)
    x = np.array([run.strains[i] * 100 for i in idx])
    e_ref = min(run.energies[i] for i in idx)
    y = np.array([(run.energies[i] - e_ref) / run.natoms * 1000 for i in idx])

    hay_gap = any(run.gaps[i] is not None for i in idx)
    hay_mag = any(run.moments[i] for i in idx)
    npan = 1 + int(hay_gap or hay_mag)
    fig, axes = qstyle.new_figure(width, journal, aspect * npan, nrows=npan,
                                  sharex=True)
    axes = np.atleast_1d(axes)
    cols = qstyle.palette(3, mono=mono)

    ax = axes[0]
    ax.plot(x, y, marker="o", ms=4, lw=st["line"], color=cols[0])
    eps0, _ = minimo(run)
    if eps0 is not None and min(x) <= eps0 * 100 <= max(x):
        ax.axvline(eps0 * 100, color=qstyle.INK_FAINT, lw=st["axis_line"],
                   dashes=[1.5, 1.5])
    ax.set_ylabel(r"$\Delta E$ (meV/átomo)")

    if npan == 2:
        ax2 = axes[1]
        if hay_gap:
            gi = [i for i in idx if run.gaps[i] is not None]
            ax2.plot([run.strains[i] * 100 for i in gi],
                     [run.gaps[i] for i in gi],
                     marker="s", ms=4, lw=st["line"], color=cols[1])
            ax2.set_ylabel("gap (eV)")
            cierre = cierre_de_gap(run)
            if cierre is not None:
                ax2.axhline(0.0, color=qstyle.INK_FAINT, lw=st["axis_line"],
                            dashes=[3.5, 2.0])
        elif hay_mag:
            mi = [i for i in idx if run.moments[i] is not None]
            ax2.plot([run.strains[i] * 100 for i in mi],
                     [run.moments[i] for i in mi],
                     marker="^", ms=4, lw=st["line"], color=cols[2])
            ax2.set_ylabel(r"$M$ ($\mu_B$/celda)")
        ax2.set_xlabel(r"deformación $\varepsilon$ (%)")
    else:
        ax.set_xlabel(r"deformación $\varepsilon$ (%)")

    written = qstyle.save(fig, outfile, formats, dpi=dpi, modulo="strain")
    plt.close(fig)
    return written
