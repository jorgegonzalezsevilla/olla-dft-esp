# Olla-DFT — command-line toolkit for Quantum ESPRESSO
# Copyright (C) 2026 Jorge Enrique González Sevilla
# SPDX-License-Identifier: AGPL-3.0-or-later
# This program is free software: you can redistribute it and/or modify it
# under the terms of the GNU Affero General Public License as published by the Free
# Software Foundation, either version 3 of the License, or (at your option)
# any later version. See the LICENSE file for details.

"""Pruebas de convergencia: ecutwfc, ecutrho y malla de puntos k.

Es lo primero que hay que hacer con un sistema nuevo y lo que más se
omite. Olla-DFT genera la serie, la corre si se le pide, y dice a partir de
qué valor la energía deja de cambiar más que el umbral.

Criterio: se compara cada punto contra el **más denso** de la serie (el
último), no contra el anterior. Comparar puntos consecutivos es el error
habitual: dos valores contiguos pueden parecerse por casualidad en mitad de
una curva que todavía no ha aplanado, y se concluye que ya convergió.
"""

from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

from qekit.core import kpoints, provenance, qeout
from qekit.core import style as qstyle
from qekit.modules import sweep
from qekit.core.errors import ErrorDeUso, FaltanDatos

KIND_LABEL = {
    "ecutwfc": "Cutoff de funciones de onda",
    "ecutrho": "Cutoff de densidad de carga",
    "kmesh": "Malla de puntos k",
}


@dataclass
class ConvergenceRun:
    kind: str
    values: list = field(default_factory=list)      # valor de cada punto
    labels: list = field(default_factory=list)      # etiqueta legible
    energies: list = field(default_factory=list)    # eV por celda (None si falló)
    natoms: int = 1
    threshold: float = 1.0                          # meV/átomo
    jobs: list = field(default_factory=list)

    def per_atom_diffs(self) -> np.ndarray:
        """|E(i) − E(referencia)| en meV/átomo, con el más denso de referencia."""
        ok = [e for e in self.energies if e is not None]
        if len(ok) < 2:
            return np.array([])
        ref = ok[-1]
        return np.array([
            abs(e - ref) * 1000.0 / self.natoms if e is not None else np.nan
            for e in self.energies
        ])

    def converged_index(self):
        """Primer punto a partir del cual todos quedan bajo el umbral."""
        d = self.per_atom_diffs()
        if d.size == 0:
            return None
        for i in range(len(d)):
            tail = d[i:]
            if np.all(np.isnan(tail) | (tail <= self.threshold)):
                return i
        return None


# ----------------------------------------------------------------------
# Preparación
# ----------------------------------------------------------------------
def prepare(atoms, kind: str, outdir: str = "convergencia",
            values: list = None, threshold: float = 1.0,
            pseudo_dir: str = None, insulator: bool = False,
            ecutwfc: float = None, ecutrho: float = None,
            kspacing: float = None, dual: float = None) -> tuple:
    """Escribe la serie de cálculos. Devuelve (ConvergenceRun, reporte)."""
    if kind not in KIND_LABEL:
        raise ErrorDeUso(
            f"tipo de convergencia desconocido '{kind}'. "
            f"Opciones: {', '.join(KIND_LABEL)}"
        )
    common = sweep.prepare_common(atoms, pseudo_dir, ecutwfc, ecutrho, insulator)
    out = Path(outdir)
    out.mkdir(parents=True, exist_ok=True)

    run = ConvergenceRun(kind=kind, natoms=len(atoms), threshold=threshold)
    report = [f"--- Convergencia: {KIND_LABEL[kind]} ---",
              f"Estructura: {atoms.get_chemical_formula()} ({len(atoms)} átomos)",
              f"Umbral: {threshold:g} meV/átomo"]
    warn = sweep.missing_pseudo_warning(common)
    if warn:
        report.append(warn)

    if kind == "ecutwfc":
        vals = values or list(range(30, 101, 10))
        d = dual if dual else common["ecutrho"] / common["ecutwfc"]
        grid = sweep.default_grid(atoms, kspacing)
        report.append(f"Malla k fija: {grid[0]}x{grid[1]}x{grid[2]}  |  "
                      f"ecutrho = {d:g} x ecutwfc")
        for v in vals:
            label = f"ecutwfc = {v:g} Ry"
            job = sweep.write_scf_job(
                atoms, common, out / f"ecutwfc_{v:g}", label, grid,
                ecutwfc=float(v), ecutrho=float(v) * d, meta={"value": float(v)},
            )
            run.values.append(float(v)); run.labels.append(label); run.jobs.append(job)

    elif kind == "ecutrho":
        base = common["ecutwfc"]
        vals = values or [base * f for f in (4, 6, 8, 10, 12)]
        grid = sweep.default_grid(atoms, kspacing)
        report.append(f"ecutwfc fijo en {base:g} Ry  |  "
                      f"malla k {grid[0]}x{grid[1]}x{grid[2]}")
        for v in vals:
            label = f"ecutrho = {v:g} Ry (dual {v / base:.1f})"
            job = sweep.write_scf_job(
                atoms, common, out / f"ecutrho_{v:g}", label, grid,
                ecutwfc=base, ecutrho=float(v), meta={"value": float(v)},
            )
            run.values.append(float(v)); run.labels.append(label); run.jobs.append(job)

    else:  # kmesh
        if values:
            grids = [tuple(int(x) for x in str(v).split("x")) if "x" in str(v)
                     else kpoints.kgrid_from_spacing(atoms, float(v))
                     for v in values]
            shown = values
        else:
            spacings = [0.40, 0.30, 0.25, 0.20, 0.15, 0.12]
            grids, shown = [], []
            for sp in spacings:
                g = kpoints.kgrid_from_spacing(atoms, sp)
                if g not in grids:
                    grids.append(g); shown.append(sp)
        report.append(f"Cutoffs fijos: ecutwfc = {common['ecutwfc']:g} Ry, "
                      f"ecutrho = {common['ecutrho']:g} Ry")
        for g, sp in zip(grids, shown):
            label = f"malla {g[0]}x{g[1]}x{g[2]}"
            job = sweep.write_scf_job(
                atoms, common, out / f"k_{g[0]}x{g[1]}x{g[2]}", label, g,
                meta={"value": float(np.prod(g)), "grid": g, "spacing": sp},
            )
            run.values.append(float(np.prod(g)))
            run.labels.append(label); run.jobs.append(job)

    sweep.write_run_script(run.jobs, out / "run.sh")
    report.append("")
    report.append(f"{len(run.jobs)} cálculos escritos en '{out.resolve()}'")
    report.append("Córrelos con --run, o a mano con ./run.sh dentro de esa carpeta.")
    return run, "\n".join(report)


# ----------------------------------------------------------------------
# Recolección y análisis
# ----------------------------------------------------------------------
def collect(run: ConvergenceRun, results: list = None) -> ConvergenceRun:
    """Rellena las energías desde los resultados (o leyendo las carpetas)."""
    run.energies = []
    if results is not None:
        by_dir = {str(r.job.directory): r for r in results}
        for job in run.jobs:
            r = by_dir.get(str(job.directory))
            run.energies.append(r.energy if (r and r.ok) else None)
        return run
    for job in run.jobs:
        try:
            run.energies.append(qeout.read_xml(str(job.directory)).total_energy)
        except Exception:
            run.energies.append(None)
    return run


def report(run: ConvergenceRun) -> str:
    lines = [f"--- Resultado de convergencia: {KIND_LABEL[run.kind]} ---"]
    d = run.per_atom_diffs()
    if d.size == 0:
        lines.append("No hay suficientes cálculos terminados para analizar.")
        return "\n".join(lines)

    lines.append(f"{'punto':>26s}  {'E (Ry/celda)':>16s}  "
                 f"{'ΔE vs. más denso':>18s}")
    for label, e, diff in zip(run.labels, run.energies, d):
        if e is None:
            lines.append(f"{label:>26s}  {'FALLÓ':>16s}")
            continue
        lines.append(f"{label:>26s}  {e / qeout.RY_EV:16.8f}  "
                     f"{diff:14.2f} meV/át")

    idx = run.converged_index()
    lines.append("")
    if idx is None:
        lines.append(
            f"NO converge dentro de {run.threshold:g} meV/átomo con los valores "
            "probados.\nExtiende la serie hacia valores más densos."
        )
    elif idx == len(run.labels) - 1:
        lines.append(
            f"Solo el último punto queda bajo {run.threshold:g} meV/átomo, así que "
            "no hay margen\npara asegurar que ahí ya aplanó: conviene extender la "
            "serie más allá."
        )
    else:
        lines.append(f"CONVERGE en: {run.labels[idx]}")
        lines.append(f"  A partir de ahí ningún punto se aparta más de "
                     f"{run.threshold:g} meV/átomo del más denso.")
        if run.kind == "ecutwfc":
            lines.append(f"  Úsalo con:  olla-dft gen estructura.cif "
                         f"--ecutwfc {run.values[idx]:g}")
        elif run.kind == "kmesh":
            g = run.jobs[idx].meta.get("grid")
            if g:
                lines.append(f"  Malla recomendada: {g[0]}x{g[1]}x{g[2]}")
    lines.append("")
    lines.append("Recuerda que la convergencia depende de la propiedad: la energía "
                 "total\nconverge antes que los esfuerzos o los fonones. Para "
                 "constantes elásticas\nconviene subir el cutoff por encima de lo "
                 "que pide la energía.")
    return "\n".join(lines)


def export(run: ConvergenceRun, outdir: str = ".") -> list:
    out = Path(outdir); out.mkdir(parents=True, exist_ok=True)
    d = run.per_atom_diffs()
    fname = out / "CONVERGENCIA.dat"
    lines = [provenance.header(
                 f"convergencia de {KIND_LABEL[run.kind]}",
                 {"umbral": f"{run.threshold:g} meV/atomo"}),
             f"# {'valor':>14s} {'E(Ry/celda)':>18s} {'dE(meV/atomo)':>16s}"]
    for v, e, diff in zip(run.values, run.energies, d if d.size else [np.nan]*len(run.values)):
        if e is None:
            continue
        lines.append(f"{v:16.6f} {e / qeout.RY_EV:18.10f} {diff:16.4f}")
    fname.write_text("\n".join(lines) + "\n")
    written = [str(fname)]
    txt = out / "CONVERGENCIA.txt"
    txt.write_text(report(run) + "\n")
    written.append(str(txt))
    return written


def plot(run: ConvergenceRun, outfile: str = "convergencia",
         formats="pdf,png", theme: str = None, size: str = None,
         family: str = None, background: str = None, palette=None,
         usetex: bool = None, width="single", journal: str = "generic",
         aspect: float = 0.75, mono: bool = False, dpi: int = None) -> list:
    """Curva de convergencia con la banda del umbral sombreada."""
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError as exc:  # pragma: no cover
        raise RuntimeError("matplotlib no está instalado.") from exc

    st = qstyle.apply(theme, size=size, family=family, background=background,
                      palette=palette, usetex=usetex, mono=mono)
    d = run.per_atom_diffs()
    xs = [v for v, e in zip(run.values, run.energies) if e is not None]
    ys = [diff for diff, e in zip(d, run.energies) if e is not None]
    if not xs:
        raise FaltanDatos("no hay puntos convergidos que graficar")

    fig, ax = qstyle.new_figure(width, journal, aspect)
    color = qstyle.palette(1, mono=mono)[0]
    ax.plot(xs, ys, marker="o", ms=4, color=color, lw=st["line"])
    ax.axhspan(0, run.threshold, color=color, alpha=0.10, lw=0)
    ax.axhline(run.threshold, color=qstyle.INK_FAINT, lw=st["axis_line"],
               dashes=[3.5, 2.0])
    ax.annotate(f"{run.threshold:g} meV/át", xy=(xs[-1], run.threshold),
                xytext=(-2, 3), textcoords="offset points", ha="right",
                fontsize=st["legend"], color=qstyle.INK_SOFT)

    idx = run.converged_index()
    if idx is not None and run.energies[idx] is not None:
        ax.axvline(run.values[idx], color=qstyle.INK_FAINT,
                   lw=st["axis_line"], dashes=[1.5, 1.5])

    xlabel = {"ecutwfc": "ecutwfc (Ry)", "ecutrho": "ecutrho (Ry)",
              "kmesh": "número de puntos k de la malla"}[run.kind]
    ax.set_xlabel(xlabel)
    ax.set_ylabel(r"$|\Delta E|$ (meV/átomo)")
    ax.set_yscale("log")
    fig.savefig  # noqa: B018  (el guardado real lo hace qstyle.save)
    written = qstyle.save(fig, outfile, formats, dpi=dpi,
                          modulo="convergencia")
    plt.close(fig)
    return written
