# Olla-DFT — command-line toolkit for Quantum ESPRESSO
# Copyright (C) 2026 Jorge Enrique González Sevilla
# SPDX-License-Identifier: GPL-3.0-or-later
# This program is free software: you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by the Free
# Software Foundation, either version 3 of the License, or (at your option)
# any later version. See the LICENSE file for details.

"""Invariantes topológicos de un Hamiltoniano de Wannier.

El módulo trabaja directamente con ``seedname_hr.dat`` (wannier90) o con
``WANNIER_hr.dat`` (Olla-DFT).  Calcula el número de Chern de un subespacio
ocupado aislado mediante la fórmula discreta, gauge-invariante, de
Fukui--Hatsugai--Suzuki y los centros híbridos asociados a lazos de Wilson.

Todas las coordenadas k son fraccionarias.  Por tanto, la fase de Fourier
es ``exp(2 pi i k.R)`` y no debe introducirse ningún factor de longitud.
"""

from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

from qekit.core import style as qstyle
from qekit.core.errors import ErrorDeUso, FaltanDatos
from qekit.modules import wannier


PLANES = {
    "xy": (0, 1, 2),
    "xz": (0, 2, 1),
    "yz": (1, 2, 0),
}


@dataclass
class TopologyRun:
    """Resultado de una sección bidimensional del modelo de Wannier."""

    model_path: str
    plane: str
    fixed: float
    grid: tuple
    occupied: int
    fermi: float = None
    direct_gap: float = float("nan")
    indirect_gap: float = float("nan")
    chern_raw: float = float("nan")
    chern: int = 0
    chern_residual: float = float("nan")
    min_overlap: float = float("nan")
    curvature: np.ndarray = None
    wilson: np.ndarray = None
    energies: np.ndarray = None
    warnings: list = field(default_factory=list)


def resolve_model(path) -> Path:
    """Resuelve un archivo ``*_hr.dat`` sin elegir uno ambiguo en silencio."""
    p = Path(path).expanduser()
    if p.is_file():
        return p.resolve()
    if not p.exists():
        raise FaltanDatos(f"no existe el modelo de Wannier '{p}'.")
    if not p.is_dir():
        raise ErrorDeUso(f"'{p}' no es un archivo ni una carpeta.")
    preferred = p / "WANNIER_hr.dat"
    if preferred.is_file():
        return preferred.resolve()
    candidates = sorted(p.glob("*_hr.dat"))
    if not candidates:
        raise FaltanDatos(
            f"en '{p}' no hay ningún WANNIER_hr.dat ni seedname_hr.dat.")
    if len(candidates) > 1:
        names = ", ".join(x.name for x in candidates)
        raise ErrorDeUso(
            f"en '{p}' hay varios modelos ({names}); indica el archivo exacto.")
    return candidates[0].resolve()


def _validate_grid(grid) -> tuple:
    try:
        n1, n2 = (int(x) for x in grid)
    except (TypeError, ValueError):
        raise ErrorDeUso("la malla topológica necesita dos enteros.") from None
    if n1 < 3 or n2 < 3:
        raise ErrorDeUso("la malla topológica debe ser de al menos 3x3.")
    return n1, n2


def kmesh(grid=(40, 40), plane="xy", fixed=0.0):
    """Malla periódica de una sección del BZ, sin repetir el extremo 1."""
    n1, n2 = _validate_grid(grid)
    if plane not in PLANES:
        raise ErrorDeUso(
            f"plano desconocido '{plane}'. Opciones: {', '.join(PLANES)}.")
    a, b, c = PLANES[plane]
    points = np.zeros((n1, n2, 3), float)
    points[..., c] = float(fixed) % 1.0
    points[..., a] = np.arange(n1)[:, None] / n1
    points[..., b] = np.arange(n2)[None, :] / n2
    return points


def _unitary_overlap(left, right):
    """Parte unitaria de ``left^dagger right`` y su menor valor singular."""
    overlap = left.conj().T @ right
    u, singular, vh = np.linalg.svd(overlap, full_matrices=False)
    return u @ vh, float(singular.min())


def invariants_from_vectors(vectors):
    """Chern FHS y lazos de Wilson de una malla de estados ocupados.

    ``vectors`` tiene forma ``(n1, n2, norb, nocc)``. La proyección polar de
    cada solape evita que errores de redondeo conviertan los enlaces en
    matrices no unitarias. El determinante hace que la curvatura dependa del
    subespacio ocupado completo y no de la gauge de cada banda.
    """
    vectors = np.asarray(vectors, complex)
    if vectors.ndim != 4 or vectors.shape[-1] < 1:
        raise ErrorDeUso(
            "los vectores ocupados deben tener forma (n1,n2,norb,nocc).")
    n1, n2, _norb, nocc = vectors.shape
    ux = np.empty((n1, n2), complex)
    uy = np.empty((n1, n2), complex)
    qx = np.empty((n1, n2, nocc, nocc), complex)
    minimum = np.inf
    for i in range(n1):
        for j in range(n2):
            qxi, sx = _unitary_overlap(vectors[i, j],
                                       vectors[(i + 1) % n1, j])
            qyi, sy = _unitary_overlap(vectors[i, j],
                                       vectors[i, (j + 1) % n2])
            qx[i, j] = qxi
            dx, dy = np.linalg.det(qxi), np.linalg.det(qyi)
            ux[i, j] = dx / abs(dx)
            uy[i, j] = dy / abs(dy)
            minimum = min(minimum, sx, sy)

    curvature = np.empty((n1, n2), float)
    for i in range(n1):
        for j in range(n2):
            loop = (ux[i, j] * uy[(i + 1) % n1, j]
                    * np.conj(ux[i, (j + 1) % n2]) * np.conj(uy[i, j]))
            curvature[i, j] = np.angle(loop)
    chern = float(curvature.sum() / (2.0 * np.pi))

    wilson = np.empty((n2, nocc), float)
    for j in range(n2):
        loop = np.eye(nocc, dtype=complex)
        for i in range(n1):
            loop = loop @ qx[i, j]
        wilson[j] = np.sort(np.mod(np.angle(np.linalg.eigvals(loop))
                                    / (2.0 * np.pi), 1.0))
    return curvature, chern, wilson, float(minimum)


def analyze(path, occupied=None, fermi=None, grid=(40, 40), plane="xy",
            fixed=0.0, gap_tol=1e-8) -> TopologyRun:
    """Analiza una sección 2D de un Hamiltoniano ``*_hr.dat``.

    La ocupación nunca se adivina: debe darse el número de bandas ocupadas o
    un nivel de Fermi que deje el mismo número de estados ocupados en toda la
    malla. Esto evita asignar un invariante de aislante a un metal.
    """
    if (occupied is None) == (fermi is None):
        raise ErrorDeUso(
            "indica exactamente una ocupación: --occupied N o --fermi EV.")
    if not np.isfinite(float(fixed)):
        raise ErrorDeUso("--fixed tiene que ser un número finito.")
    if not np.isfinite(float(gap_tol)) or float(gap_tol) <= 0:
        raise ErrorDeUso("--gap-tol tiene que ser un número positivo y finito.")
    if fermi is not None and not np.isfinite(float(fermi)):
        raise ErrorDeUso("--fermi tiene que ser un número finito.")
    n1, n2 = _validate_grid(grid)
    model = resolve_model(path)
    HR, R, deg = wannier.leer_hr(model)
    nw = HR.shape[1]
    points = kmesh((n1, n2), plane=plane, fixed=fixed)
    energy, eigvec = wannier.interpolar(
        HR, R, deg, points.reshape(-1, 3), vectores=True)
    energy = energy.reshape(n1, n2, nw)
    eigvec = eigvec.reshape(n1, n2, nw, nw)

    if occupied is not None:
        nocc = int(occupied)
        if float(occupied) != nocc:
            raise ErrorDeUso("--occupied tiene que ser un entero.")
    else:
        counts = np.sum(energy < float(fermi), axis=2)
        if counts.min() != counts.max():
            raise ErrorDeUso(
                f"el nivel de Fermi corta bandas: hay entre {counts.min()} y "
                f"{counts.max()} estados ocupados según k. El sistema es "
                "metálico en esta sección y el Chern de 'las ocupadas' no "
                "está definido.")
        nocc = int(counts.flat[0])
    if not 1 <= nocc < nw:
        raise ErrorDeUso(
            f"la ocupación debe estar entre 1 y {nw - 1}; recibí {nocc} "
            f"para un modelo de {nw} orbitales.")

    direct = float(np.min(energy[..., nocc] - energy[..., nocc - 1]))
    indirect = float(np.min(energy[..., nocc])
                     - np.max(energy[..., nocc - 1]))
    if direct <= float(gap_tol):
        raise ErrorDeUso(
            f"el subespacio ocupado no está aislado: gap directo mínimo "
            f"{direct:.3e} eV (tolerancia {gap_tol:.1e}). El número de Chern "
            "no está definido; aumenta la malla solo si esperabas un gap.")

    occupied_vectors = eigvec[..., :nocc]
    curvature, raw, wilson, minimum = invariants_from_vectors(occupied_vectors)
    integer = int(np.rint(raw))
    run = TopologyRun(
        model_path=str(model), plane=plane, fixed=float(fixed) % 1.0,
        grid=(n1, n2), occupied=nocc,
        fermi=None if fermi is None else float(fermi), direct_gap=direct,
        indirect_gap=indirect, chern_raw=raw, chern=integer,
        chern_residual=abs(raw - integer), min_overlap=minimum,
        curvature=curvature, wilson=wilson, energies=energy,
    )
    if indirect <= 0:
        run.warnings.append(
            "El gap indirecto no es positivo: el subespacio está aislado "
            "banda a banda, pero el sistema no es un aislante global.")
    if run.chern_residual > 1e-6:
        run.warnings.append(
            "El Chern discreto no cerró a un entero con precisión numérica; "
            "refina la malla y revisa la localización del modelo Wannier.")
    if minimum < 1e-6:
        run.warnings.append(
            "Hay subespacios ocupados casi ortogonales entre puntos vecinos; "
            "la malla puede ser demasiado gruesa.")
    return run


def report(run: TopologyRun) -> str:
    """Informe legible con los límites físicos explícitos."""
    occupation = (f"{run.occupied} bandas"
                  if run.fermi is None
                  else f"E < {run.fermi:g} eV ({run.occupied} bandas)")
    lines = [
        "--- Topología del modelo de Wannier ---",
        f"Modelo:       {run.model_path}",
        f"Sección BZ:   plano {run.plane}, coordenada fija {run.fixed:g}",
        f"Malla:        {run.grid[0]}x{run.grid[1]} (periódica)",
        f"Ocupación:    {occupation}",
        "",
        f"Gap directo mínimo: {run.direct_gap:.8g} eV",
        f"Gap indirecto:      {run.indirect_gap:.8g} eV",
        f"Chern discreto:     {run.chern_raw:+.12f}",
        f"Chern entero:       {run.chern:+d}",
        f"Residuo al entero:  {run.chern_residual:.3e}",
        f"Solape mínimo:      {run.min_overlap:.3e}",
        "",
        "Convención: la señal cambia al invertir la orientación del plano.",
        "Los centros de Wilson se exportan módulo 1; no se asigna un Z2 "
        "automático sin comprobar simetría de reversión temporal.",
    ]
    for warning in run.warnings:
        lines += ["", f"AVISO: {warning}"]
    return "\n".join(lines)


def export(run: TopologyRun, outdir="topology") -> list:
    """Exporta informe, curvatura de Berry y centros híbridos."""
    out = Path(outdir)
    out.mkdir(parents=True, exist_ok=True)
    written = []

    f = out / "TOPOLOGY_curvature.dat"
    n1, n2 = run.grid
    rows = [(i / n1, j / n2, run.curvature[i, j])
            for i in range(n1) for j in range(n2)]
    np.savetxt(f, rows, fmt="%14.9g",
               header="k1(frac)  k2(frac)  Berry_flux_plaquette(rad)")
    written.append(str(f))

    f = out / "TOPOLOGY_wilson.dat"
    transverse = np.arange(n2) / n2
    header = "k_transverse(frac)" + "".join(
        f"  center_{i + 1}(mod1)" for i in range(run.occupied))
    np.savetxt(f, np.column_stack([transverse, run.wilson]),
               fmt="%14.9g", header=header)
    written.append(str(f))

    f = out / "TOPOLOGY.txt"
    f.write_text(report(run) + "\n", encoding="utf-8")
    written.append(str(f))
    return written


def plot(run: TopologyRun, outfile="topology", formats="pdf,png",
         theme=None, size=None, family=None, background=None, palette=None,
         usetex=None, width="double", journal="generic", aspect=0.44,
         mono=False, dpi=None) -> list:
    """Figura final: flujo de Berry por plaqueta y centros de Wilson."""
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError as exc:                              # pragma: no cover
        raise RuntimeError("matplotlib no está instalado.") from exc

    qstyle.apply(theme, size=size, family=family, background=background,
                 palette=palette, usetex=usetex, mono=mono)
    fig, axes = plt.subplots(1, 2, figsize=qstyle.figure_size(
        width, journal=journal, aspect=aspect))
    extent = (0.0, 1.0, 0.0, 1.0)
    image = axes[0].imshow(run.curvature.T, origin="lower", extent=extent,
                           aspect="auto", cmap="RdBu_r")
    axes[0].set_xlabel("$k_1$ (fracc.)")
    axes[0].set_ylabel("$k_2$ (fracc.)")
    axes[0].set_title(f"Flujo de Berry  C = {run.chern:+d}")
    fig.colorbar(image, ax=axes[0], label="fase por plaqueta (rad)")

    transverse = np.arange(run.grid[1]) / run.grid[1]
    color = qstyle.palette(1, mono=mono)[0]
    for band in range(run.occupied):
        axes[1].scatter(transverse, run.wilson[:, band], s=8,
                        color=color, alpha=0.8)
    axes[1].set(xlabel="$k_2$ (fracc.)", ylabel="centro híbrido (mód. 1)",
                xlim=(0, 1), ylim=(0, 1), title="Lazos de Wilson")
    for label, ax in zip(("(a)", "(b)"), axes):
        qstyle.panel_label(ax, label)
        qstyle.finish_axes(ax)
    fig.tight_layout()
    written = qstyle.save(fig, outfile, formats, dpi=dpi, modulo="topology")
    plt.close(fig)
    return written
