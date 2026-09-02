# Olla-DFT — command-line toolkit for Quantum ESPRESSO
# Copyright (C) 2026 Jorge Enrique González Sevilla
# SPDX-License-Identifier: GPL-3.0-or-later
# This program is free software: you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by the Free
# Software Foundation, either version 3 of the License, or (at your option)
# any later version. See the LICENSE file for details.

"""Ecuación de estado E–V: volumen de equilibrio y módulo de bulk.

Se escala la celda en torno al volumen de partida, se calcula la energía en
cada volumen y se ajusta una ecuación de estado. De ahí salen el volumen de
equilibrio V₀, la energía mínima E₀, el módulo de bulk B₀ y su derivada B₀'.

Se ofrecen tres ecuaciones porque no dan lo mismo: Birch–Murnaghan de tercer
orden es la referencia habitual en estado sólido; Vinet describe mejor
compresiones grandes; Murnaghan es la más simple y sirve de contraste. Si
las tres coinciden, el ajuste es sólido; si no, suele faltar rango o sobran
puntos ruidosos.

Unidades: energías en eV, volúmenes en Å³, módulos en GPa
(1 eV/Å³ = 160.21766208 GPa).
"""

from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

from qekit.core import qeout
from qekit.core import provenance, structure
from qekit.core import style as qstyle
from qekit.modules import sweep
from qekit.core.errors import ErrorDeUso, FaltanDatos

EV_A3_GPA = 160.21766208


# ----------------------------------------------------------------------
# Ecuaciones de estado
# ----------------------------------------------------------------------
def birch_murnaghan(V, E0, V0, B0, Bp):
    """Birch–Murnaghan de tercer orden (B0 en eV/Å³)."""
    eta = (V0 / V) ** (2.0 / 3.0)
    return E0 + 9.0 * V0 * B0 / 16.0 * (
        (eta - 1.0) ** 3 * Bp + (eta - 1.0) ** 2 * (6.0 - 4.0 * eta)
    )


def murnaghan(V, E0, V0, B0, Bp):
    return (E0 + B0 * V / Bp * ((V0 / V) ** Bp / (Bp - 1.0) + 1.0)
            - B0 * V0 / (Bp - 1.0))


def vinet(V, E0, V0, B0, Bp):
    x = (V / V0) ** (1.0 / 3.0)
    xi = 1.5 * (Bp - 1.0)
    return (E0 + 9.0 * B0 * V0 / (xi ** 2)
            * (1.0 + (xi * (1.0 - x) - 1.0) * np.exp(xi * (1.0 - x))))


EQUATIONS = {
    "birch-murnaghan": (birch_murnaghan, "Birch–Murnaghan (3.er orden)"),
    "murnaghan": (murnaghan, "Murnaghan"),
    "vinet": (vinet, "Vinet"),
}
DEFAULT_EQ = "birch-murnaghan"


@dataclass
class EOSFit:
    equation: str = DEFAULT_EQ
    E0: float = None          # eV
    V0: float = None          # Å³
    B0: float = None          # GPa
    Bp: float = None
    a0: float = None          # Å, si la celda es cúbica
    rmse: float = None        # eV/átomo
    ok: bool = False
    message: str = ""


@dataclass
class EOSRun:
    scales: list = field(default_factory=list)     # factores lineales aplicados
    volumes: list = field(default_factory=list)    # Å³
    energies: list = field(default_factory=list)   # eV (None si falló)
    natoms: int = 1
    cubic: bool = False
    conv_ratio: float = 1.0     # V(celda convencional) / V(celda calculada)
    jobs: list = field(default_factory=list)
    fits: dict = field(default_factory=dict)

    def valid(self):
        v = [(V, E) for V, E in zip(self.volumes, self.energies) if E is not None]
        if not v:
            return np.array([]), np.array([])
        V, E = zip(*v)
        return np.array(V), np.array(E)


# ----------------------------------------------------------------------
# Preparación
# ----------------------------------------------------------------------
def prepare(atoms, outdir: str = "eos", npoints: int = 9, span: float = 0.10,
            pseudo_dir: str = None, insulator: bool = False,
            ecutwfc: float = None, ecutrho: float = None,
            kspacing: float = None, relax_ions: bool = False,
            center: float = 1.0) -> tuple:
    """Genera los cálculos a volumen fijo alrededor de la celda de entrada.

    `span` es la variación relativa de VOLUMEN a cada lado (0.10 = ±10 %).
    Con `relax_ions` se usa 'relax' en vez de 'scf', necesario cuando las
    posiciones internas no están fijadas por simetría.
    """
    if npoints < 5:
        raise ErrorDeUso("hacen falta al menos 5 puntos para un ajuste fiable")
    common = sweep.prepare_common(atoms, pseudo_dir, ecutwfc, ecutrho, insulator)
    out = Path(outdir); out.mkdir(parents=True, exist_ok=True)

    grid = sweep.default_grid(atoms, kspacing)
    V_ref = float(abs(np.linalg.det(atoms.cell.array)))
    # factores de VOLUMEN equiespaciados -> factor LINEAL = raíz cúbica
    # `center` es un factor LINEAL sobre la celda (lo que devuelve un
    # barrido previo con MLIP): centrar ahí el rango permite usar un span
    # estrecho y menos puntos, que es donde está el ahorro.
    centro_vol = float(center) ** 3
    vol_factors = np.linspace(centro_vol * (1.0 - span),
                              centro_vol * (1.0 + span), npoints)

    run = EOSRun(natoms=len(atoms))
    # Ojo: la celda PRIMITIVA de un cristal cúbico no tiene ángulos de 90°
    # (la FCC los tiene de 60°), así que mirar los parámetros de celda daría
    # "no cúbico" para el silicio. Se pregunta a la simetría, y el parámetro
    # de red se recupera del volumen de la celda CONVENCIONAL.
    try:
        ds = structure.symmetry_dataset(atoms)
        run.cubic = ds.number >= 195
        if run.cubic:
            conv = structure.conventional(atoms)
            V_conv = float(abs(np.linalg.det(conv.cell.array)))
            run.conv_ratio = V_conv / V_ref
    except Exception:
        run.cubic = False

    report = ["--- Ecuación de estado E–V ---",
              f"Estructura: {atoms.get_chemical_formula()} ({len(atoms)} átomos)",
              f"Volumen de partida: {V_ref:.4f} Å³",
              f"{npoints} puntos entre {centro_vol*(1-span)*100:.0f} % y {centro_vol*(1+span)*100:.0f} % "
              "del volumen",
              f"Malla k fija: {grid[0]}x{grid[1]}x{grid[2]}  "
              f"(la misma en todos los puntos, para que las energías sean "
              "comparables)"]
    warn = sweep.missing_pseudo_warning(common)
    if warn:
        report.append(warn)
    if relax_ions:
        report.append("Las posiciones internas se relajan en cada volumen "
                      "(calculation='relax').")

    for f in vol_factors:
        lin = f ** (1.0 / 3.0)
        scaled = atoms.copy()
        scaled.set_cell(atoms.cell.array * lin, scale_atoms=True)
        V = float(abs(np.linalg.det(scaled.cell.array)))
        label = f"V = {V:.3f} Å³ ({f * 100:.1f} %)"
        job = sweep.write_scf_job(
            scaled, common, out / f"V_{f:.4f}", label, grid,
            meta={"volume": V, "factor": f},
            calculation="relax" if relax_ions else "scf",
        )
        run.scales.append(lin); run.volumes.append(V); run.jobs.append(job)

    sweep.write_run_script(run.jobs, out / "run.sh")
    report += ["", f"{len(run.jobs)} cálculos escritos en '{out.resolve()}'",
               "Córrelos con --run, o a mano con ./run.sh dentro de esa carpeta."]
    return run, "\n".join(report)


# ----------------------------------------------------------------------
# Recolección y ajuste
# ----------------------------------------------------------------------
def collect(run: EOSRun, results: list = None) -> EOSRun:
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


def fit(run: EOSRun, equation: str = DEFAULT_EQ) -> EOSFit:
    """Ajusta una ecuación de estado a los puntos válidos."""
    from scipy.optimize import curve_fit

    if equation not in EQUATIONS:
        raise ErrorDeUso(f"ecuación desconocida '{equation}'. "
                         f"Opciones: {', '.join(EQUATIONS)}")
    func, _name = EQUATIONS[equation]
    V, E = run.valid()
    res = EOSFit(equation=equation)
    if len(V) < 4:
        res.message = (f"solo {len(V)} puntos válidos; hacen falta al menos 4 "
                       "para ajustar cuatro parámetros")
        return res

    # semilla: parábola sobre E(V), que ya da V0, E0 y B0 aproximados
    coef = np.polyfit(V, E, 2)
    V0_guess = -coef[1] / (2.0 * coef[0])
    if not (V.min() * 0.5 < V0_guess < V.max() * 1.5):
        V0_guess = float(V[np.argmin(E)])
    E0_guess = np.polyval(coef, V0_guess)
    B0_guess = max(2.0 * coef[0] * V0_guess, 1e-3)   # eV/Å³

    try:
        popt, _ = curve_fit(func, V, E,
                            p0=[E0_guess, V0_guess, B0_guess, 4.0],
                            maxfev=20000)
    except Exception as exc:
        res.message = f"el ajuste no convergió: {exc}"
        return res

    E0, V0, B0, Bp = popt
    if not (V.min() * 0.6 < V0 < V.max() * 1.4) or B0 <= 0:
        res.message = ("el ajuste dio parámetros fuera de rango; "
                       "amplía o centra mejor el barrido de volúmenes")
        return res

    residuals = E - func(V, *popt)
    res.E0, res.V0 = float(E0), float(V0)
    res.B0 = float(B0 * EV_A3_GPA)
    res.Bp = float(Bp)
    res.rmse = float(np.sqrt(np.mean(residuals ** 2)) / run.natoms)
    if run.cubic:
        # celda cúbica: el parámetro de red sale del volumen de la celda
        # CONVENCIONAL, V0 * conv_ratio (prepare() midió ese cociente con
        # spglib). Se fija aquí para que informe y exportación usen el
        # mismo número.
        res.a0 = float((res.V0 * run.conv_ratio) ** (1.0 / 3.0))
    res.ok = True
    return res


def fit_all(run: EOSRun) -> dict:
    run.fits = {eq: fit(run, eq) for eq in EQUATIONS}
    return run.fits


def report(run: EOSRun, cell_a: float = None) -> str:
    V, E = run.valid()
    lines = ["--- Ecuación de estado ---",
             f"Puntos válidos: {len(V)} de {len(run.volumes)}"]
    if len(V) == 0:
        lines.append("Ningún cálculo terminó; no hay nada que ajustar.")
        return "\n".join(lines)

    lines.append("")
    lines.append(f"{'V (Å³)':>12s} {'E (Ry/celda)':>18s} {'E (eV/átomo)':>16s}")
    for v, e in zip(run.volumes, run.energies):
        if e is None:
            lines.append(f"{v:12.4f} {'FALLÓ':>18s}")
        else:
            lines.append(f"{v:12.4f} {e / qeout.RY_EV:18.10f} "
                         f"{e / run.natoms:16.6f}")

    fits = run.fits or fit_all(run)
    lines += ["", f"{'ecuación':>26s} {'V0 (Å³)':>10s} {'B0 (GPa)':>10s} "
                  f"{'B0′':>7s} {'RMSE (meV/át)':>14s}"]
    for key, f in fits.items():
        name = EQUATIONS[key][1]
        if not f.ok:
            lines.append(f"{name:>26s}  {f.message}")
            continue
        lines.append(f"{name:>26s} {f.V0:10.4f} {f.B0:10.2f} {f.Bp:7.2f} "
                     f"{f.rmse * 1000:14.3f}")

    best = fits.get(DEFAULT_EQ)
    if best and best.ok:
        lines += ["", "Resultado (Birch–Murnaghan):",
                  f"  V₀ = {best.V0:.4f} Å³   ({best.V0 / run.natoms:.4f} Å³/átomo)",
                  f"  E₀ = {best.E0 / qeout.RY_EV:.8f} Ry",
                  f"  B₀ = {best.B0:.2f} GPa",
                  f"  B₀′ = {best.Bp:.2f}"]
        if run.cubic and best.a0 is not None:
            lines.append(f"  a₀ = {best.a0:.5f} Å   (parámetro de red cúbico)")

        ok_vals = [f.B0 for f in fits.values() if f.ok]
        if len(ok_vals) > 1:
            spread = max(ok_vals) - min(ok_vals)
            lines += ["", f"Las tres ecuaciones difieren en {spread:.1f} GPa "
                          f"({spread / best.B0 * 100:.1f} %)."]
            if spread / best.B0 > 0.05:
                lines.append("  Es bastante: suele indicar que faltan puntos o "
                             "que el rango de\n  volúmenes es muy estrecho.")
            else:
                lines.append("  Buena señal: el ajuste no depende de la "
                             "ecuación elegida.")
        V_in = (V.min() <= best.V0 <= V.max())
        if not V_in:
            lines.append("\nAVISO: V₀ cae FUERA del rango calculado. Vuelve a "
                         "correr el barrido\ncentrado en ese volumen.")
    return "\n".join(lines)


def export(run: EOSRun, outdir: str = ".", cell_a: float = None) -> list:
    out = Path(outdir); out.mkdir(parents=True, exist_ok=True)
    written = []
    fname = out / "EOS.dat"
    lines = [provenance.header("ecuación de estado",
                               {"puntos": len(run.volumes)}),
             f"# {'V(A^3)':>14s} {'E(Ry)':>20s} {'E(eV)':>18s}"]
    for v, e in zip(run.volumes, run.energies):
        if e is None:
            continue
        lines.append(f"{v:16.6f} {e / qeout.RY_EV:20.10f} {e:18.8f}")
    fname.write_text("\n".join(lines) + "\n")
    written.append(str(fname))

    txt = out / "EOS.txt"
    txt.write_text(report(run, cell_a) + "\n")
    written.append(str(txt))
    return written


def plot(run: EOSRun, outfile: str = "eos", equation: str = DEFAULT_EQ,
         formats="pdf,png", theme: str = None, size: str = None,
         family: str = None, background: str = None, palette=None,
         usetex: bool = None, width="single", journal: str = "generic",
         aspect: float = 0.80, mono: bool = False, dpi: int = None) -> list:
    """E(V) con la curva ajustada y los residuales en un panel inferior."""
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError as exc:  # pragma: no cover
        raise RuntimeError("matplotlib no está instalado.") from exc

    st = qstyle.apply(theme, size=size, family=family, background=background,
                      palette=palette, usetex=usetex, mono=mono)
    V, E = run.valid()
    if len(V) < 4:
        raise FaltanDatos("no hay puntos suficientes para graficar el ajuste")
    f = (run.fits or fit_all(run)).get(equation)
    if not (f and f.ok):
        raise FaltanDatos(f"el ajuste {equation} no es válido: "
                         f"{f.message if f else 'sin datos'}")

    figsize = qstyle.figure_size(width, journal, aspect)
    fig = plt.figure(figsize=figsize, layout="constrained")
    fig.get_layout_engine().set(w_pad=0.012, h_pad=0.012, hspace=0.0, wspace=0.0)
    gs = fig.add_gridspec(2, 1, height_ratios=[3.0, 1.0], hspace=0.06)
    ax = qstyle.finish_axes(fig.add_subplot(gs[0]))
    axr = qstyle.finish_axes(fig.add_subplot(gs[1], sharex=ax))

    color = qstyle.palette(2, mono=mono)
    func = EQUATIONS[equation][0]
    B0_ev = f.B0 / EV_A3_GPA
    Vfine = np.linspace(V.min(), V.max(), 300)
    Efit = func(Vfine, f.E0, f.V0, B0_ev, f.Bp)

    ax.plot(Vfine, (Efit - f.E0) * 1000.0 / run.natoms, color=color[0],
            lw=st["line"], zorder=2)
    ax.plot(V, (E - f.E0) * 1000.0 / run.natoms, "o", ms=4,
            color=qstyle.INK, mfc=qstyle.CURRENT.get("background", "#FFF"),
            mew=st["line"], zorder=3)
    ax.axvline(f.V0, color=qstyle.INK_FAINT, lw=st["axis_line"], dashes=[3.5, 2.0])
    ax.set_ylabel(r"$E - E_0$ (meV/átomo)")
    ax.tick_params(labelbottom=False)
    ax.annotate(
        f"$V_0$ = {f.V0:.2f} Å$^3$\n$B_0$ = {f.B0:.1f} GPa\n$B_0'$ = {f.Bp:.2f}",
        xy=(0.97, 0.95), xycoords="axes fraction", ha="right", va="top",
        fontsize=st["legend"], color=qstyle.INK,
    )

    resid = (E - func(V, f.E0, f.V0, B0_ev, f.Bp)) * 1000.0 / run.natoms
    axr.axhline(0.0, color=qstyle.INK_FAINT, lw=st["axis_line"])
    axr.plot(V, resid, "o", ms=3, color=color[1])
    axr.set_xlabel(r"$V$ (Å$^3$)")
    axr.set_ylabel("residual\n(meV/át)")
    lim = max(abs(resid).max() * 1.4, 1e-3)
    axr.set_ylim(-lim, lim)

    written = qstyle.save(fig, outfile, formats, dpi=dpi,
                          modulo="ecuación de estado")
    plt.close(fig)
    return written
