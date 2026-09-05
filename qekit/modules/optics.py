# Olla-DFT — command-line toolkit for Quantum ESPRESSO
# Copyright (C) 2026 Jorge Enrique González Sevilla
# SPDX-License-Identifier: AGPL-3.0-or-later
# This program is free software: you can redistribute it and/or modify it
# under the terms of the GNU Affero General Public License as published by the Free
# Software Foundation, either version 3 of the License, or (at your option)
# any later version. See the LICENSE file for details.

"""Propiedades ópticas con epsilon.x: ε(ω), n, k, absorción, Tauc.

Flujo: scf → nscf (malla densa, sin simetría, con bandas vacías) →
epsilon.x. De ε₁(ω) y ε₂(ω) salen todas las funciones ópticas:

    ñ = n + ik = sqrt(ε)          |ε| = sqrt(ε₁² + ε₂²)
    n = sqrt((|ε| + ε₁)/2)        k = sqrt((|ε| − ε₁)/2)
    α(ω) = 2 ω k / c              R = ((n−1)² + k²) / ((n+1)² + k²)

y la gráfica de Tauc, (αhν)^(1/r) contra hν, de cuyo tramo lineal se
extrapola el gap óptico (r = 1/2 para transición directa permitida, r = 2
para indirecta) — el mismo análisis que se hace con un espectro UV-Vis.

Limitaciones físicas que el reporte declara: epsilon.x calcula la respuesta
de partícula independiente (RPA sin campos locales ni excitones) y NO
funciona con pseudopotenciales ultrasoft/PAW — Olla-DFT revisa el tipo de cada
UPF y se niega a preparar el cálculo si no son de norma conservada. Además,
el gap óptico hereda la subestimación del funcional.
"""

import shutil
import subprocess
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

from qekit.core import provenance
from qekit.core.compat import trapezoid
from qekit.core import runner as run_mod
from qekit.core import style as qstyle
from qekit.modules import sweep
from qekit.core.errors import ErrorDeUso

HBAR_C_EV_CM = 1.9732698e-5   # ħc en eV·cm


@dataclass
class OpticsRun:
    prefix: str = ""
    outdir: Path = None
    jobs: list = field(default_factory=list)     # scf, nscf (pw.x)
    wmin: float = 0.0
    wmax: float = 20.0
    nw: int = 800
    intersmear: float = 0.10
    energies: np.ndarray = None                  # eV
    eps1: np.ndarray = None                      # promedio isótropo
    eps2: np.ndarray = None
    eps1_xyz: np.ndarray = None                  # (3, nw)
    eps2_xyz: np.ndarray = None
    scissor: float = 0.0                         # corrimiento aplicado (eV)


def prepare(atoms, outdir: str = "opticas", pseudo_dir: str = None,
            ecutwfc: float = None, ecutrho: float = None,
            kspacing: float = 0.12, insulator: bool = True,
            wmax: float = 20.0, nw: int = 800, intersmear: float = 0.10,
            nbnd_factor: float = 3.0) -> tuple:
    """Escribe scf.in, nscf.in y epsilon.in. Exige pseudos de norma conservada."""
    common = sweep.prepare_common(atoms, pseudo_dir, ecutwfc, ecutrho,
                                  insulator, tarea="optics")

    no_nc = [f"{s} ({p['type'] or 'desconocido'}: {p['filename']})"
             for s, p in common["pseudos"].items()
             if p["found"] and p["type"] != "NC"]
    if no_nc:
        raise ErrorDeUso(
            "epsilon.x solo funciona con pseudopotenciales de NORMA CONSERVADA\n"
            "(los elementos de matriz dipolares no están implementados para "
            "USPP/PAW).\nNo cumplen: " + ", ".join(no_nc) + "\n"
            "Descarga pseudos NC (por ejemplo de PseudoDojo o SG15) y apunta\n"
            "--pseudo-dir a esa carpeta."
        )

    out = Path(outdir); out.mkdir(parents=True, exist_ok=True)
    grid_scf = sweep.default_grid(atoms, None)
    grid_nscf = sweep.default_grid(atoms, kspacing)

    # bandas: bastantes vacías para cubrir la ventana de energía pedida
    from qekit.modules.inputgen import _estimate_nbnd
    nb = _estimate_nbnd(atoms, common["pseudos"])
    nbnd = int(nb * nbnd_factor) if nb else None

    run = OpticsRun(prefix=common["prefix"], outdir=out, wmax=wmax, nw=nw,
                    intersmear=intersmear)

    job_scf = sweep.write_scf_job(atoms, common, out, "scf", grid_scf)
    job_scf.input_file = "scf.in"; job_scf.output_file = "scf.out"
    if sweep.writing_inputs():
        sweep.write_input(out / "scf.in", (out / "pw.in").read_text())

    # nscf: epsilon.x necesita la malla completa sin reducción de simetría
    from qekit.modules import inputgen
    nscf_text = inputgen.build_pw_input(
        atoms=atoms, pseudos=common["pseudos"], calculation="nscf",
        prefix=common["prefix"], pseudo_dir=common["pseudo_dir"],
        ecutwfc=common["ecutwfc"], ecutrho=common["ecutrho"],
        kcard=f"K_POINTS automatic\n  {grid_nscf[0]} {grid_nscf[1]} "
              f"{grid_nscf[2]} 0 0 0\n",
        insulator=insulator, degauss=common["degauss"],
        smearing=common["smearing"], nbnd=nbnd, nosym=True,
    )
    sweep.write_input(out / "nscf.in", nscf_text)
    job_nscf = run_mod.Job(name="nscf", directory=out,
                           input_file="nscf.in", output_file="nscf.out")
    run.jobs = [job_scf, job_nscf]

    sweep.write_input(out / "epsilon.in", 
        "&INPUTPP\n"
        f"  prefix      = '{common['prefix']}'\n"
        "  outdir      = './out'\n"
        "  calculation = 'eps'\n"
        "/\n"
        "&ENERGY_GRID\n"
        "  smeartype   = 'gauss'\n"
        f"  intersmear  = {intersmear:g}\n"
        f"  wmin        = 0.0\n"
        f"  wmax        = {wmax:g}\n"
        f"  nw          = {nw}\n"
        "/\n"
    )

    report = ["--- Propiedades ópticas (epsilon.x) ---",
              f"Estructura: {atoms.get_chemical_formula()}  |  "
              f"pseudos de norma conservada verificados",
              f"Mallas k: scf {grid_scf[0]}x{grid_scf[1]}x{grid_scf[2]}, "
              f"nscf {grid_nscf[0]}x{grid_nscf[1]}x{grid_nscf[2]} "
              "(sin simetría, como exige epsilon.x)",
              f"Bandas: {nbnd or 'automáticas'}  |  ventana 0–{wmax:g} eV, "
              f"{nw} puntos, ensanchamiento {intersmear:g} eV",
              "",
              f"Archivos escritos en '{out.resolve()}': scf.in, nscf.in, "
              "epsilon.in",
              "Orden: pw.x scf -> pw.x nscf -> epsilon.x",
              "",
              "Nota física: respuesta de partícula independiente (RPA sin "
              "campos locales\nni excitones); el gap hereda la subestimación "
              "del funcional."]
    warn = sweep.missing_pseudo_warning(common)
    if warn:
        report.append(warn)
    return run, "\n".join(report)


def run_epsilon(run: OpticsRun, pw_cmd: str = None, nproc: int = None) -> str:
    """Ejecuta epsilon.x (tras scf y nscf) en la carpeta del cálculo."""
    cmd = run_mod.build_command(pw_cmd, nproc)
    exe = Path(cmd[-1]).parent / "epsilon.x" if "/" in cmd[-1] else Path("epsilon.x")
    eps_cmd = cmd[:-1] + [str(exe)]
    if not shutil.which(str(exe)) and not Path(exe).exists():
        raise FileNotFoundError(
            f"no se encontró epsilon.x junto a pw.x ('{exe}'). Compila el "
            "paquete PP de Quantum ESPRESSO (make pp)."
        )
    with open(run.outdir / "epsilon.in") as fin, \
         open(run.outdir / "epsilon.out", "w") as fout:
        proc = subprocess.run(eps_cmd, stdin=fin, stdout=fout,
                              stderr=subprocess.STDOUT, cwd=str(run.outdir))
    if proc.returncode != 0:
        raise RuntimeError("epsilon.x terminó con error; revisa epsilon.out")
    return str(run.outdir / "epsilon.out")


def collect(run: OpticsRun) -> OpticsRun:
    """Lee epsr_*.dat y epsi_*.dat de epsilon.x."""
    def read_eps(stem):
        candidates = sorted(Path(run.outdir).glob(f"{stem}_*.dat"))
        if not candidates:
            raise FileNotFoundError(
                f"no se encontró {stem}_*.dat en {run.outdir}; "
                "¿corrió epsilon.x?"
            )
        data = np.loadtxt(candidates[0], comments="#")
        return data

    dr = read_eps("epsr")
    di = read_eps("epsi")
    run.energies = dr[:, 0]
    run.eps1_xyz = dr[:, 1:4].T
    run.eps2_xyz = di[:, 1:4].T
    run.eps1 = run.eps1_xyz.mean(axis=0)
    run.eps2 = run.eps2_xyz.mean(axis=0)
    return run


# ----------------------------------------------------------------------
# Kramers-Kronig y corrimiento rígido (scissor)
# ----------------------------------------------------------------------
def kramers_kronig(E: np.ndarray, eps2: np.ndarray) -> np.ndarray:
    """ε₁(ω) a partir de ε₂(ω) por Kramers-Kronig.

        ε₁(ω) = 1 + (2/π) P ∫ ω' ε₂(ω') / (ω'² − ω²) dω'

    El polo se trata quitando el punto ω' = ω de la cuadratura (valor
    principal discreto); con la malla uniforme de epsilon.x el error es
    despreciable. Olla-DFT lo usa para el scissor: al mover ε₂ hay que
    recalcular ε₁, no basta con desplazarlo.
    """
    E = np.asarray(E, dtype=float)
    eps2 = np.asarray(eps2, dtype=float)
    out = np.empty_like(E)
    num = E * eps2
    for i, w in enumerate(E):
        den = E ** 2 - w ** 2
        f = np.divide(num, den, out=np.zeros_like(num), where=den != 0.0)
        f[i] = 0.0
        out[i] = 1.0 + (2.0 / np.pi) * trapezoid(f, E)
    return out


def scissor(run: OpticsRun, delta: float) -> OpticsRun:
    """Corrimiento rígido de las bandas de conducción en `delta` eV.

    El funcional subestima el gap; el scissor sube todas las transiciones
    en Δ sin tocar los elementos de matriz. Como ε₂ ∝ |p|²/ω²,

        ε₂'(E) = ε₂(E − Δ) · ((E − Δ)/E)²

    y ε₁ se rehace por Kramers-Kronig (desplazar ε₁ sin más rompería la
    causalidad y daría un ε₁(0) equivocado).

    Δ típico = gap experimental (o GW) − gap del cálculo. Devuelve un
    OpticsRun nuevo; el original no se toca.
    """
    import copy
    if delta == 0.0:
        return run
    new = copy.copy(run)
    E = run.energies
    with np.errstate(divide="ignore", invalid="ignore"):
        fac = np.where(E > 0, ((E - delta) / np.where(E > 0, E, 1.0)) ** 2, 0.0)
    fac = np.where(E - delta > 0, fac, 0.0)

    def shift(e2):
        e2s = np.interp(E - delta, E, e2, left=0.0, right=0.0)
        return e2s * fac

    new.eps2_xyz = np.array([shift(run.eps2_xyz[d]) for d in range(3)])
    new.eps1_xyz = np.array([kramers_kronig(E, new.eps2_xyz[d])
                             for d in range(3)])
    new.eps2 = new.eps2_xyz.mean(axis=0)
    new.eps1 = new.eps1_xyz.mean(axis=0)
    new.scissor = float(delta)
    return new


# ----------------------------------------------------------------------
# Funciones ópticas derivadas
# ----------------------------------------------------------------------
def derived(run: OpticsRun) -> dict:
    e1, e2 = run.eps1, run.eps2
    mod = np.sqrt(e1 ** 2 + e2 ** 2)
    n = np.sqrt(np.maximum((mod + e1) / 2.0, 0.0))
    k = np.sqrt(np.maximum((mod - e1) / 2.0, 0.0))
    # α = 2 ω k / c, con ħω en eV -> α en cm⁻¹: α = 2 k E / (ħc)
    alpha = 2.0 * k * run.energies / HBAR_C_EV_CM
    R = ((n - 1.0) ** 2 + k ** 2) / ((n + 1.0) ** 2 + k ** 2)
    return {"n": n, "k": k, "alpha": alpha, "R": R}


def tauc_gap(run: OpticsRun, kind: str = "direct",
             fit_window: float = 0.6, max_span: float = 1.5):
    """Gap óptico por extrapolación de Tauc.

    Se construye (αhν)^(1/r) (r = 1/2 directa, r = 2 indirecta) y se ajusta
    una recta sobre el PRIMER frente de absorción; el corte con cero es el
    gap. Buscar la pendiente máxima global no sirve: el espectro de un
    cálculo DFT tiene picos interbanda muy marcados a energías altas que
    son más empinados que el borde y darían un gap sin sentido.

    Devuelve (gap, pendiente, ventana_usada, curva_y).

    Aviso físico: epsilon.x no incluye transiciones asistidas por fonones,
    así que en un semiconductor de gap indirecto (Si, por ejemplo) ε₂ es
    cero por debajo del gap DIRECTO y el ajuste 'indirect' no reproduce el
    gap indirecto real; ahí el número que sale del ajuste 'direct' es el
    comparable.
    """
    d = derived(run)
    r = 0.5 if kind == "direct" else 2.0
    y = (d["alpha"] * run.energies) ** (1.0 / r)
    E = run.energies
    ymax = float(y.max())
    if ymax <= 0 or len(E) < 10:
        return None, None, None, y

    # suavizado ligero: el frente trae el ruido de la malla de k
    win = max(3, int(round(0.05 / (E[1] - E[0]))) | 1)
    ker = np.ones(win) / win
    ys = np.convolve(y, ker, mode="same")
    dy = np.gradient(ys, E)

    # Primer frente: desde donde y despega hasta que la pendiente vuelve a
    # caer; dentro de ese tramo, el punto más empinado.
    #
    # El umbral NO puede ser una fracción del máximo global. Si más arriba
    # hay un pico interbanda mucho más intenso que el borde, ese umbral
    # queda por encima del borde entero y la detección se va al pico —
    # exactamente el error que este ajuste existe para evitar. En su lugar
    # se mide el piso de ruido en la zona sub-gap (el arranque del
    # espectro, donde no hay absorción) y se despega de ahí.
    n_piso = max(3, len(E) // 100)
    piso_ruido = float(np.max(y[:n_piso]))
    # la escala de referencia es la MEDIANA de la curva donde hay
    # absorción, no su máximo: la mediana no la mueve un pico estrecho por
    # intenso que sea, y el máximo sí
    nz = y[y > 0]
    escala = float(np.median(nz)) if nz.size else ymax
    piso = max(2.0 * piso_ruido, 1e-3 * escala)
    above = np.where((ys > piso) & (E > 0.1))[0]
    if len(above) == 0:
        return None, None, None, y
    i_on = int(above[0])
    # Fin del frente: el primer máximo local después del despegue, o como
    # mucho `max_span` eV por encima del inicio.
    # El tope por energía es imprescindible: la curva de Tauc crece como
    # E^4 aunque alpha se aplane, así que en un espectro suave no hay
    # ningún máximo local hasta el siguiente pico interbanda — y sin tope
    # el "punto más empinado del primer frente" acaba siendo ese pico.
    tope = float(E[i_on]) + max_span
    i_end = int(np.searchsorted(E, tope))
    i_end = min(max(i_end, i_on + 3), len(E) - 1)
    for i in range(i_on + 2, i_end):
        if (ys[i] >= ys[i - 1] and ys[i] >= ys[i + 1]
                and ys[i] > 3.0 * ys[i_on]):
            i_end = i
            break
    seg = slice(i_on, max(i_end + 1, i_on + 5))
    i0 = i_on + int(np.argmax(dy[seg]))

    half = fit_window / 2.0
    sel = (E > E[i0] - half) & (E < E[i0] + half) & (y > 0)
    if sel.sum() < 4:
        return None, None, None, y
    slope, intercept = np.polyfit(E[sel], y[sel], 1)
    if slope <= 0:
        return None, None, None, y
    gap = -intercept / slope
    if not (E[0] <= gap <= E[-1]):
        return None, None, None, y
    return float(gap), float(slope), (float(E[sel].min()), float(E[sel].max())), y


def report(run: OpticsRun) -> str:
    d = derived(run)
    e1_0 = float(run.eps1[np.searchsorted(run.energies, 0.05)])
    i2 = int(np.argmax(run.eps2))
    gd, *_ = tauc_gap(run, "direct")
    gi, *_ = tauc_gap(run, "indirect")
    lines = ["--- Funciones ópticas ---"]
    if run.scissor:
        lines.append(f"Scissor aplicado: +{run.scissor:.2f} eV "
                     "(ε₁ rehecho por Kramers-Kronig)")
    lines += [f"ε₁(0) (constante dieléctrica electrónica): {e1_0:.2f}",
             f"n(0) = {float(d['n'][1]):.3f}",
             f"máximo de ε₂ en {run.energies[i2]:.2f} eV "
             f"(ε₂ = {run.eps2[i2]:.1f})",
             "",
             "Gap óptico por Tauc:",
             f"  directa   (αhν)²   : {gd:.2f} eV" if gd else
             "  directa   (αhν)²   : no se pudo ajustar",
             f"  indirecta (αhν)^½  : {gi:.2f} eV" if gi else
             "  indirecta (αhν)^½  : no se pudo ajustar",
             ""]
    if run.scissor:
        lines.append("Recuerda: sigue siendo RPA de partícula independiente "
                     "(sin campos locales\nni excitones); el scissor solo "
                     "corrige la posición del gap.")
    else:
        lines.append("Recuerda: RPA de partícula independiente y gap del "
                     "funcional. Para comparar\ncon UV-Vis conviene "
                     "'--scissor Δ' con Δ = gap experimental (o GW) − gap "
                     "del\ncálculo; Olla-DFT desplaza ε₂ y rehace ε₁ por "
                     "Kramers-Kronig.")
    return "\n".join(lines)


#: Columnas de OPTICS.dat, en orden. Quien lea el archivo (por ejemplo
#: 'tddft --compare') debe pedir la columna por nombre con `read_optics_dat`,
#: no por posición.
OPTICS_COLUMNS = ("E(eV)", "eps1", "eps2", "n", "k", "alpha(1/cm)", "R")


def export(run: OpticsRun, outdir: str = ".") -> list:
    out = Path(outdir); out.mkdir(parents=True, exist_ok=True)
    d = derived(run)
    f = out / "OPTICS.dat"
    anchos = (10, 12, 12, 10, 10, 14, 10)
    header = (provenance.header_plain(
                  "propiedades ópticas",
                  {"scissor_eV": run.scissor or 0.0,
                   "intersmear_eV": run.intersmear,
                   "promedio": "isótropo (x,y,z)"},
                  titulo="Funciones ópticas") + "\n"
              + " ".join(f"{c:>{w}s}" for c, w in zip(OPTICS_COLUMNS, anchos)))
    np.savetxt(f, np.column_stack([run.energies, run.eps1, run.eps2,
                                   d["n"], d["k"], d["alpha"], d["R"]]),
               fmt="%14.6e", header=header, comments="# ")
    return [str(f)]


def read_optics_dat(path) -> dict:
    """Lee OPTICS.dat y devuelve {nombre de columna: array}.

    Los nombres se toman de la última línea de comentario del encabezado
    (la que escribe `export`); si el archivo no la trae, se asume el orden
    de `OPTICS_COLUMNS`. Así 'tddft --compare' recibe α por su nombre y no
    la última columna, que es la reflectividad.
    """
    path = Path(path)
    nombres = None
    with open(path, encoding="utf-8", errors="ignore") as fh:
        for linea in fh:
            s = linea.strip()
            if not s.startswith("#"):
                break
            campos = s.lstrip("#").split()
            if campos and "E(eV)" in campos[0]:
                nombres = campos
    datos = np.loadtxt(path, comments="#", ndmin=2)
    if nombres is None or len(nombres) != datos.shape[1]:
        nombres = list(OPTICS_COLUMNS[:datos.shape[1]])
    if datos.shape[1] < 2:
        raise ErrorDeUso(
            f"'{path}' no parece un OPTICS.dat: tiene {datos.shape[1]} "
            "columna(s) y hacen falta al menos E y α.")
    return {n: datos[:, i] for i, n in enumerate(nombres)}


def plot(run: OpticsRun, outfile: str = "opticas", formats="pdf,png",
         theme: str = None, size: str = None, family: str = None,
         background: str = None, palette=None, usetex: bool = None,
         width="double", journal: str = "generic", aspect: float = 0.40,
         mono: bool = False, dpi: int = None, tauc_kind: str = "direct") -> list:
    """Panel (a) ε₁/ε₂, (b) absorción, (c) Tauc con el ajuste."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    st = qstyle.apply(theme, size=size, family=family, background=background,
                      palette=palette, usetex=usetex, mono=mono)
    d = derived(run)
    figsize = qstyle.figure_size(width, journal, aspect)
    fig = plt.figure(figsize=figsize, layout="constrained")
    gs = fig.add_gridspec(1, 3, wspace=0.08)
    axes = [qstyle.finish_axes(fig.add_subplot(gs[i])) for i in range(3)]
    colors = qstyle.palette(4, mono=mono)
    E = run.energies

    ax = axes[0]
    ax.plot(E, run.eps1, color=colors[0], lw=st["line"],
            label=r"$\varepsilon_1$")
    ax.plot(E, run.eps2, color=colors[1], lw=st["line"],
            **({"dashes": qstyle.dash(1)} if qstyle.use_dashes(2, "auto", mono) else {}),
            label=r"$\varepsilon_2$")
    ax.axhline(0, color=qstyle.INK_FAINT, lw=st["axis_line"])
    ax.set_xlabel(r"$\hbar\omega$ (eV)"); ax.set_ylabel(r"$\varepsilon(\omega)$")
    ax.set_xlim(0, E.max()); ax.legend()
    qstyle.panel_label(ax, "(a)")

    ax = axes[1]
    ax.plot(E, d["alpha"] / 1e5, color=colors[2], lw=st["line"])
    ax.set_xlabel(r"$\hbar\omega$ (eV)")
    ax.set_ylabel(r"$\alpha$ ($10^{5}$ cm$^{-1}$)")
    ax.set_xlim(0, E.max())
    qstyle.panel_label(ax, "(b)")

    ax = axes[2]
    gap, slope, window, y = tauc_gap(run, tauc_kind)
    exp_txt = "2" if tauc_kind == "direct" else "1/2"
    ax.plot(E, y, color=colors[0], lw=st["line"])
    if gap:
        xs = np.linspace(gap, window[1] + 0.4, 20)
        ax.plot(xs, slope * (xs - gap), color=colors[1], lw=st["line"],
                dashes=[4, 1.6])
        ax.plot([gap], [0], "o", ms=4, color=colors[1],
                mec=st["background"], mew=0.6, zorder=5)
        ax.set_xlim(max(0, gap - 1.2), min(E.max(), gap + 2.6))
        sel = (E >= ax.get_xlim()[0]) & (E <= ax.get_xlim()[1])
        ax.set_ylim(0, float(y[sel].max()) * 1.08)
        ax.annotate(rf"$E_g^{{\mathrm{{opt}}}}$ = {gap:.2f} eV",
                    xy=(0.05, 0.93), xycoords="axes fraction",
                    ha="left", va="top", fontsize=st["legend"])
    else:
        ax.set_xlim(0, E.max())
    ax.set_xlabel(r"$h\nu$ (eV)")
    ax.set_ylabel(rf"$(\alpha h\nu)^{{{exp_txt}}}$ (u. arb.)")
    ax.set_yticks([])
    qstyle.panel_label(ax, "(c)")

    written = qstyle.save(fig, outfile, formats, dpi=dpi,
                          modulo="propiedades ópticas")
    plt.close(fig)
    return written
