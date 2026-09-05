# Olla-DFT — command-line toolkit for Quantum ESPRESSO
# Copyright (C) 2026 Jorge Enrique González Sevilla
# SPDX-License-Identifier: AGPL-3.0-or-later
# This program is free software: you can redistribute it and/or modify it
# under the terms of the GNU Affero General Public License as published by the Free
# Software Foundation, either version 3 of the License, or (at your option)
# any later version. See the LICENSE file for details.

"""Densidad de carga, potencial electrostático y función trabajo (pp.x).

pp.x extrae campos escalares del cálculo terminado. Aquí se generan sus
inputs, se leen los resultados y se calculan las dos cantidades que más se
piden:

- **promedio planar** de un campo a lo largo de un eje (perfil z), con su
  promedio móvil "macroscópico" sobre una distancia típica interplanar;
- **función trabajo** Φ = V_vacío − E_F para una superficie o monocapa con
  vacío: el potencial electrostático se promedia en el plano, se lee su
  valor en la meseta del vacío y se le resta la energía de Fermi del scf.

Formatos: pp.x escribe con plot_num el campo elegido (0 densidad, 11
potencial electrostático V_H + V_ion, ...) y con output_format=6 un archivo
.cube legible. El lector de cube está aquí mismo, sin dependencias.
"""

import shutil
import subprocess
from dataclasses import dataclass
from pathlib import Path

import numpy as np

from qekit.core import qeout
from qekit.core import runner as run_mod
from qekit.core import style as qstyle
from qekit.modules import sweep
from qekit.core.errors import ErrorDeUso, FaltanDatos

RY_EV = qeout.RY_EV
BOHR = qeout.BOHR_ANG

PLOTS = {
    "density": (0, "densidad de carga (e/bohr³)"),
    "potential": (11, "potencial electrostático V_bare + V_H (Ry)"),
    "vtotal": (1, "potencial total V_bare + V_H + V_xc (Ry)"),
    "elf": (8, "función de localización electrónica (ELF)"),
    "spin": (6, "densidad de espín (up − down)"),
}


def build_pp_input(prefix: str, plot: str, cube_name: str) -> str:
    if plot not in PLOTS:
        raise ErrorDeUso(f"campo desconocido '{plot}'. Opciones: {', '.join(PLOTS)}")
    num, _desc = PLOTS[plot]
    return (
        "&INPUTPP\n"
        f"  prefix       = '{prefix}'\n"
        "  outdir       = './out'\n"
        f"  filplot      = '{cube_name}.pp'\n"
        f"  plot_num     = {num}\n"
        "/\n"
        "&PLOT\n"
        "  iflag        = 3\n"
        "  output_format = 6\n"
        f"  fileout      = '{cube_name}.cube'\n"
        "/\n"
    )


def run_pp(calc_dir, plot: str, cube_name: str, pw_cmd: str = None,
           nproc: int = None, prefix: str = None) -> Path:
    """Escribe el input de pp.x, lo ejecuta en la carpeta del cálculo y
    devuelve la ruta del .cube resultante."""
    calc_dir = Path(calc_dir)
    if prefix is None:
        prefix = qeout.read_xml(str(calc_dir)).prefix
    stem = f"pp_{plot}"
    sweep.write_input(calc_dir / f"{stem}.in",
                      build_pp_input(prefix, plot, cube_name))

    cmd = run_mod.build_command(pw_cmd, nproc)
    exe = Path(cmd[-1]).parent / "pp.x" if "/" in cmd[-1] else Path("pp.x")
    if not shutil.which(str(exe)) and not Path(exe).exists():
        raise FileNotFoundError(
            f"no se encontró pp.x junto a pw.x ('{exe}'). "
            "Compila el paquete PP de Quantum ESPRESSO (make pp)."
        )
    pp_cmd = cmd[:-1] + [str(exe)]
    with open(calc_dir / f"{stem}.in") as fin,          open(calc_dir / f"{stem}.out", "w") as fout:
        proc = subprocess.run(pp_cmd, stdin=fin, stdout=fout,
                              stderr=subprocess.STDOUT, cwd=str(calc_dir))
    out_text = (calc_dir / f"{stem}.out").read_text(errors="ignore")
    cube = calc_dir / f"{cube_name}.cube"
    if proc.returncode != 0 or "JOB DONE" not in out_text or not cube.exists():
        raise RuntimeError(run_mod.failure_message(
            "pp.x", calc_dir / f"{stem}.out", out_text))
    return cube


# ----------------------------------------------------------------------
# Lector de archivos cube
# ----------------------------------------------------------------------
@dataclass
class CubeData:
    origin: np.ndarray = None      # Å
    axes: np.ndarray = None        # (3,3) Å: paso de la malla en cada eje
    shape: tuple = None
    data: np.ndarray = None        # (n1, n2, n3)
    natoms: int = 0

    @property
    def cell(self) -> np.ndarray:
        return self.axes * np.array(self.shape)[:, None]


def read_cube(path: str) -> CubeData:
    with open(path) as fh:
        fh.readline(); fh.readline()
        parts = fh.readline().split()
        natoms = int(parts[0])
        origin = np.array([float(x) for x in parts[1:4]])
        shape, axes = [], []
        for _ in range(3):
            p = fh.readline().split()
            n = int(p[0])
            vec = np.array([float(x) for x in p[1:4]])
            # convenio cube: n > 0 => unidades atómicas (bohr)
            if n > 0:
                vec = vec * BOHR
            shape.append(abs(n)); axes.append(vec)
        for _ in range(natoms):
            fh.readline()
        # np.fromstring con separador está deprecado desde numpy 1.14
        vals = np.array(fh.read().split(), dtype=float)
    n1, n2, n3 = shape
    if vals.size != n1 * n2 * n3:
        raise FaltanDatos(f"'{path}': se esperaban {n1*n2*n3} valores y hay {vals.size}")
    return CubeData(origin=origin * BOHR, axes=np.array(axes),
                    shape=(n1, n2, n3), data=vals.reshape((n1, n2, n3)),
                    natoms=natoms)


# ----------------------------------------------------------------------
# Promedio planar y función trabajo
# ----------------------------------------------------------------------
def planar_average(cube: CubeData, axis: int = 2) -> tuple:
    """(z en Å, promedio del campo en cada plano perpendicular a `axis`)."""
    keep = axis
    other = tuple(k for k in range(3) if k != keep)
    prof = cube.data.mean(axis=other)
    step = np.linalg.norm(cube.axes[keep])
    z = np.arange(cube.shape[keep]) * step
    return z, prof


def macroscopic_average(z: np.ndarray, prof: np.ndarray, window: float) -> np.ndarray:
    """Promedio móvil periódico sobre una ventana en Å (promedio macroscópico)."""
    if window <= 0 or len(z) < 3:
        return prof
    step = z[1] - z[0]
    n = max(1, int(round(window / step)))
    kernel = np.ones(n) / n
    ext = np.concatenate([prof[-n:], prof, prof[:n]])
    sm = np.convolve(ext, kernel, mode="same")
    return sm[n:-n]


@dataclass
class WorkFunction:
    phi: float = None            # eV
    v_vacuum: float = None       # eV
    fermi: float = None          # eV
    vacuum_z: tuple = None       # tramo usado como vacío (Å)
    flatness: float = None       # variación del potencial en la meseta (eV)
    z: np.ndarray = None
    profile: np.ndarray = None   # potencial planar en eV
    axis: int = 2


def vacuum_window(cube: CubeData, positions, axis: int = 2,
                  fraction: float = 0.2) -> list:
    """Índices de plano del tramo central (`fraction`) de la región de vacío.

    La región de vacío es el hueco más ancho sin átomos a lo largo de
    `axis`, medido con periodicidad sobre las posiciones (Å) dadas. Se
    devuelve la lista de índices del perfil planar que caen en su
    `fraction` central (20 % por omisión), de modo que ni la cola de la
    densidad de una cara ni la de la otra entren en la meseta.
    """
    npts = int(cube.shape[axis])
    pos = np.asarray(positions, dtype=float).reshape(-1, 3)
    frac = ((pos - cube.origin) @ np.linalg.inv(cube.cell))[:, axis] % 1.0
    frac = np.sort(frac)
    if frac.size == 1:
        inicio, hueco = float(frac[0]), 1.0
    else:
        huecos = [(frac[i + 1] - frac[i], frac[i]) for i in range(frac.size - 1)]
        huecos.append((1.0 - frac[-1] + frac[0], frac[-1]))
        hueco, inicio = max(huecos)
    centro = (inicio + hueco / 2.0) % 1.0
    i_c = int(round(centro * npts)) % npts
    half = max(2, int(round(0.5 * fraction * hueco * npts)))
    return [(i_c + k) % npts for k in range(-half, half + 1)]


def work_function(cube: CubeData, fermi_ev: float, axis: int = 2,
                  positions=None) -> WorkFunction:
    """Φ = V_vacío − E_F a partir del cube del potencial (plot_num=11, Ry).

    El nivel de vacío es el promedio del potencial planar en la meseta de
    vacío, y la "planitud" es cuánto varía el potencial dentro de ella: si
    no es plana, el vacío es insuficiente (o hay un dipolo neto) y el
    número no es confiable.

    Dónde está la meseta depende de lo que se sepa de la estructura:

    - con `positions` (Å, las del cálculo), la meseta es el 20 % central
      del hueco más ancho sin átomos a lo largo de `axis`, que es la
      región de vacío de verdad;
    - sin posiciones, se toma una ventana de ±10 % de la celda alrededor
      del MÁXIMO del potencial planar (el vacío es la meseta alta). Es un
      criterio razonable pero ciego: con poco vacío la ventana puede
      pisar la cola del potencial atómico.
    """
    z, prof_ry = planar_average(cube, axis)
    prof = prof_ry * RY_EV
    npts = len(prof)

    if positions is not None and len(positions):
        idx = vacuum_window(cube, positions, axis=axis, fraction=0.2)
    else:
        i_vac = int(np.argmax(prof))
        half = max(2, npts // 10)
        idx = [(i_vac + k) % npts for k in range(-half, half + 1)]
    plateau = prof[idx]
    v_vac = float(plateau.mean())
    flat = float(plateau.max() - plateau.min())

    return WorkFunction(
        phi=v_vac - fermi_ev, v_vacuum=v_vac, fermi=fermi_ev,
        vacuum_z=(float(z[idx[0]]), float(z[idx[-1]])),
        flatness=flat, z=z, profile=prof, axis=axis,
    )


def report_wf(wf: WorkFunction) -> str:
    lines = ["--- Función trabajo ---",
             f"V_vacío = {wf.v_vacuum:10.4f} eV   (meseta plana a ±{wf.flatness/2:.3f} eV, "
             f"evaluada en z = {wf.vacuum_z[0]:.2f}..{wf.vacuum_z[1]:.2f} Å)",
             f"E_Fermi = {wf.fermi:10.4f} eV",
             "",
             f"Φ = V_vacío − E_F = {wf.phi:.3f} eV"]
    if wf.flatness > 0.05:
        lines.append("\nAVISO: la meseta de vacío varía más de 0.05 eV. El vacío es\n"
                     "insuficiente o hay un dipolo neto; aumenta el vacío (o usa una\n"
                     "losa simétrica) antes de confiar en este valor.")
    return "\n".join(lines)


def export_wf(wf: WorkFunction, outdir: str = ".") -> list:
    """Deja el perfil y el número en un archivo, como todos los demás módulos.

    Hasta la 0.25 este comando solo escribía la figura, así que su resultado
    no se podía cruzar contra el de nadie —en particular contra el de ESM,
    que llega a la misma función trabajo por un camino completamente
    distinto—. Un número que solo existe dentro de un PNG no es un resultado
    reutilizable.
    """
    from qekit.core import provenance
    out = Path(outdir); out.mkdir(parents=True, exist_ok=True)
    f = out / "WF.dat"
    cab = [provenance.header("función trabajo desde el potencial planar",
                             {"eje": "abc"[wf.axis]}),
           f"# Phi_eV = {wf.phi:.6f}",
           f"# V_vacio_eV = {wf.v_vacuum:.6f}",
           f"# E_Fermi_eV = {wf.fermi:.6f}",
           f"# planitud_eV = {wf.flatness:.6f}",
           f"# {'z(A)':>12s} {'V_planar(eV)':>16s}"]
    filas = [f"{zz:14.6f} {vv:16.6f}" for zz, vv in zip(wf.z, wf.profile)]
    f.write_text("\n".join(cab + filas) + "\n", encoding="utf-8")
    return [str(f)]


def plot_profile(wf: WorkFunction, outfile: str = "funcion_trabajo",
                 formats="pdf,png", theme: str = None, size: str = None,
                 family: str = None, background: str = None, palette=None,
                 usetex: bool = None, width="single", journal: str = "generic",
                 aspect: float = 0.62, mono: bool = False, dpi: int = None) -> list:
    """Perfil del potencial planar con E_F, V_vacío y Φ anotados."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    st = qstyle.apply(theme, size=size, family=family, background=background,
                      palette=palette, usetex=usetex, mono=mono)
    fig, ax = qstyle.new_figure(width, journal, aspect)
    color = qstyle.palette(2, mono=mono)

    ax.plot(wf.z, wf.profile, color=color[0], lw=st["line"])
    ax.axhline(wf.fermi, color=color[1], lw=st["line"] * 0.9, dashes=[4, 1.6])
    ax.axhline(wf.v_vacuum, color=qstyle.INK_FAINT, lw=st["axis_line"],
               dashes=[2, 1.5])
    x_txt = wf.z[-1] * 0.985
    ax.annotate(r"$E_\mathrm{F}$", xy=(x_txt, wf.fermi), ha="right",
                va="bottom", fontsize=st["legend"], color=color[1])
    ax.annotate(r"$V_\mathrm{vac}$", xy=(x_txt, wf.v_vacuum), ha="right",
                va="bottom", fontsize=st["legend"], color=qstyle.INK_SOFT)
    zmid = 0.5 * (wf.vacuum_z[0] + wf.vacuum_z[1])
    ax.annotate(
        rf"$\Phi$ = {wf.phi:.2f} eV",
        xy=(zmid, 0.5 * (wf.fermi + wf.v_vacuum)), ha="center", va="center",
        fontsize=st["legend"],
        bbox=dict(facecolor=st["background"], alpha=0.85, edgecolor="none", pad=1.2),
    )
    ax.set_xlabel(f"z ({qstyle.angstrom()})")
    ax.set_ylabel("potencial planar (eV)")
    ax.set_xlim(wf.z[0], wf.z[-1])
    written = qstyle.save(fig, outfile, formats, dpi=dpi,
                          modulo="campos (pp.x)")
    plt.close(fig)
    return written


def plot_density_profile(z, prof, outfile: str = "densidad_planar",
                         label: str = "densidad de carga",
                         formats="pdf,png", theme: str = None, size: str = None,
                         family: str = None, background: str = None, palette=None,
                         usetex: bool = None, width="single",
                         journal: str = "generic", aspect: float = 0.6,
                         mono: bool = False, dpi: int = None) -> list:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    st = qstyle.apply(theme, size=size, family=family, background=background,
                      palette=palette, usetex=usetex, mono=mono)
    fig, ax = qstyle.new_figure(width, journal, aspect)
    color = qstyle.palette(1, mono=mono)[0]
    ax.plot(z, prof, color=color, lw=st["line"])
    ax.fill_between(z, 0, prof, color=color, alpha=0.12, lw=0)
    ax.set_xlabel(f"z ({qstyle.angstrom()})")
    ax.set_ylabel(qstyle.tex_safe(label))
    ax.set_xlim(z[0], z[-1])
    written = qstyle.save(fig, outfile, formats, dpi=dpi,
                          modulo="campos (pp.x)")
    plt.close(fig)
    return written
