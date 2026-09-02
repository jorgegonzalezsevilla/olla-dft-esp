# Olla-DFT — command-line toolkit for Quantum ESPRESSO
# Copyright (C) 2026 Jorge Enrique González Sevilla
# SPDX-License-Identifier: GPL-3.0-or-later
# This program is free software: you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by the Free
# Software Foundation, either version 3 of the License, or (at your option)
# any later version. See the LICENSE file for details.

"""Fonones por DFPT: dispersión, DOS, termodinámica e IR.

Cadena completa de Quantum ESPRESSO orquestada por Olla-DFT:

    pw.x (scf muy convergido) -> ph.x (malla de q, DFPT)
        -> q2r.x (constantes de fuerza reales)
        -> matdyn.x (dispersión por el camino de alta simetría)
        -> matdyn.x (DOS de fonones en malla densa)

y de la DOS salen las funciones termodinámicas armónicas: energía de punto
cero, F(T), S(T) y C_v(T).

Modo Γ (--gamma): un solo punto q con dynmat.x, que aplica la regla de suma
acústica y reporta frecuencias y actividades IR — lo comparable con un FTIR.
Para aislantes añade epsil=.true. (tensor dieléctrico y cargas de Born, con
la separación LO–TO en Γ).

Detalles que importan y Olla-DFT impone:
- el scf lleva conv_thr = 1e-12: la DFPT deriva de la densidad y hereda su
  ruido multiplicado;
- la estructura debe estar RELAJADA (frecuencias imaginarias espurias si no);
- la regla de suma acústica (asr='simple') se aplica en q2r/matdyn/dynmat,
  y el reporte marca cualquier frecuencia imaginaria que sobreviva.
"""

import re
import shutil
import subprocess
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

from qekit.core import kpoints
from qekit.core import provenance
from qekit.core.compat import trapezoid
from qekit.core import runner as run_mod
from qekit.core import style as qstyle
from qekit.modules import sweep
from qekit.core.errors import ErrorDeUso, FaltanDatos

CM1_TO_THZ = 0.0299792458
CM1_TO_EV = 1.239841984e-4
KB_EV = 8.617333262e-5


@dataclass
class PhononRun:
    prefix: str = ""
    outdir: Path = None
    qgrid: tuple = None
    gamma_only: bool = False
    epsil: bool = False
    jobs: list = field(default_factory=list)      # solo el scf (pw.x)
    kpath = None                                   # kpoints.KPath
    # resultados
    band_q: np.ndarray = None                      # (nq, 3) crystal
    band_freqs: np.ndarray = None                  # (nq, nmodes) cm^-1
    qdist: np.ndarray = None
    labels: list = field(default_factory=list)
    dos_w: np.ndarray = None                       # cm^-1
    dos: np.ndarray = None                         # estados/cm^-1
    band_entries: list = field(default_factory=list)  # [(etiqueta, npts)]
    gamma_freqs: list = field(default_factory=list)   # [(cm-1, actividad IR o None)]
    raman: bool = False
    modes: list = field(default_factory=list)   # dicts con IR, Raman y depol


def _exe(name: str, pw_cmd: str = None, nproc: int = None) -> list:
    cmd = run_mod.build_command(pw_cmd, nproc)
    exe = Path(cmd[-1]).parent / name if "/" in cmd[-1] else Path(name)
    if not shutil.which(str(exe)) and not Path(exe).exists():
        raise FileNotFoundError(
            f"no se encontró {name} junto a pw.x. Compílalo con 'make ph'."
        )
    return cmd[:-1] + [str(exe)]


def _step_done(workdir: Path, stem: str) -> bool:
    try:
        return "JOB DONE" in (workdir / f"{stem}.out").read_text(errors="ignore")
    except OSError:
        return False


def _run_tool(tool_cmd: list, workdir: Path, stem: str):
    with open(workdir / f"{stem}.in") as fin, \
         open(workdir / f"{stem}.out", "w") as fout:
        proc = subprocess.run(tool_cmd, stdin=fin, stdout=fout,
                              stderr=subprocess.STDOUT, cwd=str(workdir))
    text = (workdir / f"{stem}.out").read_text(errors="ignore")
    if proc.returncode != 0 or "JOB DONE" not in text:
        raise RuntimeError(run_mod.failure_message(
            stem, workdir / f"{stem}.out", text))


# ----------------------------------------------------------------------
# Preparación
# ----------------------------------------------------------------------
def prepare(atoms, outdir: str = "fonones", pseudo_dir: str = None,
            ecutwfc: float = None, ecutrho: float = None,
            kspacing: float = None, insulator: bool = True,
            qgrid: tuple = None, qspacing: float = 0.6,
            gamma_only: bool = False, epsil: bool = None,
            raman: bool = False,
            tr2_ph: float = 1e-14, band_points: int = 30,
            dos_grid: tuple = (12, 12, 12),
            degauss: float = None, smearing: str = None) -> tuple:
    """Escribe scf.in, ph.in y los inputs de post-proceso.

    Sobre `insulator`: el valor por omisión de esta función es True, pero
    la CLI (`olla-dft phonons`) pasa siempre `args.insulator`, que es False
    salvo que se dé `--insulator`. O sea que desde la línea de comandos el
    scf va POR DEFECTO con smearing (seguro para metales, y también válido
    para un aislante con gap) y `epsil` NO se activa: sin `--insulator` no
    hay tensor dieléctrico, cargas de Born ni separación LO–TO en Γ. Con
    `--insulator` el scf lleva occupations='fixed' y ph.x calcula epsil.
    `--raman` activa epsil por su cuenta, porque lraman lo necesita.
    """
    from qekit.core import structure as struct_mod
    from qekit.modules import inputgen

    # marco estándar: el camino de q de la dispersión lo exige igual que
    # el de bandas electrónicas
    atoms = struct_mod.primitive(atoms)
    common = sweep.prepare_common(atoms, pseudo_dir, ecutwfc, ecutrho,
                                  insulator, tarea="fonones")
    out = Path(outdir); out.mkdir(parents=True, exist_ok=True)

    if epsil is None:
        epsil = insulator
    grid_scf = sweep.default_grid(atoms, kspacing)
    if qgrid is None and not gamma_only:
        qgrid = kpoints.kgrid_from_spacing(atoms, qspacing)
    if raman and not gamma_only:
        gamma_only = True          # el Raman de Olla-DFT se calcula en Gamma
    if raman:
        no_nc = [f"{sy} ({p['filename']})"
                 for sy, p in common["pseudos"].items()
                 if p["found"] and p["type"] != "NC"]
        if no_nc:
            raise ErrorDeUso(
                "el cálculo Raman (lraman) solo funciona con "
                "pseudopotenciales de NORMA CONSERVADA, y estos no lo "
                "son:\n  " + "\n  ".join(no_nc) +
                "\n\nCon ultrasoft o PAW, ph.x se detiene o devuelve "
                "tensores sin sentido.")
    run = PhononRun(prefix=common["prefix"], outdir=out,
                    qgrid=None if gamma_only else tuple(qgrid),
                    gamma_only=gamma_only, epsil=epsil or raman,
                    raman=raman)

    scf_text = inputgen.build_pw_input(
        atoms=atoms, pseudos=common["pseudos"], calculation="scf",
        prefix=common["prefix"], pseudo_dir=common["pseudo_dir"],
        ecutwfc=common["ecutwfc"], ecutrho=common["ecutrho"],
        kcard=f"K_POINTS automatic\n  {grid_scf[0]} {grid_scf[1]} "
              f"{grid_scf[2]} 0 0 0\n",
        insulator=insulator,
        degauss=degauss if degauss is not None else common["degauss"],
        smearing=smearing or common["smearing"], conv_thr=1e-12,
    )
    sweep.write_input(out / "scf.in", scf_text)
    job = run_mod.Job(name="scf", directory=out,
                      input_file="scf.in", output_file="scf.out")
    run.jobs = [job]

    ph = ["fonones de Olla-DFT", "&inputph",
          f"  prefix   = '{common['prefix']}'",
          "  outdir   = './out'",
          "  fildyn   = 'dyn'",
          f"  tr2_ph   = {tr2_ph:g}"]
    if epsil or raman:
        ph.append("  epsil    = .true.")
    if raman:
        # lraman activa la respuesta de tercer orden (2n+1): da los tensores
        # Raman y de ahí las actividades. Solo funciona con pseudos de norma
        # conservada, y es MUCHO más caro que un fonón normal.
        ph.append("  lraman   = .true.")
        ph.append("  trans    = .true.")
    if gamma_only:
        ph += ["/", "0.0 0.0 0.0"]
    else:
        ph += ["  ldisp    = .true.",
               f"  nq1      = {qgrid[0]}",
               f"  nq2      = {qgrid[1]}",
               f"  nq3      = {qgrid[2]}",
               "/"]
    sweep.write_input(out / "ph.in", "\n".join(ph) + "\n")

    if gamma_only:
        dyn_in = ["&input", "  fildyn = 'dyn'", "  asr    = 'simple'",
                  "  filout = 'dynmat.modos'", "/"]
        sweep.write_input(out / "dynmat.in", "\n".join(dyn_in) + "\n")
    else:
        sweep.write_input(out / "q2r.in", 
            "&input\n  fildyn = 'dyn'\n  zasr   = 'simple'\n"
            "  flfrc  = 'fuerzas.fc'\n/\n"
        )
        run.kpath = kpoints.get_kpath(atoms)
        qlines = []
        path = run.kpath.path
        entries = []
        for i, (a, b) in enumerate(path):
            entries.append((a, band_points))
            last = i == len(path) - 1
            disc = (not last) and path[i + 1][0] != b
            if last:
                entries.append((b, 1))
            elif disc:
                entries.append((b, 0))
        run.band_entries = list(entries)
        qlines.append(str(len(entries)))
        for lab, npts in entries:
            x, y, z = run.kpath.point_coords[lab]
            qlines.append(f"  {x:12.8f} {y:12.8f} {z:12.8f} {npts:4d}")
        sweep.write_input(out / "matdyn_band.in", 
            "&input\n  flfrc  = 'fuerzas.fc'\n  asr    = 'simple'\n"
            "  flfrq  = 'bandas.freq'\n  q_in_band_form   = .true.\n"
            "  q_in_cryst_coord = .true.\n/\n" + "\n".join(qlines) + "\n"
        )
        sweep.write_input(out / "matdyn_dos.in", 
            "&input\n  flfrc  = 'fuerzas.fc'\n  asr    = 'simple'\n"
            "  dos    = .true.\n  fldos  = 'fonones.dos'\n"
            f"  nk1 = {dos_grid[0]}\n  nk2 = {dos_grid[1]}\n"
            f"  nk3 = {dos_grid[2]}\n/\n"
        )

    report = ["--- Fonones (DFPT) ---",
              f"Estructura: {atoms.get_chemical_formula()} "
              f"({len(atoms)} átomos, celda primitiva estandarizada)",
              f"Malla k del scf: {grid_scf[0]}x{grid_scf[1]}x{grid_scf[2]}  |  "
              "conv_thr = 1e-12 (la DFPT lo necesita)"]
    if gamma_only:
        report.append("Modo Γ: ph.x en un solo q + dynmat.x "
                      "(frecuencias y actividades IR)")
    else:
        report.append(f"Malla de q: {qgrid[0]}x{qgrid[1]}x{qgrid[2]}  |  "
                      f"dispersión por el camino de alta simetría + DOS "
                      f"{dos_grid[0]}x{dos_grid[1]}x{dos_grid[2]}")
    if epsil:
        report.append("epsil = .true.: tensor dieléctrico y cargas de Born "
                      "(separación LO–TO)")
    report += ["", "IMPORTANTE: la estructura debe estar relajada (vc-relax) "
                    "con estos mismos\ncutoffs; si no, aparecerán frecuencias "
                    "imaginarias espurias.",
               f"Archivos escritos en '{out.resolve()}'"]
    warn = sweep.missing_pseudo_warning(common)
    if warn:
        report.append(warn)
    return run, "\n".join(report)


# ----------------------------------------------------------------------
# Ejecución de la cadena
# ----------------------------------------------------------------------
def run_chain(run: PhononRun, pw_cmd: str = None, nproc: int = None,
              verbose: bool = True):
    """ph.x -> (dynmat | q2r + matdyn x2), asumiendo el scf ya corrido."""
    import time
    steps = [("ph", "ph.x")]
    if run.gamma_only:
        steps.append(("dynmat", "dynmat.x"))
    else:
        steps += [("q2r", "q2r.x"), ("matdyn_band", "matdyn.x"),
                  ("matdyn_dos", "matdyn.x")]
    for stem, exe in steps:
        # Reanudable, con la misma filosofía del runner: un paso cuyo .out
        # ya dice JOB DONE no se repite. Importa sobre todo para ph.x, que
        # es el paso caro y el que un usuario corre a mano en el clúster.
        if _step_done(run.outdir, stem):
            if verbose:
                print(f"  {exe} ({stem}) ... ya estaba hecho")
            continue
        cmd = _exe(exe, pw_cmd, nproc)
        if verbose:
            print(f"  {exe} ({stem}) ... ", end="", flush=True)
        t0 = time.time()
        _run_tool(cmd, run.outdir, stem)
        if verbose:
            print(f"{time.time() - t0:.1f} s")


# ----------------------------------------------------------------------
# Lectura de resultados
# ----------------------------------------------------------------------
def _read_flfrq(path: Path):
    """Lee el archivo de frecuencias de matdyn (formato &plot)."""
    text = path.read_text().split()
    if not text or text[0] != "&plot":
        raise FaltanDatos(
            f"'{path}' no tiene el formato de frecuencias de matdyn.x "
            "(debe empezar con '&plot'); ¿corrió matdyn.x correctamente?"
        )
    nbnd = int(text[2].rstrip(","))
    nks = int(text[4])
    vals = [float(x) for x in text[6:] if x not in ("/",)]
    qs, freqs = [], []
    i = 0
    for _ in range(nks):
        qs.append(vals[i:i + 3]); i += 3
        freqs.append(vals[i:i + nbnd]); i += nbnd
    return np.array(qs), np.array(freqs)


_RE_TABLA = re.compile(
    r"^\s*(\d+)\s+(-?[\d.]+)\s+(-?[\d.]+)\s+([\d.eE+-]+)"
    r"(?:\s+([\d.eE+-]+)\s+([\d.eE+-]+))?\s*$")


def read_dynmat_table(path) -> list:
    """Lee la tabla '# mode [cm-1] [THz] IR [Raman depol]' de dynmat.x.

    Las columnas de Raman solo aparecen si ph.x corrió con lraman; el
    parser acepta las dos formas para no romperse con un cálculo normal.
    """
    modos = []
    en_tabla = False
    for linea in Path(path).read_text(errors="ignore").splitlines():
        if linea.lstrip().startswith("# mode"):
            en_tabla = True
            continue
        if not en_tabla:
            continue
        m = _RE_TABLA.match(linea)
        if not m:
            if modos:
                break
            continue
        g = m.groups()
        d = {"modo": int(g[0]), "omega_cm1": float(g[1]),
             "omega_thz": float(g[2]), "ir": float(g[3])}
        if g[4] is not None:
            d["raman"] = float(g[4])
            d["depol"] = float(g[5])
        modos.append(d)
    return modos


def collect(run: PhononRun) -> PhononRun:
    out = run.outdir
    if run.gamma_only:
        tabla = read_dynmat_table(out / "dynmat.out")
        if tabla:
            run.modes = tabla
            run.gamma_freqs = [(d["omega_cm1"], d.get("ir")) for d in tabla]
            return run
        text = (out / "dynmat.out").read_text(errors="ignore")
        freqs = []
        for line in text.splitlines():
            # "  1      0.00    [cm-1]   --- ..."   formato de dynmat
            parts = line.split()
            if len(parts) >= 4 and parts[0].isdigit() and "[cm-1]" in line:
                try:
                    w = float(parts[2] if parts[1] == "-" else parts[1])
                except ValueError:
                    continue
                freqs.append((w, None))
        # actividades IR (si dynmat las imprimió)
        if "IR activities" in text or "IR cross" in text:
            freqs2 = []
            grab = False
            for line in text.splitlines():
                if "mode" in line.lower() and "ir" in line.lower():
                    grab = True; continue
                if grab:
                    parts = line.split()
                    if len(parts) >= 3:
                        try:
                            freqs2.append((float(parts[1]), float(parts[-1])))
                        except ValueError:
                            if freqs2:
                                break
            if freqs2:
                run.gamma_freqs = freqs2
                return run
        run.gamma_freqs = freqs
        return run

    qs, fr = _read_flfrq(out / "bandas.freq")
    run.band_q, run.band_freqs = qs, fr
    kp = run.kpath
    if kp is not None:
        # OJO: matdyn.x escribe flfrq con q en CARTESIANAS, en unidades de
        # 2*pi/alat, aunque el input se le haya dado en cristalinas. alat es
        # |a1| de la celda usada en el scf (la primitiva estandarizada).
        cell = np.asarray(kp.primitive.cell.array, dtype=float)
        alat = np.linalg.norm(cell[0])
        qcart = qs * (2.0 * np.pi / alat)          # A^-1

        # los índices de los puntos etiquetados y de las discontinuidades
        # salen de las entradas del camino, no de coincidencias numéricas:
        # U y K, por ejemplo, son el mismo punto salvo un vector de la red
        # recíproca y una comparación de coordenadas los confunde.
        idx, labels, breaks = 0, [], []
        entries = run.band_entries or []
        for lab, npts in entries:
            labels.append((idx, kpoints.pretty_label(lab)))
            if npts == 0:
                breaks.append(idx)          # salto: no acumular distancia
            # peso n = n puntos desde este hasta el siguiente; el siguiente
            # punto etiquetado cae en idx+n. Peso 0 (salto) avanza uno solo.
            idx += max(npts, 1)

        d = np.zeros(len(qcart))
        for i in range(1, len(qcart)):
            if (i - 1) in breaks:
                d[i] = d[i - 1]             # discontinuidad: mismo x
            else:
                d[i] = d[i - 1] + np.linalg.norm(qcart[i] - qcart[i - 1])
        run.qdist = d

        # fusionar las etiquetas de una discontinuidad: "U|K"
        merged, skip = [], set()
        for j, (i, lab) in enumerate(labels):
            if j in skip:
                continue
            if j + 1 < len(labels) and labels[j + 1][0] == i + 1 and i in breaks:
                merged.append((i, f"{lab}|{labels[j + 1][1]}"))
                skip.add(j + 1)
            else:
                merged.append((i, lab))
        run.labels = [(i, lab) for i, lab in merged if i < len(qcart)]

    dos_data = np.loadtxt(out / "fonones.dos", comments="#")
    run.dos_w, run.dos = dos_data[:, 0], dos_data[:, 1]
    return run


# ----------------------------------------------------------------------
# Termodinámica armónica desde la DOS
# ----------------------------------------------------------------------

def raman_spectrum(run: PhononRun, laser_nm: float = 532.0,
                   T: float = 300.0, fwhm: float = 5.0,
                   wmin: float = 0.0, wmax: float = None,
                   npts: int = 2000):
    """Espectro Raman medible a partir de las ACTIVIDADES calculadas.

    La actividad que da dynmat.x (A^4/amu) no es lo que mide un
    espectrometro. La intensidad Stokes va como

        I(w) ~ (wL - w)^4 / w * [n(w) + 1] * A(w)

    con n(w) = 1/(exp(hbar*w/kT) - 1) el factor de Bose-Einstein y wL la
    frecuencia del laser. Los tres factores importan: el (wL-w)^4 pesa
    mucho mas los modos de baja frecuencia, y el de Bose depende de la
    temperatura. Comparar actividades crudas contra un espectro medido es
    el error tipico.

    Devuelve (w en cm-1, intensidad normalizada a 100).
    """
    if not run.modes:
        raise FaltanDatos("no hay modos con actividad Raman; corre con "
                         "--raman")
    act = [(d["omega_cm1"], d.get("raman")) for d in run.modes
           if d.get("raman") is not None]
    if not act:
        raise FaltanDatos("el calculo no trae columna Raman: ph.x tiene que "
                         "haber corrido con lraman = .true.")

    w_laser = 1.0e7 / float(laser_nm)          # nm -> cm-1
    picos = [(w, a) for w, a in act if w > 1.0]     # fuera los acusticos
    if wmax is None:
        wmax = max(w for w, _ in picos) * 1.15 if picos else 100.0
    grid = np.linspace(wmin, wmax, npts)
    total = np.zeros_like(grid)
    intens = []
    for w, a in picos:
        x = CM1_TO_EV * w / (KB_EV * T)
        nb = 1.0 / np.expm1(x) if x > 1e-8 else 0.0
        I = ((w_laser - w) ** 4) / w * (nb + 1.0) * a
        intens.append((w, I))
        g = (fwhm / 2.0) ** 2 / ((grid - w) ** 2 + (fwhm / 2.0) ** 2)
        total += I * g
    if total.max() > 0:
        total = total * (100.0 / total.max())
    return grid, total, intens


def thermodynamics(run: PhononRun, T=None, natoms: int = None):
    """ZPE, F, S, U y C_v por celda a partir de la DOS de fonones.

    La DOS de matdyn está en estados/cm⁻¹ y su integral vale 3N; se
    renormaliza a 3N exacto para absorber el error de malla.
    """
    if run.dos is None:
        raise FaltanDatos("no hay DOS de fonones; corre primero la cadena")
    if T is None:
        T = np.arange(0.0, 1001.0, 10.0)
    T = np.asarray(T, dtype=float)
    w = run.dos_w
    g = np.array(run.dos, dtype=float)
    mask = w > 1.0            # descartar el entorno de w=0 (acústicos)
    w, g = w[mask], g[mask]
    e = w * CM1_TO_EV         # eV por modo

    nmodes = trapezoid(g, w)
    if natoms:
        g = g * (3.0 * natoms / nmodes)

    zpe = trapezoid(0.5 * e * g, w)
    F, S, Cv, U = [], [], [], []
    for t in T:
        if t < 1e-6:
            F.append(zpe); U.append(zpe); S.append(0.0); Cv.append(0.0)
            continue
        x = e / (KB_EV * t)
        x = np.minimum(x, 500.0)
        occ = 1.0 / np.expm1(x)
        F.append(zpe + KB_EV * t * trapezoid(np.log1p(-np.exp(-x)) * g, w))
        U.append(trapezoid((0.5 + occ) * e * g, w))
        Cv.append(KB_EV * trapezoid(x ** 2 * np.exp(x) * occ ** 2 * g, w))
        S.append((U[-1] - F[-1]) / t)
    return {"T": T, "ZPE": float(zpe), "F": np.array(F), "U": np.array(U),
            "S": np.array(S), "Cv": np.array(Cv)}


# ----------------------------------------------------------------------
# Reporte, exportación, gráficas
# ----------------------------------------------------------------------
def report_gamma_activities(run: PhononRun) -> str:
    """Tabla de modos en Gamma con actividad IR y Raman."""
    if not run.modes:
        return ""
    tiene_raman = any("raman" in d for d in run.modes)
    lines = ["--- Modos en Gamma ---",
             f"{'modo':>5s} {'cm-1':>10s} {'THz':>9s} {'IR':>12s}"
             + (f" {'Raman':>12s} {'depol':>7s}" if tiene_raman else "")]
    for d in run.modes:
        fila = (f"{d['modo']:5d} {d['omega_cm1']:10.2f} "
                f"{d['omega_thz']:9.4f} {d['ir']:12.4f}")
        if tiene_raman:
            fila += f" {d.get('raman', float('nan')):12.4f} " \
                    f"{d.get('depol', float('nan')):7.4f}"
        lines.append(fila)
    lines.append("")
    lines.append("IR en (D/A)^2/amu; Raman en A^4/amu (unidades de QE).")

    opticos = [d for d in run.modes if d["omega_cm1"] > 1.0]
    if tiene_raman and opticos:
        ir_act = [d for d in opticos if d["ir"] > 1e-4]
        ra_act = [d for d in opticos if d.get("raman", 0.0) > 1e-4]
        if ra_act and not ir_act:
            lines.append(
                "\nTodos los modos ópticos son activos en Raman e INACTIVOS "
                "en IR: es la regla\nde exclusión mutua, y confirma que el "
                "cristal tiene centro de inversión.")
        elif ir_act and not ra_act:
            lines.append(
                "\nLos modos ópticos son activos en IR e inactivos en "
                "Raman.")
        depols = {round(d.get("depol", -1), 3) for d in ra_act}
        if depols == {0.75}:
            lines.append(
                "Factor de despolarización 0.75 en todos los modos activos: "
                "el valor exacto\nde un modo triplemente degenerado.")
    return "\n".join(lines)


def report(run: PhononRun, natoms: int = None) -> str:
    lines = ["--- Fonones ---"]
    if run.gamma_only:
        lines.append("Frecuencias en Γ (dynmat.x, con regla de suma acústica):")
        has_ir = any(a is not None for _, a in run.gamma_freqs)
        head = f"  {'modo':>4s} {'ω (cm⁻¹)':>12s} {'ω (THz)':>10s}"
        if has_ir:
            head += f" {'IR (u. arb.)':>14s}"
        lines.append(head)
        for i, (w, act) in enumerate(run.gamma_freqs, start=1):
            row = f"  {i:>4d} {w:12.2f} {w * CM1_TO_THZ:10.3f}"
            if has_ir and act is not None:
                row += f" {act:14.4f}"
            lines.append(row)
        neg = [w for w, _ in run.gamma_freqs if w < -5.0]
        if neg:
            lines.append("\nAVISO: hay frecuencias imaginarias (negativas). "
                         "O la estructura no está\nrelajada, o es inestable en Γ.")
        return "\n".join(lines)

    fr = run.band_freqs
    lines.append(f"Dispersión: {fr.shape[1]} ramas en {fr.shape[0]} puntos q")
    wmin, wmax = float(fr.min()), float(fr.max())
    lines.append(f"Rango: {wmin:.1f} a {wmax:.1f} cm⁻¹ "
                 f"({wmax * CM1_TO_THZ:.2f} THz máx.)")
    if wmin < -5.0:
        n_neg = int((fr < -5.0).sum())
        lines.append(f"\nAVISO: {n_neg} frecuencias imaginarias (mínimo "
                     f"{wmin:.1f} cm⁻¹). La estructura\nno está relajada con "
                     "estos parámetros, o es dinámicamente inestable.")
    else:
        lines.append("Sin frecuencias imaginarias: la estructura es "
                     "dinámicamente estable.")
    if run.dos is not None and natoms:
        th = thermodynamics(run, natoms=natoms)
        i300 = int(np.searchsorted(th["T"], 300.0))
        lines += ["", "Termodinámica armónica (por celda):",
                  f"  energía de punto cero: {th['ZPE'] * 1000:.2f} meV",
                  f"  a 300 K:  C_v = {th['Cv'][i300] * 1000:.3f} meV/K   "
                  f"S = {th['S'][i300] * 1000:.3f} meV/K   "
                  f"F = {th['F'][i300] * 1000:.1f} meV"]
    return "\n".join(lines)


def export(run: PhononRun, outdir: str = ".", natoms: int = None) -> list:
    out = Path(outdir); out.mkdir(parents=True, exist_ok=True)
    written = []
    if run.gamma_only:
        f = out / "FONONES_GAMMA.dat"
        lines = [provenance.header("fonones en Gamma",
                                   {"epsil": run.epsil}),
                 "# modo  omega(cm-1)  omega(THz)  IR"]
        for i, (w, a) in enumerate(run.gamma_freqs, start=1):
            lines.append(f"{i:6d} {w:12.3f} {w * CM1_TO_THZ:11.4f} "
                         f"{a if a is not None else float('nan'):12.4f}")
        f.write_text("\n".join(lines) + "\n")
        return [str(f)]

    f = out / "FONONES_BANDAS.dat"
    header = (provenance.header("dispersión de fonones",
                                {"malla_q": "x".join(map(str, run.qgrid))
                                 if run.qgrid else None},
                                titulo="Dispersion fononica")
              + "\n# q(acum, A^-1) + ramas en cm^-1")
    np.savetxt(f, np.column_stack([run.qdist, run.band_freqs]),
               fmt="%12.5f", header=header, comments="")
    written.append(str(f))
    f = out / "FONONES_DOS.dat"
    np.savetxt(f, np.column_stack([run.dos_w, run.dos]), fmt="%14.6f",
               header=provenance.header("DOS de fonones")
                      + "\n# omega(cm-1)  dos(estados/cm-1)", comments="")
    written.append(str(f))
    if natoms:
        th = thermodynamics(run, natoms=natoms)
        f = out / "FONONES_TERMO.dat"
        np.savetxt(f, np.column_stack([th["T"], th["F"], th["U"], th["S"],
                                       th["Cv"]]),
                   fmt="%14.6e",
                   header=provenance.header(
                       "termodinámica armónica",
                       {"atomos_por_celda": natoms}) + "\n"
                          f"# ZPE = {th['ZPE']:.6e} eV/celda\n"
                          "# T(K)  F(eV)  U(eV)  S(eV/K)  Cv(eV/K)",
                   comments="")
        written.append(str(f))
    return written


def has_dispersion(run: PhononRun) -> bool:
    """¿Tiene esta corrida una dispersión que se pueda dibujar?

    Los cálculos en Γ (--gamma, y también --raman, que fuerza gamma_only)
    solo tienen frecuencias en un punto q: no hay band_freqs ni qdist, y
    plot() fallaría. El despachador de la CLI decide con esto, no con la
    bandera --gamma, para no equivocarse cuando Γ se activó por otra vía.
    """
    return (not run.gamma_only and run.band_freqs is not None
            and run.qdist is not None)


def plot(run: PhononRun, outfile: str = "fonones", formats="pdf,png",
         theme: str = None, size: str = None, family: str = None,
         background: str = None, palette=None, usetex: bool = None,
         width="double", journal: str = "generic", aspect: float = 0.45,
         mono: bool = False, dpi: int = None) -> list:
    """Dispersión + DOS de fonones con eje de frecuencia compartido."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    st = qstyle.apply(theme, size=size, family=family, background=background,
                      palette=palette, usetex=usetex, mono=mono)
    figsize = qstyle.figure_size(width, journal, aspect)
    fig = plt.figure(figsize=figsize, layout="constrained")
    gs = fig.add_gridspec(1, 2, width_ratios=[2.6, 1.0], wspace=0.08)
    axb = qstyle.finish_axes(fig.add_subplot(gs[0]))
    axd = qstyle.finish_axes(fig.add_subplot(gs[1], sharey=axb))
    color = qstyle.palette(2, mono=mono)

    for b in range(run.band_freqs.shape[1]):
        axb.plot(run.qdist, run.band_freqs[:, b], color=qstyle.INK,
                 lw=st["line"])
    axb.axhline(0.0, color=qstyle.INK_FAINT, lw=st["axis_line"],
                dashes=[3.5, 2.0])
    if run.labels:
        ticks = [run.qdist[i] for i, _ in run.labels]
        axb.set_xticks(ticks)
        axb.set_xticklabels([qstyle.tex_safe(l) for _, l in run.labels])
        axb.xaxis.set_minor_locator(plt.NullLocator())
        for t in ticks[1:-1]:
            axb.axvline(t, color=qstyle.GRID, lw=st["axis_line"], zorder=0)
    axb.set_xlim(run.qdist[0], run.qdist[-1])
    axb.set_ylabel(r"$\omega$ (cm$^{-1}$)")
    qstyle.panel_label(axb, "(a)")

    axd.plot(run.dos, run.dos_w, color=color[0], lw=st["line"])
    axd.fill_betweenx(run.dos_w, 0, run.dos, color=color[0], alpha=0.12, lw=0)
    # la escala absoluta de la DOS de matdyn depende de la malla; lo que
    # se lee en la figura es la forma, así que se rotula sin números y el
    # valor numérico queda en FONONES_DOS.dat
    axd.set_xlabel(qstyle.tex_safe("DOS (estados/cm-1)").replace(
        "cm-1", "cm$^{-1}$"))
    axd.set_xticks([])
    axd.xaxis.set_minor_locator(plt.NullLocator())
    axd.tick_params(labelleft=False)
    axd.set_xlim(0, run.dos.max() * 1.1)
    qstyle.panel_label(axd, "(b)")

    ymin = min(-10.0, float(run.band_freqs.min()) * 1.05)
    ymax = float(run.band_freqs.max()) * 1.06
    axb.set_ylim(ymin, ymax)
    written = qstyle.save(fig, outfile, formats, dpi=dpi,
                          modulo="fonones (DFPT)")
    plt.close(fig)
    return written
