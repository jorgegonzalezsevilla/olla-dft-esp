# Olla-DFT — command-line toolkit for Quantum ESPRESSO
# Copyright (C) 2026 Jorge Enrique González Sevilla
# SPDX-License-Identifier: GPL-3.0-or-later
# This program is free software: you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by the Free
# Software Foundation, either version 3 of the License, or (at your option)
# any later version. See the LICENSE file for details.

"""Diagnóstico de un cálculo de pw.x: ¿sirve, y si no, por qué?

Todo lo que este módulo usa ya estaba en los archivos que QE deja: el XML
guarda si convergió, en cuántos pasos y con qué error; el stdout guarda la
historia completa de iteraciones SCF y, en un relax, cada paso iónico.
Nadie los mira, y son justo los que dicen si el resultado sirve.

LO QUE DISTINGUE ESTE MÓDULO
----------------------------
Un SCF que no converge tiene al menos dos causas con remedios OPUESTOS:

- **oscilación de carga** (charge sloshing): el error sube y baja en vez de
  bajar. Típico de losas grandes, metales y celdas con vacío. El remedio es
  mezclar MENOS (bajar mixing_beta) y usar mixing_mode='local-TF';
- **convergencia lenta monótona**: el error baja siempre, pero demasiado
  despacio. Ahí el remedio es el contrario, mezclar MÁS (subir
  mixing_beta) o simplemente dar más pasos.

Aplicar el remedio equivocado empeora el problema, así que el módulo
distingue los dos casos por la forma de la curva en vez de dar un consejo
genérico.
"""

import re
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

from qekit.core import qeout
from qekit.core import style as qstyle

_RE_ITER = re.compile(r"iteration #\s*(\d+)\s+ecut=\s*([\d.]+)\s*Ry\s+beta=\s*([\d.]+)")
_RE_ACC = re.compile(r"estimated scf accuracy\s*<\s*([\dEe.+-]+)\s*Ry")
_RE_ETOT = re.compile(r"total energy\s*=\s*(-?[\dEe.+-]+)\s*Ry")
_RE_FORCE = re.compile(r"Total force\s*=\s*([\dEe.+-]+)")
_RE_PRESS = re.compile(r"P=\s*(-?[\dEe.+-]+)")
_RE_WARN = re.compile(r"^\s*(Warning|WARNING|%%%%)", re.M)
_RE_MAXSTEP = re.compile(r"convergence NOT achieved after\s*(\d+)\s*iterations")


@dataclass
class ScfHistory:
    accuracy: list = field(default_factory=list)     # Ry, por iteración (último ciclo)
    energies: list = field(default_factory=list)     # Ry (último ciclo)
    beta: float = None
    converged: bool = None
    n_iter: int = 0
    n_ciclos: int = 1            # ciclos SCF vistos (1 en un scf, N en un relax)
    patologia: str = ""          # "oscilacion" | "lenta" | "estancada" | ""
    consejo: str = ""


@dataclass
class Trajectory:
    energies: list = field(default_factory=list)     # Ry por paso iónico
    forces: list = field(default_factory=list)       # Ry/bohr (Total force)
    pressures: list = field(default_factory=list)    # kbar
    n_steps: int = 0


@dataclass
class Diagnosis:
    result: qeout.QEResult = None
    scf: ScfHistory = None
    traj: Trajectory = None
    warnings: list = field(default_factory=list)
    problemas: list = field(default_factory=list)
    stdout_path: str = ""


def _ciclos_scf(texto: str) -> list:
    """Trozos del stdout, uno por ciclo SCF.

    En un relax cada paso iónico arranca su propio ciclo con
    'iteration #  1'; se corta ahí. Un scf da un solo trozo.
    """
    cortes = [m.start() for m in _RE_ITER.finditer(texto)
              if int(m.group(1)) == 1]
    if not cortes:
        return [texto]
    cortes.append(len(texto))
    return [texto[a:b] for a, b in zip(cortes, cortes[1:])]


def read_scf_history(stdout_path) -> ScfHistory:
    """Historia de iteraciones SCF del stdout de pw.x, y su diagnóstico.

    En un relax hay un ciclo SCF por paso iónico y son independientes:
    concatenar sus errores mezclaría el final convergido de uno con el
    arranque del siguiente y parecería una oscilación. Se clasifica solo
    el ÚLTIMO ciclo, que es el que decide si el cálculo terminó bien, y
    `n_ciclos` dice cuántos se vieron.
    """
    texto = Path(stdout_path).read_text(errors="ignore")
    h = ScfHistory()
    ciclos = _ciclos_scf(texto)
    ultimo = ciclos[-1]
    h.n_ciclos = len(ciclos)
    h.accuracy = [float(x) for x in _RE_ACC.findall(ultimo)]
    h.energies = [float(x) for x in _RE_ETOT.findall(ultimo)]
    betas = _RE_ITER.findall(texto)
    if betas:
        h.beta = float(betas[0][2])
    h.n_iter = len(h.accuracy)
    h.converged = "convergence has been achieved" in ultimo
    if _RE_MAXSTEP.search(ultimo):
        h.converged = False
    _clasificar(h)
    return h


def _clasificar(h: ScfHistory) -> None:
    """Distingue oscilacion de convergencia lenta por la FORMA de la curva.

    Se miran DOS cosas, no una: con que frecuencia sube el error y CUANTO
    sube. Unas pocas subidas enormes son mejor senal de oscilacion que
    muchas subidas minimas, asi que contar solo la frecuencia se equivoca.

    Las dos primeras iteraciones se ignoran: un salto grande al arranque es
    un transitorio normal mientras la densidad inicial se acomoda, no una
    patologia.
    """
    a = np.array(h.accuracy, dtype=float)
    if a.size < 2 or h.converged:
        return

    beta = h.beta if h.beta is not None else 0.4

    # Con muy pocas iteraciones no hay forma honesta de distinguir
    # oscilacion de lentitud: una sola subida entre tres diferencias ya da
    # 33 %, y eso es ruido, no diagnostico. Mejor decirlo que inventarlo.
    if a.size < 8:
        h.patologia = "pocos_datos"
        h.consejo = (
            f"solo {a.size} iteraciones: no alcanzan para distinguir "
            "oscilacion de carga de\nconvergencia lenta, que piden remedios "
            "opuestos. Sube electron_maxstep\n(a 100 o mas) y vuelve a "
            "mirar la curva; con el ciclo cortado tan pronto,\ncualquier "
            "diagnostico seria adivinar.")
        return

    cola = a[2:] if a.size > 5 else a          # saltar el transitorio
    difs = np.diff(cola)
    subidas = int(np.sum(difs > 0))
    frac_subidas = subidas / max(len(cola) - 1, 1)
    with np.errstate(divide="ignore", invalid="ignore"):
        razones = cola[1:] / np.maximum(cola[:-1], 1e-300)
        peor = float(np.max(razones)) if razones.size else 1.0
        decadas = np.log10(a[0] / a[-1]) if a[-1] > 0 else np.inf

    # el criterio de frecuencia solo se aplica con suficientes puntos
    frecuente = len(cola) >= 6 and frac_subidas > 0.25
    if frecuente or peor > 5.0:
        motivo = (f"el error subio en {subidas} de {len(cola)-1} "
                  "iteraciones" if frecuente
                  else f"el error se multiplico por {peor:.0f} en una "
                       "iteracion")
        h.patologia = "oscilacion"
        h.consejo = (
            f"{motivo}: es OSCILACION DE CARGA, tipica de losas,\n"
            "metales y celdas con mucho vacio. Se arregla mezclando MENOS:\n"
            f"  mixing_beta = {max(0.05, beta / 3):.2f}   "
            f"(ahora {beta:.2f})\n"
            "  mixing_mode = 'local-TF'   (pensado justo para este caso)\n"
            "  mixing_ndim = 12           (mas historia de mezcla)\n"
            "Subir mixing_beta aqui lo EMPEORA.")
    elif decadas < 3:
        h.patologia = "estancada"
        h.consejo = (
            f"en {len(a)} iteraciones el error solo bajo {decadas:.1f} "
            "ordenes de magnitud y se\nquedo plano: no es lentitud, esta "
            "estancado. Suele ser un estado magnetico o\nde ocupaciones mal "
            "planteado, o una estructura con atomos casi encima.\nRevisa "
            "starting_magnetization, el smearing y las distancias "
            "interatomicas.")
    else:
        h.patologia = "lenta"
        # Si beta ya es agresivo, subirlo mas seria justo lo contrario de
        # lo que conviene: ahi lo que falta son pasos, no mezcla.
        if beta >= 0.6:
            h.consejo = (
                f"el error baja de forma monotona y llego a {a[-1]:.1e} Ry "
                "sin alcanzar el\numbral: le faltaron PASOS, no mezcla. "
                f"mixing_beta ya esta en {beta:.2f}, que es\nagresivo; "
                "subirlo mas arriesga desestabilizarlo.\n"
                "  electron_maxstep = 300")
        else:
            h.consejo = (
                "el error baja de forma monotona pero no llego al umbral: "
                "es convergencia\nLENTA, no oscilacion. Aqui si conviene "
                "mezclar MAS o dar mas pasos:\n"
                f"  mixing_beta = {min(0.7, max(beta * 1.75, 0.3)):.2f}   "
                f"(ahora {beta:.2f})\n"
                "  electron_maxstep = 300")


def read_trajectory(stdout_path) -> Trajectory:
    """Pasos iónicos de un relax/vc-relax desde el stdout."""
    texto = Path(stdout_path).read_text(errors="ignore")
    t = Trajectory()
    # las energías de cada paso convergido llevan '!' delante
    t.energies = [float(m) for m in re.findall(
        r"^!\s+total energy\s*=\s*(-?[\dEe.+-]+)\s*Ry", texto, re.M)]
    t.forces = [float(x) for x in _RE_FORCE.findall(texto)]
    t.pressures = [float(x) for x in _RE_PRESS.findall(texto)]
    t.n_steps = len(t.energies)
    return t


def find_stdout(workdir) -> Path:
    """Busca el stdout de pw.x en una carpeta de cálculo."""
    d = Path(workdir)
    if d.is_file():
        return d
    for patron in ("*.out", "pw.out", "scf.out", "relax.out", "out.*"):
        for f in sorted(d.glob(patron)):
            try:
                head = f.read_text(errors="ignore")[:4000]
            except OSError:
                continue
            if "Program PWSCF" in head:
                return f
    return None


def diagnose(workdir, prefix: str = None) -> Diagnosis:
    """Lee XML y stdout de una carpeta y arma el diagnóstico."""
    d = Diagnosis()
    try:
        xml = qeout.find_xml(str(workdir), prefix)
        d.result = qeout.read_xml(xml)
    except (FileNotFoundError, ValueError) as exc:
        d.problemas.append(f"no se pudo leer el XML: {exc}")

    so = find_stdout(workdir)
    if so is not None:
        d.stdout_path = str(so)
        d.scf = read_scf_history(so)
        d.traj = read_trajectory(so)
        texto = so.read_text(errors="ignore")
        for linea in texto.splitlines():
            if _RE_WARN.match(linea) and linea.strip() not in d.warnings:
                d.warnings.append(linea.strip())
        if "Error in routine" in texto:
            i = texto.index("Error in routine")
            d.problemas.append(texto[i:i + 200].split("%%%")[0].strip())

    r = d.result
    if r is not None:
        if r.converged is False:
            d.problemas.append("el SCF NO convergió: el resultado no sirve")
        if r.max_force is not None and r.max_force > 0.05:
            d.problemas.append(
                f"fuerza residual máxima {r.max_force:.4f} eV/Å: la "
                "estructura no está relajada (umbral usual 0.01–0.03)")
        if r.pressure is not None and abs(r.pressure) > 1.0 and \
                r.calculation in ("vc-relax", "relax", "scf"):
            d.problemas.append(
                f"presión residual {r.pressure:+.2f} GPa: la celda no está "
                "en equilibrio con estos cutoffs")
    return d


def report(d: Diagnosis) -> str:
    r = d.result
    lines = ["--- Diagnóstico del cálculo ---"]
    if r is not None:
        lines += [f"Archivo: {r.xml_path}",
                  f"Tipo: {r.calculation or '?'}  |  "
                  f"{r.functional or '?'}  |  ecut {r.ecutwfc or '?'}/"
                  f"{r.ecutrho or '?'} Ry"]
        if r.kgrid:
            lines.append(f"Malla k: {r.kgrid[0]}x{r.kgrid[1]}x{r.kgrid[2]}"
                         f"  |  {r.nk} puntos en la ZBI  |  "
                         f"{r.n_sym or '?'} operaciones de simetría")
        estado = ("convergió" if r.converged else "NO convergió"
                  if r.converged is not None else "convergencia desconocida")
        extra = ""
        if r.n_scf_steps:
            extra = f" en {r.n_scf_steps} pasos"
        if r.scf_error is not None:
            extra += f", error final {r.scf_error:.2e} Ry"
        lines.append(f"SCF: {estado}{extra}")
        if r.max_force is not None:
            lines.append(f"Fuerza residual máxima: {r.max_force:.5f} eV/Å")
        if r.pressure is not None:
            lines.append(f"Presión residual: {r.pressure:+.3f} GPa")
        if r.total_magnetization is not None and \
                abs(r.total_magnetization) > 1e-8:
            lines.append(f"Magnetización: {r.total_magnetization:.3f} μB "
                         f"(absoluta {r.absolute_magnetization:.3f})")
        if r.wall_time:
            lines.append(f"Tiempo: {r.wall_time:.1f} s de reloj "
                         f"({r.cpu_time:.1f} s de CPU)")

    if d.traj and d.traj.n_steps > 1:
        t = d.traj
        lines += ["", f"Relajación: {t.n_steps} pasos iónicos",
                  f"  energía: {t.energies[0]:.6f} -> {t.energies[-1]:.6f} Ry "
                  f"({(t.energies[-1]-t.energies[0])*13.6057:.4f} eV)"]
        if t.forces:
            lines.append(f"  fuerza total: {t.forces[0]:.5f} -> "
                         f"{t.forces[-1]:.5f} Ry/bohr")
        subidas = sum(1 for a, b in zip(t.energies, t.energies[1:]) if b > a)
        if subidas > t.n_steps // 3:
            lines.append(
                f"  AVISO: la energía subió en {subidas} de {t.n_steps-1} "
                "pasos. Una relajación sana\n  baja casi siempre; esto "
                "sugiere una superficie de energía muy plana o un\n  paso "
                "de BFGS demasiado grande.")

    if d.scf and d.scf.n_iter:
        ciclos = ""
        if d.scf.n_ciclos > 1:
            ciclos = (f" en el último de {d.scf.n_ciclos} ciclos SCF (uno "
                      "por paso iónico; se diagnostica solo el último)")
        lines += ["", f"Historia SCF: {d.scf.n_iter} iteraciones{ciclos}"
                       f"{'' if d.scf.beta is None else f', beta = {d.scf.beta:.2f}'}"]
        if d.scf.patologia:
            lines += ["", f"PROBLEMA DE CONVERGENCIA ({d.scf.patologia}):",
                      d.scf.consejo]

    if d.problemas:
        lines += ["", "PROBLEMAS:"]
        lines += [f"  - {p}" for p in d.problemas]
    elif r is not None and r.converged:
        lines += ["", "Sin problemas detectados."]

    if d.warnings:
        lines += ["", f"Avisos de QE ({len(d.warnings)}):"]
        lines += [f"  {w}" for w in d.warnings[:6]]
    return "\n".join(lines)


def plot(d: Diagnosis, outfile: str = "diagnostico", formats="pdf,png",
         theme: str = None, family: str = None, background: str = None,
         palette=None, usetex: bool = None, width="double",
         journal: str = "generic", aspect: float = 0.42,
         mono: bool = False, dpi: int = None) -> list:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    tiene_traj = bool(d.traj and d.traj.n_steps > 1)
    n = 2 if tiene_traj else 1
    st = qstyle.apply(theme, family=family, background=background,
                      palette=palette, usetex=usetex, mono=mono)
    fig = plt.figure(figsize=qstyle.figure_size(width, journal, aspect),
                     layout="constrained")
    c = qstyle.palette(3, mono=mono)

    ax = qstyle.finish_axes(fig.add_subplot(1, n, 1))
    if d.scf and d.scf.accuracy:
        ax.semilogy(range(1, len(d.scf.accuracy) + 1), d.scf.accuracy,
                    "o-", color=c[0], lw=st["line"], ms=3)
    ax.set_xlabel("iteración SCF")
    ax.set_ylabel("precisión estimada (Ry)")
    qstyle.panel_label(ax, "(a)")

    if tiene_traj:
        ax2 = qstyle.finish_axes(fig.add_subplot(1, n, 2))
        e = np.array(d.traj.energies) * 13.605693
        ax2.plot(range(1, len(e) + 1), e - e[-1], "o-", color=c[1],
                 lw=st["line"], ms=3)
        ax2.set_xlabel("paso iónico")
        ax2.set_ylabel(r"$E - E_\mathrm{final}$ (eV)")
        qstyle.panel_label(ax2, "(b)")

    written = qstyle.save(fig, outfile, formats, dpi=dpi,
                          modulo="diagnóstico")
    plt.close(fig)
    return written
