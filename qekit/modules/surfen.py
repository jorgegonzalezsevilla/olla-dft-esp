# Olla-DFT — command-line toolkit for Quantum ESPRESSO
# Copyright (C) 2026 Jorge Enrique González Sevilla
# SPDX-License-Identifier: AGPL-3.0-or-later
# This program is free software: you can redistribute it and/or modify it
# under the terms of the GNU Affero General Public License as published by the Free
# Software Foundation, either version 3 of the License, or (at your option)
# any later version. See the LICENSE file for details.

"""Energía de superficie: cortar un cristal y ver cuánto cuesta.

    γ = (E_losa − N·E_bulto) / (2A)

La fórmula cabe en una línea y aplicarla directamente casi siempre da mal.
El problema es que E_bulto viene de OTRO cálculo, con otra celda y otra
malla de puntos k, y por bien convergido que esté queda un error residual
ε por átomo. Ese error entra multiplicado por N:

    γ(N) = γ_verdadera + N·ε / (2A)

o sea que γ no converge al engrosar la losa: DERIVA linealmente, y cuanto
más gruesa la haces, peor. Es el error clásico de este cálculo y no avisa de
nada: cada γ(N) por separado parece un número perfectamente razonable.

La salida es el ajuste lineal de Fiorentini–Methfessel: en vez de
importar E_bulto de fuera, se ajusta

    E_losa(N) = 2γA + N·E_bulto

sobre varias losas de grosor distinto. La pendiente da una energía de bulto
CONSISTENTE con los propios cálculos de losa, y la ordenada al origen da
2γA sin ningún error importado. Aquí se hacen los dos y se enseñan juntos:
si la deriva es grande, se ve.

(El método de Boettger, que toma E_bulto de la DIFERENCIA entre losas
consecutivas, es otra variante y NO es la que se hace aquí: solo el ajuste
lineal por mínimos cuadrados sobre todos los grosores.)
"""

from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

from qekit.core import provenance, qeout, structure
from qekit.core.errors import ErrorDeUso, FaltanDatos
from qekit.core import style as qstyle
from qekit.modules import builder, sweep

# 1 eV/Å² = 16.0218 J/m²
EV_A2_A_J_M2 = 16.021766


def reducir_losa(slab):
    """Celda superficial mínima, dejando el vacío intacto.

    El corte (hkl) se hace sobre la celda CONVENCIONAL, porque es la
    referencia de los índices de Miller, y eso suele dar un múltiplo de la
    celda superficial mínima: para Al(111) salen 4 átomos por plano donde
    basta 1. Como γ es por unidad de área, la celda pequeña da exactamente
    el mismo número y cuesta cuatro veces menos.

    La reducción se acepta solo si el eje c no cambia: si spglib decidiera
    reducir también a lo largo del vacío, la losa dejaría de ser una losa.
    """
    try:
        prim = structure.primitive(slab)
    except Exception:                                       # noqa: BLE001
        return slab, 1.0
    c0 = float(np.linalg.norm(slab.cell.array[2]))
    c1 = float(np.linalg.norm(prim.cell.array[2]))
    if abs(c1 - c0) > 1e-6 or len(prim) >= len(slab) or not len(prim):
        return slab, 1.0
    factor = len(slab) / len(prim)
    if abs(factor - round(factor)) > 1e-6:
        return slab, 1.0
    return prim, float(factor)


@dataclass
class GammaRun:
    miller: tuple = (1, 0, 0)
    capas: list = field(default_factory=list)      # nº de capas pedidas
    jobs: list = field(default_factory=list)
    energias: dict = field(default_factory=dict)   # capas -> eV
    natomos: dict = field(default_factory=dict)    # capas -> átomos
    convergido: dict = field(default_factory=dict)
    area: float = None                             # Å², una cara
    E_bulto: float = None                          # eV por átomo, cálculo aparte
    bulto_ok: bool = None
    natomos_bulto: int = 0
    simetrica: bool = True
    polar: bool = False
    relajado: bool = False
    vacio: float = 0.0
    avisos: list = field(default_factory=list)

    # --- resultados del ajuste ---
    gamma_ajuste: float = None      # eV/Å²
    E_bulto_ajuste: float = None    # eV por átomo, de la pendiente
    r2: float = None

    @property
    def caras(self) -> int:
        """2 en una losa simétrica; 1 si solo una cara es la que interesa."""
        return 2

    def gamma_directo(self, capas: int) -> float:
        """γ(N) usando la energía de bulto del cálculo aparte."""
        e = self.energias.get(capas)
        if e is None or self.E_bulto is None or not self.area:
            return None
        n = self.natomos[capas]
        return (e - n * self.E_bulto) / (self.caras * self.area)


def ajustar(run: GammaRun) -> GammaRun:
    """E_losa(N) = 2γA + N·E_bulto, por mínimos cuadrados sobre los átomos."""
    puntos = [(run.natomos[c], run.energias[c]) for c in run.capas
              if run.energias.get(c) is not None]
    if len(puntos) < 2:
        return run
    x = np.array([p[0] for p in puntos], dtype=float)
    y = np.array([p[1] for p in puntos], dtype=float)
    pend, orden = np.polyfit(x, y, 1)
    run.E_bulto_ajuste = float(pend)
    if run.area:
        run.gamma_ajuste = float(orden) / (run.caras * run.area)
    pred = pend * x + orden
    ss_res = float(np.sum((y - pred) ** 2))
    ss_tot = float(np.sum((y - y.mean()) ** 2))
    run.r2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else 1.0
    return run


# ----------------------------------------------------------------------
# Preparación
# ----------------------------------------------------------------------
def prepare(atoms, miller=(1, 0, 0), capas=(3, 4, 5, 6), vacuum: float = 20.0,
            outdir: str = "gamma", fijar: int = 0, relajar: bool = False,
            con_bulto: bool = True, reducir: bool = True,
            pseudo_dir: str = None,
            insulator: bool = False, ecutwfc: float = None,
            ecutrho: float = None, kspacing: float = None,
            vdw: str = None, dipolo: bool = False,
            nspin: int = 1, magnetization: dict = None) -> tuple:
    capas = sorted({int(c) for c in capas})
    if len(capas) < 2:
        raise ErrorDeUso(
            "hacen falta al menos dos grosores para poder ajustar la recta "
            "E(N); con uno solo no hay forma de separar la energía de "
            "superficie de la del bulto. Prueba --layers 3,4,5,6.")
    if min(capas) < 2:
        raise ErrorDeUso(
            f"una losa de {min(capas)} capa(s) no tiene interior: todo lo que "
            "hay es superficie, y el ajuste no significa nada. Empieza en 3.")

    losas, reducciones = {}, {}
    for n in capas:
        info = builder.surface(atoms, miller=miller, layers=n, vacuum=vacuum,
                               fix_layers=fijar)
        if reducir and not fijar:
            # con capas congeladas NO se reduce: las restricciones están
            # puestas sobre átomos concretos y la celda primitiva los
            # reordena, así que se congelarían otros
            chica, factor = reducir_losa(info.atoms)
            if factor > 1.0:
                info.atoms = chica
                reducciones[n] = factor
        losas[n] = info

    ref = losas[capas[0]]
    celda = ref.atoms.cell.array
    area = float(np.linalg.norm(np.cross(celda[0], celda[1])))

    run = GammaRun(miller=tuple(int(m) for m in miller), capas=capas,
                   area=area, simetrica=ref.simetrica, polar=ref.polar,
                   relajado=relajar, vacio=vacuum)

    # Todas las losas comparten pseudos, cutoffs y malla EN EL PLANO. La
    # malla se fija con la losa más pequeña y se reutiliza: si cada grosor
    # eligiera la suya, las energías no serían restables y el ajuste daría
    # una pendiente sin sentido.
    common = sweep.prepare_common(ref.atoms, pseudo_dir, ecutwfc, ecutrho,
                                  insulator)
    out = Path(outdir); out.mkdir(parents=True, exist_ok=True)
    grid = sweep.default_grid(ref.atoms, kspacing)
    calc = "relax" if relajar else "scf"
    extras = dict(vdw=vdw, nspin=nspin, magnetization=magnetization,
                  dipole_correction=3 if dipolo else False)

    for n in capas:
        info = losas[n]
        run.natomos[n] = len(info.atoms)
        job = sweep.write_scf_job(
            info.atoms, common, out / f"capas{n:02d}", f"{n} capas", grid,
            meta={"papel": "losa", "capas": n}, calculation=calc, **extras)
        run.jobs.append(job)
        for w in info.warnings:
            if w not in run.avisos:
                run.avisos.append(w)

    if con_bulto:
        bulto = structure.conventional(atoms)
        run.natomos_bulto = len(bulto)
        # el bulto necesita su propia malla: no tiene vacío, y usar la de la
        # losa (con 1 punto en c) daría una energía sin sentido
        grid_b = sweep.default_grid(bulto, kspacing)
        run.jobs.append(sweep.write_scf_job(
            bulto, common, out / "_bulto", "bulto", grid_b,
            meta={"papel": "bulto"}, calculation="scf",
            vdw=vdw, nspin=nspin, magnetization=magnetization))

    sweep.write_run_script(run.jobs, out / "run.sh")

    report = ["--- Energía de superficie ---",
              f"Superficie ({run.miller[0]}{run.miller[1]}{run.miller[2]}) de "
              f"{atoms.get_chemical_formula()}",
              f"Grosores: {', '.join(str(c) for c in capas)} capas  "
              f"({', '.join(str(run.natomos[c]) for c in capas)} átomos)",
              f"Área de una cara: {area:.4f} Å²   |   vacío {vacuum:g} Å",
              f"Malla k: {grid[0]}x{grid[1]}x{grid[2]}  |  "
              + ("posiciones relajadas" if relajar else "posiciones fijas")
              + (f", {fijar} capas del fondo congeladas" if fijar else "")]
    if reducciones:
        f = max(reducciones.values())
        report.append(
            f"Celda superficial reducida a la mínima: {f:.0f} veces menos "
            f"átomos por capa.\n  γ es por unidad de área, así que el número "
            f"es el mismo y el cálculo cuesta\n  {f:.0f} veces menos. Con "
            f"--no-reduce se usa el corte tal cual sale de los hkl.")
    if not ref.simetrica:
        report.append(
            "AVISO: la losa NO es simétrica, así que sus dos caras no son la "
            "misma\n  superficie. Lo que sale de dividir entre 2 es el PROMEDIO "
            "de las dos, no γ\n  de ninguna. Para una cara concreta hace falta "
            "una losa simétrica (más capas)\n  o una referencia aparte.")
    if ref.polar:
        report.append(
            "AVISO: la losa es polar (las dos caras tienen composición "
            "distinta). Además\n  del promedio, aparece un campo eléctrico a "
            "través del vacío: usa --dipole.")
    if not relajar:
        report.append(
            "Las posiciones están FIJAS: esto da la γ sin relajar, que siempre "
            "sale alta.\n  La relajación superficial la baja típicamente entre "
            "un 5 y un 20 %. Con --relax\n  se relaja (y cuesta bastante más).")
    for w in run.avisos:
        report.append(f"AVISO: {w}")
    warn = sweep.missing_pseudo_warning(common)
    if warn:
        report.append(warn)
    report += ["", f"{len(run.jobs)} cálculos escritos en '{out.resolve()}'",
               "Córrelos con --run, o a mano con ./run.sh dentro de esa carpeta."]
    return run, "\n".join(report)


# ----------------------------------------------------------------------
# Recolección
# ----------------------------------------------------------------------
def collect(run: GammaRun, results: list = None) -> GammaRun:
    por_dir = {str(r.job.directory): r for r in (results or [])}

    def _leer(job):
        r = por_dir.get(str(job.directory))
        if r is not None and r.ok and r.result is not None:
            return r.result
        try:
            return qeout.read_xml(str(job.directory))
        except Exception:                                   # noqa: BLE001
            return None

    run.energias, run.convergido = {}, {}
    for job in run.jobs:
        res = _leer(job)
        if job.meta.get("papel") == "bulto":
            if res is not None and res.total_energy is not None:
                run.E_bulto = res.total_energy / max(1, run.natomos_bulto)
                run.bulto_ok = res.converged
            continue
        n = job.meta["capas"]
        run.energias[n] = res.total_energy if res else None
        run.convergido[n] = res.converged if res else None
    return ajustar(run)


# ----------------------------------------------------------------------
# Reporte
# ----------------------------------------------------------------------
def report(run: GammaRun) -> str:
    if not run.energias:
        raise FaltanDatos(
            "no hay resultados todavía. Corre los cálculos (--run, o ./run.sh "
            "en la carpeta) y vuelve con --collect.")
    hkl = "".join(str(m) for m in run.miller)
    L = [f"--- Energía de superficie ({hkl}) ---",
         f"Área de una cara: {run.area:.4f} Å²"
         + (f"   |   E_bulto = {run.E_bulto:.6f} eV/átomo (cálculo aparte)"
            if run.E_bulto is not None else "   |   sin cálculo de bulto")]

    hay_directo = run.E_bulto is not None
    L += ["", f"  {'capas':>6s} {'átomos':>7s} {'E_losa (eV)':>15s}"
          + (f" {'γ directa (J/m²)':>18s}" if hay_directo else "")]
    L.append("  " + "-" * (48 if hay_directo else 31))
    for n in run.capas:
        e = run.energias.get(n)
        if e is None:
            L.append(f"  {n:>6d} {run.natomos[n]:>7d} {'sin resultado':>15s}")
            continue
        fila = f"  {n:>6d} {run.natomos[n]:>7d} {e:>15.6f}"
        if hay_directo:
            g = run.gamma_directo(n)
            fila += f" {g * EV_A2_A_J_M2:>18.4f}"
        if run.convergido.get(n) is False:
            fila += "   << SIN CONVERGER"
        L.append(fila)

    if hay_directo:
        gs = [run.gamma_directo(n) for n in run.capas
              if run.gamma_directo(n) is not None]
        if len(gs) >= 2:
            deriva = (gs[-1] - gs[0]) * EV_A2_A_J_M2
            L.append("")
            L.append(f"Deriva de la γ directa entre la losa más fina y la más "
                     f"gruesa: {deriva:+.4f} J/m²")
            if abs(deriva) > 0.05:
                L.append(
                    "  No converge: crece (o baja) con el grosor en vez de "
                    "estabilizarse. Es el\n  error residual de E_bulto "
                    "multiplicado por el número de átomos, no física.\n"
                    "  El valor bueno es el del ajuste, abajo.")

    if run.gamma_ajuste is not None:
        L += ["", "Ajuste lineal E_losa(N) = 2γA + N·E_bulto  "
                  "(Fiorentini–Methfessel):",
              f"  γ = {run.gamma_ajuste:.6f} eV/Å²  =  "
              f"{run.gamma_ajuste * EV_A2_A_J_M2:.4f} J/m²",
              f"  Energía de escisión (dos caras) = "
              f"{2 * run.gamma_ajuste * EV_A2_A_J_M2:.4f} J/m²",
              f"  E_bulto de la pendiente = {run.E_bulto_ajuste:.6f} eV/átomo"
              f"   (R² = {run.r2:.6f})"]
        if run.E_bulto is not None:
            d = (run.E_bulto_ajuste - run.E_bulto) * 1000
            L.append(f"  Diferencia con el cálculo de bulto aparte: "
                     f"{d:+.2f} meV/átomo")
            if abs(d) > 5:
                L.append("  Esa diferencia es la que hacía derivar la γ "
                         "directa. Viene de que la losa y\n  el bulto no "
                         "comparten malla k (no pueden: uno tiene vacío y el "
                         "otro no).")
        if run.r2 is not None and run.r2 < 0.999:
            L.append(f"  R² = {run.r2:.5f} es bajo para una recta: o falta "
                     "convergencia en algún punto,\n  o las losas finas "
                     "todavía no tienen interior de bulto. Prueba a quitar la "
                     "más fina.")
    else:
        L.append("\nNo hay puntos suficientes para el ajuste.")

    if not run.simetrica:
        L.append("\nLa losa no es simétrica: γ es el PROMEDIO de sus dos "
                 "caras, no el de una.")
    if not run.relajado:
        L.append("Sin relajar: γ sale alta. La relajación superficial la baja "
                 "entre un 5 y un 20 %.")
    sin_conv = [n for n in run.capas if run.convergido.get(n) is False]
    if sin_conv:
        L.append(f"SIN CONVERGER: {sin_conv} capas. Esos puntos sesgan el "
                 "ajuste entero.")
    return "\n".join(L)


def export(run: GammaRun, outdir: str = ".") -> list:
    out = Path(outdir); out.mkdir(parents=True, exist_ok=True)
    f = out / "GAMMA.dat"
    hkl = "".join(str(m) for m in run.miller)
    lines = [provenance.header(
        f"energia de superficie ({hkl})",
        {"area_A2": f"{run.area:.5f}",
         "gamma_eV_A2": run.gamma_ajuste,
         "gamma_J_m2": (run.gamma_ajuste * EV_A2_A_J_M2
                        if run.gamma_ajuste else None),
         "E_bulto_ajuste_eV_at": run.E_bulto_ajuste, "R2": run.r2}),
        f"# {'capas':>6s} {'atomos':>8s} {'E_losa(eV)':>18s} "
        f"{'gamma_directa(J/m2)':>21s}"]
    for n in run.capas:
        if run.energias.get(n) is None:
            continue
        g = run.gamma_directo(n)
        lines.append(f"{n:8d} {run.natomos[n]:8d} {run.energias[n]:18.8f} "
                     + (f"{g * EV_A2_A_J_M2:21.6f}" if g is not None
                        else f"{'nan':>21s}"))
    f.write_text("\n".join(lines) + "\n")
    txt = out / "GAMMA.txt"
    txt.write_text(report(run) + "\n")
    return [str(f), str(txt)]


def plot(run: GammaRun, outfile: str = "gamma", formats="pdf,png",
         theme: str = None, size: str = None, family: str = None,
         background: str = None, palette=None, usetex: bool = None,
         width="single", journal: str = "generic", aspect: float = 0.72,
         mono: bool = False, dpi: int = None) -> list:
    """γ directa contra grosor, con la del ajuste como línea horizontal."""
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError as exc:                              # pragma: no cover
        raise RuntimeError("matplotlib no está instalado.") from exc

    ns = [n for n in run.capas if run.energias.get(n) is not None]
    if len(ns) < 2:
        raise FaltanDatos("hacen falta al menos dos grosores para graficar.")
    st = qstyle.apply(theme, size=size, family=family, background=background,
                      palette=palette, usetex=usetex, mono=mono)
    fig, ax = qstyle.new_figure(width, journal, aspect)
    cols = qstyle.palette(2, mono=mono)

    directos = [(n, run.gamma_directo(n)) for n in ns]
    directos = [(n, g * EV_A2_A_J_M2) for n, g in directos if g is not None]
    if directos:
        ax.plot([d[0] for d in directos], [d[1] for d in directos],
                marker="o", ms=4, lw=st["line"], color=cols[0],
                label="directa, con E$_{bulto}$ aparte")
    if run.gamma_ajuste is not None:
        g = run.gamma_ajuste * EV_A2_A_J_M2
        ax.axhline(g, color=cols[1], lw=st["line"], dashes=[4.0, 2.0],
                   label=f"ajuste: {g:.3f} J/m$^2$")
    ax.set_xlabel("capas de la losa")
    ax.set_ylabel(r"$\gamma$ (J/m$^2$)")
    ax.set_xticks(ns)
    if directos or run.gamma_ajuste is not None:
        ax.legend(frameon=False, fontsize=st["legend"])
    written = qstyle.save(fig, outfile, formats, dpi=dpi, modulo="gamma")
    plt.close(fig)
    return written
