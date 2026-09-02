# Olla-DFT — command-line toolkit for Quantum ESPRESSO
# Copyright (C) 2026 Jorge Enrique González Sevilla
# SPDX-License-Identifier: GPL-3.0-or-later
# This program is free software: you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by the Free
# Software Foundation, either version 3 of the License, or (at your option)
# any later version. See the LICENSE file for details.

"""Energía de exfoliación de un material laminar.

Se compara la energía del cristal con la de una monocapa aislada en vacío:

    E_exf = (E_monocapa − E_bulk / N_capas) / A

con A el área en el plano de la celda. Es el costo de separar una capa del
resto, la magnitud que decide si un material es exfoliable (los laminares
tipicos rondan 0.2–0.6 J/m²; grafito ≈ 0.35 J/m² experimental).

Advertencia física seria, que el reporte repite: la cohesión entre capas es
en gran parte dispersión de van der Waals, y ni LDA ni PBE la describen —
PBE prácticamente no liga las capas y LDA liga por un error afortunado.
Para números publicables hay que activar una corrección de dispersión
(`--vdw grimme-d2` en QE ≥ 5, `--vdw grimme-d3` en QE ≥ 7.1) o un funcional
vdW. El módulo funciona sin ella, pero lo dice.
"""

from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

from qekit.core import layers as layers_mod
from qekit.core import qeout
from qekit.modules import sweep
from qekit.core.errors import ErrorDeUso

EV_A2_TO_J_M2 = 16.02176634


@dataclass
class ExfoliationRun:
    n_layers: int = 0
    area: float = None            # Å² de la celda en el plano
    natoms_bulk: int = 0
    natoms_slab: int = 0
    vacuum: float = 20.0
    vdw: str = None
    jobs: list = field(default_factory=list)   # [bulk, slab]
    E_bulk: float = None          # eV
    E_slab: float = None          # eV


def prepare(atoms, outdir: str = "exfoliacion", vacuum: float = 20.0,
            vdw: str = None, tol: float = layers_mod.DEFAULT_TOL,
            pseudo_dir: str = None, ecutwfc: float = None,
            ecutrho: float = None, kspacing: float = None,
            insulator: bool = False, relax_slab: bool = False) -> tuple:
    """Prepara los dos cálculos (bulk y monocapa) con parámetros idénticos."""
    ana = layers_mod.analyze(atoms, tol)
    if not ana.layers:
        raise ErrorDeUso(
            "no se detectaron capas en la estructura "
            f"(dimensionalidad {ana.dimensionality}). Ajusta --tol si crees "
            "que sí es laminar."
        )

    out = Path(outdir); out.mkdir(parents=True, exist_ok=True)
    common = sweep.prepare_common(atoms, pseudo_dir, ecutwfc, ecutrho, insulator)

    grid_bulk = sweep.default_grid(atoms, kspacing)
    # la monocapa lleva la misma malla en el plano y 1 punto en la normal
    grid_slab = list(grid_bulk)
    grid_slab[ana.stacking_axis] = 1
    grid_slab = tuple(grid_slab)

    slab = layers_mod.make_slab(atoms, ana, vacuum=vacuum)

    # área en el plano
    others = [k for k in range(3) if k != ana.stacking_axis]
    area = float(np.linalg.norm(np.cross(atoms.cell.array[others[0]],
                                         atoms.cell.array[others[1]])))

    run = ExfoliationRun(
        n_layers=len(ana.layers), area=area, natoms_bulk=len(atoms),
        natoms_slab=len(slab), vacuum=vacuum, vdw=vdw,
    )

    report = ["--- Energía de exfoliación ---",
              f"Estructura: {atoms.get_chemical_formula()}  |  "
              f"{run.n_layers} capa(s) por celda, "
              f"apiladas según el eje {'abc'[ana.stacking_axis]}",
              f"Espaciado basal: {ana.basal_spacing:.4f} Å  |  "
              f"área en el plano: {area:.4f} Å²",
              f"Monocapa: {run.natoms_slab} átomos con {vacuum:g} Å de vacío",
              f"Mallas k: bulk {grid_bulk[0]}x{grid_bulk[1]}x{grid_bulk[2]}, "
              f"monocapa {grid_slab[0]}x{grid_slab[1]}x{grid_slab[2]}"]
    if vdw:
        report.append(f"Corrección de dispersión: vdw_corr = '{vdw}'")
        lda_like = any(t in p["filename"].lower()
                       for p in common["pseudos"].values()
                       for t in ("pz", "lda", "pw92"))
        if lda_like and vdw.lower() in ("grimme-d2", "dft-d", "grimme-d3"):
            report.append(
                "AVISO: los pseudopotenciales parecen LDA y las correcciones de\n"
                "Grimme están parametrizadas para PBE: combinarlas cuenta la\n"
                "dispersión dos veces (LDA ya sobreliga las capas). Usa pseudos\n"
                "PBE con Grimme, o LDA sin corrección solo como referencia."
            )
    else:
        report.append(
            "SIN corrección de van der Waals: PBE apenas liga las capas y LDA\n"
            "liga por cancelación de errores. Para resultados publicables usa\n"
            "--vdw grimme-d2 (QE >= 5) o --vdw grimme-d3 (QE >= 7.1)."
        )
    warn = sweep.missing_pseudo_warning(common)
    if warn:
        report.append(warn)

    run.jobs.append(sweep.write_scf_job(
        atoms, common, out / "bulk", "bulk", grid_bulk,
        meta={"which": "bulk"}, vdw=vdw))
    run.jobs.append(sweep.write_scf_job(
        slab, common, out / "monocapa", "monocapa", grid_slab,
        meta={"which": "slab"},
        calculation="relax" if relax_slab else "scf", vdw=vdw))

    sweep.write_run_script(run.jobs, out / "run.sh")
    report += ["", f"2 cálculos escritos en '{out.resolve()}'",
               "Córrelos con --run, o a mano con ./run.sh dentro de esa carpeta."]
    return run, "\n".join(report)


def collect(run: ExfoliationRun, results: list = None) -> ExfoliationRun:
    energies = {}
    if results is not None:
        for r in results:
            if r.ok:
                energies[r.job.meta.get("which")] = r.energy
    else:
        for job in run.jobs:
            try:
                energies[job.meta.get("which")] = qeout.read_xml(
                    str(job.directory)).total_energy
            except Exception:
                pass
    run.E_bulk = energies.get("bulk")
    run.E_slab = energies.get("slab")
    return run


def report_result(run: ExfoliationRun) -> str:
    lines = ["--- Resultado de exfoliación ---"]
    if run.E_bulk is None or run.E_slab is None:
        faltan = [w for w, e in (("bulk", run.E_bulk), ("monocapa", run.E_slab))
                  if e is None]
        lines.append(f"Faltan cálculos terminados: {', '.join(faltan)}.")
        return "\n".join(lines)

    e_per_layer = run.E_bulk / run.n_layers
    diff = run.E_slab - e_per_layer          # eV por celda de monocapa
    per_area = diff / run.area               # eV/Å²
    per_atom = diff / run.natoms_slab * 1000.0

    lines += [
        f"E(bulk)     = {run.E_bulk / qeout.RY_EV:16.8f} Ry  "
        f"({run.natoms_bulk} átomos, {run.n_layers} capas)",
        f"E(monocapa) = {run.E_slab / qeout.RY_EV:16.8f} Ry  "
        f"({run.natoms_slab} átomos)",
        "",
        f"E_exf = {per_area * EV_A2_TO_J_M2:8.4f} J/m²"
        f"   = {per_area * 1000.0:8.2f} meV/Å²"
        f"   = {per_atom:8.2f} meV/átomo",
    ]
    if diff < 0:
        lines.append(
            "\nAVISO: la energía salió NEGATIVA (la monocapa sería más estable "
            "que el bulk).\nCasi siempre significa que falta la corrección vdW "
            "o que algún cálculo no está\nbien convergido."
        )
    if not run.vdw:
        lines.append(
            "\nRecordatorio: sin corrección de dispersión este número no es "
            "comparable con\nexperimento (referencia: grafito ≈ 0.35 J/m²)."
        )
    return "\n".join(lines)
