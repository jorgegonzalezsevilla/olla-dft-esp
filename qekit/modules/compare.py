# Olla-DFT — command-line toolkit for Quantum ESPRESSO
# Copyright (C) 2026 Jorge Enrique González Sevilla
# SPDX-License-Identifier: AGPL-3.0-or-later
# This program is free software: you can redistribute it and/or modify it
# under the terms of the GNU Affero General Public License as published by the Free
# Software Foundation, either version 3 of the License, or (at your option)
# any later version. See the LICENSE file for details.

"""Comparación segura de corridas de Quantum ESPRESSO.

La comparación separa tres cosas que suelen mezclarse: parámetros, métricas
y diferencias energéticas. Una energía solo se resta cuando ``audit`` dice
que las corridas son comparables y no son salidas NSCF/bandas sin energía
utilizable.
"""

from __future__ import annotations

import json
import math
from pathlib import Path

from qekit import __version__
from qekit.core import provenance
from qekit.core.errors import ErrorDeUso
from qekit.modules import audit


METRICS = (
    ("energia_eV", "energía total", "eV"),
    ("energia_por_atomo_eV", "energía por átomo", "eV/átomo"),
    ("gap_eV", "gap", "eV"),
    ("volumen_A3", "volumen", "Å³"),
    ("presion_GPa", "presión", "GPa"),
    ("fuerza_max", "fuerza máxima", "eV/Å"),
    ("magnetizacion", "magnetización", "μB/celda"),
    ("wall_s", "tiempo de pared", "s"),
)


def _finite(value):
    try:
        return value is not None and math.isfinite(float(value))
    except (TypeError, ValueError):
        return False


def _formula(result):
    try:
        from ase import Atoms
        return Atoms(symbols=result.symbols).get_chemical_formula()
    except Exception:  # noqa: BLE001
        return "".join(result.symbols or [])


def _record(run) -> dict:
    result = run.result
    if result is None:
        return {"path": run.path, "ok": False, "error": run.error}
    natoms = len(result.symbols) or 1
    energy = result.total_energy
    return {
        "path": str(Path(run.path).resolve()), "ok": True,
        "error": run.error, "origen": run.origen,
        "formula": _formula(result), "natoms": len(result.symbols),
        "calculation": result.calculation, "functional": result.functional,
        "pseudos": dict(result.pseudo_files), "ecutwfc": result.ecutwfc,
        "ecutrho": result.ecutrho,
        "kgrid": "x".join(map(str, result.kgrid)) if result.kgrid else None,
        "convergido": result.converged,
        "energia_eV": energy,
        "energia_por_atomo_eV": energy / natoms if _finite(energy) else None,
        "gap_eV": (result.lumo - result.homo
                    if _finite(result.homo) and _finite(result.lumo) else None),
        "volumen_A3": result.volume, "presion_GPa": result.pressure,
        "fuerza_max": result.max_force,
        "magnetizacion": result.total_magnetization,
        "wall_s": result.wall_time,
    }


def _reference_index(records, reference=None):
    if reference is None:
        return 0
    if isinstance(reference, int):
        if 0 <= reference < len(records):
            return reference
        raise ErrorDeUso(f"--reference debe estar entre 0 y {len(records) - 1}.")
    wanted = str(Path(reference).resolve())
    for index, record in enumerate(records):
        if record["path"] == wanted:
            return index
    raise ErrorDeUso(f"no encuentro la corrida de referencia '{reference}'.")


def compare(paths, reference=None) -> dict:
    paths = list(paths or [])
    if len(paths) < 2:
        raise ErrorDeUso("compare necesita al menos dos carpetas o XML.")
    runs = audit.collect(paths)
    checked = audit.audit(runs)
    records = [_record(run) for run in runs]
    ref_index = _reference_index(records, reference)
    ref = records[ref_index]
    comparable_energy = bool(checked["comparables"] and ref.get("ok") and
                              ref.get("convergido") is not False and
                              ref.get("calculation", "").lower() not in
                              ("nscf", "bands") and
                              _finite(ref.get("energia_eV")))
    rows = []
    for index, record in enumerate(records):
        row = dict(record)
        row["is_reference"] = index == ref_index
        for key, _label, _unit in METRICS:
            value = record.get(key)
            row[f"delta_{key}"] = (
                float(value) - float(ref[key])
                if comparable_energy and key in ("energia_eV", "energia_por_atomo_eV")
                and _finite(value) and _finite(ref.get(key)) else None)
            if key not in ("energia_eV", "energia_por_atomo_eV"):
                row[f"delta_{key}"] = (
                    float(value) - float(ref[key])
                    if _finite(value) and _finite(ref.get(key)) else None)
        rows.append(row)
    return {
        "qekit_version": __version__, "generated": provenance.fields()["generado"],
        "reference": ref["path"], "comparable_energy": comparable_energy,
        "audit": {
            "comparables": checked["comparables"],
            "difieren": [(key, [repr(v) for v in values])
                         for key, values in checked["difieren"]],
            "failed": [x.path for x in checked["fallidos"]],
            "not_converged": [x.path for x in checked["no_convergidos"]],
        },
        "runs": rows,
    }


def report(result: dict) -> str:
    lines = ["--- Comparación de corridas ---",
             f"Corridas: {len(result['runs'])}  | referencia: "
             f"{result['reference']}"]
    if result["comparable_energy"]:
        lines.append("Energías: comparables; las diferencias se calculan contra la referencia.")
    else:
        lines.append("Energías: NO se restan; faltan comparabilidad, convergencia o energía utilizable.")
    lines.append("")
    header = f"{'corrida':28s} {'fórmula':10s} {'E/átomo':>14s} {'gap':>10s} {'ΔE/át':>14s}"
    lines.append(header)
    lines.append("-" * len(header))
    for row in result["runs"]:
        path = Path(row["path"]).name[:27]
        e = row.get("energia_por_atomo_eV")
        gap = row.get("gap_eV")
        de = row.get("delta_energia_por_atomo_eV")
        e_text = f"{e:14.6f}" if _finite(e) else f"{'-':>14s}"
        gap_text = f"{gap:10.4f}" if _finite(gap) else f"{'-':>10s}"
        de_text = f"{de:14.6g}" if _finite(de) else f"{'-':>14s}"
        lines.append(f"{path:28s} {(row.get('formula') or '?'):10s} "
                     f"{e_text} {gap_text} {de_text}")
    if result["audit"]["difieren"]:
        lines += ["", "Parámetros que impiden restar energías:"]
        for key, values in result["audit"]["difieren"]:
            lines.append(f"  - {key}: {', '.join(values)}")
    if result["audit"]["not_converged"]:
        lines += ["", "Corridas no convergidas:"]
        lines.extend(f"  - {path}" for path in result["audit"]["not_converged"])
    return "\n".join(lines)


def export(result: dict, destination="comparison.json") -> Path:
    target = Path(destination)
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_text(json.dumps(result, ensure_ascii=False, indent=2) + "\n",
                      encoding="utf-8")
    return target
