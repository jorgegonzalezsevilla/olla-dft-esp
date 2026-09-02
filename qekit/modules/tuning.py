# Olla-DFT — command-line toolkit for Quantum ESPRESSO
# Copyright (C) 2026 Jorge Enrique González Sevilla
# SPDX-License-Identifier: GPL-3.0-or-later
# This program is free software: you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by the Free
# Software Foundation, either version 3 of the License, or (at your option)
# any later version. See the LICENSE file for details.

"""Recomendación adaptativa a partir de una serie de convergencia."""

from __future__ import annotations

import json
import math
from pathlib import Path

from qekit.core.errors import ErrorDeUso


def read(path) -> list:
    """Lee ``CONVERGENCIA.dat`` sin asumir que todos los puntos terminaron."""
    rows = []
    for number, line in enumerate(Path(path).read_text(encoding="utf-8",
                                                         errors="replace").splitlines(), 1):
        if not line.strip() or line.lstrip().startswith("#"):
            continue
        tokens = line.split()
        if len(tokens) < 3:
            continue
        try:
            value, energy_ry, delta = map(float, tokens[:3])
        except ValueError:
            continue
        if not all(math.isfinite(v) for v in (value, energy_ry, delta)):
            continue
        rows.append({"line": number, "value": value, "energy_ry": energy_ry,
                     "delta_mev_atom": abs(delta)})
    if not rows:
        raise ErrorDeUso(f"'{path}' no contiene filas numéricas de convergencia.")
    return rows


def analyze(path, threshold=None) -> dict:
    rows = read(path)
    if threshold is None:
        threshold = 1.0
    if threshold <= 0:
        raise ErrorDeUso("el umbral debe ser positivo.")
    index = None
    for i in range(len(rows)):
        if all(row["delta_mev_atom"] <= threshold for row in rows[i:]):
            index = i
            break
    values = [row["value"] for row in rows]
    if index is None:
        status = "extend"
        recommendation = _next_value(values)
        reason = "ningún punto mantiene toda la cola dentro del umbral"
    elif index == len(rows) - 1:
        status = "confirm"
        recommendation = _next_value(values)
        reason = "solo el último punto cumple; hace falta un punto más para confirmar"
    else:
        status = "ready"
        recommendation = rows[index]["value"]
        reason = "desde este punto toda la cola queda dentro del umbral"
    return {"file": str(Path(path).resolve()), "threshold": float(threshold),
            "rows": rows, "converged_index": index, "status": status,
            "recommended_value": recommendation, "reason": reason}


def _next_value(values):
    if len(values) < 2:
        return values[-1] * 1.25 if values[-1] > 0 else values[-1] + 1.0
    diffs = [b - a for a, b in zip(values, values[1:]) if b > a]
    if diffs:
        step = sorted(diffs)[len(diffs) // 2]
        return values[-1] + max(step, abs(values[-1]) * 0.10)
    return values[-1] * 1.25 if values[-1] > 0 else values[-1] + 1.0


def report(result: dict) -> str:
    lines = ["--- Recomendación adaptativa de convergencia ---",
             f"Archivo: {result['file']}",
             f"Umbral: {result['threshold']:g} meV/átomo"]
    index = result["converged_index"]
    if index is None:
        lines.append("Estado: EXTENDER — la serie todavía no converge.")
    elif result["status"] == "confirm":
        lines.append("Estado: CONFIRMAR — el último punto no basta como evidencia.")
    else:
        lines.append(f"Estado: LISTO — usar desde el punto {index + 1} de la serie.")
    lines.append(f"Recomendación: probar valor {result['recommended_value']:g}.")
    lines.append(f"Motivo: {result['reason']}.")
    lines.append("La propiedad energía puede converger antes que fuerzas, fonones o tensores.")
    return "\n".join(lines)


def export(result: dict, destination="CONVERGENCIA_RECOMENDACION.json") -> Path:
    target = Path(destination)
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_text(json.dumps(result, ensure_ascii=False, indent=2) + "\n",
                      encoding="utf-8")
    return target
