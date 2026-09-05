# Olla-DFT — command-line toolkit for Quantum ESPRESSO
# Copyright (C) 2026 Jorge Enrique González Sevilla
# SPDX-License-Identifier: AGPL-3.0-or-later
# This program is free software: you can redistribute it and/or modify it
# under the terms of the GNU Affero General Public License as published by the Free
# Software Foundation, either version 3 of the License, or (at your option)
# any later version. See the LICENSE file for details.

"""Validaciones estructurales y de integridad antes de gastar CPU."""

from __future__ import annotations

import math
from pathlib import Path

from qekit.modules import project


def _check(code, title, level, detail, evidence=""):
    return {"code": code, "title": title, "level": level,
            "detail": detail, "evidence": evidence}


def check(root, data) -> list:
    checks = []
    sources = data.get("sources", [])
    for source in sources:
        path = Path(source["path"])
        path = path if path.is_absolute() else root / path
        if not path.is_file():
            checks.append(_check("source.missing", "Archivo de entrada", "fail",
                                 "no existe la fuente registrada", source["path"]))
            continue
        current = project.sha256_file(path)
        if current != source.get("sha256"):
            checks.append(_check("source.changed", "Archivo de entrada", "fail",
                                 "el SHA-256 actual no coincide con el manifiesto",
                                 source["path"]))
        else:
            checks.append(_check("source.locked", "Archivo de entrada", "ok",
                                 "archivo presente y bloqueado por SHA-256",
                                 source["path"]))
        if source.get("kind") == "structure":
            checks.extend(_structure_checks(path))

    seen_outputs = {}
    for task in data.get("tasks", []):
        try:
            project._valid_command(task.get("command", ""))
            checks.append(_check("task.command", f"Comando {task.get('id', '?')}", "ok",
                                 "sintaxis Olla-DFT válida"))
        except Exception as exc:  # noqa: BLE001
            checks.append(_check("task.command", f"Comando {task.get('id', '?')}", "fail",
                                 str(exc)))
        for output in task.get("outputs", []):
            previous = seen_outputs.get(output)
            if previous:
                checks.append(_check("task.output_collision", "Salidas del workflow", "warn",
                                     "dos tareas declaran la misma salida", output))
            seen_outputs[output] = task.get("id")

    quantities = data.get("metadata", {}).get("quantities", [])
    if isinstance(quantities, dict):
        quantities = [quantities]
    for quantity in quantities:
        if not isinstance(quantity, dict):
            checks.append(_check("quantity.format", "Magnitudes", "fail",
                                 "cada magnitud debe ser un objeto con value y unit"))
            continue
        value, unit = quantity.get("value"), quantity.get("unit", "")
        try:
            valid = math.isfinite(float(value)) and bool(str(unit).strip())
        except (TypeError, ValueError):
            valid = False
        checks.append(_check("quantity.finite", str(quantity.get("name", "Magnitud")),
                             "ok" if valid else "fail",
                             "valor finito con unidad explícita" if valid else
                             "valor no finito o sin unidad"))
        if "uncertainty" in quantity:
            try:
                uncertainty_ok = (math.isfinite(float(quantity["uncertainty"])) and
                                  float(quantity["uncertainty"]) >= 0 and valid)
            except (TypeError, ValueError):
                uncertainty_ok = False
            checks.append(_check(
                "quantity.uncertainty", str(quantity.get("name", "Magnitud")),
                "ok" if uncertainty_ok else "fail",
                "incertidumbre finita, no negativa y con unidad" if uncertainty_ok
                else "la incertidumbre debe ser un número finito no negativo"))
    checks.extend(result_checks(root, data))
    if not checks:
        checks.append(_check("project.empty", "Validación avanzada", "warn",
                             "no hay fuentes ni tareas suficientes para validar"))
    return checks


def result_checks(root, data) -> list:
    """Validaciones de integridad sobre resultados ya ingeridos."""
    try:
        from qekit.modules import results
        rows = results.list_results(results.project_db(root), limit=10000)
    except Exception as exc:  # noqa: BLE001
        return [_check("results.read", "Resultados", "fail",
                        f"no se pudo leer el índice: {exc}")]
    checks = []
    for row in rows:
        label = f"Resultado {row.get('id', '?')[:10]}"
        status = row.get("status")
        metrics = row.get("metrics", {})
        energy = metrics.get("energy_total", {}).get("value")
        volume = metrics.get("volume", {}).get("value")
        gap = metrics.get("gap", {}).get("value")
        if status == "invalid":
            checks.append(_check("result.invalid", label, "fail",
                                 "la salida no pudo interpretarse",
                                 row.get("path", "")))
            continue
        if status == "not_converged":
            checks.append(_check("result.not_converged", label, "fail",
                                 "el cálculo terminó sin convergencia; no debe promocionarse",
                                 row.get("path", "")))
        elif status == "parsed_no_energy":
            checks.append(_check("result.no_energy", label, "warn",
                                 "la salida se leyó, pero no contiene energía total utilizable",
                                 row.get("calculation", "")))
        elif status == "converged" and energy is None:
            checks.append(_check("result.energy_missing", label, "fail",
                                 "figura como convergido pero carece de energía total"))
        if volume is not None and volume <= 0:
            checks.append(_check("result.volume", label, "fail",
                                 "el volumen registrado no es positivo"))
        if gap is not None and gap < -1e-9:
            checks.append(_check("result.gap", label, "fail",
                                 "el gap registrado es negativo; revisar HOMO/LUMO"))
        review = row.get("review", {})
        if review.get("status") == "rejected":
            checks.append(_check("result.review_rejected", label, "fail",
                                 "la revisión humana rechazó este resultado",
                                 review.get("note", "")))
    return checks


def _structure_checks(path: Path) -> list:
    try:
        from qekit.core import structure
        atoms = structure.load(path)
    except Exception as exc:  # noqa: BLE001
        return [_check("structure.parse", "Estructura", "fail",
                        f"no se pudo leer: {exc}", str(path))]
    checks = []
    symbols = list(atoms.get_chemical_symbols())
    positions = atoms.get_positions()
    finite = bool(len(symbols)) and all(math.isfinite(float(x))
                                        for x in positions.ravel())
    checks.append(_check("structure.geometry", "Geometría", "ok" if finite else "fail",
                         f"{len(symbols)} átomos y coordenadas finitas" if finite else
                         "la estructura no tiene átomos o contiene NaN/Inf",
                         str(path)))
    if atoms.cell.volume > 1e-9:
        checks.append(_check("structure.cell", "Celda", "ok",
                             f"volumen {atoms.cell.volume:.6g} Å³"))
    elif atoms.pbc.any():
        checks.append(_check("structure.cell", "Celda", "fail",
                             "la estructura periódica tiene volumen nulo"))
    try:
        distances = atoms.get_all_distances(mic=False)
        positive = distances[distances > 1e-8]
        minimum = float(positive.min()) if positive.size else None
    except Exception:  # noqa: BLE001
        minimum = None
    if minimum is not None and minimum < 0.5:
        checks.append(_check("structure.distance", "Distancias", "warn",
                             f"distancia mínima sospechosamente pequeña: {minimum:.4g} Å"))
    else:
        checks.append(_check("structure.distance", "Distancias", "ok",
                             "no se detectan distancias atómicas anómalas"))
    return checks


def report(checks: list) -> str:
    levels = {level: sum(c["level"] == level for c in checks)
              for level in ("ok", "warn", "fail")}
    lines = ["--- Validación avanzada del proyecto ---",
             f"OK: {levels['ok']}  avisos: {levels['warn']}  fallos: {levels['fail']}"]
    for check_item in checks:
        mark = {"ok": "OK", "warn": "AVISO", "fail": "FALLO"}[check_item["level"]]
        lines.append(f"  [{mark:5s}] {check_item['title']}: {check_item['detail']}")
        if check_item.get("evidence"):
            lines.append(f"         {check_item['evidence']}")
    return "\n".join(lines)
