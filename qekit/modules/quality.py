# Olla-DFT — command-line toolkit for Quantum ESPRESSO
# Copyright (C) 2026 Jorge Enrique González Sevilla
# SPDX-License-Identifier: GPL-3.0-or-later
# This program is free software: you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by the Free
# Software Foundation, either version 3 of the License, or (at your option)
# any later version. See the LICENSE file for details.

"""Puerta de calidad científica para un Project Hub de Olla-DFT."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

from qekit.modules import project


@dataclass(frozen=True)
class Check:
    code: str
    title: str
    level: str  # ok, warn, fail
    detail: str
    evidence: str = ""


def evaluate(root: Path, data: dict) -> dict:
    checks = []
    state = project.status(root, data)
    if data.get("sources"):
        if state["changed_sources"]:
            checks.append(Check(
                "sources.changed", "Fuentes reproducibles", "fail",
                "hay entradas que ya no coinciden con el SHA-256 registrado",
                ", ".join(state["changed_sources"])))
        else:
            checks.append(Check(
                "sources.locked", "Fuentes reproducibles", "ok",
                "las fuentes registradas conservan tamaño y SHA-256"))
    else:
        checks.append(Check("sources.missing", "Fuentes reproducibles", "fail",
                            "el proyecto no registra ninguna fuente"))

    # El bloqueo no sustituye a un contenedor ni a una receta de instalación,
    # pero sí hace visible si el equipo que reabre el proyecto cambió de
    # dependencias, Python o binarios de Quantum ESPRESSO.
    lock_path = root / project.PROJECT_DIR / "environment.lock.json"
    if lock_path.is_file():
        try:
            from qekit.modules import environment
            locked = environment.verify(root)
            checks.append(Check(
                "environment.locked", "Entorno reproducible",
                "ok" if locked["ok"] else "warn",
                "Python, dependencias y binarios coinciden" if locked["ok"] else
                "el entorno actual difiere del bloqueo guardado",
                ", ".join(locked.get("changed", []))))
        except Exception as exc:  # noqa: BLE001
            checks.append(Check("environment.unreadable", "Entorno reproducible",
                                "warn", "no se pudo verificar el bloqueo",
                                str(exc)))
    elif data.get("sources"):
        checks.append(Check(
            "environment.missing", "Entorno reproducible", "warn",
            "falta environment.lock.json; créalo antes de compartir o publicar",
            "olla-dft project environment"))

    advanced = data.get("metadata", {}).get("advanced_validation")
    if advanced:
        advanced_level = ("fail" if not advanced.get("passed") else
                          "warn" if advanced.get("warnings") else "ok")
        checks.append(Check(
            "validation.advanced", "Validación avanzada",
            advanced_level,
            "estructura, comandos, unidades y salidas revisados" if advanced_level == "ok"
            else "la validación avanzada encontró fallos" if advanced_level == "fail"
            else "validación avanzada completada con avisos",
            str(advanced.get("at", ""))))

    # La presencia del índice no es una prueba de validez física, pero sí
    # evita que un proyecto con tareas terminadas pierda silenciosamente sus
    # salidas normalizadas antes de publicarse.
    try:
        from qekit.modules import results
        indexed = results.summary(results.project_db(root))
    except Exception as exc:  # noqa: BLE001
        checks.append(Check(
            "results.unreadable", "Resultados trazables", "fail",
            "el índice de resultados existe pero no se puede leer",
            str(exc).splitlines()[0]))
        indexed = {"count": 0, "by_status": {}}
    if indexed.get("count"):
        invalid = indexed.get("by_status", {}).get("invalid", 0)
        level = "warn" if invalid else "ok"
        detail = (f"{indexed['count']} resultado(s) conservan métricas y hashes"
                  + (f"; {invalid} no se pudieron interpretar" if invalid else ""))
        checks.append(Check("results.indexed", "Resultados trazables", level, detail,
                            str(results.project_db(root))))
    elif any(x.get("status") == "succeeded" for x in data.get("tasks", [])):
        checks.append(Check(
            "results.missing", "Resultados trazables", "warn",
            "hay tareas terminadas pero todavía no se ha ingerido ningún XML"))

    counts = state["counts"]
    if counts.get("failed"):
        checks.append(Check(
            "tasks.failed", "Workflow ejecutado", "fail",
            f"{counts['failed']} tarea(s) terminaron con error"))
    elif counts.get("cancelled"):
        checks.append(Check(
            "tasks.cancelled", "Workflow ejecutado", "warn",
            f"{counts['cancelled']} tarea(s) quedaron canceladas; reanuda y revisa antes de publicar"))
    elif data.get("tasks"):
        checks.append(Check(
            "tasks.state", "Workflow ejecutado",
            "ok" if not counts.get("pending") and not counts.get("blocked")
            else "warn",
            "todas las tareas terminaron" if not counts.get("pending")
            and not counts.get("blocked") else
            "hay tareas planificadas que todavía no se han ejecutado"))
    else:
        checks.append(Check("tasks.empty", "Workflow ejecutado", "warn",
                            "todavía no hay tareas en el proyecto"))

    if any(x.get("status") == "succeeded" for x in data.get("tasks", [])):
        checks.append(Check(
            "provenance.logs", "Logs y procedencia", "ok",
            "las tareas ejecutadas tienen estado persistente; exporta el "
            "snapshot para archivarlo"))
    else:
        checks.append(Check(
            "provenance.pending", "Logs y procedencia", "warn",
            "aún no hay tareas ejecutadas para auditar"))

    # No se presenta como aprobación de publicación: solo constata si el
    # proyecto ha pasado explícitamente la suite independiente.
    selftest = data.get("metadata", {}).get("selftest", {})
    if selftest.get("passed"):
        checks.append(Check("selftest.passed", "Fórmulas independientes", "ok",
                            "el proyecto registra selftest aprobado",
                            str(selftest.get("at"))))
    else:
        checks.append(Check(
            "selftest.missing", "Fórmulas independientes", "warn",
            "no consta una ejecución aprobada de olla-dft selftest; esto no "
            "invalida una exploración, pero falta evidencia independiente"))

    fails = sum(c.level == "fail" for c in checks)
    warns = sum(c.level == "warn" for c in checks)
    score = max(0, 100 - 40 * fails - 10 * warns)
    verdict = "bloqueado" if fails else "revisar" if warns else "listo"
    return {"checks": checks, "fails": fails, "warnings": warns,
            "score": score, "verdict": verdict}


def report(result: dict) -> str:
    lines = ["--- Puerta de calidad científica ---",
             f"Veredicto: {result['verdict'].upper()}  |  puntuación orientativa: "
             f"{result['score']}/100"]
    for check in result["checks"]:
        mark = {"ok": "OK", "warn": "AVISO", "fail": "FALLO"}[check.level]
        lines.append(f"  [{mark:5s}] {check.title}: {check.detail}")
        if check.evidence:
            lines.append(f"         evidencia: {check.evidence}")
    lines.append("\nEsta puerta organiza evidencia; no sustituye revisión científica ni "
                 "autoriza publicar automáticamente.")
    return "\n".join(lines)
