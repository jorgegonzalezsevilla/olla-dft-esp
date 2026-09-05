# Olla-DFT — command-line toolkit for Quantum ESPRESSO
# Copyright (C) 2026 Jorge Enrique González Sevilla
# SPDX-License-Identifier: AGPL-3.0-or-later
# This program is free software: you can redistribute it and/or modify it
# under the terms of the GNU Affero General Public License as published by the Free
# Software Foundation, either version 3 of the License, or (at your option)
# any later version. See the LICENSE file for details.

"""Informe PDF compacto de un proyecto Olla-DFT.

El dashboard HTML sigue siendo la vista interactiva. Este módulo añade una
salida estable para adjuntar a una revisión, tesis o expediente sin depender
de un navegador ni de un servicio web.
"""

from __future__ import annotations

import textwrap
from pathlib import Path

from qekit.modules import project, quality, results


def _lines(root, data):
    state = project.status(root, data)
    gate = quality.evaluate(root, data)
    indexed = results.summary(results.project_db(root))
    counts = state["counts"]
    lines = [
        "Olla-DFT — informe reproducible del proyecto",
        "",
        f"Proyecto: {data['name']}",
        f"Carpeta: {root}",
        f"Generado por Olla-DFT: {data.get('metadata', {}).get('qekit_version', '?')}",
        "",
        f"Puerta de calidad: {gate['verdict'].upper()} ({gate['score']}/100)",
        f"Resultados indexados: {indexed['count']}",
        "Tareas: " + ", ".join(f"{key}={value}" for key, value in counts.items()),
        "",
        "Fuentes:",
    ]
    lines.extend(f"  - {item.get('path')} · SHA-256 {item.get('sha256', '')}"
                 for item in data.get("sources", []))
    if not data.get("sources"):
        lines.append("  - ninguna fuente registrada")
    lines.extend(["", "Comprobaciones:"])
    for item in gate["checks"]:
        lines.append(f"  [{item.level.upper()}] {item.title}: {item.detail}")
        if item.evidence:
            lines.append(f"      {item.evidence}")
    lines.extend(["", "Campañas:"])
    task_index = {task.get("id"): task for task in data.get("tasks", [])}
    if data.get("campaigns"):
        for campaign in data["campaigns"]:
            states = {}
            for task_id in campaign.get("tasks", []):
                status = task_index.get(task_id, {}).get("status", "missing")
                states[status] = states.get(status, 0) + 1
            state_text = ", ".join(f"{key}={value}"
                                   for key, value in sorted(states.items()))
            lines.append(f"  - {campaign['name']}: {campaign['points']} puntos ({state_text})")
    else:
        lines.append("  - ninguna campaña")
    lines.extend(["", "Nota: este informe organiza evidencia; no sustituye la revisión científica."])
    wrapped = []
    for line in lines:
        if line.startswith("  "):
            wrapped.extend(textwrap.wrap(line, width=106,
                                          subsequent_indent="      ") or [""])
        else:
            wrapped.extend(textwrap.wrap(line, width=106) or [""])
    return wrapped


def generate_pdf(root, data, destination=None) -> Path:
    """Escribe un PDF autocontenido y devuelve su ruta."""
    target = Path(destination or (Path(root) / project.PROJECT_DIR / "reports" /
                                  "project-report.pdf"))
    if not target.is_absolute():
        target = Path(root) / target
    target.parent.mkdir(parents=True, exist_ok=True)
    try:
        import matplotlib.pyplot as plt
        from matplotlib.backends.backend_pdf import PdfPages
    except ImportError as exc:
        raise RuntimeError(f"no se pudo cargar el generador PDF: {exc}") from None

    lines = _lines(Path(root).resolve(), data)
    with PdfPages(target) as pdf:
        for start in range(0, len(lines), 42):
            page = lines[start:start + 42]
            figure = plt.figure(figsize=(8.27, 11.69))
            figure.text(0.06, 0.96, "\n".join(page), va="top", ha="left",
                        family="DejaVu Sans", fontsize=9, linespacing=1.35)
            figure.text(0.06, 0.025, f"Olla-DFT · página {start // 42 + 1}",
                        color="#566573", fontsize=8)
            pdf.savefig(figure, bbox_inches="tight")
            plt.close(figure)
    return target
