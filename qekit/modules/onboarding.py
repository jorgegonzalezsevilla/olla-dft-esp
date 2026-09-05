# Olla-DFT — command-line toolkit for Quantum ESPRESSO
# Copyright (C) 2026 Jorge Enrique González Sevilla
# SPDX-License-Identifier: AGPL-3.0-or-later
# This program is free software: you can redistribute it and/or modify it
# under the terms of the GNU Affero General Public License as published by the Free
# Software Foundation, either version 3 of the License, or (at your option)
# any later version. See the LICENSE file for details.

"""Asistente de inicio para quien nunca ha usado una CLI científica."""

from __future__ import annotations

import json
from pathlib import Path

from qekit import __command_name__
from qekit.core.errors import ErrorDeUso
from qekit.modules import environment, project, validation


GOALS = (
    ("relax", "relajar posiciones y celda"),
    ("gap", "calcular bandas y band gap"),
    ("dos", "calcular DOS y PDOS"),
    ("phonons", "comprobar estabilidad con fonones"),
    ("optics", "estudiar absorción y propiedades ópticas"),
    ("scf", "obtener la energía electrónica básica"),
)
_TRANSLATION_DIR = Path(__file__).resolve().parent.parent / "data" / "i18n"


def _labels(language="es") -> dict:
    if language not in ("es", "en"):
        raise ErrorDeUso("language debe ser es o en")
    target = _TRANSLATION_DIR / f"onboarding_{language}.json"
    try:
        return json.loads(target.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise ErrorDeUso(f"no se pudo cargar el idioma {language}: {exc}") from None


def _ask(prompt, default="", input_fn=None):
    input_fn = input_fn or input
    suffix = f" [{default}]" if default else ""
    answer = input_fn(f"{prompt}{suffix}: ").strip()
    return answer or default


def _goal_from_answer(answer):
    text = str(answer or "").strip().lower()
    if text.isdigit() and 1 <= int(text) <= len(GOALS):
        return GOALS[int(text) - 1][0]
    aliases = {"relajacion": "relax", "relajación": "relax", "relaxation": "relax",
               "bandas": "gap", "bands": "gap", "band gap": "gap",
               "gap": "gap", "pdos": "dos", "fonones": "phonons",
               "phonons": "phonons", "optica": "optics", "óptica": "optics",
               "optics": "optics", "optical": "optics", "energia": "scf",
               "energía": "scf", "energy": "scf", "scf": "scf"}
    return aliases.get(text, text if text in {x[0] for x in GOALS} else None)


def guide(project_path=".", structure_path=None, goal=None, name=None,
          interactive=True, validate=True, input_fn=None, language="es") -> dict:
    """Inicializa o abre un proyecto y deja un workflow revisable."""
    labels = _labels(language)
    localized_goals = tuple((key, labels[f"goal_{key}"]) for key, _ in GOALS)
    project_path = Path(project_path).expanduser().resolve()
    created = False
    try:
        root, data = project.load(project_path)
    except ErrorDeUso:
        # Solo crear cuando realmente no existe un manifiesto. Un proyecto
        # corrupto debe conservar su error original, no parecer uno nuevo.
        try:
            project.manifest_path(project_path)
        except ErrorDeUso:
            root, data = project.init(project_path, name=name)
            created = True
        else:
            raise

    if structure_path is None and not data.get("sources") and interactive:
        structure_path = _ask(labels["structure_path"],
                               input_fn=input_fn)
    if structure_path:
        project.add_source(root, data, structure_path)
    if not data.get("sources"):
        raise ErrorDeUso("el inicio guiado necesita una estructura; indica --structure.")

    selected_goal = _goal_from_answer(goal)
    if selected_goal is None and interactive:
        print("\n" + labels["question"])
        for index, (_key, label) in enumerate(localized_goals, 1):
            print(f"  {index}) {label}")
        answer = _ask(labels["option"], "1", input_fn=input_fn)
        selected_goal = _goal_from_answer(answer) or "relax"
    if selected_goal is None:
        selected_goal = "scf"

    tasks = project.plan(root, data, selected_goal)
    lock = root / project.PROJECT_DIR / environment.LOCK_NAME
    if not lock.is_file():
        environment.write(root)
        data.setdefault("metadata", {})["environment_lock"] = str(lock.relative_to(root))
    if validate:
        checks = validation.check(root, data)
        fails = sum(item["level"] == "fail" for item in checks)
        data.setdefault("metadata", {})["onboarding"] = {
            "created": created, "goal": selected_goal,
            "validated": True, "fails": fails, "at": project._now(),
        }
    project.save(root, data)
    return {"root": root, "data": data, "created": created,
            "goal": selected_goal, "tasks": tasks,
            "validation": validation.check(root, data) if validate else [],
            "language": language}


def report(result: dict) -> str:
    root, data = result["root"], result["data"]
    labels = _labels(result.get("language", "es"))
    lines = ["--- " + labels["title"] + " ---",
             f"{labels['project']}: {data['name']}",
             f"{labels['folder']}: {root}",
             f"{labels['goal']}: {result['goal']}",
             f"{labels['tasks_prepared']}: {len(result['tasks'])}"]
    checks = result.get("validation", [])
    if checks:
        fails = sum(item["level"] == "fail" for item in checks)
        lines.append(f"{labels['validation']}: {len(checks) - fails} "
                     f"{labels['correct']}, {fails} {labels['failures']}")
    lines += ["", labels["next"],
              f"  {__command_name__} project status --project {root}       # {labels['status']}",
              f"  {__command_name__} project validate --project {root} --advanced  # {labels['validate']}",
              f"  {__command_name__} project run --project {root}       # {labels['simulation']}",
              f"  {__command_name__} project run --project {root} --execute  # {labels['run']}",
              f"  {__command_name__} project dashboard --project {root}   # {labels['web']}"]
    return "\n".join(lines)
