# Olla-DFT — command-line toolkit for Quantum ESPRESSO
# Copyright (C) 2026 Jorge Enrique González Sevilla
# SPDX-License-Identifier: GPL-3.0-or-later
# This program is free software: you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by the Free
# Software Foundation, either version 3 of the License, or (at your option)
# any later version. See the LICENSE file for details.

"""Project Hub de Olla-DFT: proyectos, workflows reanudables y procedencia.

El proyecto es un manifiesto JSON pequeño dentro de ``.qekit``. Los archivos
grandes siguen donde el usuario los puso; se guardan ruta relativa, tamaño y
SHA-256 para saber si cambiaron. Las tareas son comandos de Olla-DFT generados o
revisados por el usuario y se ejecutan como listas de argumentos, nunca con
``shell=True``.
"""

from __future__ import annotations

import hashlib
import json
import os
import shlex
import subprocess
import sys
import tempfile
from concurrent.futures import ThreadPoolExecutor
from contextlib import contextmanager
from datetime import datetime, timezone
from pathlib import Path

from qekit import __command_name__, __version__
from qekit.core.errors import ErrorDeUso


SCHEMA_VERSION = 2
LEGACY_SCHEMA_VERSIONS = (1,)
PROJECT_DIR = ".qekit"
MANIFEST_NAME = "project.json"
CANCEL_NAME = "CANCEL.json"


def _now() -> str:
    return datetime.now(timezone.utc).replace(microsecond=0).isoformat()


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as fh:
        for chunk in iter(lambda: fh.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def manifest_path(path=".") -> Path:
    """Encuentra el manifiesto desde un directorio, manifiesto o subcarpeta."""
    raw = Path(path).expanduser()
    if raw.is_file():
        if raw.name != MANIFEST_NAME:
            raise ErrorDeUso(f"'{raw}' no es un manifiesto {MANIFEST_NAME}.")
        return raw.resolve()
    if raw.name == PROJECT_DIR and raw.is_dir():
        candidate = raw / MANIFEST_NAME
        if candidate.is_file():
            return candidate.resolve()
    candidate = raw / PROJECT_DIR / MANIFEST_NAME
    if candidate.is_file():
        return candidate.resolve()
    if raw.is_dir():
        for parent in (raw.resolve(), *raw.resolve().parents):
            candidate = parent / PROJECT_DIR / MANIFEST_NAME
            if candidate.is_file():
                return candidate
    raise ErrorDeUso(
        f"no encuentro un proyecto Olla-DFT desde '{raw}'. "
        "Inicializa uno con 'olla-dft project init'.")


def _migrate(data: dict) -> tuple:
    """Actualiza manifiestos antiguos sin perder tareas ni procedencia."""
    version = data.get("schema_version")
    if version == SCHEMA_VERSION:
        return data, False
    if version not in LEGACY_SCHEMA_VERSIONS:
        raise ErrorDeUso(
            f"se esperaba esquema {SCHEMA_VERSION} (o legado "
            f"{', '.join(map(str, LEGACY_SCHEMA_VERSIONS))})")
    data.setdefault("metadata", {})
    data.setdefault("campaigns", [])
    data.setdefault("sources", [])
    data.setdefault("tasks", [])
    data["metadata"].setdefault("migrations", []).append({
        "from": version, "to": SCHEMA_VERSION, "at": _now(),
    })
    data["schema_version"] = SCHEMA_VERSION
    return data, True


@contextmanager
def _manifest_lock(directory: Path):
    """Serializa escrituras del manifiesto dentro de un mismo host."""
    lock_path = directory / ".project.lock"
    handle = lock_path.open("a+", encoding="utf-8")
    locked = False
    try:
        try:
            import fcntl
            fcntl.flock(handle.fileno(), fcntl.LOCK_EX)
            locked = True
        except (ImportError, OSError):
            # En Windows el replace atómico evita archivos parciales; el
            # lock queda como marcador documentado para futuras herramientas.
            pass
        yield
    finally:
        if locked:
            try:
                fcntl.flock(handle.fileno(), fcntl.LOCK_UN)
            except OSError:
                pass
        handle.close()


def load(path=".") -> tuple:
    file = manifest_path(path)
    try:
        data = json.loads(file.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise ErrorDeUso(f"no se pudo leer el proyecto {file}: {exc}") from None
    if not isinstance(data, dict):
        raise ErrorDeUso(f"manifiesto incompatible en {file}: no es un objeto JSON.")
    try:
        data, migrated = _migrate(data)
    except ErrorDeUso as exc:
        raise ErrorDeUso(f"manifiesto incompatible en {file}: {exc}") from None
    root = file.parent.parent.resolve()
    data.setdefault("sources", [])
    data.setdefault("tasks", [])
    data.setdefault("metadata", {})
    data.setdefault("campaigns", [])
    if migrated:
        save(root, data)
    return root, data


def save(root: Path, data: dict) -> Path:
    root = Path(root).resolve()
    directory = root / PROJECT_DIR
    directory.mkdir(parents=True, exist_ok=True)
    data["updated"] = _now()
    target = directory / MANIFEST_NAME
    # El reemplazo es atómico y el temporal pertenece exclusivamente a Olla-DFT.
    with _manifest_lock(directory):
        fd, name = tempfile.mkstemp(prefix="project.", suffix=".tmp", dir=str(directory))
        try:
            with os.fdopen(fd, "w", encoding="utf-8", newline="\n") as fh:
                json.dump(data, fh, ensure_ascii=False, indent=2, sort_keys=True)
                fh.write("\n")
            Path(name).replace(target)
        except Exception:
            Path(name).unlink(missing_ok=True)
            raise
    return target


def _cancel_path(root: Path, cancel_file=None) -> Path:
    target = Path(cancel_file).expanduser() if cancel_file else (
        Path(root) / PROJECT_DIR / CANCEL_NAME)
    return target if target.is_absolute() else Path(root) / target


def cancel(root: Path, reason="", cancel_file=None) -> Path:
    """Solicita detener tareas nuevas o pendientes de forma cooperativa."""
    target = _cancel_path(Path(root).resolve(), cancel_file)
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_text(json.dumps({"requested": _now(), "reason": str(reason or "")},
                                 ensure_ascii=False, indent=2) + "\n",
                      encoding="utf-8")
    return target


def resume(root: Path, data: dict, cancel_file=None) -> int:
    """Retira la solicitud de cancelación y deja sus tareas listas para reanudar."""
    target = _cancel_path(Path(root).resolve(), cancel_file)
    target.unlink(missing_ok=True)
    changed = 0
    for task in data.get("tasks", []):
        if task.get("status") == "cancelled":
            task["status"] = "pending"
            task["resumed_at"] = _now()
            changed += 1
    data.setdefault("metadata", {})["last_resume"] = _now()
    save(Path(root), data)
    return changed


def init(directory=".", name=None) -> tuple:
    root = Path(directory).expanduser().resolve()
    root.mkdir(parents=True, exist_ok=True)
    target = root / PROJECT_DIR / MANIFEST_NAME
    if target.exists():
        raise ErrorDeUso(f"ya existe el proyecto {target}; no se sobrescribe.")
    project_name = name or root.name or "olla-dft-project"
    data = {
        "schema_version": SCHEMA_VERSION,
        "id": _slug(project_name),
        "name": project_name,
        "created": _now(),
        "updated": _now(),
        "metadata": {"qekit_version": __version__,
                      "manifest_schema": SCHEMA_VERSION},
        "sources": [],
        "tasks": [],
        "campaigns": [],
    }
    save(root, data)
    for directory_name in ("artifacts", "reports", "logs"):
        (root / PROJECT_DIR / directory_name).mkdir(exist_ok=True)
    # Cada proyecto nuevo arranca con una referencia reproducible del entorno.
    # El import local evita que el módulo de captura forme un ciclo de carga.
    try:
        from qekit.modules import environment
        lock = environment.write(root)
        data["metadata"]["environment_lock"] = _relative(root, lock)
        save(root, data)
    except OSError:
        # El proyecto sigue siendo utilizable en un medio muy restringido;
        # `project environment` permite capturarlo más tarde.
        pass
    return root, data


def _slug(text: str) -> str:
    out = "".join(c.lower() if c.isalnum() else "-" for c in str(text))
    return "-".join(x for x in out.split("-") if x) or "olla-dft-project"


def _relative(root: Path, path: Path) -> str:
    try:
        return path.resolve().relative_to(root.resolve()).as_posix()
    except ValueError:
        return str(path.resolve())


def add_source(root: Path, data: dict, source) -> dict:
    path = Path(source).expanduser().resolve()
    if not path.is_file():
        raise ErrorDeUso(f"no existe el archivo de entrada '{source}'.")
    rel = _relative(root, path)
    existing = next((x for x in data["sources"] if x.get("path") == rel), None)
    record = {
        "path": rel,
        "kind": _source_kind(path),
        "size": path.stat().st_size,
        "sha256": sha256_file(path),
        "added": _now(),
    }
    if existing:
        existing.update(record)
        return existing
    data["sources"].append(record)
    return record


def _source_kind(path: Path) -> str:
    suffix = path.suffix.lower()
    if suffix in (".cif", ".vasp", ".poscar", ".xyz", ".xsf"):
        return "structure"
    if suffix in (".in", ".pwi", ".pw"):
        return "qe-input"
    return "file"


def source_path(root: Path, data: dict) -> Path:
    structures = [x for x in data["sources"] if x.get("kind") == "structure"]
    candidates = structures or data["sources"]
    if not candidates:
        raise ErrorDeUso(
            "el proyecto no tiene archivos. Añade una estructura con "
            "'olla-dft project add estructura.cif'.")
    path = Path(candidates[0]["path"])
    return path if path.is_absolute() else root / path


def _task(task_id, label, command, depends=(), outputs=()):
    return {
        "id": task_id,
        "label": label,
        "command": command,
        "depends_on": list(depends),
        "outputs": list(outputs),
        "status": "pending",
        "returncode": None,
        "started": None,
        "ended": None,
        "log": None,
        "output_hashes": {},
    }


def plan(root: Path, data: dict, goal: str, commands=(), task_prefix=None) -> list:
    """Añade un plan idempotente a partir de una intención humana."""
    goal_text = str(goal or "scf").strip().lower()
    explicit = [str(x).strip() for x in (commands or ()) if str(x).strip()]
    tasks = []
    if explicit:
        previous = []
        prefix = _slug(task_prefix or "custom")
        for i, command in enumerate(explicit, 1):
            task_id = f"{prefix}-{i}"
            tasks.append(_task(task_id, f"Tarea personalizada {i}", command,
                               previous[-1:] if previous else []))
            previous.append(task_id)
    else:
        source = _relative(root, source_path(root, data))
        if any(k in goal_text for k in ("dos", "densidad", "pdos")):
            tasks = [
                _task("info", "Revisar estructura", f"olla-dft info {source}"),
                _task("gen-dos", "Preparar SCF + NSCF + DOS",
                      f"olla-dft gen {source} --preset dos --outdir artifacts/dos",
                      ["info"], ["artifacts/dos"]),
                _task("dos", "Analizar DOS/PDOS",
                      "olla-dft dos artifacts/dos --outdir reports/dos",
                      ["gen-dos"], ["reports/dos"]),
            ]
        elif any(k in goal_text for k in ("opt", "epsilon", "absorb", "tauc")):
            tasks = [
                _task("info", "Revisar estructura", f"olla-dft info {source}"),
                _task("gen-optics", "Preparar cálculo óptico",
                      f"olla-dft optics {source} --outdir artifacts/optics",
                      ["info"], ["artifacts/optics"]),
            ]
        elif any(k in goal_text for k in ("phonon", "fonon", "vibr", "raman")):
            tasks = [
                _task("info", "Revisar estructura", f"olla-dft info {source}"),
                _task("phonons", "Preparar fonones DFPT",
                      f"olla-dft phonons {source} --outdir artifacts/phonons",
                      ["info"], ["artifacts/phonons"]),
            ]
        elif any(k in goal_text for k in ("band", "topolog", "chern", "gap")):
            tasks = [
                _task("info", "Revisar estructura", f"olla-dft info {source}"),
                _task("gen-bands", "Preparar bandas",
                      f"olla-dft gen {source} --preset bands --outdir artifacts/bands",
                      ["info"], ["artifacts/bands"]),
                _task("bands", "Analizar bandas y gap",
                      "olla-dft bands artifacts/bands --outdir reports/bands",
                      ["gen-bands"], ["reports/bands"]),
            ]
        elif any(k in goal_text for k in ("relax", "relaj", "optim")):
            tasks = [
                _task("info", "Revisar estructura", f"olla-dft info {source}"),
                _task("gen-relax", "Preparar relajación de posiciones y celda",
                      f"olla-dft gen {source} --preset relax --outdir artifacts/relax",
                      ["info"], ["artifacts/relax"]),
            ]
        else:
            tasks = [
                _task("info", "Revisar estructura", f"olla-dft info {source}"),
                _task("gen-scf", "Preparar SCF",
                      f"olla-dft gen {source} --preset scf --outdir artifacts/scf",
                      ["info"], ["artifacts/scf"]),
            ]

    by_id = {x["id"]: x for x in data["tasks"]}
    new_ids = {x["id"] for x in tasks}
    merged = []
    for task in tasks:
        old = by_id.get(task["id"])
        if old and old.get("command") == task["command"]:
            preserved = dict(task)
            preserved.update(old)
            task = preserved
        merged.append(task)
    data["tasks"] = [x for x in data["tasks"] if x["id"] not in new_ids]
    data["tasks"].extend(merged)
    data["metadata"]["last_goal"] = goal_text
    return tasks


def _toposort(tasks: list) -> list:
    by_id = {t["id"]: t for t in tasks}
    order, temporary, permanent = [], set(), set()

    def visit(task_id):
        if task_id in permanent:
            return
        if task_id in temporary:
            raise ErrorDeUso(f"dependencia circular en la tarea '{task_id}'.")
        if task_id not in by_id:
            raise ErrorDeUso(f"la tarea depende de '{task_id}', que no existe.")
        temporary.add(task_id)
        for dependency in by_id[task_id].get("depends_on", []):
            visit(dependency)
        temporary.remove(task_id)
        permanent.add(task_id)
        order.append(by_id[task_id])

    for task in tasks:
        visit(task["id"])
    return order


def _valid_command(command: str) -> list:
    if any(x in command for x in (";", "|", "&", ">", "<", "`", "$", "\n")):
        raise ErrorDeUso(f"tarea no segura, contiene shell: {command}")
    try:
        tokens = shlex.split(command, posix=(os.name != "nt"))
    except ValueError as exc:
        raise ErrorDeUso(f"tarea ilegible '{command}': {exc}") from None
    if len(tokens) < 2 or tokens[0].lower() != __command_name__:
        raise ErrorDeUso(
            f"cada tarea debe comenzar por '{__command_name__}'.")
    if tokens[1] == "project":
        raise ErrorDeUso(
            "una tarea no puede invocar 'project' recursivamente; "
            "usa comandos científicos concretos.")
    return tokens


def _task_fingerprint(root: Path, data: dict, task: dict) -> str:
    """Huella de entradas relevantes para evitar reutilizar una salida vieja."""
    sources = []
    for source in sorted(data.get("sources", []), key=lambda item: item.get("path", "")):
        path = Path(source.get("path", ""))
        path = path if path.is_absolute() else root / path
        try:
            current = sha256_file(path) if path.is_file() else "missing"
        except OSError:
            current = "unreadable"
        sources.append((source.get("path"), current))
    payload = {"command": task.get("command", ""), "sources": sources,
               "qekit_version": __version__}
    return hashlib.sha256(json.dumps(payload, sort_keys=True,
                                     ensure_ascii=False).encode("utf-8")).hexdigest()


def _outputs_match(root: Path, task: dict) -> bool:
    """Comprueba los hashes guardados sin exigirlos a tareas antiguas."""
    expected = task.get("output_hashes") or {}
    if not expected:
        return True
    for raw, digest in expected.items():
        path = root / raw
        if not path.is_file():
            return False
        try:
            if sha256_file(path) != digest:
                return False
        except OSError:
            return False
    return True


def _execute_task(root: Path, data: dict, task: dict,
                  python_executable=None, retries=0, timeout=None,
                  cancel_file=None) -> tuple:
    """Ejecuta una tarea aislada; no guarda el manifiesto desde el worker."""
    _valid_command(task["command"])
    if _cancel_requested(root, cancel_file):
        task["status"] = "cancelled"
        task["progress"] = 0
        task["cancelled_at"] = _now()
        return task, None, "cancelada por solicitud del usuario"

    for dependency in task.get("depends_on", []):
        dep = next(x for x in data["tasks"] if x["id"] == dependency)
        if dep.get("status") != "succeeded":
            task["status"] = "blocked"
            return task, None, f"bloqueada por {dependency}"

    task["status"] = "running"
    task["progress"] = 0
    task["started"] = _now()
    task["input_fingerprint"] = _task_fingerprint(root, data, task)
    log_rel = Path(PROJECT_DIR) / "logs" / f"{task['id']}.log"
    log_path = root / log_rel
    log_path.parent.mkdir(parents=True, exist_ok=True)
    argv = [python_executable or sys.executable, "-m", "qekit.cli"]
    argv.extend(_valid_command(task["command"])[1:])
    attempts = []
    final_output = ""
    completed = None
    for attempt in range(1, retries + 2):
        attempt_started = _now()
        try:
            completed = subprocess.run(
                argv, cwd=str(root), capture_output=True, text=True,
                encoding="utf-8", errors="replace", check=False,
                timeout=timeout)
            output = completed.stdout + ("\n" + completed.stderr
                                         if completed.stderr else "")
            code = completed.returncode
        except subprocess.TimeoutExpired as exc:
            output = f"timeout después de {timeout:g} s"
            if exc.stdout:
                output += "\n" + str(exc.stdout)
            if exc.stderr:
                output += "\n" + str(exc.stderr)
            code = 124
        except OSError as exc:
            output, code = str(exc), 127
        final_output = output
        attempts.append({"number": attempt, "started": attempt_started,
                         "ended": _now(), "returncode": code,
                         "detail": output[-1000:]})
        if code == 0:
            break
    task["attempts"] = attempts
    task["retry_count"] = max(0, len(attempts) - 1)
    task["returncode"] = attempts[-1]["returncode"]
    task["status"] = ("succeeded" if task["returncode"] == 0 else "failed")
    task["progress"] = 100
    log_path.write_text("\n\n".join(
        f"=== intento {item['number']} · código {item['returncode']} ===\n"
        f"{item['detail']}" for item in attempts), encoding="utf-8")
    task["log"] = log_rel.as_posix()
    task["output_hashes"] = _hash_outputs(root, task.get("outputs", []))
    task["result_ids"] = []
    if task.get("status") == "succeeded":
        # La ingestión es oportunista: un generador que solo escribió
        # pw.in no es todavía un resultado. La acción explícita
        # ``project ingest`` cubre cálculos corridos fuera del workflow.
        try:
            from qekit.modules import results as results_mod
            candidates = []
            for raw in task.get("outputs", []):
                candidate = Path(raw)
                candidate = (candidate if candidate.is_absolute()
                             else root / candidate)
                if candidate.is_file() and candidate.suffix.lower() == ".xml":
                    candidates.append(candidate)
                elif candidate.is_dir() and any(
                        child.is_file() and child.suffix.lower() == ".xml"
                        for child in candidate.rglob("*.xml")):
                    candidates.append(candidate)
            if candidates:
                ingested = results_mod.ingest(
                    candidates, results_mod.project_db(root), tag=task["id"])
                task["result_ids"] = ingested["ids"]
        except Exception as exc:  # noqa: BLE001
            task["result_ingest_warning"] = str(exc)
    result = (task, task["returncode"], final_output[-1000:])
    task["ended"] = _now()
    return result


def run(root: Path, data: dict, execute=False, python_executable=None,
        task_ids=None, force=False, parallel=1, retries=0, timeout=None,
        cancel_file=None) -> list:
    """Ejecuta el DAG, con paralelismo opcional para tareas independientes."""
    try:
        parallel = max(1, int(parallel))
    except (TypeError, ValueError):
        raise ErrorDeUso("parallel debe ser un entero positivo.") from None
    try:
        retries = int(retries)
    except (TypeError, ValueError):
        raise ErrorDeUso("retries debe ser un entero no negativo.") from None
    if retries < 0:
        raise ErrorDeUso("retries debe ser un entero no negativo.")
    if timeout is not None:
        try:
            timeout = float(timeout)
        except (TypeError, ValueError):
            raise ErrorDeUso("timeout debe ser un número positivo.") from None
        if timeout <= 0:
            raise ErrorDeUso("timeout debe ser un número positivo.")
    ordered = _toposort(data["tasks"])
    data.setdefault("metadata", {})["last_run"] = {
        "at": _now(), "execute": bool(execute), "parallel": parallel,
        "retries": retries, "timeout": timeout,
    }
    selected = set(task_ids) if task_ids is not None else None
    selected_tasks = [task for task in ordered
                      if selected is None or task.get("id") in selected]
    results = []
    runnable = []
    for task in selected_tasks:
        if task.get("status") == "running":
            task["status"] = "pending"
            task["recovered_from_running"] = True
        if task.get("status") == "succeeded":
            cached = (not force and (not task.get("input_fingerprint") or
                                     task.get("input_fingerprint") ==
                                     _task_fingerprint(root, data, task)) and
                      _outputs_match(root, task))
            if cached:
                results.append((task, 0, "omitida: caché válida"))
                continue
            task["status"] = "pending"
            task["cache_invalidated"] = True
        if not execute:
            results.append((task, None, "pendiente; usa --execute para correrla"))
        else:
            runnable.append(task)

    if not execute:
        save(root, data)
        return results
    if parallel > 1 and any(task.get("depends_on") for task in runnable):
        raise ErrorDeUso(
            "--parallel solo puede usarse con tareas independientes; "
            "ejecuta el DAG normal o usa una campaña.")

    if parallel > 1 and len(runnable) > 1:
        with ThreadPoolExecutor(max_workers=parallel) as pool:
            results.extend(pool.map(
                lambda task: _execute_task(root, data, task, python_executable,
                                           retries=retries, timeout=timeout,
                                           cancel_file=cancel_file),
                runnable))
        save(root, data)
        return results

    for task in runnable:
        result = _execute_task(root, data, task, python_executable,
                               retries=retries, timeout=timeout,
                               cancel_file=cancel_file)
        results.append(result)
        save(root, data)
        if task.get("status") == "failed":
            break
    return results


def _cancel_requested(root: Path, cancel_file=None) -> bool:
    return _cancel_path(Path(root).resolve(), cancel_file).is_file()


def _hash_outputs(root: Path, outputs) -> dict:
    result = {}
    for raw in outputs:
        path = root / raw
        if path.is_file():
            result[str(raw)] = sha256_file(path)
        elif path.is_dir():
            for child in sorted(path.rglob("*")):
                if child.is_file():
                    result[_relative(root, child)] = sha256_file(child)
    return result


def status(root: Path, data: dict) -> dict:
    counts = {key: 0 for key in ("pending", "running", "succeeded", "failed",
                                 "blocked", "cancelled")}
    for task in data["tasks"]:
        key = task.get("status", "pending")
        counts[key] = counts.get(key, 0) + 1
    total = len(data["tasks"])
    complete = sum(counts.get(key, 0) for key in ("succeeded", "failed", "cancelled"))
    changed = []
    for source in data["sources"]:
        path = Path(source["path"])
        path = path if path.is_absolute() else root / path
        if not path.is_file() or path.stat().st_size != source.get("size") \
                or sha256_file(path) != source.get("sha256"):
            changed.append(source["path"])
    return {"name": data["name"], "root": str(root), "counts": counts,
            "sources": len(data["sources"]), "changed_sources": changed,
            "goal": data.get("metadata", {}).get("last_goal"),
            "progress": int(100 * complete / total) if total else 0}


def report_status(root: Path, data: dict) -> str:
    state = status(root, data)
    c = state["counts"]
    lines = [f"Proyecto: {state['name']}", f"Raíz: {state['root']}",
             f"Fuentes: {state['sources']}",
             f"Progreso: {state['progress']}%",
             "Tareas: " + ", ".join(f"{k}={v}" for k, v in c.items())]
    if state["goal"]:
        lines.append(f"Objetivo: {state['goal']}")
    if state["changed_sources"]:
        lines += ["AVISO — cambiaron fuentes registradas:"]
        lines.extend(f"  - {x}" for x in state["changed_sources"])
    if data["tasks"]:
        lines.append("\nFlujo:")
        for task in _toposort(data["tasks"]):
            lines.append(f"  [{task.get('status', 'pending'):9s}] "
                         f"{task['id']:16s} {task['label']}")
            lines.append(f"             {task['command']}")
    return "\n".join(lines)


def export_snapshot(root: Path, data: dict, destination=None) -> Path:
    target = Path(destination or (root / PROJECT_DIR / "reports" / "provenance.json"))
    if not target.is_absolute():
        target = root / target
    target.parent.mkdir(parents=True, exist_ok=True)
    snapshot = json.loads(json.dumps(data))
    snapshot["exported"] = _now()
    snapshot["root"] = str(root)
    snapshot["source_state"] = status(root, data)
    target.write_text(json.dumps(snapshot, ensure_ascii=False, indent=2) + "\n",
                      encoding="utf-8")
    return target


def diff(left, right) -> dict:
    """Compara dos manifiestos/snapshots sin modificar ninguno."""
    def as_data(value):
        if isinstance(value, dict):
            return value
        path = Path(value)
        if path.is_dir():
            try:
                _root, data = load(path)
                return data
            except ErrorDeUso:
                raise ErrorDeUso(
                    f"'{path}' no contiene un proyecto Olla-DFT; indica una "
                    "carpeta con .qekit/project.json o un snapshot JSON.") from None
        try:
            return json.loads(path.read_text(encoding="utf-8"))
        except (OSError, json.JSONDecodeError) as exc:
            raise ErrorDeUso(f"no se pudo leer el snapshot '{path}': {exc}") from None

    a, b = as_data(left), as_data(right)
    def keyed(items, key):
        return {str(item.get(key)): item for item in items if isinstance(item, dict)}

    a_sources, b_sources = keyed(a.get("sources", []), "path"), keyed(b.get("sources", []), "path")
    a_tasks, b_tasks = keyed(a.get("tasks", []), "id"), keyed(b.get("tasks", []), "id")
    source_changes = []
    for key in sorted(set(a_sources) | set(b_sources)):
        if key not in a_sources:
            source_changes.append({"path": key, "change": "added"})
        elif key not in b_sources:
            source_changes.append({"path": key, "change": "removed"})
        elif a_sources[key].get("sha256") != b_sources[key].get("sha256"):
            source_changes.append({"path": key, "change": "hash_changed"})
    task_changes = []
    for key in sorted(set(a_tasks) | set(b_tasks)):
        if key not in a_tasks:
            task_changes.append({"id": key, "change": "added",
                                 "status": b_tasks[key].get("status")})
        elif key not in b_tasks:
            task_changes.append({"id": key, "change": "removed",
                                 "status": a_tasks[key].get("status")})
        elif (a_tasks[key].get("status") != b_tasks[key].get("status") or
              a_tasks[key].get("command") != b_tasks[key].get("command")):
            task_changes.append({"id": key, "change": "modified",
                                 "from": a_tasks[key].get("status"),
                                 "to": b_tasks[key].get("status")})
    return {
        "left": a.get("name") or a.get("root"),
        "right": b.get("name") or b.get("root"),
        "source_changes": source_changes,
        "task_changes": task_changes,
        "campaigns_left": len(a.get("campaigns", [])),
        "campaigns_right": len(b.get("campaigns", [])),
    }


def diff_report(result: dict) -> str:
    lines = ["--- Diferencia de snapshots ---",
             f"Izquierda: {result['left']}", f"Derecha: {result['right']}"]
    sources, tasks = result["source_changes"], result["task_changes"]
    lines.append(f"Fuentes cambiadas: {len(sources)}  | tareas cambiadas: {len(tasks)}")
    for item in sources:
        lines.append(f"  fuente [{item['change']:12s}] {item['path']}")
    for item in tasks:
        lines.append(f"  tarea  [{item['change']:12s}] {item['id']}"
                     + (f" ({item.get('from')} -> {item.get('to')})"
                        if item["change"] == "modified" else ""))
    lines.append(f"Campañas: {result['campaigns_left']} -> {result['campaigns_right']}")
    return "\n".join(lines)
