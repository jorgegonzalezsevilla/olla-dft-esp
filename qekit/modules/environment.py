# Olla-DFT — command-line toolkit for Quantum ESPRESSO
# Copyright (C) 2026 Jorge Enrique González Sevilla
# SPDX-License-Identifier: AGPL-3.0-or-later
# This program is free software: you can redistribute it and/or modify it
# under the terms of the GNU Affero General Public License as published by the Free
# Software Foundation, either version 3 of the License, or (at your option)
# any later version. See the LICENSE file for details.

"""Bloqueo ligero del entorno para reproducibilidad local."""

from __future__ import annotations

import hashlib
import importlib.metadata
import json
import os
import platform
import shutil
import sys
import tempfile
from datetime import datetime, timezone
from pathlib import Path

from qekit import __version__
from qekit.core.errors import ErrorDeUso


LOCK_NAME = "environment.lock.json"
LOCK_VERSION = 1
PACKAGES = ("numpy", "ase", "spglib", "seekpath", "matplotlib", "scipy")
BINARIES = ("pw.x", "ph.x", "q2r.x", "matdyn.x", "bands.x", "dos.x",
            "projwfc.x", "epsilon.x", "pp.x", "mpirun")


def _now():
    return datetime.now(timezone.utc).replace(microsecond=0).isoformat()


def _hash(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def capture() -> dict:
    packages = {}
    for name in PACKAGES:
        try:
            packages[name] = importlib.metadata.version(name)
        except importlib.metadata.PackageNotFoundError:
            packages[name] = None
    binaries = {}
    for name in BINARIES:
        path = shutil.which(name)
        if not path:
            binaries[name] = None
            continue
        resolved = Path(path).resolve()
        try:
            binaries[name] = {"path": str(resolved), "size": resolved.stat().st_size,
                              "sha256": _hash(resolved)}
        except OSError:
            binaries[name] = {"path": str(resolved)}
    return {
        "schema_version": LOCK_VERSION,
        "qekit_version": __version__,
        "python": sys.version,
        "python_executable": str(Path(sys.executable).resolve()),
        "platform": platform.platform(),
        "machine": platform.machine(),
        "packages": packages,
        "binaries": binaries,
        "environment": {key: os.environ.get(key) for key in
                        ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "CUDA_VISIBLE_DEVICES")
                        if os.environ.get(key) is not None},
    }


def path(root, destination=None) -> Path:
    target = Path(destination or (Path(root) / ".qekit" / LOCK_NAME))
    return target if target.is_absolute() else Path(root) / target


def write(root, destination=None) -> Path:
    target = path(root, destination)
    target.parent.mkdir(parents=True, exist_ok=True)
    data = capture()
    data["captured"] = _now()
    fd, name = tempfile.mkstemp(prefix="environment.", suffix=".tmp",
                                 dir=str(target.parent))
    try:
        with os.fdopen(fd, "w", encoding="utf-8", newline="\n") as handle:
            json.dump(data, handle, ensure_ascii=False, indent=2, sort_keys=True)
            handle.write("\n")
        Path(name).replace(target)
    except Exception:
        Path(name).unlink(missing_ok=True)
        raise
    return target


def verify(root, lock=None) -> dict:
    target = path(root, lock)
    if not target.is_file():
        raise ErrorDeUso(f"no existe el bloqueo de entorno '{target}'.")
    try:
        expected = json.loads(target.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise ErrorDeUso(f"no se pudo leer '{target}': {exc}") from None
    current = capture()
    changes = []
    for key in ("qekit_version", "python", "python_executable", "platform",
                "machine", "packages", "binaries", "environment"):
        if expected.get(key) != current.get(key):
            changes.append(key)
    return {"path": str(target.resolve()), "locked_at": expected.get("captured"),
            "ok": not changes, "changed": changes,
            "expected": expected, "current": current}


def report(result: dict) -> str:
    lines = ["--- Entorno reproducible ---", f"Bloqueo: {result['path']}"]
    if result.get("locked_at"):
        lines.append(f"Capturado: {result['locked_at']}")
    if result.get("ok"):
        lines.append("OK: Python, dependencias, binarios y variables coinciden.")
    else:
        lines.append("CAMBIOS: " + ", ".join(result.get("changed", [])))
        lines.append("Vuelve a capturar el entorno si el cambio fue intencional.")
    return "\n".join(lines)
