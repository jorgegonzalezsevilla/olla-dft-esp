# Olla-DFT — command-line toolkit for Quantum ESPRESSO
# Copyright (C) 2026 Jorge Enrique González Sevilla
# SPDX-License-Identifier: GPL-3.0-or-later
# This program is free software: you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by the Free
# Software Foundation, either version 3 of the License, or (at your option)
# any later version. See the LICENSE file for details.

"""Motor local de resultados normalizados y trazables.

Un archivo XML de Quantum ESPRESSO contiene muchos números, pero por sí solo
no dice qué entrada produjo el número ni si cambió desde la última lectura.
Este módulo convierte cada lectura en un registro inmutable dentro de una
SQLite local del proyecto. La identidad incluye la ruta, el hash de todos los
archivos observados y la etiqueta del cálculo: volver a ingerir es idempotente
y una salida modificada crea una nueva evidencia, no pisa la anterior.

La base es deliberadamente pequeña y local. No sustituye la auditoría de
compatibilidad: conserva sus advertencias y evita presentar un resultado no
convergido como publicable.
"""

from __future__ import annotations

import hashlib
import json
import math
import sqlite3
from datetime import datetime, timezone
from pathlib import Path

from qekit import __version__
from qekit.core.errors import ErrorDeUso
from qekit.modules import audit


DB_NAME = "results.sqlite3"
SCHEMA_VERSION = 2

SCHEMA = """
CREATE TABLE IF NOT EXISTS results (
    id              TEXT PRIMARY KEY,
    path            TEXT NOT NULL,
    status          TEXT NOT NULL,
    formula         TEXT,
    natoms          INTEGER,
    calculation     TEXT,
    origin          TEXT,
    converged       INTEGER,
    metrics_json    TEXT NOT NULL,
    provenance_json TEXT NOT NULL,
    tag             TEXT,
    ingested        TEXT NOT NULL,
    qekit_version   TEXT NOT NULL
);
CREATE INDEX IF NOT EXISTS idx_results_formula ON results(formula);
CREATE INDEX IF NOT EXISTS idx_results_calculation ON results(calculation);
CREATE INDEX IF NOT EXISTS idx_results_status ON results(status);
CREATE INDEX IF NOT EXISTS idx_results_ingested ON results(ingested);
CREATE TABLE IF NOT EXISTS schema_meta (
    name            TEXT PRIMARY KEY,
    version         INTEGER NOT NULL,
    migrated        TEXT NOT NULL
);
"""


def _now() -> str:
    return datetime.now(timezone.utc).replace(microsecond=0).isoformat()


def project_db(root) -> Path:
    """Ruta estable de la base de resultados de un proyecto."""
    return Path(root).resolve() / ".qekit" / DB_NAME


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _observed_files(path: Path) -> dict:
    """Devuelve hashes deterministas sin seguir symlinks fuera del cálculo."""
    if path.is_file():
        return {path.name: _sha256(path)}
    if not path.is_dir():
        return {}
    result = {}
    root = path.resolve()
    for child in sorted(path.rglob("*")):
        if not child.is_file() or child.is_symlink():
            continue
        try:
            rel = child.resolve().relative_to(root).as_posix()
        except ValueError:
            continue
        result[rel] = _sha256(child)
    return result


def _contains_xml(path: Path) -> bool:
    """Indica si una ruta parece contener una salida parseable de QE."""
    if path.is_file():
        return path.suffix.lower() == ".xml"
    return path.is_dir() and any(
        child.is_file() and child.suffix.lower() == ".xml"
        for child in path.rglob("*.xml"))


def _finite(value) -> bool:
    try:
        return value is not None and math.isfinite(float(value))
    except (TypeError, ValueError):
        return False


def _formula(result) -> str:
    try:
        from ase import Atoms
        return Atoms(symbols=result.symbols).get_chemical_formula()
    except Exception:  # noqa: BLE001
        return "".join(result.symbols or [])


def _metric(value, unit: str, uncertainty=None) -> dict:
    metric = {"value": float(value) if _finite(value) else None, "unit": unit}
    if uncertainty is not None:
        metric["uncertainty"] = (float(uncertainty)
                                  if _finite(uncertainty) and float(uncertainty) >= 0
                                  else None)
    return metric


def _record(run, tag=None) -> dict:
    path = Path(run.path).resolve()
    files = _observed_files(path)
    result = run.result
    observed = json.dumps(files, ensure_ascii=False, sort_keys=True)
    identity = hashlib.sha256(
        json.dumps({"path": str(path), "files": files, "tag": tag},
                   ensure_ascii=False, sort_keys=True).encode("utf-8")
    ).hexdigest()[:32]

    if result is None:
        return {
            "id": identity, "path": str(path), "status": "invalid",
            "formula": None, "natoms": None, "calculation": None,
            "origin": getattr(run, "origen", "dft"), "converged": None,
            "metrics": {},
            "provenance": {"files": files, "observed": observed,
                           "error": run.error, "qekit_version": __version__},
            "tag": tag,
        }

    natoms = len(result.symbols) or 1
    energy = result.total_energy
    gap = (result.lumo - result.homo
           if _finite(result.homo) and _finite(result.lumo) else None)
    calculation = (result.calculation or "").lower() or None
    if result.converged is False:
        status = "not_converged"
    elif calculation in ("nscf", "bands") and not _finite(energy):
        status = "parsed_no_energy"
    else:
        status = "converged" if result.converged is True else "parsed"
    metrics = {
        "energy_total": _metric(energy, "eV"),
        "energy_per_atom": _metric(energy / natoms if _finite(energy) else None,
                                    "eV/atom"),
        "gap": _metric(gap, "eV"),
        "volume": _metric(result.volume, "angstrom^3"),
        "pressure": _metric(result.pressure, "GPa"),
        "max_force": _metric(result.max_force, "eV/angstrom"),
        "magnetization": _metric(result.total_magnetization, "mu_B/cell"),
        "wall_time": _metric(result.wall_time, "s"),
    }
    provenance = {
        "files": files,
        "observed": observed,
        "fingerprint": [str(value) for value in result.fingerprint],
        "calculation": result.calculation,
        "functional": result.functional,
        "pseudos": dict(result.pseudo_files),
        "qekit_version": __version__,
    }
    return {
        "id": identity, "path": str(path), "status": status,
        "formula": _formula(result), "natoms": len(result.symbols),
        "calculation": calculation,
        "origin": getattr(run, "origen", "dft"),
        "converged": result.converged, "metrics": metrics,
        "provenance": provenance, "tag": tag,
    }


def _connect(db_path) -> sqlite3.Connection:
    target = Path(db_path)
    target.parent.mkdir(parents=True, exist_ok=True)
    connection = sqlite3.connect(str(target))
    connection.execute("PRAGMA journal_mode=WAL")
    connection.executescript(SCHEMA)
    _migrate(connection)
    return connection


def _migrate(connection: sqlite3.Connection) -> None:
    """Migra el índice local y deja constancia de su versión."""
    row = connection.execute(
        "SELECT version FROM schema_meta WHERE name = 'results'"
    ).fetchone()
    current = int(row[0]) if row else 1
    columns = {item[1] for item in connection.execute("PRAGMA table_info(results)")}
    if current < 2 and "review_json" not in columns:
        connection.execute(
            "ALTER TABLE results ADD COLUMN review_json TEXT NOT NULL DEFAULT '{}'"
        )
    connection.execute(
        "INSERT INTO schema_meta(name, version, migrated) VALUES('results', ?, ?) "
        "ON CONFLICT(name) DO UPDATE SET version=excluded.version, "
        "migrated=excluded.migrated",
        (SCHEMA_VERSION, _now()))
    connection.commit()


def ingest(paths, db_path, tag=None) -> dict:
    """Ingiere carpetas/XML y devuelve estadísticas más los IDs nuevos."""
    paths = [Path(path).expanduser() for path in (paths or [])]
    if not paths:
        raise ErrorDeUso("results ingest necesita al menos una carpeta o XML.")
    candidates = [path for path in paths if _contains_xml(path)]
    skipped = len(paths) - len(candidates)
    runs = audit.collect(candidates)
    records = [_record(run, tag=tag) for run in runs]
    connection = _connect(db_path)
    inserted, existing = 0, 0
    ids = []
    try:
        for record in records:
            values = (
                record["id"], record["path"], record["status"],
                record["formula"], record["natoms"], record["calculation"],
                record["origin"],
                None if record["converged"] is None else int(record["converged"]),
                json.dumps(record["metrics"], ensure_ascii=False, sort_keys=True),
                json.dumps(record["provenance"], ensure_ascii=False, sort_keys=True),
                record["tag"], _now(), __version__,
            )
            cursor = connection.execute(
                "INSERT OR IGNORE INTO results "
                "(id,path,status,formula,natoms,calculation,origin,converged,"
                "metrics_json,provenance_json,tag,ingested,qekit_version) "
                "VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?)", values)
            if cursor.rowcount:
                inserted += 1
                ids.append(record["id"])
            else:
                existing += 1
        connection.commit()
    finally:
        connection.close()
    return {"read": len(records), "skipped_without_xml": skipped,
            "inserted": inserted, "existing": existing,
            "invalid": sum(r["status"] == "invalid" for r in records),
            "ids": ids, "db": str(Path(db_path).resolve())}


def ingest_project(root, data, paths=None, tag=None, db_path=None) -> dict:
    """Ingiere y actualiza el manifiesto con el último estado del índice."""
    root = Path(root).resolve()
    if paths is None:
        paths = [root / ".qekit" / "artifacts"]
    database = Path(db_path).expanduser() if db_path else project_db(root)
    result = ingest(paths, database, tag=tag)
    try:
        database_label = str(database.resolve().relative_to(root))
    except ValueError:
        database_label = str(database.resolve())
    data.setdefault("metadata", {})["results"] = {
        "db": database_label,
        "last_ingest": _now(),
        "read": result["read"],
        "inserted": result["inserted"],
        "existing": result["existing"],
        "invalid": result["invalid"],
    }
    from qekit.modules import project
    project.save(root, data)
    return result


def _row(row) -> dict:
    result = dict(row)
    result["metrics"] = json.loads(result.pop("metrics_json"))
    result["provenance"] = json.loads(result.pop("provenance_json"))
    review = result.pop("review_json", "{}")
    try:
        result["review"] = json.loads(review or "{}")
    except (TypeError, json.JSONDecodeError):
        result["review"] = {}
    result["converged"] = (None if result["converged"] is None
                            else bool(result["converged"]))
    return result


def list_results(db_path, formula=None, calculation=None, status=None,
                 limit=100) -> list:
    """Lista resultados con filtros parametrizados."""
    try:
        limit = max(1, min(int(limit), 10000))
    except (TypeError, ValueError):
        raise ErrorDeUso("--limit debe ser un entero positivo.") from None
    if not Path(db_path).exists():
        return []
    clauses, values = [], []
    if formula:
        clauses.append("formula LIKE ?")
        values.append(f"%{formula}%")
    if calculation:
        clauses.append("calculation = ?")
        values.append(str(calculation).lower())
    if status:
        clauses.append("status = ?")
        values.append(status)
    where = " WHERE " + " AND ".join(clauses) if clauses else ""
    connection = sqlite3.connect(str(db_path))
    connection.row_factory = sqlite3.Row
    try:
        rows = connection.execute(
            "SELECT * FROM results" + where +
            " ORDER BY ingested DESC, path LIMIT ?", values + [limit]).fetchall()
        return [_row(row) for row in rows]
    finally:
        connection.close()


def get(db_path, result_id) -> dict:
    if not Path(db_path).exists():
        raise ErrorDeUso(f"no existe el índice de resultados '{db_path}'.")
    connection = sqlite3.connect(str(db_path))
    connection.row_factory = sqlite3.Row
    try:
        row = connection.execute("SELECT * FROM results WHERE id = ?",
                                 (str(result_id),)).fetchone()
    finally:
        connection.close()
    if row is None:
        raise ErrorDeUso(f"no encuentro el resultado '{result_id}'.")
    return _row(row)


def review(db_path, result_id, status, note="") -> dict:
    """Añade una revisión humana sin modificar la evidencia ingerida."""
    allowed = {"unreviewed", "accepted", "rejected"}
    if status not in allowed:
        raise ErrorDeUso("--review-status debe ser unreviewed, accepted o rejected.")
    connection = _connect(db_path)
    try:
        payload = {"status": status, "note": str(note or "").strip(),
                   "at": _now()}
        cursor = connection.execute(
            "UPDATE results SET review_json = ? WHERE id = ?",
            (json.dumps(payload, ensure_ascii=False, sort_keys=True),
             str(result_id)))
        if cursor.rowcount == 0:
            raise ErrorDeUso(f"no encuentro el resultado '{result_id}'.")
        connection.commit()
    finally:
        connection.close()
    return get(db_path, result_id)


def summary(db_path) -> dict:
    if not Path(db_path).exists():
        return {"db": str(Path(db_path)), "count": 0, "by_status": {},
                "results": []}
    connection = sqlite3.connect(str(db_path))
    try:
        count = connection.execute("SELECT COUNT(*) FROM results").fetchone()[0]
        by_status = {row[0]: row[1] for row in connection.execute(
            "SELECT status, COUNT(*) FROM results GROUP BY status")}
    finally:
        connection.close()
    return {"db": str(Path(db_path)), "count": count, "by_status": by_status,
            "results": list_results(db_path, limit=20)}


def report(rows: list, db_path=None) -> str:
    lines = ["--- Resultados normalizados ---"]
    if db_path:
        lines.append(f"Índice: {Path(db_path).resolve()}")
    if not rows:
        lines.append("No hay resultados ingeridos.")
        return "\n".join(lines)
    lines.append(f"Registros mostrados: {len(rows)}")
    header = f"{'id':10s} {'fórmula':12s} {'tipo':10s} {'estado':16s} {'E/átomo':>14s}"
    lines += [header, "-" * len(header)]
    for row in rows:
        energy = row.get("metrics", {}).get("energy_per_atom", {}).get("value")
        energy_text = f"{energy:14.6f}" if _finite(energy) else f"{'-':>14s}"
        lines.append(f"{row['id'][:10]:10s} {(row.get('formula') or '?'):12s} "
                     f"{(row.get('calculation') or '?'):10s} "
                     f"{row['status']:16s} {energy_text}")
    return "\n".join(lines)


def export(db_path, destination) -> Path:
    target = Path(destination)
    target.parent.mkdir(parents=True, exist_ok=True)
    rows = list_results(db_path, limit=10000)
    payload = {"schema_version": SCHEMA_VERSION, "qekit_version": __version__,
               "generated": _now(), "db": str(Path(db_path).resolve()),
               "count": len(rows), "results": rows}
    target.write_text(json.dumps(payload, ensure_ascii=False, indent=2) + "\n",
                      encoding="utf-8")
    return target
