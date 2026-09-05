# Olla-DFT — command-line toolkit for Quantum ESPRESSO
# Copyright (C) 2026 Jorge Enrique González Sevilla
# SPDX-License-Identifier: AGPL-3.0-or-later
# This program is free software: you can redistribute it and/or modify it
# under the terms of the GNU Affero General Public License as published by the Free
# Software Foundation, either version 3 of the License, or (at your option)
# any later version. See the LICENSE file for details.

"""Auditoría de consistencia entre cálculos, y base de datos local.

EL ERROR QUE ESTO EVITA
-----------------------
Restar energías totales de dos cálculos de Quantum ESPRESSO solo es válido
si comparten funcional, pseudopotenciales, cutoffs y tratamiento de las
ocupaciones. Si no, la diferencia no significa nada — y NO hay ningún
aviso: QE devuelve dos números perfectamente formados y la resta da un
resultado plausible pero falso.

Eso invalida, en silencio, ecuaciones de estado, energías de exfoliación,
energías de formación de defectos y cascos convexos. Es de los errores más
caros que se pueden cometer con DFT porque no deja rastro.

La auditoría lee un conjunto de cálculos y agrupa por "huella": si sale más
de un grupo, dice exactamente qué parámetro difiere.

SOBRE LA MALLA DE K
-------------------
La malla se trata aparte, como AVISO y no como incompatibilidad. Un bulk y
una losa necesitan mallas distintas por construcción, así que exigir que
coincidan sería falso. Lo que sí se compara es la DENSIDAD de puntos k por
unidad de volumen recíproco, que es la cantidad que debería ser parecida.
"""

import json
import sqlite3
from dataclasses import dataclass
from pathlib import Path

import numpy as np

from qekit import __version__
from qekit.core import provenance, qeout
from qekit.core.errors import ErrorDeUso

ETIQUETAS = {
    "origen": "origen del cálculo (DFT o potencial aprendido)",
    "functional": "funcional",
    "pseudos": "pseudopotenciales",
    "ecutwfc": "ecutwfc (Ry)",
    "ecutrho": "ecutrho (Ry)",
    "smearing": "tipo de smearing",
    "degauss": "degauss (Ry)",
    "occupations": "ocupaciones",
    "nspin": "nspin",
}


@dataclass
class RunInfo:
    path: str = ""
    result: qeout.QEResult = None
    error: str = ""
    origen: str = "dft"          # "dft" o "mlip:<modelo>"
    mlip: dict = None            # marca de procedencia, si la hay

    @property
    def ok(self) -> bool:
        return self.result is not None and not self.error


def _campos(x) -> dict:
    r = x.result if hasattr(x, "result") else x
    origen = getattr(x, "origen", "dft") if hasattr(x, "result") else "dft"
    return {
        "origen": origen,
        "functional": r.functional or "?",
        "pseudos": tuple(sorted(r.pseudo_files.items())),
        "ecutwfc": r.ecutwfc,
        "ecutrho": r.ecutrho,
        "smearing": r.smearing or "",
        "degauss": r.degauss,
        "occupations": r.occupations_kind or "",
        "nspin": r.nspin,
    }


def collect(paths) -> list:
    """Lee todos los cálculos de una lista de carpetas o archivos XML."""
    runs = []
    for raw in paths:
        p = Path(raw)
        candidatos = [p]
        if p.is_dir():
            hijos = [d for d in sorted(p.iterdir()) if d.is_dir()]
            # si la carpeta no tiene XML propio pero sus hijas sí, se
            # audita el conjunto de hijas (el caso de un barrido)
            try:
                qeout.find_xml(str(p))
            except (FileNotFoundError, ValueError):
                if hijos:
                    candidatos = hijos
        for c in candidatos:
            info = RunInfo(path=str(c))
            try:
                from qekit.modules import mlip as _mlip
                marca = _mlip.read_provenance(c)
                if marca:
                    info.mlip = marca
                    info.origen = marca.get("origen", "mlip:?")
            except Exception:                          # noqa: BLE001
                pass
            try:
                info.result = qeout.read_xml(str(c))
            except Exception as exc:            # noqa: BLE001
                info.error = str(exc).splitlines()[0]
            runs.append(info)
    return runs


def kdensity(r: qeout.QEResult):
    """Puntos k por Å⁻³ de zona de Brillouin (comparable entre celdas)."""
    if not r.kgrid or r.volume in (None, 0):
        return None
    n = r.kgrid[0] * r.kgrid[1] * r.kgrid[2]
    v_bz = (2 * np.pi) ** 3 / r.volume
    return n / v_bz


def audit(runs: list) -> dict:
    """Agrupa por huella y describe las diferencias."""
    buenos = [x for x in runs if x.ok]
    grupos = {}
    for x in buenos:
        # el origen entra en la huella: una energia de MLIP y una de DFT
        # jamas deben caer en el mismo grupo
        grupos.setdefault((x.origen,) + x.result.fingerprint, []).append(x)

    difieren = []
    if len(grupos) > 1:
        campos = [_campos(x) for x in buenos]
        for clave in ETIQUETAS:
            valores = {c[clave] for c in campos}
            if len(valores) > 1:
                difieren.append((clave, valores))

    return {
        "runs": runs,
        "buenos": buenos,
        "fallidos": [x for x in runs if not x.ok],
        # Solo tiene sentido exigir convergencia SCF a los cálculos que
        # tienen ciclo SCF. Un nscf o un bands parte de una densidad ya
        # convergida y no vuelve a hacer el ciclo: marcarlo como "no
        # convergió" sería un falso positivo. Además su energía total no es
        # una cantidad utilizable, y eso se avisa aparte.
        "no_convergidos": [x for x in buenos
                           if x.result.converged is False
                           and (x.result.calculation or "").lower()
                           in ("scf", "relax", "vc-relax", "md", "vc-md")],
        "sin_energia": [x for x in buenos
                        if (x.result.calculation or "").lower()
                        in ("nscf", "bands")],
        "grupos": grupos,
        "difieren": difieren,
        "comparables": len(grupos) <= 1,
    }


def report(a: dict) -> str:
    lines = ["--- Auditoría de consistencia ---",
             f"Cálculos leídos: {len(a['runs'])}  |  "
             f"válidos: {len(a['buenos'])}  |  "
             f"ilegibles: {len(a['fallidos'])}"]

    for x in a["fallidos"]:
        lines.append(f"  no se pudo leer {x.path}: {x.error}")

    if a["no_convergidos"]:
        lines += ["", f"NO CONVERGIERON ({len(a['no_convergidos'])}) — sus "
                      "energías no sirven:"]
        for x in a["no_convergidos"]:
            lines.append(f"  {x.path}")

    if a.get("sin_energia"):
        lines += ["", f"Sin energía utilizable ({len(a['sin_energia'])}): "
                      "son nscf o bands, que parten de una",
                  "densidad ya convergida y no dan una energía total "
                  "comparable."]
        for x in a["sin_energia"]:
            lines.append(f"  {x.path}  ({x.result.calculation})")

    lines.append("")
    if not a["buenos"]:
        lines.append("No hay nada que auditar.")
        return "\n".join(lines)

    if a["comparables"]:
        lines.append("COMPARABLES: los cálculos comparten funcional, "
                     "pseudos, cutoffs y ocupaciones.")
        lines.append("Restar sus energías totales es legítimo.")
    else:
        lines.append(f"NO COMPARABLES: hay {len(a['grupos'])} "
                     "configuraciones distintas entre estos cálculos.")
        lines.append("Restar sus energías totales daría un número sin "
                     "significado físico.")
        lines.append("")
        lines.append("Difieren en:")
        for clave, valores in a["difieren"]:
            lines.append(f"  {ETIQUETAS[clave]}:")
            for v in sorted(valores, key=repr):
                if clave == "pseudos":
                    v = ", ".join(f"{k}={n}" for k, n in v) or "(ninguno)"
                lines.append(f"    - {v}")

    # malla k: aviso, no error
    dens = [(x.path, kdensity(x.result)) for x in a["buenos"]]
    dens = [(p, d) for p, d in dens if d is not None]
    if len(dens) > 1:
        vals = [d for _, d in dens]
        if max(vals) / min(vals) > 2.0:
            lines += ["",
                      "AVISO — densidad de puntos k muy dispar "
                      f"(de {min(vals):.1f} a {max(vals):.1f} puntos/Å⁻³):"]
            for p, d in sorted(dens, key=lambda t: t[1]):
                lines.append(f"    {d:8.1f}   {p}")
            lines.append(
                "  No es una incompatibilidad por sí sola —un bulk y una "
                "losa necesitan\n  mallas distintas— pero si estos cálculos "
                "se van a restar entre sí,\n  la malla debería estar "
                "convergida en todos por igual.")
    return "\n".join(lines)


# ----------------------------------------------------------------------
# Base de datos local
# ----------------------------------------------------------------------
ESQUEMA = """
CREATE TABLE IF NOT EXISTS calculos (
    ruta            TEXT PRIMARY KEY,
    origen          TEXT,
    formula         TEXT,
    natoms          INTEGER,
    calculation     TEXT,
    functional      TEXT,
    ecutwfc         REAL,
    ecutrho         REAL,
    kgrid           TEXT,
    kdensity        REAL,
    nspin           INTEGER,
    smearing        TEXT,
    degauss         REAL,
    pseudos         TEXT,
    energia_eV      REAL,
    energia_por_atomo_eV REAL,
    volumen_A3      REAL,
    presion_GPa     REAL,
    fuerza_max      REAL,
    gap_eV          REAL,
    magnetizacion   REAL,
    convergido      INTEGER,
    n_scf           INTEGER,
    nk              INTEGER,
    nbnd            INTEGER,
    n_bfgs          INTEGER,
    wall_s          REAL,
    huella          TEXT,
    qekit_version   TEXT,
    registrado      TEXT
);
"""


def _fila(x: RunInfo) -> dict:
    r = x.result
    from ase import Atoms
    try:
        formula = Atoms(symbols=r.symbols).get_chemical_formula()
    except Exception:                                   # noqa: BLE001
        formula = "".join(r.symbols[:6])
    nat = len(r.symbols) or 1
    gap = None
    if r.homo is not None and r.lumo is not None:
        gap = r.lumo - r.homo
    return {
        "ruta": str(Path(x.path).resolve()),
        "origen": getattr(x, "origen", "dft"),
        "formula": formula,
        "natoms": nat,
        "calculation": r.calculation,
        "functional": r.functional,
        "ecutwfc": r.ecutwfc,
        "ecutrho": r.ecutrho,
        "kgrid": "x".join(map(str, r.kgrid)) if r.kgrid else None,
        "kdensity": kdensity(r),
        "nspin": r.nspin,
        "smearing": r.smearing,
        "degauss": r.degauss,
        "pseudos": json.dumps(r.pseudo_files, sort_keys=True),
        "energia_eV": r.total_energy,
        "energia_por_atomo_eV": (r.total_energy / nat
                                 if r.total_energy is not None else None),
        "volumen_A3": r.volume,
        "presion_GPa": r.pressure,
        "fuerza_max": r.max_force,
        "gap_eV": gap,
        "magnetizacion": r.total_magnetization,
        "convergido": None if r.converged is None else int(r.converged),
        "n_scf": r.n_scf_steps,
        # Puntos k que USÓ pw.x, no los de la malla. La reducción por
        # simetría llega a ser de 48 a 1, y además pw.x y spglib no siempre
        # coinciden (visto en silicio: 12 operaciones frente a 48, 85 puntos
        # frente a 35). Guardar el número real es lo único que hace
        # comparables la calibración del coste y su predicción.
        # getattr y no r.nk: por aquí pasan también resultados de MLIP, que
        # no tienen puntos k ni bandas porque no resuelven Kohn-Sham. Que un
        # relax con MACE no se pueda indexar por eso sería absurdo.
        "nk": getattr(r, "nk", None) or None,
        # Bandas que USÓ pw.x. Si la calibración del coste supone un número
        # y la predicción lee otro del input, el modelo sale sesgado por su
        # cociente y nada lo delata: los dos lados parecen razonables.
        "nbnd": getattr(r, "nbnd", None) or None,
        # Pasos iónicos de una relajación. Sin esto no se puede saber si un
        # relax tardó por ser grande o por haber dado treinta pasos, y el
        # estimador de coste se queda corto justo en los cálculos largos.
        "n_bfgs": getattr(r, "n_bfgs_steps", None) or None,
        "wall_s": r.wall_time,
        "huella": json.dumps([getattr(x, "origen", "dft")]
                             + [str(v) for v in r.fingerprint]),
        "qekit_version": __version__,
        "registrado": provenance.fields()["generado"],
    }


# Columnas añadidas después de la primera versión del esquema. CREATE TABLE
# IF NOT EXISTS no toca una tabla que ya existe, así que una base creada con
# una versión anterior se queda sin ellas y los INSERT fallan enteros.
_COLUMNAS_NUEVAS = {"nk": "INTEGER", "nbnd": "INTEGER",
                    "n_bfgs": "INTEGER"}


def _migrar(con) -> None:
    """Añade a una base antigua las columnas que le falten."""
    try:
        tiene = {f[1] for f in con.execute("PRAGMA table_info(calculos)")}
    except sqlite3.Error:
        return
    for nombre, tipo in _COLUMNAS_NUEVAS.items():
        if nombre not in tiene:
            try:
                con.execute(f"ALTER TABLE calculos ADD COLUMN {nombre} {tipo}")
            except sqlite3.Error:
                pass


def index(runs: list, db_path="olla-dft.db") -> tuple:
    """Registra (o actualiza) los cálculos en la base de datos local."""
    con = sqlite3.connect(str(db_path))
    con.execute(ESQUEMA)
    _migrar(con)
    nuevos = act = 0
    for x in runs:
        if not x.ok:
            continue
        fila = _fila(x)
        existe = con.execute("SELECT 1 FROM calculos WHERE ruta = ?",
                             (fila["ruta"],)).fetchone()
        cols = ", ".join(fila)
        marks = ", ".join("?" for _ in fila)
        con.execute(f"INSERT OR REPLACE INTO calculos ({cols}) "
                    f"VALUES ({marks})", tuple(fila.values()))
        if existe:
            act += 1
        else:
            nuevos += 1
    con.commit()
    con.close()
    return nuevos, act


def query(sql: str, db_path="olla-dft.db") -> list:
    """Consulta libre sobre la base (solo lectura)."""
    if not Path(db_path).exists():
        raise FileNotFoundError(
            f"no existe la base '{db_path}'. Créala con "
            "'olla-dft db carpeta/ --index'.")
    bajo = sql.strip().lower()
    if not bajo.startswith("select"):
        raise ErrorDeUso("solo se admiten consultas SELECT: la base es un "
                         "índice de resultados, no se edita a mano.")
    con = sqlite3.connect(str(db_path))
    con.row_factory = sqlite3.Row
    try:
        filas = [dict(f) for f in con.execute(sql).fetchall()]
    except sqlite3.OperationalError as exc:
        # "no such column: prefix" a secas obliga a ir a leer el esquema.
        # Como el esquema está aquí mismo, se listan las columnas reales.
        cols = columnas(db_path, con=con)
        raise ErrorDeUso(
            f"{exc}.\nColumnas de 'calculos': {', '.join(cols)}"
        ) from None
    finally:
        con.close()
    return filas


def columnas(db_path="olla-dft.db", tabla: str = "calculos", con=None) -> list:
    """Nombres de columna de una tabla de la base."""
    propio = con is None
    if propio:
        con = sqlite3.connect(str(db_path))
    try:
        return [f[1] for f in con.execute(f"PRAGMA table_info({tabla})")]
    except sqlite3.Error:
        return []
    finally:
        if propio:
            con.close()


def search(db_path="olla-dft.db", formula=None, calculation=None,
           gap_min=None, gap_max=None, limit=100) -> list:
    """Búsqueda parametrizada para no obligar al usuario a escribir SQL."""
    if not Path(db_path).exists():
        raise FileNotFoundError(
            f"no existe la base '{db_path}'. Créala con 'olla-dft db carpeta/'.")
    try:
        limit = max(1, min(int(limit), 10000))
    except (TypeError, ValueError):
        raise ErrorDeUso("--limit debe ser un entero positivo.") from None
    try:
        if gap_min is not None:
            gap_min = float(gap_min)
        if gap_max is not None:
            gap_max = float(gap_max)
    except (TypeError, ValueError):
        raise ErrorDeUso("los límites de gap deben ser numéricos.") from None
    clauses, values = [], []
    if formula:
        clauses.append("formula LIKE ?")
        values.append(f"%{formula}%")
    if calculation:
        clauses.append("calculation = ?")
        values.append(calculation)
    if gap_min is not None:
        clauses.append("gap_eV >= ?")
        values.append(gap_min)
    if gap_max is not None:
        clauses.append("gap_eV <= ?")
        values.append(gap_max)
    where = " WHERE " + " AND ".join(clauses) if clauses else ""
    con = sqlite3.connect(str(db_path))
    con.row_factory = sqlite3.Row
    try:
        rows = con.execute(
            "SELECT formula, calculation, functional, gap_eV, "
            "energia_por_atomo_eV, convergido, ruta FROM calculos" + where +
            " ORDER BY formula, gap_eV LIMIT ?", values + [limit]).fetchall()
        return [dict(row) for row in rows]
    finally:
        con.close()


def summary(db_path="olla-dft.db") -> str:
    filas = query("SELECT formula, calculation, functional, ecutwfc, kgrid, "
                  "energia_por_atomo_eV, convergido, ruta FROM calculos "
                  "ORDER BY formula, energia_por_atomo_eV", db_path)
    if not filas:
        return "La base está vacía."
    lines = [f"--- Base de cálculos ({len(filas)} registros) ---",
             f"{'fórmula':>10s} {'tipo':>10s} {'func':>6s} {'ecut':>6s} "
             f"{'malla':>9s} {'E/átomo(eV)':>13s} {'conv':>5s}"]
    for f in filas:
        conv = {1: "sí", 0: "NO"}.get(f["convergido"], "?")
        e = f["energia_por_atomo_eV"]
        lines.append(
            f"{(f['formula'] or '?'):>10s} {(f['calculation'] or '?'):>10s} "
            f"{(f['functional'] or '?'):>6s} "
            f"{(f['ecutwfc'] or 0):>6.0f} {(f['kgrid'] or '?'):>9s} "
            f"{(e if e is not None else float('nan')):>13.4f} {conv:>5s}")
    return "\n".join(lines)


def export_json(db_path="olla-dft.db", out="calculos.json") -> str:
    filas = query("SELECT * FROM calculos", db_path)
    doc = {"qekit_version": __version__,
           "generado": provenance.fields()["generado"],
           "n": len(filas), "calculos": filas}
    Path(out).write_text(json.dumps(doc, ensure_ascii=False, indent=2))
    return str(out)
