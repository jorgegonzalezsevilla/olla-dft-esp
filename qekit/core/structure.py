# Olla-DFT — command-line toolkit for Quantum ESPRESSO
# Copyright (C) 2026 Jorge Enrique González Sevilla
# SPDX-License-Identifier: GPL-3.0-or-later
# This program is free software: you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by the Free
# Software Foundation, either version 3 of the License, or (at your option)
# any later version. See the LICENSE file for details.

"""Lectura de estructuras y análisis de simetría (ASE + spglib).

Formatos de entrada soportados (autodetectados por ASE):
CIF, POSCAR/CONTCAR (VASP), entradas/salidas de pw.x, XYZ (con celda), etc.
"""

from pathlib import Path

import spglib
from ase import Atoms
from qekit.core.errors import ErrorDeUso


def _ase_read(*args, **kwargs):
    """ase.io se importa al usarlo: cuesta ~0.2 s y no todo comando lee archivos."""
    from ase.io import read
    return read(*args, **kwargs)


def _ase_write(*args, **kwargs):
    from ase.io import write
    return write(*args, **kwargs)

SYMPREC = 1e-4  # tolerancia de simetría (Å)


def load(filename: str) -> Atoms:
    """Lee una estructura desde un archivo y valida que tenga celda."""
    path = Path(filename)
    if not path.exists():
        raise FileNotFoundError(f"no existe el archivo '{filename}'")
    name = path.name.upper()
    try:
        if name.startswith(("POSCAR", "CONTCAR")):
            atoms = _ase_read(path, format="vasp")
        else:
            atoms = _ase_read(path)
    except Exception as exc:                            # noqa: BLE001
        # ASE adivina el formato por el contenido y, si se le pasa un .dat
        # o una tabla, revienta con un error interno que no dice nada al
        # usuario. Aquí sí es un error de uso: el archivo no es una estructura.
        raise ErrorDeUso(
            f"no se pudo leer '{filename}' como estructura "
            f"({type(exc).__name__}: {exc}). Se esperaba un CIF, POSCAR, "
            "XYZ con celda o un input de pw.x.") from None
    if isinstance(atoms, list):
        atoms = atoms[-1]
    if atoms is None or atoms.cell.volume < 1e-6:
        raise ErrorDeUso(f"'{filename}' no contiene una celda unitaria válida")
    return atoms


def to_spglib_cell(atoms: Atoms):
    """Convierte un objeto Atoms de ASE a la tupla que espera spglib."""
    return (
        atoms.cell.array,
        atoms.get_scaled_positions(),
        atoms.get_atomic_numbers(),
    )


def from_spglib_cell(cell) -> Atoms:
    """Convierte la tupla (celda, posiciones, números) de spglib a Atoms."""
    lattice, positions, numbers = cell
    return Atoms(
        numbers=numbers,
        scaled_positions=positions,
        cell=lattice,
        pbc=True,
    )


def symmetry_dataset(atoms: Atoms, symprec: float = SYMPREC):
    """Dataset completo de simetría de spglib (grupo espacial, etc.)."""
    dataset = spglib.get_symmetry_dataset(to_spglib_cell(atoms), symprec=symprec)
    if dataset is None:
        raise RuntimeError("spglib no pudo determinar la simetría de la estructura")
    return dataset


def primitive(atoms: Atoms, symprec: float = SYMPREC) -> Atoms:
    """Celda primitiva estandarizada."""
    cell = spglib.standardize_cell(
        to_spglib_cell(atoms), to_primitive=True, symprec=symprec
    )
    if cell is None:
        raise RuntimeError("spglib no pudo estandarizar la celda")
    return from_spglib_cell(cell)


def conventional(atoms: Atoms, symprec: float = SYMPREC) -> Atoms:
    """Celda convencional estandarizada."""
    cell = spglib.standardize_cell(
        to_spglib_cell(atoms), to_primitive=False, symprec=symprec
    )
    if cell is None:
        raise RuntimeError("spglib no pudo estandarizar la celda")
    return from_spglib_cell(cell)


def supercell(atoms: Atoms, nx: int, ny: int, nz: int) -> Atoms:
    """Supercelda nx × ny × nz."""
    if min(nx, ny, nz) < 1:
        raise ErrorDeUso("los factores de la supercelda deben ser >= 1")
    return atoms.repeat((nx, ny, nz))


def conserva_fijos(outfile: str) -> bool:
    """¿El formato de `outfile` guarda qué átomos están congelados?

    POSCAR/CONTCAR/.vasp lo hacen ("Selective dynamics") y ASE lo lee de
    vuelta como FixAtoms; el CIF y el XYZ lo pierden.
    """
    path = Path(outfile)
    name = path.name.upper()
    return (name.startswith(("POSCAR", "CONTCAR"))
            or path.suffix.lower() == ".vasp")


def convert(atoms: Atoms, outfile: str) -> str:
    """Escribe la estructura en el formato que indique la extensión.

    .cif → CIF | POSCAR/.vasp → VASP | .xyz → XYZ (extendido, con celda)
    Solo POSCAR/.vasp conserva los átomos congelados (ver `conserva_fijos`).
    """
    path = Path(outfile)
    # crear la carpeta si hace falta: escribir a una ruta nueva es lo normal
    # (olla-dft surface ... -o resultados/losa.cif) y fallar ahí es gratuito
    if path.parent and not path.parent.exists():
        path.parent.mkdir(parents=True, exist_ok=True)
    name = path.name.upper()
    if name.startswith(("POSCAR", "CONTCAR")) or path.suffix.lower() == ".vasp":
        _ase_write(path, atoms, format="vasp", direct=True, sort=True)
    else:
        _ase_write(path, atoms)
    return str(path)


def info_text(atoms: Atoms, symprec: float = SYMPREC) -> str:
    """Resumen legible de la estructura: celda, composición y simetría."""
    ds = symmetry_dataset(atoms, symprec)
    a, b, c, alpha, beta, gamma = atoms.cell.cellpar()
    prim = primitive(atoms, symprec)

    symbols = atoms.get_chemical_symbols()
    composition = {}
    for s in symbols:
        composition[s] = composition.get(s, 0) + 1
    comp_str = " ".join(f"{el}{n}" for el, n in sorted(composition.items()))

    lines = [
        "--- Estructura ---",
        f"Fórmula:            {atoms.get_chemical_formula()}   ({comp_str})",
        f"Número de átomos:   {len(atoms)}",
        f"Volumen:            {atoms.cell.volume:.4f} Å³",
        "",
        "Parámetros de red:",
        f"  a = {a:.5f} Å   b = {b:.5f} Å   c = {c:.5f} Å",
        f"  α = {alpha:.3f}°   β = {beta:.3f}°   γ = {gamma:.3f}°",
        "",
        "--- Simetría (spglib) ---",
        f"Grupo espacial:     {ds.international} (N.º {ds.number})",
        f"Símbolo de Hall:    {ds.hall}",
        f"Grupo puntual:      {ds.pointgroup}",
        f"Átomos en la celda primitiva: {len(prim)}",
    ]

    wyckoffs = sorted(set(ds.wyckoffs))
    lines.append(f"Posiciones de Wyckoff: {' '.join(wyckoffs)}")

    lines += [
        "",
        "Vectores de celda (Å):",
    ]
    for vec in atoms.cell.array:
        lines.append(f"  {vec[0]:12.6f} {vec[1]:12.6f} {vec[2]:12.6f}")
    return "\n".join(lines)
