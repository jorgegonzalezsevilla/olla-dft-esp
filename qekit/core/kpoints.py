# Olla-DFT — command-line toolkit for Quantum ESPRESSO
# Copyright (C) 2026 Jorge Enrique González Sevilla
# SPDX-License-Identifier: GPL-3.0-or-later
# This program is free software: you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by the Free
# Software Foundation, either version 3 of the License, or (at your option)
# any later version. See the LICENSE file for details.

"""Mallas de k-points y caminos de alta simetría.

- Malla uniforme centrada en Γ (K_POINTS automatic, sin desplazamiento)
  a partir de un espaciado en el espacio recíproco (equivalente al KSPACING
  de VASP que usa VASPKIT). No es una malla Monkhorst-Pack desplazada: con
  n par la malla de MP no contiene Γ y esta sí.
- K-path de alta simetría con seekpath (convención de Hinuma et al.,
  Comput. Mater. Sci. 128, 140 (2017) — la misma que usa Materials Cloud).
"""

from dataclasses import dataclass, field

import numpy as np
import seekpath
from ase import Atoms

from qekit.core import structure

# Niveles de densidad estilo VASPKIT (espaciado en Å^-1, con factor 2π)
KSPACING_LEVELS = {
    "gamma": None,     # solo el punto Γ
    "coarse": 0.30,
    "medium": 0.20,
    "fine": 0.15,
    "very-fine": 0.10,
}


# Vacío (Å) a partir del cual una dirección se considera no periódica de
# verdad: ahí la banda no dispersa y un segundo k-point solo cuesta tiempo.
VACIO_MINIMO = 8.0


def direcciones_con_vacio(atoms: Atoms, minimo: float = VACIO_MINIMO) -> list:
    """Ejes (0, 1, 2) cuyo hueco entre átomos supera `minimo` angstrom."""
    from qekit.modules.inputgen import hueco_vacio
    ejes = []
    for eje in range(3):
        try:
            _, _, hueco_A = hueco_vacio(atoms, eje)
        except Exception:                                   # noqa: BLE001
            continue
        if hueco_A >= minimo:
            ejes.append(eje)
    return ejes


def kgrid_from_spacing(atoms: Atoms, kspacing: float,
                       respetar_vacio: bool = True) -> tuple:
    """Malla (n1, n2, n3) tal que el espaciado entre k-points sea <= kspacing.

    kspacing en Å^-1 e incluye el factor 2π: n_i = ceil(|b_i| / kspacing),
    donde b_i son los vectores recíprocos con 2π.

    Con `respetar_vacio` (por defecto), una dirección que solo contiene vacío
    recibe un único k-point. Una losa de 30 Å de celda con kspacing 0.20 caía
    en 2 k-points a lo largo del vacío: el doble de trabajo para muestrear una
    banda que no dispersa, porque las réplicas ni se ven.
    """
    if kspacing is None or kspacing <= 0:
        return (1, 1, 1)
    recip = 2.0 * np.pi * atoms.cell.reciprocal().array  # filas = b1, b2, b3
    norms = np.linalg.norm(recip, axis=1)
    grid = np.maximum(1, np.ceil(norms / kspacing).astype(int))
    if respetar_vacio:
        for eje in direcciones_con_vacio(atoms):
            grid[eje] = 1
    return tuple(int(n) for n in grid)


@dataclass
class KPath:
    """Resultado del análisis con seekpath."""

    primitive: Atoms                 # celda primitiva estandarizada (¡usarla en el cálculo!)
    path: list                       # lista de pares (etiqueta_inicial, etiqueta_final)
    point_coords: dict               # etiqueta -> coordenadas fraccionarias
    spacegroup: str = ""
    spacegroup_number: int = 0
    cell_changed: bool = True        # la celda primitiva difiere de la de entrada
    segments: list = field(default_factory=list)


def pretty_label(label: str) -> str:
    """GAMMA → Γ, DELTA_0 → Δ₀... para mostrar en pantalla y en gráficas."""
    greek = {"GAMMA": "Γ", "DELTA": "Δ", "SIGMA": "Σ", "LAMBDA": "Λ"}
    base = label
    sub = ""
    if "_" in label:
        base, sub = label.split("_", 1)
    base = greek.get(base, base)
    return base + sub


def get_kpath(atoms: Atoms, symprec: float = structure.SYMPREC) -> KPath:
    """Obtiene el k-path estándar para la estructura dada.

    IMPORTANTE: las coordenadas de los k-points están referidas a la celda
    primitiva estandarizada que devuelve seekpath, no a la celda original.
    Los inputs de bandas deben generarse con esa celda primitiva.
    """
    res = seekpath.get_path(structure.to_spglib_cell(atoms), symprec=symprec)
    prim = structure.from_spglib_cell(
        (res["primitive_lattice"], res["primitive_positions"], res["primitive_types"])
    )
    same_natoms = len(prim) == len(atoms)
    same_cell = same_natoms and np.allclose(
        prim.cell.array, atoms.cell.array, atol=1e-5
    )
    return KPath(
        primitive=prim,
        path=list(res["path"]),
        point_coords=dict(res["point_coords"]),
        spacegroup=res.get("spacegroup_international", ""),
        spacegroup_number=res.get("spacegroup_number", 0),
        cell_changed=not same_cell,
    )


def kpath_card(kpath: KPath, points_per_segment: int = 20) -> tuple:
    """Construye la tarjeta K_POINTS crystal_b para pw.x.

    Devuelve (texto_de_la_tarjeta, lista_de_etiquetas) donde la lista
    contiene tuplas (índice, etiqueta, coordenadas) para post-proceso.

    Convención de pesos en crystal_b: el peso de cada punto es el número
    de puntos entre ese punto y el siguiente; un peso 0 marca una
    discontinuidad (salto directo al siguiente punto).
    """
    entries = []  # (label, coords, weight)
    path = kpath.path
    for i, (start, end) in enumerate(path):
        entries.append((start, kpath.point_coords[start], points_per_segment))
        is_last = i == len(path) - 1
        discontinuous = (not is_last) and (path[i + 1][0] != end)
        if is_last:
            entries.append((end, kpath.point_coords[end], 1))
        elif discontinuous:
            entries.append((end, kpath.point_coords[end], 0))
        # si es continuo, el punto final es el inicial del siguiente segmento

    lines = ["K_POINTS crystal_b", f"{len(entries)}"]
    labels = []
    for idx, (label, coords, weight) in enumerate(entries):
        x, y, z = coords
        lines.append(f"  {x:12.8f} {y:12.8f} {z:12.8f} {weight:4d}")
        labels.append((idx, pretty_label(label), tuple(coords)))
    return "\n".join(lines) + "\n", labels


def kpath_text(kpath: KPath) -> str:
    """Resumen legible del k-path para mostrar al usuario."""
    lines = [
        "--- Camino de alta simetría (seekpath) ---",
        f"Grupo espacial: {kpath.spacegroup} (N.º {kpath.spacegroup_number})",
        "",
        "Camino:",
    ]
    # compactar: G-X-U | K-G-L ...
    chunks = []
    current = []
    for i, (start, end) in enumerate(kpath.path):
        if not current:
            current = [start, end]
        elif current[-1] == start:
            current.append(end)
        else:
            chunks.append(current)
            current = [start, end]
    if current:
        chunks.append(current)
    for chunk in chunks:
        lines.append("  " + " — ".join(pretty_label(p) for p in chunk))
    lines += ["", "Puntos especiales (coordenadas fraccionarias de la celda primitiva):"]
    seen = []
    for start, end in kpath.path:
        for label in (start, end):
            if label not in seen:
                seen.append(label)
    for label in seen:
        x, y, z = kpath.point_coords[label]
        lines.append(f"  {pretty_label(label):8s} {x:10.6f} {y:10.6f} {z:10.6f}")
    if kpath.cell_changed:
        lines += [
            "",
            "AVISO: el k-path está referido a la celda primitiva estandarizada,",
            "que difiere de la celda de entrada. Usa esa celda primitiva en el",
            "cálculo de bandas (Olla-DFT lo hace automáticamente con 'olla-dft gen').",
        ]
    return "\n".join(lines)
