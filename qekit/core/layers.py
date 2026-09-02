# Olla-DFT — command-line toolkit for Quantum ESPRESSO
# Copyright (C) 2026 Jorge Enrique González Sevilla
# SPDX-License-Identifier: GPL-3.0-or-later
# This program is free software: you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by the Free
# Software Foundation, either version 3 of the License, or (at your option)
# any later version. See the LICENSE file for details.

"""Detección de capas en materiales laminares.

El criterio no es geométrico sino de enlace: se construye la red de enlaces
covalentes (suma de radios covalentes más una tolerancia) respetando la
periodicidad, se separan las componentes conexas y se determina en cuántas
direcciones es periódica cada una. Una componente periódica en exactamente
dos direcciones es una capa; en tres, un armazón 3D; en una, una cadena; en
ninguna, una molécula.

La dimensionalidad de una componente se obtiene del rango del retículo de
sus "vectores de cierre": al recorrer la red de enlaces asignando a cada
átomo un desplazamiento de celda respecto a un átomo raíz, cada enlace que
no cuadra con esa asignación aporta un vector entero no nulo; el subespacio
que generan esos vectores es exactamente el conjunto de direcciones en que
la componente se repite.
"""

from dataclasses import dataclass, field

import numpy as np
from ase import Atoms
from ase.data import covalent_radii
from ase.neighborlist import neighbor_list
from qekit.core.errors import FaltanDatos

DEFAULT_TOL = 0.45   # Å añadidos a la suma de radios covalentes


@dataclass
class Layer:
    indices: list                 # átomos de la capa (índices en la celda)
    formula: str = ""
    center: float = 0.0           # posición del centro a lo largo del apilado (Å)
    thickness: float = 0.0        # grosor atómico (Å, sin radios de vdW)


@dataclass
class LayerAnalysis:
    dimensionality: str = ""      # "2D", "3D", "1D", "0D" o mezcla ("2D+0D")
    n_components: int = 0
    components_dim: list = field(default_factory=list)
    layers: list = field(default_factory=list)      # solo componentes 2D
    stacking_axis: int = None     # índice del vector de celda de apilamiento
    normal: np.ndarray = None     # normal unitaria a las capas (cartesiana)
    basal_spacing: float = None   # distancia entre capas equivalentes (Å)
    gap: float = None             # hueco entre superficies atómicas (Å)
    period: float = None          # altura de la celda a lo largo de la normal
    tol: float = DEFAULT_TOL


def bonds(atoms: Atoms, tol: float = DEFAULT_TOL):
    """Lista de enlaces (i, j, S) con S el vector de celda del vecino."""
    numbers = atoms.get_atomic_numbers()
    cutoffs = [covalent_radii[z] + tol / 2.0 for z in numbers]
    i, j, S = neighbor_list("ijS", atoms, cutoffs)
    return i, j, S


def _components_and_dim(natoms: int, bi, bj, bS):
    """Componentes conexas y dimensionalidad de cada una.

    Devuelve (etiqueta_por_átomo, [(átomos, dim, base_de_offsets), ...]).
    """
    from collections import defaultdict, deque

    adj = defaultdict(list)
    for a, b, S in zip(bi, bj, bS):
        adj[a].append((b, tuple(S)))

    label = [-1] * natoms
    comps = []
    for start in range(natoms):
        if label[start] != -1:
            continue
        cid = len(comps)
        offset = {start: np.zeros(3, dtype=int)}
        label[start] = cid
        closure = []
        queue = deque([start])
        members = [start]
        while queue:
            a = queue.popleft()
            for b, S in adj[a]:
                S = np.array(S, dtype=int)
                if b not in offset:
                    offset[b] = offset[a] + S
                    label[b] = cid
                    members.append(b)
                    queue.append(b)
                else:
                    d = offset[a] + S - offset[b]
                    if np.any(d):
                        closure.append(d)
        if closure:
            M = np.array(closure)
            dim = int(np.linalg.matrix_rank(M))
        else:
            M = np.zeros((0, 3), dtype=int)
            dim = 0
        comps.append((members, dim, M))
    return label, comps


def analyze(atoms: Atoms, tol: float = DEFAULT_TOL) -> LayerAnalysis:
    """Análisis completo de capas de una estructura periódica."""
    res = LayerAnalysis(tol=tol)
    bi, bj, bS = bonds(atoms, tol)
    label, comps = _components_and_dim(len(atoms), bi, bj, bS)

    res.n_components = len(comps)
    res.components_dim = [dim for _, dim, _ in comps]
    dims = sorted(set(res.components_dim), reverse=True)
    res.dimensionality = "+".join(f"{d}D" for d in dims) if dims else "0D"

    twod = [(members, M) for members, dim, M in comps if dim == 2]
    if not twod:
        return res

    # --- eje de apilamiento -------------------------------------------------
    # Los vectores de cierre de una capa generan un plano en coordenadas
    # fraccionarias; la dirección fraccionaria fuera de ese plano es la de
    # apilamiento. Se toma el vector de celda con mayor componente fuera
    # del plano.
    M = twod[0][1].astype(float)
    _u, _s, vh = np.linalg.svd(M)
    frac_normal = vh[-1]                      # normal al plano, en frac
    axis = int(np.argmax(np.abs(frac_normal)))
    res.stacking_axis = axis

    # Normal cartesiana: perpendicular a los dos vectores de celda en el plano
    others = [k for k in range(3) if k != axis]
    n = np.cross(atoms.cell.array[others[0]], atoms.cell.array[others[1]])
    n /= np.linalg.norm(n)
    if np.dot(n, atoms.cell.array[axis]) < 0:
        n = -n
    res.normal = n
    res.period = float(abs(np.dot(atoms.cell.array[axis], n)))

    # --- capas ordenadas a lo largo de la normal ---------------------------
    # Para medir posiciones sin cortes de celda, cada capa se reconstruye
    # contigua usando los desplazamientos del recorrido.
    pos = atoms.get_positions()
    layers = []
    for members, _M in twod:
        # reconstrucción contigua: BFS local con offsets cartesianos
        from collections import deque
        adj = {}
        mset = set(members)
        for a, b, S in zip(bi, bj, bS):
            if a in mset and b in mset:
                adj.setdefault(a, []).append((b, S))
        root = members[0]
        fixed = {root: pos[root]}
        queue = deque([root])
        while queue:
            a = queue.popleft()
            for b, S in adj.get(a, []):
                if b not in fixed:
                    fixed[b] = fixed[a] + (pos[b] + S @ atoms.cell.array - pos[a])
                    queue.append(b)
        zs = np.array([np.dot(fixed[m], n) for m in members])
        sub = atoms[members]
        layers.append(Layer(
            indices=list(members),
            formula=sub.get_chemical_formula(),
            center=float(zs.mean()),
            thickness=float(zs.max() - zs.min()),
        ))
    layers.sort(key=lambda L: L.center)
    res.layers = layers

    # --- espaciados --------------------------------------------------------
    nlay = len(layers)
    res.basal_spacing = res.period / nlay if nlay else None
    if nlay >= 1:
        # hueco: de la cara superior de una capa a la inferior de la imagen
        # siguiente (la propia capa desplazada un espaciado si solo hay una)
        tops = [L.center + L.thickness / 2 for L in layers]
        bottoms = [L.center - L.thickness / 2 for L in layers]
        bottoms.append(bottoms[0] + res.period)
        gaps = [bottoms[k + 1] - tops[k] for k in range(nlay)]
        res.gap = float(min(gaps))
    return res


def report(atoms: Atoms, res: LayerAnalysis, wavelength: float = 1.5406,
           radiation: str = None) -> str:
    """Reporte legible; incluye dónde caería la reflexión basal (00l).

    `radiation` es el nombre de la radiación ('Cu Kα', 'Mo Kα'...) que
    acompaña a λ en el rótulo; sin él solo se imprime el valor de λ.
    """
    lines = ["--- Análisis de capas ---",
             f"Fórmula: {atoms.get_chemical_formula()}  |  "
             f"tolerancia de enlace: +{res.tol:g} Å sobre radios covalentes",
             f"Componentes conexas: {res.n_components}  |  "
             f"dimensionalidad: {res.dimensionality}"]
    if not res.layers:
        lines.append("")
        lines.append("No se detectaron capas (ninguna componente es periódica "
                     "en exactamente 2 direcciones).")
        if 3 in res.components_dim:
            lines.append("La estructura es un armazón 3D con esta tolerancia; "
                         "puedes probar con --tol menor.")
        return "\n".join(lines)

    axis_name = "abc"[res.stacking_axis]
    lines += [
        "",
        f"Capas detectadas: {len(res.layers)}  "
        f"(apiladas a lo largo del eje {axis_name})",
    ]
    for k, L in enumerate(res.layers, start=1):
        lines.append(f"  capa {k}: {L.formula:12s} grosor {L.thickness:6.3f} Å  "
                     f"centro en {L.center:8.3f} Å")
    lines += [
        "",
        f"Espaciado basal d = {res.basal_spacing:.4f} Å",
        f"Hueco interlaminar (entre superficies atómicas): {res.gap:.4f} Å",
        f"Periodo de apilamiento: {res.period:.4f} Å "
        f"({len(res.layers)} capa(s) por celda)",
    ]
    # posición del pico basal en un difractograma
    lam = wavelength
    lines.append("")
    # el nombre de la radiación viene de fuera: aquí solo se conoce λ y no
    # hay que suponer que es Cu Kα si el usuario pidió otra
    rad = f", {radiation}" if radiation else ""
    lines.append(f"Reflexiones basales esperadas (λ = {lam:.4f} Å{rad}):")
    for order in (1, 2, 3):
        d = res.basal_spacing / order
        s = lam / (2.0 * d)
        if s >= 1.0:
            break
        tt = np.degrees(2.0 * np.arcsin(s))
        lines.append(f"  d = {d:7.4f} Å  ->  2θ = {tt:6.2f}°")
    lines.append("")
    lines.append("El espaciado basal es la distancia entre capas equivalentes "
                 "(centro a centro);\nel hueco interlaminar resta el grosor "
                 "atómico de la capa, sin radios de van der Waals.")
    return "\n".join(lines)


def make_slab(atoms: Atoms, res: LayerAnalysis, layer_index: int = 0,
              vacuum: float = 20.0) -> Atoms:
    """Monocapa aislada con vacío, para energía de exfoliación o superficies.

    Toma una capa detectada, la reconstruye contigua y la coloca centrada en
    una celda con `vacuum` Å de vacío total a lo largo de la normal.
    """
    if not res.layers:
        raise FaltanDatos("no hay capas detectadas de las que construir la monocapa")
    L = res.layers[layer_index]
    slab = atoms[L.indices]
    axis = res.stacking_axis
    n = res.normal

    # Desenrollar la capa a lo largo del apilamiento ANTES de tocar la celda:
    # si la capa cruza la frontera periódica (parte de los átomos en z ≈ 0 y
    # parte en z ≈ c), con las posiciones envueltas el "grosor" sería casi
    # c entero y el centrado la partiría en dos. Imagen mínima respecto del
    # primer átomo de la capa, en fraccionarias del eje de apilamiento.
    frac = slab.get_scaled_positions(wrap=False)
    f = frac[:, axis]
    frac[:, axis] = f - np.round(f - f[0])
    slab.set_scaled_positions(frac)

    # celda: los dos vectores en el plano se conservan; el tercero se
    # sustituye por la normal con la altura capa+vacío
    height = L.thickness + vacuum
    new_cell = np.array(atoms.cell.array, dtype=float)
    new_cell[axis] = n * height
    slab.set_cell(new_cell, scale_atoms=False)

    # centrar la capa
    z = slab.get_positions() @ n
    shift = (height / 2.0) - (z.max() + z.min()) / 2.0
    slab.translate(n * shift)
    slab.wrap()
    return slab
