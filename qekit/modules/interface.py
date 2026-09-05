# Olla-DFT — command-line toolkit for Quantum ESPRESSO
# Copyright (C) 2026 Jorge Enrique González Sevilla
# SPDX-License-Identifier: AGPL-3.0-or-later
# This program is free software: you can redistribute it and/or modify it
# under the terms of the GNU Affero General Public License as published by the Free
# Software Foundation, either version 3 of the License, or (at your option)
# any later version. See the LICENSE file for details.

"""Heteroestructuras: apilar dos materiales con la menor deformación posible.

EL PROBLEMA
-----------
Poner grafeno sobre hBN, o un óxido sobre un sustrato, exige que las dos
redes quepan en la MISMA celda periódica. Casi nunca encajan: hay que
buscar una supercelda de cada una cuyos vectores en el plano se parezcan
lo bastante, y estirar un poco la más pequeña.

Ese "un poco" es el número que decide si el cálculo significa algo. Una
deformación del 1 % cambia poco; una del 8 % está cambiando el material
antes de empezar. Este módulo busca la mejor coincidencia, la cuantifica
y la reporta — no la esconde.

CÓMO BUSCA
----------
Sobre los vectores en el plano de cada material se prueban todas las
combinaciones enteras

    A' = m11*a1 + m12*a2      B' = n11*b1 + n12*b2
    A'' = m21*a1 + m22*a2     B'' = n21*b1 + n22*b2

con |mij|, |nij| <= max_index, y se queda con las que dan celdas
parecidas en longitud y ángulo. La deformación se mide sobre la matriz
que lleva la red B a la red A:

    epsilon = B'^-1 A' - I

y se reporta el MAYOR |ε_ij| de la matriz, no un promedio: una
deformación de 0 % en una dirección y 6 % en la otra no es "3 %".

DECISIONES QUE HAY QUE TOMAR A MANO
-----------------------------------
- **Quién se deforma.** Por omisión el material 2 (el de arriba), porque
  lo normal es que el 1 sea el sustrato. Con `--strain both` se reparte
  a medias, ponderado por el área.
- **La distancia entre capas.** No se calcula aquí: hay que relajarla, o
  barrerla. Olla-DFT pone una distancia inicial razonable a partir de los
  radios de van der Waals y avisa de que es un punto de partida.
- **El registro (el desplazamiento lateral).** Dos capas pueden apilarse
  de varias formas no equivalentes (AA, AB...). El módulo genera el
  apilamiento sin desplazar y ofrece desplazamientos; cuál es el estable
  hay que calcularlo.
"""

from dataclasses import dataclass, field
from itertools import product

import numpy as np

from qekit.core import provenance
from qekit.core.errors import ErrorDeUso

#: Radios de van der Waals (Å) para estimar la separación inicial.
R_VDW = {
    "H": 1.20, "He": 1.40, "Li": 1.82, "Be": 1.53, "B": 1.92, "C": 1.70,
    "N": 1.55, "O": 1.52, "F": 1.47, "Ne": 1.54, "Na": 2.27, "Mg": 1.73,
    "Al": 1.84, "Si": 2.10, "P": 1.80, "S": 1.80, "Cl": 1.75, "Ar": 1.88,
    "K": 2.75, "Ca": 2.31, "Ni": 1.63, "Cu": 1.40, "Zn": 1.39, "Ga": 1.87,
    "Ge": 2.11, "As": 1.85, "Se": 1.90, "Br": 1.85, "Kr": 2.02, "Mo": 2.10,
    "Ag": 1.72, "Cd": 1.58, "In": 1.93, "Sn": 2.17, "Sb": 2.06, "Te": 2.06,
    "I": 1.98, "W": 2.10, "Au": 1.66, "Pb": 2.02,
}
R_VDW_DEFECTO = 2.0


@dataclass
class Coincidencia:
    """Una supercelda común candidata."""
    M: np.ndarray = None          # 2x2 enteros para el material 1
    N: np.ndarray = None          # 2x2 enteros para el material 2
    celda: np.ndarray = None      # 2x2 vectores en el plano (Å) resultante
    deformacion: np.ndarray = None    # 2x2, epsilon
    eps_max: float = None         # componente mayor en valor absoluto
    natoms: int = 0
    area: float = None
    n1: int = 0                   # celdas del material 1
    n2: int = 0

    @property
    def eps_pct(self) -> float:
        return 100.0 * self.eps_max


@dataclass
class Heteroestructura:
    atoms: object = None
    coincidencia: Coincidencia = None
    separacion: float = None
    vacio: float = None
    formula1: str = ""
    formula2: str = ""
    candidatas: list = field(default_factory=list)
    avisos: list = field(default_factory=list)


# ----------------------------------------------------------------------
# Búsqueda de la supercelda común
# ----------------------------------------------------------------------
def _plano(atoms) -> np.ndarray:
    """Los dos vectores de red en el plano (se asume que c es la normal)."""
    return np.array(atoms.get_cell())[:2, :2]


def _deformacion(A: np.ndarray, B: np.ndarray) -> np.ndarray:
    """Matriz de deformación que lleva B a A: epsilon = B^-1 A - I."""
    try:
        return np.linalg.solve(B, A) - np.eye(2)
    except np.linalg.LinAlgError:
        return np.full((2, 2), np.inf)


def _celdas_candidatas(vectores: np.ndarray, max_index: int) -> dict:
    """Todas las superceldas enteras, agrupadas por su determinante.

    Agrupar por determinante es lo que hace la búsqueda viable: dos redes
    solo pueden coincidir si sus áreas coinciden, y el área es
    det(M) * area_de_la_celda. Sin esto hay que probar (2n+1)^8
    combinaciones — con max_index=3 son 5.7 millones y tarda medio minuto.
    """
    rango = range(-max_index, max_index + 1)
    por_det = {}
    for m in product(rango, repeat=4):
        M = np.array(m, dtype=int).reshape(2, 2)
        det = int(round(np.linalg.det(M)))
        if det <= 0:
            continue
        por_det.setdefault(det, []).append((M, M @ vectores))
    return por_det


def reducir_2d(A: np.ndarray) -> np.ndarray:
    """Base reducida de Lagrange-Gauss de una red 2D.

    Hace falta para no reportar la misma red seis veces. Una celda de
    (2.46, 4.26 Å) con 30 grados y otra de (2.46, 2.46 Å) con 60 son la
    MISMA red hexagonal escrita con otra base; sin reducir, la búsqueda
    devuelve todas como si fueran alternativas distintas.
    """
    v = np.array(A, dtype=float).copy()
    for _ in range(50):
        if np.dot(v[0], v[0]) > np.dot(v[1], v[1]):
            v = v[::-1].copy()
        n = np.dot(v[0], v[0])
        if n < 1e-12:
            break
        m = int(round(np.dot(v[0], v[1]) / n))
        if m == 0:
            break
        v[1] = v[1] - m * v[0]
    return v


def _forma(A: np.ndarray) -> tuple:
    """(|a1|, |a2|, ángulo) de la red REDUCIDA: única salvo simetría."""
    v = reducir_2d(A)
    l1 = float(np.linalg.norm(v[0]))
    l2 = float(np.linalg.norm(v[1]))
    cos = float(np.dot(v[0], v[1]) / max(l1 * l2, 1e-12))
    corto, largo = sorted((l1, l2))
    return (round(corto, 3), round(largo, 3), round(abs(cos), 3))


def _simplicidad(M: np.ndarray) -> tuple:
    """Cuánto de 'fea' es una transformación entera; menor es mejor.

    Entre las muchas M que dan la misma supercelda, la que se reporta
    debería ser la que un humano escribiría: [[1,0],[0,1]] antes que
    [[-3,-2],[-1,-1]].
    """
    M = np.asarray(M)
    return (int(np.abs(M).sum()), int(np.abs(M).max()),
            int((M < 0).sum()), int(np.count_nonzero(M)))


def buscar(atoms1, atoms2, max_index: int = 4, tol: float = 0.05,
           max_atoms: int = 200, n_mejores: int = 10) -> list:
    """Superceldas comunes ordenadas por deformación y luego por tamaño.

    Se ordena primero por deformación porque una celda pequeña con 7 % de
    deformación no es un buen resultado: es un material distinto.
    """
    a = _plano(atoms1)
    b = _plano(atoms2)
    if abs(np.linalg.det(a)) < 1e-8 or abs(np.linalg.det(b)) < 1e-8:
        raise ErrorDeUso(
            "alguna de las dos celdas es degenerada en el plano ab. Las "
            "estructuras tienen que ser losas con el vacío a lo largo de c "
            "('olla-dft surface' o 'olla-dft layers --slab' las dejan así).")

    n_at1, n_at2 = len(atoms1), len(atoms2)
    area1 = abs(np.linalg.det(a))
    area2 = abs(np.linalg.det(b))
    cand1 = _celdas_candidatas(a, max_index)
    cand2 = _celdas_candidatas(b, max_index)

    mejores = {}
    for d1, lista1 in cand1.items():
        if d1 * n_at1 >= max_atoms:
            continue
        for d2, lista2 in cand2.items():
            if d1 * n_at1 + d2 * n_at2 > max_atoms:
                continue
            # las áreas tienen que coincidir dentro de la tolerancia
            if abs(d1 * area1 - d2 * area2) > 2 * tol * d1 * area1:
                continue
            for M, A in lista1:
                for N, B in lista2:
                    eps = _deformacion(A, B)
                    emax = float(np.abs(eps).max())
                    if emax > tol:
                        continue
                    clave = (d1, d2, _forma(A), round(emax, 4))
                    previo = mejores.get(clave)
                    puntos = _simplicidad(M) + _simplicidad(N)
                    if previo is None or puntos < previo[0]:
                        mejores[clave] = (puntos, Coincidencia(
                            M=M, N=N, celda=A, deformacion=eps, eps_max=emax,
                            natoms=d1 * n_at1 + d2 * n_at2,
                            area=abs(np.linalg.det(A)), n1=d1, n2=d2))

    fuera = [c for _, c in mejores.values()]
    fuera.sort(key=lambda c: (round(c.eps_max, 4), c.natoms,
                              _simplicidad(c.M) + _simplicidad(c.N)))
    return fuera[:n_mejores]


# ----------------------------------------------------------------------
# Construcción
# ----------------------------------------------------------------------
def separacion_vdw(atoms1, atoms2) -> float:
    """Separación inicial entre las capas, de los radios de van der Waals.

    Es un PUNTO DE PARTIDA para relajar, no un resultado: la distancia de
    equilibrio depende del enlace, y con funcionales sin dispersión sale
    demasiado grande.
    """
    s1 = atoms1.get_chemical_symbols()
    s2 = atoms2.get_chemical_symbols()
    z1 = atoms1.get_positions()[:, 2]
    z2 = atoms2.get_positions()[:, 2]
    arriba = [s for s, z in zip(s1, z1) if z > z1.max() - 0.5]
    abajo = [s for s, z in zip(s2, z2) if z < z2.min() + 0.5]
    r1 = max(R_VDW.get(s, R_VDW_DEFECTO) for s in (arriba or s1))
    r2 = max(R_VDW.get(s, R_VDW_DEFECTO) for s in (abajo or s2))
    return 0.85 * (r1 + r2)


def construir(atoms1, atoms2, coincidencia: Coincidencia,
              separacion: float = None, vacio: float = 20.0,
              deformar: str = "second", desplazamiento=(0.0, 0.0)) -> object:
    """Arma la heteroestructura a partir de una coincidencia."""
    from ase import Atoms

    a = _plano(atoms1)
    b = _plano(atoms2)
    A = coincidencia.M @ a
    B = coincidencia.N @ b

    if deformar == "second":
        objetivo1, objetivo2 = A, A
    elif deformar == "first":
        objetivo1, objetivo2 = B, B
    elif deformar == "both":
        # promedio ponderado por el área, que es lo que hace que las dos
        # deformaciones queden del mismo orden
        w = coincidencia.n1 * abs(np.linalg.det(a))
        v = coincidencia.n2 * abs(np.linalg.det(b))
        objetivo1 = objetivo2 = (w * A + v * B) / (w + v)
    else:
        raise ErrorDeUso(
            f"--strain '{deformar}' desconocido. Opciones: first, second, "
            "both.")

    s1 = _supercelda_deformada(atoms1, coincidencia.M, objetivo1)
    s2 = _supercelda_deformada(atoms2, coincidencia.N, objetivo2)

    if separacion is None:
        separacion = separacion_vdw(atoms1, atoms2)

    z1 = s1.get_positions()[:, 2]
    z2 = s2.get_positions()[:, 2]
    desplaza_z = (z1.max() - z2.min()) + separacion
    pos2 = s2.get_positions()
    pos2[:, 2] += desplaza_z
    # desplazamiento lateral en coordenadas de la celda común
    plano = np.array(s1.get_cell())[:2, :2]
    pos2[:, :2] += desplazamiento[0] * plano[0] + desplazamiento[1] * plano[1]

    simbolos = list(s1.get_chemical_symbols()) + list(s2.get_chemical_symbols())
    posiciones = np.vstack([s1.get_positions(), pos2])
    altura = posiciones[:, 2].max() - posiciones[:, 2].min()
    celda = np.array(s1.get_cell())
    celda[2] = [0.0, 0.0, altura + vacio]
    fuera = Atoms(symbols=simbolos, positions=posiciones, cell=celda,
                  pbc=[True, True, True])
    fuera.center(axis=2)
    return fuera


def _supercelda_deformada(atoms, M: np.ndarray, objetivo: np.ndarray):
    """Supercelda M del material, con sus vectores del plano llevados a
    `objetivo` y las posiciones arrastradas con la celda."""
    from ase.build import make_supercell

    P = np.eye(3)
    P[:2, :2] = M
    sc = make_supercell(atoms, P)
    celda = np.array(sc.get_cell())
    frac = sc.get_scaled_positions()
    celda[:2, :2] = objetivo
    sc.set_cell(celda)
    sc.set_scaled_positions(frac)
    return sc


def emparejar(atoms1, atoms2, max_index: int = 4, tol: float = 0.05,
              max_atoms: int = 200, indice: int = 0,
              separacion: float = None, vacio: float = 20.0,
              deformar: str = "second",
              desplazamiento=(0.0, 0.0)) -> Heteroestructura:
    """Busca, elige y construye, en un solo paso."""
    candidatas = buscar(atoms1, atoms2, max_index=max_index, tol=tol,
                        max_atoms=max_atoms)
    if not candidatas:
        raise ErrorDeUso(
            f"no hay ninguna supercelda común con menos de {tol * 100:.0f} % "
            f"de deformación y menos de {max_atoms} átomos.\n"
            "Opciones: subir --tol (y aceptar más deformación), subir "
            "--max-atoms, o subir --max-index para buscar celdas más "
            "giradas.")
    if indice >= len(candidatas):
        raise ErrorDeUso(
            f"solo hay {len(candidatas)} candidatas; pediste la {indice}.")
    c = candidatas[indice]
    sep = separacion if separacion is not None else separacion_vdw(atoms1,
                                                                   atoms2)
    het = Heteroestructura(
        atoms=construir(atoms1, atoms2, c, separacion=sep, vacio=vacio,
                        deformar=deformar, desplazamiento=desplazamiento),
        coincidencia=c, separacion=sep, vacio=vacio,
        formula1=atoms1.get_chemical_formula(),
        formula2=atoms2.get_chemical_formula(),
        candidatas=candidatas)

    if c.eps_pct > 3.0:
        het.avisos.append(
            f"La deformación es del {c.eps_pct:.1f} %. Por encima de ~3 % "
            "no se está\nmodelando el material sino una versión estirada "
            "de él: las bandas, el gap\ny las constantes elásticas cambian. "
            "Busca otra coincidencia (--index) o\nacepta una celda más "
            "grande (--max-atoms).")
    if separacion is None:
        het.avisos.append(
            f"La separación de {sep:.2f} Å sale de los radios de van der "
            "Waals: es un punto\nde partida para relajar, NO un resultado. "
            "Y con un funcional sin corrección\nde dispersión la distancia "
            "de equilibrio saldrá demasiado grande — usa\n'olla-dft gen "
            "--vdw grimme-d3' o equivalente.")
    het.avisos.append(
        "El REGISTRO (cómo se alinean lateralmente las dos capas) no está "
        "optimizado.\nDos apilamientos distintos pueden diferir en decenas "
        "de meV por átomo; para\nsaber cuál es el estable hay que barrer "
        "--shift y comparar energías.")
    return het


# ----------------------------------------------------------------------
# Reporte
# ----------------------------------------------------------------------
def report(het: Heteroestructura, n_candidatas: int = 6) -> str:
    c = het.coincidencia
    lines = ["--- Heteroestructura ---",
             f"Material 1 (abajo): {het.formula1}",
             f"Material 2 (arriba): {het.formula2}",
             "",
             f"Supercelda elegida: {c.n1} celda(s) del 1 y {c.n2} del 2  "
             f"->  {c.natoms} átomos",
             f"Área en el plano: {c.area:.2f} Å²",
             f"Separación inicial: {het.separacion:.2f} Å   "
             f"vacío: {het.vacio:.1f} Å",
             "",
             "Deformación (matriz epsilon, en %):"]
    for fila in c.deformacion:
        lines.append("   " + "  ".join(f"{100 * x:+7.3f}" for x in fila))
    lines.append(f"  mayor componente: {c.eps_pct:.2f} %")
    lines += ["",
              "Transformaciones enteras:",
              f"  material 1:  M = [[{c.M[0,0]:2d} {c.M[0,1]:2d}] "
              f"[{c.M[1,0]:2d} {c.M[1,1]:2d}]]",
              f"  material 2:  N = [[{c.N[0,0]:2d} {c.N[0,1]:2d}] "
              f"[{c.N[1,0]:2d} {c.N[1,1]:2d}]]"]

    if len(het.candidatas) > 1:
        lines += ["", f"Otras candidatas (de {len(het.candidatas)}):",
                  f"  {'#':>2s} {'átomos':>7s} {'deformación':>12s} "
                  f"{'n1':>4s} {'n2':>4s}"]
        for i, cc in enumerate(het.candidatas[:n_candidatas]):
            marca = " <-" if cc is c else ""
            lines.append(f"  {i:2d} {cc.natoms:7d} {cc.eps_pct:11.2f} % "
                         f"{cc.n1:4d} {cc.n2:4d}{marca}")
        lines.append("  Se eligen con --index.")

    for a in het.avisos:
        lines += ["", a]
    return "\n".join(lines)


def export(het: Heteroestructura, outdir: str = ".",
           nombre: str = "heteroestructura") -> list:
    from pathlib import Path
    from ase.io import write

    out = Path(outdir); out.mkdir(parents=True, exist_ok=True)
    cif = out / f"{nombre}.cif"
    write(str(cif), het.atoms)
    txt = out / f"{nombre}.txt"
    txt.write_text(provenance.header_plain(
        "heteroestructura",
        {"deformacion_pct": round(het.coincidencia.eps_pct, 3),
         "separacion_A": round(het.separacion, 3),
         "natoms": het.coincidencia.natoms},
        titulo="Emparejamiento de redes") + "\n" + report(het) + "\n")
    return [str(cif), str(txt)]
