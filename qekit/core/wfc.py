# Olla-DFT — command-line toolkit for Quantum ESPRESSO
# Copyright (C) 2026 Jorge Enrique González Sevilla
# SPDX-License-Identifier: GPL-3.0-or-later
# This program is free software: you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by the Free
# Software Foundation, either version 3 of the License, or (at your option)
# any later version. See the LICENSE file for details.

"""Lector de los archivos binarios de función de onda de Quantum ESPRESSO.

QE 6.x escribe cada punto k en `<outdir>/<prefix>.save/wfc<N>.dat`, en
formato Fortran sin formato (secuencial). No es un formato documentado
para el usuario, pero es simple y estable: cada registro va envuelto
entre dos marcas de 4 bytes con su longitud.

    registro 1:  ik(i4)  xk(3, r8)  ispin(i4)  gamma_only(l4)  scalef(r8)
    registro 2:  ngw(i4)  igwx(i4)  npol(i4)  nbnd(i4)
    registro 3:  b1(3,r8) b2(3,r8) b3(3,r8)
    registro 4:  mill(3, igwx) (i4)
    registros 5..: un registro por banda, con npol*igwx complejos (c16)

`xk` y los `b` vienen en unidades de 2*pi/alat, y los índices de Miller
son enteros sobre esa base recíproca. Esto es lo único que hacía falta
para desdoblar bandas, y es la razón por la que ese módulo no existía.
"""

from dataclasses import dataclass
from pathlib import Path

import numpy as np


@dataclass
class WfcK:
    """Las funciones de onda de un punto k."""
    ik: int = 0
    xk: np.ndarray = None            # (3,) en 2*pi/alat
    ispin: int = 1
    gamma_only: bool = False
    scalef: float = 1.0
    ngw: int = 0
    igwx: int = 0                    # número de ondas planas
    npol: int = 1
    nbnd: int = 0
    b: np.ndarray = None             # (3,3) vectores recíprocos, 2*pi/alat
    mill: np.ndarray = None          # (igwx, 3) índices de Miller
    coef: np.ndarray = None          # (nbnd, npol*igwx) complejos


class _Lector:
    """Lector de registros Fortran sin formato."""

    def __init__(self, ruta):
        self.f = open(ruta, "rb")

    def registro(self) -> bytes:
        cab = self.f.read(4)
        if len(cab) < 4:
            raise EOFError("fin de archivo inesperado")
        n = int(np.frombuffer(cab, dtype="<i4")[0])
        datos = self.f.read(n)
        self.f.read(4)               # marca de cierre
        return datos

    def cerrar(self):
        self.f.close()

    def __enter__(self):
        return self

    def __exit__(self, *a):
        self.cerrar()


def leer_wfc(ruta, bandas=None) -> WfcK:
    """Lee un wfc<N>.dat. `bandas` limita qué bandas se cargan en memoria.

    Cargar las 200 bandas de una supercelda grande puede ser un gigabyte;
    con `bandas` se leen solo las que interesan y las demás se saltan sin
    materializarlas.
    """
    ruta = Path(ruta)
    w = WfcK()
    with _Lector(ruta) as r:
        d = r.registro()
        w.ik = int(np.frombuffer(d[0:4], dtype="<i4")[0])
        w.xk = np.frombuffer(d[4:28], dtype="<f8").copy()
        w.ispin = int(np.frombuffer(d[28:32], dtype="<i4")[0])
        w.gamma_only = bool(np.frombuffer(d[32:36], dtype="<i4")[0])
        w.scalef = float(np.frombuffer(d[36:44], dtype="<f8")[0])

        d = r.registro()
        w.ngw, w.igwx, w.npol, w.nbnd = [
            int(x) for x in np.frombuffer(d[:16], dtype="<i4")]

        d = r.registro()
        w.b = np.frombuffer(d[:72], dtype="<f8").reshape(3, 3).copy()

        d = r.registro()
        w.mill = np.frombuffer(d[:12 * w.igwx],
                               dtype="<i4").reshape(w.igwx, 3).copy()

        pedidas = range(w.nbnd) if bandas is None else set(bandas)
        n = w.npol * w.igwx
        fuera = []
        for ib in range(w.nbnd):
            d = r.registro()
            if ib in pedidas:
                fuera.append(np.frombuffer(d[:16 * n],
                                           dtype="<c16").copy())
        w.coef = np.array(fuera) if fuera else None
    return w


_ESPINES = {"up": "up", "arriba": "up", "1": "up",
            "dw": "dw", "down": "dw", "abajo": "dw", "2": "dw"}


def _numero_k(f: Path) -> int:
    return int("".join(c for c in f.stem if c.isdigit()) or 0)


def es_lsda(save_dir) -> bool:
    """¿El cálculo guardó las funciones de onda por canal de espín?

    Con lsda pw.x no escribe wfc<N>.dat sino wfcup<N>.dat y wfcdw<N>.dat.
    """
    p = Path(save_dir)
    return any(p.glob("wfcup*.dat")) or any(p.glob("wfcdw*.dat"))


def buscar_wfc(save_dir, spin: str = "up") -> list:
    """Los wfc*.dat de UN canal de espín, ordenados por número de k.

    `spin` es "up" (por omisión) o "dw". En un cálculo con lsda pw.x escribe
    wfcup<N>.dat y wfcdw<N>.dat; el glob `wfc*.dat` a secas los cogía todos
    y, ordenados por número, los intercalaba (up1, dw1, up2, dw2, ...), así
    que la mitad de los puntos k salían del canal equivocado. Aquí nunca se
    mezclan: si hay archivos por canal se devuelve solo el pedido; si no los
    hay (cálculo sin polarización de espín) se devuelven los wfc<N>.dat y
    `spin` no cambia nada.
    """
    p = Path(save_dir)
    canal = _ESPINES.get(str(spin).lower())
    if canal is None:
        raise ValueError(f"canal de espín desconocido: {spin!r} "
                         f"(usa 'up' o 'dw').")
    if es_lsda(p):
        archivos = list(p.glob(f"wfc{canal}*.dat"))
    else:
        archivos = [f for f in p.glob("wfc*.dat")
                    if f.stem[3:].isdigit()]
    return sorted(archivos, key=_numero_k)
