# Olla-DFT — command-line toolkit for Quantum ESPRESSO
# Copyright (C) 2026 Jorge Enrique González Sevilla
# SPDX-License-Identifier: AGPL-3.0-or-later
# This program is free software: you can redistribute it and/or modify it
# under the terms of the GNU Affero General Public License as published by the Free
# Software Foundation, either version 3 of the License, or (at your option)
# any later version. See the LICENSE file for details.

"""Cargas atómicas y diferencia de densidad.

Tres análisis que se hacen sobre la densidad ya calculada, sin volver a
correr pw.x:

**Cargas de Löwdin** (de projwfc.x). Proyección de los estados sobre los
orbitales atómicos, ortogonalizada. Es barata y estable, y da la carga por
átomo y por orbital — útil para ver estados de oxidación y transferencia
de carga.

**Cargas de Bader** (sobre la rejilla de densidad). Divide el espacio en
cuencas de atracción del gradiente de la densidad: cada punto se asigna al
máximo al que lleva la subida más empinada. No depende de una base de
orbitales, que es su ventaja sobre Löwdin y Mulliken.

**Diferencia de densidad de carga.** rho(AB) - rho(A) - rho(B) sobre la
misma rejilla: dónde se acumula y dónde se vacía carga al formarse el
enlace. Es la figura estándar de un estudio de adsorción.

Unidades: pp.x escribe la densidad (plot_num=0) en e/bohr³ también cuando
el .cube lleva la rejilla en bohr, y `fields.read_cube` pasa la rejilla a
Å. Aquí se integra por eso con el volumen del vóxel en bohr³; los
volúmenes de cuenca se reportan en Å³.

Aviso sobre Bader: esta implementación es "on-grid" (subida más empinada
entre vecinos de la rejilla). Es la variante simple del método de Henkelman
y hereda su sesgo: las cargas dependen algo de la densidad de la malla, y
el error típico es de unas centésimas de electrón. Para números finos hace
falta la variante near-grid del código `bader` original. El reporte lo dice
y avisa si la carga total no se conserva.
"""

import re
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

from qekit.core import style as qstyle
from qekit.modules import fields
from qekit.core.errors import ErrorDeUso, FaltanDatos

BOHR = fields.BOHR

#: Unidades en que viene la densidad del .cube. pp.x escribe e/bohr³
#: (unidades atómicas, el convenio del formato cube).
DENSITY_UNITS = ("e/bohr3", "e/A3")


def _voxel_volume(cube: fields.CubeData, density_units: str) -> tuple:
    """(dV para integrar la densidad, dV en Å³) de un vóxel del cube.

    Las aristas de `cube.axes` están en Å; si la densidad está en e/bohr³
    el volumen con el que se integra tiene que estar en bohr³ para que
    ρ·dV sea un número de electrones.
    """
    if density_units not in DENSITY_UNITS:
        raise ErrorDeUso(
            f"unidades de densidad '{density_units}' desconocidas; "
            f"opciones: {', '.join(DENSITY_UNITS)}")
    celda = cube.axes * np.array(cube.shape)[:, None]
    dv_A3 = abs(np.linalg.det(celda)) / float(np.prod(cube.shape))
    dv = dv_A3 / BOHR ** 3 if density_units == "e/bohr3" else dv_A3
    return dv, dv_A3


def valence_from_pseudos(symbols, pseudo_dir) -> np.ndarray:
    """Electrones de valencia por átomo, leídos de los UPF de `pseudo_dir`.

    Devuelve None (sin reventar) si la carpeta no existe o algún elemento
    no tiene UPF con `z_valence` legible: la columna 'neta' queda en n/d y
    la CLI avisa.
    """
    if not pseudo_dir:
        return None
    pdir = Path(pseudo_dir).expanduser()
    if not pdir.is_dir():
        return None
    from qekit.core import pseudo as ps
    try:
        info = ps.resolve(list(symbols), str(pdir))
    except Exception:                                    # noqa: BLE001
        return None
    zval = {s: d.get("z_valence") for s, d in info.items()}
    if any(zval.get(s) is None for s in symbols):
        return None
    return np.array([float(zval[s]) for s in symbols])


# ----------------------------------------------------------------------
# Cargas de Löwdin desde projwfc.x
# ----------------------------------------------------------------------
_RE_LOWDIN = re.compile(
    r"Atom\s*#\s*(\d+):\s*total charge\s*=\s*([\d.]+)", re.IGNORECASE)
_RE_SPILL = re.compile(r"Spilling Parameter:\s*([\d.]+)", re.IGNORECASE)


@dataclass
class LowdinResult:
    charges: np.ndarray = None       # carga total por átomo (electrones)
    symbols: list = field(default_factory=list)
    valence: np.ndarray = None       # electrones de valencia del pseudo
    spilling: float = None


def read_lowdin(projwfc_out, symbols=None, valence=None) -> LowdinResult:
    """Lee las cargas de Löwdin de la salida de projwfc.x."""
    txt = Path(projwfc_out).read_text(errors="ignore")
    pares = _RE_LOWDIN.findall(txt)
    if not pares:
        raise FaltanDatos(
            f"no hay cargas de Löwdin en '{projwfc_out}'. projwfc.x las "
            "imprime solo si corrió completo; revisa su salida.")
    orden = sorted((int(i), float(q)) for i, q in pares)
    res = LowdinResult(charges=np.array([q for _, q in orden]))
    m = _RE_SPILL.search(txt)
    if m:
        res.spilling = float(m.group(1))
    if symbols is not None:
        res.symbols = list(symbols)
    if valence is not None:
        res.valence = np.asarray(valence, dtype=float)
    return res


def report_lowdin(res: LowdinResult) -> str:
    lines = ["--- Cargas de Löwdin ---"]
    if res.spilling is not None:
        lines.append(f"Parámetro de derrame (spilling): {res.spilling:.4f}")
        if res.spilling > 0.05:
            lines.append(
                "  AVISO: por encima de ~0.05 la base atómica no describe "
                "bien los estados;\n  las cargas proyectadas pierden "
                "significado.")
    lines += ["", f"{'átomo':>8s} {'especie':>8s} {'carga(e)':>10s} "
                  f"{'neta':>8s}"]
    for i, q in enumerate(res.charges):
        sym = res.symbols[i] if i < len(res.symbols) else "?"
        if res.valence is not None and i < len(res.valence):
            neta = res.valence[i] - q
            lines.append(f"{i+1:8d} {sym:>8s} {q:10.4f} {neta:+8.3f}")
        else:
            lines.append(f"{i+1:8d} {sym:>8s} {q:10.4f} {'n/d':>8s}")
    lines += ["",
              "La 'neta' es Z_valencia − carga proyectada: positiva = el "
              "átomo cedió carga.",
              "Löwdin depende de la base de orbitales del pseudo; es útil "
              "para COMPARAR\nátomos entre sí, no como carga absoluta."]
    return "\n".join(lines)


# ----------------------------------------------------------------------
# Bader on-grid
# ----------------------------------------------------------------------
@dataclass
class BaderResult:
    charges: np.ndarray = None       # electrones en la cuenca de cada átomo
    volumes: np.ndarray = None       # A^3
    symbols: list = field(default_factory=list)
    valence: np.ndarray = None
    total: float = None
    total_grid: float = None


def bader(cube: fields.CubeData, positions: np.ndarray, symbols=None,
          valence=None, density_units: str = "e/bohr3") -> BaderResult:
    """Cargas de Bader por subida más empinada sobre la rejilla.

    Cada punto de la rejilla se mueve al vecino (de los 26) con mayor
    densidad, hasta llegar a un máximo local; los puntos que terminan en el
    mismo máximo forman una cuenca. Cada cuenca se asigna al átomo más
    cercano a su máximo.

    `density_units` dice en qué unidades viene `cube.data`: "e/bohr3" (lo
    que escribe pp.x) o "e/A3". La suma de las cargas de cuenca debe dar el
    número de electrones de valencia de la celda.
    """
    rho = np.asarray(cube.data, dtype=float)
    n1, n2, n3 = rho.shape
    npts = n1 * n2 * n3

    # desplazamientos a los 26 vecinos y su longitud real (celda no cúbica)
    offs, longs = [], []
    for a in (-1, 0, 1):
        for b in (-1, 0, 1):
            for c in (-1, 0, 1):
                if a == b == c == 0:
                    continue
                d = a * cube.axes[0] + b * cube.axes[1] + c * cube.axes[2]
                offs.append((a, b, c))
                longs.append(np.linalg.norm(d))
    longs = np.array(longs)

    # pendiente hacia cada vecino: (rho_vecino - rho) / distancia
    mejor = np.zeros(rho.shape, dtype=np.int8)
    mejor_pend = np.zeros(rho.shape)
    for k, ((a, b, c), L) in enumerate(zip(offs, longs)):
        vec = np.roll(np.roll(np.roll(rho, -a, 0), -b, 1), -c, 2)
        pend = (vec - rho) / L
        mask = pend > mejor_pend
        mejor_pend[mask] = pend[mask]
        mejor[mask] = k + 1          # 0 = ya es máximo

    # seguir el camino hasta el máximo, con compresión de caminos
    idx = np.arange(npts).reshape(rho.shape)
    destino = np.full(npts, -1, dtype=np.int64)
    salto = np.empty(npts, dtype=np.int64)
    plano = mejor.ravel()
    for k, (a, b, c) in enumerate(offs, start=1):
        m = plano == k
        if not m.any():
            continue
        vecino_idx = np.roll(np.roll(np.roll(idx, -a, 0), -b, 1), -c, 2)
        salto[m] = vecino_idx.ravel()[m]
    salto[plano == 0] = np.flatnonzero(plano == 0)   # un máximo se apunta a sí

    # iterar salto[salto[...]] hasta punto fijo (doblado de caminos)
    destino = salto.copy()
    for _ in range(64):
        nuevo = destino[destino]
        if np.array_equal(nuevo, destino):
            break
        destino = nuevo

    maximos = np.unique(destino)
    # posición cartesiana de cada máximo
    mi = np.array(np.unravel_index(maximos, rho.shape)).T
    pos_max = (mi[:, 0:1] * cube.axes[0] + mi[:, 1:2] * cube.axes[1]
               + mi[:, 2:3] * cube.axes[2]) + cube.origin

    pos = np.asarray(positions, dtype=float)
    # asignar cada máximo al átomo más cercano (con imágenes periódicas)
    celda = np.array([cube.axes[0] * n1, cube.axes[1] * n2,
                      cube.axes[2] * n3])
    inv = np.linalg.inv(celda)
    d = pos_max[:, None, :] - pos[None, :, :]
    f = d @ inv
    f -= np.round(f)
    d = f @ celda
    atomo_de_max = np.argmin(np.linalg.norm(d, axis=2), axis=1)

    mapa = np.zeros(destino.max() + 1, dtype=np.int64)
    mapa[maximos] = atomo_de_max
    atomo_de_punto = mapa[destino]

    # dv en las unidades de la densidad (bohr³ para pp.x) para que ρ·dV
    # sean electrones; dv_A3 solo para reportar volúmenes en Å³
    dv, dv_A3 = _voxel_volume(cube, density_units)
    rho_flat = rho.ravel()
    nat = len(pos)
    cargas = np.zeros(nat)
    vols = np.zeros(nat)
    for i in range(nat):
        m = atomo_de_punto == i
        cargas[i] = float(rho_flat[m].sum()) * dv
        vols[i] = float(m.sum()) * dv_A3

    res = BaderResult(charges=cargas, volumes=vols, total=float(cargas.sum()),
                      total_grid=float(rho_flat.sum()) * dv)
    if symbols is not None:
        res.symbols = list(symbols)
    if valence is not None:
        res.valence = np.asarray(valence, dtype=float)
    return res


def report_bader(res: BaderResult) -> str:
    lines = ["--- Cargas de Bader (on-grid) ---",
             f"Carga total en las cuencas: {res.total:.4f} e   "
             f"(integral de la rejilla: {res.total_grid:.4f} e)"]
    err = abs(res.total - res.total_grid)
    if err > 1e-3:
        lines.append(f"  AVISO: se perdieron {err:.4f} e al repartir en "
                     "cuencas; la malla del cube es demasiado gruesa.")
    if res.valence is not None and len(res.valence) == len(res.charges):
        z_tot = float(np.sum(res.valence))
        lines.append(f"Electrones de valencia según los UPF: {z_tot:.4f} e")
        if abs(res.total_grid - z_tot) > 0.05 * max(z_tot, 1.0):
            lines.append(
                f"  AVISO: la integral de la rejilla ({res.total_grid:.3f} e) "
                f"no coincide con Z_valencia total ({z_tot:.3f} e).\n"
                "  Revisa que el cube sea la densidad de valencia completa "
                "(plot_num=0) y que\n  los UPF de --pseudo-dir sean los del "
                "cálculo.")
    lines += ["", f"{'átomo':>8s} {'especie':>8s} {'carga(e)':>10s} "
                  f"{'neta':>8s} {'volumen(Å³)':>12s}"]
    for i, (q, v) in enumerate(zip(res.charges, res.volumes)):
        sym = res.symbols[i] if i < len(res.symbols) else "?"
        if res.valence is not None and i < len(res.valence):
            lines.append(f"{i+1:8d} {sym:>8s} {q:10.4f} "
                         f"{res.valence[i]-q:+8.3f} {v:12.3f}")
        else:
            lines.append(f"{i+1:8d} {sym:>8s} {q:10.4f} {'n/d':>8s} "
                         f"{v:12.3f}")
    lines += ["",
              "Método on-grid: la asignación sigue la subida más empinada "
              "entre vecinos de\nla rejilla. Hereda el sesgo de malla del "
              "método (centésimas de electrón);\npara números finos usa la "
              "variante near-grid del código `bader` de Henkelman.",
              "OJO: la densidad debe incluir la carga de valencia completa "
              "(plot_num=0).\nCon pseudopotenciales, la carga de Bader se "
              "compara contra Z_valencia, no Z."]
    return "\n".join(lines)


# ----------------------------------------------------------------------
# Diferencia de densidad de carga
# ----------------------------------------------------------------------
def difference(total: fields.CubeData, partes: list) -> fields.CubeData:
    """rho(total) - suma de rho(partes), sobre la misma rejilla."""
    d = np.array(total.data, dtype=float)
    for p in partes:
        if tuple(p.shape) != tuple(total.shape):
            raise ErrorDeUso(
                f"las rejillas no coinciden: {total.shape} contra {p.shape}. "
                "Las tres densidades tienen que salir de cálculos con la "
                "MISMA celda, la misma malla FFT y los mismos cutoffs — si "
                "no, la resta no significa nada.")
        d -= np.asarray(p.data, dtype=float)
    return fields.CubeData(origin=total.origin, axes=total.axes,
                           shape=total.shape, data=d, natoms=total.natoms)


def report_difference(cube: fields.CubeData, axis: int = 2,
                      density_units: str = "e/bohr3") -> str:
    z, prof = fields.planar_average(cube, axis=axis)
    # mismo criterio que en bader(): la densidad de pp.x está en e/bohr³
    dv, _dv_A3 = _voxel_volume(cube, density_units)
    d = np.asarray(cube.data, dtype=float)
    neto = float(d.sum()) * dv
    acum = float(d[d > 0].sum()) * dv
    return "\n".join([
        "--- Diferencia de densidad de carga ---",
        f"Carga neta transferida: {neto:+.4f} e   "
        f"(debería ser ~0 si las partes suman el total)",
        f"Carga acumulada (regiones positivas): {acum:.4f} e",
        f"Máximo del perfil planar: {prof.max():+.5f}  en "
        f"z = {z[int(np.argmax(prof))]:.2f} Å",
        f"Mínimo del perfil planar: {prof.min():+.5f}  en "
        f"z = {z[int(np.argmin(prof))]:.2f} Å",
        "",
        "Positivo = se acumula carga al formarse el sistema; negativo = se "
        "vacía.",
    ])


def plot_difference(cube: fields.CubeData, outfile: str = "diferencia_carga",
                    axis: int = 2, formats="pdf,png", theme: str = None,
                    family: str = None, background: str = None, palette=None,
                    usetex: bool = None, width="single",
                    journal: str = "generic", aspect: float = 0.62,
                    mono: bool = False, dpi: int = None) -> list:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    st = qstyle.apply(theme, family=family, background=background,
                      palette=palette, usetex=usetex, mono=mono)
    fig = plt.figure(figsize=qstyle.figure_size(width, journal, aspect),
                     layout="constrained")
    ax = qstyle.finish_axes(fig.add_subplot(111))
    z, prof = fields.planar_average(cube, axis=axis)
    c = qstyle.palette(2, mono=mono)
    ax.plot(z, prof, color=c[0], lw=st["line"])
    ax.fill_between(z, 0, prof, where=prof > 0, color=c[0], alpha=0.18, lw=0)
    ax.fill_between(z, 0, prof, where=prof < 0, color=c[1], alpha=0.18, lw=0)
    ax.axhline(0.0, color=qstyle.INK_FAINT, lw=st["axis_line"])
    ax.set_xlabel(f"{'xyz'[axis]} (Å)")
    ax.set_ylabel(r"$\Delta\rho$ promediado en el plano (e/bohr$^3$)")
    ax.set_xlim(z.min(), z.max())
    written = qstyle.save(fig, outfile, formats, dpi=dpi,
                          modulo="diferencia de densidad")
    plt.close(fig)
    return written
