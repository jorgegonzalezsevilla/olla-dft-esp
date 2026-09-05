# Olla-DFT — command-line toolkit for Quantum ESPRESSO
# Copyright (C) 2026 Jorge Enrique González Sevilla
# SPDX-License-Identifier: AGPL-3.0-or-later
# This program is free software: you can redistribute it and/or modify it
# under the terms of the GNU Affero General Public License as published by the Free
# Software Foundation, either version 3 of the License, or (at your option)
# any later version. See the LICENSE file for details.

"""Energías de formación, casco convexo y estabilidad de fases.

A partir de un conjunto de cálculos ya hechos —nada nuevo que correr— sale
la pregunta que de verdad importa de un material: ¿es estable, o se
descompone en otra cosa?

    E_f(por átomo) = [E(compuesto) - SUMA_i n_i * mu_i] / N

con mu_i la energía por átomo de la referencia elemental. El CASCO CONVEXO
inferior de E_f contra la composición marca las fases estables; la
distancia vertical de un punto al casco (E_hull) mide cuánta energía por
átomo gana descomponiéndose. E_hull = 0 significa estable; por encima de
unos 25 meV/átomo, la fase difícilmente se sintetiza como fase única.

REQUISITOS QUE NO SON NEGOCIABLES
---------------------------------
Todas las energías tienen que venir de cálculos con el MISMO funcional,
pseudos, cutoffs y tratamiento de ocupaciones. Restar energías de cálculos
distintos da un número sin sentido y sin ningún aviso, así que este módulo
usa la auditoría de `olla-dft audit` y SE NIEGA a construir el casco si el
conjunto no es homogéneo.

Y una limitación honesta: esto es energía a 0 K sin punto cero ni entropía.
Un compuesto puede salir por encima del casco y sintetizarse igual si lo
estabiliza la entropía a la temperatura de reacción.
"""

from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

from qekit.core import provenance
from qekit.core import style as qstyle
from qekit.core.errors import ErrorDeUso


@dataclass
class Fase:
    nombre: str = ""
    ruta: str = ""
    conteo: dict = field(default_factory=dict)   # símbolo -> nº de átomos
    energia: float = None       # eV, energía total de la celda
    natoms: int = 0
    x: np.ndarray = None        # fracciones molares, en el orden elementos
    e_form: float = None        # eV/átomo
    e_hull: float = None        # eV/átomo por encima del casco
    en_casco: bool = False


@dataclass
class HullResult:
    elementos: list = field(default_factory=list)
    fases: list = field(default_factory=list)
    referencias: dict = field(default_factory=dict)   # símbolo -> eV/átomo
    faltan_ref: list = field(default_factory=list)
    warnings: list = field(default_factory=list)


def _conteo(simbolos) -> dict:
    d = {}
    for s in simbolos:
        d[s] = d.get(s, 0) + 1
    return d


def from_runs(runs, elementos=None) -> HullResult:
    """Construye las fases a partir de una lista de RunInfo de `audit`."""
    from ase import Atoms

    res = HullResult()
    origenes = {getattr(x, "origen", "dft") for x in runs if x.ok}
    if len(origenes) > 1:
        res.warnings.append(
            "el conjunto mezcla energias de DFT y de potenciales "
            "aprendidos: " + ", ".join(sorted(origenes)) +
            ".\nSon superficies de energia distintas; un casco convexo "
            "construido con ambas\nno significa nada.")
        return res
    fases = []
    for x in runs:
        if not x.ok:
            continue
        r = x.result
        if (r.calculation or "").lower() in ("nscf", "bands"):
            continue
        if r.total_energy is None or r.converged is False:
            continue
        c = _conteo(r.symbols)
        nat = sum(c.values())
        if nat == 0:
            continue
        try:
            nombre = Atoms(symbols=r.symbols).get_chemical_formula()
        except Exception:                              # noqa: BLE001
            nombre = "".join(sorted(c))
        fases.append(Fase(nombre=nombre, ruta=x.path, conteo=c,
                          energia=r.total_energy, natoms=nat))

    if not fases:
        res.warnings.append("no hay cálculos con energía utilizable")
        return res

    els = elementos or sorted({s for f in fases for s in f.conteo})
    res.elementos = list(els)

    # referencias elementales: la energía por átomo MÁS BAJA de cada
    # elemento puro que haya en el conjunto
    for el in els:
        puros = [f for f in fases
                 if set(f.conteo) == {el}]
        if puros:
            res.referencias[el] = min(f.energia / f.natoms for f in puros)
        else:
            res.faltan_ref.append(el)

    if res.faltan_ref:
        res.warnings.append(
            "faltan las referencias elementales de: "
            + ", ".join(res.faltan_ref)
            + ".\nSin ellas no se puede definir una energía de formación: "
            "hay que calcular\ncada elemento puro en su fase estable, con "
            "los mismos parámetros.")
        res.fases = fases
        return res

    for f in fases:
        ref = sum(n * res.referencias[s] for s, n in f.conteo.items())
        f.e_form = (f.energia - ref) / f.natoms
        f.x = np.array([f.conteo.get(el, 0) / f.natoms for el in els])

    res.fases = fases
    _casco(res)
    return res


def from_table(filas, elementos=None) -> HullResult:
    """Casco a partir de una tabla explícita de (fórmula, energía total).

    `filas` es una lista de (nombre, {simbolo: n}, energia_eV). Sirve para
    mezclar energías que no salieron de Olla-DFT, y es lo que usan las
    pruebas: con datos sintéticos el casco tiene una respuesta exacta
    conocida y se puede verificar de verdad.
    """
    res = HullResult()
    fases = []
    for nombre, conteo, energia in filas:
        nat = sum(conteo.values())
        fases.append(Fase(nombre=nombre, conteo=dict(conteo),
                          energia=float(energia), natoms=nat))
    els = elementos or sorted({s for f in fases for s in f.conteo})
    res.elementos = list(els)
    for el in els:
        puros = [f for f in fases if set(f.conteo) == {el}]
        if puros:
            res.referencias[el] = min(f.energia / f.natoms for f in puros)
        else:
            res.faltan_ref.append(el)
    if res.faltan_ref:
        res.warnings.append("faltan referencias elementales: "
                            + ", ".join(res.faltan_ref))
        res.fases = fases
        return res
    for f in fases:
        ref = sum(n * res.referencias[s] for s, n in f.conteo.items())
        f.e_form = (f.energia - ref) / f.natoms
        f.x = np.array([f.conteo.get(el, 0) / f.natoms for el in els])
    res.fases = fases
    _casco(res)
    return res


def _casco(res: HullResult) -> None:
    """Casco convexo inferior de E_f contra composición."""
    els = res.elementos
    fases = res.fases
    if len(els) == 1:
        for f in fases:
            f.e_hull = f.e_form
            f.en_casco = abs(f.e_form) < 1e-9
        return

    # coordenadas independientes: las primeras n-1 fracciones
    pts = np.array([np.append(f.x[:-1], f.e_form) for f in fases])

    if len(els) == 2:
        # en binario basta recorrer la envolvente inferior ordenando por x
        orden = np.argsort(pts[:, 0])
        casco = []
        for i in orden:
            while len(casco) >= 2:
                (x1, y1), (x2, y2) = pts[casco[-2]][:2], pts[casco[-1]][:2]
                x3, y3 = pts[i][:2]
                # quitar el punto de en medio si no dobla hacia abajo
                if (x2 - x1) * (y3 - y1) - (y2 - y1) * (x3 - x1) <= 0:
                    casco.pop()
                else:
                    break
            casco.append(int(i))
        xs = pts[casco, 0]
        ys = pts[casco, 1]
        for f, p in zip(fases, pts):
            f.e_hull = float(p[1] - np.interp(p[0], xs, ys))
            f.en_casco = f.e_hull < 1e-9
        return

    # ternario o más: casco convexo n-dimensional, quedándose con las
    # facetas cuya normal apunta hacia abajo en el eje de energía
    try:
        from scipy.spatial import ConvexHull, Delaunay
    except ImportError:                                # pragma: no cover
        res.warnings.append("scipy no está disponible: no se puede calcular "
                            "el casco de un sistema de 3+ elementos")
        return
    try:
        ch = ConvexHull(pts)
    except Exception as exc:                           # noqa: BLE001
        res.warnings.append(f"no se pudo construir el casco: {exc}")
        return
    abajo = [s for s, eq in zip(ch.simplices, ch.equations) if eq[-2] < -1e-12]
    if not abajo:
        res.warnings.append("el casco no tiene facetas inferiores; "
                            "¿faltan fases?")
        return
    verts = sorted({int(i) for s in abajo for i in s})
    tri = Delaunay(pts[verts][:, :-1])
    for f, p in zip(fases, pts):
        s = tri.find_simplex(p[:-1])
        if s < 0:
            f.e_hull = None
            continue
        idx = [verts[j] for j in tri.simplices[s]]
        A = np.vstack([pts[idx][:, :-1].T, np.ones(len(idx))])
        try:
            w = np.linalg.solve(A, np.append(p[:-1], 1.0))
        except np.linalg.LinAlgError:
            f.e_hull = None
            continue
        f.e_hull = float(p[-1] - np.dot(w, pts[idx][:, -1]))
        f.en_casco = f.e_hull < 1e-9


def report(res: HullResult, umbral: float = 0.025) -> str:
    lines = ["--- Estabilidad de fases (casco convexo) ---",
             f"Elementos: {', '.join(res.elementos) or '?'}  |  "
             f"fases: {len(res.fases)}"]
    if res.referencias:
        lines.append("Referencias elementales (eV/átomo): " + "  ".join(
            f"{k} = {v:.4f}" for k, v in res.referencias.items()))
    if res.faltan_ref:
        lines += ["", "AVISO: " + res.warnings[0] if res.warnings else ""]
        return "\n".join(lines)

    lines += ["", f"{'fase':>14s} {'E_f (eV/át)':>13s} "
                  f"{'E_hull (eV/át)':>15s}  estado"]
    for f in sorted(res.fases, key=lambda z: (tuple(z.x), z.e_form or 0)):
        eh = f.e_hull
        if eh is None:
            estado, ehs = "fuera del dominio", "     n/d"
        elif f.en_casco:
            estado, ehs = "ESTABLE (en el casco)", f"{0.0:15.4f}"
        elif eh <= umbral:
            estado = f"metaestable (< {umbral*1000:.0f} meV/át)"
            ehs = f"{eh:15.4f}"
        else:
            estado, ehs = "inestable", f"{eh:15.4f}"
        lines.append(f"{f.nombre:>14s} {f.e_form:13.4f} {ehs}  {estado}")

    lines += ["",
              "E_hull es cuánta energía por átomo gana la fase "
              "descomponiéndose en las\nfases del casco. Cero = estable.",
              "",
              "Esto es energía a 0 K, sin punto cero ni entropía: una fase "
              "por encima del\ncasco puede sintetizarse igual si la "
              "estabiliza la entropía a la temperatura\nde reacción. Y todo "
              "depende de que las energías vengan de cálculos con los\n"
              "mismos parámetros — pásalas antes por 'olla-dft audit'."]
    for w in res.warnings:
        lines.append(f"\nAVISO: {w}")
    return "\n".join(lines)


def export(res: HullResult, outdir: str = ".") -> list:
    out = Path(outdir); out.mkdir(parents=True, exist_ok=True)
    f = out / "CASCO_CONVEXO.dat"
    lines = [provenance.header("casco convexo",
                               {"elementos": ",".join(res.elementos)},
                               titulo="Estabilidad de fases"),
             f"# {'fase':>14s} {'natoms':>7s} {'E_total(eV)':>14s} "
             f"{'E_f(eV/at)':>12s} {'E_hull(eV/at)':>14s}  ruta"]
    for x in res.fases:
        ef = x.e_form if x.e_form is not None else float("nan")
        eh = x.e_hull if x.e_hull is not None else float("nan")
        lines.append(f"{x.nombre:>16s} {x.natoms:7d} {x.energia:14.6f} "
                     f"{ef:12.5f} {eh:14.5f}  {x.ruta}")
    f.write_text("\n".join(lines) + "\n")
    return [str(f)]


def plot(res: HullResult, outfile: str = "casco", formats="pdf,png",
         theme: str = None, family: str = None, background: str = None,
         palette=None, usetex: bool = None, width="single",
         journal: str = "generic", aspect: float = 0.7,
         mono: bool = False, dpi: int = None) -> list:
    """Diagrama de casco convexo (solo binarios)."""
    if len(res.elementos) != 2:
        raise ErrorDeUso("la gráfica del casco solo está implementada para "
                         "sistemas binarios; para ternarios usa los datos "
                         "de CASCO_CONVEXO.dat")
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    st = qstyle.apply(theme, family=family, background=background,
                      palette=palette, usetex=usetex, mono=mono)
    fig = plt.figure(figsize=qstyle.figure_size(width, journal, aspect),
                     layout="constrained")
    ax = qstyle.finish_axes(fig.add_subplot(111))
    c = qstyle.palette(3, mono=mono)

    x = np.array([f.x[1] for f in res.fases])
    y = np.array([f.e_form for f in res.fases])
    est = np.array([bool(f.en_casco) for f in res.fases])

    orden = np.argsort(x[est])
    ax.plot(x[est][orden], y[est][orden], "-", color=c[0], lw=st["line"],
            zorder=2)
    ax.plot(x[~est], y[~est], "o", ms=4, color=c[1], mec=st["background"],
            mew=0.5, label="inestable", zorder=3)
    ax.plot(x[est], y[est], "o", ms=5, color=c[0], mec=st["background"],
            mew=0.5, label="estable", zorder=4)
    ax.axhline(0.0, color=qstyle.INK_FAINT, lw=st["axis_line"])
    for f in res.fases:
        if f.en_casco and 0 < f.x[1] < 1:
            ax.annotate(qstyle.tex_safe(f.nombre),
                        xy=(f.x[1], f.e_form), xytext=(0, -11),
                        textcoords="offset points", ha="center",
                        fontsize=st["legend"])
    # margen abajo para que las etiquetas de las fases del casco no
    # choquen contra el eje
    lo, hi = float(min(y.min(), 0.0)), float(max(y.max(), 0.0))
    rango = max(hi - lo, 1e-6)
    ax.set_ylim(lo - 0.18 * rango, hi + 0.06 * rango)
    ax.set_xlabel(f"fracción de {qstyle.tex_safe(res.elementos[1])}")
    ax.set_ylabel(r"$E_\mathrm{f}$ (eV/átomo)")
    ax.set_xlim(-0.02, 1.02)
    ax.legend()
    written = qstyle.save(fig, outfile, formats, dpi=dpi,
                          modulo="casco convexo")
    plt.close(fig)
    return written
