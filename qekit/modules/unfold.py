# Olla-DFT — command-line toolkit for Quantum ESPRESSO
# Copyright (C) 2026 Jorge Enrique González Sevilla
# SPDX-License-Identifier: AGPL-3.0-or-later
# This program is free software: you can redistribute it and/or modify it
# under the terms of the GNU Affero General Public License as published by the Free
# Software Foundation, either version 3 of the License, or (at your option)
# any later version. See the LICENSE file for details.

"""Desdoblamiento de bandas: recuperar la dispersión de una supercelda.

EL PROBLEMA
-----------
Una supercelda de N celdas primitivas tiene una zona de Brillouin N veces
más pequeña, así que su estructura de bandas sale PLEGADA: donde la celda
primitiva tenía una banda, la supercelda tiene N ramas amontonadas. Mirar
esa maraña no dice nada.

Eso importa justo cuando la supercelda hace falta: una aleación, un
dopante, un defecto, una capa sobre otra. En esos casos uno quiere saber
qué le pasó a la banda del material original — si el dopante abrió un
gap, si el defecto metió un estado, cuánto se difuminó la banda.

CÓMO SE DESHACE
---------------
Cada banda de la supercelda se escribe como una suma de ondas planas.
Cada onda plana tiene un vector de onda bien definido q = K + G. La
pregunta "¿cuánto de esta banda de supercelda pertenece al punto k de la
celda primitiva?" tiene una respuesta exacta:

    P(k) = suma de |C(G)|^2 sobre las G tales que K + G ≡ k
           (módulo la red recíproca PRIMITIVA)

Ese peso espectral es lo que se dibuja. Si la supercelda es perfecta
—sin defecto ni desorden—, cada banda desdoblada tiene peso 1 en su k y
0 en los demás, y se recupera exactamente la banda primitiva. Cuanto más
rompe la periodicidad el defecto, más se reparte el peso: la banda se
difumina. Esa difuminación ES el resultado físico.

LO QUE HACE FALTA
-----------------
Las funciones de onda. Un cálculo con `disk_io='nowf'` o `'low'` no las
guarda, y entonces esto no se puede hacer: hay que repetirlo. Olla-DFT lo
comprueba antes de nada y lo dice.

También hace falta la MATRIZ de la supercelda: los enteros M tales que
A = M·a. Olla-DFT la deduce de las dos celdas y verifica que salgan enteros;
si no salen, es que las dos estructuras no están relacionadas por una
supercelda y el desdoblamiento no tiene sentido.
"""

from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

from qekit.core import provenance, wfc
from qekit.core import style as qstyle
from qekit.core.errors import ErrorDeUso

TOL_ENTERO = 1e-4


@dataclass
class Desdoblado:
    kpath: np.ndarray = None        # (nk, 3) fraccionarios de la primitiva
    distancias: np.ndarray = None   # (nk,) para el eje x
    energias: np.ndarray = None     # (nk, nbnd) eV
    pesos: np.ndarray = None        # (nk, nbnd) en [0, 1]
    etiquetas: list = field(default_factory=list)
    e_fermi: float = None
    M: np.ndarray = None            # matriz entera de la supercelda
    ncel: int = 0
    avisos: list = field(default_factory=list)
    spin: str = None                # canal desdoblado si el cálculo es lsda


def aviso_lsda(spin: str, otro: str) -> str:
    """El texto que va al reporte cuando el cálculo es de espín polarizado."""
    return (f"AVISO: el cálculo es de espín polarizado (lsda) y aquí solo se "
            f"ha desdoblado el\ncanal '{spin}' (wfc{spin}<N>.dat y sus "
            f"energías). El otro canal no se mezcla ni se\nsuma: para verlo "
            f"repite el desdoblamiento con --spin {otro}.")


# ----------------------------------------------------------------------
# Relación entre las dos celdas
# ----------------------------------------------------------------------
def matriz_supercelda(celda_sc, celda_prim) -> np.ndarray:
    """M entera tal que A_sc = M · a_prim.

    Se comprueba que salga entera de verdad. Si no, las dos estructuras
    no están relacionadas por una supercelda —por ejemplo porque una está
    relajada y la otra no— y desdoblar sobre ellas da basura.
    """
    A = np.asarray(celda_sc, dtype=float)
    a = np.asarray(celda_prim, dtype=float)
    M = A @ np.linalg.inv(a)
    Mi = np.round(M)
    err = float(np.abs(M - Mi).max())
    if err > 1e-3:
        # Segunda oportunidad: las dos celdas pueden describir la misma
        # relacion pero estar ORIENTADAS distinto — pasa siempre que una
        # de las dos viene de un CIF, porque el formato reorienta la celda
        # a su forma canonica. La relacion A = M*a es invariante bajo una
        # rotacion comun, asi que se busca M sobre las METRICAS, que no
        # dependen de la orientacion.
        Mi = _m_por_metricas(A, a)
        if Mi is not None:
            return Mi
        raise ErrorDeUso(
            "la celda de la supercelda no es un múltiplo entero de la "
            f"primitiva (error {err:.4f}).\n"
            "Puede que una esté relajada y la otra no, o que la primitiva no "
            "sea la que\ncorresponde. El desdoblamiento necesita la relación "
            "exacta A = M·a.\n"
            f"M calculada:\n{np.array2string(M, precision=4)}")
    return Mi.astype(int)



def _m_por_metricas(A, a):
    """Busca M entera con G_sc = M G_p M^T, sin depender de la orientacion.

    El tensor metrico G = X X^T no cambia al rotar la celda, asi que
    comparar metricas encuentra la relacion aunque las dos estructuras
    esten escritas en ejes distintos (un CIF siempre reorienta la celda).

    Se busca FILA POR FILA: primero los vectores enteros con la longitud
    correcta, y solo entre esos se prueban las combinaciones. Barrer las
    nueve componentes a la vez son decenas de millones de casos y no
    termina nunca.
    """
    from itertools import product

    G_sc = A @ A.T
    G_p = a @ a.T
    det_sc = float(np.linalg.det(G_sc))
    det_p = float(np.linalg.det(G_p))
    if det_p <= 0 or det_sc <= 0:
        return None
    ncel = np.sqrt(det_sc / det_p)
    if abs(ncel - round(ncel)) > 1e-3:
        return None            # ni siquiera es un multiplo entero del volumen
    ncel = int(round(ncel))
    lim = int(np.ceil(ncel ** (1.0 / 3.0))) + 2
    rango = range(-lim, lim + 1)

    tol = 1e-3 * max(1.0, abs(G_sc).max())
    candidatos = [[], [], []]
    for v in product(rango, repeat=3):
        w = np.array(v, dtype=int)
        if not w.any():
            continue
        norma = float(w @ G_p @ w)
        for i in range(3):
            if abs(norma - G_sc[i, i]) < tol:
                candidatos[i].append(w)

    for v0 in candidatos[0]:
        for v1 in candidatos[1]:
            if abs(v0 @ G_p @ v1 - G_sc[0, 1]) > tol:
                continue
            for v2 in candidatos[2]:
                if abs(v0 @ G_p @ v2 - G_sc[0, 2]) > tol:
                    continue
                if abs(v1 @ G_p @ v2 - G_sc[1, 2]) > tol:
                    continue
                M = np.array([v0, v1, v2], dtype=int)
                if int(round(np.linalg.det(M))) == ncel:
                    return M
    return None


def pesos_de_k(w, M, k_sc_frac, k_prim_frac) -> np.ndarray:
    """Peso espectral de cada banda de la supercelda en un k primitivo.

    Se hace TODO en coordenadas fraccionarias, sin tocar alat ni los
    vectores reciprocos cartesianos. Los indices de Miller ya son las
    coordenadas de G en la base reciproca de la supercelda, y como
    b_sc = M^-T b_prim, las coordenadas de G en la base primitiva son
    simplemente

        g = mill @ M^-T

    Una onda plana pertenece al k primitivo pedido si g, corrida por el
    desplazamiento entero que separa ese k del propio k de la supercelda,
    tiene las tres componentes ENTERAS. Es exacto y no depende de en qué
    unidades escribiera QE nada.
    """
    if w.coef is None:
        raise ErrorDeUso("el archivo de funciones de onda no trae "
                         "coeficientes.")
    Minv_T = np.linalg.inv(M).T
    # desplazamiento entero entre el k pedido y el de la supercelda
    m0 = np.asarray(k_prim_frac, dtype=float) @ np.asarray(M, dtype=float).T \
        - np.asarray(k_sc_frac, dtype=float)
    if np.abs(m0 - np.round(m0)).max() > TOL_ENTERO:
        # ese k primitivo no se pliega sobre este k de supercelda
        return np.zeros(w.coef.shape[0])
    g = (w.mill - np.round(m0)[None, :]) @ Minv_T
    pertenece = np.all(np.abs(g - np.round(g)) < TOL_ENTERO, axis=1)

    c = w.coef
    if w.npol > 1:
        c = c.reshape(c.shape[0], w.npol, w.igwx)
        p2 = np.sum(np.abs(c) ** 2, axis=1)
    else:
        p2 = np.abs(c) ** 2
    total = p2.sum(axis=1)
    peso = p2[:, pertenece].sum(axis=1)
    return np.divide(peso, total, out=np.zeros_like(peso),
                     where=total > 0)


# ----------------------------------------------------------------------
# Todo el camino
# ----------------------------------------------------------------------
def desdoblar(calc_dir, celda_primitiva, kpath_frac=None,
              prefix: str = None, bandas=None,
              etiquetas=None, spin: str = "up") -> Desdoblado:
    """Desdobla las bandas de una supercelda sobre la ZB primitiva.

    `kpath_frac` son los puntos k EN COORDENADAS DE LA PRIMITIVA. Si no
    se dan, se usan los propios k del cálculo interpretados en la
    primitiva, que es lo que corresponde cuando el camino se generó ya
    plegado.

    `spin` elige el canal ("up" o "dw") en un cálculo con lsda. Se desdobla
    UN solo canal por llamada: las funciones de onda y las energías salen
    del mismo, y el reporte lo avisa.
    """
    from qekit.core import qeout

    p = Path(calc_dir)
    res = qeout.read_xml(str(p))
    save = _carpeta_save(p, prefix)
    archivos = wfc.buscar_wfc(save, spin=spin)
    lsda = wfc.es_lsda(save) or res.nspin == 2
    # índice del canal en res.eigenvalues (nspin, nk, nbnd)
    ispin = 1 if str(spin).lower() in ("dw", "down", "abajo", "2") else 0
    spin = "dw" if ispin else "up"
    if not archivos:
        if lsda:
            raise ErrorDeUso(
                f"no hay ningún wfc{spin}*.dat en {save}.\n"
                "El cálculo es de espín polarizado (lsda) pero falta el "
                f"canal '{spin}': revisa\nque pw.x terminó y guardó las "
                "funciones de onda de los dos canales.")
        raise ErrorDeUso(
            f"no hay ningún wfc*.dat en {save}.\n"
            "El cálculo no guardó las funciones de onda: eso pasa con "
            "disk_io='nowf' o\n'low'. Hay que repetirlo con disk_io='medium' "
            "o 'high' — el desdoblamiento\nnecesita los coeficientes de las "
            "ondas planas, no solo las energías.")

    celda_sc = np.array(res.cell, dtype=float)
    a_prim = np.array(celda_primitiva, dtype=float)
    M = matriz_supercelda(celda_sc, a_prim)
    ncel = int(round(abs(np.linalg.det(M))))
    # La primitiva se REDERIVA de la supercelda: a = M^-1 A. Asi las dos
    # quedan en el MISMO sistema de ejes aunque la que dio el usuario
    # viniera rotada (un CIF siempre reorienta la celda), y las
    # coordenadas fraccionarias del espacio reciproco salen bien.
    a_prim = np.linalg.inv(M) @ celda_sc

    # QE escribe xk y b en unidades de 2*pi/alat, con alat = |a1| de la
    # SUPERCELDA. La primitiva hay que expresarla en esas mismas unidades.
    k_sc = np.array(res.kpoints_frac, dtype=float)
    kfrac = k_sc.copy() if kpath_frac is None \
        else np.array(kpath_frac, dtype=float)
    if kpath_frac is None:
        # Los k del calculo estan en coordenadas fraccionarias de la
        # SUPERCELDA. Como B_sc = M^-T b_prim, pasar a coordenadas de la
        # primitiva es k_prim = k_sc M^-T, NO k_sc M: con la version
        # equivocada solo gamma sale bien, porque cero por cualquier cosa
        # sigue siendo cero.
        kfrac = kfrac @ np.linalg.inv(M).T

    energias, pesos = [], []
    n = min(len(archivos), len(kfrac))
    for i in range(n):
        w = wfc.leer_wfc(archivos[i], bandas=bandas)
        pesos.append(pesos_de_k(w, M, k_sc[i], kfrac[i]))
        if res.eigenvalues.ndim == 3:
            eig = res.eigenvalues[min(ispin, res.eigenvalues.shape[0] - 1)]
        else:
            eig = res.eigenvalues
        e = np.array(eig[i], dtype=float)
        energias.append(e if bandas is None else e[list(bandas)])

    d = Desdoblado(kpath=kfrac[:n], energias=np.array(energias),
                   pesos=np.array(pesos), M=M, ncel=ncel,
                   e_fermi=res.fermi,
                   etiquetas=list(etiquetas or []))
    d.distancias = _distancias(kfrac[:n], a_prim)
    d.spin = spin if lsda else None
    if lsda:
        otro = "dw" if ispin == 0 else "up"
        d.avisos.append(aviso_lsda(spin, otro))

    suma = d.pesos.sum(axis=1)
    if np.any(suma > len(d.pesos[0]) * 1.001):
        d.avisos.append("Algún peso salió mayor que 1: revisa la tolerancia.")
    frac_alta = float(np.mean(d.pesos > 0.9))
    if frac_alta > 0.95:
        d.avisos.append(
            "Casi todos los pesos valen 1: la supercelda parece PERFECTA "
            "(sin defecto\nni desorden). En ese caso el desdoblamiento "
            "reproduce exactamente las bandas\nprimitivas — que es la "
            "comprobación de que funciona, pero no un resultado nuevo.")
    return d


def _carpeta_save(p: Path, prefix: str = None) -> Path:
    for base in (p, p / "out"):
        if not base.is_dir():
            continue
        saves = sorted(base.glob("*.save"))
        if prefix:
            saves = [s for s in saves if s.stem == prefix] or saves
        if saves:
            return saves[0]
    raise ErrorDeUso(f"no se encontró ninguna carpeta .save bajo {p}.")


def _distancias(kfrac, a_prim) -> np.ndarray:
    """Distancia acumulada a lo largo del camino, en el espacio recíproco."""
    b_prim = 2.0 * np.pi * np.linalg.inv(np.asarray(a_prim, dtype=float).T)
    kcart = np.asarray(kfrac, dtype=float) @ b_prim
    d = np.zeros(len(kcart))
    if len(kcart) > 1:
        paso = np.linalg.norm(np.diff(kcart, axis=0), axis=1)
        # un salto grande es un cambio de rama del camino, no distancia
        paso[paso > 5 * np.median(paso[paso > 0]) if np.any(paso > 0)
             else False] = 0.0
        d[1:] = np.cumsum(paso)
    return d


# ----------------------------------------------------------------------
# Reporte y figura
# ----------------------------------------------------------------------
def report(d: Desdoblado) -> str:
    lines = ["--- Desdoblamiento de bandas ---",
             f"Supercelda: {d.ncel} celdas primitivas",
             f"M =\n{np.array2string(d.M, prefix='    ')}",
             f"Puntos k: {len(d.kpath)}   bandas: {d.energias.shape[1]}"]
    if d.e_fermi is not None:
        lines.append(f"Nivel de Fermi: {d.e_fermi:.4f} eV")

    p = d.pesos
    lines += ["", "Distribución del peso espectral:",
              f"  peso medio            {p.mean():.4f}",
              f"  fracción con peso > 0.9  {np.mean(p > 0.9) * 100:5.1f} %",
              f"  fracción con peso < 0.1  {np.mean(p < 0.1) * 100:5.1f} %"]
    lines += ["",
              "Un peso de 1 quiere decir que ese estado de la supercelda ES "
              "un estado de la\ncelda primitiva en ese k. Un peso repartido "
              "quiere decir que la periodicidad\nprimitiva está rota — y eso, "
              "no la posición de la banda, es lo que informa\nsobre el "
              "defecto o el desorden."]
    for a in d.avisos:
        lines += ["", a]
    return "\n".join(lines)


def export(d: Desdoblado, outdir: str = ".") -> list:
    out = Path(outdir); out.mkdir(parents=True, exist_ok=True)
    f = out / "UNFOLD.dat"
    filas = []
    for i in range(len(d.distancias)):
        for j in range(d.energias.shape[1]):
            e = d.energias[i, j] - (d.e_fermi or 0.0)
            filas.append([d.distancias[i], e, d.pesos[i, j]])
    np.savetxt(f, np.array(filas), fmt="%14.6f",
               header=provenance.header_plain(
                   "desdoblamiento de bandas",
                   {"celdas": d.ncel, "E_fermi_eV": d.e_fermi},
                   titulo="Peso espectral desdoblado") +
               "\n     distancia          E-EF(eV)          peso",
               comments="# ")
    txt = out / "UNFOLD.txt"
    txt.write_text(report(d) + "\n")
    return [str(f), str(txt)]


def plot(d: Desdoblado, outfile: str = "unfold", formats="pdf,png",
         emin: float = -6.0, emax: float = 6.0, theme: str = None,
         family: str = None, background: str = None, palette=None,
         usetex: bool = None, width="single", journal: str = "generic",
         mono: bool = False, dpi: int = None, escala: float = 60.0) -> list:
    try:
        import matplotlib
        matplotlib.use("Agg")
    except ImportError as exc:                          # pragma: no cover
        raise RuntimeError("matplotlib no está instalado.") from exc

    st = qstyle.apply(theme, family=family, background=background,
                      palette=palette, usetex=usetex, mono=mono)
    fig, ax = qstyle.new_figure(width, journal, 0.85)
    colores = qstyle.palette(4, mono=mono)

    ef = d.e_fermi or 0.0
    X, Y, S = [], [], []
    for i in range(len(d.distancias)):
        for j in range(d.energias.shape[1]):
            e = d.energias[i, j] - ef
            if emin <= e <= emax and d.pesos[i, j] > 0.005:
                X.append(d.distancias[i]); Y.append(e)
                S.append(d.pesos[i, j])
    S = np.array(S)
    # el tamaño del punto ES el peso: es la forma honesta de dibujarlo,
    # porque una banda difuminada tiene que VERSE difuminada
    ax.scatter(X, Y, s=escala * S, c=[colores[0]], alpha=0.55,
               linewidths=0, zorder=2)
    ax.axhline(0.0, color=qstyle.INK_FAINT, lw=st["axis_line"],
               dashes=[3.5, 2.0])
    ax.set_ylabel(r"$E - E_\mathrm{F}$ (eV)")
    ax.set_ylim(emin, emax)
    ax.set_xlim(d.distancias[0], d.distancias[-1])
    if d.etiquetas:
        pos = [t[0] for t in d.etiquetas]
        nom = [t[1] for t in d.etiquetas]
        ax.set_xticks(pos); ax.set_xticklabels(nom)
        for t in pos:
            ax.axvline(t, color=qstyle.GRID, lw=st["axis_line"], zorder=0)
    else:
        ax.set_xlabel("camino en el espacio recíproco")
    return qstyle.save(fig, outfile, formats, dpi=dpi)
