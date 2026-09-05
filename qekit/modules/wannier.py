# Olla-DFT — command-line toolkit for Quantum ESPRESSO
# Copyright (C) 2026 Jorge Enrique González Sevilla
# SPDX-License-Identifier: AGPL-3.0-or-later
# This program is free software: you can redistribute it and/or modify it
# under the terms of the GNU Affero General Public License as published by the Free
# Software Foundation, either version 3 of the License, or (at your option)
# any later version. See the LICENSE file for details.

"""Funciones de Wannier: bajar la estructura de bandas a un modelo pequeño.

Un cálculo de DFT te da ε_n(k) en los puntos k que calculaste y en ningún
otro. Si quieres la banda en un punto intermedio hay que volver a correr
pw.x. Eso es carísimo para las cosas que de verdad necesitan una malla fina:
transporte, superficies de Fermi, densidades de estados suaves, masas
efectivas por diferencias finitas.

Las funciones de Wannier resuelven eso. Son la transformada de Fourier de
los estados de Bloch a una base localizada en el espacio real:

    |R n⟩ = (V/(2π)³) ∫ dk e^{-ik·R} Σ_m U_mn(k) |ψ_mk⟩

Con esa base el hamiltoniano se vuelve una matriz pequeña H_mn(R) que decae
con |R|. Y una vez que la tienes, la banda en CUALQUIER k sale de
diagonalizar una matriz de num_wann×num_wann:

    H(k) = Σ_R e^{ik·R} H(R) / deg(R)

que en un portátil son microsegundos por punto k. Ahí está todo el negocio:
se paga una vez la malla gruesa de DFT y luego se interpola gratis.

Lo delicado es U(k). Si eliges la matriz identidad —o sea, si transformas
directamente las energías propias— la interpolación sale MAL, porque la fase
de cada |ψ_nk⟩ que devuelve un diagonalizador es arbitraria y cambia de un
punto k al siguiente: la función que transformas no es suave, y su
transformada de Fourier no decae. Ese es el error que hace que "interpolar
bandas" parezca fácil y dé basura. Aquí se compara explícitamente con y sin
gauge para que se vea.

La U buena se obtiene proyectando sobre orbitales localizados de prueba
(orbitales atómicos, o híbridos sp3 en un enlace) y ortonormalizando:

    A_mn(k) = ⟨ψ_mk|g_n⟩          U(k) = A (A†A)^{-1/2}

Esto es la "gauge de proyección" de Marzari-Vanderbilt, que es también el
punto de partida de wannier90 antes de minimizar la dispersión. No es
maximalmente localizada, pero para interpolar bandas suele bastar, y tiene
una ventaja grande para quien no tiene un clúster: **no hace falta instalar
wannier90**. Los A_mn y los solapes M los calcula pw2wannier90.x, que viene
con Quantum ESPRESSO, y el resto se hace aquí.

Si además tienes wannier90 instalado, `leer_hr()` lee su seedname_hr.dat y
el resto del módulo funciona igual, con funciones ya maximalmente
localizadas. Los dos caminos acaban en el mismo H(R).

Qué se comprueba antes de creerse nada:

  - La malla de puntos k tiene que ser COMPLETA (nosym, noinv). Con la malla
    reducida por simetría los solapes M^{k,b} no existen y pw2wannier90 se
    queja o, peor, da números sin sentido.
  - Las capas de vectores b tienen que cumplir Σ_b w_b b_α b_β = δ_αβ. Si no,
    la derivada en k que hay detrás de todo esto está mal y los centros de
    Wannier salen desplazados. Se verifica numéricamente y se reporta el
    residuo.
  - La banda interpolada tiene que reproducir EXACTAMENTE la de DFT en los
    puntos de la malla (es una interpolación, no un ajuste). Se mide.
  - Y tiene que reproducirla APROXIMADAMENTE en puntos que no estaban en la
    malla. Eso es lo único que dice si el modelo sirve, y para saberlo hay
    que correr un cálculo de bandas de verdad y comparar.
"""

from dataclasses import dataclass, field
from pathlib import Path
import itertools
import re

import numpy as np

from qekit.core.errors import ErrorDeUso, FaltanDatos
from qekit.core import style as qstyle

BOHR = 0.529177210903

# Tabla 3.1/3.2 del manual de wannier90: (l, mr) de cada orbital de prueba.
ORBITALES = {
    "s": [(0, 1)],
    "p": [(1, 1), (1, 2), (1, 3)],
    "pz": [(1, 1)], "px": [(1, 2)], "py": [(1, 3)],
    "d": [(2, 1), (2, 2), (2, 3), (2, 4), (2, 5)],
    "dz2": [(2, 1)], "dxz": [(2, 2)], "dyz": [(2, 3)],
    "dx2-y2": [(2, 4)], "dxy": [(2, 5)],
    "f": [(3, m) for m in range(1, 8)],
    "sp": [(-1, 1), (-1, 2)],
    "sp2": [(-2, 1), (-2, 2), (-2, 3)],
    "sp3": [(-3, 1), (-3, 2), (-3, 3), (-3, 4)],
    "sp3d": [(-4, m) for m in range(1, 6)],
    "sp3d2": [(-5, m) for m in range(1, 7)],
}

# Tolerancia con la que se exige Σ_b w_b b_a b_b = δ_ab.
TOL_COMPLETITUD = 1e-5
# Peso por debajo del cual una capa se considera inútil y se descarta.
TOL_PESO = 1e-8
# Diferencia máxima (eV) que se tolera al reproducir la malla de partida.
TOL_EXACTA = 1e-6


# ----------------------------------------------------------------------
# Geometría: la malla completa y las capas de vectores b
# ----------------------------------------------------------------------
def malla_completa(n):
    """Malla uniforme centrada en Γ (sin desplazamiento), SIN reducir por
    simetría.

    El orden es el mismo que usa Quantum ESPRESSO al escribir K_POINTS
    crystal: el índice 3 (el último) corre más rápido. Hay que respetarlo
    porque el archivo .nnkp identifica a los vecinos por su número de orden,
    no por sus coordenadas.
    """
    n1, n2, n3 = (int(x) for x in n)
    if min(n1, n2, n3) < 1:
        raise ErrorDeUso("la malla tiene que ser de al menos 1×1×1.")
    return np.array([[i / n1, j / n2, k / n3]
                     for i in range(n1) for j in range(n2) for k in range(n3)])


def _candidatos(bg, n, nmax=5):
    """Vectores b = (h/n)·bg agrupados por longitud."""
    n = np.asarray(n, float)
    grupos = {}
    for h in itertools.product(*[range(-nmax, nmax + 1)] * 3):
        if h == (0, 0, 0):
            continue
        frac = np.array(h, float) / n
        b = frac @ bg
        grupos.setdefault(round(float(np.linalg.norm(b)), 8), []).append(
            (np.array(h, int), b))
    return [(d, grupos[d]) for d in sorted(grupos)]


def capas_b(bg, n, nmax=5, max_capas=12):
    """Capas de vecinos y sus pesos, tal que Σ_b w_b b_α b_β = δ_αβ.

    Es la condición de completitud de Marzari-Vanderbilt. Sin ella, la
    diferencia finita ∇_k que define el centro de Wannier no es una
    derivada: da un número, pero no el que se busca.

    Se van añadiendo capas de vecinos por orden de distancia y se resuelve
    el sistema de 6 ecuaciones (la parte simétrica de δ_αβ) por mínimos
    cuadrados. Una capa que no aporta rango se descarta —pasa en celdas muy
    anisótropas, donde la primera capa es un plano y no define la tercera
    dirección— y una capa cuyo peso sale nulo también, porque solo añade
    trabajo a pw2wannier90.
    """
    bg = np.asarray(bg, float)
    todas = _candidatos(bg, n, nmax)
    q = np.array([1., 1., 1., 0., 0., 0.])
    elegidas, cols = [], []
    for d, lst in todas[:max_capas]:
        col = np.zeros(6)
        for _, b in lst:
            col += np.array([b[0] * b[0], b[1] * b[1], b[2] * b[2],
                             b[0] * b[1], b[1] * b[2], b[2] * b[0]])
        M = np.array(cols + [col]).T
        s = np.linalg.svd(M, compute_uv=False)
        if s[-1] < 1e-8 * s[0]:      # no añade rango: la capa es redundante
            continue
        cols.append(col)
        elegidas.append((d, lst))
        w, *_ = np.linalg.lstsq(np.array(cols).T, q, rcond=None)
        if np.linalg.norm(np.array(cols).T @ w - q) < TOL_COMPLETITUD:
            # quitar capas de peso nulo y volver a resolver
            keep = [i for i, wi in enumerate(w) if abs(wi) > TOL_PESO]
            if len(keep) < len(w):
                cols = [cols[i] for i in keep]
                elegidas = [elegidas[i] for i in keep]
                w, *_ = np.linalg.lstsq(np.array(cols).T, q, rcond=None)
            return elegidas, w
    raise FaltanDatos(
        "no encuentro un conjunto de capas de vecinos que cumpla la condición "
        "de completitud\ncon esta malla. Suele pasar con mallas muy "
        "anisótropas (por ejemplo 8×8×1 en\nuna celda con mucho vacío): usa "
        "una malla más parecida en las tres direcciones,\no aumenta max_capas.")


def residuo_completitud(capas, pesos):
    """‖Σ_b w_b b⊗b − 1‖_∞. Debe ser cero a precisión de máquina."""
    S = np.zeros((3, 3))
    for w, (_, lst) in zip(pesos, capas):
        for _, b in lst:
            S += w * np.outer(b, b)
    return float(np.abs(S - np.eye(3)).max())


def vecinos(kpts, capas):
    """Para cada k y cada b: qué punto de la malla es k+b y con qué G.

    Devuelve (indices, G) con forma (nk, nnb) y (nk, nnb, 3). Los índices
    son base 1, como los quiere el .nnkp.
    """
    nk = len(kpts)
    clave = {tuple(np.round(k % 1.0, 6)): i for i, k in enumerate(kpts)}
    hs = [h for _, lst in capas for h, _ in lst]
    n = np.array([len(set(np.round(kpts[:, i], 8))) for i in range(3)])
    idx = np.zeros((nk, len(hs)), int)
    G = np.zeros((nk, len(hs), 3), int)
    for ik, k in enumerate(kpts):
        for ib, h in enumerate(hs):
            kb = k + h / n
            red = np.round(kb % 1.0, 6)
            j = clave.get(tuple(red))
            if j is None:                       # empate por redondeo
                d = np.abs((kpts - kb + 0.5) % 1.0 - 0.5).sum(axis=1)
                j = int(np.argmin(d))
            idx[ik, ib] = j + 1
            G[ik, ib] = np.round(kb - kpts[j]).astype(int)
    return idx, G, np.array([b for _, lst in capas for _, b in lst])


def pesos_por_b(capas, pesos):
    """El peso w_b repetido, un valor por vector b (no por capa)."""
    return np.array([w for w, (_, lst) in zip(pesos, capas) for _ in lst])


# ----------------------------------------------------------------------
# Orbitales de prueba
# ----------------------------------------------------------------------
@dataclass
class Proyeccion:
    """Un orbital de prueba: dónde está centrado y de qué tipo es."""
    centro: tuple          # coordenadas fraccionarias
    l: int
    mr: int
    zona: float = 1.0
    r: int = 1
    zaxis: tuple = (0., 0., 1.)
    xaxis: tuple = (1., 0., 0.)
    etiqueta: str = ""


def _orbitales(nombre):
    n = nombre.strip().lower()
    if n in ORBITALES:
        return n, ORBITALES[n]
    raise ErrorDeUso(
        f"no conozco el orbital de prueba '{nombre}'. Los que hay son: "
        f"{', '.join(sorted(ORBITALES))}.")


def proyecciones(spec, atoms, zona=1.0):
    """Traduce 'Si:sp3' o 'f=0.25,0.25,0.25:s' a una lista de Proyeccion.

    Se admiten varias separadas por ';'. El centro puede darse por elemento
    (una proyección por cada átomo de ese elemento) o por coordenadas
    fraccionarias explícitas, que es lo que hace falta para poner un
    orbital en el centro de un enlace —el sitio natural de una función de
    Wannier en un semiconductor tetraédrico—.
    """
    if not spec or str(spec).strip().lower() == "auto":
        return _proyecciones_auto(atoms, zona)
    fuera = []
    simb = atoms.get_chemical_symbols()
    frac = atoms.get_scaled_positions()
    for trozo in str(spec).split(";"):
        trozo = trozo.strip()
        if not trozo:
            continue
        if ":" not in trozo:
            raise ErrorDeUso(
                f"la proyección '{trozo}' no tiene la forma sitio:orbital. "
                f"Ejemplos: Si:sp3, O:p, f=0.5,0.5,0.5:s.")
        sitio, orb = trozo.split(":", 1)
        nom, pares = _orbitales(orb)
        sitio = sitio.strip()
        if sitio.lower().startswith("f="):
            try:
                c = tuple(float(x) for x in sitio[2:].split(","))
            except ValueError:
                raise ErrorDeUso(
                    f"no entiendo las coordenadas de '{sitio}'; se escriben "
                    f"f=x,y,z en fraccionarias.") from None
            if len(c) != 3:
                raise ErrorDeUso(f"'{sitio}' necesita tres coordenadas.")
            centros = [(c, f"f={c[0]:g},{c[1]:g},{c[2]:g}")]
        else:
            centros = [(tuple(frac[i]), f"{sitio}{i + 1}")
                       for i, s in enumerate(simb) if s == sitio]
            if not centros:
                raise ErrorDeUso(
                    f"no hay ningún átomo de {sitio} en la estructura; hay "
                    f"{', '.join(sorted(set(simb)))}.")
        for c, et in centros:
            for l, mr in pares:
                fuera.append(Proyeccion(centro=c, l=l, mr=mr, zona=zona,
                                        etiqueta=f"{et}:{nom}"))
    return fuera


def _proyecciones_auto(atoms, zona=1.0):
    """s y p en cada átomo. Es la apuesta segura, no la buena.

    Cubre la mayoría de los sólidos sp, pero no acierta con metales de
    transición (faltan las d) ni con enlaces muy covalentes, donde las
    funciones de Wannier viven en el enlace y no en el átomo. Por eso
    siempre avisa.
    """
    fuera = []
    frac = atoms.get_scaled_positions()
    for i, (s, c) in enumerate(zip(atoms.get_chemical_symbols(), frac)):
        for l, mr in ORBITALES["s"] + ORBITALES["p"]:
            fuera.append(Proyeccion(centro=tuple(c), l=l, mr=mr, zona=zona,
                                    etiqueta=f"{s}{i + 1}:{'sp'[min(l, 1)]}"))
    return fuera


# ----------------------------------------------------------------------
# El archivo .nnkp, que normalmente escribe 'wannier90.x -pp'
# ----------------------------------------------------------------------
def escribir_nnkp(ruta, atoms, n, proys, excluir=(), nmax=5):
    """Escribe el .nnkp sin necesidad de wannier90.

    pw2wannier90.x solo lee de este archivo la red, los puntos k, los
    orbitales de prueba y la lista de vecinos. Todo eso se puede calcular
    aquí, así que no hace falta tener wannier90 instalado para obtener los
    solapes M y las proyecciones A. Es la diferencia entre poder hacer esto
    en un portátil y no poder.

    El orden de los bloques importa: pw2wannier90 los busca en secuencia y
    no rebobina el archivo entre uno y otro.
    """
    cell = np.array(atoms.cell.array, float)             # Å
    bg = 2 * np.pi * np.linalg.inv(cell).T               # Å⁻¹, con el 2π
    kpts = malla_completa(n)
    capas, pesos = capas_b(bg, n, nmax=nmax)
    idx, G, bvec = vecinos(kpts, capas)
    res = residuo_completitud(capas, pesos)

    L = ["File written by Olla-DFT: entrada para pw2wannier90.x", "",
         "calc_only_A  :  F", "", "begin real_lattice"]
    L += [f"  {v[0]:14.8f}{v[1]:14.8f}{v[2]:14.8f}" for v in cell]
    L += ["end real_lattice", "", "begin recip_lattice"]
    L += [f"  {v[0]:14.8f}{v[1]:14.8f}{v[2]:14.8f}" for v in bg]
    L += ["end recip_lattice", "", "begin kpoints", f"{len(kpts):8d}"]
    L += [f"  {k[0]:14.8f}{k[1]:14.8f}{k[2]:14.8f}" for k in kpts]
    L += ["end kpoints", "", "begin projections", f"{len(proys):8d}"]
    for p in proys:
        L.append(f"  {p.centro[0]:12.8f}{p.centro[1]:12.8f}"
                 f"{p.centro[2]:12.8f}{p.l:4d}{p.mr:3d}{p.r:3d}")
        L.append(f"  {p.zaxis[0]:10.5f}{p.zaxis[1]:10.5f}{p.zaxis[2]:10.5f}"
                 f"{p.xaxis[0]:10.5f}{p.xaxis[1]:10.5f}{p.xaxis[2]:10.5f}"
                 f"{p.zona:12.5f}")
    L += ["end projections", "", "begin nnkpts", f"{idx.shape[1]:8d}"]
    for ik in range(len(kpts)):
        for ib in range(idx.shape[1]):
            g = G[ik, ib]
            L.append(f"{ik + 1:6d}{idx[ik, ib]:6d}"
                     f"{g[0]:5d}{g[1]:5d}{g[2]:5d}")
    L += ["end nnkpts", "", "begin exclude_bands", f"{len(excluir):8d}"]
    L += [f"{int(b):6d}" for b in excluir]
    L += ["end exclude_bands", ""]
    Path(ruta).write_text("\n".join(L), encoding="utf-8")
    return {"nk": len(kpts), "nnb": int(idx.shape[1]), "capas": len(capas),
            "pesos": pesos, "residuo": res, "b": bvec,
            "wb": pesos_por_b(capas, pesos), "kpts": kpts,
            "idx": idx, "G": G}


# ----------------------------------------------------------------------
# Lo que escribe pw2wannier90.x
# ----------------------------------------------------------------------
def leer_eig(ruta):
    """seedname.eig → matriz (nk, nb) de energías en eV."""
    dat = np.loadtxt(ruta)
    if dat.ndim == 1:
        dat = dat[None, :]
    nb, nk = int(dat[:, 0].max()), int(dat[:, 1].max())
    E = np.zeros((nk, nb))
    for ib, ik, e in dat:
        E[int(ik) - 1, int(ib) - 1] = e
    return E


def leer_amn(ruta):
    """seedname.amn → A(nk, nb, nw) = ⟨ψ_mk|g_n⟩."""
    with open(ruta, encoding="utf-8", errors="replace") as f:
        f.readline()
        nb, nk, nw = (int(x) for x in f.readline().split()[:3])
        A = np.zeros((nk, nb, nw), complex)
        for linea in f:
            p = linea.split()
            if len(p) < 5:
                continue
            m, n, ik = int(p[0]), int(p[1]), int(p[2])
            A[ik - 1, m - 1, n - 1] = float(p[3]) + 1j * float(p[4])
    return A


def leer_mmn(ruta):
    """seedname.mmn → M(nk, nnb, nb, nb) y la lista de vecinos.

    M^{k,b}_{mn} = ⟨u_mk|u_{n,k+b}⟩. Es el único ingrediente que hace falta
    para los centros y la dispersión de las funciones de Wannier: toda la
    fase de Berry sale de aquí.
    """
    with open(ruta, encoding="utf-8", errors="replace") as f:
        f.readline()
        nb, nk, nnb = (int(x) for x in f.readline().split()[:3])
        M = np.zeros((nk, nnb, nb, nb), complex)
        vec = np.zeros((nk, nnb, 5), int)
        cont = np.zeros(nk, int)
        for _ in range(nk * nnb):
            cab = f.readline().split()
            ik = int(cab[0]) - 1
            ib = cont[ik]; cont[ik] += 1
            vec[ik, ib] = [int(cab[0]), int(cab[1]),
                           int(cab[2]), int(cab[3]), int(cab[4])]
            blk = np.empty((nb * nb, 2))
            for i in range(nb * nb):
                p = f.readline().split()
                blk[i] = (float(p[0]), float(p[1]))
            z = blk[:, 0] + 1j * blk[:, 1]
            # en el archivo m corre más rápido: la columna n es el índice lento
            M[ik, ib] = z.reshape(nb, nb, order="F")
    return M, vec


def leer_hr(ruta):
    """seedname_hr.dat de wannier90 → (H(R), R, degeneraciones).

    Si tienes wannier90 instalado y ya minimizaste la dispersión, esto entra
    por aquí y el resto del módulo es idéntico. Las energías del archivo
    están en eV.
    """
    lin = Path(ruta).read_text(encoding="utf-8", errors="replace").split("\n")
    nw = int(lin[1].split()[0])
    nr = int(lin[2].split()[0])
    deg, i = [], 3
    while len(deg) < nr:
        deg += [int(x) for x in lin[i].split()]
        i += 1
    deg = np.array(deg[:nr], float)
    R = np.zeros((nr, 3), int)
    H = np.zeros((nr, nw, nw), complex)
    for j in range(nr * nw * nw):
        p = lin[i + j].split()
        ir = j // (nw * nw)
        R[ir] = [int(p[0]), int(p[1]), int(p[2])]
        H[ir, int(p[3]) - 1, int(p[4]) - 1] = float(p[5]) + 1j * float(p[6])
    return H, R, deg


# ----------------------------------------------------------------------
# El gauge: de bandas de Bloch a funciones localizadas
# ----------------------------------------------------------------------
def gauge_proyeccion(A):
    """U = A (A†A)^{-1/2}, ortonormalización de Löwdin.

    A_mn(k) = ⟨ψ_mk|g_n⟩ mide cuánto se parece cada banda a cada orbital de
    prueba. Löwdin es la matriz unitaria más parecida a A en norma de
    Frobenius, así que U(k) es "la proyección, pero ortonormal". Como los
    g_n son los mismos en todos los k, U(k) hereda su suavidad, que es lo
    único que hace falta para que H(R) decaiga.

    Vale igual si hay más bandas que funciones de Wannier: entonces U es
    rectangular (nb×nw) y esto es la desenredada por proyección, sin
    ventanas. Es lo que wannier90 usa como punto de partida.

    Se avisa si algún valor singular de A es pequeño: eso significa que uno
    de los orbitales de prueba no tiene con qué solaparse, y la
    ortonormalización lo va a amplificar hasta convertirlo en ruido.
    """
    A = np.asarray(A, complex)
    U = np.empty((A.shape[0], A.shape[1], A.shape[2]), complex)
    sv_min = np.inf
    for ik in range(A.shape[0]):
        u, s, vh = np.linalg.svd(A[ik], full_matrices=False)
        sv_min = min(sv_min, float(s.min()))
        U[ik] = u @ vh                     # = A (A†A)^{-1/2}
    return U, sv_min


def ventana(E, rango):
    """Máscara de las bandas que caen dentro de una ventana de energía."""
    if rango is None:
        return np.ones(np.asarray(E).shape, bool)
    lo, hi = (float(rango[0]), float(rango[1]))
    if hi <= lo:
        raise ErrorDeUso(
            f"la ventana [{lo}, {hi}] está al revés o es vacía: el segundo "
            f"número tiene que ser mayor.")
    E = np.asarray(E, float)
    return (E >= lo) & (E <= hi)


def omega_I(M, U, idx, wb):
    """Ω_I = (1/N) Σ_kb w_b [J − ‖U†(k) M U(k+b)‖²_F].

    Es la parte de la dispersión que NO depende del gauge: mide cuánto se
    parece el subespacio elegido en k al de sus vecinos. Minimizarla es
    exactamente el problema del desenredado, y es un problema distinto —y
    anterior— al de localizar dentro de un subespacio ya fijo.
    """
    nk, nnb = M.shape[0], M.shape[1]
    J = U[0].shape[1] if isinstance(U, (list, tuple)) else U.shape[2]
    tot = 0.0
    for ik in range(nk):
        Uk = U[ik]
        for ib in range(nnb):
            M2 = Uk.conj().T @ M[ik, ib] @ U[idx[ik, ib] - 1]
            tot += wb[ib] * (J - float(np.sum(np.abs(M2) ** 2)))
    return tot / nk


def gauge_desenredo(M, A, idx, wb, E=None, exterior=None, congelada=None,
                    pasos=200, mezcla=0.5, tol=1e-10, traza=None):
    """Desenredado de Souza-Marzari-Vanderbilt (PRB 65, 035109).

    Cuando las bandas que interesan están ENREDADAS con otras —los estados
    de conducción de casi cualquier sólido, o todo el espectro de un metal—
    no existe un grupo aislado que wannierizar. Hay que elegir, en cada
    punto k, el subespacio de dimensión J que mejor se conecta con el de sus
    vecinos: el que minimiza Ω_I.

    La proyección sola no hace eso. Elige el subespacio que más se parece a
    los orbitales de prueba en CADA k por separado, sin mirar a los vecinos,
    y por eso la interpolación sale mal en cuanto hay enredo.

    El algoritmo es el punto fijo de SMV: en cada k se construye

        Z(k) = Σ_b w_b M^{k,b} P(k+b) M^{k,b}†

    y se toman los J autovectores de mayor autovalor. Si hay ventana
    congelada, esos estados se conservan exactos y solo se eligen los
    J − N_congelados restantes, en el complemento ortogonal.

    Ω_I baja monótonamente. Si no lo hace, es que la mezcla es demasiado
    agresiva; se reduce sola.
    """
    M = np.asarray(M, complex)
    A = np.asarray(A, complex)
    nk, nnb, nb = M.shape[0], M.shape[1], M.shape[2]
    J = A.shape[2]

    if E is None:
        mask = np.ones((nk, nb), bool)
    else:
        mask = ventana(E, exterior)
    if mask.shape != (nk, nb):
        raise ErrorDeUso(
            f"las energías no encajan con los solapes: {mask.shape} contra "
            f"({nk}, {nb}).")
    n_win = mask.sum(axis=1)
    if (n_win < J).any():
        k_malo = int(np.argmin(n_win))
        raise ErrorDeUso(
            f"la ventana exterior deja solo {int(n_win[k_malo])} bandas en "
            f"el punto k {k_malo + 1}, y hacen falta al menos {J}.\n"
            f"Ensánchala, o pide menos funciones de Wannier.")
    frozen = (ventana(E, congelada) & mask if congelada is not None
              else np.zeros((nk, nb), bool))
    n_fro = frozen.sum(axis=1)
    if (n_fro > J).any():
        k_malo = int(np.argmax(n_fro))
        raise ErrorDeUso(
            f"la ventana congelada mete {int(n_fro[k_malo])} bandas en el "
            f"punto k {k_malo + 1}, más que las {J} funciones de Wannier "
            f"que pides.\nEstrecha la ventana congelada.")

    # punto de partida: la proyección, restringida a la ventana
    U = []
    for ik in range(nk):
        w = np.where(mask[ik])[0]
        Ak = A[ik][w, :]
        u, _s, vh = np.linalg.svd(Ak, full_matrices=False)
        Uk = np.zeros((nb, J), complex)
        Uk[w, :] = u @ vh
        U.append(Uk)
    U = np.array(U)

    hist = [omega_I(M, U, idx, wb)]
    Zprev = None
    for it in range(int(pasos)):
        Znew = []
        for ik in range(nk):
            Z = np.zeros((nb, nb), complex)
            for ib in range(nnb):
                Uv = U[idx[ik, ib] - 1]
                Mb = M[ik, ib] @ Uv
                Z += wb[ib] * (Mb @ Mb.conj().T)
            Znew.append(Z)
        Znew = np.array(Znew)
        Z = Znew if Zprev is None else (mezcla * Znew
                                        + (1.0 - mezcla) * Zprev)
        Unew = np.empty_like(U)
        for ik in range(nk):
            w = np.where(mask[ik])[0]
            f = np.where(frozen[ik])[0]
            Uk = np.zeros((nb, J), complex)
            libres = J - len(f)
            if len(f):
                Uk[f, np.arange(len(f))] = 1.0
            if libres > 0:
                Zk = Z[ik][np.ix_(w, w)].copy()
                if len(f):
                    # proyectar fuera lo congelado: esos estados ya están
                    pos = {b: i for i, b in enumerate(w)}
                    P = np.zeros((len(w), len(w)), complex)
                    for b in f:
                        P[pos[b], pos[b]] = 1.0
                    Q = np.eye(len(w)) - P
                    Zk = Q @ Zk @ Q
                val, vec = np.linalg.eigh(0.5 * (Zk + Zk.conj().T))
                orden = np.argsort(val)[::-1][:libres]
                Uk[np.ix_(w, np.arange(len(f), J))] = vec[:, orden]
            Unew[ik] = Uk
        oi = omega_I(M, Unew, idx, wb)
        if oi > hist[-1] + 1e-12:
            mezcla *= 0.5                      # se pasó: paso más corto
            if mezcla < 1e-4:
                break
            continue
        U, Zprev = Unew, Z
        hist.append(oi)
        if abs(hist[-2] - hist[-1]) < tol:
            break
    # El subespacio ya es el óptimo, pero el gauge DENTRO de él es el de los
    # autovectores de Z, que es arbitrario y no es suave. Sin este paso, Ω
    # arranca en decenas de Å² y el descenso posterior no lo arregla: se
    # midió Ω = 51.8 Å² en el silicio, con la valencia peor que sin
    # desenredar. Volver a proyectar sobre los mismos orbitales de prueba,
    # ya dentro del subespacio, devuelve un gauge suave de partida.
    for ik in range(nk):
        Ad = U[ik].conj().T @ A[ik]                 # (J, J)
        u, _sv, vh = np.linalg.svd(Ad, full_matrices=False)
        U[ik] = U[ik] @ (u @ vh)

    if traza:
        np.savetxt(traza, np.column_stack([np.arange(len(hist)), hist]),
                   header="iteracion  Omega_I(A^2)")
    return U, np.array(hist), mask, frozen


def hamiltoniano_k(E, U):
    """H(k) = U†(k) diag(ε(k)) U(k), en la base de Wannier."""
    E = np.asarray(E, float)
    nk, nw = U.shape[0], U.shape[2]
    H = np.empty((nk, nw, nw), complex)
    for ik in range(nk):
        H[ik] = U[ik].conj().T @ (E[ik, :U.shape[1], None] * U[ik])
    return H


def celda_wigner_seitz(cell, n, tol=1e-5):
    """Vectores R de la superred con su degeneración.

    Un R que cae justo en la frontera de la celda de Wigner-Seitz de la
    superred es equidistante de varias imágenes, y hay que repartir su peso
    entre ellas o la interpolación deja de reproducir los puntos de partida.
    Esa degeneración es la diferencia entre una interpolación exacta y una
    que se equivoca en decenas de meV.
    """
    cell = np.asarray(cell, float)
    n = np.asarray(n, int)
    (cell.T * n).T                       # vectores de la superred
    despl = np.array(list(itertools.product(*[range(-2, 3)] * 3)))
    i0 = int(np.where((despl == 0).all(axis=1))[0][0])
    Rs, degs = [], []
    for R in itertools.product(*[range(-int(m), int(m) + 1) for m in n]):
        d = np.array(R, float) - despl * n
        dist = np.einsum("ij,ij->i", d @ cell, d @ cell)
        dmin = dist.min()
        if dist[i0] - dmin < tol:
            Rs.append(R)
            degs.append(int((dist - dmin < tol).sum()))
    return np.array(Rs, int), np.array(degs, float)


def a_reales(Hk, kpts, R):
    """H(R) = (1/N_k) Σ_k e^{-2πi k·R} H(k)."""
    f = np.exp(-2j * np.pi * (np.asarray(kpts) @ np.asarray(R).T))  # (nk, nR)
    return np.einsum("kr,kmn->rmn", f, Hk) / len(kpts)


def interpolar(HR, R, deg, kpts, vectores=False):
    """Diagonaliza H(k) = Σ_R e^{2πi k·R} H(R)/deg(R) en los k que le des."""
    kpts = np.atleast_2d(np.asarray(kpts, float))
    f = np.exp(2j * np.pi * (kpts @ np.asarray(R).T)) / np.asarray(deg)
    Hk = np.einsum("kr,rmn->kmn", f, HR)
    Hk = 0.5 * (Hk + np.conj(np.transpose(Hk, (0, 2, 1))))
    if vectores:
        w, v = np.linalg.eigh(Hk)
        return w, v
    return np.linalg.eigvalsh(Hk)


def decaimiento(HR, R, cell):
    """max|H_mn(R)| frente a |R|. Es el diagnóstico de si esto va a servir.

    Si el máximo no ha bajado un par de órdenes de magnitud al llegar al
    borde de la superred, la base no está localizada y la interpolación va
    a inventar estructura entre los puntos de la malla.
    """
    d = np.linalg.norm(np.asarray(R, float) @ np.asarray(cell, float), axis=1)
    amp = np.abs(HR).max(axis=(1, 2))
    orden = np.argsort(d)
    return d[orden], amp[orden]


# ----------------------------------------------------------------------
# Centros, dispersión y minimización (Marzari-Vanderbilt)
# ----------------------------------------------------------------------
def _m_gauge(M, U, idx):
    """M en la base de Wannier: M^W = U†(k) M^{k,b} U(k+b)."""
    nk, nnb = M.shape[0], M.shape[1]
    nw = U.shape[2]
    MW = np.empty((nk, nnb, nw, nw), complex)
    for ik in range(nk):
        for ib in range(nnb):
            MW[ik, ib] = U[ik].conj().T @ M[ik, ib] @ U[idx[ik, ib] - 1]
    return MW


@dataclass
class Dispersion:
    """El funcional de dispersión troceado como lo trocea wannier90."""
    centros: np.ndarray          # (nw, 3) en Å
    spreads: np.ndarray          # (nw,) en Å²
    omega: float
    omega_I: float               # invariante: no baja al minimizar
    omega_D: float               # diagonal
    omega_OD: float              # fuera de la diagonal
    deriva_I: float = 0.0        # cuánto se movió Ω_I al minimizar (debe ser ~0)


def dispersion(MW, b, wb):
    """Centros y dispersión de cada función de Wannier.

    Marzari-Vanderbilt (PRB 56, 12847), ecuaciones 31 y 34-36. Todo sale de
    la fase de los elementos diagonales de M: Im ln M_nn es literalmente
    una fase de Berry, y su derivada respecto de k es la posición.

    Es el mismo objeto que aparece en la polarización eléctrica de
    Berry-King-Smith. Por eso el centro de Wannier tiene unidades de
    longitud y no de "algo parecido a una posición": lo es.
    """
    nk, nnb, nw, _ = MW.shape
    diag = np.einsum("kbnn->kbn", MW)
    fase = np.angle(diag)                                # Im ln M_nn
    # r̄_n = −(1/N) Σ_kb w_b b · Im ln M_nn
    r = -np.einsum("b,bx,kbn->nx", wb, b, fase) / nk
    r2 = np.einsum("b,kbn->n", wb, 1.0 - np.abs(diag) ** 2 + fase ** 2) / nk
    spreads = r2 - np.einsum("nx,nx->n", r, r)
    absM2 = np.abs(MW) ** 2
    om_I = float(np.einsum("b,kb->", wb, nw - absM2.sum(axis=(2, 3))) / nk)
    om_OD = float(np.einsum("b,kb->", wb,
                            absM2.sum(axis=(2, 3))
                            - (np.abs(diag) ** 2).sum(axis=2)) / nk)
    br = b @ r.T                                          # (nnb, nw)
    om_D = float(np.einsum("b,kbn->", wb, (-fase - br[None, :, :]) ** 2) / nk)
    return Dispersion(centros=r, spreads=spreads, omega=float(spreads.sum()),
                      omega_I=om_I, omega_D=om_D, omega_OD=om_OD)


def _gradiente(MW, b, wb):
    """dΩ/dW en cada k, con W una rotación antihermítica de U(k).

    Marzari-Vanderbilt ec. 52-57:

        G^(k) = −(4/N_k) Σ_b w_b [ A(R^{k,b}) − S(T^{k,b}) ]

    con R_mn = M_mn M*_nn,  T_mn = (M_mn/M_nn)·q_n,  q_n = Im ln M_nn + b·r̄_n,
    A(B) = (B − B†)/2 y S(B) = (B + B†)/2i.

    El 1/N_k y el signo NO son cosmética: se fijaron comparando contra la
    derivada numérica de Ω, que da exactamente −1/N_k veces la expresión sin
    normalizar. Con el factor equivocado el paso sale 64 veces demasiado
    largo en una malla 4×4×4 y el descenso da tumbos en vez de bajar.
    """
    nk, nnb, nw, _ = MW.shape
    diag = np.einsum("kbnn->kbn", MW)
    fase = np.angle(diag)
    r = -np.einsum("b,bx,kbn->nx", wb, b, fase) / nk
    q = fase + (b @ r.T)[None, :, :]
    G = np.zeros((nk, nw, nw), complex)
    for ib in range(nnb):
        Mb, dn = MW[:, ib], diag[:, ib]
        Rm = Mb * dn.conj()[:, None, :]
        # M_nn puede acercarse a cero con mallas gruesas; sin este suelo el
        # gradiente explota.
        seguro = np.where(np.abs(dn) < 1e-10, 1e-10 + 0j, dn)
        Tm = (Mb / seguro[:, None, :]) * q[:, ib][:, None, :]
        A = 0.5 * (Rm - np.conj(np.transpose(Rm, (0, 2, 1))))
        S = 0.5 / 1j * (Tm + np.conj(np.transpose(Tm, (0, 2, 1))))
        G += 4.0 * wb[ib] * (A - S)
    return -G / nk


def _rotar(U, W_):
    """U ← U·exp(W), con W antihermítica, por diagonalización exacta."""
    Un = np.empty_like(U)
    for ik in range(U.shape[0]):
        w_, v_ = np.linalg.eigh(1j * W_[ik])
        Un[ik] = U[ik] @ ((v_ * np.exp(-1j * w_)) @ v_.conj().T)
    return Un


def minimizar(M, U, idx, b, wb, pasos=500, alfa=2.0, tol=1e-10,
              retroceso=12, traza=None):
    """Descenso por gradiente sobre Ω: de la proyección a maximalmente localizada.

    La dirección de descenso es −G. El paso de prueba es el de wannier90,
    Δt = α/(4 Σ_b w_b), y si el paso sube Ω se parte por la mitad hasta que
    baje (hasta `retroceso` veces). Esa búsqueda de línea es lo que hace que
    esto funcione también cuando el punto de partida es malo; sin ella, un
    gauge inicial desordenado hace que el método se pasee sin bajar.

    Solo bajan Ω_D y Ω_OD. Ω_I es invariante de gauge y tiene que quedarse
    quieto: se comprueba y se devuelve, porque si se mueve es que hay un
    error, no que el cálculo esté "casi convergido".
    """
    U = np.array(U, complex, copy=True)
    dt0 = float(alfa) / (4.0 * float(np.sum(wb)))
    d = dispersion(_m_gauge(M, U, idx), b, wb)
    hist = [d.omega]
    omI0 = d.omega_I
    for _ in range(int(pasos)):
        MW = _m_gauge(M, U, idx)
        G = _gradiente(MW, b, wb)
        dt = dt0
        for _ in range(int(retroceso)):
            Un = _rotar(U, -dt * G)
            dn = dispersion(_m_gauge(M, Un, idx), b, wb)
            if dn.omega <= hist[-1]:
                break
            dt *= 0.5
        else:
            break                       # ni partiendo el paso baja: parar
        U, d = Un, dn
        hist.append(d.omega)
        if abs(hist[-2] - hist[-1]) < tol:
            break
    if traza:
        np.savetxt(traza, np.column_stack([np.arange(len(hist)), hist]),
                   header="iteracion  Omega(A^2)")
    d.deriva_I = abs(d.omega_I - omI0)
    return U, d, np.array(hist)


# ----------------------------------------------------------------------
# El .win, para quien sí tenga wannier90
# ----------------------------------------------------------------------
def escribir_win(ruta, atoms, n, proys, nbnd, nwann, camino=None,
                 iteraciones=200, excluir=()):
    """Archivo de entrada de wannier90, por si prefieres su minimización.

    Olla-DFT no lo necesita —hace la proyección y la minimización él mismo—
    pero escribirlo cuesta nada y deja la puerta abierta: si tienes
    wannier90 instalado puedes correr `wannier90.x -pp seed`, `pw2wannier90.x`
    y `wannier90.x seed`, y luego leer el seedname_hr.dat con `leer_hr()`.
    Los dos caminos desembocan en el mismo H(R).
    """
    cell = np.array(atoms.cell.array, float)
    L = [f"num_wann        = {int(nwann)}",
         f"num_bands       = {int(nbnd)}",
         f"num_iter        = {int(iteraciones)}",
         "dis_num_iter    = 0" if nbnd == nwann else "dis_num_iter    = 200",
         "write_hr        = .true.",
         "write_xyz       = .true.", ""]
    if excluir:
        L.append("exclude_bands   = " + ",".join(str(int(b)) for b in excluir))
        L.append("")
    L += ["begin unit_cell_cart", "ang"]
    L += [f"  {v[0]:14.8f}{v[1]:14.8f}{v[2]:14.8f}" for v in cell]
    L += ["end unit_cell_cart", "", "begin atoms_frac"]
    for s, f in zip(atoms.get_chemical_symbols(), atoms.get_scaled_positions()):
        L.append(f"  {s:3s}{f[0]:14.8f}{f[1]:14.8f}{f[2]:14.8f}")
    L += ["end atoms_frac", "", "begin projections"]
    for p in proys:
        L.append(f"  f={p.centro[0]:.6f},{p.centro[1]:.6f},{p.centro[2]:.6f}"
                 f":l={p.l},mr={p.mr}")
    L += ["end projections", ""]
    if camino:
        L += ["bands_plot      = .true.", "begin kpoint_path"]
        for (a, ka), (b, kb) in zip(camino[:-1], camino[1:]):
            L.append(f"  {a} {ka[0]:.5f} {ka[1]:.5f} {ka[2]:.5f}   "
                     f"{b} {kb[0]:.5f} {kb[1]:.5f} {kb[2]:.5f}")
        L += ["end kpoint_path", ""]
    L += [f"mp_grid : {int(n[0])} {int(n[1])} {int(n[2])}", "", "begin kpoints"]
    L += [f"  {k[0]:14.8f}{k[1]:14.8f}{k[2]:14.8f}" for k in malla_completa(n)]
    L += ["end kpoints", ""]
    Path(ruta).write_text("\n".join(L), encoding="utf-8")
    return ruta


# ----------------------------------------------------------------------
# Preparar los tres pasos
# ----------------------------------------------------------------------
@dataclass
class WannierRun:
    """Todo lo que sale de un cálculo de funciones de Wannier."""
    formula: str = ""
    malla: tuple = (4, 4, 4)
    nk: int = 0
    nnb: int = 0
    capas: int = 0
    residuo: float = 0.0
    nbnd: int = 0
    nwann: int = 0
    excluir: tuple = ()
    proyecciones: list = field(default_factory=list)
    fuente: str = "proyección"          # o "wannier90 (_hr.dat)"
    sv_min: float = float("nan")
    exterior: tuple = None              # ventana de desenredado
    congelada: tuple = None
    omega_I_hist: np.ndarray = None     # bajada de Ω_I al desenredar
    n_congeladas: tuple = None          # mínimo y máximo por punto k
    disp0: object = None                # dispersión antes de minimizar
    disp: object = None                 # y después
    hist: np.ndarray = None
    HR: np.ndarray = None
    R: np.ndarray = None
    deg: np.ndarray = None
    cell: np.ndarray = None
    error_malla: float = float("nan")   # debe ser ~0: es interpolación
    camino: list = field(default_factory=list)
    k_camino: np.ndarray = None
    E_wann: np.ndarray = None
    E_dft: np.ndarray = None
    E_sin_gauge: np.ndarray = None
    fermi: float = None
    avisos: list = field(default_factory=list)


def prepare(atoms, outdir: str = "wannier", malla=(4, 4, 4), proy="auto",
            nbnd: int = None, excluir=(), pseudo_dir: str = None,
            ecutwfc: float = None, ecutrho: float = None, kgrid_scf=None,
            insulator: bool = False, camino=None, iteraciones: int = 200,
            zona: float = 1.0):
    """Escribe el scf, el nscf de malla completa, el .nnkp y el pw2wannier90.

    Los dos detalles que arruinan este cálculo si se hacen mal:

    1. **El nscf va con nosym y noinv.** wannier90 y pw2wannier90 necesitan
       TODOS los puntos de la malla, no la cuña irreducible. Con la malla
       reducida el número de puntos no coincide con el del .nnkp y
       pw2wannier90 se para; si por casualidad coincidiera, los solapes
       serían entre estados que no son vecinos y no lo diría nadie.
    2. **nbnd tiene que ser al menos num_wann.** Si pides ocho funciones de
       Wannier de cuatro bandas, la proyección no tiene de dónde salir.
    """
    from qekit.modules import inputgen, sweep

    proys = proyecciones(proy, atoms, zona=zona)
    nw = len(proys)
    excluir = tuple(int(b) for b in excluir)
    nb_total = int(nbnd) if nbnd else nw + len(excluir)
    if nb_total - len(excluir) < nw:
        raise ErrorDeUso(
            f"pides {nw} funciones de Wannier pero solo quedan "
            f"{nb_total - len(excluir)} bandas después de excluir "
            f"{len(excluir)}. Sube --bands o quita proyecciones.")

    common = sweep.prepare_common(atoms, pseudo_dir, ecutwfc, ecutrho,
                                  insulator)
    out = Path(outdir); out.mkdir(parents=True, exist_ok=True)
    k_scf = tuple(kgrid_scf or sweep.default_grid(atoms, None))
    kfull = malla_completa(malla)
    kcard = ("K_POINTS crystal\n" + f"{len(kfull)}\n" +
             "".join(f" {a:.10f} {b:.10f} {c:.10f} {1.0 / len(kfull):.10f}\n"
                     for a, b, c in kfull))
    for nombre, calc, card, extra in (
            ("1_scf", "scf",
             f"K_POINTS automatic\n  {k_scf[0]} {k_scf[1]} {k_scf[2]} "
             "0 0 0\n", {}),
            ("2_nscf", "nscf", kcard, {"nosym": True})):
        txt = inputgen.build_pw_input(
            atoms=atoms, pseudos=common["pseudos"], calculation=calc,
            prefix=common["prefix"], pseudo_dir=common["pseudo_dir"],
            ecutwfc=common["ecutwfc"], ecutrho=common["ecutrho"],
            kcard=card, insulator=common["insulator"],
            degauss=common["degauss"], smearing=common["smearing"],
            conv_thr=1e-10, nbnd=nb_total if calc == "nscf" else None, **extra)
        sweep.write_input(out / f"{nombre}.in", txt)

    seed = common["prefix"]
    info = escribir_nnkp(out / f"{seed}.nnkp", atoms, malla, proys,
                         excluir=excluir)
    sweep.write_input(out / "3_pw2wan.in",
                      "&inputpp\n"
                      "  outdir    = './out'\n"
                      f"  prefix    = '{seed}'\n"
                      f"  seedname  = '{seed}'\n"
                      "  write_amn = .true.\n"
                      "  write_mmn = .true.\n"
                      "  write_unk = .false.\n/\n")
    escribir_win(out / f"{seed}.win", atoms, malla, proys, nb_total, nw,
                 camino=camino, iteraciones=iteraciones, excluir=excluir)

    # 4º paso: las bandas de DFT sobre un camino, que es la ÚNICA forma de
    # saber si el modelo interpolado sirve. Sin esto solo se puede comprobar
    # que reproduce los puntos que ya conocía, que es trivial.
    ks, _, _, _ = camino_denso(atoms, 30)
    txt = inputgen.build_pw_input(
        atoms=atoms, pseudos=common["pseudos"], calculation="bands",
        prefix=common["prefix"], pseudo_dir=common["pseudo_dir"],
        ecutwfc=common["ecutwfc"], ecutrho=common["ecutrho"],
        kcard="K_POINTS crystal\n" + f"{len(ks)}\n" +
              "".join(f" {a:.10f} {b:.10f} {c:.10f} 1.0\n" for a, b, c in ks),
        insulator=common["insulator"], degauss=common["degauss"],
        smearing=common["smearing"], nbnd=nb_total)
    txt = txt.replace("outdir", "outdir", 1)
    sweep.write_input(out / "4_bands.in", txt.replace("'./out'", "'./out_bandas'"))

    from qekit.core import plataforma
    plataforma.escribir_par_de_guiones(out, [
        ("pw.x", "1_scf.in", "1_scf.out"),
        ("pw.x", "2_nscf.in", "2_nscf.out"),
        ("pw2wannier90.x", "3_pw2wan.in", "3_pw2wan.out"),
    ])

    run = WannierRun(formula=atoms.get_chemical_formula(),
                     malla=tuple(int(m) for m in malla), nk=info["nk"],
                     nnb=info["nnb"], capas=info["capas"],
                     residuo=info["residuo"], nbnd=nb_total, nwann=nw,
                     excluir=excluir,
                     proyecciones=[p.etiqueta or f"l={p.l},mr={p.mr}"
                                   for p in proys],
                     cell=np.array(atoms.cell.array, float))
    if str(proy).strip().lower() in ("", "auto"):
        run.avisos.append(
            "Proyecciones automáticas (s y p en cada átomo). Es una apuesta "
            "razonable en\n  sólidos sp y una mala idea en metales de "
            "transición (faltan las d) y en enlaces\n  muy covalentes, donde "
            "la función de Wannier vive en el enlace y no en el átomo.\n"
            "  Si la dispersión sale grande o la interpolación mala, empieza "
            "por aquí.")
    rep = [f"--- Funciones de Wannier: {run.formula} ---",
           f"Malla completa: {malla[0]}×{malla[1]}×{malla[2]} = {info['nk']} "
           f"puntos k (sin reducir por simetría)",
           f"Vecinos por punto k: {info['nnb']} en {info['capas']} capa(s); "
           f"residuo de completitud {info['residuo']:.1e}",
           f"Bandas: {nb_total}   Funciones de Wannier: {nw}"
           + (f"   Excluidas: {','.join(map(str, excluir))}" if excluir else ""),
           "",
           f"Archivos en '{out.resolve()}':",
           "  1_scf.in      scf normal",
           "  2_nscf.in     malla COMPLETA, con nosym y noinv",
           f"  {seed}.nnkp   vecinos y orbitales de prueba (lo que "
           "normalmente escribe wannier90.x -pp)",
           "  3_pw2wan.in   pw2wannier90.x: escribe .amn, .mmn y .eig",
           "  4_bands.in    bandas de DFT sobre el camino, para comparar",
           f"  {seed}.win    por si prefieres correr wannier90 tú mismo",
           "",
           "  bash correr.sh    (o los tres a mano, en ese orden)",
           "",
           "Luego:  olla-dft wannier <estructura> --collect -o "
           f"{out}",
           "",
           "No hace falta tener wannier90 instalado: el .nnkp lo escribe "
           "Olla-DFT y la\nlocalización se hace aquí. Si lo tienes y prefieres "
           "su minimización, corre\nwannier90.x y Olla-DFT leerá su "
           "seedname_hr.dat."]
    warn = sweep.missing_pseudo_warning(common)
    if warn:
        rep.append(warn)
    return run, common, "\n".join(rep)


def camino_denso(atoms, puntos_por_tramo: int = 30):
    """Camino de alta simetría muestreado, con etiquetas y abscisa."""
    from qekit.core import kpoints as kp
    kpath = kp.get_kpath(atoms)
    cell = np.array(atoms.cell.array, float)
    bg = 2 * np.pi * np.linalg.inv(cell).T
    ks, etiquetas, x = [], [], [0.0]
    for i, (a, b) in enumerate(kpath.path):
        ka = np.array(kpath.point_coords[a], float)
        kb = np.array(kpath.point_coords[b], float)
        if not ks:
            ks.append(ka); etiquetas.append((0, kp.pretty_label(a)))
        elif np.linalg.norm(ka - ks[-1]) > 1e-8:
            etiquetas.append((len(ks) - 1, etiquetas[-1][1] + "|"
                              + kp.pretty_label(a)))
            ks.append(ka); x.append(x[-1])
        for t in np.linspace(0, 1, puntos_por_tramo + 1)[1:]:
            k = ka + t * (kb - ka)
            x.append(x[-1] + np.linalg.norm((k - ks[-1]) @ bg))
            ks.append(k)
        etiquetas.append((len(ks) - 1, kp.pretty_label(b)))
    # quitar etiquetas repetidas seguidas
    limpio = []
    for idx, et in etiquetas:
        if limpio and limpio[-1][0] == idx:
            limpio[-1] = (idx, et if "|" in et else limpio[-1][1])
        else:
            limpio.append((idx, et))
    return np.array(ks), np.array(x[:len(ks)]), limpio, kpath


def _leer_nnkp(ruta):
    """Recupera celda, malla y bandas excluidas del .nnkp que se escribió."""
    txt = Path(ruta).read_text(encoding="utf-8", errors="replace")

    def bloque(nombre):
        m = re.search(rf"begin\s+{nombre}\s*\n(.*?)\n\s*end\s+{nombre}",
                      txt, re.S | re.I)
        return m.group(1).strip().split("\n") if m else []

    cell = np.array([[float(x) for x in l.split()]
                     for l in bloque("real_lattice")])
    kl = bloque("kpoints")
    kpts = np.array([[float(x) for x in l.split()] for l in kl[1:]])
    n = tuple(int(round(1.0 / min(v for v in np.unique(np.round(kpts[:, i], 8))
                                  if v > 1e-8))) if
              len(np.unique(np.round(kpts[:, i], 8))) > 1 else 1
              for i in range(3))
    ex = bloque("exclude_bands")
    excluir = tuple(int(x) for x in ex[1:]) if len(ex) > 1 else ()
    return cell, n, kpts, excluir


def collect(outdir: str = "wannier", minimizar_=True, pasos: int = 500,
            bandas_dft=None, puntos_por_tramo: int = 30, atoms=None,
            exterior=None, congelada=None):
    """Lee lo que dejó pw2wannier90 (o wannier90) y construye el modelo.

    El orden importa: primero se comprueba que la interpolación reproduce
    EXACTAMENTE los puntos de la malla —si no, hay un error de índices o de
    degeneraciones y todo lo demás sobra— y solo después se mira si
    reproduce los puntos que no estaban.
    """
    out = Path(outdir)
    nnkp = sorted(out.glob("*.nnkp"))
    # OJO: el H(R) que exporta Olla-DFT tiene el mismo nombre que el de
    # wannier90 a propósito (para que sea intercambiable), así que hay que
    # distinguirlos por la cabecera. Sin esto, un segundo --collect leía su
    # propio resultado en vez de recalcularlo, y todos los diagnósticos
    # —dispersión, centros, exactitud— desaparecían del informe sin avisar.
    hr = [f for f in sorted(out.glob("*_hr.dat"))
          if "Olla-DFT" not in f.read_text(errors="replace").split("\n")[0]]
    if not nnkp and not hr:
        raise FaltanDatos(
            f"en {out} no hay ni un .nnkp ni un seedname_hr.dat. Corre antes "
            f"`olla-dft wannier ... -o {out}` y luego los tres pasos de "
            f"correr.sh.")
    run = WannierRun()

    if nnkp:
        seed = nnkp[0].stem
        cell, n, kpts, excluir = _leer_nnkp(nnkp[0])
        run.malla, run.excluir, run.nk, run.cell = n, excluir, len(kpts), cell
    else:
        seed = hr[0].name.replace("_hr.dat", "")
        cell = np.array(atoms.cell.array, float) if atoms is not None else None
        n, kpts = None, None

    if hr:
        # wannier90 ya hizo el trabajo: nos quedamos con su H(R)
        run.HR, run.R, run.deg = leer_hr(hr[0])
        run.fuente = f"wannier90 ({hr[0].name})"
        run.nwann = run.HR.shape[1]
        if cell is None:
            raise FaltanDatos(
                "para usar un _hr.dat hace falta la estructura: pásala como "
                "primer argumento.")
        run.cell = cell
        if n is None:
            n = tuple(int(abs(run.R[:, i]).max()) + 1 for i in range(3))
            run.malla = n
    else:
        for suf in (".amn", ".mmn", ".eig"):
            if not (out / (seed + suf)).exists():
                raise FaltanDatos(
                    f"falta {seed}{suf}. Es lo que escribe pw2wannier90.x: "
                    f"corre el paso 3 (`pw2wannier90.x -in 3_pw2wan.in`) "
                    f"dentro de {out}.")
        E = leer_eig(out / (seed + ".eig"))
        A = leer_amn(out / (seed + ".amn"))
        M, _ = leer_mmn(out / (seed + ".mmn"))
        run.nbnd, run.nwann = A.shape[1], A.shape[2]
        bg = 2 * np.pi * np.linalg.inv(cell).T
        capas, pes = capas_b(bg, n)
        idx, _, bvec = vecinos(kpts, capas)
        wb = pesos_por_b(capas, pes)
        run.nnb, run.capas = len(bvec), len(capas)
        run.residuo = residuo_completitud(capas, pes)
        run.exterior, run.congelada = exterior, congelada
        if exterior is not None or congelada is not None or A.shape[1] > A.shape[2]:
            # hay bandas de sobra: el subespacio hay que ELEGIRLO, no solo
            # proyectar. Con más bandas que funciones de Wannier, proyectar
            # sin más elige en cada k por separado y la interpolación sale
            # mal en cuanto las bandas están enredadas.
            U, run.omega_I_hist, _mask, fro = gauge_desenredo(
                M, A, idx, wb, E=E, exterior=exterior, congelada=congelada,
                traza=str(out / "omega_I.dat"))
            nf = fro.sum(axis=1)
            run.n_congeladas = (int(nf.min()), int(nf.max()))
            run.sv_min = float("nan")
        else:
            fro = np.zeros(A.shape[:2], bool)
            U, run.sv_min = gauge_proyeccion(A)
        run.disp0 = dispersion(_m_gauge(M, U, idx), bvec, wb)
        if minimizar_:
            U, run.disp, run.hist = minimizar(M, U, idx, bvec, wb, pasos=pasos,
                                              traza=str(out / "omega.dat"))
            run.fuente = "proyección + minimización (Olla-DFT)"
        else:
            run.disp, run.hist = run.disp0, np.array([run.disp0.omega])
            run.fuente = "proyección (sin minimizar)"
        run.R, run.deg = celda_wigner_seitz(cell, n)
        run.HR = a_reales(hamiltoniano_k(E, U), kpts, run.R)
        Ei = np.sort(interpolar(run.HR, run.R, run.deg, kpts), axis=1)
        if run.omega_I_hist is None:
            # manifold aislado: la interpolación reproduce TODAS las bandas
            run.error_malla = float(np.abs(Ei - np.sort(E, axis=1)).max())
        elif run.n_congeladas and run.n_congeladas[0] == run.n_congeladas[1] \
                and run.n_congeladas[0] > 0:
            # con desenredado solo las bandas CONGELADAS tienen que salir
            # exactas: el resto del subespacio es una mezcla elegida por
            # suavidad, no un conjunto de bandas
            nf = run.n_congeladas[0]
            Ec = np.sort(np.where(fro, E, np.inf), axis=1)[:, :nf]
            run.error_malla = float(np.abs(Ei[:, :nf] - Ec).max())
        else:
            # sin ventana congelada no hay nada que tenga que salir exacto
            run.error_malla = float("nan")
            run.avisos.append(
                "Sin ventana congelada no hay ninguna banda que la "
                "interpolación tenga que\n  reproducir exactamente: el "
                "subespacio se eligió por suavidad, no por bandas. Si "
                "quieres\n  que la valencia salga exacta, pásala en "
                "--frozen.")
        # control negativo: lo mismo sin gauge, para poder enseñar la diferencia
        if run.omega_I_hist is None:
            Uid = np.tile(np.eye(A.shape[1], A.shape[2], dtype=complex),
                          (len(kpts), 1, 1))
            run._HR_sin = a_reales(hamiltoniano_k(E, Uid), kpts, run.R)

    if atoms is not None:
        run.formula = atoms.get_chemical_formula()
        ks, x, et, kpath = camino_denso(atoms, puntos_por_tramo)
        if kpath.cell_changed:
            run.avisos.append(
                "La celda que usaste no es la primitiva estándar de seekpath, "
                "así que las\n  etiquetas del camino pueden no corresponder a "
                "los puntos que nombran. Pasa la\n  celda primitiva "
                "(`olla-dft prim`) si quieres el camino canónico.")
        run.k_camino, run.camino = ks, et
        run._x = x
        run.E_wann = np.sort(interpolar(run.HR, run.R, run.deg, ks), axis=1)
        if getattr(run, "_HR_sin", None) is not None:
            run.E_sin_gauge = np.sort(
                interpolar(run._HR_sin, run.R, run.deg, ks), axis=1)

    if bandas_dft:
        from qekit.core import qeout
        res = qeout.read_xml(str(bandas_dft))
        Eb = np.array(res.eigenvalues)[0]
        kb = np.array(res.kpoints_frac)
        run.fermi = res.fermi
        if run.k_camino is None or len(kb) != len(run.k_camino):
            run.k_camino = kb
            run._x = np.concatenate([[0.0], np.cumsum(np.linalg.norm(
                np.diff(kb, axis=0) @ (2 * np.pi * np.linalg.inv(run.cell).T),
                axis=1))])
            run.camino = []
            run.E_wann = np.sort(
                interpolar(run.HR, run.R, run.deg, kb), axis=1)
            if getattr(run, "_HR_sin", None) is not None:
                run.E_sin_gauge = np.sort(
                    interpolar(run._HR_sin, run.R, run.deg, kb), axis=1)
        nb = run.E_wann.shape[1]
        run.E_dft = np.sort(Eb[:, :nb], axis=1)
    return run


def asignar(centros, atoms, d_enlace=(0.5, 3.2)):
    """¿Dónde vive cada función de Wannier: sobre un átomo o sobre un enlace?

    ``d_enlace`` es la ventana de distancias (Å) dentro de la cual dos
    átomos se consideran enlazados y su punto medio entra como candidato.

    Es el diagnóstico más rápido de si el modelo tiene sentido. En un
    semiconductor covalente los centros deben caer en los enlaces; en un
    óxido iónico, sobre el oxígeno; en un metal de transición, sobre el
    metal. Un centro que no cae en ningún sitio reconocible casi siempre
    significa que las proyecciones estaban mal elegidas.
    """
    cell = np.array(atoms.cell.array, float)
    pos = atoms.get_positions()
    sim = atoms.get_chemical_symbols()
    inv = np.linalg.inv(cell)
    sitios = [(p, f"{s}{i + 1}") for i, (p, s) in enumerate(zip(pos, sim))]
    # puntos medios de los pares más cercanos (candidatos a enlace)
    for i in range(len(pos)):
        for j in range(len(pos)):
            for T in itertools.product((-1, 0, 1), repeat=3):
                pj = pos[j] + np.array(T) @ cell
                d = np.linalg.norm(pj - pos[i])
                if d_enlace[0] < d < d_enlace[1] and i <= j:
                    sitios.append(((pos[i] + pj) / 2,
                                   f"{sim[i]}{i + 1}–{sim[j]}{j + 1}"))
    fuera = []
    for c in np.atleast_2d(centros):
        mejor, dmin = "?", np.inf
        for p, et in sitios:
            df = (c - p) @ inv
            df -= np.round(df)
            d = float(np.linalg.norm(df @ cell))
            if d < dmin:
                dmin, mejor = d, et
        fuera.append((mejor, dmin))
    return fuera


# Unidades de la DOS interpolada, tal como van a la cabecera del archivo.
# Sin factor de espín: una función de Wannier = un estado, integra a num_wann.
DOS_UNIDADES = ("estados/eV/celda, sin factor de espín: integra a num_wann "
                "(x2 para comparar con dos.x sin espín)")


def dos_interpolada(run, malla=24, sigma: float = 0.05, npuntos: int = 601,
                    rango=None):
    """Densidad de estados en una malla densa, gratis, desde H(R).

    Es el rendimiento inmediato de haber wannierizado: la malla de DFT era
    de 4³ u 8³ puntos y aquí se usan 24³ = 13 824 sin volver a tocar pw.x.
    Solo vale dentro del rango de energías que cubren las funciones de
    Wannier; fuera de ahí la curva es una invención y se recorta.

    Normalización (ver DOS_UNIDADES): estados por eV y por celda, SIN el
    factor 2 de espín: cada función de Wannier aporta exactamente un
    estado, así que ∫ DOS dE = num_wann. Es lo que hace wannier90 con su
    propia DOS y lo que comprueba el test de integración. Para comparar
    con la DOS de dos.x de un cálculo sin espín (que sí lleva el 2 de
    ocupación) hay que multiplicar por 2 — y por eso la cabecera del archivo
    lo dice en vez de poner solo "estados/eV/celda".
    """
    n = (int(malla),) * 3 if np.isscalar(malla) else tuple(int(m) for m in malla)
    ks = malla_completa(n)
    E = interpolar(run.HR, run.R, run.deg, ks)
    lo, hi = (float(E.min()), float(E.max())) if rango is None else rango
    ejes = np.linspace(lo - 5 * sigma, hi + 5 * sigma, int(npuntos))
    dos = np.zeros_like(ejes)
    pref = 1.0 / (sigma * np.sqrt(2 * np.pi) * len(ks))
    for e in E.ravel():
        dos += np.exp(-0.5 * ((ejes - e) / sigma) ** 2)
    return ejes, dos * pref


def report(run, atoms=None) -> str:
    L = [f"--- Funciones de Wannier: {run.formula or ''} ---".strip(),
         f"Fuente: {run.fuente}",
         f"Malla {run.malla[0]}×{run.malla[1]}×{run.malla[2]} = {run.nk} "
         f"puntos k   |   {run.nwann} funciones de Wannier"
         + (f" de {run.nbnd} bandas" if run.nbnd else "")]
    if run.excluir:
        L.append(f"Bandas excluidas: {','.join(map(str, run.excluir))}")
    if run.nnb:
        L.append(f"Vecinos: {run.nnb} en {run.capas} capa(s)   "
                 f"|   Σ w_b b⊗b − 1 = {run.residuo:.1e}"
                 + ("  ✓" if run.residuo < TOL_COMPLETITUD else "  ← MAL"))
    if run.omega_I_hist is not None:
        L += ["", "Desenredado (Souza-Marzari-Vanderbilt):",
              "  ventana exterior: "
              + (f"[{run.exterior[0]:g}, {run.exterior[1]:g}] eV"
                 if run.exterior else "todas las bandas del cálculo"),
              "  ventana congelada: "
              + (f"[{run.congelada[0]:g}, {run.congelada[1]:g}] eV — "
                 f"{run.n_congeladas[0]}–{run.n_congeladas[1]} bandas por "
                 f"punto k, reproducidas exactas"
                 if run.congelada else "ninguna"),
              f"  Ω_I: {run.omega_I_hist[0]:.4f} → "
              f"{run.omega_I_hist[-1]:.4f} Å² en "
              f"{len(run.omega_I_hist) - 1} pasos",
              "  Ω_I mide cuánto se parece el subespacio elegido al de sus "
              "vecinos en k. Es lo",
              "  único que el desenredado puede bajar, y a partir de aquí ya "
              "no se mueve."]
    if run.disp is not None:
        d, d0 = run.disp, run.disp0
        L += ["", "Dispersión (Å²):",
              f"  Ω total      {d.omega:10.4f}   ({d.omega / run.nwann:.4f} "
              f"por función)",
              f"  Ω_I          {d.omega_I:10.4f}   invariante de gauge: no "
              f"baja al minimizar",
              f"  Ω_D          {d.omega_D:10.4f}",
              f"  Ω_OD         {d.omega_OD:10.4f}",
              f"  suma         {d.omega_I + d.omega_D + d.omega_OD:10.4f}   "
              f"(tiene que ser Ω)"]
        if d0 is not None and run.hist is not None and len(run.hist) > 1:
            L.append(f"  minimización: {d0.omega:.4f} → {d.omega:.4f} Å² en "
                     f"{len(run.hist) - 1} pasos; Ω_I se movió "
                     f"{getattr(d, 'deriva_I', 0.0):.1e} Å²")
    if run.disp is not None:
        L += ["", "Centros y dispersión de cada función:"]
        asign = asignar(run.disp.centros, atoms) if atoms is not None else None
        for i, (c, s) in enumerate(zip(run.disp.centros, run.disp.spreads)):
            extra = ""
            if asign:
                sitio, dd = asign[i]
                extra = f"   ← {sitio} (a {dd:.3f} Å)"
            L.append(f"  {i + 1:2d}  ({c[0]:8.4f},{c[1]:8.4f},{c[2]:8.4f}) Å"
                     f"   Ω_n = {s:7.4f} Å²{extra}")
    if run.HR is not None and run.cell is not None:
        d, a = decaimiento(run.HR, run.R, run.cell)
        L += ["", "Decaimiento de H(R):",
              f"  |R| = 0            max|H| = {a[0]:9.3e} eV",
              f"  |R| = {d[-1]:6.2f} Å      max|H| = {a[-1]:9.3e} eV"
              f"   (razón {a[-1] / max(a[0], 1e-30):.1e})"]
        if a[-1] > 0.05 * a[0]:
            run.avisos.append(
                "H(R) apenas ha decaído al borde de la superred: la base no "
                "está localizada.\n  Interpolar con esto INVENTA estructura "
                "entre los puntos de la malla. Prueba otras\n  proyecciones, "
                "o una malla más densa.")
    if np.isfinite(run.error_malla):
        ok = run.error_malla < TOL_EXACTA
        etiqueta = ("Reproduce las bandas CONGELADAS en la malla"
                    if run.omega_I_hist is not None
                    else "Reproduce la malla de partida")
        L += ["", f"{etiqueta}: "
                  f"max|ΔE| = {run.error_malla:.1e} eV"
                  + ("  ✓ (es interpolación, tiene que ser exacto)"
                     if ok else "  ← MAL: hay un error de índices")]
    if run.E_dft is not None and run.E_wann is not None:
        d = run.E_wann - run.E_dft
        L += ["", "Contra las bandas de DFT en puntos que NO estaban en la "
                  "malla:",
              f"  máximo {np.abs(d).max() * 1000:8.1f} meV      "
              f"rms {np.sqrt((d ** 2).mean()) * 1000:6.1f} meV"]
        if run.E_sin_gauge is not None:
            d0 = run.E_sin_gauge - run.E_dft
            L.append("  sin gauge (transformando las energías propias "
                     "directamente):")
            L.append(f"  máximo {np.abs(d0).max() * 1000:8.1f} meV      "
                     f"rms {np.sqrt((d0 ** 2).mean()) * 1000:6.1f} meV"
                     f"   ← el gauge es {np.abs(d0).max() / max(np.abs(d).max(), 1e-12):.1f}× mejor")
    else:
        run.avisos.append(
            "No has comparado con bandas de DFT. Que la interpolación "
            "reproduzca la malla es\n  trivial —es interpolación—; lo único "
            "que dice si el modelo sirve es compararla en\n  puntos que no "
            "estaban. Corre 4_bands.in y vuelve con --dft-bands.")
    if run.sv_min == run.sv_min and run.sv_min < 0.2:
        run.avisos.append(
            f"El valor singular más pequeño de A es {run.sv_min:.3f}: alguno "
            f"de los orbitales de\n  prueba casi no solapa con ninguna banda, "
            f"y la ortonormalización lo va a\n  amplificar hasta convertirlo "
            f"en ruido. Cambia esa proyección.")
    for a in run.avisos:
        L += ["", f"AVISO: {a}"]
    return "\n".join(L)


def export(run, outdir: str = "wannier") -> list:
    """Deja el modelo en texto: H(R), centros, bandas interpoladas y el informe."""
    out = Path(outdir); out.mkdir(parents=True, exist_ok=True)
    escritos = []

    if run.HR is not None:
        f = out / "WANNIER_hr.dat"
        nw, nr = run.HR.shape[1], len(run.R)
        L = [f"# H(R) escrito por Olla-DFT  ({run.fuente})", f"{nw:12d}",
             f"{nr:12d}"]
        for i in range(0, nr, 15):
            L.append("".join(f"{int(d):5d}" for d in run.deg[i:i + 15]))
        for ir in range(nr):
            for n in range(nw):
                for m in range(nw):
                    z = run.HR[ir, m, n]
                    L.append(f"{run.R[ir, 0]:5d}{run.R[ir, 1]:5d}"
                             f"{run.R[ir, 2]:5d}{m + 1:5d}{n + 1:5d}"
                             f"{z.real:14.6f}{z.imag:14.6f}")
        f.write_text("\n".join(L) + "\n", encoding="utf-8")
        escritos.append(str(f))

    if run.disp is not None:
        f = out / "WANNIER_centros.dat"
        L = ["#  n   x(A)        y(A)        z(A)        Omega_n(A^2)"]
        for i, (c, s) in enumerate(zip(run.disp.centros, run.disp.spreads)):
            L.append(f"{i + 1:4d} {c[0]:12.6f}{c[1]:12.6f}{c[2]:12.6f}"
                     f"{s:14.6f}")
        L.append(f"# Omega = {run.disp.omega:.6f}   Omega_I = "
                 f"{run.disp.omega_I:.6f}   Omega_D = {run.disp.omega_D:.6f}"
                 f"   Omega_OD = {run.disp.omega_OD:.6f}")
        f.write_text("\n".join(L) + "\n", encoding="utf-8")
        escritos.append(str(f))

    if run.E_wann is not None:
        f = out / "WANNIER_bandas.dat"
        x = getattr(run, "_x", np.arange(len(run.E_wann)))
        cab = "#  x(1/A)" + "".join(f"   wann{i + 1:02d}"
                                    for i in range(run.E_wann.shape[1]))
        cols = [x, *run.E_wann.T]
        if run.E_dft is not None:
            cab += "".join(f"    dft{i + 1:02d}"
                           for i in range(run.E_dft.shape[1]))
            cols += list(run.E_dft.T)
        np.savetxt(f, np.column_stack(cols), header=cab[1:], fmt="%12.6f")
        escritos.append(str(f))

    f = out / "WANNIER.txt"
    f.write_text(report(run) + "\n", encoding="utf-8")
    escritos.append(str(f))
    return escritos


def plot(run, outfile: str = "wannier", formats="pdf,png", theme: str = None,
         size: str = None, family: str = None, background: str = None,
         palette=None, usetex: bool = None, width="single",
         journal: str = "generic", aspect: float = 0.78, mono: bool = False,
         dpi: int = None) -> list:
    """Bandas interpoladas sobre las de DFT. Si no coinciden, se ve."""
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError as exc:                              # pragma: no cover
        raise RuntimeError("matplotlib no está instalado.") from exc
    if run.E_wann is None:
        raise FaltanDatos("no hay bandas interpoladas que dibujar.")

    st = qstyle.apply(theme, size=size, family=family, background=background,
                      palette=palette, usetex=usetex, mono=mono)
    fig, ax = qstyle.new_figure(width, journal, aspect)
    cols = qstyle.palette(3, mono=mono)
    x = getattr(run, "_x", np.arange(len(run.E_wann)))
    cero = run.fermi if run.fermi is not None else 0.0

    if run.E_dft is not None:
        for i in range(run.E_dft.shape[1]):
            ax.plot(x, run.E_dft[:, i] - cero, lw=st["line"] * 2.2,
                    color=cols[1], alpha=0.35, solid_capstyle="round",
                    label="DFT" if i == 0 else None, zorder=1)
    for i in range(run.E_wann.shape[1]):
        ax.plot(x, run.E_wann[:, i] - cero, lw=st["line"], color=cols[0],
                label="Wannier" if i == 0 else None, zorder=3)
    for idx, et in run.camino or []:
        if idx < len(x):
            ax.axvline(x[idx], color=st.get("grid", "0.8"), lw=0.5, zorder=0)
    if run.camino:
        ax.set_xticks([x[i] for i, _ in run.camino if i < len(x)])
        ax.set_xticklabels([e for i, e in run.camino if i < len(x)])
    ax.set_xlim(float(x[0]), float(x[-1]))
    ax.set_ylabel("E − E$_F$ (eV)" if run.fermi is not None else "E (eV)")
    if run.E_dft is not None:
        ax.legend(frameon=False, fontsize=st["legend"], loc="best")
    escritos = qstyle.save(fig, outfile, formats, dpi=dpi, modulo="wannier")
    plt.close(fig)

    if run.hist is not None and len(run.hist) > 2:
        fig2, ax2 = qstyle.new_figure(width, journal, aspect)
        ax2.plot(np.arange(len(run.hist)), run.hist, lw=st["line"],
                 color=cols[0])
        if run.disp is not None:
            ax2.axhline(run.disp.omega_I, color=cols[1], lw=st["line"],
                        dashes=[4.0, 2.0],
                        label=f"$\\Omega_I$ = {run.disp.omega_I:.3f} "
                              f"{qstyle.angstrom()}$^2$ (invariante)")
            ax2.legend(frameon=False, fontsize=st["legend"])
        ax2.set_xlabel("iteración")
        ax2.set_ylabel(f"$\\Omega$ ({qstyle.angstrom()}$^2$)")
        escritos += qstyle.save(fig2, str(outfile) + "_omega", formats,
                                dpi=dpi, modulo="wannier")
        plt.close(fig2)
    return escritos


def correr(outdir: str = "wannier", pw_cmd: str = None, nproc: int = None,
           pw2wan_cmd: str = None, con_bandas: bool = True,
           timeout: float = None, verbose: bool = True):
    """Lanza los pasos en orden: scf → nscf → pw2wannier90 → bandas.

    No se paraleliza nada porque cada paso necesita el anterior: el nscf
    lee la densidad del scf y pw2wannier90 lee las funciones de onda del
    nscf. Se para en el primero que falle, con el mensaje de QE delante,
    porque seguir con los siguientes solo produce errores derivados que
    despistan.
    """
    import shutil
    import subprocess
    from qekit.core import runner as run_mod

    out = Path(outdir)
    base = run_mod.build_command(pw_cmd, nproc)
    if pw2wan_cmd:
        p2w = base[:-1] + pw2wan_cmd.split() if len(base) > 1 else \
            pw2wan_cmd.split()
    else:
        p2w = base[:-1] + [str(Path(base[-1]).with_name("pw2wannier90.x"))]
    if not shutil.which(p2w[-1]) and not Path(p2w[-1]).exists():
        raise FaltanDatos(
            f"no encuentro pw2wannier90.x en {p2w[-1]}. Viene con Quantum "
            f"ESPRESSO (make pp).\nSi está en otro sitio, pásalo con "
            f"--pw2wan-cmd /ruta/a/pw2wannier90.x.")

    pasos = [("1_scf", base), ("2_nscf", base), ("3_pw2wan", p2w)]
    if con_bandas and (out / "4_bands.in").exists():
        pasos.append(("4_bands", base))
    hechos = []
    for nombre, cmd in pasos:
        entrada, salida = out / f"{nombre}.in", out / f"{nombre}.out"
        if not entrada.exists():
            continue
        if nombre == "4_bands" and not (out / "out_bandas").exists():
            shutil.copytree(out / "out", out / "out_bandas")
        if verbose:
            print(f"  {nombre} ...", end="", flush=True)
        with open(entrada) as fi, open(salida, "w") as fo:
            r = subprocess.run(cmd, stdin=fi, stdout=fo,
                               stderr=subprocess.STDOUT, cwd=str(out),
                               timeout=timeout)
        txt = salida.read_text(errors="replace")
        ok = r.returncode == 0 and "JOB DONE" in txt
        if verbose:
            print("  ok" if ok else "  FALLÓ")
        hechos.append((nombre, ok))
        if not ok:
            cola = "\n".join(txt.strip().split("\n")[-15:])
            raise FaltanDatos(
                f"{nombre} falló (código {r.returncode}). Últimas líneas de "
                f"{salida.name}:\n\n{cola}")
    return hechos
