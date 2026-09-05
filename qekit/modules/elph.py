# Olla-DFT — command-line toolkit for Quantum ESPRESSO
# Copyright (C) 2026 Jorge Enrique González Sevilla
# SPDX-License-Identifier: AGPL-3.0-or-later
# This program is free software: you can redistribute it and/or modify it
# under the terms of the GNU Affero General Public License as published by the Free
# Software Foundation, either version 3 of the License, or (at your option)
# any later version. See the LICENSE file for details.

"""Acoplamiento electrón-fonón: lambda, alpha^2F, Tc y un tau de verdad.

POR QUÉ IMPORTA
---------------
El módulo de transporte de Olla-DFT usa CRTA: supone que el tiempo de
relajación tau es el MISMO para todos los estados y a todas las energías.
Es la aproximación más burda que existe, y la razón por la que la
conductividad sale siempre como "sigma/tau" — un número sin unidades
útiles hasta que uno inventa un tau.

El acoplamiento electrón-fonón da un tau de verdad. Este módulo NO lo
enchufa solo en el de transporte: `olla-dft transport` sigue reportando
sigma/tau y kappa_e/tau, y eres tú quien multiplica esas columnas de
TRANSPORTE.dat por el tau(T) de la tabla de este informe (o de ELPH.txt)
para tener sigma en S/m. La secuencia (es la de la receta "termoelectrico"
de `olla-dft recipes`, resumida):

    olla-dft transport Si.cif --grid 24x24x24 --run -o trans
    olla-dft transport Si.cif --collect -o trans       # sigma/tau, S, ...
    olla-dft elph Si.cif -o elph --qgrid 2x2x2         # lambda y tau(T)
    sigma(T) = [sigma/tau](T) * tau(T)                 # a mano

De la misma DFPT que ya se corre para fonones sale, con
`electron_phonon='interpolated'`, la anchura de línea de cada modo, y de
ahí:

    lambda      = 2 * integral de alpha^2F(w)/w dw
    1/tau(T)    = 2*pi*lambda*k_B*T/hbar        (régimen de alta T)
    omega_log   = exp( (2/lambda) * integral ln(w) alpha^2F(w)/w dw )
    T_c         = (omega_log/1.2) * exp[ -1.04(1+lambda) /
                                         (lambda - mu*(1+0.62 lambda)) ]

La última es Allen-Dynes, y da la temperatura crítica de un
superconductor convencional.

LO QUE HAY QUE SABER
--------------------
1. **La malla de k tiene que ser MUY densa.** El acoplamiento se calcula
   sobre la superficie de Fermi; con la malla que basta para la energía
   total, lambda sale mal por factores. La regla práctica es una malla de
   k al menos cuatro veces más fina en cada dirección que la de q, y
   comprobar que lambda no se mueva al refinarla.

2. **El ensanchamiento no es un detalle.** ph.x calcula lambda para una
   serie de ensanchamientos (`el_ph_nsigma`). El valor bueno es el del
   PLATÓ: donde lambda deja de depender del ensanchamiento. Si no hay
   plató, la malla de k es insuficiente. Olla-DFT reporta la serie entera y
   señala el plató, en vez de dar un solo número.

3. **mu\\* es empírico.** El pseudopotencial de Coulomb vale entre 0.10 y
   0.16 y NO se calcula aquí. Tc depende mucho de él, así que Olla-DFT da
   Tc para un rango, no un solo valor.

4. **1/tau = 2 pi lambda k_B T / hbar vale a alta temperatura**, por
   encima de la temperatura de Debye. Por debajo, la dispersión se
   congela y esa fórmula sobreestima muchísimo. Olla-DFT lo dice y marca el
   régimen.
"""

import re
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

from qekit.core import provenance
from qekit.core import style as qstyle
from qekit.core.compat import trapezoid
from qekit.core.errors import ErrorDeUso

HBAR_EVS = 6.582119569e-16       # eV*s
KB_EV = 8.617333262e-5           # eV/K
THZ_K = 47.9924                  # 1 THz en K (h*nu/k_B)
THZ_CM1 = 33.35641

#: Cuanto mas fina que la del scf se pide la malla de k del nscf.
FACTOR_NSCF = 2.0


@dataclass
class ElPhRun:
    sigmas: np.ndarray = None       # ensanchamientos (Ry)
    lambdas: np.ndarray = None      # lambda por ensanchamiento
    dos_ef: np.ndarray = None       # N(EF) por ensanchamiento
    omega_log: np.ndarray = None    # K, por ensanchamiento
    Tc: np.ndarray = None           # K, por ensanchamiento
    mustar: float = 0.10
    freq: np.ndarray = None         # THz, malla de alpha^2F
    a2F: np.ndarray = None
    omega_2: float = None           # K, media cuadrática de alpha^2F (factor f2)
    i_plato: int = None
    nq: int = 0
    Tc_fuente: str = None           # de dónde sale la columna Tc(K)
    avisos: list = field(default_factory=list)

    @property
    def lam(self) -> float:
        if self.lambdas is None:
            return float("nan")
        i = self.i_plato if self.i_plato is not None else len(self.lambdas) - 1
        return float(self.lambdas[i])

    @property
    def wlog(self) -> float:
        if self.omega_log is None:
            return float("nan")
        i = self.i_plato if self.i_plato is not None else len(self.omega_log) - 1
        return float(self.omega_log[i])


# ----------------------------------------------------------------------
# Fórmulas
# ----------------------------------------------------------------------
def lambda_de_a2F(freq_thz, a2F) -> float:
    """lambda = 2 * int a2F(w)/w dw."""
    w = np.asarray(freq_thz, dtype=float)
    f = np.asarray(a2F, dtype=float)
    sel = w > 1e-9
    return float(2.0 * trapezoid(f[sel] / w[sel], w[sel]))


def omega_log_de_a2F(freq_thz, a2F, lam: float = None) -> float:
    """omega_log en K: la media logarítmica que pesa Allen-Dynes."""
    w = np.asarray(freq_thz, dtype=float)
    f = np.asarray(a2F, dtype=float)
    sel = w > 1e-9
    if lam is None:
        lam = lambda_de_a2F(w, f)
    if lam <= 0:
        return float("nan")
    integ = trapezoid(np.log(w[sel]) * f[sel] / w[sel], w[sel])
    return float(np.exp(2.0 / lam * integ) * THZ_K)


def omega_2(freq_thz, a2F, lam: float = None) -> float:
    """Media cuadrática de la frecuencia, en K. Entra en el factor f2."""
    w = np.asarray(freq_thz, dtype=float)
    f = np.asarray(a2F, dtype=float)
    sel = w > 1e-9
    if lam is None:
        lam = lambda_de_a2F(w, f)
    if lam <= 0:
        return float("nan")
    return float(np.sqrt(2.0 / lam * trapezoid(f[sel] * w[sel], w[sel]))
                 * THZ_K)


def factores_correccion(lam: float, mustar: float,
                        omega_log_K: float = None,
                        omega_2_K: float = None) -> tuple:
    """Los factores f1 y f2 de Allen-Dynes para acoplamiento fuerte.

    La fórmula desnuda SUBESTIMA Tc cuando lambda pasa de 1: para el
    plomo da 6.0 K contra 7.2 K experimentales. f1 corrige el
    acoplamiento fuerte y f2 la forma del espectro. Con lambda < 1 los
    dos valen casi 1 y no cambian nada.
    """
    if not np.isfinite(lam) or lam <= 0:
        return 1.0, 1.0
    L1 = 2.46 * (1.0 + 3.8 * mustar)
    f1 = (1.0 + (lam / L1) ** 1.5) ** (1.0 / 3.0)
    f2 = 1.0
    if omega_log_K and omega_2_K and np.isfinite(omega_2_K) \
            and omega_log_K > 0:
        r = omega_2_K / omega_log_K
        L2 = 1.82 * (1.0 + 6.3 * mustar) * r
        f2 = 1.0 + ((r - 1.0) * lam ** 2) / (lam ** 2 + L2 ** 2)
    return float(f1), float(f2)


def allen_dynes(lam: float, omega_log_K: float, mustar: float = 0.10,
                omega_2_K: float = None, correcciones: bool = True) -> float:
    """Tc de Allen-Dynes, en K.

    Devuelve 0 si el denominador se vuelve no positivo: eso quiere decir
    que el acoplamiento no basta para superconducir con ese mu*, no que
    la fórmula falle.
    """
    if not np.isfinite(lam) or not np.isfinite(omega_log_K) or lam <= 0:
        return float("nan")
    den = lam - mustar * (1.0 + 0.62 * lam)
    if den <= 0:
        return 0.0
    f1, f2 = factores_correccion(lam, mustar, omega_log_K, omega_2_K) \
        if correcciones else (1.0, 1.0)
    return float(f1 * f2 * omega_log_K / 1.2 *
                 np.exp(-1.04 * (1.0 + lam) / den))


def tau_elph(lam: float, T) -> np.ndarray:
    """Tiempo de relajación por fonones, en segundos.

        1/tau = 2*pi*lambda*k_B*T/hbar

    Vale por encima de la temperatura de Debye. Por debajo sobreestima la
    dispersión, porque los modos que dispersan se van congelando.
    """
    T = np.atleast_1d(np.asarray(T, dtype=float))
    with np.errstate(divide="ignore"):
        tau = HBAR_EVS / (2.0 * np.pi * lam * KB_EV * T)
    return tau


# ----------------------------------------------------------------------
# Lectura de la salida de ph.x y de lambda.x
# ----------------------------------------------------------------------
def leer_elph_ph(path) -> ElPhRun:
    """Lee lambda y N(EF) por ensanchamiento de la salida de ph.x."""
    p = Path(path)
    if p.is_dir():
        cand = [p / "ph.out"] + sorted(p.glob("*.out"))
        p = next((c for c in cand if c.exists()), None)
        if p is None:
            raise ErrorDeUso(f"no hay ninguna salida de ph.x en {path}.")
    texto = p.read_text(errors="ignore")
    if "electron-phonon" not in texto.lower() and "lambda" not in texto:
        raise ErrorDeUso(
            f"{p.name} no trae acoplamiento electrón-fonón. Hace falta correr "
            "ph.x con\n  electron_phonon = 'interpolated'\ny una malla de k "
            "mucho más densa que la del scf normal.")

    # "Gaussian Broadening: 0.005 Ry, ngauss= 1" ... "DOS = ... states/spin/Ry"
    sig = [float(x) for x in re.findall(
        r"Gaussian Broadening:\s*([\d.]+)\s*Ry", texto)]
    dos = [float(x) for x in re.findall(
        r"DOS =\s*([\d.]+)\s*states/spin/Ry", texto)]
    run = ElPhRun()
    if sig:
        n = len(dict.fromkeys(sig))
        run.sigmas = np.array(sorted(dict.fromkeys(sig)))
        run.dos_ef = np.array(dos[:n]) if len(dos) >= n else None
    return run


def leer_a2F(path) -> tuple:
    """Lee alpha^2F(w) de un archivo a2F.dos* que escribe lambda.x."""
    p = Path(path)
    datos = np.loadtxt(p, comments="#")
    if datos.ndim != 2 or datos.shape[1] < 2:
        raise ErrorDeUso(f"{p.name} no tiene la forma de un alpha^2F.")
    return datos[:, 0], datos[:, 1]


def build_lambda_input(qpuntos, pesos, archivos, emax_thz: float = 20.0,
                       degaussq: float = 0.12, ngaussq: int = 1,
                       mustar: float = 0.10) -> str:
    """Arma el input de lambda.x.

    Escribirlo a mano es el paso que hace que la gente abandone este
    cálculo: hay que copiar los q-puntos y sus pesos de la salida de ph.x
    en el mismo orden que los archivos elph, y un desorden ahí da un
    lambda equivocado sin ningún error.
    """
    if not (len(qpuntos) == len(pesos) == len(archivos)):
        raise ErrorDeUso(
            f"no cuadran los tamaños: {len(qpuntos)} q-puntos, {len(pesos)} "
            f"pesos y {len(archivos)} archivos. Los tres tienen que ir en el "
            "MISMO orden.")
    lineas = [f"{emax_thz}  {degaussq}  {ngaussq}", f"{len(qpuntos)}"]
    for q, w in zip(qpuntos, pesos):
        lineas.append(f"  {q[0]:12.7f} {q[1]:12.7f} {q[2]:12.7f} {w:10.4f}")
    lineas += list(archivos)
    lineas.append(f"{mustar}")
    return "\n".join(lineas) + "\n"


def _mustar_de_lambda_in(p: Path) -> float:
    """mu* que se le dio a lambda.x: la última línea numérica de lambda.in.

    Si no hay lambda.in al lado se devuelve None y se usa el mu* por
    omisión de ElPhRun.
    """
    for nombre in ("lambda.in", p.stem + ".in"):
        cand = p.with_name(nombre)
        if not cand.exists():
            continue
        lineas = [l.strip() for l in cand.read_text(errors="ignore")
                  .splitlines() if l.strip()]
        if lineas:
            try:
                mu = float(lineas[-1])
            except ValueError:
                return None
            if 0.0 <= mu < 1.0:
                return mu
    return None


def _tabla_tc_lambda_out(texto: str):
    """La tabla final de lambda.x: 'lambda  omega_log  T_c', una fila por
    ensanchamiento. Devuelve (lambda, omega_log, Tc) o None si no está."""
    m = re.search(r"lambda\s+omega_log\s+T_c\s*\n((?:[ \t]*\S+[ \t]+\S+"
                  r"[ \t]+\S+[ \t]*\n?)+)", texto)
    if not m:
        return None
    filas = []
    for linea in m.group(1).splitlines():
        partes = linea.split()
        if len(partes) != 3:
            continue
        try:
            filas.append([float(x) for x in partes])
        except ValueError:
            filas.append([float("nan")] * 3)
    if not filas:
        return None
    d = np.array(filas)
    return d[:, 0], d[:, 1], d[:, 2]


def leer_lambda_out(path) -> ElPhRun:
    """Lee la salida de lambda.x: lambda, omega_log y Tc por ensanchamiento.

    Prefiere el `lambda.dat` que lambda.x escribe al lado, porque es una
    tabla con encabezado y no depende del formato del texto. En mallas de
    q gruesas, lambda.x deja <log w> y Tc en NaN: eso no es un fallo, es
    que su integral de alpha^2F sale mal muestreada. Olla-DFT lo detecta y
    recalcula omega_log por su cuenta si tiene el alpha^2F.

    Tc: lambda.x la imprime en la tabla final de lambda.out ('lambda
    omega_log T_c', Allen-Dynes sin factores de corrección y con el mu*
    del input). Si esa tabla no está, se calcula aquí fila a fila con
    `allen_dynes(lambda_i, omega_log_i, mu*)`, y `run.Tc_fuente` dice cuál
    de las dos cosas se hizo para que el informe lo cuente.
    """
    p = Path(path)
    run = ElPhRun()
    mu = _mustar_de_lambda_in(p)
    if mu is not None:
        run.mustar = mu
    dat = p.with_name("lambda.dat")
    if dat.exists():
        d = np.genfromtxt(dat, comments="#")
        if d.ndim == 1:
            d = d.reshape(1, -1)
        if d.shape[1] >= 5:
            run.sigmas, run.lambdas = d[:, 0], d[:, 1]
            run.omega_log, run.dos_ef = d[:, 3], d[:, 4]
    if run.lambdas is None:
        texto = p.read_text(errors="ignore")
        filas = re.findall(
            r"lambda\s*=\s*([\d.]+).*?<log w>=\s*(\S+)\s*K\s+"
            r"N\(Ef\)=\s*([\d.]+)\s+at degauss=\s*([\d.]+)", texto)
        if filas:
            def _f(x):
                try:
                    return float(x)
                except ValueError:
                    return float("nan")
            d = np.array([[_f(v) for v in f] for f in filas])
            run.lambdas, run.omega_log = d[:, 0], d[:, 1]
            run.dos_ef, run.sigmas = d[:, 2], d[:, 3]

    # Tc por ensanchamiento: de la tabla final de lambda.out si la hay y
    # cuadra en tamaño; si no, Allen-Dynes aquí mismo con el mu* de lambda.in
    if run.lambdas is not None:
        tabla = None
        if p.exists():
            tabla = _tabla_tc_lambda_out(p.read_text(errors="ignore"))
        if tabla is not None and len(tabla[2]) == len(run.lambdas):
            run.Tc = tabla[2]
            run.Tc_fuente = f"lambda.x (Allen-Dynes con mu* = {run.mustar:.2f})"
            if run.omega_log is None or not np.any(np.isfinite(run.omega_log)):
                # la tabla trae omega_log en K: si el .dat no lo tenía, vale
                if np.any(np.isfinite(tabla[1])):
                    run.omega_log = tabla[1]
        else:
            run.Tc = np.array([
                allen_dynes(float(l), float(w), run.mustar, correcciones=False)
                for l, w in zip(run.lambdas,
                                run.omega_log if run.omega_log is not None
                                else np.full(len(run.lambdas), np.nan))])
            run.Tc_fuente = (f"Allen-Dynes calculada por Olla-DFT con mu* = "
                             f"{run.mustar:.2f} (lambda.out no trae la tabla "
                             "de T_c)")
    if run.omega_log is not None and not np.any(np.isfinite(run.omega_log)):
        run.avisos.append(
            "lambda.x dejo omega_log en NaN. Pasa cuando la malla de q "
            "es gruesa: su integral de alpha^2F sale mal muestreada. "
            "lambda sigue siendo utilizable, pero para Tc hace falta "
            "omega_log: calculalo desde el alpha^2F si tienes el archivo "
            "a2F.dos, o refina la malla de q.")
    return run


def plato(lambdas, tol: float = 0.05) -> int:
    """Índice del PLATÓ: donde lambda deja de depender del ensanchamiento.

    Sin plató, el número que se reporte es arbitrario. Se busca el tramo
    más largo de valores consecutivos que no varíen más de `tol` en
    términos relativos, y se devuelve su punto medio.
    """
    x = np.asarray(lambdas, dtype=float)
    if x.size < 3:
        return None
    mejor, mejor_len = None, 0
    i = 0
    while i < len(x):
        j = i
        while j + 1 < len(x) and abs(x[j + 1] - x[i]) <= tol * abs(x[i]):
            j += 1
        if j - i + 1 > mejor_len:
            mejor_len, mejor = j - i + 1, (i + j) // 2
        i = j + 1
    return mejor if mejor_len >= 3 else None


# ----------------------------------------------------------------------
# Reporte
# ----------------------------------------------------------------------
def report(run: ElPhRun, T_debye: float = None,
           mus=(0.10, 0.13, 0.16)) -> str:
    lines = ["--- Acoplamiento electrón-fonón ---"]
    if run.nq:
        lines.append(f"Puntos q: {run.nq}")

    if run.lambdas is not None:
        lines += ["", "Dependencia del ENSANCHAMIENTO (el valor bueno es el "
                  "del plató):",
                  f"  {'sigma(Ry)':>10s} {'lambda':>8s} {'N(EF)':>10s} "
                  f"{'w_log(K)':>10s} {'Tc(K)':>8s}"]
        for i in range(len(run.lambdas)):
            s = run.sigmas[i] if run.sigmas is not None and \
                i < len(run.sigmas) else float("nan")
            d = run.dos_ef[i] if run.dos_ef is not None else float("nan")
            w = run.omega_log[i] if run.omega_log is not None else float("nan")
            t = run.Tc[i] if run.Tc is not None else float("nan")
            marca = "  <- plató" if i == run.i_plato else ""
            lines.append(f"  {s:10.4f} {run.lambdas[i]:8.4f} {d:10.4f} "
                         f"{w:10.2f} {t:8.3f}{marca}")
        if run.Tc_fuente:
            lines.append(f"  Tc(K) de la tabla: {run.Tc_fuente}.")
        if run.i_plato is None:
            lines += ["",
                      "NO hay plató: lambda cambia con el ensanchamiento a lo "
                      "largo de toda la\nserie. Eso quiere decir que la malla "
                      "de k es insuficiente — el acoplamiento\nse calcula "
                      "sobre la superficie de Fermi y necesita una malla mucho "
                      "más\nfina que la de la energía total. Cualquier lambda "
                      "que se reporte de aquí es\narbitrario."]

    lam, wlog = run.lam, run.wlog
    if np.isfinite(lam):
        lines += ["", f"lambda   = {lam:.4f}",
                  f"omega_log = {wlog:.1f} K"]
        regimen = ("acoplamiento débil" if lam < 0.5 else
                   "acoplamiento intermedio" if lam < 1.0 else
                   "acoplamiento fuerte")
        lines.append(f"  ({regimen})")

        w2 = run.omega_2 if run.omega_2 and np.isfinite(run.omega_2) else None
        if w2:
            lines.append(f"omega_2  = {w2:.1f} K  (entra en el factor de forma f2)")
        lines += ["", "Temperatura crítica (Allen-Dynes), según mu*:"]
        for mu in mus:
            lines.append(f"  mu* = {mu:.2f}   Tc = "
                         f"{allen_dynes(lam, wlog, mu, omega_2_K=w2):7.3f} K")
        lines.append("  mu* es empírico (0.10-0.16) y NO se calcula aquí; por "
                     "eso se da el rango.")
        if lam < 0.3:
            lines.append("  Con lambda tan bajo, Tc es prácticamente cero: "
                         "este material no es\n  un superconductor "
                         "convencional apreciable.")

        lines += ["", "Tiempo de relajación por fonones (1/tau = "
                  "2*pi*lambda*k_B*T/hbar):",
                  f"  {'T (K)':>8s} {'tau (fs)':>12s}"]
        for T in (100.0, 300.0, 500.0, 800.0):
            tau = float(tau_elph(lam, T)[0])
            lines.append(f"  {T:8.0f} {tau * 1e15:12.3f}")
        lines.append("  Este tau es el que le falta a la CRTA de 'olla-dft "
                     "transport'. No se aplica\n  solo: multiplica las "
                     "columnas sigma/tau y kappa_e/tau de TRANSPORTE.dat "
                     "por\n  el tau(T) de esta tabla y sigma sale en S/m.")
        if T_debye:
            lines.append(f"  La fórmula vale por encima de la temperatura de "
                         f"Debye ({T_debye:.0f} K).\n  Por debajo "
                         "sobreestima la dispersión: los modos que dispersan "
                         "se congelan.")
        else:
            lines.append("  Vale por ENCIMA de la temperatura de Debye; por "
                         "debajo sobreestima la\n  dispersión. Calcúlala con "
                         "'olla-dft derived' y compara.")

    if run.a2F is not None:
        lines += ["", "alpha^2F(w):",
                  f"  máximo en {run.freq[int(np.argmax(run.a2F))]:.2f} THz "
                  f"({run.freq[int(np.argmax(run.a2F))] * THZ_CM1:.0f} cm⁻¹)"]

    for a in run.avisos:
        lines += ["", a]
    return "\n".join(lines)


def export(run: ElPhRun, outdir: str = ".") -> list:
    out = Path(outdir); out.mkdir(parents=True, exist_ok=True)
    escritos = []
    cab = provenance.header_plain(
        "acoplamiento electron-fonon",
        {"lambda": None if not np.isfinite(run.lam) else round(run.lam, 5),
         "omega_log_K": None if not np.isfinite(run.wlog)
         else round(run.wlog, 2),
         "omega_2_K": None if not (run.omega_2 and np.isfinite(run.omega_2))
         else round(run.omega_2, 2), "nq": run.nq or None},
        titulo="Electron-fonon")
    if run.lambdas is not None:
        f = out / "ELPH.dat"
        cols, nombres = [], []
        for arr, nombre in ((run.sigmas, "sigma(Ry)"),
                            (run.lambdas, "lambda"),
                            (run.dos_ef, "N(EF)"),
                            (run.omega_log, "w_log(K)"),
                            (run.Tc, "Tc(K)")):
            if arr is not None and len(arr) == len(run.lambdas):
                cols.append(arr); nombres.append(nombre)
        np.savetxt(f, np.column_stack(cols), fmt="%14.6f",
                   header=cab + "\n" + "  ".join(f"{n:>14s}" for n in nombres),
                   comments="# ")
        escritos.append(str(f))
    if run.a2F is not None:
        f = out / "A2F.dat"
        np.savetxt(f, np.column_stack([run.freq, run.a2F]), fmt="%14.6f",
                   header=cab + "\n   frecuencia(THz)      alpha^2F",
                   comments="# ")
        escritos.append(str(f))
    txt = out / "ELPH.txt"
    txt.write_text(report(run) + "\n")
    escritos.append(str(txt))
    return escritos


def plot(run: ElPhRun, outfile: str = "elph", formats="pdf,png",
         theme: str = None, family: str = None, background: str = None,
         palette=None, usetex: bool = None, width="single",
         journal: str = "generic", mono: bool = False,
         dpi: int = None) -> list:
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError as exc:                          # pragma: no cover
        raise RuntimeError("matplotlib no está instalado.") from exc

    qstyle.apply(theme, family=family, background=background,
                 palette=palette, usetex=usetex, mono=mono)
    tiene_a2F = run.a2F is not None
    n = 2 if tiene_a2F and run.lambdas is not None else 1
    size = qstyle.figure_size(width, journal, 0.75 if n == 1 else 0.42)
    fig, axes = plt.subplots(1, n, figsize=size, layout="constrained")
    axes = np.atleast_1d(axes)
    colores = qstyle.palette(4, mono=mono)

    i = 0
    if tiene_a2F:
        ax = axes[i]; i += 1
        ax.fill_between(run.freq, 0, run.a2F, alpha=0.28, color=colores[0],
                        lw=0)
        ax.plot(run.freq, run.a2F, lw=1.3, color=colores[0])
        ax.set_xlabel(r"$\omega$ (THz)")
        ax.set_ylabel(r"$\alpha^2 F(\omega)$")
        ax.set_xlim(0, run.freq[-1]); ax.set_ylim(bottom=0)
    if run.lambdas is not None and i < n:
        ax = axes[i]
        x = run.sigmas if run.sigmas is not None else \
            np.arange(len(run.lambdas))
        ax.plot(x, run.lambdas, "o-", ms=4, lw=1.2, color=colores[1])
        if run.i_plato is not None:
            ax.plot(x[run.i_plato], run.lambdas[run.i_plato], "o", ms=8,
                    mfc="none", mew=1.5, color=colores[2], label="plató")
            ax.legend(frameon=False, fontsize="small")
        ax.set_xlabel("ensanchamiento (Ry)")
        ax.set_ylabel(r"$\lambda$")
        ax.set_ylim(bottom=0)
    return qstyle.save(fig, outfile, formats, dpi=dpi)


# ----------------------------------------------------------------------
# Preparación de los TRES pasos
# ----------------------------------------------------------------------
def prepare(atoms, outdir: str = "elph", qgrid=(2, 2, 2), kgrid_scf=None,
            kgrid_nscf=None, pseudo_dir: str = None, ecutwfc: float = None,
            ecutrho: float = None, degauss: float = 0.02,
            smearing: str = "methfessel-paxton", nsigma: int = 10,
            sigma_paso: float = 0.005, tr2: float = 1e-12) -> tuple:
    """Escribe los TRES inputs del cálculo electrón-fonón.

    Son tres, no dos, y ese es el detalle que hace que la gente abandone:

      1. scf con malla de k normal
      2. **nscf con malla de k MUY densa y `la2F = .true.`** -> deja el
         archivo a2Fsave que ph.x necesita para interpolar
      3. ph.x con `electron_phonon = 'interpolated'`

    Sin el paso 2, ph.x llega hasta "electron-phonon interaction ..." y
    se muere sin decir por qué. Olla-DFT escribe los tres y los numera.
    """
    from qekit.modules import inputgen, sweep

    common = sweep.prepare_common(atoms, pseudo_dir, ecutwfc, ecutrho,
                                  insulator=False)
    out = Path(outdir); out.mkdir(parents=True, exist_ok=True)
    k_scf = tuple(kgrid_scf or sweep.default_grid(atoms, None))
    # la malla del nscf tiene que ser un MÚLTIPLO entero de la de q, y
    # mucho más fina: el acoplamiento vive en la superficie de Fermi
    # La malla del nscf tiene que ser un MULTIPLO entero de la de q (si no,
    # los q pedidos no caen sobre puntos de la malla y la interpolacion no
    # existe) y bastante mas fina que la del scf.
    k_nscf = tuple(kgrid_nscf or tuple(
        int(q) * max(2, int(np.ceil(FACTOR_NSCF * k / max(q, 1))))
        for k, q in zip(k_scf, qgrid)))

    for nombre, calc, malla, extra_sys in (
            ("1_scf", "scf", k_scf, {}),
            ("2_nscf", "nscf", k_nscf, {"la2F": True})):
        txt = inputgen.build_pw_input(
            atoms=atoms, pseudos=common["pseudos"], calculation=calc,
            prefix=common["prefix"], pseudo_dir=common["pseudo_dir"],
            ecutwfc=common["ecutwfc"], ecutrho=common["ecutrho"],
            kcard=f"K_POINTS automatic\n  {malla[0]} {malla[1]} {malla[2]} "
                  "0 0 0\n",
            insulator=False, degauss=degauss, smearing=smearing,
            conv_thr=1e-10)
        if extra_sys.get("la2F"):
            txt = re.sub(r"(&SYSTEM\n)", r"\1  la2F            = .true.\n",
                         txt, count=1)
        sweep.write_input(out / f"{nombre}.in", txt)

    ph = ["Acoplamiento electron-fonon", " &INPUTPH",
          f"   prefix = '{common['prefix']}',", "   outdir = './out',",
          f"   tr2_ph = {tr2:.1e},".replace("e-", "d-"),
          "   fildyn = 'dyn',", "   ldisp = .true.,",
          f"   nq1 = {qgrid[0]}, nq2 = {qgrid[1]}, nq3 = {qgrid[2]},",
          "   electron_phonon = 'interpolated',",
          f"   el_ph_sigma = {sigma_paso}, el_ph_nsigma = {nsigma},",
          "   fildvscf = 'dvscf',", " /"]
    sweep.write_input(out / "3_ph.in", "\n".join(ph) + "\n")

    nk = int(np.prod(k_nscf))
    rep = ["--- Acoplamiento electron-fonon ---",
           f"Estructura: {atoms.get_chemical_formula()} "
           f"({len(atoms)} átomos)",
           f"Malla de q: {qgrid[0]}x{qgrid[1]}x{qgrid[2]}",
           f"Malla de k del scf:  {k_scf[0]}x{k_scf[1]}x{k_scf[2]}",
           f"Malla de k del nscf: {k_nscf[0]}x{k_nscf[1]}x{k_nscf[2]}  "
           f"({nk} puntos antes de simetría)",
           f"Ensanchamientos: {nsigma} de {sigma_paso} Ry",
           "",
           f"Archivos en '{out.resolve()}':",
           "  1_scf.in    scf normal",
           "  2_nscf.in   malla densa CON la2F = .true.",
           "  3_ph.in     ph.x con electron_phonon = 'interpolated'",
           "",
           "Orden, y los tres son obligatorios:",
           "  pw.x -in 1_scf.in  &&  pw.x -in 2_nscf.in  &&  ph.x -in 3_ph.in",
           "",
           "El paso 2 es el que se olvida. Sin el archivo a2Fsave que "
           "deja ese nscf, ph.x llega hasta 'electron-phonon "
           "interaction ...' y se muere sin decir por que.",
           "",
           "Esto es CARO: el acoplamiento vive en la superficie de Fermi "
           "y necesita una malla de k mucho mas fina que cualquier otra "
           "cosa. Si lambda cambia al refinarla, todavia no esta "
           "convergido."]
    warn = sweep.missing_pseudo_warning(common)
    if warn:
        rep.append(warn)
    return common, "\n".join(rep)
