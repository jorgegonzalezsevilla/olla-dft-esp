# Olla-DFT — command-line toolkit for Quantum ESPRESSO
# Copyright (C) 2026 Jorge Enrique González Sevilla
# SPDX-License-Identifier: GPL-3.0-or-later
# This program is free software: you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by the Free
# Software Foundation, either version 3 of the License, or (at your option)
# any later version. See the LICENSE file for details.

"""Aproximación cuasi-armónica: expansión térmica y a(T).

La aproximación armónica no da expansión térmica: si las frecuencias no
dependen del volumen, el mínimo de la energía libre está siempre en el
mismo sitio y el sólido no se dilata nunca. La QHA arregla eso de la
manera mínima: mantiene los modos armónicos pero deja que sus frecuencias
dependan del VOLUMEN.

    F(V,T) = E_estatica(V) + F_vib(V,T)
    F_vib  = ZPE + k_B T INT g(w) ln(1 - exp(-hbar w / k_B T)) dw

Para cada T se minimiza F(V,T) en V. De V(T) salen:

- el coeficiente de expansión térmica  alpha = (1/V)(dV/dT)
- el parámetro de red a(T), comparable con dilatometría o DRX a temperatura
- el parámetro de Grüneisen  gamma = -dln(w)/dln(V)
- C_p = C_v + alpha^2 * B * V * T

QUÉ HACE FALTA Y QUÉ CUESTA
---------------------------
Fonones en VARIOS volúmenes (5 suele bastar) más la curva E(V) de la EOS.
Con DFPT eso son cinco cálculos de fonones, que es caro; por eso Olla-DFT
admite calcular las frecuencias con un potencial aprendido para explorar
primero, y dejar el DFT para el resultado final.

LÍMITES HONESTOS
----------------
La QHA no es anarmónica: cada modo sigue siendo armónico, solo cambia con
el volumen. Eso funciona bien hasta aproximadamente la mitad de la
temperatura de fusión, y falla cerca de ella. Tampoco describe
transiciones de fase ni modos blandos.

Y un caso que sí reproduce y conviene mirar: el silicio tiene expansión
térmica NEGATIVA por debajo de unos 120 K, porque las ramas transversales
acústicas tienen Grüneisen negativo. Si una implementación de QHA no
recupera ese signo, está mal.
"""

from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

from qekit.core import provenance
from qekit.core import style as qstyle

KB_EV = 8.617333262e-5
CM1_EV = 1.239841984e-4
EV_A3_GPA = 160.21766208


@dataclass
class QHAResult:
    volumenes: np.ndarray = None          # A^3
    energias: np.ndarray = None           # eV (estatica)
    frecuencias: list = field(default_factory=list)   # cm^-1 por volumen
    T: np.ndarray = None                  # K
    V_T: np.ndarray = None                # A^3
    a_T: np.ndarray = None                # A (si es cubico)
    # True si a_T es el parámetro de red CONVENCIONAL (hizo falta la
    # estructura para saber cuántas celdas primitivas caben en la
    # convencional); False si es solo la raíz cúbica del volumen primitivo
    a_convencional: bool = False
    alpha: np.ndarray = None              # 1/K, volumetrico
    B_T: np.ndarray = None                # GPa
    Cv: np.ndarray = None                 # meV/K por celda
    Cp: np.ndarray = None
    gruneisen: float = None
    natoms: int = 1
    cubico: bool = False
    avisos: list = field(default_factory=list)


def f_vib(frecuencias_cm1, T: float) -> float:
    """Energía libre vibracional armónica de un conjunto de modos, en eV.

    Los modos con w <= 0 (acústicos en Gamma, o imaginarios) se descartan;
    si hay imaginarias de verdad se avisa aparte, porque entonces la QHA no
    aplica.
    """
    w = np.asarray(frecuencias_cm1, dtype=float)
    w = w[w > 1.0]
    if w.size == 0:
        return 0.0
    e = w * CM1_EV
    zpe = 0.5 * float(np.sum(e))
    if T <= 0:
        return zpe
    x = e / (KB_EV * T)
    x = np.clip(x, 1e-12, 300.0)
    return float(zpe + KB_EV * T * np.sum(np.log1p(-np.exp(-x))))


def cv_modos(frecuencias_cm1, T: float) -> float:
    """C_v armónica de un conjunto de modos, en meV/K."""
    w = np.asarray(frecuencias_cm1, dtype=float)
    w = w[w > 1.0]
    if w.size == 0 or T <= 0:
        return 0.0
    x = (w * CM1_EV) / (KB_EV * T)
    x = np.clip(x, 1e-9, 300.0)
    occ = 1.0 / np.expm1(x)
    return float(KB_EV * np.sum(x ** 2 * np.exp(x) * occ ** 2) * 1000.0)


def factor_convencional(atoms) -> float:
    """Cuántas celdas primitivas caben en la convencional: V_conv / V_prim.

    Se cuenta por átomos (spglib estandariza las dos celdas), que es exacto:
    4 en fcc y diamante, 2 en bcc, 1 en cúbica simple. Es lo que hace falta
    para convertir el volumen por celda primitiva de la QHA en el parámetro
    de red convencional que se compara con dilatometría o DRX.
    """
    from qekit.core import structure as struct_mod
    n_prim = len(struct_mod.primitive(atoms))
    n_conv = len(struct_mod.conventional(atoms))
    return n_conv / n_prim


def es_cubico(atoms) -> bool:
    """¿Es cúbica la estructura (grupos espaciales 195–230)?"""
    from qekit.core import structure as struct_mod
    ds = struct_mod.symmetry_dataset(atoms)
    numero = ds.number if hasattr(ds, "number") else ds["number"]
    return int(numero) >= 195


def run(volumenes, energias, frecuencias, T=None, natoms: int = 1,
        cubico: bool = False, celdas_por_modo: int = 1,
        factor_conv: float = None) -> QHAResult:
    """Minimiza F(V,T) para cada T y deriva las propiedades térmicas.

    `frecuencias` es una lista con las frecuencias (cm^-1) de cada volumen.
    `celdas_por_modo` es cuántas celdas primitivas contiene la supercelda
    con que se calcularon esos modos: la energía libre se divide entre ese
    número para dejarlo todo POR CELDA, que es como está E(V).

    `factor_conv` es V_conv/V_prim (lo da `factor_convencional(atoms)`).
    Con él, y si `cubico`, a(T) es el parámetro de red CONVENCIONAL. Sin
    él, a(T) es solo la raíz cúbica del volumen primitivo, que en fcc/bcc/
    diamante NO es el parámetro de red: el resultado lo dice y se avisa.
    """
    V = np.asarray(volumenes, dtype=float)
    E = np.asarray(energias, dtype=float)
    if T is None:
        T = np.arange(0.0, 1001.0, 10.0)
    T = np.asarray(T, dtype=float)
    res = QHAResult(volumenes=V, energias=E, frecuencias=list(frecuencias),
                    T=T, natoms=natoms, cubico=cubico)

    if len(V) < 4:
        res.avisos.append(
            f"solo {len(V)} volúmenes: hacen falta al menos 4 para ajustar "
            "F(V,T) con sentido.\nCon menos, el mínimo sale de una "
            "extrapolación.")
    for i, f in enumerate(frecuencias):
        w = np.asarray(f, dtype=float)
        if np.any(w < -5.0):
            res.avisos.append(
                f"el volumen {V[i]:.2f} Å³ tiene frecuencias imaginarias: "
                "ahí la estructura no\nestá en un mínimo y la QHA no aplica. "
                "Quita ese punto o relaja mejor.")

    V_T, B_T, Cv = [], [], []
    for t in T:
        F = np.array([E[i] + f_vib(frecuencias[i], t) / celdas_por_modo
                      for i in range(len(V))])
        # parábola local alrededor del mínimo muestreado
        j = int(np.argmin(F))
        lo, hi = max(0, j - 2), min(len(V), j + 3)
        if hi - lo < 3:
            lo, hi = 0, len(V)
        coef = np.polyfit(V[lo:hi], F[lo:hi], 2)
        if coef[0] <= 0:
            V_T.append(V[j]); B_T.append(np.nan)
        else:
            v0 = -coef[1] / (2.0 * coef[0])
            V_T.append(float(np.clip(v0, V.min(), V.max())))
            B_T.append(float(2.0 * coef[0] * v0 * EV_A3_GPA))
        cv = np.interp(V_T[-1], V,
                       [cv_modos(f, t) / celdas_por_modo
                        for f in frecuencias])
        Cv.append(float(cv))

    res.V_T = np.array(V_T)
    res.B_T = np.array(B_T)
    res.Cv = np.array(Cv)
    # alpha = (1/V) dV/dT. Con una sola temperatura no hay derivada que
    # tomar: se devuelve NaN y se avisa, en vez de reventar.
    if len(T) < 2:
        res.alpha = np.full_like(res.V_T, np.nan)
        res.avisos.append(
            "una sola temperatura: la expansion termica es una DERIVADA "
            "respecto de T y no\nse puede calcular con un punto. Pasa un "
            "rango de temperaturas.")
    else:
        res.alpha = np.gradient(res.V_T, T) / res.V_T
    # C_p = C_v + alpha^2 B V T   (unidades: meV/K)
    with np.errstate(invalid="ignore"):
        extra = (res.alpha ** 2 * res.B_T / EV_A3_GPA * res.V_T * T) * 1000.0
    res.Cp = res.Cv + np.nan_to_num(extra)
    if cubico:
        if factor_conv is not None and factor_conv > 0:
            res.a_T = (res.V_T * float(factor_conv)) ** (1.0 / 3.0)
            res.a_convencional = True
        else:
            res.a_T = res.V_T ** (1.0 / 3.0)
            res.a_convencional = False
            res.avisos.append(
                "a(T) es la raíz cúbica del volumen de la celda PRIMITIVA. "
                "En fcc, bcc o\ndiamante eso NO es el parámetro de red "
                "convencional (difieren en 4^(1/3) o\n2^(1/3)). Pasa la "
                "estructura con --structure para que se convierta.")

    # Gruneisen medio: -dln<w>/dlnV
    medias = []
    for f in frecuencias:
        w = np.asarray(f, dtype=float)
        w = w[w > 1.0]
        medias.append(float(np.mean(w)) if w.size else np.nan)
    medias = np.array(medias)
    ok = np.isfinite(medias) & (medias > 0)
    if ok.sum() >= 3:
        p = np.polyfit(np.log(V[ok]), np.log(medias[ok]), 1)
        res.gruneisen = float(-p[0])
    return res


def report(res: QHAResult, T_ref: float = 300.0) -> str:
    i = int(np.argmin(np.abs(res.T - T_ref)))
    lines = ["--- Aproximación cuasi-armónica ---",
             f"{len(res.volumenes)} volúmenes entre "
             f"{res.volumenes.min():.2f} y {res.volumenes.max():.2f} Å³"]
    if res.gruneisen is not None:
        lines.append(f"Parámetro de Grüneisen medio: {res.gruneisen:.3f}")
    lines += ["",
              f"A {res.T[i]:.0f} K:",
              f"  volumen de equilibrio: {res.V_T[i]:.3f} Å³"]
    if np.isfinite(res.alpha[i]):
        lines += [f"  expansión térmica volumétrica: "
                  f"{res.alpha[i]*1e6:.2f} × 10⁻⁶ K⁻¹",
                  f"  lineal (α/3): {res.alpha[i]/3*1e6:.2f} × 10⁻⁶ K⁻¹"]
    lines += [
              f"  C_v = {res.Cv[i]:.4f} meV/K   C_p = {res.Cp[i]:.4f} meV/K",
              f"  B(T) = {res.B_T[i]:.1f} GPa"]
    if res.a_T is not None:
        if res.a_convencional:
            lines.append(f"  parámetro de red (celda convencional): "
                         f"{res.a_T[i]:.4f} Å")
        else:
            lines.append(f"  V_prim^(1/3) (NO es el parámetro de red "
                         f"convencional): {res.a_T[i]:.4f} Å")

    finita = np.isfinite(res.alpha)
    neg = res.T[finita & (res.alpha < 0) & (res.T > 1)]
    if neg.size:
        lines += ["",
                  f"Expansión térmica NEGATIVA por debajo de "
                  f"{neg.max():.0f} K.",
                  "No es un error: en el silicio y otros con estructura "
                  "tipo diamante ocurre de\nverdad, y viene de que las "
                  "ramas transversales acústicas tienen Grüneisen\n"
                  "negativo. Si la implementación no lo reprodujera, sería "
                  "señal de que algo\nfalla."]
    for a in res.avisos:
        lines.append(f"\nAVISO: {a}")
    lines += ["",
              "La QHA deja que las frecuencias dependan del volumen, pero "
              "cada modo sigue\nsiendo armónico. Vale hasta ~la mitad de la "
              "temperatura de fusión; cerca de\nella hace falta "
              "anarmonicidad explícita."]
    return "\n".join(lines)


def export(res: QHAResult, outdir: str = ".") -> list:
    out = Path(outdir); out.mkdir(parents=True, exist_ok=True)
    f = out / "QHA.dat"
    cab = provenance.header_plain(
        "cuasi-armónica",
        {"n_volumenes": len(res.volumenes),
         "gruneisen": round(res.gruneisen, 4) if res.gruneisen else None},
        titulo="Aproximacion cuasi-armonica")
    cols = [res.T, res.V_T, res.alpha, res.Cv, res.Cp, res.B_T]
    nombres = f"{'T(K)':>12s} {'V(A^3)':>14s} {'alpha(1/K)':>14s} " \
              f"{'Cv(meV/K)':>14s} {'Cp(meV/K)':>14s} {'B(GPa)':>12s}"
    if res.a_T is not None:
        cols.append(res.a_T)
        # el nombre de la columna dice qué es: parámetro de red convencional
        # o solo la raíz cúbica del volumen primitivo
        nombres += f" {'a_conv(A)' if res.a_convencional else 'Vprim^1/3(A)':>12s}"
    np.savetxt(f, np.column_stack(cols), fmt="%14.6e",
               header=cab + "\n" + nombres, comments="# ")
    return [str(f)]


def plot(res: QHAResult, outfile: str = "qha", formats="pdf,png",
         theme: str = None, family: str = None, background: str = None,
         palette=None, usetex: bool = None, width="double",
         journal: str = "generic", aspect: float = 0.38,
         mono: bool = False, dpi: int = None) -> list:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    st = qstyle.apply(theme, family=family, background=background,
                      palette=palette, usetex=usetex, mono=mono)
    fig = plt.figure(figsize=qstyle.figure_size(width, journal, aspect),
                     layout="constrained")
    c = qstyle.palette(3, mono=mono)
    ax = [qstyle.finish_axes(fig.add_subplot(1, 3, i + 1)) for i in range(3)]

    ax[0].plot(res.T, res.V_T, color=c[0], lw=st["line"])
    ax[0].set_xlabel("T (K)"); ax[0].set_ylabel(r"$V$ (Å$^3$)")
    qstyle.panel_label(ax[0], "(a)")

    ax[1].plot(res.T, res.alpha * 1e6, color=c[1], lw=st["line"])
    ax[1].axhline(0.0, color=qstyle.INK_FAINT, lw=st["axis_line"])
    ax[1].set_xlabel("T (K)")
    ax[1].set_ylabel(r"$\alpha$ ($10^{-6}$ K$^{-1}$)")
    qstyle.panel_label(ax[1], "(b)")

    ax[2].plot(res.T, res.Cv, color=c[0], lw=st["line"], label=r"$C_v$")
    ax[2].plot(res.T, res.Cp, color=c[2], lw=st["line"],
               dashes=[4, 1.6], label=r"$C_p$")
    ax[2].axhline(3.0 * res.natoms * KB_EV * 1000.0,
                  color=qstyle.INK_FAINT, lw=st["axis_line"],
                  dashes=[2, 2])
    ax[2].set_xlabel("T (K)"); ax[2].set_ylabel("meV/K por celda")
    ax[2].legend()
    qstyle.panel_label(ax[2], "(c)")

    written = qstyle.save(fig, outfile, formats, dpi=dpi,
                          modulo="cuasi-armónica")
    plt.close(fig)
    return written
