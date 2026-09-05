# Olla-DFT — command-line toolkit for Quantum ESPRESSO
# Copyright (C) 2026 Jorge Enrique González Sevilla
# SPDX-License-Identifier: AGPL-3.0-or-later
# This program is free software: you can redistribute it and/or modify it
# under the terms of the GNU Affero General Public License as published by the Free
# Software Foundation, either version 3 of the License, or (at your option)
# any later version. See the LICENSE file for details.

"""Fonones a temperatura electrónica: ¿se estabiliza el modo imaginario?

Un modo con frecuencia imaginaria dice que la estructura no está en un
mínimo. En muchos materiales eso no es un error del cálculo sino física: la
fase de alta simetría es inestable a T = 0 y se estabiliza al calentar, que
es como funcionan las ondas de densidad de carga y muchas transiciones
estructurales. Para verlo hay que repetir los fonones a varias temperaturas
electrónicas y mirar cómo se mueve el modo.

**El punto que hay que tener claro**: en pw.x la temperatura electrónica se
mete a través del ensanchamiento, y solo el ensanchamiento de FERMI-DIRAC es
una temperatura de verdad,

    degauss = k_B · T,      k_B = 6.33362e-6 Ry/K

Los otros (gaussiano, Methfessel-Paxton, Marzari-Vanderbilt "cold") son
trucos numéricos para converger antes la integral sobre la superficie de
Fermi: dan una ocupación que NO es la de Fermi-Dirac y su anchura no
corresponde a ninguna temperatura. Un barrido de degauss con smearing cold
produce una gráfica con pinta de "frecuencia contra temperatura" que no lo
es. Por eso este módulo impone fermi-dirac y lo dice.

Comprobado en silicio (Γ, LDA): el modo óptico sale a 507 cm⁻¹ a 300 K
(experimental 520) y se ablanda hasta 464 cm⁻¹ a 6000 K. No es un artefacto:
a esa temperatura k_B·T = 0.52 eV supera el gap LDA, se excitan portadores
a través de él y los enlaces se debilitan. Es el mismo ablandamiento que se
mide tras excitar silicio con un láser intenso.
"""

from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

from qekit.core import provenance
from qekit.core.errors import ErrorDeUso, FaltanDatos
from qekit.core import style as qstyle

# k_B en Ry/K
KB_RY = 6.333621e-6
# por debajo de esto una frecuencia imaginaria es ruido numérico, no física
UMBRAL_IMAGINARIO = 10.0        # cm^-1


@dataclass
class BarridoT:
    temperaturas: list = field(default_factory=list)     # K
    carpetas: list = field(default_factory=list)
    frecuencias: dict = field(default_factory=dict)      # T -> array cm^-1
    gamma_only: bool = True
    smearing: str = "fermi-dirac"
    avisos: list = field(default_factory=list)

    def imaginarias(self, T) -> np.ndarray:
        f = self.frecuencias.get(T)
        if f is None:
            return np.array([])
        return f[f < -UMBRAL_IMAGINARIO]

    def modo_blando(self, T):
        """Frecuencia CON SIGNO del modo más blando, sin contar los acústicos.

        Se sigue esta y no "la peor imaginaria, o cero si no hay": con esa
        definición la curva se aplana en cuanto la estructura se estabiliza,
        y la temperatura de estabilización salía siempre pegada al punto del
        barrido donde desaparecía el modo, no entre dos puntos.
        """
        f = self.frecuencias.get(T)
        if f is None or not len(f):
            return None
        util = f[np.abs(f) > UMBRAL_IMAGINARIO]
        if not len(util):
            return None
        return float(np.min(util))

    @property
    def con_datos(self) -> list:
        return [T for T in self.temperaturas if self.frecuencias.get(T) is not None]


def degauss_de_T(T: float) -> float:
    """Ensanchamiento de Fermi-Dirac en Ry que corresponde a T kelvin."""
    if T <= 0:
        raise ErrorDeUso(
            f"la temperatura electrónica tiene que ser positiva; recibí {T}.")
    return KB_RY * float(T)


def T_de_degauss(degauss: float) -> float:
    return float(degauss) / KB_RY


def prepare(atoms, temperaturas, outdir: str = "fonones_T",
            gamma_only: bool = True, **kw) -> tuple:
    """Un cálculo de fonones completo por cada temperatura."""
    from qekit.modules import phonons

    temperaturas = sorted({float(T) for T in temperaturas})
    if len(temperaturas) < 2:
        raise ErrorDeUso(
            "un barrido de temperatura necesita al menos dos valores; con uno "
            "solo no hay nada que comparar. Prueba --tscan 300,1000,3000.")
    if max(temperaturas) > 20000:
        raise ErrorDeUso(
            f"{max(temperaturas):g} K de temperatura ELECTRÓNICA es enorme: "
            "el ensanchamiento se come varios eV y las ocupaciones dejan de "
            "parecerse a nada. Los estudios de ondas de densidad de carga "
            "llegan a unos 6000 K.")

    out = Path(outdir); out.mkdir(parents=True, exist_ok=True)
    run = BarridoT(temperaturas=temperaturas, gamma_only=gamma_only)
    partes = []
    for T in temperaturas:
        sub = out / f"T{int(round(T)):05d}"
        _, rep = phonons.prepare(
            atoms, outdir=str(sub), gamma_only=gamma_only,
            insulator=False,                 # una temperatura electrónica
            degauss=degauss_de_T(T),         # exige ocupaciones con smearing
            smearing="fermi-dirac", **kw)
        run.carpetas.append(str(sub))
        partes.append((T, rep))

    report = ["--- Fonones a temperatura electrónica ---",
              "Temperaturas: "
              + ", ".join(f"{T:g} K" for T in temperaturas),
              "Ensanchamiento: fermi-dirac, degauss = k_B·T = "
              + ", ".join(f"{degauss_de_T(T):.5f}" for T in temperaturas)
              + " Ry",
              "",
              "Se impone smearing='fermi-dirac' a propósito: es el único cuya "
              "anchura ES una\n  temperatura. Con gaussiano o cold el "
              "ensanchamiento es un truco numérico y\n  la curva "
              "frecuencia-contra-temperatura no significaría nada.",
              "",
              f"{len(temperaturas)} cálculos de fonones en "
              f"'{out.resolve()}'.",
              "Cada uno es la cadena entera (scf, ph.x y post-proceso): "
              "esto cuesta N veces\n  un cálculo de fonones normal."]
    if not gamma_only:
        report.append(
            "AVISO: barrido de temperatura con malla de q completa. Es lo "
            "correcto si el\n  modo blando NO está en Γ (una onda de densidad "
            "de carga casi nunca lo está),\n  pero multiplica un cálculo ya "
            "caro por el número de temperaturas.")
    return run, "\n".join(report)


def collect(run: BarridoT) -> BarridoT:
    from qekit.modules import phonons

    for T, carpeta in zip(run.temperaturas, run.carpetas):
        try:
            pr = phonons.PhononRun(prefix="", outdir=Path(carpeta),
                                   gamma_only=run.gamma_only)
            pr = phonons.collect(pr)
        except Exception:                                   # noqa: BLE001
            run.frecuencias[T] = None
            continue
        if run.gamma_only and pr.gamma_freqs:
            f = np.array([g[0] for g in pr.gamma_freqs], dtype=float)
        elif pr.band_freqs is not None and len(pr.band_freqs):
            f = np.asarray(pr.band_freqs, dtype=float).ravel()
        else:
            f = None
        run.frecuencias[T] = f
    return run


def temperatura_de_estabilizacion(run: BarridoT) -> float:
    """T a la que el modo blando cruza el cero, interpolando linealmente."""
    Ts = run.con_datos
    if len(Ts) < 2:
        return None
    blandos = [(T, run.modo_blando(T)) for T in Ts]
    blandos = [(T, w) for T, w in blandos if w is not None]
    for (T1, w1), (T2, w2) in zip(blandos, blandos[1:]):
        if w1 < 0.0 <= w2:
            if w2 == w1:
                return float(T2)
            return float(T1 + (0.0 - w1) * (T2 - T1) / (w2 - w1))
    return None


def monotono(run: BarridoT) -> bool:
    """¿El número de modos imaginarios solo baja al subir la temperatura?"""
    cuentas = [len(run.imaginarias(T)) for T in run.con_datos]
    return all(b <= a for a, b in zip(cuentas, cuentas[1:]))


def report(run: BarridoT) -> str:
    Ts = run.con_datos
    if not Ts:
        raise FaltanDatos(
            "no hay frecuencias que leer todavía. Corre las cadenas de "
            "fonones y vuelve con --collect.")
    L = ["--- Fonones contra temperatura electrónica ---",
         f"Ensanchamiento fermi-dirac; {'solo Γ' if run.gamma_only else 'malla de q completa'}",
         "",
         f"  {'T (K)':>8s} {'degauss (Ry)':>13s} {'modos imag.':>12s} "
         f"{'peor (cm⁻¹)':>13s} {'menor real':>12s}"]
    L.append("  " + "-" * 62)
    for T in run.temperaturas:
        f = run.frecuencias.get(T)
        if f is None:
            L.append(f"  {T:>8.0f} {degauss_de_T(T):>13.5f} "
                     f"{'sin resultado':>12s}")
            continue
        im = run.imaginarias(T)
        reales = f[f >= -UMBRAL_IMAGINARIO]
        L.append(f"  {T:>8.0f} {degauss_de_T(T):>13.5f} {len(im):>12d} "
                 + (f"{im.min():>13.2f}" if len(im) else f"{'—':>13s}")
                 + (f" {reales.min():>12.2f}" if len(reales) else
                    f" {'—':>12s}"))

    T_est = temperatura_de_estabilizacion(run)
    L.append("")
    if not monotono(run):
        L.append("El número de modos imaginarios NO baja de forma monótona "
                 "con la temperatura.\n  Suele querer decir que falta "
                 "convergencia en la malla de k: con ensanchamiento\n  "
                 "pequeño hace falta malla fina, y si no, el ruido se "
                 "confunde con el modo blando.\n  Cualquier temperatura de "
                 "estabilización sacada de aquí sería inventada.")
    elif T_est is not None:
        L += [f"La estructura se estabiliza alrededor de {T_est:.0f} K: por "
              f"debajo hay modos\n  imaginarios y por encima no.",
              "  Es la firma de una transición estructural o de una onda de "
              "densidad de carga.\n  El número sale de interpolar entre dos "
              "puntos del barrido, así que su\n  precisión es la del paso que "
              "hayas usado."]
    elif all(len(run.imaginarias(T)) == 0 for T in Ts):
        L.append("No hay modos imaginarios a ninguna temperatura: la "
                 "estructura es estable\n  en todo el rango barrido.")
    elif all(len(run.imaginarias(T)) > 0 for T in Ts):
        peor_alta = run.imaginarias(Ts[-1])
        L += [f"Sigue habiendo {len(peor_alta)} modo(s) imaginario(s) a "
              f"{Ts[-1]:.0f} K, el peor en\n  {peor_alta.min():.1f} cm⁻¹. "
              "O la inestabilidad no es de origen electrónico (y no se cura\n"
              "  calentando), o hace falta subir más la temperatura, o la "
              "celda es demasiado\n  pequeña para contener la distorsión que "
              "el material quiere hacer."]
    else:
        L.append("Hay modos imaginarios en parte del rango, pero el barrido "
                 "no llega a cruzar\n  el cero. Añade temperaturas más altas.")

    L += ["",
          "Recordatorio: esto es temperatura ELECTRÓNICA. Los iones siguen "
          "estando quietos\n  en sus posiciones de equilibrio; no hay "
          "movimiento térmico ni dilatación. Para\n  eso hacen falta dinámica "
          "molecular (olla-dft gen -p md) o la cuasi-armónica (olla-dft qha)."]
    return "\n".join(L)


def export(run: BarridoT, outdir: str = ".") -> list:
    out = Path(outdir); out.mkdir(parents=True, exist_ok=True)
    f = out / "FONONES_T.dat"
    lines = [provenance.header(
        "fonones contra temperatura electronica",
        {"smearing": "fermi-dirac", "gamma_only": run.gamma_only,
         "T_estabilizacion_K": temperatura_de_estabilizacion(run)}),
        f"# {'T(K)':>10s} {'degauss(Ry)':>14s} {'n_imag':>8s} "
        f"{'peor(cm-1)':>13s}"]
    for T in run.temperaturas:
        fr = run.frecuencias.get(T)
        if fr is None:
            continue
        im = run.imaginarias(T)
        lines.append(f"{T:12.2f} {degauss_de_T(T):14.6f} {len(im):8d} "
                     + (f"{im.min():13.3f}" if len(im) else f"{'0.000':>13s}"))
    f.write_text("\n".join(lines) + "\n")
    txt = out / "FONONES_T.txt"
    txt.write_text(report(run) + "\n")
    return [str(f), str(txt)]


def plot(run: BarridoT, outfile: str = "fonones_T", formats="pdf,png",
         theme: str = None, size: str = None, family: str = None,
         background: str = None, palette=None, usetex: bool = None,
         width="single", journal: str = "generic", aspect: float = 0.75,
         mono: bool = False, dpi: int = None) -> list:
    """Frecuencias contra temperatura, con la zona imaginaria sombreada."""
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError as exc:                              # pragma: no cover
        raise RuntimeError("matplotlib no está instalado.") from exc
    Ts = run.con_datos
    if len(Ts) < 2:
        raise FaltanDatos("hacen falta al menos dos temperaturas con datos.")

    st = qstyle.apply(theme, size=size, family=family, background=background,
                      palette=palette, usetex=usetex, mono=mono)
    fig, ax = qstyle.new_figure(width, journal, aspect)
    cols = qstyle.palette(2, mono=mono)

    nmax = max(len(run.frecuencias[T]) for T in Ts)
    for i in range(nmax):
        y = [run.frecuencias[T][i] if i < len(run.frecuencias[T]) else np.nan
             for T in Ts]
        blando = any(v < -UMBRAL_IMAGINARIO for v in y if v == v)
        ax.plot(Ts, y, marker="o", ms=3.5, lw=st["line"],
                color=cols[0] if blando else qstyle.INK_FAINT,
                zorder=3 if blando else 1)
    ax.axhline(0.0, color=qstyle.INK_FAINT, lw=st["axis_line"],
               dashes=[3.5, 2.0])
    lo = ax.get_ylim()[0]
    ax.axhspan(lo, 0.0, color=cols[0], alpha=0.07, lw=0)
    ax.annotate("imaginarias", xy=(Ts[0], lo), xytext=(3, 4),
                textcoords="offset points", fontsize=st["legend"],
                color=qstyle.INK_SOFT)
    T_est = temperatura_de_estabilizacion(run)
    if T_est is not None:
        ax.axvline(T_est, color=qstyle.INK_FAINT, lw=st["axis_line"],
                   dashes=[1.5, 1.5])
        ax.annotate(f"{T_est:.0f} K", xy=(T_est, 0.0), xytext=(4, 6),
                    textcoords="offset points", fontsize=st["legend"],
                    color=qstyle.INK_SOFT)
    ax.set_xlabel("temperatura electrónica (K)")
    ax.set_ylabel(r"frecuencia (cm$^{-1}$)")
    written = qstyle.save(fig, outfile, formats, dpi=dpi, modulo="fonones_T")
    plt.close(fig)
    return written
