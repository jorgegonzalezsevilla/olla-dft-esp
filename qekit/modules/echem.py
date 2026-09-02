# Olla-DFT — command-line toolkit for Quantum ESPRESSO
# Copyright (C) 2026 Jorge Enrique González Sevilla
# SPDX-License-Identifier: GPL-3.0-or-later
# This program is free software: you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by the Free
# Software Foundation, either version 3 of the License, or (at your option)
# any later version. See the LICENSE file for details.

"""Electrodo de hidrógeno computacional: HER, OER y el diagrama de Pourbaix.

El truco del CHE (Nørskov 2004) es no calcular nunca un protón solvatado.
En el equilibrio del electrodo de hidrógeno estándar, a U = 0 V y pH = 0,

    ½ H₂(g)  ⇌  H⁺ + e⁻

tienen la MISMA energía libre. Así que cada vez que un paso libera un par
(H⁺ + e⁻), se le puede poner la energía de ½H₂, que sí se calcula. El
potencial y el pH entran después, como términos que se suman:

    ΔG(U, pH) = ΔG(0, 0) − n·e·U − n·k_B·T·ln(10)·pH

con n el número de pares transferidos y U medido frente al electrodo
estándar de hidrógeno (SHE): el término de pH es justo lo que convierte
la escala SHE en la RHE, U_RHE = U_SHE + 0.0592·pH, así que en la escala
RHE los ΔG no dependen del pH. A pH 0 las dos escalas coinciden. De ahí
salen las dos magnitudes que se citan:

  potencial limitante U_L = max(ΔG_i)/(n·e)   el potencial (vs RHE) al que
                                              TODOS los pasos se vuelven
                                              cuesta abajo
  sobrepotencial      η   = U_L − U_eq        lo que hay que aplicar de más

η se devuelve CON SIGNO: positivo quiere decir que al potencial de
equilibrio el paso limitante sigue cuesta arriba (lo normal; con los
perfiles de aquí nunca sale negativo porque el paso peor está siempre
por encima del promedio). Un η negativo solo puede aparecer con un
`dG_total` distinto del experimental y significaría que al equilibrio ya
todo va cuesta abajo.

**Lo que el CHE no incluye**, y conviene no olvidar: no hay barreras
cinéticas (solo termodinámica de los intermedios), no hay disolvente
explícito, no hay campo eléctrico de la doble capa, y la entropía del
adsorbato se toma como la de un sólido. Es un descriptor excelente para
comparar catalizadores entre sí y una predicción mediocre de una corriente.
"""

import math
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

from qekit.core import provenance
from qekit.core.errors import ErrorDeUso, FaltanDatos
from qekit.core import style as qstyle

KB_EV = 8.617333262e-5

# Energía libre de la reacción global 2H₂O -> O₂ + 2H₂, a 298 K y 1 bar.
# Se toma del EXPERIMENTO a propósito: el O₂ triplete es un caso conocido de
# error grande de los funcionales semilocales (unos 0.4 eV), así que anclar
# aquí el total es más fiable que calcularlo.
DG_AGUA_TOTAL = 4.92          # eV, para los cuatro pasos de la OER
U_EQ_OER = 1.229              # V frente al RHE
U_EQ_HER = 0.0

# Correcciones térmicas típicas (ZPE − TΔS) a 298 K de los intermedios
# sobre un metal de transición, en eV. Vienen de la literatura estándar
# (Nørskov y col.) y sirven de valor por omisión cuando no se han calculado
# las frecuencias; el reporte dice cuáles se usaron.
CORRECCIONES = {"H": 0.24, "OH": 0.35, "O": 0.05, "OOH": 0.40}

PASOS_OER = [
    ("H₂O + * → OH* + H⁺ + e⁻", "OH"),
    ("OH* → O* + H⁺ + e⁻", "O"),
    ("O* + H₂O → OOH* + H⁺ + e⁻", "OOH"),
    ("OOH* → * + O₂ + H⁺ + e⁻", None),
]


@dataclass
class Echem:
    reaccion: str = "her"
    energias: dict = field(default_factory=dict)     # intermedio -> E_ads (eV)
    correcciones: dict = field(default_factory=dict)
    T: float = 298.15
    pH: float = 0.0
    U: float = 0.0
    pasos: list = field(default_factory=list)        # [(nombre, dG en eV)]
    avisos: list = field(default_factory=list)

    @property
    def dG_H(self):
        """El descriptor de la HER: ΔG del hidrógeno adsorbido."""
        if "H" not in self.energias:
            return None
        return self.energias["H"] + self.correcciones.get("H", 0.0)

    def dG(self, U: float = None, pH: float = None) -> list:
        """ΔG de cada paso a un potencial (V vs SHE) y pH dados.

        El término de pH convierte U a la escala RHE, así que
        dG(U, pH) == dG(u_rhe(U, pH), 0).
        """
        U = self.U if U is None else U
        pH = self.pH if pH is None else pH
        desp = U + KB_EV * self.T * math.log(10.0) * pH
        return [(nombre, g - desp) for nombre, g in self.pasos]

    @property
    def limitante(self):
        """(nombre, ΔG) del paso que manda: el de mayor ΔG a U = 0."""
        if not self.pasos:
            return None, None
        return max(self.pasos, key=lambda t: t[1])

    @property
    def U_limitante(self):
        _, g = self.limitante
        return None if g is None else g          # un electrón por paso

    @property
    def sobrepotencial(self):
        """η = U_L − U_eq, con signo (positivo = cuesta arriba en U_eq)."""
        u = self.U_limitante
        if u is None:
            return None
        eq = U_EQ_OER if self.reaccion == "oer" else U_EQ_HER
        return u - eq

    def U_rhe(self, U: float = None, pH: float = None) -> float:
        """Potencial frente al RHE que corresponde a (U vs SHE, pH)."""
        U = self.U if U is None else U
        pH = self.pH if pH is None else pH
        return u_rhe(U, pH, self.T)


def u_rhe(U_she: float, pH: float, T: float = 298.15) -> float:
    """U_RHE = U_SHE + k_B·T·ln(10)·pH  (0.0592·pH a 298 K)."""
    return float(U_she) + KB_EV * T * math.log(10.0) * float(pH)


# ----------------------------------------------------------------------
# Construcción de los perfiles
# ----------------------------------------------------------------------
def her(E_ads_H: float, correccion: float = None, T: float = 298.15) -> Echem:
    """Perfil de la reacción de evolución de hidrógeno.

    Dos pasos, y el descriptor es que ΔG_H* esté cerca de cero: si el
    hidrógeno se pega demasiado poco no llega a adsorberse, y si se pega
    demasiado no se suelta. Es la cumbre del volcán de Nørskov, y en Pt(111)
    vale −0.09 eV.
    """
    c = CORRECCIONES["H"] if correccion is None else float(correccion)
    e = Echem(reaccion="her", T=T,
              energias={"H": float(E_ads_H)}, correcciones={"H": c})
    g = e.dG_H
    e.pasos = [("H⁺ + e⁻ + * → H*", g),
               ("H* → ½H₂ + *", -g)]
    if correccion is None:
        e.avisos.append(
            "La corrección térmica de H* (ZPE − TΔS = +0.24 eV) es el valor "
            "estándar de\n  la literatura, no una calculada para tu "
            "superficie. Para calcularla:\n  frecuencias del adsorbato y "
            "'olla-dft thermochem'.")
    return e


def oer(energias: dict, correcciones: dict = None,
        T: float = 298.15, dG_total: float = DG_AGUA_TOTAL) -> Echem:
    """Perfil de la evolución de oxígeno, con sus cuatro pasos.

    Las energías van referidas al agua y a la superficie limpia:
        E_ads(OH)  de  * + H₂O → OH* + ½H₂
        E_ads(O)   de  * + H₂O → O*  + H₂
        E_ads(OOH) de  * + 2H₂O → OOH* + 3/2 H₂
    que es lo que sale de restar medias moléculas de H₂ con el truco del CHE.

    El cuarto paso NO se calcula: se obtiene por diferencia con el total
    experimental de 4.92 eV. Es deliberado, porque calcular el O₂ triplete
    con un funcional semilocal se equivoca en unos 0.4 eV, y ese error
    entraría entero en el sobrepotencial.
    """
    faltan = [k for k in ("OH", "O", "OOH") if k not in energias]
    if faltan:
        raise ErrorDeUso(
            f"para la OER hacen falta las energías de OH, O y OOH; faltan "
            f"{', '.join(faltan)}. Sácalas de 'olla-dft adsorb' con cada "
            f"adsorbato, referidas al agua.")
    corr = dict(CORRECCIONES)
    corr.update(correcciones or {})
    e = Echem(reaccion="oer", T=T, energias=dict(energias), correcciones=corr)

    g = {k: float(energias[k]) + corr.get(k, 0.0) for k in ("OH", "O", "OOH")}
    dg1 = g["OH"]
    dg2 = g["O"] - g["OH"]
    dg3 = g["OOH"] - g["O"]
    dg4 = dG_total - (dg1 + dg2 + dg3)
    e.pasos = [(PASOS_OER[0][0], dg1), (PASOS_OER[1][0], dg2),
               (PASOS_OER[2][0], dg3), (PASOS_OER[3][0], dg4)]
    if dg4 < 0:
        e.avisos.append(
            f"El cuarto paso sale NEGATIVO ({dg4:+.2f} eV) por diferencia con "
            f"los 4.92 eV\n  totales. Quiere decir que la suma de los tres "
            f"calculados ya se pasa: o hay\n  un error en las referencias, o "
            f"tu superficie liga los intermedios muchísimo.")
    if correcciones is None:
        e.avisos.append(
            "Correcciones térmicas por omisión (OH 0.35, O 0.05, OOH 0.40 eV), "
            "de la\n  literatura estándar sobre metales de transición. Si tu "
            "superficie es otra cosa,\n  calcula las frecuencias y pásalas.")
    return e


# ----------------------------------------------------------------------
# Relación de escala de la OER
# ----------------------------------------------------------------------
#: ΔG(OOH*) − ΔG(OH*) universal (eV): 3.2 ± 0.2 en casi cualquier óxido.
ESCALA_OOH_OH = 3.2


def escala_ooh_oh(e: Echem) -> float:
    """ΔG(OOH*) − ΔG(OH*) del perfil, con sus correcciones térmicas.

    Es el número que se compara con la relación de escala universal
    (`ESCALA_OOH_OH`): si se sale mucho de 3.2 eV, o la superficie es
    especial o hay un error en las referencias.
    """
    if e.reaccion != "oer" or "OOH" not in e.energias or "OH" not in e.energias:
        raise FaltanDatos("la relación de escala OOH−OH solo tiene sentido "
                          "en un perfil de OER con OH y OOH.")
    return float(e.energias["OOH"] + e.correcciones.get("OOH", 0.0)
                 - e.energias["OH"] - e.correcciones.get("OH", 0.0))


def sobrepotencial_minimo_escala(delta: float = ESCALA_OOH_OH,
                                 dG_total: float = DG_AGUA_TOTAL) -> float:
    """η mínimo que impone la relación de escala: (Δ/2 − ΔG_total/4) V.

    Si OOH* y OH* están separados por Δ fijo, los pasos 2 y 3 (OH*→O* y
    O*→OOH*) suman Δ, y el peor de los dos no puede bajar de Δ/2; frente
    al equilibrio ΔG_total/4 = 1.23 V eso deja ~0.37 V con Δ = 3.2 eV.
    """
    return float(delta / 2.0 - dG_total / 4.0)


# ----------------------------------------------------------------------
# Reporte
# ----------------------------------------------------------------------
def report(e: Echem) -> str:
    nombre = {"her": "evolución de hidrógeno (HER)",
              "oer": "evolución de oxígeno (OER)"}[e.reaccion]
    L = [f"--- Electrodo de hidrógeno computacional: {nombre} ---",
         f"T = {e.T:.2f} K   |   U = {e.U:.2f} V vs SHE   |   pH = {e.pH:g}"
         + (f"   (= {e.U_rhe():.2f} V vs RHE)" if e.pH else ""),
         ""]
    L.append(f"  {'paso':44s} {'ΔG(0 V)':>10s} {'ΔG(U,pH)':>11s}")
    L.append("  " + "-" * 68)
    for (nom, g0), (_, gu) in zip(e.pasos, e.dG(e.U, e.pH)):
        L.append(f"  {nom:44s} {g0:>10.3f} {gu:>11.3f}")

    paso, gmax = e.limitante
    L += ["", f"Paso limitante: {paso}   (ΔG = {gmax:+.3f} eV)",
          f"Potencial limitante U_L = {e.U_limitante:+.3f} V vs RHE"]
    eq = U_EQ_OER if e.reaccion == "oer" else U_EQ_HER
    L.append(f"Sobrepotencial η = U_L − {eq:.3f} = "
             f"{e.sobrepotencial:+.3f} V   (positivo = en U_eq el paso "
             "limitante sigue cuesta arriba)")

    if e.reaccion == "her":
        g = e.dG_H
        L += ["", f"Descriptor ΔG_H* = {g:+.3f} eV"]
        if abs(g) < 0.10:
            L.append("  Muy cerca de cero: está en la cumbre del volcán, como "
                     "el Pt (−0.09 eV).")
        elif g < 0:
            L.append("  Negativo: el hidrógeno se pega demasiado y cuesta "
                     "soltarlo. La rama\n  izquierda del volcán; el paso "
                     "lento es la desorción.")
        else:
            L.append("  Positivo: el hidrógeno apenas se adsorbe. La rama "
                     "derecha del volcán;\n  el paso lento es la adsorción.")
        L.append("  Regla práctica: |ΔG_H*| < 0.2 eV es un buen catalizador "
                 "de HER.")
    else:
        eta = e.sobrepotencial
        if eta < 0.3:
            L.append("  η < 0.3 V es excelente; los mejores óxidos medidos "
                     "rondan 0.3 V.")
        elif eta > 0.8:
            L.append("  η > 0.8 V: poco prometedor. Mira si algún paso está "
                     "descompensado.")
        # relación de escala OOH-OH
        d = escala_ooh_oh(e)
        eta_min = sobrepotencial_minimo_escala()
        L += ["", f"Diferencia ΔG(OOH*) − ΔG(OH*) = {d:.3f} eV",
              f"  La relación de escala universal la fija en {ESCALA_OOH_OH} "
              "± 0.2 eV para casi cualquier\n  superficie, y de ahí sale el "
              f"límite de ~{eta_min:.2f} V en el sobrepotencial de la OER.\n"
              f"  Si tu número se sale mucho de {ESCALA_OOH_OH}, o has "
              "encontrado algo interesante o hay un\n  error en las "
              "referencias."]
        if abs(d - 3.2) > 0.5:
            L.append(f"  El tuyo está a {abs(d - 3.2):.2f} eV de 3.2: "
                     "revísalo antes de celebrarlo.")

    L += ["", "El CHE es termodinámica de intermedios: NO hay barreras "
              "cinéticas, ni\n  disolvente explícito, ni doble capa. Compara "
              "catalizadores muy bien y\n  predice corrientes regular."]
    for a in e.avisos:
        L.append(f"\nAVISO: {a}")
    return "\n".join(L)


def pourbaix(e: Echem, U=None, pH=None) -> dict:
    """ΔG del paso limitante en la rejilla (U, pH)."""
    U = np.linspace(-0.5, 2.0, 121) if U is None else np.asarray(U, float)
    pH = np.linspace(0.0, 14.0, 57) if pH is None else np.asarray(pH, float)
    g0 = np.array([g for _, g in e.pasos])
    # ΔG_i(U,pH) = ΔG_i(0,0) − eU − k_B T ln10 pH, igual para todos los pasos
    desp = U[None, :] + KB_EV * e.T * math.log(10.0) * pH[:, None]
    lim = g0.max() - desp                       # (npH, nU)
    return {"U": U, "pH": pH, "dG_limitante": lim}


def export(e: Echem, outdir: str = ".") -> list:
    out = Path(outdir); out.mkdir(parents=True, exist_ok=True)
    f = out / "ECHEM.dat"
    lines = [provenance.header(
        f"electrodo de hidrogeno computacional ({e.reaccion})",
        {"T_K": e.T, "U_V": e.U, "pH": e.pH,
         "U_limitante_V": e.U_limitante,
         "sobrepotencial_V": e.sobrepotencial}),
        f"# {'paso':46s} {'dG0(eV)':>10s}"]
    for nom, g in e.pasos:
        lines.append(f"  {nom:46s} {g:10.5f}")
    f.write_text("\n".join(lines) + "\n")
    txt = out / "ECHEM.txt"
    txt.write_text(report(e) + "\n")
    return [str(f), str(txt)]


def plot(e: Echem, outfile: str = "echem", formats="pdf,png",
         potenciales=None, theme: str = None, size: str = None,
         family: str = None, background: str = None, palette=None,
         usetex: bool = None, width="single", journal: str = "generic",
         aspect: float = 0.72, mono: bool = False, dpi: int = None) -> list:
    """Diagrama de energía libre en escalera, a varios potenciales."""
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError as exc:                              # pragma: no cover
        raise RuntimeError("matplotlib no está instalado.") from exc
    if not e.pasos:
        raise FaltanDatos("no hay pasos que graficar.")

    if potenciales is None:
        eq = U_EQ_OER if e.reaccion == "oer" else U_EQ_HER
        potenciales = [0.0, eq, e.U_limitante]
    potenciales = sorted(set(round(float(u), 3) for u in potenciales))

    st = qstyle.apply(theme, size=size, family=family, background=background,
                      palette=palette, usetex=usetex, mono=mono)
    fig, ax = qstyle.new_figure(width, journal, aspect)
    cols = qstyle.palette(max(len(potenciales), 3), mono=mono)

    n = len(e.pasos)
    for k, u in enumerate(potenciales):
        gs = [g for _, g in e.dG(u, e.pH)]
        acum = np.concatenate([[0.0], np.cumsum(gs)])
        x, y = [], []
        for i, val in enumerate(acum):
            x += [i - 0.35, i + 0.35]
            y += [val, val]
        ax.plot(x, y, lw=st["line"] * 1.6, color=cols[k % len(cols)],
                solid_capstyle="butt", label=f"U = {u:.2f} V")
        for i in range(n):
            ax.plot([i + 0.35, i + 1 - 0.35], [acum[i], acum[i + 1]],
                    lw=st["line"] * 0.8, color=cols[k % len(cols)],
                    dashes=[2.5, 2.0])
    ax.set_xticks(range(n + 1))
    ax.set_xlabel("coordenada de reacción")
    ax.set_ylabel(r"$\Delta G$ acumulada (eV)")
    ax.axhline(0.0, color=qstyle.INK_FAINT, lw=st["axis_line"],
               dashes=[3.5, 2.0])
    ax.legend(frameon=False, fontsize=st["legend"])
    written = qstyle.save(fig, outfile, formats, dpi=dpi, modulo="echem")
    plt.close(fig)
    return written
