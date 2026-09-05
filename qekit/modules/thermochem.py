# Olla-DFT — command-line toolkit for Quantum ESPRESSO
# Copyright (C) 2026 Jorge Enrique González Sevilla
# SPDX-License-Identifier: AGPL-3.0-or-later
# This program is free software: you can redistribute it and/or modify it
# under the terms of the GNU Affero General Public License as published by the Free
# Software Foundation, either version 3 of the License, or (at your option)
# any later version. See the LICENSE file for details.

"""Termoquímica: de una energía DFT a una energía libre comparable.

EL PROBLEMA
-----------
Una energía de adsorción o una barrera de DFT es una diferencia de
energías ELECTRÓNICAS a 0 K. Lo que se mide en un laboratorio es una
energía libre a la temperatura y la presión del experimento. Entre las
dos hay tres cosas:

    G(T,p) = E_DFT + ZPE + [H(T) - H(0)] - T*S(T)

- **ZPE**, la energía de punto cero: incluso a 0 K los modos vibran.
  Para un enlace X-H son varias décimas de eV, y como el número de modos
  cambia al adsorber, NO se cancela en la diferencia.
- **La corrección entálpica**, que puebla los modos a temperatura T.
- **-T*S**, el término entrópico. Es el que más pesa y el que más se
  olvida: una molécula en fase gas pierde casi toda su entropía
  traslacional y rotacional al adsorberse, y a 500 K eso vale del orden
  de 1 eV. Ignorarlo puede cambiar el signo de una energía de adsorción.

LOS DOS REGÍMENES
-----------------
- **Sólido / adsorbato**: todos los grados de libertad son vibracionales.
  Se usa el oscilador armónico para los 3N modos.
- **Gas ideal**: 3 traslaciones, 2 o 3 rotaciones y 3N-6 (o 3N-5)
  vibraciones. Las traslaciones y rotaciones se tratan con sus fórmulas
  clásicas, que necesitan la masa, los momentos de inercia y el número
  de simetría del grupo puntual.

MODOS IMAGINARIOS Y MODOS BLANDOS
---------------------------------
Una frecuencia imaginaria en un mínimo quiere decir que no es un mínimo.
En un estado de transición TIENE que haber exactamente una, y esa NO
entra en la suma. Olla-DFT cuenta las imaginarias y se comporta según lo que
se le diga que es la estructura.

Los modos muy blandos (por debajo de unos 100 cm⁻¹) son el otro problema:
la entropía vibracional de un modo diverge como -ln(w) cuando w->0, así
que un modo de 5 cm⁻¹ mal calculado mete un error enorme. La práctica
habitual es subirlos a un piso; Olla-DFT lo hace si se le pide y dice
cuántos modos tocó.
"""

from dataclasses import dataclass, field

import numpy as np

from qekit.core.errors import ErrorDeUso

# constantes
H_EVS = 4.135667696e-15        # h en eV*s
KB_EV = 8.617333262e-5         # k_B en eV/K
C_CM = 2.99792458e10           # c en cm/s
AMU_KG = 1.66053906660e-27
HBAR_JS = 1.054571817e-34
KB_J = 1.380649e-23
EV_J = 1.602176634e-19
NA = 6.02214076e23

#: Piso al que se suben los modos blandos, en cm^-1.
PISO_BLANDO = 100.0


def cm1_a_eV(nu):
    """cm^-1 -> eV."""
    return np.asarray(nu, dtype=float) * C_CM * H_EVS


@dataclass
class Termoquimica:
    ZPE: float = None            # eV
    H_corr: float = None         # eV, H(T)-H(0) sin ZPE
    S: float = None              # eV/K
    TS: float = None             # eV
    G_corr: float = None         # eV, ZPE + H_corr - T*S
    Cv: float = None             # eV/K
    T: float = None
    p: float = None              # Pa
    fase: str = ""
    n_imaginarias: int = 0
    n_subidos: int = 0
    S_trans: float = None
    S_rot: float = None
    S_vib: float = None
    S_elec: float = None
    avisos: list = field(default_factory=list)


# ----------------------------------------------------------------------
# Vibracional
# ----------------------------------------------------------------------
def limpiar_frecuencias(nu, fase: str = "solido", piso: float = None,
                        tol_imaginaria: float = 1.0) -> tuple:
    """Quita las que no cuentan y sube las blandas.

    Devuelve (frecuencias usables en cm^-1, n_imaginarias, n_subidas,
    avisos).

    Una frecuencia imaginaria se representa como negativa. En un mínimo no
    debería haber ninguna; en un estado de transición, exactamente una. En
    los dos casos se EXCLUYE de las sumas termodinámicas, porque la
    fórmula del oscilador armónico no está definida para ella.
    """
    nu = np.asarray(nu, dtype=float)
    avisos = []
    imag = nu < -tol_imaginaria
    n_imag = int(imag.sum())
    reales = nu[~imag]
    # las casi-cero (traslaciones y rotaciones residuales) tampoco cuentan
    reales = reales[reales > tol_imaginaria]

    if fase == "transicion":
        if n_imag == 0:
            avisos.append(
                "Se declaró un ESTADO DE TRANSICIÓN pero no hay ninguna "
                "frecuencia imaginaria.\nUn estado de transición es un punto "
                "de silla de primer orden: tiene que\ntener exactamente una. "
                "O la estructura no es el estado de transición, o el\ncálculo "
                "de fonones no está convergido.")
        elif n_imag > 1:
            avisos.append(
                f"Hay {n_imag} frecuencias imaginarias en un estado de "
                "transición; debería\nhaber exactamente una. Con más de una "
                "es un punto de silla de orden mayor,\nno un estado de "
                "transición.")
    elif n_imag > 0:
        avisos.append(
            f"Hay {n_imag} frecuencia(s) imaginaria(s) en algo declarado como "
            f"MÍNIMO\n(la mayor: {abs(nu[imag]).max():.1f}i cm⁻¹). Eso quiere "
            "decir que la estructura NO\nes un mínimo: relaja mejor antes de "
            "calcular termoquímica sobre ella.\nSe excluyen de las sumas, "
            "pero el número que salga no describe un estado\nestable.")

    n_sub = 0
    if piso:
        blandos = reales < piso
        n_sub = int(blandos.sum())
        reales = np.where(blandos, piso, reales)
        if n_sub:
            avisos.append(
                f"Se subieron {n_sub} modo(s) por debajo de {piso:.0f} cm⁻¹ "
                f"hasta ese piso.\nLa entropía vibracional de un modo diverge "
                "como -ln(w) cuando w tiende a 0,\nasí que un modo blando mal "
                "calculado domina el resultado. Subirlos es la\npráctica "
                "habitual, pero es una CORRECCIÓN, no un cálculo: dilo si "
                "publicas\nestos números.")
    return reales, n_imag, n_sub, avisos


def zpe(nu_cm1) -> float:
    """Energía de punto cero: la mitad de la suma de los cuantos."""
    return float(0.5 * np.sum(cm1_a_eV(nu_cm1)))


def _x(nu_cm1, T):
    return cm1_a_eV(nu_cm1) / (KB_EV * T)


def H_vib(nu_cm1, T: float) -> float:
    """H(T)-H(0) vibracional (sin ZPE), en eV."""
    if T <= 0:
        return 0.0
    e = cm1_a_eV(nu_cm1)
    x = _x_acotado(nu_cm1, T)
    return float(np.sum(e / np.expm1(x)))


#: Por encima de este x = h*nu/(k_B*T) el modo esta congelado y sus
#: formulas desbordan; su contribucion es cero a la precision de un float.
X_MAX = 500.0


def _x_acotado(nu_cm1, T):
    """x = h*nu/(k_B*T), acotado para que exp(x) no desborde.

    A T = 1 K un modo de 1000 cm-1 da x ~ 1440: exp(1440) es infinito y
    numpy avisa. El modo esta congelado y aporta cero, asi que acotar no
    cambia el resultado — solo quita un aviso que asusta sin motivo.
    """
    return np.minimum(_x(nu_cm1, T), X_MAX)


def S_vib(nu_cm1, T: float) -> float:
    """Entropía vibracional en eV/K."""
    if T <= 0:
        return 0.0
    x = _x_acotado(nu_cm1, T)
    return float(KB_EV * np.sum(x / np.expm1(x) - np.log1p(-np.exp(-x))))


def Cv_vib(nu_cm1, T: float) -> float:
    if T <= 0:
        return 0.0
    x = _x_acotado(nu_cm1, T)
    ex = np.exp(x)
    return float(KB_EV * np.sum(x ** 2 * ex / np.expm1(x) ** 2))


# ----------------------------------------------------------------------
# Gas ideal: traslación y rotación
# ----------------------------------------------------------------------
def S_traslacional(masa_amu: float, T: float, p: float = 101325.0) -> float:
    """Sackur-Tetrode, en eV/K.

    Depende de la PRESIÓN: es la razón por la que una energía de adsorción
    "a 1 bar" y otra "a la presión parcial del reactor" no son el mismo
    número.
    """
    m = masa_amu * AMU_KG
    V = KB_J * T / p                      # volumen por molécula
    lam = np.sqrt(2 * np.pi * m * KB_J * T) / (2 * np.pi * HBAR_JS)
    S_J = KB_J * (np.log(lam ** 3 * V) + 2.5)
    return float(S_J / EV_J)


def momentos_inercia(atoms) -> np.ndarray:
    """Momentos principales de inercia en amu*Å²."""
    m = atoms.get_masses()
    pos = atoms.get_positions() - np.average(atoms.get_positions(), axis=0,
                                             weights=m)
    I = np.zeros((3, 3))
    for mi, r in zip(m, pos):
        I += mi * (np.dot(r, r) * np.eye(3) - np.outer(r, r))
    return np.sort(np.linalg.eigvalsh(I))


def S_rotacional(atoms, T: float, simetria: int = 1) -> float:
    """Entropía rotacional del rotor rígido clásico, en eV/K.

    `simetria` es el número de simetría del grupo puntual: 2 para H2O y
    O2, 3 para NH3, 12 para CH4 y benceno. Olvidarlo sobreestima la
    entropía en k_B*ln(sigma), que para el metano son 0.06 eV a 300 K.
    """
    if len(atoms) < 2:
        return 0.0
    I = momentos_inercia(atoms) * AMU_KG * 1e-20      # kg*m^2
    lineal = I[0] < 1e-3 * I[2]
    if lineal:
        Ib = I[2]
        theta = HBAR_JS ** 2 / (2 * Ib * KB_J)
        S_J = KB_J * (np.log(T / (simetria * theta)) + 1.0)
    else:
        thetas = HBAR_JS ** 2 / (2 * I * KB_J)
        S_J = KB_J * (0.5 * np.log(np.pi * T ** 3 /
                                   (simetria ** 2 * np.prod(thetas))) + 1.5)
    return float(S_J / EV_J)


# ----------------------------------------------------------------------
# Todo junto
# ----------------------------------------------------------------------
def corregir(nu_cm1, T: float = 298.15, fase: str = "solido",
             atoms=None, p: float = 101325.0, simetria: int = 1,
             multiplicidad: int = 1, piso: float = None) -> Termoquimica:
    """Correcciones térmicas completas para una estructura.

    `fase`:
      - 'solido' / 'adsorbato': solo vibraciones (3N modos)
      - 'gas': traslaciones + rotaciones + vibraciones
      - 'transicion': como sólido, pero exigiendo una imaginaria
    """
    if fase not in ("solido", "adsorbato", "gas", "transicion"):
        raise ErrorDeUso(
            f"fase '{fase}' desconocida. Opciones: solido, adsorbato, gas, "
            "transicion.")
    if T <= 0:
        raise ErrorDeUso("la temperatura tiene que ser positiva.")

    reales, n_imag, n_sub, avisos = limpiar_frecuencias(
        nu_cm1, fase=fase, piso=piso)
    if len(reales) == 0:
        raise ErrorDeUso(
            "no queda ninguna frecuencia real utilizable. Revisa el cálculo "
            "de fonones.")

    tq = Termoquimica(T=T, p=p, fase=fase, n_imaginarias=n_imag,
                      n_subidos=n_sub, avisos=avisos)
    tq.ZPE = zpe(reales)
    tq.S_vib = S_vib(reales, T)
    tq.Cv = Cv_vib(reales, T)
    H = H_vib(reales, T)
    S = tq.S_vib

    if fase == "gas":
        if atoms is None:
            raise ErrorDeUso(
                "para la fase gas hace falta la estructura (masa y momentos "
                "de inercia).")
        tq.S_trans = S_traslacional(float(sum(atoms.get_masses())), T, p)
        tq.S_rot = S_rotacional(atoms, T, simetria)
        S += tq.S_trans + tq.S_rot
        # traslación 3/2 kT, rotación 3/2 kT (o kT si es lineal), + pV = kT
        lineal = len(atoms) > 1 and \
            momentos_inercia(atoms)[0] < 1e-3 * momentos_inercia(atoms)[2]
        n_rot = 1.0 if lineal else 1.5
        if len(atoms) == 1:
            n_rot = 0.0
        H += (1.5 + n_rot + 1.0) * KB_EV * T
        n_esperados = 3 * len(atoms) - (5 if lineal else 6)
        if len(reales) != n_esperados and len(atoms) > 1:
            tq.avisos.append(
                f"Para una molécula {'lineal' if lineal else 'no lineal'} de "
                f"{len(atoms)} átomos se esperan\n{n_esperados} modos "
                f"vibracionales y hay {len(reales)}. Si sobran, seguramente "
                "son\ntraslaciones y rotaciones residuales que no se "
                "separaron: cuentan doble\ncon los términos traslacional y "
                "rotacional.")

    if multiplicidad > 1:
        tq.S_elec = KB_EV * np.log(multiplicidad)
        S += tq.S_elec

    tq.S = S
    tq.H_corr = H
    tq.TS = T * S
    tq.G_corr = tq.ZPE + H - T * S
    return tq


def report(tq: Termoquimica, E_dft: float = None) -> str:
    lines = ["--- Correcciones termoquímicas ---",
             f"Fase: {tq.fase}    T = {tq.T:.2f} K"
             + (f"    p = {tq.p / 1e5:.4g} bar" if tq.fase == "gas" else ""),
             ""]
    lines.append(f"  ZPE (punto cero)          {tq.ZPE:+10.4f} eV")
    lines.append(f"  H(T) - H(0)               {tq.H_corr:+10.4f} eV")
    lines.append(f"  -T*S                      {-tq.TS:+10.4f} eV")
    lines.append(f"  {'-' * 40}")
    lines.append(f"  G(T) - E_DFT              {tq.G_corr:+10.4f} eV")
    if E_dft is not None:
        lines += ["",
                  f"  E_DFT                     {E_dft:+10.4f} eV",
                  f"  G(T)                      {E_dft + tq.G_corr:+10.4f} eV"]
    lines += ["", "Desglose de la entropía (eV/K y su contribución -T*S):"]
    for nombre, val in (("vibracional", tq.S_vib),
                        ("traslacional", tq.S_trans),
                        ("rotacional", tq.S_rot),
                        ("electrónica", tq.S_elec)):
        if val is None:
            continue
        lines.append(f"  {nombre:15s} {val:12.3e}   {-tq.T * val:+8.4f} eV")
    lines.append(f"  {'total':15s} {tq.S:12.3e}   {-tq.TS:+8.4f} eV")
    lines.append("")
    lines.append(f"  C_v vibracional           {tq.Cv:.4e} eV/K  "
                 f"({tq.Cv * NA * EV_J:.2f} J/(mol·K))")

    if tq.fase == "gas":
        lines += ["",
                  "La entropía traslacional depende de la PRESIÓN y la "
                  "rotacional del número\nde simetría. Cambiar de 1 bar a la "
                  "presión parcial de un reactor mueve\nG en décimas de eV; "
                  "olvidar el número de simetría del metano (12) la\nmueve "
                  "0.06 eV a 300 K."]
    for a in tq.avisos:
        lines += ["", a]
    return "\n".join(lines)


def adsorcion(E_slab_ads: float, E_slab: float, E_gas: float,
              tq_ads: Termoquimica = None, tq_gas: Termoquimica = None,
              n: int = 1) -> dict:
    """Energía de adsorción, electrónica y libre.

        E_ads = E(slab+ads) - E(slab) - n*E(gas)

    Negativa quiere decir favorable. La versión libre suma las
    correcciones del adsorbato y RESTA las del gas: la molécula pierde su
    entropía traslacional y rotacional al pegarse, y ese término suele
    ser el que decide.
    """
    E = E_slab_ads - E_slab - n * E_gas
    fuera = {"E_ads": E}
    if tq_ads is not None and tq_gas is not None:
        fuera["G_ads"] = E + tq_ads.G_corr - n * tq_gas.G_corr
        fuera["dZPE"] = tq_ads.ZPE - n * tq_gas.ZPE
        fuera["dTS"] = tq_ads.TS - n * tq_gas.TS
        fuera["T"] = tq_ads.T
    return fuera
