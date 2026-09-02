# Olla-DFT — command-line toolkit for Quantum ESPRESSO
# Copyright (C) 2026 Jorge Enrique González Sevilla
# SPDX-License-Identifier: GPL-3.0-or-later
# This program is free software: you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by the Free
# Software Foundation, either version 3 of the License, or (at your option)
# any later version. See the LICENSE file for details.

"""Cantidades derivadas de las constantes elásticas y de los fonones.

Post-proceso puro: de las Cij que ya calcula `olla-dft elastic` salen las
velocidades del sonido, la temperatura de Debye, el parámetro de Grüneisen
y una estimación de la conductividad térmica de red. No cuesta ningún
cálculo nuevo, y son justo las cantidades que permiten CRUZAR el módulo
elástico contra el de fonones — dos rutas físicamente independientes al
mismo número.

    v_l = sqrt((B + 4G/3)/rho)      longitudinal
    v_t = sqrt(G/rho)               transversal
    v_m = [ (1/3)(2/v_t^3 + 1/v_l^3) ]^(-1/3)    promedio
    theta_D = (hbar/k_B) [6 pi^2 n]^(1/3) v_m

CUIDADO CON LA TEMPERATURA DE DEBYE
-----------------------------------
No hay UNA temperatura de Debye: hay varias definiciones que dan números
distintos para el mismo material.

- la **elástica** sale de las velocidades del sonido: es el límite de baja
  temperatura, donde solo cuentan las ramas acústicas;
- la que sale de la **DOS de fonones** por momentos usa TODO el espectro,
  ópticas incluidas.

Comparar las dos es útil —si difieren mucho, algo va mal— pero esperar que
coincidan al 1 % es un error: describen cosas distintas. Olla-DFT reporta las
dos, dice cuál es cuál, y usa una tolerancia amplia al cruzarlas.
"""

from dataclasses import dataclass, field

import numpy as np

HBAR = 1.054571817e-34          # J*s
KB = 1.380649e-23               # J/K
AMU = 1.66053906660e-27         # kg
NA = 6.02214076e23
# Factor para pasar una pendiente dnu/dq [cm^-1 / A^-1] a m/s.
#   w = 2*pi*c*nu   con c = 2.99792458e10 cm/s y nu en cm^-1
#   v = dw/dq = 2*pi*c*(dnu/dq)  -> (cm/s)*(cm^-1*A) = A/s -> *1e-10 = m/s
#   => v[m/s] = 2*pi * 2.99792458e10 * 1e-10 * (dnu/dq) = 18.836*(dnu/dq)
CM1_A_A_MS = 2.0 * np.pi * 2.99792458e10 * 1e-10


@dataclass
class Termoelastico:
    rho: float = None            # kg/m^3
    v_l: float = None            # m/s
    v_t: float = None
    v_m: float = None
    theta_D: float = None        # K, definición elástica
    gruneisen: float = None      # de la razón de Poisson
    poisson: float = None
    kappa_slack: float = None    # W/(m*K) a la temperatura T
    T: float = 300.0             # K, temperatura de la kappa de Slack
    natoms: int = 0
    volumen: float = None        # A^3
    avisos: list = field(default_factory=list)


def density(masas_amu, volumen_A3: float) -> float:
    """Densidad en kg/m^3 a partir de las masas y el volumen de celda."""
    return float(np.sum(masas_amu)) * AMU / (float(volumen_A3) * 1e-30)


def sound_velocities(B_GPa: float, G_GPa: float, rho: float) -> tuple:
    """(v_l, v_t, v_m) en m/s. B y G en GPa, rho en kg/m^3."""
    if G_GPa <= 0 or rho <= 0:
        return (None, None, None)
    B, G = B_GPa * 1e9, G_GPa * 1e9
    v_l = float(np.sqrt((B + 4.0 * G / 3.0) / rho))
    v_t = float(np.sqrt(G / rho))
    v_m = float((1.0 / 3.0 * (2.0 / v_t ** 3 + 1.0 / v_l ** 3)) ** (-1.0 / 3.0))
    return (v_l, v_t, v_m)


def debye_from_velocity(v_m: float, natoms: int, volumen_A3: float) -> float:
    """Temperatura de Debye elástica (límite acústico, baja temperatura)."""
    if not v_m:
        return None
    n = natoms / (volumen_A3 * 1e-30)          # átomos por m^3
    return float(HBAR / KB * (6.0 * np.pi ** 2 * n) ** (1.0 / 3.0) * v_m)


def debye_from_dos(omega_cm1, dos, natoms: int) -> float:
    """Temperatura de Debye por el segundo momento de la DOS de fonones.

    theta_D = (hbar/k_B) sqrt(5/3 <w^2>), con <w^2> pesado por la DOS
    normalizada a 3N. Usa TODO el espectro, no solo las acústicas: por eso
    no tiene por qué coincidir con la elástica.
    """
    from qekit.core.compat import trapezoid

    w = np.asarray(omega_cm1, dtype=float)
    g = np.asarray(dos, dtype=float)
    m = w > 1.0
    w, g = w[m], g[m]
    if w.size < 3:
        return None
    norm = trapezoid(g, w)
    if norm <= 0:
        return None
    w2 = trapezoid(g * w ** 2, w) / norm            # <w^2> en cm^-2
    # cm^-1 -> rad/s : w_rad = 2 pi c w_cm  (c en cm/s)
    w2_rad = w2 * (2.0 * np.pi * 2.99792458e10) ** 2
    return float(HBAR / KB * np.sqrt(5.0 / 3.0 * w2_rad))


def poisson_ratio(B_GPa: float, G_GPa: float) -> float:
    if G_GPa <= 0:
        return None
    return float((3.0 * B_GPa - 2.0 * G_GPa) / (2.0 * (3.0 * B_GPa + G_GPa)))


def gruneisen_from_poisson(nu: float) -> float:
    """Grüneisen a partir de la razón de Poisson (relación de Belomestnykh).

    Es una CORRELACIÓN empírica, no una derivación: sirve para tener el
    orden de magnitud cuando no hay fonones a varios volúmenes. El valor
    riguroso sale de la QHA.
    """
    if nu is None or nu >= 0.5:
        return None
    return float(3.0 * (1.0 + nu) / (2.0 * (2.0 - 3.0 * nu)))


def slack(theta_D: float, gamma: float, masa_media_amu: float,
          natoms: int, volumen_A3: float, T: float = 300.0) -> float:
    """Conductividad térmica de red por el modelo de Slack (W/m/K).

        kappa = A * M_avg * theta_D^3 * delta / (gamma^2 * n^(2/3) * T)

    con delta el volumen por átomo elevado a 1/3. Es una ESTIMACIÓN de
    orden de magnitud: el prefactor A es empírico y el modelo supone
    dispersión de tres fonones dominante y cristal simple.
    """
    if not theta_D or not gamma:
        return None
    # delta va en ANGSTROM: el prefactor empirico A = 3.1e-6 esta ajustado
    # para M en amu, theta en K y delta en A. Pasarlo a metros mata el
    # resultado por diez ordenes de magnitud.
    delta = (volumen_A3 / natoms) ** (1.0 / 3.0)                # A
    A = 3.1e-6 / (1.0 - 0.514 / gamma + 0.228 / gamma ** 2)
    return float(A * masa_media_amu * theta_D ** 3 * delta
                 / (gamma ** 2 * natoms ** (2.0 / 3.0) * T))


def analyze(B_GPa: float, G_GPa: float, masas_amu, volumen_A3: float,
            natoms: int = None, T: float = 300.0) -> Termoelastico:
    natoms = natoms or len(masas_amu)
    r = Termoelastico(natoms=natoms, volumen=volumen_A3, T=float(T))
    r.rho = density(masas_amu, volumen_A3)
    r.v_l, r.v_t, r.v_m = sound_velocities(B_GPa, G_GPa, r.rho)
    r.theta_D = debye_from_velocity(r.v_m, natoms, volumen_A3)
    r.poisson = poisson_ratio(B_GPa, G_GPa)
    r.gruneisen = gruneisen_from_poisson(r.poisson)
    r.kappa_slack = slack(r.theta_D, r.gruneisen,
                          float(np.mean(masas_amu)), natoms, volumen_A3, T)
    if r.poisson is not None and r.poisson < 0:
        r.avisos.append(
            f"razón de Poisson negativa ({r.poisson:.3f}): es posible "
            "(materiales auxéticos) pero raro;\nrevisa las Cij antes de "
            "creerlo.")
    return r


def is_cubic_tensor(C: np.ndarray, rtol: float = 0.05,
                    atol_GPa: float = 2.0) -> bool:
    """¿Tiene el tensor Cij la forma de un cristal cúbico?

    Tres constantes independientes: C11 = C22 = C33, C12 = C13 = C23,
    C44 = C55 = C66, y ceros fuera de esos bloques. Se admite un `rtol`
    relativo (o `atol_GPa` absoluto, lo que sea mayor) porque las Cij de un
    barrido de deformaciones llevan ruido numérico. Sirve para no imprimir
    v_L[100] = sqrt(C11/rho) como si fuera cúbico un tensor que no lo es.
    """
    C = np.asarray(C, dtype=float)
    if C.shape != (6, 6) or not np.all(np.isfinite(C)):
        return False
    S = 0.5 * (C + C.T)

    def iguales(vals):
        vals = np.asarray(vals, dtype=float)
        ref = float(np.mean(np.abs(vals)))
        return float(np.ptp(vals)) <= max(rtol * ref, atol_GPa)

    diag = [S[0, 0], S[1, 1], S[2, 2]]
    offd = [S[0, 1], S[0, 2], S[1, 2]]
    shear = [S[3, 3], S[4, 4], S[5, 5]]
    if not (iguales(diag) and iguales(offd) and iguales(shear)):
        return False
    # acoplamientos normal-cortante y cortante-cortante cruzados: nulos
    escala = max(abs(S[0, 0]), abs(S[3, 3]), 1.0)
    resto = np.array(S)
    resto[:3, :3] = 0.0
    resto[3, 3] = resto[4, 4] = resto[5, 5] = 0.0
    return bool(np.max(np.abs(resto)) <= max(rtol * escala, atol_GPa))


def cubic_directional(C: np.ndarray, rho: float) -> dict:
    """Velocidades del sonido a lo largo de [100] en un cristal cubico.

    Son las que hay que comparar contra la pendiente de las ramas
    acusticas del camino Gamma->X, que va justo en esa direccion. Los
    promedios de Voigt-Reuss-Hill son isotropos (policristal) y NO son
    lo mismo: compararlos contra una rama concreta mezcla dos cosas.

        v_L[100] = sqrt(C11/rho)      v_T[100] = sqrt(C44/rho)
    """
    if rho <= 0:
        return {}
    c11 = float(C[0, 0]) * 1e9
    c44 = float(C[3, 3]) * 1e9
    if c11 <= 0 or c44 <= 0:
        return {}
    return {"v_l_100": float(np.sqrt(c11 / rho)),
            "v_t_100": float(np.sqrt(c44 / rho))}


def acoustic_velocities(qdist, freqs, n_puntos: int = 4) -> dict:
    """Velocidades del sonido desde la PENDIENTE de las ramas acústicas.

    Ajusta una recta por el origen a las tres ramas más bajas cerca de Γ.
    Es la ruta independiente con la que se cruza el módulo elástico.

    AVISO: la pendiente en q->0 es justo lo que peor interpola una malla
    de q gruesa. Con 2x2x2 las ramas transversales de un semiconductor
    salen mal, y el cruce lo delata.
    """
    q = np.asarray(qdist, dtype=float)
    f = np.asarray(freqs, dtype=float)
    n = min(int(n_puntos), len(q) - 1)
    if n < 2:
        return {}
    qq = q[1:n + 1]
    out = {}
    vs = []
    for b in range(min(3, f.shape[1])):
        ff = f[1:n + 1, b]
        # recta forzada por el origen: v = sum(q*f)/sum(q^2)
        pend = float(np.sum(qq * ff) / np.sum(qq ** 2))
        vs.append(pend * CM1_A_A_MS)
    vs = sorted(vs)
    out["v_ramas"] = vs
    out["v_t1"], out["v_t2"], out["v_l"] = vs[0], vs[1], vs[2]
    out["n_puntos"] = n
    out["qmax_A-1"] = float(qq[-1])
    return out


def export(r: Termoelastico, outdir: str = ".") -> list:
    """Deja los números en un .dat, no solo en la pantalla.

    Hasta la 0.26 este comando solo imprimía. Un número que solo existe en
    la terminal no lo puede leer `crosscheck` para cruzarlo, ni `datasheet`
    para meterlo en la ficha, ni tú dentro de tres semanas. Lo detectó
    escribir las recetas: la receta decía «deja DERIVED.dat» y no era
    verdad.
    """
    from pathlib import Path as _P
    from qekit.core import provenance
    out = _P(outdir); out.mkdir(parents=True, exist_ok=True)
    f = out / "DERIVED.dat"
    filas = [("densidad", r.rho, "kg/m3"),
             ("v_longitudinal", r.v_l, "m/s"),
             ("v_transversal", r.v_t, "m/s"),
             ("v_media", r.v_m, "m/s"),
             ("theta_Debye_elastica", r.theta_D, "K"),
             ("razon_de_Poisson", r.poisson, ""),
             ("gruneisen", r.gruneisen, ""),
             (f"kappa_Slack_{r.T:g}K", r.kappa_slack, "W/m/K")]
    L = [provenance.header("derivadas termoelásticas",
                           {"atomos": r.natoms,
                            "volumen_A3": f"{r.volumen:.4f}"
                            if r.volumen else "?"}),
         f"# {'magnitud':<24s} {'valor':>16s}  unidad"]
    for nom, val, uni in filas:
        if val is None:
            continue
        L.append(f"  {nom:<24s} {val:16.6f}  {uni}")
    f.write_text("\n".join(L) + "\n", encoding="utf-8")
    return [str(f)]


def report(r: Termoelastico) -> str:
    lines = ["--- Propiedades termoelásticas ---",
             f"Densidad: {r.rho:.1f} kg/m³  "
             f"({r.natoms} átomos en {r.volumen:.2f} Å³)"]
    if r.v_l:
        lines += [f"Velocidades del sonido:  longitudinal {r.v_l:.0f} m/s  |  "
                  f"transversal {r.v_t:.0f} m/s",
                  f"                         promedio     {r.v_m:.0f} m/s"]
    if r.theta_D:
        lines.append(f"Temperatura de Debye (elástica): {r.theta_D:.0f} K")
    if r.poisson is not None:
        lines.append(f"Razón de Poisson: {r.poisson:.4f}")
    if r.gruneisen:
        lines.append(f"Grüneisen (correlación de Poisson): {r.gruneisen:.2f}")
    if r.kappa_slack:
        lines.append(f"Conductividad térmica de red (Slack, {r.T:g} K): "
                     f"{r.kappa_slack:.1f} W/(m·K)")
    for a in r.avisos:
        lines.append(f"\nAVISO: {a}")
    lines += ["",
              "La Debye de aquí es la ELÁSTICA: sale de las velocidades del "
              "sonido y describe\nel límite de baja temperatura. La que "
              "sale de la DOS de fonones usa todo el\nespectro y da otro "
              "número; no son la misma cantidad.",
              "El Grüneisen viene de una correlación empírica con la razón "
              "de Poisson y la\nconductividad de Slack es una estimación de "
              "orden de magnitud, no un valor\npara reportar sin más."]
    return "\n".join(lines)
