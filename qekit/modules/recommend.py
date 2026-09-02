# Olla-DFT — command-line toolkit for Quantum ESPRESSO
# Copyright (C) 2026 Jorge Enrique González Sevilla
# SPDX-License-Identifier: GPL-3.0-or-later
# This program is free software: you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by the Free
# Software Foundation, either version 3 of the License, or (at your option)
# any later version. See the LICENSE file for details.

"""Recomendaciones a partir de TU propio historial de cálculos.

La base de `olla-dft db` acumula, sin que cueste nada extra, lo que hace falta
para no repetir el mismo tanteo: qué cutoffs convergieron, qué mezcla
funcionó, cuánto tardó cada cosa y en qué sistemas.

POR QUÉ ESTO NO ES UNA RED NEURONAL
-----------------------------------
Con unas decenas de cálculos —que es lo que va a haber— un modelo aprendido
sobreajusta y no se puede auditar. Aquí se usa lo que sí funciona con pocos
datos: buscar los cálculos PARECIDOS al que quieres hacer (mismos
elementos, sistema parecido) y mirar qué les funcionó, diciendo siempre
cuántos casos respaldan cada número.

Una recomendación con un solo caso detrás se marca como tal. Es la
diferencia entre "esto suele funcionar" y "esto funcionó una vez".

Y una cosa que NO hace: inventar cutoffs. Si no hay historial del elemento,
lo dice y remite a los cutoffs que declara el propio pseudopotencial, que
es un dato y no una predicción.
"""

from dataclasses import dataclass

import numpy as np


@dataclass
class Sugerencia:
    campo: str = ""
    valor: object = None
    n_casos: int = 0
    rango: tuple = None
    razon: str = ""
    confianza: str = "baja"      # "baja" | "media" | "alta"


def _confianza(n: int) -> str:
    if n >= 8:
        return "alta"
    if n >= 3:
        return "media"
    return "baja"


def similares(filas: list, elementos, natoms: int = None,
              tol_natoms: float = 2.0) -> list:
    """Cálculos del historial parecidos al que se quiere hacer.

    "Parecido" = comparte al menos un elemento; se prioriza compartir
    todos. No se usa una métrica sofisticada a propósito: con pocos datos,
    cualquier cosa más elaborada da una falsa sensación de precisión.
    """
    els = set(elementos)
    out = []
    for f in filas:
        if not f.get("convergido"):
            continue
        formula = f.get("formula") or ""
        f_els = set(_elementos_de(formula))
        if not (f_els & els):
            continue
        puntaje = len(f_els & els) / max(len(f_els | els), 1)
        if natoms and f.get("natoms"):
            razon = f["natoms"] / natoms
            if razon > tol_natoms or razon < 1.0 / tol_natoms:
                puntaje *= 0.5
        out.append((puntaje, f))
    return [f for _p, f in sorted(out, key=lambda t: -t[0])]


def _elementos_de(formula: str) -> list:
    import re
    return re.findall(r"[A-Z][a-z]?", formula or "")


def sugerir(filas: list, elementos, natoms: int = None,
            es_losa: bool = False) -> list:
    """Sugerencias de parámetros para un cálculo nuevo."""
    sug = []
    vecinos = similares(filas, elementos, natoms)
    if not vecinos:
        sug.append(Sugerencia(
            campo="(sin historial)", n_casos=0, confianza="baja",
            razon="No hay cálculos previos con estos elementos. Usa los "
                  "cutoffs que declara\nel propio pseudopotencial (Olla-DFT "
                  "los lee del UPF) o la tabla SSSP: eso es un\ndato "
                  "medido, no una predicción, y siempre le gana a una "
                  "extrapolación."))
        return sug

    ecuts = [f["ecutwfc"] for f in vecinos if f.get("ecutwfc")]
    if ecuts:
        v = float(np.max(ecuts))
        sug.append(Sugerencia(
            campo="ecutwfc", valor=v, n_casos=len(ecuts),
            rango=(float(np.min(ecuts)), float(np.max(ecuts))),
            confianza=_confianza(len(ecuts)),
            razon=f"el MÁXIMO de {len(ecuts)} cálculos convergidos con "
                  f"estos elementos (rango {min(ecuts):.0f}–{max(ecuts):.0f} "
                  "Ry). Se toma el máximo, no la media: un cutoff bajo que "
                  "funcionó en un\nsistema no garantiza nada en otro."))

    duales = [(f["ecutrho"] / f["ecutwfc"]) for f in vecinos
              if f.get("ecutrho") and f.get("ecutwfc")]
    if duales:
        sug.append(Sugerencia(
            campo="dual (ecutrho/ecutwfc)", valor=float(np.max(duales)),
            n_casos=len(duales), confianza=_confianza(len(duales)),
            razon="el dual que usaron los cálculos previos; depende del "
                  "tipo de pseudo\n(4 para norma conservada, 8-12 para "
                  "ultrasoft y PAW)."))

    dens = [f["kdensity"] for f in vecinos if f.get("kdensity")]
    if dens:
        sug.append(Sugerencia(
            campo="densidad de k (puntos/Å⁻³)", valor=float(np.median(dens)),
            n_casos=len(dens), rango=(float(np.min(dens)),
                                      float(np.max(dens))),
            confianza=_confianza(len(dens)),
            razon="mediana de los cálculos convergidos. La densidad, no el "
                  "número de puntos,\nes lo comparable entre celdas de "
                  "tamaño distinto."))

    pasos = [f["n_scf"] for f in vecinos if f.get("n_scf")]
    if pasos and np.median(pasos) > 40:
        sug.append(Sugerencia(
            campo="electron_maxstep", valor=300, n_casos=len(pasos),
            confianza=_confianza(len(pasos)),
            razon=f"tus cálculos con estos elementos necesitaron una "
                  f"mediana de {np.median(pasos):.0f} pasos SCF:\nel "
                  "máximo por defecto se queda corto."))

    if es_losa:
        sug.append(Sugerencia(
            campo="mixing_beta", valor=0.3, n_casos=0, confianza="baja",
            razon="es una losa con vacío: son las que más oscilación de "
                  "carga dan. Empezar\ncon mixing_beta bajo y "
                  "mixing_mode='local-TF' ahorra reintentos. Esto no sale "
                  "de\ntu historial, es una regla general."))
    return sug


def report(sug: list, elementos, n_historial: int = 0) -> str:
    lines = ["--- Sugerencias desde tu historial ---",
             f"Elementos: {', '.join(elementos)}  |  "
             f"cálculos en la base: {n_historial}", ""]
    if not sug or sug[0].campo == "(sin historial)":
        lines.append(sug[0].razon if sug else "Sin datos.")
        return "\n".join(lines)

    for s in sug:
        val = s.valor
        if isinstance(val, float):
            val = f"{val:.4g}"
        # la confianza 'baja' cubre 1 y 2 casos: no se dice "un solo caso"
        # cuando hay dos
        marca = {"alta": "", "media": "  (pocos casos)",
                 "baja": ("  (UN SOLO CASO: tómalo como indicio)"
                          if s.n_casos == 1 else
                          f"  (SOLO {s.n_casos} CASOS: tómalo como indicio)")
                 }[s.confianza]
        if s.n_casos == 0:
            marca = "  (regla general, no de tu historial)"
        lines.append(f"  {s.campo}: {val}"
                     f"   [{s.n_casos} caso{'s' if s.n_casos != 1 else ''}]"
                     f"{marca}")
        for l in s.razon.splitlines():
            lines.append(f"      {l}")
        lines.append("")
    lines.append("Estas sugerencias salen de lo que YA te funcionó, no de un "
                 "modelo entrenado.\nNo sustituyen a una prueba de "
                 "convergencia: 'olla-dft converge' sigue siendo\nla forma de "
                 "saberlo de verdad para un sistema nuevo.")
    return "\n".join(lines)
