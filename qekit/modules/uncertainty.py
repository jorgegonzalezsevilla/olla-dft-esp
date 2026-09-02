# Olla-DFT — command-line toolkit for Quantum ESPRESSO
# Copyright (C) 2026 Jorge Enrique González Sevilla
# SPDX-License-Identifier: GPL-3.0-or-later
# This program is free software: you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by the Free
# Software Foundation, either version 3 of the License, or (at your option)
# any later version. See the LICENSE file for details.

"""Utilidades pequeñas para propagar incertidumbres experimentales o numéricas.

Son deliberadamente transparentes: usa derivadas centradas y asume entradas
independientes. Olla-DFT nunca inventa una incertidumbre si el cálculo no la
reporta; estas funciones solo hacen explícita una que el usuario proporciona.
"""

from __future__ import annotations

import math

from qekit.core.errors import ErrorDeUso


def propagate(function, values, uncertainties, relative_step=1e-6) -> dict:
    """Propaga sigma mediante cuadratura de derivadas parciales centradas."""
    values = [float(value) for value in values]
    sigmas = [float(value) for value in uncertainties]
    if len(values) != len(sigmas) or not values:
        raise ErrorDeUso("values y uncertainties deben tener la misma longitud no vacía.")
    if any(not math.isfinite(x) or x < 0 for x in sigmas):
        raise ErrorDeUso("las incertidumbres deben ser finitas y no negativas.")
    try:
        central = float(function(values))
    except Exception as exc:  # noqa: BLE001
        raise ErrorDeUso(f"no se pudo evaluar la función: {exc}") from None
    if not math.isfinite(central):
        raise ErrorDeUso("la función produjo un valor no finito.")
    derivatives = []
    variance = 0.0
    for index, (value, sigma) in enumerate(zip(values, sigmas)):
        step = max(abs(value) * float(relative_step), float(relative_step))
        plus, minus = values[:], values[:]
        plus[index] += step
        minus[index] -= step
        try:
            derivative = (float(function(plus)) - float(function(minus))) / (2 * step)
        except Exception as exc:  # noqa: BLE001
            raise ErrorDeUso(f"no se pudo derivar la entrada {index}: {exc}") from None
        if not math.isfinite(derivative):
            raise ErrorDeUso(f"la derivada de la entrada {index} no es finita.")
        derivatives.append(derivative)
        variance += (derivative * sigma) ** 2
    return {"value": central, "uncertainty": math.sqrt(variance),
            "derivatives": derivatives, "assumption": "entradas independientes"}


def weighted_mean(values, uncertainties) -> dict:
    """Media ponderada con incertidumbres sigma no nulas."""
    values = [float(value) for value in values]
    sigmas = [float(value) for value in uncertainties]
    if len(values) != len(sigmas) or not values:
        raise ErrorDeUso("values y uncertainties deben tener la misma longitud no vacía.")
    if any(not math.isfinite(x) or x <= 0 for x in sigmas):
        raise ErrorDeUso("las incertidumbres deben ser finitas y positivas.")
    weights = [1.0 / sigma ** 2 for sigma in sigmas]
    total = sum(weights)
    return {"value": sum(w * value for w, value in zip(weights, values)) / total,
            "uncertainty": math.sqrt(1.0 / total), "weights": weights}
