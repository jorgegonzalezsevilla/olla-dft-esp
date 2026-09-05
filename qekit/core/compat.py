# Olla-DFT — command-line toolkit for Quantum ESPRESSO
# Copyright (C) 2026 Jorge Enrique González Sevilla
# SPDX-License-Identifier: AGPL-3.0-or-later
# This program is free software: you can redistribute it and/or modify it
# under the terms of the GNU Affero General Public License as published by the Free
# Software Foundation, either version 3 of the License, or (at your option)
# any later version. See the LICENSE file for details.

"""Compatibilidad entre versiones de las dependencias.

Aquí viven los parches de nombres que cambian entre versiones, para que el
resto del código no tenga que preguntarse con qué numpy está corriendo.
"""

import numpy as np

# numpy 2.0 renombró trapz -> trapezoid y ELIMINÓ el nombre viejo; el nombre
# nuevo no existe en numpy 1.x. Sin esto, Olla-DFT funciona con la versión que
# tenga quien lo escribió y falla con AttributeError en la otra — justo el
# tipo de error que no aparece en la máquina de desarrollo.
trapezoid = getattr(np, "trapezoid", None)
if trapezoid is None:                       # numpy < 2.0
    trapezoid = np.trapz                    # noqa: NPY201

__all__ = ["trapezoid"]
