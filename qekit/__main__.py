# Olla-DFT — command-line toolkit for Quantum ESPRESSO
# Copyright (C) 2026 Jorge Enrique González Sevilla
# SPDX-License-Identifier: AGPL-3.0-or-later
# This program is free software: you can redistribute it and/or modify it
# under the terms of the GNU Affero General Public License as published by the Free
# Software Foundation, either version 3 of the License, or (at your option)
# any later version. See the LICENSE file for details.

import sys

from qekit.cli import main

if __name__ == "__main__":
    # sys.exit y no una llamada suelta: sin esto `python -m qekit` devuelve
    # siempre 0 y un error de uso pasa por exito en cualquier script que
    # encadene comandos con && o revise $?.
    sys.exit(main())
