# Olla-DFT — command-line toolkit for Quantum ESPRESSO
# Copyright (C) 2026 Jorge Enrique González Sevilla
# SPDX-License-Identifier: GPL-3.0-or-later
# This program is free software: you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by the Free
# Software Foundation, either version 3 of the License, or (at your option)
# any later version. See the LICENSE file for details.

"""Errores de Olla-DFT, separados por a quién le toca arreglarlos.

La distinción importa para el registro de incidencias. Si cada bandera mal
escrita queda archivada como "falla del programa", el registro se llena de
ruido y deja de servir para lo que sirve: ver qué se rompe de verdad.

- `ErrorDeUso`  -> el comando está mal escrito, o el dato no encaja con lo
  que el comando necesita. El programa hizo lo correcto: avisar. Se anota
  como tipo "uso" —sin traza y sin alarma— porque la estadística de qué
  mensajes de uso se repiten SÍ vale: es justo donde la interfaz confunde.
- cualquier otra excepción -> falla del programa. Esa sí se archiva
  completa, con traza y versiones.

Un `ErrorDeUso` sale con código 2, igual que argparse cuando rechaza una
bandera; una falla del programa sale con 1.
"""


class ErrorDeUso(ValueError):
    """El comando o sus datos no encajan; el mensaje ya explica qué hacer.

    Es una subclase de ValueError a propósito: el código que ya atrapaba
    ValueError —incluidas las pruebas— sigue funcionando igual.
    """


class FaltanDatos(ErrorDeUso):
    """Falta un resultado previo: hay que correr otro paso antes."""
