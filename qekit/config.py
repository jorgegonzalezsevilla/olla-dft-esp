# Olla-DFT — command-line toolkit for Quantum ESPRESSO
# Copyright (C) 2026 Jorge Enrique González Sevilla
# SPDX-License-Identifier: GPL-3.0-or-later
# This program is free software: you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by the Free
# Software Foundation, either version 3 of the License, or (at your option)
# any later version. See the LICENSE file for details.

"""Configuración persistente de Olla-DFT.

Se guarda en la carpeta de configuración de cada sistema (en Linux,
~/.config/olla-dft/config.ini) y define valores que el usuario no quiere
repetir en cada comando (ruta de pseudopotenciales, cutoffs por defecto,
idioma, etc.).
"""

import configparser
import shutil
from pathlib import Path

from qekit.core import i18n, plataforma

CONFIG_DIR = plataforma.dir_config()
CONFIG_FILE = CONFIG_DIR / "config.ini"


def _migrar_si_hace_falta():
    """Trae la configuración de una versión anterior si la hay.

    Hasta la 0.27 se guardaba en ~/.config/qekit; de la 0.28 a la 0.35 en la
    carpeta ``QEkit`` de cada sistema. Al cambiar de carpeta, la
    configuración de quien ya usaba el programa se quedaría atrás sin decir
    nada: sus pseudos y sus cutoffs dejarían de encontrarse y parecería que
    el programa se ha olvidado de todo. Se copia una vez, sin borrar el
    original.
    """
    if CONFIG_FILE.exists():
        return False
    for vieja_dir in plataforma.dirs_config_heredados():
        vieja = vieja_dir / "config.ini"
        if not vieja.exists() or vieja_dir == CONFIG_DIR:
            continue
        try:
            CONFIG_DIR.mkdir(parents=True, exist_ok=True)
            shutil.copy2(vieja, CONFIG_FILE)
            plantillas = vieja_dir / "templates"
            if plantillas.is_dir() and not (CONFIG_DIR / "templates").exists():
                shutil.copytree(plantillas, CONFIG_DIR / "templates")
            return True
        except OSError:
            return False
    return False

# Valores por defecto (se usan si el usuario no ha configurado nada).
DEFAULTS = {
    "pseudo_dir": str(Path.home() / "pseudos"),
    "ecutwfc": "60.0",       # Ry — revisa la tabla SSSP para tu sistema
    "dual": "8",             # ecutrho = dual * ecutwfc (8 para USPP/PAW, 4 para NC)
    "kspacing": "0.20",      # Å^-1 (incluye el factor 2*pi) — malla scf
    "kspacing_nscf": "0.12", # Å^-1 — malla densa para nscf/DOS
    "band_points": "20",     # puntos por segmento del k-path
    "degauss": "0.01",       # Ry — smearing
    "smearing": "cold",
    "nproc": "4",            # procesos MPI para el script run.sh y para --run
    # El ejecutable y el lanzador dependen del sistema: en Windows es pw.exe
    # y no hay mpirun (hay mpiexec, de MS-MPI). Se detectan al vuelo en vez
    # de dejar un valor fijo que solo acierta en Linux.
    "pw_cmd": plataforma.nombres_ejecutable("pw")[0],
    "mpi_cmd": plataforma.lanzador_mpi() or "",
    "language": i18n.DEFAULT_LANGUAGE,   # idioma de la interfaz: es o en
}

VALID_KEYS = sorted(DEFAULTS.keys())


def load() -> dict:
    """Lee la configuración; devuelve DEFAULTS combinado con lo guardado."""
    _migrar_si_hace_falta()
    values = dict(DEFAULTS)
    if CONFIG_FILE.exists():
        parser = configparser.ConfigParser()
        parser.read(CONFIG_FILE)
        if parser.has_section("qekit"):
            for key, val in parser.items("qekit"):
                values[key] = val
    return values


def save(values: dict) -> None:
    """Guarda solo las claves válidas en el archivo de configuración."""
    parser = configparser.ConfigParser()
    parser["qekit"] = {k: str(v) for k, v in values.items() if k in DEFAULTS}
    CONFIG_DIR.mkdir(parents=True, exist_ok=True)
    with open(CONFIG_FILE, "w") as fh:
        parser.write(fh)


def set_value(key: str, value: str) -> None:
    if key not in DEFAULTS:
        raise KeyError(
            f"clave desconocida '{key}'. Claves válidas: {', '.join(VALID_KEYS)}"
        )
    if key == "language" and value not in ("es", "en"):
        raise KeyError("language admite 'es' o 'en'")
    values = load()
    values[key] = value
    save(values)


def show() -> str:
    values = load()
    lines = [f"Archivo de configuración: {CONFIG_FILE}"]
    if not CONFIG_FILE.exists():
        lines.append("(aún no existe — se muestran los valores por defecto)")
    for key in VALID_KEYS:
        lines.append(f"  {key:14s} = {values[key]}")
    return "\n".join(lines)
