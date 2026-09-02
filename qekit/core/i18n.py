# Olla-DFT — command-line toolkit for Quantum ESPRESSO
# Copyright (C) 2026 Jorge Enrique González Sevilla
# SPDX-License-Identifier: GPL-3.0-or-later
# This program is free software: you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by the Free
# Software Foundation, either version 3 of the License, or (at your option)
# any later version. See the LICENSE file for details.

"""Idioma de la interfaz: español (por defecto) o inglés.

Olla-DFT se escribió en español y ese sigue siendo el idioma de los
informes científicos. Lo que sí se traduce es la *interfaz*: la ayuda de
cada comando y de cada bandera, el menú interactivo, el inicio guiado, el
dashboard y la referencia HTML.

El idioma se decide, en este orden:

  1. la bandera global ``--language en`` (se acepta en cualquier posición);
  2. la variable de entorno ``OLLA_DFT_LANG``;
  3. la clave ``language`` de la configuración (``olla-dft config set language en``);
  4. español.

Las traducciones viven en ``qekit/data/i18n/cli_en.json`` (ayuda de la CLI)
y en los demás archivos de esa carpeta. Este módulo no traduce mensajes en
tiempo de ejecución: solo textos de ayuda e interfaz.

Los catálogos de datos (las recetas, las metas del asistente) no se
duplican por idioma: ``translate_data`` recorre la estructura original y
devuelve una copia con cada cadena pasada por una tabla ``{es: en}`` que
vive en ``qekit/data/i18n/<nombre>.json`` (``load_table``).
"""

import dataclasses
import json
import os
from functools import lru_cache
from pathlib import Path

LANGUAGES = ("es", "en")
DEFAULT_LANGUAGE = "es"
ENV_VAR = "OLLA_DFT_LANG"

_I18N_DIR = Path(__file__).resolve().parent.parent / "data" / "i18n"
_current = None


def set_language(language) -> str:
    """Fija el idioma del proceso; devuelve el que quedó activo."""
    global _current
    if language is None:
        _current = None
        return get_language()
    if language not in LANGUAGES:
        raise ValueError(f"idioma no admitido: {language!r} (usa es o en)")
    _current = language
    return _current


def get_language() -> str:
    """El idioma activo, siguiendo el orden de prioridad del módulo."""
    if _current in LANGUAGES:
        return _current
    env = os.environ.get(ENV_VAR, "").strip().lower()
    if env in LANGUAGES:
        return env
    try:
        from qekit import config as qcfg
        cfg = qcfg.load().get("language", "").strip().lower()
        if cfg in LANGUAGES:
            return cfg
    except Exception:                                  # noqa: BLE001
        pass
    return DEFAULT_LANGUAGE


def extract_language(argv):
    """Quita ``--language X`` / ``--language=X`` de argv, esté donde esté.

    Devuelve (argv_limpio, idioma_o_None). Argparse solo aceptaría la bandera
    delante del subcomando; la gente la escribe donde le sale.
    """
    limpio, idioma = [], None
    i = 0
    while i < len(argv):
        tok = argv[i]
        if tok == "--language" and i + 1 < len(argv):
            idioma = argv[i + 1]
            i += 2
            continue
        if tok.startswith("--language="):
            idioma = tok.split("=", 1)[1]
            i += 1
            continue
        limpio.append(tok)
        i += 1
    if idioma is not None and idioma not in LANGUAGES:
        raise ValueError(idioma)
    return limpio, idioma


@lru_cache(maxsize=None)
def _load(name: str) -> dict:
    path = _I18N_DIR / name
    if not path.exists():
        return {}
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError):
        return {}


def default_help(flag: str) -> str:
    """Ayuda en español para banderas que se repiten en muchos comandos."""
    return _load("cli_es.json").get("defaults", {}).get(flag, "")


def translate(text: str, language=None) -> str:
    """Traduce un texto de ayuda al idioma activo; si no hay, lo deja igual."""
    language = language or get_language()
    if language == "es" or not text:
        return text
    table = _load("cli_en.json")
    return table.get("help", {}).get(text) or text


def ui(key: str, language=None) -> str:
    """Rótulos fijos de argparse (títulos de secciones, 'mostrar esta ayuda')."""
    language = language or get_language()
    es = _load("cli_es.json").get("ui", {})
    if language == "es":
        return es.get(key, key)
    return _load("cli_en.json").get("ui", {}).get(key, es.get(key, key))


def load_table(name: str) -> dict:
    """Carga ``qekit/data/i18n/<name>.json`` (con caché).

    Devuelve el diccionario entero del archivo. Por convención lleva una
    clave ``"strings"`` con el mapa ``{texto_es: texto_en}`` que usa
    ``translate_data`` y, si hace falta para la búsqueda, ``"keywords"``.
    Si el archivo no existe, devuelve ``{}``: la interfaz sigue en español.
    """
    return _load(f"{name}.json")


def translate_data(obj, table):
    """Copia profunda de ``obj`` con cada ``str`` pasado por ``table``.

    Recorre dataclasses, listas, tuplas (siguen siendo tuplas) y
    diccionarios (claves y valores). Lo que no es ``str`` se deja tal cual,
    y una cadena que no esté en la tabla se copia sin cambios: por eso los
    comandos, los nombres de archivo y las salidas reales de los programas
    no necesitan marcarse, basta con no incluirlos en la tabla. Aplicarlo
    dos veces no cambia nada mientras ningún texto traducido sea a su vez
    una clave de la tabla con otro valor.
    """
    if isinstance(obj, str):
        return table.get(obj, obj)
    if isinstance(obj, list):
        return [translate_data(x, table) for x in obj]
    if isinstance(obj, tuple):
        return tuple(translate_data(x, table) for x in obj)
    if isinstance(obj, dict):
        return {translate_data(k, table): translate_data(v, table)
                for k, v in obj.items()}
    if dataclasses.is_dataclass(obj) and not isinstance(obj, type):
        cambios = {f.name: translate_data(getattr(obj, f.name), table)
                   for f in dataclasses.fields(obj)}
        return dataclasses.replace(obj, **cambios)
    return obj
