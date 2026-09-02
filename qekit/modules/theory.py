# Olla-DFT — command-line toolkit for Quantum ESPRESSO
# Copyright (C) 2026 Jorge Enrique González Sevilla
# SPDX-License-Identifier: GPL-3.0-or-later
# This program is free software: you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by the Free
# Software Foundation, either version 3 of the License, or (at your option)
# any later version. See the LICENSE file for details.

"""Fundamento físico de cada comando, consultable desde la terminal.

La misma documentación que se publica en ``docs/TEORIA.md`` y
``docs/THEORY.md`` viaja dentro del paquete (``qekit/data/theory``) para
poder leerla sin conexión: ``olla-dft teoria eos`` imprime qué responde el
comando, el fundamento explicado sin jerga, las fórmulas que el código
implementa, de qué módulo salen y de dónde sale cada dato.

Los archivos están escritos a mano, en Markdown, uno por área y por idioma.
Cada sección empieza por ``### `olla-dft <comando>` — título``; ese
encabezado es el índice.
"""

import re
from dataclasses import dataclass
from pathlib import Path

from qekit import __command_name__
from qekit.core.errors import ErrorDeUso

THEORY_DIR = Path(__file__).resolve().parent.parent / "data" / "theory"

#: Orden de las áreas en el documento completo.
AREAS = ("electronica", "mecanica", "espectros")

_HEADING = re.compile(r"^### `" + re.escape(__command_name__) + r" ([\w-]+)`"
                      r"(?: (?:y|and) `" + re.escape(__command_name__) + r" ([\w-]+)`)?"
                      r"\s+[—-]\s+(.*)$")


@dataclass
class Seccion:
    comandos: tuple      # ("eos",) o ("audit", "db")
    titulo: str
    area: str
    texto: str           # el Markdown completo de la sección, con su encabezado


def _archivo(area: str, language: str) -> Path:
    if language not in ("es", "en"):
        raise ErrorDeUso("language debe ser es o en")
    return THEORY_DIR / f"{area}.{language}.md"


def area_titulo(area: str, language: str = "es") -> str:
    """El encabezado de nivel 2 del archivo de un área."""
    for line in _archivo(area, language).read_text(encoding="utf-8").splitlines():
        if line.startswith("## "):
            return line[3:].strip()
    return area


def secciones(language: str = "es", areas=AREAS) -> list:
    """Todas las secciones, en el orden del documento."""
    fuera = []
    for area in areas:
        path = _archivo(area, language)
        if not path.exists():
            continue
        actual, cabecera = None, None
        for line in path.read_text(encoding="utf-8").splitlines(keepends=True):
            m = _HEADING.match(line.rstrip("\n"))
            if m:
                if actual is not None:
                    fuera.append(Seccion(cabecera[0], cabecera[1], area,
                                         "".join(actual).rstrip() + "\n"))
                cmds = tuple(c for c in (m.group(1), m.group(2)) if c)
                cabecera = (cmds, m.group(3).strip())
                actual = [line]
            elif actual is not None:
                actual.append(line)
        if actual is not None:
            fuera.append(Seccion(cabecera[0], cabecera[1], area,
                                 "".join(actual).rstrip() + "\n"))
    return fuera


def buscar(comando: str, language: str = "es"):
    """La sección de un comando, o None si no está documentado."""
    comando = comando.strip().lower()
    for sec in secciones(language):
        if comando in sec.comandos:
            return sec
    return None


def indice(language: str = "es") -> str:
    """Lista legible de todo lo documentado, agrupado por área."""
    intro = {"es": f"Fundamento físico por comando. Uso:  {__command_name__} "
                   "teoria COMANDO   (o --all para todo, -o para guardarlo)",
             "en": f"Physics behind each command. Usage:  {__command_name__} "
                   "theory COMMAND   (or --all for everything, -o to save it)"}
    lines = [intro[language], ""]
    for area in AREAS:
        lines.append(area_titulo(area, language))
        for sec in secciones(language, (area,)):
            nombre = " / ".join(sec.comandos)
            lines.append(f"  {nombre:14s} {sec.titulo}")
        lines.append("")
    return "\n".join(lines).rstrip()


def _quitar_markdown(texto: str) -> str:
    """Deja el Markdown legible en una terminal sin renderizarlo del todo."""
    fuera = []
    for line in texto.splitlines():
        if line.startswith("### "):
            line = line[4:].replace("`", "")
            fuera += [line, "=" * min(len(line), 78)]
            continue
        line = re.sub(r"\*\*(.+?)\*\*", r"\1", line)
        fuera.append(line)
    return "\n".join(fuera)


def texto(comando: str = None, language: str = "es", crudo: bool = False) -> str:
    """El texto que imprime el comando: una sección, o el índice."""
    if not comando:
        return indice(language)
    sec = buscar(comando, language)
    if sec is None:
        conocidos = sorted({c for s in secciones(language) for c in s.comandos})
        raise ErrorDeUso(
            f"no hay fundamento escrito para '{comando}'. Los que hay: "
            + ", ".join(conocidos))
    return sec.texto if crudo else _quitar_markdown(sec.texto)


def documento(language: str = "es") -> str:
    """El documento completo de un idioma, tal como se publica en docs/."""
    titulo = {"es": "# Fundamento físico de Olla-DFT",
              "en": "# The physics behind Olla-DFT"}
    partes = [titulo[language], ""]
    for area in AREAS:
        partes.append(_archivo(area, language).read_text(encoding="utf-8").rstrip())
        partes.append("")
    return "\n".join(partes).rstrip() + "\n"
