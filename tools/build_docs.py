#!/usr/bin/env python3
# Olla-DFT — command-line toolkit for Quantum ESPRESSO
# Copyright (C) 2026 Jorge Enrique González Sevilla
# SPDX-License-Identifier: AGPL-3.0-or-later
# This program is free software: you can redistribute it and/or modify it
# under the terms of the GNU Affero General Public License as published by the Free
# Software Foundation, either version 3 of the License, or (at your option)
# any later version. See the LICENSE file for details.

"""Genera la documentación derivada del código, en el idioma del repositorio.

    docs/TEORIA.md   (o THEORY.md)    <- qekit/data/theory/*.md
    docs/COMANDOS.md (o COMMANDS.md)  <- el árbol de argparse (olla-dft --help)

El idioma lo decide ``qekit.core.i18n.DEFAULT_LANGUAGE``: el repositorio en
español publica TEORIA.md y COMANDOS.md; el repositorio en inglés, THEORY.md
y COMMANDS.md. Con ``--all`` se generan los dos idiomas.

Correr después de tocar la teoría o de añadir/renombrar un comando:

    python tools/build_docs.py

tests/test_teoria.py y tests/test_docs.py comprueban que estén al día.
"""

import argparse
import sys
from pathlib import Path

# Build this checkout, even when another release is installed in site-packages.
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from qekit import __command_name__, __version__
from qekit.cli import COMMAND_GROUPS, build_parser, _menu_labels
from qekit.core.i18n import DEFAULT_LANGUAGE
from qekit.modules import theory

RAIZ = Path(__file__).resolve().parent.parent
DOCS = RAIZ / "docs"

TEXTOS = {
    "es": {
        "titulo": "# Referencia de comandos de Olla-DFT",
        "intro": (f"Los {{n}} subcomandos de `{__command_name__}`, agrupados por área, con "
                  "sus opciones. Generado del propio código con "
                  "`python tools/build_docs.py`; la misma información sale en la "
                  f"terminal con `{__command_name__} COMANDO --help` y, navegable, con "
                  f"`{__command_name__} docs`."),
        "uso": "Uso", "argumentos": "Argumentos", "opciones": "Opciones",
        "teoria": "Fundamento físico", "ninguno": "(ninguno)",
        "col": "| Opción | Descripción |\n|---|---|",
        "default": "default",
        "indice": "## Índice",
    },
    "en": {
        "titulo": "# Olla-DFT command reference",
        "intro": (f"The {{n}} `{__command_name__}` subcommands, grouped by area, with "
                  "their options. Generated from the code itself with "
                  "`python tools/build_docs.py`; the same information is printed by "
                  f"`{__command_name__} COMMAND --help --language en` and, as a browsable page, by "
                  f"`{__command_name__} docs --language en`."),
        "uso": "Usage", "argumentos": "Arguments", "opciones": "Options",
        "teoria": "Physics", "ninguno": "(none)",
        "col": "| Option | Description |\n|---|---|",
        "default": "default",
        "indice": "## Index",
    },
}


def _esc(text: str) -> str:
    return (text or "").replace("|", "\\|").replace("\n", " ")


def _opcion(action) -> str:
    if action.option_strings:
        nombre = ", ".join(action.option_strings)
        if action.metavar or (action.nargs != 0 and action.type is not None
                              and not action.choices):
            nombre += f" {action.metavar or action.dest.upper()}"
        elif action.choices:
            nombre += " {" + ",".join(str(c) for c in action.choices) + "}"
        return f"`{nombre}`"
    return f"`{action.metavar or action.dest}`"


def comandos_md(language: str) -> str:
    t = TEXTOS[language]
    parser = build_parser(language)
    sub = parser._subparsers._group_actions[0]
    resumen = {ca.dest.split(" ")[0]: ca.help for ca in sub._choices_actions}
    labels = _menu_labels(language)
    con_teoria = {c for s in theory.secciones(language) for c in s.comandos}
    n = sum(len(cmds) for _, cmds in COMMAND_GROUPS)
    out = [t["titulo"], "", t["intro"].format(n=n), "", t["indice"], ""]
    for title, cmds in COMMAND_GROUPS:
        out.append(f"- **{labels['catalog_groups'].get(title, title)}**: "
                   + ", ".join(f"[`{c}`](#{c})" for c in cmds))
    out.append("")
    for title, cmds in COMMAND_GROUPS:
        out += [f"## {labels['catalog_groups'].get(title, title)}", ""]
        for name in cmds:
            p = sub.choices[name]
            out += [f"### `{name}`", "", f"{resumen.get(name, '')}", ""]
            uso = p.format_usage().replace("usage: ", "").replace("Usage: ", "")
            out += [f"**{t['uso']}:** `{' '.join(uso.split())}`", ""]
            pos = [a for a in p._actions if not a.option_strings]
            opts = [a for a in p._actions
                    if a.option_strings and a.dest != "help"
                    and a.help is not argparse.SUPPRESS]
            if pos:
                out += [f"**{t['argumentos']}:**", ""]
                for a in pos:
                    ch = (" {" + ",".join(str(c) for c in a.choices) + "}"
                          if a.choices else "")
                    out.append(f"- {_opcion(a)}{ch} — {_esc(a.help)}")
                out.append("")
            if opts:
                out += [f"**{t['opciones']}:**", "", t["col"]]
                for a in opts:
                    d = a.default
                    extra = ""
                    if d not in (None, False, "", argparse.SUPPRESS) and \
                            f"{t['default']}" not in (a.help or ""):
                        extra = f" ({t['default']}: `{d}`)"
                    out.append(f"| {_opcion(a)} | {_esc(a.help)}{extra} |")
                out.append("")
            if name in con_teoria:
                doc = "TEORIA.md" if language == "es" else "THEORY.md"
                out += [f"**{t['teoria']}:** [`{__command_name__} teoria {name}`]({doc})", ""]
    out.append(f"---\n\n*Olla-DFT {__version__}*")
    return "\n".join(out) + "\n"


NOMBRES = {"es": ("TEORIA.md", "COMANDOS.md"), "en": ("THEORY.md", "COMMANDS.md")}


def salidas(language: str) -> dict:
    teoria, comandos = NOMBRES[language]
    return {DOCS / teoria: theory.documento(language),
            DOCS / comandos: comandos_md(language)}


def main(argv=None) -> int:
    import sys
    argv = sys.argv[1:] if argv is None else argv
    idiomas = ("es", "en") if "--all" in argv else (DEFAULT_LANGUAGE,)
    DOCS.mkdir(parents=True, exist_ok=True)
    todas = {}
    for lang in idiomas:
        todas.update(salidas(lang))
    for path, texto in todas.items():
        path.write_text(texto, encoding="utf-8")
        print(f"{path.relative_to(RAIZ)}: {path.stat().st_size // 1024} KB")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
