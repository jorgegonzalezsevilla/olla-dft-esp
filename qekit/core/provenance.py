# Olla-DFT — command-line toolkit for Quantum ESPRESSO
# Copyright (C) 2026 Jorge Enrique González Sevilla
# SPDX-License-Identifier: GPL-3.0-or-later
# This program is free software: you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by the Free
# Software Foundation, either version 3 of the License, or (at your option)
# any later version. See the LICENSE file for details.

"""Procedencia: de dónde salió cada número.

Sin repositorio público no hay commit que citar, así que el propio archivo
tiene que decir con qué versión de Olla-DFT y con qué parámetros se generó. Es
lo que permite, meses después, rastrear una cifra de una tesis o un artículo
hasta el cálculo que la produjo — y darse cuenta de que una figura vieja se
hizo con una versión que tenía un error.

El encabezado va como comentarios '#' al principio de cada .dat, y en las
figuras va en los metadatos del PNG/PDF (invisible al imprimir, legible con
`pdfinfo`, `exiftool` o `olla-dft info --figura`).
"""

import datetime as _dt
import shlex
import sys

from qekit import __command_name__, __product_name__, __version__

_ARGV = None            # se fija una vez al entrar por la CLI


def record_argv(argv=None):
    """Guarda la línea de comandos con que se invocó Olla-DFT."""
    global _ARGV
    _ARGV = list(argv if argv is not None else sys.argv)


def command_line() -> str:
    """La invocación como texto, o '' si se usó como biblioteca."""
    if not _ARGV:
        return ""
    args = _ARGV[1:] if _ARGV and _ARGV[0].endswith(
        ("cli.py", "qekit", __command_name__)) else _ARGV
    return __command_name__ + " " + " ".join(shlex.quote(str(a)) for a in args)


def _timestamp() -> str:
    return _dt.datetime.now(_dt.timezone.utc).strftime("%Y-%m-%d %H:%M UTC")


def fields(modulo: str = "", params: dict = None, extra: dict = None) -> dict:
    """Diccionario de procedencia (base para encabezado y metadatos)."""
    d = {
        "qekit_version": __version__,
        "generado": _timestamp(),
    }
    if modulo:
        d["modulo"] = modulo
    cmd = command_line()
    if cmd:
        d["comando"] = cmd
    for src in (params, extra):
        if src:
            for k, v in src.items():
                if v is not None and v != "":
                    d[str(k)] = v
    return d


def header(modulo: str = "", params: dict = None, titulo: str = "") -> str:
    """Bloque de comentarios '#' para encabezar un archivo de datos.

    No lleva salto final: quien lo use añade el resto del encabezado.
    """
    d = fields(modulo, params)
    lines = []
    if titulo:
        lines.append(f"# {titulo}")
    lines.append(f"# {__product_name__} {d['qekit_version']} — generado {d['generado']}")
    if "comando" in d:
        lines.append(f"# comando: {d['comando']}")
    resto = {k: v for k, v in d.items()
             if k not in ("qekit_version", "generado", "comando", "modulo")}
    if resto:
        partes = [f"{k} = {v}" for k, v in resto.items()]
        # partir en líneas de ~76 caracteres para que no se desborde
        linea = "# parámetros:"
        for parte in partes:
            if len(linea) + len(parte) + 2 > 76:
                lines.append(linea)
                linea = "#   " + parte
            else:
                linea += (" " if linea.endswith(":") else "; ") + parte
        lines.append(linea)
    return "\n".join(lines)


def header_plain(modulo: str = "", params: dict = None,
                 titulo: str = "") -> str:
    """Igual que header() pero SIN el '#' inicial en cada línea.

    Para np.savetxt(..., comments="# "), que ya lo antepone; mezclarlos
    produce '## '.
    """
    txt = header(modulo, params, titulo)
    return "\n".join(l[2:] if l.startswith("# ") else l.lstrip("#")
                     for l in txt.splitlines())


def figure_metadata(modulo: str = "", params: dict = None,
                    fmt: str = "pdf") -> dict:
    """Metadatos para savefig(); las claves válidas dependen del formato.

    PDF usa las claves del estándar (Title, Author, Subject, Keywords,
    Creator, Producer). PNG usa palabras clave tEXt, donde 'Subject' NO es
    estándar y el escritor la descarta en silencio — ahí el detalle va en
    'Description', que sí lo es.
    """
    d = fields(modulo, params)
    resto = {k: v for k, v in d.items()
             if k not in ("qekit_version", "generado", "comando", "modulo")}
    detalle = "; ".join(f"{k}={v}" for k, v in resto.items())
    meta = {
        "Creator": f"{__product_name__} {d['qekit_version']}",
        "Producer": f"{__product_name__} {d['qekit_version']}",
        "Title": (modulo or f"figura {__product_name__}"),
    }
    if d.get("comando"):
        meta["Author"] = d["comando"]
    if fmt.lower() == "png":
        meta["Creation Time"] = d["generado"]
        meta["Software"] = f"{__product_name__} {d['qekit_version']}"
        if detalle:
            meta["Description"] = detalle
    else:
        meta["Keywords"] = d["generado"]
        if detalle:
            meta["Subject"] = detalle
    return meta
