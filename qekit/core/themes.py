# Olla-DFT — command-line toolkit for Quantum ESPRESSO
# Copyright (C) 2026 Jorge Enrique González Sevilla
# SPDX-License-Identifier: AGPL-3.0-or-later
# This program is free software: you can redistribute it and/or modify it
# under the terms of the GNU Affero General Public License as published by the Free
# Software Foundation, either version 3 of the License, or (at your option)
# any later version. See the LICENSE file for details.

"""Plantillas visuales de Olla-DFT.

Una plantilla agrupa todo lo que define el aspecto de una figura — fondo,
tintas, paleta, ejes, tipografía — para que cambiar de estilo sea una sola
bandera y no una docena de ajustes sueltos.

Las plantillas se combinan con el preset de tamaño (`paper`, `presentation`,
`poster`), que solo controla la escala tipográfica y los grosores. Así, la
misma plantilla sirve para un artículo y para una diapositiva sin duplicar
definiciones.

Paletas: los colores no se eligieron a ojo. Se verificaron con un validador
que simula protanopia y deuteranopia y mide la separación en OKLab; los
números están anotados en cada entrada.
"""

import json
from pathlib import Path
from qekit.core.errors import ErrorDeUso

from qekit.core import plataforma

#: Plantillas del usuario, en la carpeta que corresponde a cada sistema.
USER_DIR = plataforma.dir_config() / "templates"

# ----------------------------------------------------------------------
# Paletas
# ----------------------------------------------------------------------
PALETTES = {
    # Okabe & Ito (2008), estándar para figuras accesibles con daltonismo.
    # Reordenada: las primeras cuatro son las de mayor separación
    # (ΔE OKLab = 11.0 bajo deuteranopia, contra el objetivo de 8).
    "okabe-ito": [
        "#0072B2", "#D55E00", "#009E73", "#E69F00",
        "#CC79A7", "#56B4E9", "#000000", "#F0E442",
    ],
    # Versión para fondo oscuro. No es la clara invertida: se recalcularon
    # los pasos sobre los MISMOS tonos, dentro de la banda de luminosidad
    # que exige un fondo oscuro. Las cuatro primeras pasan todas las
    # comprobaciones (ΔE 8.9 deuteranopia, 18.4 visión normal, contraste
    # >= 3:1); de la quinta en adelante el patrón de línea es obligatorio,
    # y por eso el tema oscuro activa los guiones antes.
    "okabe-ito-dark": [
        "#027BBF", "#C95907", "#0EA276", "#C91290",
        "#069CA3", "#9C8A08", "#8FA6B8", "#E0C64A",
    ],
    # Escala de grises para impresión sin color; la identidad la llevan
    # los patrones de línea, no el tono.
    "grayscale": ["#1a1a1a", "#4d4d4d", "#808080", "#a6a6a6"],
}

# ----------------------------------------------------------------------
# Plantillas
# ----------------------------------------------------------------------
# Claves de cada plantilla:
#   family      sans | serif | latex        (latex = Computer Modern)
#   usetex      True para renderizar con LaTeX de verdad (requiere TeX)
#   background  color de fondo de la figura y de los ejes
#   ink         color de texto, ejes y marcas
#   palette     nombre en PALETTES o lista de colores
#   spines      box (marco completo) | lr (izquierda e inferior)
#   ticks       in | out
#   all_sides   marcas también arriba y a la derecha
#   grid        False | 'x' | 'y' | 'both'
#   size        preset de escala por defecto
#   dashes      auto | always | never
#   fill_alpha  opacidad del relleno bajo la DOS total

THEMES = {
    "journal": dict(
        description="Artículo científico: fondo blanco, sans, marco completo.",
        family="sans", usetex=False, background="#FFFFFF",
        ink="#1a1a1a", ink_soft="#4d4d4d", ink_faint="#8c8c8c", grid_color="#d9d9d9",
        palette="okabe-ito", spines="box", ticks="in", all_sides=True,
        grid=False, size="paper", dashes="auto", fill_alpha=0.10,
    ),
    "latex": dict(
        description="Tipografía Computer Modern, igual que un documento LaTeX, "
                    "sin necesidad de tener LaTeX instalado.",
        family="latex", usetex=False, background="#FFFFFF",
        ink="#1a1a1a", ink_soft="#4d4d4d", ink_faint="#8c8c8c", grid_color="#d9d9d9",
        palette="okabe-ito", spines="box", ticks="in", all_sides=True,
        grid=False, size="paper", dashes="auto", fill_alpha=0.10,
    ),
    "latex-true": dict(
        description="Renderizado con LaTeX real: la figura usa exactamente el "
                    "mismo motor tipográfico que tu manuscrito. Requiere una "
                    "instalación de TeX con dvipng.",
        family="latex", usetex=True, background="#FFFFFF",
        ink="#1a1a1a", ink_soft="#4d4d4d", ink_faint="#8c8c8c", grid_color="#d9d9d9",
        palette="okabe-ito", spines="box", ticks="in", all_sides=True,
        grid=False, size="paper", dashes="auto", fill_alpha=0.10,
    ),
    "minimal": dict(
        description="Solo ejes izquierdo e inferior, rejilla tenue. Aire limpio "
                    "para informes y tesis.",
        family="sans", usetex=False, background="#FFFFFF",
        ink="#262626", ink_soft="#595959", ink_faint="#a6a6a6", grid_color="#e8e8e8",
        palette="okabe-ito", spines="lr", ticks="out", all_sides=False,
        grid="both", size="paper", dashes="auto", fill_alpha=0.09,
    ),
    "dark": dict(
        description="Fondo oscuro para diapositivas: paleta propia validada "
                    "contra la superficie oscura, no la clara invertida.",
        family="sans", usetex=False, background="#1a1a19",
        ink="#f2f2f0", ink_soft="#c4c4c0", ink_faint="#8a8a86", grid_color="#3a3a38",
        palette="okabe-ito-dark", spines="box", ticks="in", all_sides=True,
        grid=False, size="presentation", dashes="auto", fill_alpha=0.16,
    ),
    "slides": dict(
        description="Diapositivas en claro: tipografía y trazos grandes, fondo "
                    "hueso que cansa menos que el blanco puro proyectado.",
        family="sans", usetex=False, background="#FBFBF9",
        ink="#1a1a1a", ink_soft="#4d4d4d", ink_faint="#8c8c8c", grid_color="#e0e0dc",
        palette="okabe-ito", spines="lr", ticks="out", all_sides=False,
        grid="y", size="presentation", dashes="auto", fill_alpha=0.12,
    ),
    "poster": dict(
        description="Cartel: todo a mayor escala, para leerse a un metro.",
        family="sans", usetex=False, background="#FFFFFF",
        ink="#1a1a1a", ink_soft="#4d4d4d", ink_faint="#8c8c8c", grid_color="#dcdcdc",
        palette="okabe-ito", spines="box", ticks="in", all_sides=True,
        grid=False, size="poster", dashes="auto", fill_alpha=0.12,
    ),
    "mono": dict(
        description="Monocromo para revistas que cobran el color: la identidad "
                    "de cada serie la lleva el patrón de línea.",
        family="sans", usetex=False, background="#FFFFFF",
        ink="#000000", ink_soft="#4d4d4d", ink_faint="#8c8c8c", grid_color="#d9d9d9",
        palette="grayscale", spines="box", ticks="in", all_sides=True,
        grid=False, size="paper", dashes="always", fill_alpha=0.08,
    ),
    "mono-latex": dict(
        description="Monocromo con tipografía Computer Modern.",
        family="latex", usetex=False, background="#FFFFFF",
        ink="#000000", ink_soft="#4d4d4d", ink_faint="#8c8c8c", grid_color="#d9d9d9",
        palette="grayscale", spines="box", ticks="in", all_sides=True,
        grid=False, size="paper", dashes="always", fill_alpha=0.08,
    ),
}

DEFAULT = "journal"


# ----------------------------------------------------------------------
# Carga y combinación
# ----------------------------------------------------------------------
def user_templates() -> dict:
    """Plantillas del usuario, en <config>/templates/*.json.

    La carpeta depende del sistema: `olla-dft sistema` dice cuál es la tuya.
    """
    found = {}
    if USER_DIR.is_dir():
        for f in sorted(USER_DIR.glob("*.json")):
            try:
                data = json.loads(f.read_text())
            except (OSError, json.JSONDecodeError):
                continue
            if isinstance(data, dict):
                found[f.stem] = data
    return found


def names() -> list:
    """Todas las plantillas disponibles: las de fábrica y las del usuario."""
    return sorted(set(THEMES) | set(user_templates()))


def resolve_palette(value, mono: bool = False) -> list:
    """Nombre de paleta o lista de colores -> lista de colores."""
    if mono:
        return list(PALETTES["grayscale"])
    if value is None:
        return list(PALETTES["okabe-ito"])
    if isinstance(value, (list, tuple)):
        return list(value)
    text = str(value).strip()
    if text in PALETTES:
        return list(PALETTES[text])
    if "," in text or text.startswith("#"):
        colors = [c.strip() for c in text.split(",") if c.strip()]
        bad = [c for c in colors if not c.startswith("#")]
        if bad:
            raise ErrorDeUso(
                f"colores no reconocidos: {', '.join(bad)}. "
                "Usa valores hexadecimales (#0072B2) separados por coma."
            )
        return colors
    raise ErrorDeUso(
        f"paleta desconocida '{value}'. Disponibles: {', '.join(sorted(PALETTES))}, "
        "o una lista de colores hexadecimales separados por coma."
    )


def load(name=None, **overrides) -> dict:
    """Devuelve una plantilla completa, aplicando las modificaciones dadas.

    `name` puede ser el nombre de una plantilla de fábrica, el de una del
    usuario, o la ruta a un archivo JSON. Los `overrides` (por ejemplo
    background='#000000') se aplican encima.
    """
    name = name or DEFAULT
    base = dict(THEMES[DEFAULT])

    if isinstance(name, dict):
        base.update(name)
    else:
        text = str(name)
        path = Path(text)
        if path.suffix.lower() == ".json" and path.is_file():
            base.update(json.loads(path.read_text()))
        elif text in THEMES:
            base.update(THEMES[text])
        else:
            users = user_templates()
            if text in users:
                # una plantilla de usuario puede heredar de otra con "extends"
                data = users[text]
                parent = data.get("extends")
                if parent and parent in THEMES:
                    base.update(THEMES[parent])
                base.update({k: v for k, v in data.items() if k != "extends"})
            else:
                # La confusión típica: '-t nature'. `nature` existe, pero es
                # una REVISTA (ancho de columna y tipografía del formato),
                # no una plantilla (colores y estilo). Decirlo aquí ahorra
                # ir a buscar la diferencia en el manual.
                extra = ""
                try:
                    from qekit.core.style import JOURNALS
                    if text in JOURNALS:
                        extra = (f"\n'{text}' existe, pero es una REVISTA, no "
                                 f"una plantilla: usa --journal {text}.")
                except Exception:                      # noqa: BLE001
                    pass
                if not extra and text in PALETTES:
                    extra = (f"\n'{text}' existe, pero es una PALETA: usa "
                             f"--palette {text}.")
                raise ErrorDeUso(
                    f"plantilla desconocida '{text}'. Disponibles: "
                    f"{', '.join(names())}" + extra
                )

    for key, value in overrides.items():
        if value is not None:
            base[key] = value
    return base


def describe(name: str) -> str:
    """Descripción legible de una plantilla."""
    t = load(name)
    lines = [f"Plantilla: {name}"]
    if t.get("description"):
        lines.append(f"  {t['description']}")
    colors = resolve_palette(t.get("palette"))
    lines += [
        f"  tipografía   : {t.get('family')}"
        + ("  (LaTeX real)" if t.get("usetex") else ""),
        f"  fondo        : {t.get('background')}",
        f"  tinta        : {t.get('ink')}",
        f"  paleta       : {t.get('palette')}",
        f"                 {' '.join(colors[:6])}",
        f"  ejes         : {t.get('spines')} · marcas {t.get('ticks')}"
        + (" en los cuatro lados" if t.get("all_sides") else ""),
        f"  rejilla      : {t.get('grid')}",
        f"  escala       : {t.get('size')}",
        f"  guiones      : {t.get('dashes')}",
    ]
    return "\n".join(lines)


def export(name: str, path: str = None) -> str:
    """Guarda una plantilla como JSON para que el usuario la modifique."""
    data = dict(load(name))
    data["extends"] = name if name in THEMES else DEFAULT
    target = Path(path) if path else USER_DIR / f"{name}-copia.json"
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_text(json.dumps(data, indent=2, ensure_ascii=False) + "\n")
    return str(target)
