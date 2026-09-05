# Olla-DFT — command-line toolkit for Quantum ESPRESSO
# Copyright (C) 2026 Jorge Enrique González Sevilla
# SPDX-License-Identifier: AGPL-3.0-or-later
# This program is free software: you can redistribute it and/or modify it
# under the terms of the GNU Affero General Public License as published by the Free
# Software Foundation, either version 3 of the License, or (at your option)
# any later version. See the LICENSE file for details.

"""Estilo de figuras para publicación científica.

El problema que resuelve este módulo: una figura hecha con los valores por
defecto de matplotlib mide unas 6 x 4.6 pulgadas. Al insertarla en un
manuscrito hay que reducirla a 8.6 cm de ancho, y esa reducción encoge
también la tipografía: lo que se veía a 10 pt termina impreso a 4 pt, por
debajo del mínimo que aceptan las revistas.

La solución es construir la figura al tamaño físico final desde el inicio,
con las fuentes en puntos reales, y guardarla sin reescalar. Eso es lo que
hace `figure_size` junto con `apply`.

Referencias de anchos de columna (mm) tomadas de las guías de autor de cada
editorial; se pueden ajustar con `--width` si tu revista pide otra medida.
"""

from pathlib import Path

# El parser consulta los presets incluso en comandos sin figuras. Cargar
# Matplotlib y descubrir fuentes solo dentro de las funciones de dibujo.
from qekit.core import themes
from qekit.core.errors import ErrorDeUso

MM_PER_INCH = 25.4

# ----------------------------------------------------------------------
# Anchos de columna por editorial (mm): (una columna, 1.5 columnas, doble)
# ----------------------------------------------------------------------
JOURNALS = {
    "generic": (86.0, 140.0, 178.0),
    "elsevier": (90.0, 140.0, 190.0),
    "aps": (86.0, 129.0, 172.0),        # Physical Review
    "nature": (89.0, 120.0, 183.0),
    "acs": (84.6, 122.0, 177.8),
    "rsc": (83.0, 110.0, 171.0),
    "iop": (86.0, 129.0, 176.0),
    "wiley": (80.0, 120.0, 166.0),
}
WIDTH_KEYS = {"single": 0, "onehalf": 1, "double": 2}

# Patrones de línea como codificación secundaria: permiten distinguir las
# series en impresión a blanco y negro y refuerzan la separación de color.
DASHES = [
    (None, None),          # sólida
    (4.0, 1.6),
    (1.2, 1.2),
    (6.0, 1.6, 1.2, 1.6),
    (3.0, 1.4, 1.2, 1.4, 1.2, 1.4),
    (8.0, 2.0),
    (2.4, 1.2, 0.8, 1.2),
    (5.0, 1.2, 2.0, 1.2),
]

# Tintas neutras: el texto y los ejes nunca llevan el color de una serie.
INK = "#1a1a1a"
INK_SOFT = "#4d4d4d"
INK_FAINT = "#8c8c8c"
GRID = "#d9d9d9"

# ----------------------------------------------------------------------
# Presets de estilo. Los tamaños están en puntos tipográficos reales,
# tal como se verán impresos.
# ----------------------------------------------------------------------
STYLES = {
    "paper": dict(
        base=8.0, axes_label=9.0, tick=8.0, legend=7.5, title=9.0,
        line=1.0, axis_line=0.7, tick_len=3.0, tick_len_minor=1.8,
    ),
    "presentation": dict(
        base=11.0, axes_label=12.0, tick=11.0, legend=10.5, title=13.0,
        line=1.8, axis_line=1.1, tick_len=4.0, tick_len_minor=2.4,
    ),
    "poster": dict(
        base=14.0, axes_label=16.0, tick=14.0, legend=13.0, title=17.0,
        line=2.4, axis_line=1.4, tick_len=5.0, tick_len_minor=3.0,
    ),
}

FAMILIES = {
    # Se prueban en orden; se usa la primera instalada. Liberation Sans es
    # métricamente compatible con Arial, y Liberation Serif con Times.
    "sans": (
        ["Arial", "Helvetica", "Liberation Sans", "TeX Gyre Heros",
         "Nimbus Sans", "DejaVu Sans"],
        "stixsans",
    ),
    "serif": (
        ["Times New Roman", "Liberation Serif", "TeX Gyre Termes",
         "Nimbus Roman", "STIXGeneral", "DejaVu Serif"],
        "stix",
    ),
    # Computer Modern: la tipografía por omisión de LaTeX. Latin Modern es
    # su versión moderna y viene con cualquier distribución de TeX; si no
    # está, matplotlib trae sus propias cmr10.
    "latex": (
        ["Latin Modern Roman", "CMU Serif", "Computer Modern Roman",
         "cmr10", "DejaVu Serif"],
        "cm",
    ),
}

# Estado del tema activo. `apply` lo actualiza y el resto del código lee
# de aquí, de modo que cambiar de plantilla cambia toda la figura.
CURRENT = dict(themes.THEMES[themes.DEFAULT])
_PALETTE = list(themes.PALETTES["okabe-ito"])
USETEX = False


# ----------------------------------------------------------------------
# Texto compatible con y sin LaTeX
# ----------------------------------------------------------------------
# Con text.usetex, matplotlib pasa el texto a LaTeX tal cual, y ahí los
# caracteres griegos o el ångström en modo texto provocan un error de
# compilación. Estas tablas los traducen a las macros equivalentes.
_TEX_MAP = {
    "Γ": r"$\Gamma$", "Δ": r"$\Delta$", "Σ": r"$\Sigma$", "Λ": r"$\Lambda$",
    "Θ": r"$\Theta$", "Ω": r"$\Omega$", "Φ": r"$\Phi$", "Ψ": r"$\Psi$",
    "Ξ": r"$\Xi$", "Π": r"$\Pi$",
    "Å": r"\AA{}", "↑": r"$\uparrow$", "↓": r"$\downarrow$",
    "−": "-", "—": "---", "·": r"$\cdot$", "×": r"$\times$",
    "|": r"$|$", "<": r"$<$", ">": r"$>$", "~": r"$\sim$",
}

# Caracteres que LaTeX reserva y hay que escapar en modo texto.
_TEX_ESCAPE = set("#%&_{}$")


def tex_safe(text: str) -> str:
    """Adapta una etiqueta de texto plano al motor activo.

    Sin LaTeX se devuelve tal cual (matplotlib dibuja bien el Unicode). Con
    LaTeX hay que hacer dos cosas: escapar los caracteres que LaTeX reserva
    y traducir los símbolos que en modo texto no existen o se dibujan mal.
    El caso que más muerde es la barra de las discontinuidades del camino:
    un '|' suelto en modo texto sale como raya larga, y 'U|K' se leería
    'U—K'.
    """
    if not USETEX or not text:
        return text
    out = []
    for ch in str(text):
        if ch in _TEX_ESCAPE:
            out.append("\\" + ch)
        elif ch in _TEX_MAP:
            out.append(_TEX_MAP[ch])
        else:
            out.append(ch)
    return "".join(out)


def angstrom() -> str:
    """El símbolo de ångström en la forma que entienda el motor activo."""
    return r"\AA" if USETEX else "Å"


def available_font(candidates: list) -> str:
    """Primera fuente de la lista que esté instalada en el sistema."""
    from matplotlib import font_manager

    installed = {f.name for f in font_manager.fontManager.ttflist}
    for name in candidates:
        if name in installed:
            return name
    return candidates[-1]


def apply(theme=None, size: str = None, family: str = None,
          background: str = None, palette=None, usetex: bool = None,
          mono: bool = False, **overrides) -> dict:
    """Aplica una plantilla completa a matplotlib.

    `theme` es el nombre de una plantilla (o un dict, o la ruta a un JSON).
    El resto de argumentos la modifican puntualmente, de modo que se puede
    partir de una plantilla y cambiarle solo el fondo o la paleta.

    Devuelve un diccionario con la escala tipográfica y los datos del tema,
    que es lo que consultan las funciones de graficado.
    """
    import matplotlib

    global CURRENT, _PALETTE, INK, INK_SOFT, INK_FAINT, GRID, USETEX

    t = themes.load(theme, family=family, background=background,
                    palette=palette, usetex=usetex, size=size, **overrides)
    if mono:
        t["palette"] = "grayscale"
        t["dashes"] = "always"
        t["ink"] = "#000000"

    size_key = t.get("size", "paper")
    if size_key not in STYLES:
        raise ErrorDeUso(
            f"escala desconocida '{size_key}'. Opciones: {', '.join(STYLES)}"
        )
    fam = t.get("family", "sans")
    if fam not in FAMILIES:
        raise ErrorDeUso(
            f"familia desconocida '{fam}'. Opciones: {', '.join(FAMILIES)}"
        )

    s = STYLES[size_key]
    candidates, mathfont = FAMILIES[fam]
    chosen = available_font(candidates)
    rc_family = "sans-serif" if fam == "sans" else "serif"

    CURRENT = t
    _PALETTE = themes.resolve_palette(t.get("palette"))
    INK = t.get("ink", "#1a1a1a")
    INK_SOFT = t.get("ink_soft", "#4d4d4d")
    INK_FAINT = t.get("ink_faint", "#8c8c8c")
    GRID = t.get("grid_color", "#d9d9d9")
    bg = t.get("background", "#FFFFFF")
    USETEX = bool(t.get("usetex", False))

    matplotlib.rcParams.update({
        # --- tipografía ---
        "font.family": rc_family,
        f"font.{rc_family}": [chosen] + candidates,
        "font.size": s["base"],
        "axes.labelsize": s["axes_label"],
        "axes.titlesize": s["title"],
        "xtick.labelsize": s["tick"],
        "ytick.labelsize": s["tick"],
        "legend.fontsize": s["legend"],
        "axes.unicode_minus": True,      # signo menos tipográfico, no guion

        # --- incrustación de fuentes ---
        # Las revistas rechazan las fuentes Type 3; el tipo 42 (TrueType)
        # incrusta la fuente completa y es el que piden Elsevier, APS e IEEE.
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "svg.fonttype": "path",

        # --- colores ---
        "figure.facecolor": bg,
        "axes.facecolor": bg,
        "savefig.facecolor": bg,
        "axes.edgecolor": INK,
        "axes.labelcolor": INK,
        "axes.titlecolor": INK,
        "text.color": INK,
        "xtick.color": INK,
        "ytick.color": INK,
        "xtick.labelcolor": INK,
        "ytick.labelcolor": INK,
        "grid.color": GRID,

        # --- ejes ---
        "axes.linewidth": s["axis_line"],
        "axes.axisbelow": True,
        "grid.linewidth": s["axis_line"] * 0.7,

        # --- líneas ---
        "lines.linewidth": s["line"],
        "lines.markersize": 3.5,
        "lines.solid_capstyle": "round",

        # --- leyenda ---
        "legend.frameon": False,
        "legend.handlelength": 1.8,
        "legend.handletextpad": 0.5,
        "legend.labelspacing": 0.25,
        "legend.columnspacing": 1.0,
        "legend.borderaxespad": 0.4,
        "legend.labelcolor": INK,

        # --- guardado ---
        "savefig.dpi": 600,          # por encima del mínimo de 300 dpi
        "savefig.transparent": False,
    })

    # --- marcas de escala ---
    direction = t.get("ticks", "in")
    all_sides = bool(t.get("all_sides", True))
    matplotlib.rcParams.update({
        "xtick.direction": direction,
        "ytick.direction": direction,
        "xtick.top": all_sides,
        "ytick.right": all_sides,
        "xtick.minor.visible": True,
        "ytick.minor.visible": True,
        "xtick.major.size": s["tick_len"],
        "ytick.major.size": s["tick_len"],
        "xtick.minor.size": s["tick_len_minor"],
        "ytick.minor.size": s["tick_len_minor"],
        "xtick.major.width": s["axis_line"],
        "ytick.major.width": s["axis_line"],
        "xtick.minor.width": s["axis_line"] * 0.8,
        "ytick.minor.width": s["axis_line"] * 0.8,
    })

    # --- rejilla ---
    grid_mode = t.get("grid", False)
    matplotlib.rcParams["axes.grid"] = bool(grid_mode)
    if grid_mode:
        matplotlib.rcParams["axes.grid.axis"] = (
            "both" if grid_mode is True else str(grid_mode)
        )

    # --- texto matemático ---
    if USETEX:
        preamble = t.get("latex_preamble", r"\usepackage{amsmath}\usepackage{amssymb}")
        if fam == "sans":
            preamble += r"\usepackage{helvet}\renewcommand{\familydefault}{\sfdefault}"
        matplotlib.rcParams.update({
            "text.usetex": True,
            "text.latex.preamble": preamble,
        })
    else:
        matplotlib.rcParams["text.usetex"] = False
        if fam == "latex":
            # Computer Modern sin depender de una instalación de TeX
            matplotlib.rcParams["mathtext.fontset"] = "cm"
        else:
            # El texto matemático debe usar la MISMA familia que el resto: con
            # los conjuntos predefinidos ('stixsans') matplotlib mezcla glifos
            # de STIX serif, y la figura acaba con dos tipografías distintas.
            matplotlib.rcParams.update({
                "mathtext.fontset": "custom",
                "mathtext.rm": chosen,
                "mathtext.it": f"{chosen}:italic",
                "mathtext.bf": f"{chosen}:bold",
                "mathtext.sf": chosen,
                "mathtext.cal": f"{chosen}:italic",
                "mathtext.tt": "DejaVu Sans Mono",
                "mathtext.default": "it",
            })

    return dict(s, font=chosen, family=fam, theme=t, background=bg,
                usetex=USETEX)


def finish_axes(ax):
    """Aplica al eje lo que no cabe en rcParams: modo de marco y rejilla.

    Con spines='lr' se ocultan el marco superior y el derecho, que es el
    aspecto habitual de informes y diapositivas; 'box' deja el marco
    completo, que es la convención en revistas de física.
    """
    import matplotlib

    mode = CURRENT.get("spines", "box")
    if mode in ("lr", "left-bottom", "open"):
        for side in ("top", "right"):
            ax.spines[side].set_visible(False)
        ax.tick_params(top=False, right=False, which="both")
    grid_mode = CURRENT.get("grid", False)
    if grid_mode:
        axis = "both" if grid_mode is True else str(grid_mode)
        ax.grid(True, axis=axis, color=GRID,
                lw=matplotlib.rcParams["grid.linewidth"], zorder=0)
        ax.set_axisbelow(True)
    return ax


# ----------------------------------------------------------------------
# Tamaño físico
# ----------------------------------------------------------------------
def width_mm(width="single", journal: str = "generic") -> float:
    """Ancho en milímetros a partir de un nombre o un número."""
    if isinstance(width, (int, float)):
        return float(width)
    w = str(width).strip().lower()
    if w.replace(".", "", 1).isdigit():
        return float(w)
    if journal not in JOURNALS:
        raise ErrorDeUso(
            f"revista desconocida '{journal}'. Opciones: {', '.join(JOURNALS)}"
        )
    if w not in WIDTH_KEYS:
        raise ErrorDeUso(
            f"ancho desconocido '{width}'. Usa single/onehalf/double "
            "o un número en mm."
        )
    return JOURNALS[journal][WIDTH_KEYS[w]]


def figure_size(width="single", journal: str = "generic",
                aspect: float = 0.75, height_mm: float = None) -> tuple:
    """Tamaño de figura en pulgadas, al ancho físico real de la columna.

    `aspect` es alto/ancho. 0.75 es un rectángulo cómodo; para figuras
    anchas de dos columnas suele convenir 0.42–0.5.
    """
    w_mm = width_mm(width, journal)
    h_mm = height_mm if height_mm else w_mm * aspect
    return (w_mm / MM_PER_INCH, h_mm / MM_PER_INCH)


# ----------------------------------------------------------------------
# Colores y patrones
# ----------------------------------------------------------------------
def palette(n: int = None, mono: bool = False) -> list:
    """Colores de serie del tema activo, en orden fijo (nunca reciclado
    arbitrariamente: el color sigue a la serie, no a su posición)."""
    base = list(themes.PALETTES["grayscale"]) if mono else list(_PALETTE)
    if not base:
        base = list(themes.PALETTES["okabe-ito"])
    if n is None:
        return base
    if n <= len(base):
        return base[:n]
    # Más series de las que la paleta admite con seguridad: se repiten los
    # colores, pero el patrón de línea las sigue distinguiendo.
    return [base[i % len(base)] for i in range(n)]


def dash(index: int):
    """Patrón de guiones para la serie `index` (None = línea sólida)."""
    on_off = DASHES[index % len(DASHES)]
    return None if on_off[0] is None else list(on_off)


def use_dashes(n_series: int, mode: str = "auto", mono: bool = False) -> bool:
    """¿Conviene distinguir también por patrón de línea?

    'auto' los activa a partir de cinco series, que es donde la paleta deja
    de garantizar separación bajo daltonismo, y siempre en monocromo. Si la
    plantilla activa pide otra cosa ('always'/'never'), manda la plantilla.
    """
    if mono:
        return True
    if mode == "auto":
        mode = CURRENT.get("dashes", "auto")
    if mode == "always":
        return True
    if mode == "never":
        return False
    return n_series >= 5


def style_line(index: int, n_series: int, dash_mode: str = "auto",
               mono: bool = False) -> dict:
    """kwargs de color y guiones para la serie `index`."""
    kw = {"color": palette(max(n_series, index + 1), mono=mono)[index]}
    if use_dashes(n_series, dash_mode, mono):
        d = dash(index)
        if d:
            kw["dashes"] = d
    return kw


# ----------------------------------------------------------------------
# Etiquetas de panel y guardado
# ----------------------------------------------------------------------
def panel_label(ax, text: str, loc: str = "upper left", pad: float = 0.03,
                weight: str = "bold"):
    """Etiqueta (a), (b)... para figuras de varios paneles."""
    import matplotlib

    positions = {
        "upper left": (pad, 1 - pad, "left", "top"),
        "upper right": (1 - pad, 1 - pad, "right", "top"),
        "lower left": (pad, pad, "left", "bottom"),
        "lower right": (1 - pad, pad, "right", "bottom"),
    }
    x, y, ha, va = positions.get(loc, positions["upper left"])
    return ax.text(
        x, y, text, transform=ax.transAxes, ha=ha, va=va,
        fontweight=weight, color=INK,
        fontsize=matplotlib.rcParams["axes.labelsize"],
        zorder=10,
    )


def parse_formats(fmt) -> list:
    """'pdf,png' -> ['pdf', 'png']."""
    if isinstance(fmt, (list, tuple)):
        return list(fmt)
    return [f.strip().lower() for f in str(fmt).split(",") if f.strip()]


def save(fig, outbase: str, formats="pdf,png", dpi: int = None,
         modulo: str = "", params: dict = None) -> list:
    """Guarda la figura en varios formatos conservando el tamaño físico.

    No se usa bbox_inches='tight': recortar el margen cambiaría el ancho
    final y la figura ya no mediría exactamente el ancho de columna pedido.
    El ajuste interno lo hace el layout 'constrained'.

    La procedencia (versión de Olla-DFT, fecha, comando y parámetros) se
    incrusta en los METADATOS del archivo, no sobre la imagen: no estorba en
    la revista y sigue ahí para saber qué versión hizo la figura.
    """
    import matplotlib

    base = Path(outbase)
    if base.suffix.lower().lstrip(".") in ("pdf", "png", "svg", "eps", "tif", "tiff"):
        base = base.with_suffix("")
    base.parent.mkdir(parents=True, exist_ok=True)

    from qekit.core import provenance

    written = []
    for ext in parse_formats(formats):
        target = base.with_suffix(f".{ext}")
        kwargs = {}
        if ext in ("png", "tif", "tiff"):
            kwargs["dpi"] = dpi or matplotlib.rcParams["savefig.dpi"]
        if ext in ("pdf", "png", "svg", "eps"):
            # cada backend acepta un juego distinto de claves; si alguna no
            # le gusta, la figura importa más que el metadato
            try:
                fig.savefig(
                    target,
                    metadata=provenance.figure_metadata(modulo, params, ext),
                    **kwargs)
            except (TypeError, ValueError):
                fig.savefig(target, **kwargs)
        else:
            fig.savefig(target, **kwargs)
        written.append(str(target))
    return written


def new_figure(width="single", journal: str = "generic", aspect: float = 0.75,
               height_mm: float = None, **kwargs):
    """Crea figura y ejes al tamaño físico correcto, con layout constrained."""
    import matplotlib.pyplot as plt

    size = figure_size(width, journal, aspect, height_mm)
    fig, ax = plt.subplots(figsize=size, layout="constrained", **kwargs)
    fig.get_layout_engine().set(w_pad=0.012, h_pad=0.012, hspace=0.0, wspace=0.0)
    finish_axes(ax)
    return fig, ax
