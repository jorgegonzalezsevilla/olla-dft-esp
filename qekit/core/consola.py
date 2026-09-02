# Olla-DFT — command-line toolkit for Quantum ESPRESSO
# Copyright (C) 2026 Jorge Enrique González Sevilla
# SPDX-License-Identifier: GPL-3.0-or-later
# This program is free software: you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by the Free
# Software Foundation, either version 3 of the License, or (at your option)
# any later version. See the LICENSE file for details.

"""Que la salida se vea igual en cualquier terminal, sin reventar.

EL PROBLEMA
-----------
Los informes de Olla-DFT están llenos de caracteres que la física necesita: Å,
α, ε, Ω, ħ, ², ←, →, ①, ✓. En Linux y macOS eso es UTF-8 y no pasa nada. En
la consola heredada de Windows la página de códigos por defecto es cp1252 y
`print` lanza UnicodeEncodeError. No es un problema estético: el comando
MUERE a media salida, con código 1 y una traza, y el usuario cree que el
cálculo falló.

Se comprobó: con la salida forzada a cp1252, `olla-dft info`, `olla-dft selftest` y
`olla-dft recetas` reventaban las tres. La única que sobrevivía era la que no
imprimía ni una Å.

LA SOLUCIÓN, EN DOS ESCALONES
-----------------------------
1. **Intentar UTF-8.** Windows 10 en adelante lo admite, y Python permite
   reconfigurar el flujo en marcha. Si sale bien, no se pierde nada: se ve
   la Å de verdad.
2. **Si no se puede, transliterar.** Cuando el destino no admite un carácter
   —una consola vieja, una redirección a un archivo con otra codificación—,
   en vez de morir se sustituye por su equivalente ASCII: Å → A, α → alpha,
   → → ->, ① → (1). La salida pierde belleza y conserva TODA la información,
   que es exactamente el orden de prioridades correcto.

Lo que NO se hace es poner `errors="replace"`. Eso convertiría cada Å en un
'?' y dejaría informes con "el parametro de red es 5.43 ?", que es peor que
no imprimir nada porque parece un dato corrupto.
"""

import sys

# Sustituciones para cuando el terminal no admite el carácter. La regla es
# que se pueda LEER en voz alta y signifique lo mismo, no que sea bonito.
TRANSLITERACION = {
    # unidades y constantes
    "Å": "A", "å": "a", "Ω": "Ohm", "ω": "omega", "µ": "u", "μ": "u",
    "ħ": "hbar", "ε": "eps", "ν": "nu", "ρ": "rho", "σ": "sigma",
    "τ": "tau", "λ": "lambda", "γ": "gamma", "Γ": "Gamma", "α": "alpha",
    "β": "beta", "δ": "delta", "Δ": "Delta", "θ": "theta", "Θ": "Theta",
    "κ": "kappa", "π": "pi", "Φ": "Phi", "φ": "phi", "χ": "chi",
    "Ψ": "Psi", "ψ": "psi", "Σ": "Sum", "Π": "Prod", "∫": "int",
    "∂": "d", "∇": "grad", "∞": "inf", "√": "sqrt", "≈": "~=",
    "≠": "!=", "≤": "<=", "≥": ">=", "±": "+/-", "×": "x", "·": ".",
    "°": "deg", "‰": "por mil",
    # superíndices y subíndices
    "⁰": "0", "¹": "1", "²": "2", "³": "3", "⁴": "4", "⁵": "5", "⁶": "6",
    "⁷": "7", "⁸": "8", "⁹": "9", "⁻": "-", "⁺": "+",
    "₀": "0", "₁": "1", "₂": "2", "₃": "3", "₄": "4", "₅": "5", "₆": "6",
    "₇": "7", "₈": "8", "₉": "9", "ₑ": "e", "ₐ": "a",
    "Λ": "Lambda", "Ξ": "Xi", "ξ": "xi", "ζ": "zeta", "η": "eta",
    "⇌": "<=>", "⇄": "<=>", "∝": "prop a", "≡": "==", "∼": "~",
    "⟨": "<", "⟩": ">", "⊗": "(x)", "⊕": "(+)", "∈": "en", "∀": "para todo",
    "≫": ">>", "≪": "<<", "∅": "vacio", "⌀": "vacio",
    # el menos tipográfico (U+2212) NO es el guion ASCII, y es el que más se
    # cuela al copiar de un PDF o al escribir −0.5 en un editor decente
    "−": "-", "‐": "-", "‑": "-", "‒": "-", "\u00a0": " ", "\u2009": " ",
    "⁄": "/",
    # flechas y marcas
    "→": "->", "←": "<-", "↑": "^", "↓": "v", "↔": "<->", "⇒": "=>",
    "✓": "ok", "✔": "ok", "✗": "X", "✘": "X", "●": "*", "·": ".",
    "—": "--", "–": "-", "“": '"', "”": '"', "‘": "'", "’": "'",
    "«": '"', "»": '"', "…": "...",
    # cajas y viñetas
    "│": "|", "─": "-", "├": "|-", "└": "\\-", "┌": ",-", "┐": "-.",
    "┘": "-'", "┬": "-,", "┴": "-'", "┼": "-|-", "║": "|", "═": "=",
    "▪": "-", "▸": ">", "▹": ">", "•": "-", "↳": "->",
    # castellano: solo hacen falta cuando el destino es ASCII puro, pero
    # sin ellas "energía" sale como "energ\xeda" y no se puede leer
    "á": "a", "é": "e", "í": "i", "ó": "o", "ú": "u", "ü": "u", "ñ": "n",
    "Á": "A", "É": "E", "Í": "I", "Ó": "O", "Ú": "U", "Ü": "U", "Ñ": "N",
    "¿": "?", "¡": "!", "º": "o", "ª": "a",
    # números encerrados
    "①": "(1)", "②": "(2)", "③": "(3)", "④": "(4)", "⑤": "(5)",
    "⑥": "(6)", "⑦": "(7)", "⑧": "(8)", "⑨": "(9)", "⑩": "(10)",
}

_TABLA = str.maketrans({k: v for k, v in TRANSLITERACION.items()})

# se rellena en preparar(): (codificación efectiva, si hubo que transliterar)
ESTADO = {"codificacion": None, "translitera": False, "forzado_utf8": False}


def transliterar(texto: str) -> str:
    """Cambia lo que no es ASCII por su equivalente legible."""
    return str(texto).translate(_TABLA)


def _admite(flujo, muestra="Å α → ① ✓ ε²") -> bool:
    cod = getattr(flujo, "encoding", None)
    if not cod:
        return False
    try:
        muestra.encode(cod)
        return True
    except (UnicodeEncodeError, LookupError):
        return False


class _Transliterando:
    """Envoltorio de un flujo que sustituye lo que su codificación no admite.

    Solo toca lo que haga falta: si la línea es ASCII pura —la mayoría— pasa
    tal cual. Y si aun después de transliterar queda algo que no cabe, se
    escribe con `backslashreplace` en vez de reventar: se ve feo, pero el
    comando termina y el resto del informe llega.
    """

    def __init__(self, flujo):
        self._f = flujo

    def write(self, texto):
        try:
            return self._f.write(texto)
        except UnicodeEncodeError:
            pass
        limpio = transliterar(texto)
        try:
            return self._f.write(limpio)
        except UnicodeEncodeError:
            cod = getattr(self._f, "encoding", "ascii") or "ascii"
            return self._f.write(
                limpio.encode(cod, "backslashreplace").decode(cod))

    def __getattr__(self, nombre):
        return getattr(self._f, nombre)


def preparar(forzar_ascii: bool = False) -> dict:
    """Deja stdout y stderr en condiciones. Se llama una vez, al arrancar.

    Devuelve el estado para que `olla-dft sistema` pueda contarlo.
    """
    if forzar_ascii:
        sys.stdout = _Transliterando(_SoloAscii(sys.stdout))
        sys.stderr = _Transliterando(_SoloAscii(sys.stderr))
        ESTADO.update(codificacion="ascii (forzado)", translitera=True)
        return dict(ESTADO)

    for nombre in ("stdout", "stderr"):
        flujo = getattr(sys, nombre)
        if _admite(flujo):
            continue
        # 1) intentar hablar UTF-8: Windows 10+ lo admite
        try:
            flujo.reconfigure(encoding="utf-8")
            ESTADO["forzado_utf8"] = True
            continue
        except (AttributeError, ValueError, OSError):
            pass
        # 2) no se puede: transliterar antes de escribir
        setattr(sys, nombre, _Transliterando(flujo))
        ESTADO["translitera"] = True

    ESTADO["codificacion"] = getattr(sys.stdout, "encoding", None) or "?"
    return dict(ESTADO)


class _SoloAscii:
    """Flujo que finge no admitir nada fuera de ASCII. Para --ascii y pruebas."""

    encoding = "ascii"

    def __init__(self, flujo):
        self._f = flujo

    def write(self, texto):
        texto.encode("ascii")          # lanza si no es ASCII: lo que queremos
        return self._f.write(texto)

    def __getattr__(self, nombre):
        return getattr(self._f, nombre)
