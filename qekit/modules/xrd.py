# Olla-DFT — command-line toolkit for Quantum ESPRESSO
# Copyright (C) 2026 Jorge Enrique González Sevilla
# SPDX-License-Identifier: GPL-3.0-or-later
# This program is free software: you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by the Free
# Software Foundation, either version 3 of the License, or (at your option)
# any later version. See the LICENSE file for details.

"""Difracción de polvos simulada a partir de la estructura cristalina.

Para materiales laminares es la conexión directa entre el cálculo y el
experimento: la posición del pico basal mide el espaciado interlaminar, y
comparar el patrón simulado de la estructura relajada contra el
difractograma medido dice si el modelo estructural es el correcto.

Física implementada:
- posiciones por la ley de Bragg sobre todos los hkl alcanzables;
- intensidades |F(hkl)|² con factores de dispersión atómica analíticos
  (f(s) = Z − 41.78214·s²·Σ aᵢ·e^(−bᵢ·s²), s = sinθ/λ; coeficientes del
  archivo de datos tomado de pymatgen, licencia MIT — mismos valores de las
  International Tables), corregidas por Lorentz–polarización
  (1+cos²2θ)/(sin²θ·cosθ);
- la multiplicidad sale sola de enumerar todos los hkl y agrupar por 2θ;
- perfil pseudo-Voigt, con anchura fija o dada por tamaño de cristalito
  (Scherrer), que es lo que ensancha los picos de un laminar real;
- factor de temperatura global opcional (Debye–Waller isotrópico).

Verificado numéricamente contra la implementación independiente de
pymatgen sobre Si y grafito (posiciones e intensidades relativas).
"""

import json
from dataclasses import dataclass, field
from importlib import resources
from pathlib import Path

import numpy as np
from qekit.core import provenance
from ase import Atoms

from qekit.core import style as qstyle
from qekit.core.errors import ErrorDeUso

WAVELENGTHS = {
    "CuKa": 1.54184, "CuKa1": 1.54056, "CuKa2": 1.54439,
    "MoKa": 0.71073, "MoKa1": 0.70930,
    "CoKa": 1.79026, "CoKa1": 1.78896,
    "FeKa": 1.93735, "CrKa": 2.29100, "AgKa": 0.56087,
}
SCHERRER_K = 0.9


def _load_scattering():
    with resources.files("qekit.data").joinpath(
        "atomic_scattering_params.json"
    ).open() as fh:
        return json.load(fh)


_SCATTERING = None


def scattering_params(symbol: str):
    global _SCATTERING
    if _SCATTERING is None:
        _SCATTERING = _load_scattering()
    if symbol not in _SCATTERING:
        raise ErrorDeUso(
            f"no hay factores de dispersión para '{symbol}' en la tabla."
        )
    return np.array(_SCATTERING[symbol])


def wavelength_name(w) -> str:
    """Nombre legible de la radiación: 'CuKa' -> 'Cu Kα', 'MoKa1' -> 'Mo Kα1'.

    Si lo que llega es una longitud de onda numérica (o una clave que no
    está en la tabla) se devuelve 'λ dada', para que ningún rótulo diga
    'Cu Kα' cuando el usuario pidió otra cosa.
    """
    if isinstance(w, (int, float)):
        return "λ dada"
    key = str(w).strip()
    if key in WAVELENGTHS:
        elemento, resto = key[:-3] if key[-1].isdigit() else key[:-2], \
            key[-3:] if key[-1].isdigit() else key[-2:]
        letra = "α" if resto[1] == "a" else "β"
        return f"{elemento} K{letra}{resto[2:]}"
    return "λ dada"


def wavelength_value(w) -> float:
    if isinstance(w, (int, float)):
        return float(w)
    key = str(w).strip()
    if key in WAVELENGTHS:
        return WAVELENGTHS[key]
    try:
        return float(key)
    except ValueError:
        raise ErrorDeUso(
            f"radiación desconocida '{w}'. Opciones: "
            f"{', '.join(WAVELENGTHS)}, o la longitud de onda en Å."
        )


def _friedel(hkl):
    """Orienta (hkl) con su compañero de Friedel: primer índice no nulo > 0."""
    for x in hkl:
        if x != 0:
            return hkl if x > 0 else tuple(-v for v in hkl)
    return hkl


@dataclass
class Peak:
    two_theta: float
    d: float
    intensity: float          # normalizada a 100
    hkls: list = field(default_factory=list)   # [(h,k,l), ...] representativos

    @property
    def label(self) -> str:
        # Elegir el representante legible NO es tomar valores absolutos: en
        # una red hexagonal (1,-1,0) es equivalente a (100) pero (110) es
        # otra reflexión distinta. De los índices fusionados en el pico
        # (todos con el mismo 2θ), se prefiere el que no tiene componentes
        # negativas; si ninguno, el de menos negativas.
        best = min(self.hkls, key=lambda t: (sum(1 for v in t if v < 0),
                                             [-v for v in t]))
        def fmt(x):
            return str(x) if x >= 0 else f"-{abs(x)}"
        h, k, l = best
        return f"({fmt(h)}{fmt(k)}{fmt(l)})"


@dataclass
class Pattern:
    wavelength: float
    two_theta: np.ndarray = None       # perfil continuo
    intensity: np.ndarray = None
    peaks: list = field(default_factory=list)
    fwhm: float = None                 # anchura usada (°) en el punto medio
    size_nm: float = None
    basis: str = "conventional"        # celda en que se indexaron los hkl
    formula: str = ""


# ----------------------------------------------------------------------
# Cálculo del patrón
# ----------------------------------------------------------------------
def compute(atoms: Atoms, wavelength="CuKa", two_theta_range=(5.0, 70.0),
            b_iso: float = 0.0, min_intensity: float = 0.1,
            merge_tol: float = 0.02, basis: str = "conventional") -> Pattern:
    """Patrón de picos (posiciones e intensidades integradas).

    `basis` decide en qué celda se indexan los hkl:

    - "conventional" (por defecto): celda convencional estandarizada. Es la
      base de las fichas PDF y de la literatura — el (220) del silicio se
      llama (220). Las posiciones e intensidades no cambian: es el mismo
      cristal, y las extinciones salen solas del factor de estructura.
    - "input": la celda tal como viene del archivo. Si esa celda es la
      primitiva, los índices NO coincidirán con los de la ficha PDF (el
      (220) cúbico se vuelve (211) en la base primitiva FCC), así que se
      usa solo cuando se quiere indexar en esa base concreta.
    """
    lam = wavelength_value(wavelength)
    basis = (basis or "conventional").lower()
    if basis not in ("conventional", "input"):
        raise ErrorDeUso("basis debe ser 'conventional' o 'input'")
    if basis == "conventional":
        from qekit.core import structure as _struct
        try:
            atoms = _struct.conventional(atoms)
        except Exception:
            basis = "input"       # sin simetría utilizable, se deja igual
    tt_min, tt_max = two_theta_range
    if tt_min <= 0:
        tt_min = 0.5

    cell = atoms.cell.array
    recip = np.linalg.inv(cell).T          # filas b_i SIN 2π (|g| = 1/d)
    frac = atoms.get_scaled_positions()
    symbols = atoms.get_chemical_symbols()
    coeffs = {s: scattering_params(s) for s in set(symbols)}
    from ase.data import atomic_numbers
    zs = {s: atomic_numbers[s] for s in set(symbols)}

    g_max = 2.0 * np.sin(np.radians(tt_max / 2.0)) / lam
    g_min = 2.0 * np.sin(np.radians(tt_min / 2.0)) / lam

    # rangos de hkl que caben en la esfera |g| <= g_max
    lengths = np.linalg.norm(recip, axis=1)
    hmax = np.maximum(1, np.ceil(g_max / lengths).astype(int) + 1)
    hs = np.arange(-hmax[0], hmax[0] + 1)
    ks = np.arange(-hmax[1], hmax[1] + 1)
    ls = np.arange(-hmax[2], hmax[2] + 1)
    H, K, L = np.meshgrid(hs, ks, ls, indexing="ij")
    hkl = np.stack([H.ravel(), K.ravel(), L.ravel()], axis=1)
    hkl = hkl[np.any(hkl != 0, axis=1)]

    g = hkl @ recip
    gn = np.linalg.norm(g, axis=1)
    sel = (gn >= max(g_min, 1e-8)) & (gn <= g_max)
    hkl, gn = hkl[sel], gn[sel]

    s2 = (gn / 2.0) ** 2                     # (sinθ/λ)²
    theta = np.arcsin(lam * gn / 2.0)
    two_theta = np.degrees(2.0 * theta)

    # factor de estructura
    phases = np.exp(2j * np.pi * (hkl @ frac.T))     # (nrefl, nat)
    F = np.zeros(len(hkl), dtype=complex)
    for sym in set(symbols):
        mask = np.array([sy == sym for sy in symbols])
        a_b = coeffs[sym]
        fs = zs[sym] - 41.78214 * s2 * np.sum(
            a_b[:, 0] * np.exp(-a_b[:, 1] * s2[:, None]), axis=1
        )
        if b_iso:
            fs = fs * np.exp(-b_iso * s2)
        F += fs * phases[:, mask].sum(axis=1)

    lp = (1.0 + np.cos(2.0 * theta) ** 2) / (np.sin(theta) ** 2 * np.cos(theta))
    I = np.abs(F) ** 2 * lp

    # descartar las reflexiones extinguidas por simetría (F = 0): entran en
    # la fusión con peso cero y, si dos coinciden en 2θ, el promedio pesado
    # se vuelve 0/0. Además no son picos: no deben aparecer en la tabla.
    imax_raw = float(I.max()) if len(I) else 0.0
    if imax_raw > 0:
        keep = I > 1e-8 * imax_raw
        hkl, gn, two_theta, I = hkl[keep], gn[keep], two_theta[keep], I[keep]

    # agrupar reflexiones por 2θ (la multiplicidad queda incluida)
    order = np.argsort(two_theta)
    peaks = []
    for idx in order:
        tt = float(two_theta[idx])
        if peaks and abs(tt - peaks[-1].two_theta) < merge_tol:
            p = peaks[-1]
            wsum = p.intensity + I[idx]
            if wsum > 0:
                p.two_theta = (p.two_theta * p.intensity
                               + tt * I[idx]) / wsum
            p.intensity = wsum
            cand = _friedel(tuple(int(x) for x in hkl[idx]))
            if len(p.hkls) < 24 and cand not in p.hkls:
                p.hkls.append(cand)
        else:
            peaks.append(Peak(two_theta=tt, d=float(1.0 / gn[idx]),
                              intensity=float(I[idx]),
                              hkls=[_friedel(tuple(int(x) for x in hkl[idx]))]))

    if not peaks:
        return Pattern(wavelength=lam, peaks=[], basis=basis,
                       formula=atoms.get_chemical_formula())
    imax = max(p.intensity for p in peaks)
    for p in peaks:
        p.intensity *= 100.0 / imax
    peaks = [p for p in peaks if p.intensity >= min_intensity]
    return Pattern(wavelength=lam, peaks=peaks, basis=basis,
                   formula=atoms.get_chemical_formula())


def broaden(pattern: Pattern, two_theta_range=(5.0, 70.0), step: float = 0.02,
            fwhm: float = 0.15, size_nm: float = None, eta: float = 0.5) -> Pattern:
    """Convierte los picos en un perfil continuo pseudo-Voigt.

    Con `size_nm`, la anchura de cada pico sale de la ecuación de Scherrer
    (K = 0.9), que es el mecanismo dominante en cristalitos pequeños; si no,
    se usa una anchura instrumental constante `fwhm` (en grados 2θ).
    """
    lam = pattern.wavelength
    x = np.arange(two_theta_range[0], two_theta_range[1] + step, step)
    y = np.zeros_like(x)
    for p in pattern.peaks:
        if size_nm:
            theta = np.radians(p.two_theta / 2.0)
            w = np.degrees(SCHERRER_K * lam / (size_nm * 10.0 * np.cos(theta)))
        else:
            w = fwhm
        sigma = w / (2.0 * np.sqrt(2.0 * np.log(2.0)))
        gauss = np.exp(-0.5 * ((x - p.two_theta) / sigma) ** 2)
        lorentz = 1.0 / (1.0 + ((x - p.two_theta) / (w / 2.0)) ** 2)
        y += p.intensity * ((1.0 - eta) * gauss + eta * lorentz)
    if y.max() > 0:
        y *= 100.0 / y.max()
    pattern.two_theta, pattern.intensity = x, y
    pattern.fwhm, pattern.size_nm = fwhm, size_nm
    return pattern


# ----------------------------------------------------------------------
# Datos experimentales para comparar
# ----------------------------------------------------------------------
def read_experimental(path: str):
    """Lee un difractograma medido: dos columnas (2θ, intensidad).

    Acepta .xy, .dat, .txt y .csv con comentarios (#) y separadores por
    espacio, tabulador o coma. La intensidad se normaliza a 100.
    """
    raw = Path(path).read_text(errors="ignore")
    rows = []
    for line in raw.splitlines():
        line = line.strip().replace(",", " ")
        if not line or line[0] in "#!;/'":
            continue
        parts = line.split()
        if len(parts) < 2:
            continue
        try:
            rows.append((float(parts[0]), float(parts[1])))
        except ValueError:
            continue
    if len(rows) < 10:
        raise ErrorDeUso(
            f"'{path}' no parece un difractograma de dos columnas (2θ, I)."
        )
    arr = np.array(rows)
    arr = arr[np.argsort(arr[:, 0])]
    y = arr[:, 1] - arr[:, 1].min()
    if y.max() > 0:
        y *= 100.0 / y.max()
    return arr[:, 0], y


# ----------------------------------------------------------------------
# Reporte, exportación y gráfica
# ----------------------------------------------------------------------
def report(pattern: Pattern, top: int = 12) -> str:
    basis_txt = ("celda convencional (índices como en las fichas PDF)"
                 if pattern.basis == "conventional"
                 else "celda de entrada (¡los hkl NO son los de la ficha "
                      "PDF si es primitiva!)")
    lines = ["--- Difracción de polvos simulada ---",
             f"λ = {pattern.wavelength:.5f} Å  |  "
             f"{pattern.formula}  |  indexado en la {basis_txt}",
             "",
             f"{'2θ (°)':>9s} {'d (Å)':>9s} {'I rel':>7s}   hkl"]
    for p in sorted(pattern.peaks, key=lambda q: -q.intensity)[:top]:
        lines.append(f"{p.two_theta:9.3f} {p.d:9.4f} {p.intensity:7.1f}   "
                     f"{p.label}")
    lines.append("")
    lines.append(f"{len(pattern.peaks)} reflexiones en total "
                 "(tabla completa en XRD_HKL.dat).")
    return "\n".join(lines)


def export(pattern: Pattern, outdir: str = ".") -> list:
    out = Path(outdir); out.mkdir(parents=True, exist_ok=True)
    written = []
    if pattern.two_theta is not None:
        f = out / "XRD.dat"
        np.savetxt(f, np.column_stack([pattern.two_theta, pattern.intensity]),
                   fmt="%12.4f",
                   header=provenance.header_plain(
                       "difracción de polvos",
                       {"lambda_A": f"{pattern.wavelength:.5f}",
                        "base_hkl": pattern.basis,
                        "formula": pattern.formula},
                       titulo="Difractograma simulado")
                   + f"\n{'2theta':>10s} {'I':>12s}", comments="# ")
        written.append(str(f))
    f = out / "XRD_HKL.dat"
    lines = [provenance.header(
                 "reflexiones",
                 {"lambda_A": f"{pattern.wavelength:.5f}",
                  "base_hkl": pattern.basis, "formula": pattern.formula},
                 titulo="Reflexiones"),
             f"# {'2theta':>9s} {'d(A)':>10s} {'I':>9s}  hkl"]
    for p in sorted(pattern.peaks, key=lambda q: q.two_theta):
        lines.append(f"{p.two_theta:11.4f} {p.d:10.5f} {p.intensity:9.2f}  "
                     f"{p.label}")
    f.write_text("\n".join(lines) + "\n")
    written.append(str(f))
    return written


def plot(pattern: Pattern, outfile: str = "xrd", exp=None,
         label_top: int = 8, formats="pdf,png", theme: str = None,
         size: str = None, family: str = None, background: str = None,
         palette=None, usetex: bool = None, width="single",
         journal: str = "generic", aspect: float = 0.62, mono: bool = False,
         dpi: int = None, exp_label: str = "experimental") -> list:
    """Difractograma con etiquetas hkl; con `exp`, superpone el medido."""
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError as exc:  # pragma: no cover
        raise RuntimeError("matplotlib no está instalado.") from exc

    st = qstyle.apply(theme, size=size, family=family, background=background,
                      palette=palette, usetex=usetex, mono=mono)
    fig, ax = qstyle.new_figure(width, journal, aspect)
    colors = qstyle.palette(2, mono=mono)

    if exp is not None:
        xe, ye = exp
        ax.plot(xe, ye + 105.0, color=colors[1], lw=st["line"] * 0.9,
                label=qstyle.tex_safe(exp_label))
        ax.plot(pattern.two_theta, pattern.intensity, color=colors[0],
                lw=st["line"], label="simulado")
        ax.set_ylim(0, 215)
        ax.legend(loc="upper right")
        ax.set_yticks([])
    else:
        ax.plot(pattern.two_theta, pattern.intensity, color=colors[0],
                lw=st["line"])
        ax.set_ylim(0, 118)

    # etiquetas hkl sobre los picos más intensos
    labeled = sorted(pattern.peaks, key=lambda p: -p.intensity)[:label_top]
    x0, x1 = pattern.two_theta.min(), pattern.two_theta.max()
    labeled = [p for p in labeled if x0 + 0.8 < p.two_theta < x1 - 0.8]
    for p in labeled:
        ax.annotate(
            qstyle.tex_safe(p.label), xy=(p.two_theta, min(p.intensity + 2, 113)),
            ha="center", va="bottom", fontsize=st["legend"] * 0.85,
            color=qstyle.INK_SOFT, rotation=90,
        )

    ax.set_xlim(pattern.two_theta.min(), pattern.two_theta.max())
    ax.set_xlabel(r"$2\theta$ (°)" if not qstyle.USETEX else r"$2\theta$ (grados)")
    ax.set_ylabel("intensidad (u. arb.)")
    ax.tick_params(axis="y", which="both", left=False, right=False,
                   labelleft=False)
    written = qstyle.save(fig, outfile, formats, dpi=dpi,
                          modulo="difracción de polvos")
    plt.close(fig)
    return written
