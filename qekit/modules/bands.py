# Olla-DFT — command-line toolkit for Quantum ESPRESSO
# Copyright (C) 2026 Jorge Enrique González Sevilla
# SPDX-License-Identifier: GPL-3.0-or-later
# This program is free software: you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by the Free
# Software Foundation, either version 3 of the License, or (at your option)
# any later version. See the LICENSE file for details.

"""Post-proceso de estructura de bandas.

Lee el resultado de un cálculo `calculation='bands'` (o nscf) desde el XML de
pw.x, construye el camino en el espacio recíproco, analiza el band gap y
exporta los datos listos para graficar en cualquier programa (OriginPro,
gnuplot, Excel) además de producir la gráfica directamente.
"""

import re
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

from qekit.core import qeout
from qekit.core import provenance
from qekit.core import style as qstyle
from qekit.core.errors import ErrorDeUso, FaltanDatos

# Tolerancia para considerar que una banda cruza el nivel de referencia (eV)
CROSS_TOL = 1e-6


@dataclass
class BandStructure:
    result: qeout.QEResult
    kdist: np.ndarray = None          # (nk,) distancia acumulada en Å⁻¹
    labels: list = field(default_factory=list)   # [(índice, etiqueta), ...]
    breaks: list = field(default_factory=list)   # índices donde el camino salta

    @property
    def energies(self) -> np.ndarray:
        """(nspin, nk, nbnd) en eV."""
        return self.result.eigenvalues


def _build_kdist(kcart: np.ndarray, breaks: set) -> np.ndarray:
    """Distancia acumulada a lo largo del camino, con saltos de longitud cero."""
    dist = np.zeros(len(kcart))
    for i in range(1, len(kcart)):
        if i in breaks:
            step = 0.0
        else:
            step = np.linalg.norm(kcart[i] - kcart[i - 1])
        dist[i] = dist[i - 1] + step
    return dist


def _detect_breaks(matched: list) -> set:
    """Discontinuidades del camino: dos puntos especiales en índices seguidos.

    Cuando seekpath corta el camino (por ejemplo U|K), pw.x genera los dos
    puntos uno tras otro; ahí la distancia no debe acumularse.
    """
    breaks = set()
    for (i1, _), (i2, _) in zip(matched, matched[1:]):
        if i2 == i1 + 1:
            breaks.add(i2)
    return breaks


def _merge_break_labels(matched: list, breaks: set) -> list:
    """Une las etiquetas de una discontinuidad en una sola: 'U|K'."""
    merged = []
    skip = set()
    for pos, (idx, label) in enumerate(matched):
        if idx in skip:
            continue
        if pos + 1 < len(matched):
            nxt_idx, nxt_label = matched[pos + 1]
            if nxt_idx in breaks:
                merged.append((idx, f"{label}|{nxt_label}"))
                skip.add(nxt_idx)
                continue
        merged.append((idx, label))
    return merged


def load(
    path: str = ".",
    prefix: str = None,
    kpath_file: str = None,
    bands_input: str = None,
) -> BandStructure:
    """Carga una estructura de bandas desde la carpeta de un cálculo.

    Las etiquetas de alta simetría se toman de KPATH.txt (que escribe
    `olla-dft gen`); si no existe, se intenta leer la tarjeta K_POINTS
    crystal_b de bands.in. Sin ninguna de las dos, se grafica sin etiquetas.
    """
    res = qeout.read_xml(path, prefix)
    base = Path(path if Path(path).is_dir() else Path(path).parent)

    raw_labels = []
    kf = Path(kpath_file) if kpath_file else base / "KPATH.txt"
    if kf.is_file():
        raw_labels = qeout.read_kpath_labels(kf)
    else:
        bi = Path(bands_input) if bands_input else base / "bands.in"
        if bi.is_file():
            pts = qeout.read_crystal_b_card(bi)
            raw_labels = [(lab or "", c) for lab, c in pts if lab]

    matched = (
        qeout.match_labels_to_kpoints(res.kpoints_frac, raw_labels)
        if raw_labels
        else []
    )
    breaks = _detect_breaks(matched)
    kdist = _build_kdist(res.kpoints_cart, breaks)
    labels = _merge_break_labels(matched, breaks)

    return BandStructure(result=res, kdist=kdist, labels=labels, breaks=sorted(breaks))


# ----------------------------------------------------------------------
# Análisis del band gap
# ----------------------------------------------------------------------
@dataclass
class GapInfo:
    is_metal: bool = False
    vbm: float = None
    cbm: float = None
    gap: float = None                 # gap fundamental (indirecto o directo)
    direct_gap: float = None          # menor gap vertical
    vbm_kindex: int = None
    cbm_kindex: int = None
    direct_kindex: int = None
    vbm_label: str = ""
    cbm_label: str = ""
    direct_label: str = ""
    is_direct: bool = False
    vbm_band: int = None
    spin: int = 0
    fermi: float = None


def _label_for(kindex: int, labels: list, kfrac: np.ndarray) -> str:
    """Etiqueta del punto k, o sus coordenadas si no es de alta simetría."""
    for idx, lab in labels:
        if idx == kindex:
            return lab
    x, y, z = kfrac[kindex]
    return f"({x:.3f}, {y:.3f}, {z:.3f})"


def analyze_gap(bs: BandStructure, spin: int = 0) -> GapInfo:
    """Determina VBM, CBM, gap fundamental y gap directo de un canal de espín."""
    res = bs.result
    energies = res.eigenvalues[spin]          # (nk, nbnd)
    info = GapInfo(spin=spin, fermi=res.fermi)

    if energies.ndim != 2 or not energies.size or not np.isfinite(energies).all():
        raise FaltanDatos("eigenvalores vacíos o no finitos; no se puede determinar el gap")

    # Nivel de referencia para separar bandas ocupadas de vacías
    ref = res.fermi
    if ref is None:
        ref = res.homo
    if ref is None:
        # último recurso: contar electrones
        n_occ = int(round(res.nelec / (1 if res.noncolin else 2))) if res.nspin == 1 else None
        if n_occ and n_occ >= res.nbnd and res.occupations_kind == "fixed":
            info.vbm = float(energies.max())
            return info
        if n_occ and 0 < n_occ < res.nbnd:
            ref = 0.5 * (energies[:, n_occ - 1].max() + energies[:, n_occ].min())
        else:
            ref = float(np.median(energies))

    # ¿alguna banda cruza el nivel de referencia? -> metal
    band_min = energies.min(axis=0)
    band_max = energies.max(axis=0)
    crossing = (band_min < ref - CROSS_TOL) & (band_max > ref + CROSS_TOL)
    if np.any(crossing):
        info.is_metal = True
        return info

    # bandas completamente por debajo / por encima de la referencia
    below = band_max <= ref + CROSS_TOL
    above = band_min > ref - CROSS_TOL
    if np.any(below) and not np.any(above):
        # No empty bands is insufficient evidence for a gap, not a metal.
        info.vbm = float(energies[:, below].max())
        return info
    if not np.any(below):
        raise FaltanDatos("no hay bandas ocupadas para determinar el gap")

    vb_index = int(np.max(np.where(below)[0]))
    cb_index = int(np.min(np.where(above)[0]))
    if cb_index <= vb_index:
        cb_index = vb_index + 1
    if cb_index >= res.nbnd:
        # no hay bandas de conducción calculadas
        info.is_metal = False
        info.vbm = float(energies[:, vb_index].max())
        return info

    vb = energies[:, vb_index]
    cb = energies[:, cb_index]

    info.vbm_band = vb_index
    info.vbm_kindex = int(np.argmax(vb))
    info.cbm_kindex = int(np.argmin(cb))
    info.vbm = float(vb[info.vbm_kindex])
    info.cbm = float(cb[info.cbm_kindex])
    info.gap = info.cbm - info.vbm

    vertical = cb - vb
    info.direct_kindex = int(np.argmin(vertical))
    info.direct_gap = float(vertical[info.direct_kindex])
    info.is_direct = info.vbm_kindex == info.cbm_kindex

    kfrac = res.kpoints_frac
    info.vbm_label = _label_for(info.vbm_kindex, bs.labels, kfrac)
    info.cbm_label = _label_for(info.cbm_kindex, bs.labels, kfrac)
    info.direct_label = _label_for(info.direct_kindex, bs.labels, kfrac)
    return info


def gap_report(bs: BandStructure) -> str:
    """Reporte legible del análisis de gap (todos los canales de espín)."""
    res = bs.result
    lines = ["--- Análisis de band gap ---"]
    lines.append(f"Archivo XML: {res.xml_path}")
    lines.append(f"Cálculo: {res.calculation or '?'}  |  bandas: {res.nbnd}  |  "
                 f"puntos k: {res.nk}  |  electrones: {res.nelec:g}")
    if res.fermi is not None:
        lines.append(f"Energía de Fermi: {res.fermi:.4f} eV")
    if res.converged is False:
        lines.append("ADVERTENCIA: el cálculo no convergió; el gap no está validado.")
    lines.append("")

    for spin in range(res.nspin):
        if res.nspin == 2:
            lines.append(f"[Canal de espín {'up' if spin == 0 else 'down'}]")
        info = analyze_gap(bs, spin)
        if info.is_metal:
            lines.append("  Sistema METÁLICO: hay bandas que cruzan el nivel de Fermi.")
            lines.append("")
            continue
        if info.gap is None:
            lines.append("  No hay bandas de conducción en el cálculo "
                         "(aumenta nbnd para obtener el gap).")
            lines.append("")
            continue
        tipo = "DIRECTO" if info.is_direct else "INDIRECTO"
        lines.append(f"  Band gap fundamental: {info.gap:.4f} eV  ({tipo})")
        lines.append(f"    VBM = {info.vbm:.4f} eV  en  {info.vbm_label}")
        lines.append(f"    CBM = {info.cbm:.4f} eV  en  {info.cbm_label}")
        if not info.is_direct:
            lines.append(f"  Gap directo mínimo: {info.direct_gap:.4f} eV  "
                         f"en  {info.direct_label}")
        lines.append("")

    tiene_gap = any(
        not analyze_gap(bs, s).is_metal and analyze_gap(bs, s).gap is not None
        for s in range(res.nspin)
    )
    if tiene_gap:
        lines.append("Recuerda: los funcionales GGA/LDA subestiman el gap "
                     "sistemáticamente (típicamente 30–50 %).")
    return "\n".join(lines).rstrip()


# ----------------------------------------------------------------------
# Exportación de datos
# ----------------------------------------------------------------------
def reference_energy(bs: BandStructure, ref: str = "auto") -> tuple:
    """Energía de referencia para desplazar el cero.

    'auto'  -> VBM si es semiconductor/aislante, Fermi si es metal
    'fermi' -> energía de Fermi
    'vbm'   -> máximo de la banda de valencia
    'none'  -> sin desplazar
    Devuelve (valor, descripción).
    """
    res = bs.result
    if ref == "none":
        return 0.0, "sin desplazar (energías absolutas)"
    if ref == "fermi":
        if res.fermi is None:
            return 0.0, "sin desplazar (no hay energía de Fermi en el XML)"
        return res.fermi, "energía de Fermi"
    info = analyze_gap(bs, 0)
    if ref == "vbm":
        if info.vbm is None:
            return (res.fermi or 0.0), "energía de Fermi (no se determinó el VBM)"
        return info.vbm, "VBM"
    # auto
    if info.is_metal or info.vbm is None:
        return (res.fermi or 0.0), "energía de Fermi (sistema metálico)"
    return info.vbm, "VBM (sistema con gap)"


def export(bs: BandStructure, outdir: str = ".", ref: str = "auto") -> list:
    """Escribe BAND.dat y KLABELS. Devuelve la lista de archivos creados."""
    out = Path(outdir)
    out.mkdir(parents=True, exist_ok=True)
    res = bs.result
    shift, ref_desc = reference_energy(bs, ref)
    written = []

    for spin in range(res.nspin):
        energies = res.eigenvalues[spin] - shift
        suffix = "" if res.nspin == 1 else ("_up" if spin == 0 else "_dw")
        fname = out / f"BAND{suffix}.dat"
        header = [
            provenance.header("bandas", {"origen_energias": ref_desc},
                              titulo="Estructura de bandas"),
            "# Columnas: k(Ang^-1)  " +
            "  ".join(f"banda_{i + 1}" for i in range(res.nbnd)),
        ]
        table = np.column_stack([bs.kdist, energies])
        np.savetxt(fname, table, fmt="%14.8f", header="\n".join(header), comments="")
        written.append(str(fname))

    if bs.labels:
        fname = out / "KLABELS.dat"
        lines = ["# etiqueta   k(Ang^-1)   indice_k"]
        for idx, lab in bs.labels:
            lines.append(f"{lab:12s} {bs.kdist[idx]:14.8f} {idx:8d}")
        Path(fname).write_text("\n".join(lines) + "\n")
        written.append(str(fname))

    fname = out / "BAND_GAP.txt"
    Path(fname).write_text(gap_report(bs) + "\n")
    written.append(str(fname))
    return written


# ----------------------------------------------------------------------
# Gráfica
# ----------------------------------------------------------------------

def CURRENT_BG():
    """Color de fondo del tema activo (para bordes y rellenos de marcador)."""
    return qstyle.CURRENT.get("background", "#FFFFFF")


def _mark_extrema(ax, bs, info, shift, mono: bool):
    """Marca VBM y CBM. La forma distingue uno de otro incluso sin color,
    que es lo que hace legible la figura impresa en blanco y negro."""
    accent = qstyle.palette(4, mono=mono)
    vbm_color = qstyle.INK if mono else accent[2]
    cbm_color = qstyle.INK if mono else accent[3]
    ax.plot(bs.kdist[info.vbm_kindex], info.vbm - shift, marker="o",
            color=vbm_color, mfc=vbm_color, ms=4.0, mec=CURRENT_BG(), mew=0.6,
            ls="none", zorder=6)
    ax.plot(bs.kdist[info.cbm_kindex], info.cbm - shift, marker="s",
            color=cbm_color, mfc=CURRENT_BG(), ms=4.0, mec=cbm_color, mew=1.0,
            ls="none", zorder=6)

def plot(
    bs: BandStructure,
    outfile: str = "bandas",
    ref: str = "auto",
    emin: float = -6.0,
    emax: float = 6.0,
    formats="pdf,png",
    theme: str = None,
    size: str = None,
    family: str = None,
    background: str = None,
    palette=None,
    usetex: bool = None,
    width="single",
    journal: str = "generic",
    aspect: float = 0.88,
    mono: bool = False,
    title: str = None,
    gap_label: bool = False,
    mark_extrema: bool = True,
    panel: str = None,
    dpi: int = None,
    fat: np.ndarray = None,
    fat_label: str = None,
    fat_scale: float = 55.0,
) -> list:
    """Grafica la estructura de bandas al tamaño físico de una columna.

    Con `fat` (un arreglo (nspin, nk, nbnd) de pesos) se dibujan fatbands:
    las bandas quedan en gris fino y encima va un punto por cada (k, banda)
    de tamaño proporcional al peso del orbital elegido. Se ve de un vistazo
    de dónde sale cada banda, que es para lo que sirven.

    Devuelve la lista de archivos escritos (uno por formato).
    """
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError as exc:  # pragma: no cover
        raise RuntimeError(
            "matplotlib no está instalado. Instálalo con:\n"
            "  pip install matplotlib --break-system-packages"
        ) from exc

    st = qstyle.apply(theme, size=size, family=family, background=background,
                      palette=palette, usetex=usetex, mono=mono)
    res = bs.result
    shift, _ref_desc = reference_energy(bs, ref)

    fig, ax = qstyle.new_figure(width, journal, aspect)

    # Una sola serie no necesita color de identidad: se dibuja en tinta.
    # Con espín hay dos series y ahí sí entra la paleta.
    if res.nspin == 1:
        colors = [qstyle.INK]
    else:
        colors = qstyle.palette(2, mono=mono)

    segments = []
    start = 0
    for b in bs.breaks:
        segments.append((start, b))
        start = b
    segments.append((start, res.nk))

    fatcols = qstyle.palette(2, mono=mono) if fat is not None else None
    for spin in range(res.nspin):
        energies = res.eigenvalues[spin] - shift
        if fat is not None:
            # las bandas de fondo, finas y apagadas: lo que se lee es el peso
            kw = {"color": qstyle.INK_FAINT, "lw": st["line"] * 0.6}
        else:
            kw = {"color": colors[spin], "lw": st["line"]}
            if res.nspin == 2 and spin == 1:
                kw["dashes"] = qstyle.dash(1)
        for ib in range(res.nbnd):
            for s0, s1 in segments:
                if s1 - s0 < 2:
                    continue
                label = None
                if ib == 0 and s0 == 0 and res.nspin == 2 and fat is None:
                    label = qstyle.tex_safe(
                        "espín ↑" if spin == 0 else "espín ↓")
                ax.plot(bs.kdist[s0:s1], energies[s0:s1, ib], label=label, **kw)
        if fat is not None:
            w = np.clip(np.asarray(fat)[min(spin, fat.shape[0] - 1)], 0.0, None)
            nb = min(res.nbnd, w.shape[1])
            x = np.repeat(bs.kdist[:, None], nb, axis=1).ravel()
            y = energies[:, :nb].ravel()
            ss = (w[:, :nb].ravel() ** 1.0) * fat_scale
            visible = ss > 0.35
            etiqueta = None
            if spin == 0 and fat_label:
                etiqueta = qstyle.tex_safe(fat_label)
            elif spin == 1 and fat_label:
                etiqueta = qstyle.tex_safe(fat_label + " ↓")
            ax.scatter(x[visible], y[visible], s=ss[visible],
                       color=fatcols[min(spin, 1)], alpha=0.65,
                       linewidths=0, zorder=3, label=etiqueta)

    # Nivel de referencia
    ax.axhline(0.0, color=qstyle.INK_FAINT, lw=st["axis_line"], dashes=[3.5, 2.0],
               zorder=1)

    # Eje k: las marcas menores no significan nada sobre un camino de alta
    # simetría, así que se desactivan y solo quedan los puntos especiales.
    if bs.labels:
        ticks = [bs.kdist[i] for i, _ in bs.labels]
        ax.set_xticks(ticks)
        ax.set_xticklabels([qstyle.tex_safe(lab) for _, lab in bs.labels])
        ax.tick_params(axis="x", which="minor", top=False, bottom=False)
        ax.xaxis.set_minor_locator(plt.NullLocator())
        for t in ticks[1:-1]:
            ax.axvline(t, color=qstyle.GRID, lw=st["axis_line"], zorder=0)
    else:
        ax.set_xlabel(f"$k$ ({qstyle.angstrom()}$^{{-1}}$)")
    ax.set_xlim(bs.kdist[0], bs.kdist[-1])

    ax.set_ylabel(r"$E - E_\mathrm{F}$ (eV)" if _ref_desc.startswith("energía")
                  else r"$E - E_\mathrm{VBM}$ (eV)")
    ax.set_ylim(emin, emax)

    info = analyze_gap(bs, 0)
    if mark_extrema and not info.is_metal and info.gap is not None:
        _mark_extrema(ax, bs, info, shift, mono)
        if gap_label:
            tipo = "directo" if info.is_direct else "indirecto"
            ax.annotate(
                f"$E_\\mathrm{{g}}$ = {info.gap:.2f} eV ({tipo})",
                xy=(0.5, 0.02), xycoords="axes fraction",
                ha="center", va="bottom", fontsize=st["legend"],
                color=qstyle.INK,
                bbox=dict(facecolor=st["background"], alpha=0.85,
                          edgecolor="none", pad=1.0),
            )

    if res.nspin == 2 or fat_label:
        ax.legend(loc="upper right", frameon=False, fontsize=st["legend"],
                  scatterpoints=1, markerscale=1.0)
    if title:
        ax.set_title(title)
    if panel:
        qstyle.panel_label(ax, panel)

    written = qstyle.save(fig, outfile, formats, dpi=dpi,
                          modulo="bandas")
    plt.close(fig)
    return written


# ----------------------------------------------------------------------
# Fatbands: cuánto pesa cada orbital en cada banda
# ----------------------------------------------------------------------
_RE_ESTADO = re.compile(
    r"state #\s*(\d+):\s*atom\s*(\d+)\s*\(\s*([A-Za-z]{1,2})\s*\).*?"
    r"l=\s*(\d+)")
_RE_K = re.compile(r"^\s*k\s*=\s*([-\d.]+)\s+([-\d.]+)\s+([-\d.]+)")
_RE_BANDA = re.compile(r"^====\s*e\(\s*(\d+)\)\s*=\s*([-\d.]+)")
_RE_PESO = re.compile(r"([-\d.]+)\*\[#\s*(\d+)\]")

_LETRA = {0: "s", 1: "p", 2: "d", 3: "f"}


@dataclass
class Proyecciones:
    """Peso de cada estado atómico en cada banda, punto k a punto k."""

    estados: list = field(default_factory=list)   # [{atomo, elemento, l, orb}]
    pesos: np.ndarray = None                      # (nspin, nk, nbnd, nestados)
    kpuntos: np.ndarray = None                    # (nk, 3)
    nspin: int = 1
    fuente: str = ""

    @property
    def elementos(self) -> list:
        vistos = []
        for e in self.estados:
            if e["elemento"] not in vistos:
                vistos.append(e["elemento"])
        return vistos

    @property
    def etiquetas(self) -> list:
        """Combinaciones elemento-orbital disponibles, como 'Ni-d'."""
        vistos = []
        for e in self.estados:
            t = f"{e['elemento']}-{e['orb']}"
            if t not in vistos:
                vistos.append(t)
        return vistos


def leer_proyecciones(path) -> Proyecciones:
    """Lee las proyecciones banda a banda de la salida de projwfc.x.

    projwfc.x escribe la composición de CADA banda en CADA punto k en su
    salida de texto (los bloques `psi = 0.498*[#  1]+...`). Es lo que hace
    falta para las fatbands, y es distinto de los archivos .pdos, que ya
    vienen sumados sobre los puntos k y por eso solo sirven para la DOS.
    """
    p = Path(path)
    if p.is_dir():
        for nombre in ("projwfc.out", "proj.out", "projwfc_bands.out"):
            if (p / nombre).exists():
                p = p / nombre
                break
        else:
            raise FileNotFoundError(
                f"no encontré la salida de projwfc.x en '{path}'. Debe "
                f"llamarse projwfc.out y venir del cálculo de BANDAS (el mismo "
                f"camino de k), no del nscf de la DOS.")
    texto = p.read_text(errors="ignore")

    estados = []
    for m in _RE_ESTADO.finditer(texto):
        n, atomo, el, l = m.groups()
        estados.append({"n": int(n), "atomo": int(atomo),
                        "elemento": el.strip().capitalize(),
                        "l": int(l), "orb": _LETRA.get(int(l), f"l{l}")})
    if not estados:
        raise ErrorDeUso(
            f"'{p}' no tiene la lista de estados atómicos de projwfc.x. "
            f"¿Es de verdad una salida de projwfc.x?")
    nest = max(e["n"] for e in estados)

    bloques, actual, spin = [], None, 0
    nspin = 1
    for linea in texto.splitlines():
        if "SPIN UP" in linea.upper():
            spin, nspin = 0, 2
            continue
        if "SPIN DOWN" in linea.upper():
            spin, nspin = 1, 2
            continue
        mk = _RE_K.match(linea)
        if mk:
            actual = {"k": [float(x) for x in mk.groups()], "spin": spin,
                      "bandas": {}}
            bloques.append(actual)
            banda = None
            continue
        mb = _RE_BANDA.match(linea.strip())
        if mb and actual is not None:
            banda = int(mb.group(1))
            actual["bandas"][banda] = {}
            continue
        if actual is not None and banda is not None and "[#" in linea:
            for peso, idx in _RE_PESO.findall(linea):
                actual["bandas"][banda][int(idx)] = (
                    actual["bandas"][banda].get(int(idx), 0.0) + float(peso))

    if not bloques:
        raise ErrorDeUso(
            f"'{p}' no trae los bloques 'psi = ...' con la composición de cada "
            f"banda. projwfc.x solo los escribe si NO se le pasa filproj a "
            f"secas; comprueba que la salida esté completa.")

    por_spin = {}
    for b in bloques:
        por_spin.setdefault(b["spin"], []).append(b)
    nk = max(len(v) for v in por_spin.values())
    nbnd = max((max(b["bandas"]) if b["bandas"] else 0) for b in bloques)

    pesos = np.zeros((nspin, nk, nbnd, nest), dtype=np.float32)
    kpts = np.zeros((nk, 3))
    for s, lista in por_spin.items():
        for ik, b in enumerate(lista):
            if s == 0:
                kpts[ik] = b["k"]
            for ib, comp in b["bandas"].items():
                for idx, w in comp.items():
                    if 1 <= idx <= nest:
                        pesos[s, ik, ib - 1, idx - 1] = w
    return Proyecciones(estados=estados, pesos=pesos, kpuntos=kpts,
                        nspin=nspin, fuente=str(p))


def peso_de(proy: Proyecciones, selector: str) -> np.ndarray:
    """Peso total del selector en cada banda. Devuelve (nspin, nk, nbnd).

    El selector es 'Ni', 'Ni-d', 'd' o 'atomo:3'. Sumar sobre los estados
    que encajan y no normalizar: el peso que falta hasta 1 es la parte de la
    función de onda que no cae en ninguna esfera atómica, y esconderlo daría
    una imagen más limpia y menos cierta.
    """
    sel = str(selector).strip()
    if sel.lower().startswith("atomo:") or sel.lower().startswith("atom:"):
        try:
            n = int(sel.split(":", 1)[1])
        except ValueError:
            raise ErrorDeUso(f"'{selector}' debería ser atomo:N.") from None
        idx = [e["n"] - 1 for e in proy.estados if e["atomo"] == n]
        if not idx:
            raise ErrorDeUso(
                f"no hay ningún estado del átomo {n}; los átomos van de 1 a "
                f"{max(e['atomo'] for e in proy.estados)}.")
        return proy.pesos[:, :, :, idx].sum(axis=3)

    partes = sel.replace("_", "-").split("-")
    elemento = partes[0].capitalize() if partes[0] else None
    orbital = partes[1].lower() if len(partes) > 1 else None
    if elemento and elemento.lower() in _LETRA.values() and orbital is None:
        elemento, orbital = None, partes[0].lower()

    idx = [e["n"] - 1 for e in proy.estados
           if (elemento is None or e["elemento"] == elemento)
           and (orbital is None or e["orb"] == orbital)]
    if not idx:
        raise ErrorDeUso(
            f"'{selector}' no encaja con nada. Lo que hay: "
            f"{', '.join(proy.etiquetas)}.")
    return proy.pesos[:, :, :, idx].sum(axis=3)


def comprobar_compatibilidad(bs: BandStructure, proy: Proyecciones) -> None:
    """Se niega a mezclar proyecciones de un cálculo con bandas de otro.

    Es el fallo silencioso de este módulo: si projwfc.x corrió sobre el nscf
    de la DOS y las bandas vienen del camino de alta simetría, los dos
    arreglos tienen forma parecida, la gráfica sale y los colores no
    significan nada.
    """
    e = bs.energies
    if e is None:
        raise FaltanDatos("no hay bandas que decorar.")
    nk_b, nb_b = e.shape[1], e.shape[2]
    nk_p, nb_p = proy.pesos.shape[1], proy.pesos.shape[2]
    if nk_b != nk_p:
        raise ErrorDeUso(
            f"las bandas tienen {nk_b} puntos k y las proyecciones {nk_p}. "
            f"No son del mismo cálculo: projwfc.x tiene que correr sobre el "
            f"cálculo de BANDAS (el del camino de alta simetría), no sobre el "
            f"nscf de la DOS.")
    if nb_p < nb_b:
        raise ErrorDeUso(
            f"las proyecciones solo llegan a la banda {nb_p} y hay {nb_b}. "
            f"Vuelve a correr projwfc.x sobre el mismo cálculo.")


def report_fat(proy: Proyecciones, selector: str, bs: BandStructure = None,
               shift: float = 0.0) -> str:
    w = peso_de(proy, selector)
    L = [f"--- Proyección '{selector}' sobre las bandas ---",
         f"Fuente: {proy.fuente}",
         f"{len(proy.estados)} estados atómicos, "
         f"{proy.pesos.shape[1]} puntos k, {proy.pesos.shape[2]} bandas"
         + ("  (dos canales de espín)" if proy.nspin == 2 else ""),
         f"Disponibles: {', '.join(proy.etiquetas)}",
         "",
         f"Peso medio del selector sobre todas las bandas: {w.mean():.4f}",
         f"Peso máximo en una banda: {w.max():.4f}"]
    total = proy.pesos.sum(axis=3)
    perdido = 1.0 - float(total.mean())
    if perdido > 0.10:
        L.append(f"\nDe media, un {perdido * 100:.0f} % de cada función de "
                 f"onda NO cae dentro de ninguna\n  esfera atómica. En un "
                 f"material abierto o con mucho vacío es normal, pero quiere "
                 f"decir\n  que los pesos de la gráfica no suman uno.")
    if bs is not None and bs.energies is not None:
        e = bs.energies[0] - shift
        pesada = np.unravel_index(int(np.argmax(w[0])), w[0].shape)
        L.append(f"\nLa banda con más peso es la {pesada[1] + 1} en el punto k "
                 f"{pesada[0] + 1}, a {e[pesada]:+.3f} eV.")
    return "\n".join(L)
