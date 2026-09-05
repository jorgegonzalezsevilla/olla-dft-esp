# Olla-DFT — command-line toolkit for Quantum ESPRESSO
# Copyright (C) 2026 Jorge Enrique González Sevilla
# SPDX-License-Identifier: AGPL-3.0-or-later
# This program is free software: you can redistribute it and/or modify it
# under the terms of the GNU Affero General Public License as published by the Free
# Software Foundation, either version 3 of the License, or (at your option)
# any later version. See the LICENSE file for details.

"""Figura combinada de estructura de bandas y densidad de estados.

Es la figura estándar en artículos de estructura electrónica: las bandas a
la izquierda y la DOS girada a la derecha, compartiendo el eje de energía.
Los dos paneles usan forzosamente el mismo cero, que se toma del análisis
de bandas.
"""


from qekit.core import style as qstyle
from qekit.modules import bands as bands_mod
from qekit.modules import dos as dos_mod


def plot(
    bs: "bands_mod.BandStructure",
    dd: "dos_mod.DOSData",
    outfile: str = "bandas_dos",
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
    width="double",
    journal: str = "generic",
    aspect: float = 0.46,
    mono: bool = False,
    dash_mode: str = "auto",
    dos_mode: str = "orbital",
    ratio: float = 2.6,
    title: str = None,
    gap_label: bool = False,
    mark_extrema: bool = True,
    panel_labels: bool = True,
    dpi: int = None,
) -> list:
    """Bandas + DOS con eje de energía compartido, al ancho físico pedido."""
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

    # El cero lo fija el análisis de bandas y la DOS lo hereda: si cada
    # panel eligiera su propia referencia, las dos mitades de la figura no
    # coincidirían en energía.
    shift, ref_desc = bands_mod.reference_energy(bs, ref)

    figsize = qstyle.figure_size(width, journal, aspect)
    fig = plt.figure(figsize=figsize, layout="constrained")
    fig.get_layout_engine().set(w_pad=0.012, h_pad=0.012, hspace=0.0, wspace=0.0)
    gs = fig.add_gridspec(1, 2, width_ratios=[ratio, 1.0], wspace=0.03)
    axb = qstyle.finish_axes(fig.add_subplot(gs[0]))
    axd = qstyle.finish_axes(fig.add_subplot(gs[1], sharey=axb))

    # ---------------- panel (a): bandas ----------------
    if res.nspin == 1:
        band_colors = [qstyle.INK]
    else:
        band_colors = qstyle.palette(2, mono=mono)

    segments = []
    start = 0
    for b in bs.breaks:
        segments.append((start, b))
        start = b
    segments.append((start, res.nk))

    for spin in range(res.nspin):
        energies = res.eigenvalues[spin] - shift
        kw = {"color": band_colors[spin], "lw": st["line"]}
        if res.nspin == 2 and spin == 1:
            kw["dashes"] = qstyle.dash(1)
        for ib in range(res.nbnd):
            for s0, s1 in segments:
                if s1 - s0 < 2:
                    continue
                lab = None
                if ib == 0 and s0 == 0 and res.nspin == 2:
                    lab = qstyle.tex_safe("espín ↑" if spin == 0 else "espín ↓")
                axb.plot(bs.kdist[s0:s1], energies[s0:s1, ib], label=lab, **kw)

    axb.axhline(0.0, color=qstyle.INK_FAINT, lw=st["axis_line"],
                dashes=[3.5, 2.0], zorder=1)

    if bs.labels:
        ticks = [bs.kdist[i] for i, _ in bs.labels]
        axb.set_xticks(ticks)
        axb.set_xticklabels([qstyle.tex_safe(lab) for _, lab in bs.labels])
        axb.tick_params(axis="x", which="minor", top=False, bottom=False)
        axb.xaxis.set_minor_locator(plt.NullLocator())
        for t in ticks[1:-1]:
            axb.axvline(t, color=qstyle.GRID, lw=st["axis_line"], zorder=0)
    else:
        axb.set_xlabel(f"$k$ ({qstyle.angstrom()}$^{{-1}}$)")
    axb.set_xlim(bs.kdist[0], bs.kdist[-1])
    axb.set_ylim(emin, emax)
    axb.set_ylabel(r"$E - E_\mathrm{F}$ (eV)" if ref_desc.startswith("energía")
                   else r"$E - E_\mathrm{VBM}$ (eV)")

    info = bands_mod.analyze_gap(bs, 0)
    if mark_extrema and not info.is_metal and info.gap is not None:
        bands_mod._mark_extrema(axb, bs, info, shift, mono)
        if gap_label:
            tipo = "directo" if info.is_direct else "indirecto"
            axb.annotate(
                f"$E_\\mathrm{{g}}$ = {info.gap:.2f} eV ({tipo})",
                xy=(0.5, 0.02), xycoords="axes fraction", ha="center",
                va="bottom", fontsize=st["legend"], color=qstyle.INK,
                bbox=dict(facecolor=st["background"], alpha=0.85,
                          edgecolor="none", pad=1.0),
            )
    if res.nspin == 2:
        axb.legend(loc="upper right")

    # ---------------- panel (b): DOS ----------------
    dos_mod.draw(axd, dd, shift, st, mode=dos_mode, mono=mono,
                 dash_mode=dash_mode, vertical=True, emin=emin, emax=emax)
    axd.set_xlabel("DOS (estados eV$^{-1}$)")
    axd.set_ylim(emin, emax)
    axd.tick_params(labelleft=False)
    axd.xaxis.set_major_locator(plt.MaxNLocator(nbins=3, prune="lower"))

    n_series = len(dos_mod.series_list(dd, dos_mode)) + (
        1 if dd.total is not None else 0)
    if n_series >= 2:
        axd.legend(loc="upper right")

    if panel_labels:
        qstyle.panel_label(axb, "(a)")
        qstyle.panel_label(axd, "(b)")
    if title:
        fig.suptitle(title, fontsize=st["title"])

    written = qstyle.save(fig, outfile, formats, dpi=dpi,
                          modulo="bandas+DOS")
    plt.close(fig)
    return written
