# Olla-DFT — command-line toolkit for Quantum ESPRESSO
# Copyright (C) 2026 Jorge Enrique González Sevilla
# SPDX-License-Identifier: AGPL-3.0-or-later
# This program is free software: you can redistribute it and/or modify it
# under the terms of the GNU Affero General Public License as published by the Free
# Software Foundation, either version 3 of the License, or (at your option)
# any later version. See the LICENSE file for details.

"""Post-proceso de densidad de estados (DOS) y DOS proyectada (PDOS).

Lee la salida de `dos.x` (archivo `<fildos>`) y la de `projwfc.x`
(archivos `<filpdos>.pdos_atm#N(El)_wfc#M(l)`), suma las proyecciones por
elemento y por orbital, y exporta tablas listas para graficar además de
producir la gráfica.
"""

import math
import re
from collections import OrderedDict
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

from qekit.core import qeout
from qekit.core import provenance
from qekit.core.errors import ErrorDeUso
from qekit.core import style as qstyle

L_OF_LETTER = {"s": 0, "p": 1, "d": 2, "f": 3}
ORBITAL_ORDER = ["s", "p", "d", "f"]

# pdos_atm#1(Si)_wfc#2(p)   /   ...#2(p_j1.5) con acoplamiento espín-órbita
_RE_PDOS = re.compile(
    r"pdos_atm#(\d+)\(([A-Za-z]+)\)_wfc#(\d+)\(([A-Za-z])[^)]*\)"
)
_RE_EFERMI = re.compile(r"EFermi\s*=\s*([-\d.eEdD+]+)")


@dataclass
class DOSData:
    energies: np.ndarray = None            # (ne,) en eV, absolutas
    total: np.ndarray = None               # (nspin, ne) DOS total
    integrated: np.ndarray = None          # (ne,) DOS integrada, si existe
    projected: "OrderedDict" = field(default_factory=OrderedDict)
    # projected: (elemento, orbital) -> (nspin, ne)
    fermi: float = None                    # eV
    nspin: int = 1
    source_files: list = field(default_factory=list)
    avisos: list = field(default_factory=list)   # se imprimen en el reporte

    @property
    def elements(self) -> list:
        seen = []
        for el, _ in self.projected:
            if el not in seen:
                seen.append(el)
        return seen

    def by_element(self) -> "OrderedDict":
        """PDOS sumada por elemento (todos sus orbitales)."""
        out = OrderedDict()
        for (el, _orb), data in self.projected.items():
            out[el] = out.get(el, 0.0) + data
        return out


# ----------------------------------------------------------------------
# Lectura
# ----------------------------------------------------------------------
def _read_table(path: Path) -> tuple:
    """Devuelve (datos, texto_del_encabezado)."""
    header = ""
    with open(path, errors="ignore") as fh:
        first = fh.readline()
        if first.lstrip().startswith("#"):
            header = first
    data = np.loadtxt(path, comments="#")
    if data.ndim == 1:
        data = data.reshape(1, -1)
    return data, header


def read_dos_file(path: str) -> tuple:
    """Lee el archivo de dos.x -> (energías, dos(nspin,ne), integrada, fermi)."""
    data, header = _read_table(Path(path))
    energies = data[:, 0]
    fermi = None
    m = _RE_EFERMI.search(header)
    if m:
        try:
            fermi = float(m.group(1).replace("D", "E"))
        except ValueError:
            fermi = None

    ncols = data.shape[1]
    if ncols >= 4:
        # E, dosup, dosdw, int  -> polarizado en espín
        total = np.vstack([data[:, 1], data[:, 2]])
        integrated = data[:, 3]
    elif ncols == 3:
        total = data[:, 1][None, :]
        integrated = data[:, 2]
    else:
        total = data[:, 1][None, :]
        integrated = None
    return energies, total, integrated, fermi


def read_pdos_file(path: Path) -> tuple:
    """Lee un archivo pdos_atm#... -> (energías, ldos(nspin,ne), nspin).

    Se usa la columna `ldos`, que ya es la suma sobre los números cuánticos
    magnéticos del orbital. Con polarización de espín hay dos columnas ldos.
    """
    data, _ = _read_table(path)
    energies = data[:, 0]
    ncols = data.shape[1]

    m = _RE_PDOS.search(path.name)
    letter = m.group(4).lower() if m else "s"
    l = L_OF_LETTER.get(letter, 0)
    nm = 2 * l + 1

    # sin espín: E + ldos + (2l+1) columnas pdos
    if ncols == 1 + 1 + nm:
        return energies, data[:, 1][None, :], 1
    # con espín: E + ldosup + ldosdw + 2*(2l+1)
    if ncols == 1 + 2 + 2 * nm:
        return energies, np.vstack([data[:, 1], data[:, 2]]), 2
    # respaldo: asumir sin espín si el número de columnas es impar
    if ncols >= 3 and (ncols - 2) % 2 == 0 and ncols > 2 + nm:
        return energies, np.vstack([data[:, 1], data[:, 2]]), 2
    return energies, data[:, 1][None, :], 1


def _interp_to(target_e: np.ndarray, source_e: np.ndarray,
               arr: np.ndarray) -> np.ndarray:
    """Interpola cada canal de espín de `arr` sobre la malla `target_e`."""
    return np.vstack([
        np.interp(target_e, source_e, arr[s], left=0.0, right=0.0)
        for s in range(arr.shape[0])
    ])


def load(
    path: str = ".",
    prefix: str = None,
    dos_file: str = None,
    pdos_prefix: str = None,
) -> DOSData:
    """Carga DOS y PDOS desde la carpeta de un cálculo.

    Busca automáticamente `<prefix>.dos` y los archivos
    `<prefix>.pdos.pdos_atm#...` que generan dos.x y projwfc.x.
    """
    base = Path(path if Path(path).is_dir() else Path(path).parent)
    result = DOSData()

    # --- energía de Fermi desde el XML (más confiable) ---
    try:
        qe = qeout.read_xml(str(base), prefix)
        result.fermi = qe.fermi
        if not prefix:
            prefix = qe.prefix
    except (FileNotFoundError, ValueError):
        qe = None

    # --- DOS total ---
    dos_path = None
    if dos_file:
        dos_path = Path(dos_file)
    else:
        candidates = sorted(base.glob(f"{prefix}.dos")) if prefix else []
        if not candidates:
            candidates = sorted(base.glob("*.dos"))
        if candidates:
            dos_path = candidates[0]

    if dos_path and dos_path.is_file():
        energies, total, integrated, fermi = read_dos_file(dos_path)
        result.energies = energies
        result.total = total
        result.integrated = integrated
        result.nspin = total.shape[0]
        if result.fermi is None:
            result.fermi = fermi
        result.source_files.append(str(dos_path))

    # --- PDOS ---
    pattern = f"{pdos_prefix}*pdos_atm#*" if pdos_prefix else "*pdos_atm#*"
    pdos_files = sorted(base.glob(pattern))
    if pdos_files:
        accum = OrderedDict()
        energies_p = None
        nspin_p = 1
        saltados = []          # (archivo, puntos) con otra malla de energía
        for f in pdos_files:
            m = _RE_PDOS.search(f.name)
            if not m:
                continue
            element = m.group(2)
            orbital = m.group(4).lower()
            e, ldos, nspin = read_pdos_file(f)
            if energies_p is None:
                energies_p = e
                nspin_p = nspin
            elif len(e) != len(energies_p):
                # Otra malla de energía: no se puede sumar punto a punto.
                # Se salta, pero se dice, porque una PDOS a la que le falta
                # un orbital entero es un resultado silenciosamente falso.
                saltados.append((f.name, len(e)))
                continue
            key = (element, orbital)
            if key in accum:
                accum[key] = accum[key] + ldos
            else:
                accum[key] = ldos
            result.source_files.append(str(f))

        # ordenar: por elemento (orden de aparición) y orbital s,p,d,f
        def sort_key(item):
            (el, orb) = item[0]
            elements = []
            for (e2, _o) in accum:
                if e2 not in elements:
                    elements.append(e2)
            return (elements.index(el), ORBITAL_ORDER.index(orb)
                    if orb in ORBITAL_ORDER else 9)

        result.projected = OrderedDict(sorted(accum.items(), key=sort_key))
        if saltados:
            lista = "\n".join(f"    {nombre}  ({n} puntos)"
                              for nombre, n in saltados)
            result.avisos.append(
                f"se han SALTADO {len(saltados)} archivo(s) de projwfc.x "
                f"cuya malla de energía no coincide\n  con la del primero "
                f"({len(energies_p)} puntos), así que la PDOS está "
                f"incompleta:\n{lista}\n  Casi siempre es que hay archivos "
                f"de dos corridas de projwfc.x mezclados en la\n  misma "
                f"carpeta (distinto Emin/Emax/DeltaE). Aparta los viejos o "
                f"vuelve a correr projwfc.x.")
        if result.energies is None:
            result.energies = energies_p
            result.nspin = nspin_p
            # sin dos.x, la DOS total es la suma de las proyecciones
            if result.projected:
                result.total = sum(result.projected.values())
        elif energies_p is not None and not np.array_equal(
            energies_p, result.energies
        ):
            # dos.x y projwfc.x pueden escribir mallas de energía distintas:
            # se interpolan las proyecciones sobre la malla de la DOS total
            result.projected = OrderedDict(
                (key, _interp_to(result.energies, energies_p, arr))
                for key, arr in result.projected.items()
            )

    if result.energies is None:
        raise FileNotFoundError(
            f"no se encontraron archivos de DOS ni PDOS en '{base}'.\n"
            "Ejecuta primero dos.x (genera <prefix>.dos) y/o projwfc.x "
            "(genera <prefix>.pdos.pdos_atm#...)."
        )
    return result


# ----------------------------------------------------------------------
# Referencia de energía y exportación
# ----------------------------------------------------------------------
def reference_energy(data: DOSData, ref: str = "auto") -> tuple:
    """Cero de energías: 'auto'/'fermi' -> Fermi, 'none' -> absolutas."""
    if ref == "none" or data.fermi is None:
        if ref != "none" and data.fermi is None:
            return 0.0, "sin desplazar (no se encontró la energía de Fermi)"
        return 0.0, "sin desplazar (energías absolutas)"
    return data.fermi, "energía de Fermi"


def export(data: DOSData, outdir: str = ".", ref: str = "auto") -> list:
    """Escribe DOS.dat y PDOS.dat. Devuelve la lista de archivos creados."""
    out = Path(outdir)
    out.mkdir(parents=True, exist_ok=True)
    shift, ref_desc = reference_energy(data, ref)
    energies = data.energies - shift
    written = []
    spin_tag = ["up", "dw"]

    # --- DOS total ---
    if data.total is not None:
        cols = [energies]
        names = ["E(eV)"]
        for s in range(data.total.shape[0]):
            cols.append(data.total[s])
            names.append("DOS" if data.total.shape[0] == 1 else f"DOS_{spin_tag[s]}")
        if data.integrated is not None:
            cols.append(data.integrated)
            names.append("DOS_integrada")
        fname = out / "DOS.dat"
        header = (
            provenance.header("DOS total",
                              {"origen_energias": ref_desc}) + "\n"
            "# " + "  ".join(f"{n:>14s}" for n in names)
        )
        np.savetxt(fname, np.column_stack(cols), fmt="%16.8f",
                   header=header, comments="")
        written.append(str(fname))

    # --- PDOS por elemento y orbital ---
    if data.projected:
        cols = [energies]
        names = ["E(eV)"]
        for (el, orb), arr in data.projected.items():
            for s in range(arr.shape[0]):
                cols.append(arr[s])
                tag = f"{el}_{orb}"
                if arr.shape[0] == 2:
                    tag += f"_{spin_tag[s]}"
                names.append(tag)
        # sumas por elemento
        for el, arr in data.by_element().items():
            for s in range(arr.shape[0]):
                cols.append(arr[s])
                tag = f"{el}_total"
                if arr.shape[0] == 2:
                    tag += f"_{spin_tag[s]}"
                names.append(tag)

        fname = out / "PDOS.dat"
        header = (
            "# Densidad de estados proyectada — Olla-DFT\n"
            f"# Origen de energías: {ref_desc} (E = 0)\n"
            "# " + "  ".join(f"{n:>14s}" for n in names)
        )
        np.savetxt(fname, np.column_stack(cols), fmt="%16.8f",
                   header=header, comments="")
        written.append(str(fname))
    return written


def report(data: DOSData, ref: str = "auto") -> str:
    """Resumen legible de lo que se cargó."""
    shift, ref_desc = reference_energy(data, ref)
    lines = ["--- DOS / PDOS ---"]
    lines.append(f"Puntos de energía: {len(data.energies)}  "
                 f"({data.energies.min():.2f} a {data.energies.max():.2f} eV)")
    if data.fermi is not None:
        lines.append(f"Energía de Fermi: {data.fermi:.4f} eV")
    lines.append(f"Origen de energías para exportar/graficar: {ref_desc}")
    lines.append(f"Canales de espín: {data.nspin}")
    if data.projected:
        lines.append("")
        lines.append("Proyecciones encontradas:")
        for el in data.elements:
            orbs = [o for (e, o) in data.projected if e == el]
            lines.append(f"  {el:3s} -> {', '.join(orbs)}")
    else:
        lines.append("")
        lines.append("Sin PDOS (no se encontraron archivos de projwfc.x).")

    # DOS en el nivel de Fermi: útil para distinguir metal de semiconductor
    if data.fermi is not None and data.total is not None:
        idx = int(np.argmin(np.abs(data.energies - data.fermi)))
        dos_ef = float(np.sum(data.total[:, idx]))
        lines.append("")
        lines.append(f"DOS en E_F: {dos_ef:.4f} estados/eV")
        if dos_ef < 1e-3:
            lines.append("  -> compatible con un sistema con gap "
                         "(semiconductor o aislante).")
        else:
            lines.append("  -> compatible con un sistema metálico.")
    for aviso in data.avisos:
        lines += ["", f"AVISO: {aviso}"]
    return "\n".join(lines)


# ----------------------------------------------------------------------
# Gráfica
# ----------------------------------------------------------------------
def series_list(data: DOSData, mode: str = "orbital") -> list:
    """[(etiqueta, arreglo), ...] según cómo se quiera descomponer la PDOS."""
    if mode == "element":
        return [(el, arr) for el, arr in data.by_element().items()]
    if mode == "orbital":
        return [(f"{el} {orb}", arr) for (el, orb), arr in data.projected.items()]
    return []


def plot(
    data: DOSData,
    outfile: str = "dos",
    ref: str = "auto",
    emin: float = -10.0,
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
    aspect: float = 0.75,
    mono: bool = False,
    dash_mode: str = "auto",
    mode: str = "orbital",
    title: str = None,
    panel: str = None,
    vertical: bool = False,
    dpi: int = None,
) -> list:
    """Grafica DOS total y proyectada al tamaño físico de una columna.

    `vertical=True` intercambia los ejes (energía en el eje y), útil cuando
    la figura acompaña a un panel de bandas.
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
    shift, _desc = reference_energy(data, ref)

    fig, ax = qstyle.new_figure(width, journal, aspect)
    draw(ax, data, shift, st, mode=mode, mono=mono, dash_mode=dash_mode,
         vertical=vertical, emin=emin, emax=emax)

    elabel = (r"$E - E_\mathrm{F}$ (eV)" if _desc.startswith("energía")
              else r"$E$ (eV)")
    dlabel = "DOS (estados eV$^{-1}$)"
    if vertical:
        ax.set_ylabel(elabel)
        ax.set_xlabel(dlabel)
        ax.set_ylim(emin, emax)
    else:
        ax.set_xlabel(elabel)
        ax.set_ylabel(dlabel)
        ax.set_xlim(emin, emax)

    n = len(series_list(data, mode)) + (1 if data.total is not None else 0)
    if n >= 2:
        ax.legend(loc="best", ncol=1 if n <= 4 else 2)
    if title:
        ax.set_title(title)
    if panel:
        qstyle.panel_label(ax, panel)

    written = qstyle.save(fig, outfile, formats, dpi=dpi,
                          modulo="DOS/PDOS")
    plt.close(fig)
    return written


def draw(ax, data: DOSData, shift: float, st: dict, mode: str = "orbital",
         mono: bool = False, dash_mode: str = "auto", vertical: bool = False,
         emin: float = None, emax: float = None):
    """Dibuja las curvas de DOS sobre unos ejes ya creados.

    Se separa de `plot` para que la figura combinada pueda reutilizarla sin
    duplicar la lógica de colores, relleno y espín.
    """
    import numpy as np

    energies = data.energies - shift
    series = series_list(data, mode)
    colors = qstyle.palette(max(len(series), 1), mono=mono)
    dashes_on = qstyle.use_dashes(len(series), dash_mode, mono)

    def _curve(values, label, color, lw, dash=None, fill=False):
        for s in range(values.shape[0]):
            sign = 1.0 if s == 0 else -1.0
            kw = {"color": color, "lw": lw}
            if dash:
                kw["dashes"] = dash
            lab = label if s == 0 else None
            if vertical:
                ax.plot(sign * values[s], energies, label=lab, **kw)
                if fill:
                    ax.fill_betweenx(energies, 0, sign * values[s],
                                     color=color, alpha=0.10, lw=0)
            else:
                ax.plot(energies, sign * values[s], label=lab, **kw)
                if fill:
                    ax.fill_between(energies, 0, sign * values[s],
                                    color=color, alpha=0.10, lw=0)

    if data.total is not None:
        _curve(data.total, "Total", qstyle.INK, st["line"], fill=True)
    for i, (label, arr) in enumerate(series):
        _curve(arr, label, colors[i], st["line"] * 0.95,
               dash=qstyle.dash(i + 1) if dashes_on else None)

    # Línea del nivel de referencia
    zero_kw = dict(color=qstyle.INK_FAINT, lw=st["axis_line"],
                   dashes=[3.5, 2.0], zorder=1)
    if vertical:
        ax.axhline(0.0, **zero_kw)
    else:
        ax.axvline(0.0, **zero_kw)

    # Con espín, el canal minoritario se dibuja reflejado
    if data.nspin == 2:
        if vertical:
            ax.axvline(0.0, color=qstyle.INK_FAINT, lw=st["axis_line"] * 0.8,
                       zorder=1)
        else:
            ax.axhline(0.0, color=qstyle.INK_FAINT, lw=st["axis_line"] * 0.8,
                       zorder=1)

    # Escala de la DOS ajustada a la ventana de energía visible
    if emin is not None and emax is not None and data.total is not None:
        mask = (energies >= emin) & (energies <= emax)
        if np.any(mask):
            dmax = float(np.max(data.total[:, mask])) * 1.12
            if dmax > 0:
                lo = -dmax if data.nspin == 2 else 0.0
                if vertical:
                    ax.set_xlim(lo, dmax)
                else:
                    ax.set_ylim(lo, dmax)
    if data.nspin == 2:
        for s, txt in ((0.965, "↑"), (0.035, "↓")):
            txt = qstyle.tex_safe(txt)
            if vertical:
                ax.text(0.965 if s > 0.5 else 0.035, 0.5, txt,
                        transform=ax.transAxes, ha="center", va="center",
                        fontsize=st["legend"], color=qstyle.INK_SOFT)
            else:
                ax.text(0.985, s, txt, transform=ax.transAxes, ha="right",
                        va="center", fontsize=st["legend"], color=qstyle.INK_SOFT)


# ----------------------------------------------------------------------
# Momentos de una banda proyectada (centro de banda d y compañía)
# ----------------------------------------------------------------------
def momentos(data: DOSData, elemento: str, orbital: str = "d",
             emax: float = None) -> dict:
    """Centro, anchura y llenado de una banda proyectada.

        ε_c = ∫ E ρ(E) dE / ∫ ρ(E) dE          (respecto al nivel de Fermi)
        W   = sqrt(∫ (E−ε_c)² ρ dE / ∫ ρ dE)

    El centro de la banda d es EL descriptor de la catálisis en metales de
    transición: cuanto más arriba está respecto al Fermi, más fuerte se
    adsorbe. Se compara directamente con la energía de adsorción que calcula
    `olla-dft adsorb`.

    Dónde se corta la integral importa y no hay un convenio único. Aquí se
    integra todo el rango disponible y se dice cuál es; si la PDOS no ha
    decaído en el extremo superior, la banda está cortada y el centro sale
    desplazado hacia abajo. El reporte lo comprueba y avisa, porque es un
    sesgo silencioso de varias décimas de eV.
    """
    orb = str(orbital).lower()
    if orb not in L_OF_LETTER:
        raise ErrorDeUso(
            f"orbital '{orbital}' desconocido; son "
            f"{', '.join(ORBITAL_ORDER)}.")
    clave = (elemento, orb)
    if clave not in data.projected:
        disp = sorted({f"{el}-{o}" for el, o in data.projected})
        raise ErrorDeUso(
            f"no hay PDOS de {elemento}-{orb} en este cálculo. "
            f"Lo que hay: {', '.join(disp) if disp else 'nada'}. "
            f"Hace falta projwfc.x con la proyección por orbital.")
    if data.fermi is None:
        raise ErrorDeUso(
            "no se encontró la energía de Fermi, y el centro de banda se mide "
            "respecto a ella. Revisa la salida de projwfc.x.")

    e = np.asarray(data.energies, dtype=float) - data.fermi
    rho = np.asarray(data.projected[clave], dtype=float)   # (nspin, ne)
    if emax is not None:
        m = e <= emax
        e, rho = e[m], rho[:, m]

    fuera = {"elemento": elemento, "orbital": orb,
             "rango": (float(e[0]), float(e[-1])), "canales": []}
    for i in range(rho.shape[0]):
        r = rho[i]
        n0 = float(np.trapezoid(r, e)) if hasattr(np, "trapezoid") \
            else float(np.trapz(r, e))
        if n0 <= 0:
            fuera["canales"].append(None)
            continue
        _int = np.trapezoid if hasattr(np, "trapezoid") else np.trapz
        centro = float(_int(e * r, e)) / n0
        ancho = math.sqrt(max(0.0, float(_int((e - centro) ** 2 * r, e)) / n0))
        ocup = e <= 0.0
        llenado = float(_int(r[ocup], e[ocup])) / n0 if ocup.any() else 0.0
        # ¿decayó la banda al final del rango? Si no, está cortada.
        cola = float(np.max(r[-max(3, len(r) // 50):]))
        pico = float(np.max(r)) or 1.0
        fuera["canales"].append(
            {"centro": centro, "ancho": ancho, "llenado": llenado,
             "estados": n0, "cola_relativa": cola / pico})
    validos = [c for c in fuera["canales"] if c]
    if validos:
        pesos = [c["estados"] for c in validos]
        total = sum(pesos)
        fuera["centro"] = sum(c["centro"] * w for c, w in zip(validos, pesos)) / total
        fuera["ancho"] = sum(c["ancho"] * w for c, w in zip(validos, pesos)) / total
        fuera["llenado"] = sum(c["llenado"] * w for c, w in zip(validos, pesos)) / total
        fuera["estados"] = total
        fuera["cola_relativa"] = max(c["cola_relativa"] for c in validos)
    return fuera


def report_momentos(m: dict) -> str:
    el, orb = m["elemento"], m["orbital"]
    if "centro" not in m:
        return f"No hay peso de {el}-{orb} que integrar."
    L = [f"--- Banda {orb} de {el} ---",
         f"Rango integrado: {m['rango'][0]:+.2f} a {m['rango'][1]:+.2f} eV "
         f"respecto al Fermi",
         "",
         f"  centro   ε_{orb} = {m['centro']:+.4f} eV",
         f"  anchura  W    = {m['ancho']:.4f} eV",
         f"  llenado       = {m['llenado'] * 100:.1f} %  "
         f"({m['estados']:.2f} estados en total)"]
    canales = [c for c in m["canales"] if c]
    if len(canales) == 2:
        L += ["", "  Por canal de espín:",
              f"    up  ε = {canales[0]['centro']:+.4f} eV   "
              f"W = {canales[0]['ancho']:.4f}   "
              f"llenado {canales[0]['llenado'] * 100:.1f} %",
              f"    down ε = {canales[1]['centro']:+.4f} eV   "
              f"W = {canales[1]['ancho']:.4f}   "
              f"llenado {canales[1]['llenado'] * 100:.1f} %",
              f"    desdoblamiento de intercambio = "
              f"{canales[0]['centro'] - canales[1]['centro']:+.4f} eV"]
    if m.get("cola_relativa", 0) > 0.05:
        L += ["",
              f"AVISO: al final del rango todavía queda un "
              f"{m['cola_relativa'] * 100:.0f} % del pico de PDOS. La banda "
              f"está\n  CORTADA por arriba, así que el centro sale más bajo de "
              f"lo que debería.\n  Vuelve a correr projwfc.x con un Emax mayor "
              f"(o nscf con más bandas)."]
    if orb == "d":
        L += ["",
              "  El centro de banda d se compara con la energía de adsorción: "
              "cuanto más\n  arriba (menos negativo), más fuerte adsorbe la "
              "superficie. Es una\n  correlación empírica dentro de una misma "
              "familia de metales, no una ley."]
    return "\n".join(L)
