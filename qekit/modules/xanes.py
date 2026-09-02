# Olla-DFT — command-line toolkit for Quantum ESPRESSO
# Copyright (C) 2026 Jorge Enrique González Sevilla
# SPDX-License-Identifier: GPL-3.0-or-later
# This program is free software: you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by the Free
# Software Foundation, either version 3 of the License, or (at your option)
# any later version. See the LICENSE file for details.

"""XANES / NEXAFS: absorción de rayos X cerca del borde, con `xspectra.x`.

QUÉ MIDE Y QUÉ CALCULA ESTO
---------------------------
Un experimento XANES arranca un electrón de un nivel de core y lo manda a
los estados vacíos. El espectro es, en esencia, la densidad de estados
desocupados PROYECTADA sobre el sitio del átomo absorbedor y filtrada por
la regla de selección dipolar: en un borde K (1s inicial) solo se ven los
estados finales de carácter p.

Por eso XANES dice cosas que ni la DOS ni la DRX dicen: es local (solo
ese átomo), es selectivo por elemento, y responde al estado de oxidación
y a la geometría de coordinación.

`xspectra.x` lo calcula con el método de Lanczos y fracciones continuas.
No necesita calcular estados vacíos explícitamente —solo la densidad de
carga del scf— y por eso es rápido comparado con lo que uno esperaría.

LAS TRES COSAS QUE HAY QUE HACER BIEN
-------------------------------------
1. **El hueco de core.** El átomo que absorbe se pone como una especie
   APARTE, con el pseudopotencial de hueco de core (`olla-dft corehole`), y
   el sistema lleva `tot_charge = +1` porque falta ese electrón. Sin esto
   se calcula el espectro del estado fundamental, que no es lo que mide
   el experimento: el hueco atrae los estados vacíos y corre el borde.

2. **La supercelda.** Con condiciones periódicas, el hueco de core ve sus
   propias imágenes. Hay que separarlas lo suficiente; Olla-DFT mide la
   distancia mínima entre imágenes del absorbedor y avisa si es corta.
   No hay un número universal, pero por debajo de ~8 Å el espectro suele
   depender del tamaño de la celda.

3. **La polarización.** `xepsilon` es la dirección del campo eléctrico.
   En un cristal anisótropo el espectro DEPENDE de ella, igual que en el
   experimento con luz polarizada. Para comparar con una muestra en polvo
   hay que promediar tres direcciones ortogonales — `--average` lo hace.

LO QUE ESTE MÓDULO NO DA
------------------------
La energía ABSOLUTA del borde. El eje que sale es relativo al nivel de
Fermi, no la energía de fotón del experimento (1839 eV para el borde K
del silicio). Para alinearlo hace falta la energía de enlace del nivel de
core, que no sale de aquí. Lo normal —y lo honesto— es alinear el borde
calculado con el experimental y comparar la FORMA, que es donde está la
información química.
"""

import re
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

from qekit.core import provenance
from qekit.core import style as qstyle
from qekit.core.errors import ErrorDeUso
from qekit.modules import sweep

#: Direcciones ortogonales para el promedio de polvo.
EJES = {"x": (1.0, 0.0, 0.0), "y": (0.0, 1.0, 0.0), "z": (0.0, 0.0, 1.0)}

#: Distancia mínima (Å) entre imágenes del absorbedor por debajo de la cual
#: el espectro suele depender del tamaño de la supercelda.
DIST_MINIMA = 8.0

#: Bordes que entiende xspectra.x (su variable `edge`). Los bordes M de
#: `atomconf.BORDES` sirven para generar el pseudo con hueco, pero
#: xspectra.x no los calcula.
BORDES_XSPECTRA = ("K", "L1", "L2", "L3", "L23")
#: Con qué --edge de 'olla-dft corehole' se genera el hueco de cada borde
#: (L2 y L3 comparten el hueco 2p).
BORDE_COREHOLE = {"K": "K", "L1": "L1", "L2": "L23", "L3": "L23",
                  "L23": "L23"}


def validar_borde(borde: str) -> str:
    """Devuelve el borde en mayúsculas o lanza ErrorDeUso si xspectra.x no
    lo admite (en particular, los bordes M)."""
    b = str(borde or "").strip().upper()
    if b in BORDES_XSPECTRA:
        return b
    if b.startswith("M"):
        raise ErrorDeUso(
            f"borde '{borde}': xspectra.x solo calcula bordes K y L "
            f"({', '.join(BORDES_XSPECTRA)}); los bordes M no están "
            "implementados en QE, aunque 'olla-dft corehole' pueda generar "
            "el pseudo con ese hueco.")
    raise ErrorDeUso(
        f"borde '{borde}' desconocido para XANES. Opciones: "
        f"{', '.join(BORDES_XSPECTRA)}.")


@dataclass
class XanesRun:
    elemento: str = ""
    borde: str = ""
    sitio: int = 0
    energias: np.ndarray = None       # eV, relativas al nivel de Fermi
    sigma: np.ndarray = None          # sección eficaz (unidades arbitrarias)
    componentes: dict = field(default_factory=dict)   # dirección -> sigma
    xgamma: float = None
    distancia_imagen: float = None
    natoms: int = 0
    avisos: list = field(default_factory=list)


# ----------------------------------------------------------------------
# Preparación
# ----------------------------------------------------------------------

def etiqueta_excitada(elemento: str) -> str:
    """Nombre de la especie excitada, dentro del límite de QE.

    Las etiquetas de especie en Quantum ESPRESSO son de TRES caracteres
    (`CHARACTER(LEN=3) :: atm`). Un 'Si_h' se trunca a 'Si_' en silencio y
    después pw.x se queja de que la especie de ATOMIC_POSITIONS no existe,
    sin decir que el problema es la longitud.
    """
    base = elemento.strip()
    return (base + "h")[:3]

def distancia_imagen_minima(atoms, indice: int) -> float:
    """Distancia del absorbedor a su imagen periódica más cercana."""
    cell = np.array(atoms.get_cell())
    mejor = np.inf
    for i in (-1, 0, 1):
        for j in (-1, 0, 1):
            for k in (-1, 0, 1):
                if (i, j, k) == (0, 0, 0):
                    continue
                mejor = min(mejor, float(np.linalg.norm(
                    i * cell[0] + j * cell[1] + k * cell[2])))
    return mejor


def build_xspectra_input(prefix: str, xiabs: int, epsilon, filecore: str,
                         r_paw: float = 3.0, borde: str = "K",
                         xemin: float = -10.0, xemax: float = 30.0,
                         xgamma: float = 0.8, xnepoint: int = 1000,
                         save_file: str = "xanes.sav",
                         xanes_file: str = "xanes.dat",
                         kgrid=(4, 4, 4), xniter: int = 2000,
                         solo_plot: bool = False) -> str:
    borde = validar_borde(borde)
    e = tuple(float(v) for v in epsilon)
    lineas = [
        " &input_xspectra",
        "    calculation='xanes_dipole',",
        f"    edge='{borde}',",
        f"    prefix='{prefix}',",
        "    outdir='./out',",
        f"    xonly_plot=.{'true' if solo_plot else 'false'}.,",
        f"    xniter={xniter},",
        "    xcheck_conv=10,",
        f"    xepsilon(1)={e[0]},",
        f"    xepsilon(2)={e[1]},",
        f"    xepsilon(3)={e[2]},",
        "    xcoordcrys=.false.,",
        f"    xiabs={xiabs},",
        f"    x_save_file='{save_file}',",
        "    xerror=0.001,",
        " /",
        " &plot",
        # xanes_file va en &plot, NO en &input_xspectra: ponerlo en el
        # primero hace que xspectra.x falle al leer el namelist entero
        # con un error que no nombra la variable culpable.
        f"    xanes_file='{xanes_file}',",
        f"    xnepoint={xnepoint},",
        f"    xgamma={xgamma},",
        f"    xemin={xemin},",
        f"    xemax={xemax},",
        "    terminator=.true.,",
        "    cut_occ_states=.true.,",
        " /",
        " &pseudos",
        f"    filecore='{filecore}',",
        f"    r_paw(1)={r_paw},",
        " /",
        " &cut_occ",
        "    cut_desmooth=0.1,",
        "    cut_stepl=0.01,",
        " /",
        f"{kgrid[0]} {kgrid[1]} {kgrid[2]} 1 1 1",
    ]
    return "\n".join(lineas) + "\n"


def prepare(atoms, elemento: str, core_hole_upf: str, outdir: str = "xanes",
            sitio: int = 0, borde: str = "K", pseudo_dir: str = None,
            ecutwfc: float = None, ecutrho: float = None,
            kspacing: float = None, insulator: bool = False,
            polarizacion=(1.0, 0.0, 0.0), promedio: bool = False,
            xemin: float = -10.0, xemax: float = 30.0, xgamma: float = 0.8,
            r_paw: float = 3.0, kgrid=None) -> tuple:
    """Escribe scf.in (con el absorbedor como especie aparte) y xspectra.in."""
    from qekit.modules import corehole, inputgen

    borde = validar_borde(borde)
    simbolos = list(atoms.get_chemical_symbols())
    if elemento not in simbolos:
        raise ErrorDeUso(
            f"'{elemento}' no está en la estructura "
            f"({', '.join(dict.fromkeys(simbolos))}).")
    indices = [i for i, s in enumerate(simbolos) if s == elemento]
    if sitio >= len(indices):
        raise ErrorDeUso(
            f"solo hay {len(indices)} átomos de {elemento}; pediste el "
            f"sitio {sitio} (se cuentan desde 0).")
    iabs = indices[sitio]

    out = Path(outdir); out.mkdir(parents=True, exist_ok=True)

    # El absorbedor pasa a ser su propia especie, y va PRIMERO para que
    # xiabs=1 y no haya que contar tipos a mano.
    marcado = atoms.copy()
    orden = [iabs] + [i for i in range(len(atoms)) if i != iabs]
    marcado = marcado[orden]
    etiqueta = etiqueta_excitada(elemento)

    common = sweep.prepare_common(
        marcado, pseudo_dir, ecutwfc, ecutrho, insulator,
        exclude_pseudos=[core_hole_upf], tarea="xanes")
    grid = kgrid or sweep.default_grid(marcado, kspacing)

    # El primer átomo se escribe con la especie excitada; el resto normal.
    scf = inputgen.build_pw_input(
        atoms=marcado, pseudos=common["pseudos"], calculation="scf",
        prefix=common["prefix"], pseudo_dir=common["pseudo_dir"],
        ecutwfc=common["ecutwfc"], ecutrho=common["ecutrho"],
        kcard=f"K_POINTS automatic\n  {grid[0]} {grid[1]} {grid[2]} 0 0 0\n",
        insulator=insulator, degauss=common["degauss"],
        smearing=common["smearing"], tot_charge=1.0)
    scf = _marcar_absorbedor(scf, elemento, etiqueta,
                             Path(core_hole_upf).name)
    sweep.write_input(out / "scf.in", scf)

    # la función de onda de core sale del propio UPF
    filecore = corehole.core_wfc(core_hole_upf, out / f"{elemento}.wfc")

    direcciones = (list(EJES.items()) if promedio
                   else [("pol", tuple(polarizacion))])
    for nombre, eps in direcciones:
        sweep.write_input(
            out / f"xspectra_{nombre}.in",
            build_xspectra_input(
                common["prefix"], xiabs=1, epsilon=eps,
                filecore=Path(filecore).name, r_paw=r_paw, borde=borde,
                xemin=xemin, xemax=xemax, xgamma=xgamma,
                save_file=f"xanes_{nombre}.sav",
                xanes_file=f"xanes_{nombre}.dat", kgrid=grid))

    dmin = distancia_imagen_minima(marcado, 0)
    rep = ["--- XANES / NEXAFS ---",
           f"Estructura: {atoms.get_chemical_formula()} ({len(atoms)} átomos)",
           f"Absorbedor: {elemento} (sitio {sitio}), borde {borde}",
           f"Especie excitada: {etiqueta} con {Path(core_hole_upf).name}",
           "Carga total: +1 (falta el electrón de core)",
           f"Distancia mínima entre imágenes del absorbedor: {dmin:.2f} Å",
           ""]
    if dmin < DIST_MINIMA:
        rep += [f"AVISO: {dmin:.1f} Å es poco. Con condiciones periódicas el "
                "hueco de core ve sus\npropias imágenes y el espectro puede "
                "depender del tamaño de la celda. Haz\nuna supercelda "
                "('olla-dft supercell') y comprueba que el espectro no cambie.",
                ""]
    rep += [f"Archivos en '{out.resolve()}':",
            "  scf.in           el scf con el hueco de core",
            f"  {Path(filecore).name:16s} función de onda de core (del UPF)"]
    for nombre, _ in direcciones:
        rep.append(f"  xspectra_{nombre}.in{'':<3s}"
                   f"{'promedio de polvo' if promedio else 'polarización'}")
    rep += ["",
            "Orden:  pw.x -in scf.in   ->   xspectra.x -in xspectra_*.in",
            "",
            "El eje de energía sale RELATIVO al nivel de Fermi, no en energía",
            "de fotón. Para comparar con un experimento se alinea el borde y",
            "se compara la FORMA."]
    warn = sweep.missing_pseudo_warning(common)
    if warn:
        rep.append(warn)
    return common, "\n".join(rep)


def _marcar_absorbedor(scf: str, elemento: str, etiqueta: str,
                       upf: str) -> str:
    """Convierte el PRIMER átomo del elemento en la especie excitada.

    Se hace sobre el texto ya generado en vez de tocar el constructor de
    inputs porque aquí la especie excitada SÍ tiene un átomo — al revés
    que en XPS, donde existe solo como tipo al que apuntar.
    """
    lineas = scf.split("\n")
    fuera, en_especies, en_pos, hecho = [], False, False, False
    for ln in lineas:
        if ln.startswith("ATOMIC_SPECIES"):
            en_especies, en_pos = True, False
            fuera.append(ln)
            continue
        if ln.startswith("ATOMIC_POSITIONS"):
            en_especies, en_pos = False, True
            fuera.append(ln)
            continue
        if ln.startswith(("CELL_PARAMETERS", "K_POINTS")):
            en_especies = en_pos = False
        if en_especies and ln.strip().startswith(elemento + " "):
            masa = ln.split()[1]
            fuera.append(f"  {etiqueta:3s} {float(masa):10.4f}  {upf}")
            fuera.append(ln)
            continue
        if en_pos and not hecho and ln.strip().startswith(elemento + " "):
            fuera.append(ln.replace(elemento, etiqueta, 1))
            hecho = True
            continue
        fuera.append(ln)
    texto = "\n".join(fuera)
    return re.sub(r"(ntyp\s*=\s*)(\d+)",
                  lambda m: m.group(1) + str(int(m.group(2)) + 1), texto,
                  count=1)


# ----------------------------------------------------------------------
# Lectura
# ----------------------------------------------------------------------


def _leer_dat(path) -> tuple:
    d = np.loadtxt(path, comments="#")
    if d.ndim == 1:
        d = d.reshape(1, -1)
    return d[:, 0], d[:, 1]


def collect(path, elemento: str = "", borde: str = "K") -> XanesRun:
    """Lee los xanes_*.dat de la carpeta y promedia si hay varios."""
    p = Path(path)
    archivos = sorted(p.glob("xanes_*.dat"))
    if not archivos:
        raise ErrorDeUso(
            f"no hay ningún xanes_*.dat en {p}. Corre primero:\n"
            "  pw.x -in scf.in  &&  xspectra.x -in xspectra_pol.in")
    run = XanesRun(elemento=elemento, borde=borde)
    for f in archivos:
        nombre = f.stem.replace("xanes_", "")
        e, s = _leer_dat(f)
        run.componentes[nombre] = s
        if run.energias is None:
            run.energias = e
        elif len(e) != len(run.energias):
            raise ErrorDeUso(
                f"{f.name} tiene {len(e)} puntos y los otros "
                f"{len(run.energias)}: no se pueden promediar. Corre las tres "
                "direcciones con los mismos xemin/xemax/xnepoint.")
    comps = list(run.componentes.values())
    run.sigma = np.mean(comps, axis=0) if len(comps) > 1 else comps[0]
    if len(comps) == 3:
        run.avisos.append(
            "Promedio de las tres direcciones ortogonales: es lo que "
            "corresponde a una muestra en polvo o a un cristal cúbico.")
    elif len(comps) == 1 and "pol" in run.componentes:
        run.avisos.append(
            "UNA sola polarización. En un cristal anisótropo el espectro "
            "depende de la dirección; para comparar con un polvo hace falta "
            "el promedio de tres direcciones (--average).")
    m = re.search(r"Broadening parameter \(in eV\):\s*([\d.]+)",
                  archivos[0].read_text(errors="ignore"))
    if m:
        run.xgamma = float(m.group(1))
    return run


def onset(run: XanesRun, fraccion: float = 0.5) -> float:
    """Energía del borde: el primer punto donde σ alcanza `fraccion` (50 %
    por omisión) de su máximo GLOBAL.

    Es una definición operativa, no la del experimento (que usa la
    derivada máxima o el primer punto de inflexión). Sirve para comparar
    espectros entre sí, no para dar una energía absoluta. Si un pre-borde
    débil precede a la línea blanca no se cuenta como borde, porque el
    umbral se mide contra el máximo de todo el espectro.
    """
    if run.sigma is None or run.energias is None:
        return float("nan")
    s = np.asarray(run.sigma)
    if s.max() <= 0:
        return float("nan")
    # primer punto que supera la fracción del máximo global
    umbral = fraccion * s.max()
    idx = np.argmax(s >= umbral)
    return float(run.energias[idx])


def report(run: XanesRun) -> str:
    lines = ["--- XANES / NEXAFS ---"]
    if run.elemento:
        lines.append(f"Absorbedor: {run.elemento}, borde {run.borde}")
    if run.energias is not None:
        lines.append(f"Rango: {run.energias[0]:.1f} a {run.energias[-1]:.1f} "
                     f"eV respecto del nivel de Fermi "
                     f"({len(run.energias)} puntos)")
    if run.xgamma is not None:
        lines.append(f"Ensanchamiento: {run.xgamma:.2f} eV")
    if run.componentes:
        lines.append(f"Polarizaciones: {', '.join(run.componentes)}")
    lines += ["", f"Borde (σ al 50 % del máximo): {onset(run):+.2f} eV "
              "respecto de E_F"]
    if run.sigma is not None:
        i = int(np.argmax(run.sigma))
        lines.append(f"Máximo principal: {run.energias[i]:+.2f} eV")
        picos = _picos(run)
        if picos:
            lines.append("Estructuras (máximos locales, en eV desde E_F):")
            for e, alto in picos[:6]:
                lines.append(f"  {e:+8.2f}   intensidad relativa {alto:.2f}")
    if len(run.componentes) == 3:
        aniso = _anisotropia(run)
        lines += ["", f"Anisotropía entre direcciones: {aniso * 100:.1f} % "
                  "del máximo"]
        if aniso > 0.1:
            lines.append(
                "  El espectro DEPENDE de la polarización: un experimento con "
                "luz\n  polarizada sobre monocristal vería espectros "
                "distintos según la orientación.")
    for a in run.avisos:
        lines += ["", a]
    lines += ["",
              "El eje es relativo al nivel de Fermi. Para comparar con un "
              "experimento hay\nque alinear el borde: la información química "
              "está en la FORMA, no en la\nposición absoluta, que este "
              "método no da."]
    return "\n".join(lines)


def _picos(run: XanesRun, prominencia: float = 0.05) -> list:
    s = np.asarray(run.sigma, dtype=float)
    e = np.asarray(run.energias, dtype=float)
    if s.size < 3 or s.max() <= 0:
        return []
    sn = s / s.max()
    idx = [i for i in range(1, len(sn) - 1)
           if sn[i] > sn[i - 1] and sn[i] >= sn[i + 1] and sn[i] > prominencia]
    return sorted([(float(e[i]), float(sn[i])) for i in idx],
                  key=lambda t: -t[1])


def _anisotropia(run: XanesRun) -> float:
    comps = np.array(list(run.componentes.values()), dtype=float)
    if comps.shape[0] < 2 or comps.max() <= 0:
        return 0.0
    return float(np.ptp(comps, axis=0).max() / comps.max())


def export(run: XanesRun, outdir: str = ".") -> list:
    out = Path(outdir); out.mkdir(parents=True, exist_ok=True)
    f = out / "XANES.dat"
    cols = [run.energias, run.sigma]
    nombres = ["E-EF(eV)", "sigma"]
    for nombre, s in run.componentes.items():
        if len(run.componentes) > 1:
            cols.append(s); nombres.append(f"sigma_{nombre}")
    np.savetxt(f, np.column_stack(cols), fmt="%14.6f",
               header=provenance.header_plain(
                   "XANES", {"elemento": run.elemento, "borde": run.borde,
                             "xgamma_eV": run.xgamma},
                   titulo="Absorcion de rayos X cerca del borde") +
               "\n" + "  ".join(f"{n:>14s}" for n in nombres),
               comments="# ")
    txt = out / "XANES.txt"
    txt.write_text(report(run) + "\n")
    return [str(f), str(txt)]


def plot(run: XanesRun, outfile: str = "xanes", formats="pdf,png",
         theme: str = None, family: str = None, background: str = None,
         palette=None, usetex: bool = None, width="single",
         journal: str = "generic", mono: bool = False,
         dpi: int = None) -> list:
    try:
        import matplotlib
        matplotlib.use("Agg")
    except ImportError as exc:                          # pragma: no cover
        raise RuntimeError("matplotlib no está instalado.") from exc

    st = qstyle.apply(theme, family=family, background=background,
                      palette=palette, usetex=usetex, mono=mono)
    fig, ax = qstyle.new_figure(width, journal, 0.75)
    colores = qstyle.palette(4, mono=mono)

    if len(run.componentes) > 1:
        for i, (nombre, s) in enumerate(run.componentes.items()):
            kw = qstyle.style_line(i + 1, 4, mono=mono)
            ax.plot(run.energias, s, lw=0.9, alpha=0.75,
                    label=rf"$\varepsilon \parallel {nombre}$", **kw)
    ax.plot(run.energias, run.sigma, lw=1.6, color=colores[0],
            label="promedio" if len(run.componentes) > 1 else "σ(E)")
    ax.axvline(0.0, color=qstyle.INK_FAINT, lw=st["axis_line"],
               dashes=[3.5, 2.0])
    ax.set_xlabel(r"$E - E_\mathrm{F}$ (eV)")
    ax.set_ylabel(r"$\sigma$ (u. arb.)")
    ax.set_xlim(run.energias[0], run.energias[-1])
    ax.set_ylim(bottom=0)
    if run.elemento:
        ax.set_title(f"XANES borde {run.borde} de {run.elemento}")
    if len(run.componentes) > 1:
        ax.legend(frameon=False)
    return qstyle.save(fig, outfile, formats, dpi=dpi)
