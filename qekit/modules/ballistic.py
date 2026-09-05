# Olla-DFT — command-line toolkit for Quantum ESPRESSO
# Copyright (C) 2026 Jorge Enrique González Sevilla
# SPDX-License-Identifier: AGPL-3.0-or-later
# This program is free software: you can redistribute it and/or modify it
# under the terms of the GNU Affero General Public License as published by the Free
# Software Foundation, either version 3 of the License, or (at your option)
# any later version. See the LICENSE file for details.

"""Transporte balístico: conductancia con `pwcond.x`.

QUÉ ES Y EN QUÉ SE DIFERENCIA DE `olla-dft transport`
--------------------------------------------------
`olla-dft transport` calcula transporte DIFUSIVO: electrones que se dispersan
muchas veces mientras cruzan el material, descritos por la ecuación de
Boltzmann. Es lo que corresponde a un cristal macroscópico.

Esto es lo contrario. En un nanocontacto, una molécula entre dos
electrodos o un nanohilo corto, el electrón cruza SIN dispersarse. Ahí no
hay conductividad: hay CONDUCTANCIA, y la da la fórmula de Landauer:

    G = G0 * T(E_F)        con G0 = 2e^2/h = 7.748e-5 S = 1/(12.906 kOhm)

T es la probabilidad de que un electrón que entra por la izquierda salga
por la derecha, sumada sobre todos los canales abiertos. Como T <= número
de canales, la conductancia viene CUANTIZADA en escalones de G0 — y ver
esos escalones es la comprobación de que el cálculo está bien.

LA GEOMETRÍA IMPORTA MÁS QUE NADA
---------------------------------
`pwcond.x` no acepta cualquier estructura. Necesita:

- el transporte a lo largo de **z**, siempre;
- un electrodo IZQUIERDO periódico en z, con su propio cálculo scf;
- opcionalmente una región de dispersión (la molécula, el defecto);
- un electrodo DERECHO, que en Olla-DFT es siempre el MISMO que el
  izquierdo (`ikind=1`). pwcond.x admite electrodos distintos (`ikind=2`,
  con `prefixr` y `bdr`), pero Olla-DFT no lo prepara: haría falta un
  tercer scf y la comprobación de que las tres celdas empalman.

Las regiones tienen que tener la MISMA celda en el plano xy, y los límites
(`bdl`, `bds`) se dan en unidades de alat a lo largo de z. Para pwcond.x el
electrodo ocupa de z = 0 a z = bdl y la región de dispersión de z = 0 a
z = bds, cada una en SU celda: bdl y bds son las LONGITUDES de esas celdas a
lo largo de z (celldm(3) si alat es |a1|), no la altura del último átomo.
Poner mal esos límites es la causa número uno de resultados sin sentido, y
Olla-DFT los calcula de la geometría en vez de dejarlos al ojo.

LOS DOS MODOS
-------------
- `ikind=0`: solo la **estructura de bandas compleja** del electrodo. Da
  el número de canales abiertos a cada energía, que es la COTA SUPERIOR de
  la conductancia. Es barato y es lo primero que hay que mirar.
- `ikind=1`: la conductancia de verdad, con la región de dispersión en
  medio y el mismo electrodo a los dos lados. Mucho más caro.
"""

import re
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

from qekit.core import provenance
from qekit.core import style as qstyle
from qekit.core.errors import ErrorDeUso
from qekit.modules import sweep

#: Cuanto de conductancia 2e^2/h, en siemens.
G0 = 7.748091729e-5
#: Su inverso, en ohmios.
R0 = 1.0 / G0


@dataclass
class CondRun:
    energias: np.ndarray = None      # eV respecto de E_F
    transmision: np.ndarray = None   # T(E), adimensional
    canales: np.ndarray = None       # nº de canales abiertos por energía
    ikind: int = None
    G_fermi: float = None            # en unidades de G0
    avisos: list = field(default_factory=list)

    @property
    def G_siemens(self) -> float:
        return (self.G_fermi or float("nan")) * G0

    @property
    def R_ohm(self) -> float:
        g = self.G_fermi
        return R0 / g if g else float("inf")


# ----------------------------------------------------------------------
# Geometría
# ----------------------------------------------------------------------
def limites_z(atoms, margen: float = 0.0) -> tuple:
    """(z mínimo, z máximo) de los átomos en unidades de alat.

    pwcond.x pide los límites de cada región en unidades de alat a lo
    largo de z, y alat es |a1|. Calcularlos de la geometría evita el error
    más común del módulo: cortar la región por donde no toca.
    """
    celda = np.array(atoms.get_cell(), dtype=float)
    alat = float(np.linalg.norm(celda[0]))
    z = atoms.get_positions()[:, 2]
    return (float(z.min() - margen) / alat, float(z.max() + margen) / alat)


def longitud_z(atoms) -> float:
    """Longitud de la celda a lo largo de z, en unidades de alat (= |a1|).

    Es lo que pwcond.x espera en `bdl` (electrodo) y `bds` (región de
    dispersión): la frontera derecha de cada región es el final de SU celda,
    no la altura del último átomo. Con un átomo en z = 0 y otro a mitad de
    celda, la altura máxima atómica sería medio periodo y pwcond.x
    empalmaría las regiones donde no toca.
    """
    celda = np.array(atoms.get_cell(), dtype=float)
    alat = float(np.linalg.norm(celda[0]))
    return float(np.linalg.norm(celda[2])) / alat


def comprobar_geometria(atoms) -> list:
    """Lo que tiene que cumplir una estructura para pwcond.x."""
    problemas = []
    celda = np.array(atoms.get_cell(), dtype=float)
    # el eje de transporte es z y tiene que ser ortogonal al plano
    if abs(celda[2, 0]) > 1e-6 or abs(celda[2, 1]) > 1e-6:
        problemas.append(
            "el tercer vector de red no es paralelo a z. pwcond.x transporta "
            "SIEMPRE a lo\nlargo de z, y la celda tiene que estar orientada "
            "así.")
    if abs(celda[0, 2]) > 1e-6 or abs(celda[1, 2]) > 1e-6:
        problemas.append(
            "los vectores del plano tienen componente z. La celda tiene que "
            "ser tetragonal\nu ortorrómbica con z separado.")
    z = atoms.get_positions()[:, 2]
    largo = float(np.linalg.norm(celda[2]))
    if largo <= 0:
        problemas.append("la celda no tiene extensión en z.")
    elif (z.max() - z.min()) > 0.98 * largo:
        problemas.append(
            "los átomos llenan la celda entera en z. Para un electrodo "
            "periódico eso está\nbien; para una región de dispersión hace "
            "falta dejar sitio a los electrodos.")
    return problemas


def build_cond_input(prefixl: str, ikind: int = 0, prefixs: str = None,
                     prefixr: str = None, outdir: str = "./out",
                     energia0: float = 3.0, denergia: float = -0.1,
                     nenergia: int = 61, bdl: float = None,
                     bds: tuple = None, bdr: tuple = None,
                     band_file: str = "bands", tran_file: str = "trans.dat",
                     kpuntos=((0.0, 0.0, 1.0),), ewind: float = 1.0,
                     epsproj: float = 1e-3, nz1: int = 3,
                     ecut2d: float = None) -> str:
    lineas = [" &inputcond", f"    outdir='{outdir}'",
              f"    prefixl='{prefixl}'"]
    if prefixs:
        lineas.append(f"    prefixs='{prefixs}'")
    if prefixr:
        lineas.append(f"    prefixr='{prefixr}'")
    lineas += [f"    band_file='{band_file}'"]
    if ikind > 0:
        lineas.append(f"    tran_file='{tran_file}'")
    lineas += [f"    ikind={ikind}",
               f"    energy0={energia0}d0",
               f"    denergy={denergia}d0",
               f"    ewind={ewind}d0",
               f"    epsproj={epsproj:.1e}".replace("e-0", "d-0"),
               f"    nz1={nz1}"]
    if bdl is not None:
        lineas.append(f"    bdl={bdl:.6f}")
    if bds:
        lineas.append(f"    bds={bds[0]:.6f}")
    if bdr:
        lineas.append(f"    bdr={bdr[0]:.6f}")
    if ecut2d:
        lineas.append(f"    ecut2d={ecut2d}")
    lineas.append(" /")
    lineas.append(f"    {len(kpuntos)}")
    for k in kpuntos:
        lineas.append(f"    {k[0]} {k[1]} {k[2]}")
    lineas.append(f"    {nenergia}")
    return "\n".join(lineas) + "\n"


def prepare(electrodo, outdir: str = "balistico", dispersor=None,
            ikind: int = None, emin: float = -3.0, emax: float = 3.0,
            npuntos: int = 61, pseudo_dir: str = None,
            ecutwfc: float = None, ecutrho: float = None,
            kspacing: float = None, kpuntos=None,
            nz1: int = 3) -> tuple:
    """Escribe los scf de las regiones y el input de pwcond.x."""
    from qekit.modules import inputgen

    problemas = comprobar_geometria(electrodo)
    if dispersor is not None:
        problemas += comprobar_geometria(dispersor)
        c1 = np.array(electrodo.get_cell())[:2, :2]
        c2 = np.array(dispersor.get_cell())[:2, :2]
        if not np.allclose(c1, c2, atol=1e-4):
            problemas.append(
                "el electrodo y la región de dispersión NO tienen la misma "
                "celda en el plano xy.\npwcond.x empalma las dos regiones "
                "por ahí: si no coinciden, no hay empalme.")
    if problemas:
        raise ErrorDeUso("la geometría no sirve para pwcond.x:\n\n" +
                         "\n\n".join("  " + p for p in problemas))

    if ikind is None:
        ikind = 1 if dispersor is not None else 0
    if ikind == 2:
        raise ErrorDeUso(
            "ikind=2 (electrodos izquierdo y derecho DISTINTOS) no está "
            "implementado: Olla-DFT\nsolo prepara el caso de electrodos "
            "iguales (ikind=1). Para dos electrodos distintos\nhay que "
            "escribir a mano el tercer scf y 'prefixr' y 'bdr' en cond.in.")
    if ikind == 1 and dispersor is None:
        raise ErrorDeUso(
            "ikind=1 pide una región de dispersión: pásala con --scatterer, "
            "o usa ikind=0\npara ver solo las bandas complejas del electrodo.")

    out = Path(outdir); out.mkdir(parents=True, exist_ok=True)
    regiones = [("electrodo", electrodo)]
    if dispersor is not None:
        regiones.append(("dispersor", dispersor))

    prefijos = {}
    for nombre, at in regiones:
        common = sweep.prepare_common(at, pseudo_dir, ecutwfc, ecutrho,
                                      insulator=False, prefix=nombre[:6])
        grid = sweep.default_grid(at, kspacing)
        txt = inputgen.build_pw_input(
            atoms=at, pseudos=common["pseudos"], calculation="scf",
            prefix=common["prefix"], pseudo_dir=common["pseudo_dir"],
            ecutwfc=common["ecutwfc"], ecutrho=common["ecutrho"],
            kcard=f"K_POINTS automatic\n  {grid[0]} {grid[1]} {grid[2]} "
                  "0 0 0\n",
            insulator=False, degauss=common["degauss"],
            smearing=common["smearing"])
        sweep.write_input(out / f"scf_{nombre}.in", txt)
        prefijos[nombre] = common["prefix"]

    # el barrido de energia se hace de arriba abajo, como en los ejemplos
    paso = -(emax - emin) / max(npuntos - 1, 1)
    cond = build_cond_input(
        prefixl=prefijos["electrodo"], ikind=ikind,
        prefixs=prefijos.get("dispersor"),
        energia0=emax, denergia=paso, nenergia=npuntos,
        # fronteras derechas de cada región = longitud de SU celda en z
        bdl=None if dispersor is None else longitud_z(electrodo),
        bds=None if dispersor is None else (longitud_z(dispersor),),
        kpuntos=kpuntos or ((0.0, 0.0, 1.0),), nz1=nz1)
    sweep.write_input(out / "cond.in", cond)

    rep = ["--- Transporte balístico (pwcond.x) ---",
           f"Electrodo: {electrodo.get_chemical_formula()} "
           f"({len(electrodo)} átomos)"]
    if dispersor is not None:
        rep.append(f"Región de dispersión: "
                   f"{dispersor.get_chemical_formula()} "
                   f"({len(dispersor)} átomos)")
    rep += [f"Modo: ikind = {ikind}"
            + ("  (solo bandas complejas: el número de canales, que es la "
               "cota superior)" if ikind == 0 else
               "  (conductancia con región de dispersión)"),
            f"Ventana: {emin} a {emax} eV respecto de E_F, "
            f"{npuntos} puntos",
            "",
            f"Archivos en '{out.resolve()}':"]
    for nombre, _ in regiones:
        rep.append(f"  scf_{nombre}.in")
    rep += ["  cond.in", "",
            "Orden:  " + "  &&  ".join(
                f"pw.x -in scf_{n}.in" for n, _ in regiones)
            + "  &&  pwcond.x -in cond.in",
            ""]
    if ikind == 0:
        rep += ["Con ikind=0 NO sale la conductancia: sale el número de "
                "canales abiertos a\ncada energía, que es cuánto PODRÍA "
                "transmitir como mucho. Es barato y es\nlo primero que hay "
                "que mirar: si a E_F hay dos canales, la conductancia\nno "
                "puede pasar de 2 G0.", ""]
    rep += [f"La conductancia sale en unidades de G0 = 2e²/h = "
            f"{G0:.4e} S,\nque equivale a una resistencia de "
            f"{R0 / 1000:.3f} kΩ por canal perfecto.",
            "",
            "Esto es transporte BALÍSTICO: vale para un nanocontacto o una "
            "molécula entre\nelectrodos, no para un cristal macroscópico. "
            "Para eso está 'olla-dft transport'."]
    return {"prefijos": prefijos}, "\n".join(rep)


# ----------------------------------------------------------------------
# Lectura
# ----------------------------------------------------------------------


def collect(path) -> CondRun:
    p = Path(path)
    run = CondRun()

    trans = sorted(p.glob("trans*.dat")) + sorted(p.glob("*.tran"))
    if trans:
        d = np.loadtxt(trans[0], comments="#")
        if d.ndim == 1:
            d = d.reshape(1, -1)
        run.energias, run.transmision = d[:, 0], d[:, 1]
        run.ikind = 1

    salida = sorted(p.glob("cond*.out")) + sorted(p.glob("*.cond.out"))
    if salida:
        texto = salida[0].read_text(errors="ignore")
        m = re.search(r"ikind\s*=\s*(\d+)", texto)
        if m:
            run.ikind = int(m.group(1))
        if run.transmision is None:
            filas = re.findall(r"T_tot\s+(-?[\d.]+)\s+([\dEe.+-]+)", texto)
            if filas:
                d = np.array([[float(a), float(b)] for a, b in filas])
                run.energias, run.transmision = d[:, 0], d[:, 1]
        # canales abiertos: "Nchannels of the left tip = N", una vez por
        # energia y punto k. Si hay varios k por energia se toma el maximo.
        pares = re.findall(
            r"---\s+E-Ef\s*=\s*(-?[\d.]+).*?Nchannels of the left tip\s*=\s*"
            r"(\d+)", texto, re.S)
        if pares:
            d = np.array([[float(a), float(b)] for a, b in pares])
            energias = np.unique(d[:, 0])
            canales = np.array([d[d[:, 0] == e, 1].max() for e in energias])
            orden = np.argsort(energias)
            if run.energias is None:
                run.energias = energias[orden]
            run.canales = (canales[orden]
                           if len(energias) == len(run.energias) else None)

    if run.ikind is None:
        # pwcond.x no repite ikind en su salida; esta en el input, que
        # normalmente esta al lado. Sin el, el reporte no sabria si lo que
        # tiene delante es una conductancia o solo canales abiertos.
        entrada = p / "cond.in"
        if entrada.exists():
            m = re.search(r"ikind\s*=\s*(\d+)",
                          entrada.read_text(errors="ignore"))
            if m:
                run.ikind = int(m.group(1))

    if run.energias is None:
        raise ErrorDeUso(
            f"no se pudo leer ningun resultado de pwcond.x en {p}.\n"
            "Busca trans*.dat (conductancia) o la salida cond.out (bandas "
            "complejas).")

    if run.transmision is not None:
        i = int(np.argmin(np.abs(run.energias)))
        run.G_fermi = float(run.transmision[i])
    _avisar(run)
    return run


def _avisar(run: CondRun) -> None:
    if run.transmision is not None and run.canales is not None and \
            len(run.canales) == len(run.transmision):
        exceso = run.transmision - run.canales
        if np.any(exceso > 0.01):
            run.avisos.append(
                "La transmisión supera el número de canales abiertos en "
                "alguna energía.\nEso es imposible: T <= N por construcción. "
                "Revisa que los límites bdl/bds\nsepan dónde acaba cada "
                "región y que las celdas del plano coincidan.")
    if run.transmision is not None and np.any(run.transmision < -1e-6):
        run.avisos.append(
            "Hay transmisiones NEGATIVAS. Es señal de que el cálculo no "
            "convergió o de que\nla geometría de las regiones está mal "
            "cortada.")
    if run.ikind == 0:
        run.avisos.append(
            "Modo ikind=0: esto NO es la conductancia. Es el número de "
            "canales abiertos,\nque acota la conductancia por arriba. Para "
            "el valor real hace falta la región\nde dispersión y ikind=1 o 2.")


def report(run: CondRun) -> str:
    lines = ["--- Transporte balístico ---",
             f"Modo: ikind = {run.ikind}"]
    if run.energias is not None:
        lines.append(f"Ventana: {run.energias.min():.2f} a "
                     f"{run.energias.max():.2f} eV respecto de E_F "
                     f"({len(run.energias)} puntos)")

    if run.G_fermi is not None:
        lines += ["",
                  "Conductancia en el nivel de Fermi:",
                  f"  T(E_F)  = {run.G_fermi:.4f}",
                  f"  G       = {run.G_fermi:.4f} G0 = "
                  f"{run.G_siemens:.4e} S",
                  f"  R       = {run.R_ohm / 1000:.3f} kΩ"]
        cerca = round(run.G_fermi)
        if cerca >= 1 and abs(run.G_fermi - cerca) < 0.08:
            lines.append(
                f"  T está muy cerca de {cerca} entero: son {cerca} canal(es) "
                "transmitiendo casi\n  perfectamente. Eso es la "
                "cuantización de la conductancia, y verla es la\n  mejor "
                "señal de que el cálculo está bien planteado.")
        elif run.G_fermi < 0.1:
            lines.append(
                "  T casi cero: el contacto está cerrado a esa energía. Con "
                "una molécula en\n  medio es lo normal si su gap cae sobre "
                "E_F.")

    if run.canales is not None:
        lines += ["", "Canales abiertos (cota superior de T):",
                  f"  {'E-Ef (eV)':>11s} {'canales':>8s}"
                  + ("  {:>8s}".format("T") if run.transmision is not None
                     else "")]
        n = len(run.energias)
        idx = range(n) if n <= 12 else \
            [int(round(x)) for x in np.linspace(0, n - 1, 12)]
        for i in idx:
            fila = f"  {run.energias[i]:11.3f} {run.canales[i]:8.0f}"
            if run.transmision is not None and i < len(run.transmision):
                fila += f"  {run.transmision[i]:8.4f}"
            lines.append(fila)

    for a in run.avisos:
        lines += ["", a]
    lines += ["",
              f"G0 = 2e²/h = {G0:.4e} S; un canal perfecto son "
              f"{R0 / 1000:.3f} kΩ.",
              "Transporte BALÍSTICO: el electrón cruza sin dispersarse. "
              "Para un cristal\nmacroscópico, que es difusivo, está "
              "'olla-dft transport'."]
    return "\n".join(lines)


def export(run: CondRun, outdir: str = ".") -> list:
    out = Path(outdir); out.mkdir(parents=True, exist_ok=True)
    f = out / "BALISTICO.dat"
    cols, nombres = [run.energias], ["E-EF(eV)"]
    if run.transmision is not None:
        cols.append(run.transmision); nombres.append("T")
    if run.canales is not None and len(run.canales) == len(run.energias):
        cols.append(run.canales); nombres.append("canales")
    np.savetxt(f, np.column_stack(cols), fmt="%14.6f",
               header=provenance.header_plain(
                   "transporte balistico",
                   {"ikind": run.ikind,
                    "G_en_G0": None if run.G_fermi is None
                    else round(run.G_fermi, 5)},
                   titulo="Conductancia de Landauer") + "\n" +
               "  ".join(f"{n:>14s}" for n in nombres), comments="# ")
    txt = out / "BALISTICO.txt"
    txt.write_text(report(run) + "\n")
    return [str(f), str(txt)]


def plot(run: CondRun, outfile: str = "balistico", formats="pdf,png",
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

    if run.canales is not None and len(run.canales) == len(run.energias):
        ax.step(run.energias, run.canales, where="mid", lw=1.0,
                color=qstyle.INK_FAINT, dashes=[3, 2],
                label="canales abiertos")
    if run.transmision is not None:
        ax.plot(run.energias, run.transmision, lw=1.5, color=colores[0],
                label="T(E)")
    ax.axvline(0.0, color=qstyle.INK_FAINT, lw=st["axis_line"],
               dashes=[3.5, 2.0])
    ax.set_xlabel(r"$E - E_\mathrm{F}$ (eV)")
    ax.set_ylabel(r"$T$  ($G/G_0$)")
    ax.set_ylim(bottom=0)
    ax.legend(frameon=False, fontsize="small")
    return qstyle.save(fig, outfile, formats, dpi=dpi)
