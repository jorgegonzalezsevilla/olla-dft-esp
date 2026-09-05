# Olla-DFT — command-line toolkit for Quantum ESPRESSO
# Copyright (C) 2026 Jorge Enrique González Sevilla
# SPDX-License-Identifier: AGPL-3.0-or-later
# This program is free software: you can redistribute it and/or modify it
# under the terms of the GNU Affero General Public License as published by the Free
# Software Foundation, either version 3 of the License, or (at your option)
# any later version. See the LICENSE file for details.

"""U de Hubbard por respuesta lineal, con `hp.x`.

EL PROBLEMA QUE RESUELVE
------------------------
DFT+U necesita un número, U, que corrige la autointeracción de los
electrones localizados (los d de un metal de transición, los f de una
tierra rara). Sin él, un óxido como el NiO sale metálico cuando en el
laboratorio es un aislante de varios eV.

Lo que hace casi todo el mundo es copiar el U de un artículo. Eso tiene
dos problemas: el U depende del SISTEMA (no es una propiedad del
elemento) y depende del ESQUEMA DE PROYECCIÓN con el que se aplique. Un
U de 4.6 eV medido con proyectores 'atomic' no es el mismo número que
uno con 'ortho-atomic'. Copiar entre esquemas distintos es un error
silencioso.

`hp.x` lo calcula para tu sistema por respuesta lineal:

    U = (chi0^-1 - chi^-1)_II

donde chi0 es la respuesta de las ocupaciones a una perturbación del
potencial SIN dejar que el sistema se reajuste, y chi la respuesta
completa. La diferencia es la parte "espuria" que U tiene que cancelar.

LA AUTOCONSISTENCIA
-------------------
Aquí está la trampa que hace que mucha gente reporte un U mal calculado.
El U que sale DEPENDE del U que usaste en el scf de partida. Un cálculo
con U = 0 da un primer U que NO es el autoconsistente: hay que volver a
correr el scf con ese U, recalcular, y repetir hasta que deje de moverse.

Olla-DFT hace ese ciclo (`ciclo()`), reporta cuántas vueltas hizo y cuánto
se movió el número. Si alguien solo quiere una vuelta, se lo permite —
pero el reporte dice claramente que eso es un U DE PRIMERA ITERACIÓN.

LO QUE HAY QUE SABER ANTES DE CREÉRSELO
---------------------------------------
- La malla de q (`--qgrid`) es una supercelda encubierta: nq=2x2x2 son
  8 celdas. Con nq=1x1x1 la perturbación ve sus imágenes y el U sale
  mal. Hay que converger esa malla como cualquier otra.
- El `U_projection_type` del scf tiene que ser el MISMO con el que
  después uses el U. Olla-DFT lo fija en 'ortho-atomic' y lo escribe en la
  procedencia.
- El scf de partida necesita `conv_thr` muy apretado (1e-15): la
  respuesta lineal deriva ocupaciones, y derivar amplifica el ruido.
- `hp.x` arranca con un U diminuto (1e-8) en vez de cero: es como se
  activa la maquinaria de proyección sin cambiar la física.
"""

import re
import shutil
import subprocess
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

from qekit.core import provenance, runner as run_mod
from qekit.core.errors import ErrorDeUso
from qekit.modules import sweep

#: U inicial (eV) para activar la proyección sin cambiar la física.
U_SEMILLA = 1.0e-8

#: Elementos donde DFT+U es habitual, con el orbital que se corrige.
ORBITAL_HUBBARD = {
    **{el: "3d" for el in "Sc Ti V Cr Mn Fe Co Ni Cu Zn".split()},
    **{el: "4d" for el in "Y Zr Nb Mo Tc Ru Rh Pd Ag Cd".split()},
    **{el: "5d" for el in "Hf Ta W Re Os Ir Pt Au Hg".split()},
    **{el: "4f" for el in ("La Ce Pr Nd Pm Sm Eu Gd Tb Dy Ho Er Tm Yb "
                           "Lu").split()},
}


@dataclass
class SitioU:
    sitio: int = 0
    tipo: int = 0
    etiqueta: str = ""
    espin: int = 1
    U: float = None            # eV


@dataclass
class ParV:
    """Un parámetro V entre dos sitios vecinos."""

    i: int = 0                 # índice del primer átomo (numeración de hp.x)
    el_i: str = ""
    j: int = 0                 # índice del segundo, ya en la SUPERCELDA
    el_j: str = ""
    distancia: float = None    # bohr
    V: float = None            # eV

    @property
    def es_sitio(self) -> bool:
        """V(i,i) es en realidad la U del sitio."""
        return self.i == self.j


@dataclass
class HubbardRun:
    sitios: list = field(default_factory=list)
    qgrid: tuple = None
    proyeccion: str = "ortho-atomic"
    U_entrada: dict = field(default_factory=dict)    # etiqueta -> U usado
    iteraciones: list = field(default_factory=list)  # [{etiqueta: U}, ...]
    convergido: bool = False
    tol: float = None
    v_pares: list = field(default_factory=list)      # [ParV, ...]
    supercelda_v: tuple = None
    avisos: list = field(default_factory=list)

    @property
    def U(self) -> dict:
        """U por etiqueta de especie (promedio si hay varios sitios)."""
        por = {}
        for s in self.sitios:
            por.setdefault(s.etiqueta, []).append(s.U)
        return {k: float(np.mean(v)) for k, v in por.items()}


# hp.x escribe los V en una tabla de vecinos:
#     Atom 1     Atom 2     Distance (Bohr)   Hubbard V (eV)
#       1 Ni       1 Ni        0.000000          5.0104
#       1 Ni      19 O         3.947603          0.7521
# El segundo índice está en la numeración de la SUPERCELDA que hp.x monta
# para los vecinos, no en la de la celda de entrada: hay que arrastrarlo tal
# cual a la tarjeta HUBBARD o los pares apuntarán a átomos equivocados.
_RE_V = re.compile(
    r"^\s*(\d+)\s+([A-Za-z][a-z]?)\s+(\d+)\s+([A-Za-z][a-z]?)\s+"
    r"(-?\d+\.\d+)\s+(-?\d+\.\d+)\s*$")
_RE_SUPER_V = re.compile(r"supercell\s+(\d+)\s*x\s*(\d+)\s*x\s*(\d+)",
                         re.IGNORECASE)


# ----------------------------------------------------------------------
# Lectura
# ----------------------------------------------------------------------
_RE_SITIO = re.compile(
    r"^\s*(\d+)\s+(\d+)\s+(\S+)\s+(\d+)\s+(\d+)\s+(\S+)\s+"
    r"(-?\d+\.\d+)\s*$")


def leer_parametros(path) -> list:
    """Lee el <prefix>.Hubbard_parameters.dat que escribe hp.x."""
    p = Path(path)
    if not p.exists():
        raise ErrorDeUso(
            f"no existe {p}. hp.x escribe ese archivo al terminar; si no "
            "está, revisa la salida de hp.x.")
    sitios, dentro = [], False
    for linea in p.read_text(errors="ignore").splitlines():
        if "Hubbard U parameters" in linea:
            dentro = True
            continue
        if dentro and "chi0 matrix" in linea:
            break
        m = _RE_SITIO.match(linea)
        if dentro and m:
            sitios.append(SitioU(sitio=int(m.group(1)), tipo=int(m.group(2)),
                                 etiqueta=m.group(3), espin=int(m.group(4)),
                                 U=float(m.group(7))))
    if not sitios:
        raise ErrorDeUso(
            f"{p.name} no trae ninguna U. Si hp.x terminó bien, revisa que "
            "el scf llevara lda_plus_u=.true. y un Hubbard_U de arranque.")
    return sitios


def leer_v(path) -> tuple:
    """Lee los V intersitio del <prefix>.Hubbard_parameters.dat de hp.x.

    Devuelve (lista de ParV, supercelda). La lista viene vacía si hp.x no
    calculó V, que es lo normal salvo que el scf pidiera DFT+U+V.

    Por qué importa: hp.x devuelve los V junto a las U en el mismo archivo,
    y hasta ahora Olla-DFT leía solo las U y tiraba el resto. La diferencia
    entre DFT+U y DFT+U+V no es cosmética: en óxidos de metales de
    transición con enlace covalente apreciable, el término intersitio
    cambia el gap y hasta el orden de las bandas.
    """
    p = Path(path)
    if not p.exists():
        raise ErrorDeUso(
            f"no existe {p}. hp.x escribe ese archivo al terminar.")
    texto = p.read_text(errors="ignore")
    if "Hubbard V parameters" not in texto:
        return [], None

    trozo = texto.split("Hubbard V parameters", 1)[1]
    ms = _RE_SUPER_V.search(trozo[:400])
    super_v = tuple(int(x) for x in ms.groups()) if ms else None

    pares = []
    for linea in trozo.splitlines():
        if "=---" in linea:
            break
        m = _RE_V.match(linea)
        if m:
            i, el_i, j, el_j, dist, v = m.groups()
            pares.append(ParV(i=int(i), el_i=el_i, j=int(j), el_j=el_j,
                              distancia=float(dist), V=float(v)))
    return pares, super_v


def tarjeta_hubbard(sitios: list, pares: list = None,
                    proyeccion: str = "ortho-atomic",
                    umbral_v: float = 0.01) -> str:
    """Tarjeta HUBBARD de QE >= 7.1, con las U y los V.

    ADVERTENCIA de alcance: esta sintaxis es de QE 7.1 en adelante. En QE
    <= 7.0 los V se metían con Hubbard_V(na,nb,k) en &SYSTEM y con
    lda_plus_u_kind=2, que es otra cosa; Olla-DFT no escribe esa forma. El
    generador de esta tarjeta está probado contra la sintaxis documentada,
    no contra una corrida de QE 7.1, porque el QE de esta máquina es 6.6.
    """
    lineas = [f"HUBBARD ({proyeccion})"]
    vistos = set()
    for s in sitios:
        if s.etiqueta in vistos or s.U is None:
            continue
        vistos.add(s.etiqueta)
        lineas.append(f"  U {s.etiqueta}-{ORBITAL_HUBBARD.get(s.etiqueta, '3d')} "
                      f"{s.U:.4f}")
    for pv in (pares or []):
        if pv.es_sitio or pv.V is None or abs(pv.V) < umbral_v:
            continue
        oi = ORBITAL_HUBBARD.get(pv.el_i, "3d")
        oj = ORBITAL_HUBBARD.get(pv.el_j, "2p")
        lineas.append(f"  V {pv.el_i}-{oi} {pv.el_j}-{oj} "
                      f"{pv.i} {pv.j} {pv.V:.4f}")
    return "\n".join(lineas) + "\n"


def report_v(pares: list, supercelda=None, umbral: float = 0.01) -> str:
    if not pares:
        return ("No hay parámetros V en la salida de hp.x.\n"
                "  hp.x solo los calcula si el scf pidió DFT+U+V "
                "(lda_plus_u_kind=2 en QE <= 7.0,\n  o una tarjeta HUBBARD "
                "con términos V en QE >= 7.1). Con un scf de U a secas, "
                "hp.x\n  devuelve solo las U y eso es lo que hay.")
    sitio = [p for p in pares if p.es_sitio]
    inter = [p for p in pares if not p.es_sitio and abs(p.V) >= umbral]
    L = ["--- Parámetros V intersitio ---"]
    if supercelda:
        L.append(f"Vecinos buscados en una supercelda "
                 f"{supercelda[0]}x{supercelda[1]}x{supercelda[2]}")
    if sitio:
        L += ["", "En el sitio (equivalen a la U):"]
        for p in sitio:
            L.append(f"  {p.i:>3d} {p.el_i:<3s}                       "
                     f"{p.V:8.4f} eV")
    L += ["", f"Entre vecinos (por encima de {umbral:g} eV):",
          f"  {'i':>4s} {'':4s} {'j':>4s} {'':4s} {'d (bohr)':>10s} "
          f"{'V (eV)':>9s}"]
    for p in sorted(inter, key=lambda x: -abs(x.V)):
        L.append(f"  {p.i:>4d} {p.el_i:<4s} {p.j:>4d} {p.el_j:<4s} "
                 f"{p.distancia:>10.4f} {p.V:>9.4f}")
    descartados = len([p for p in pares
                       if not p.es_sitio and abs(p.V) < umbral])
    if descartados:
        L.append(f"  ({descartados} pares por debajo de {umbral:g} eV, "
                 f"no se listan)")
    if inter:
        mayor = max(inter, key=lambda x: abs(x.V))
        L += ["",
              f"El V mayor es {mayor.V:.4f} eV entre {mayor.el_i}{mayor.i} y "
              f"{mayor.el_j}{mayor.j}, a {mayor.distancia:.3f} bohr.",
              "  Un V del orden de una décima de eV ya cambia el gap de un "
              "óxido de transición;\n  por eso el nivel U+V no es lo mismo "
              "que U y no se pueden comparar entre sí."]
    return "\n".join(L)


# ----------------------------------------------------------------------
# Preparación
# ----------------------------------------------------------------------
def build_hp_input(prefix: str, qgrid=(2, 2, 2), conv_thr_chi: float = 1e-8,
                   iverbosity: int = 2, solo_contar: bool = False) -> str:
    lineas = [" &inputhp",
              f"   prefix = '{prefix}',",
              "   outdir = './out',",
              f"   nq1 = {qgrid[0]}, nq2 = {qgrid[1]}, nq3 = {qgrid[2]},",
              f"   conv_thr_chi = {conv_thr_chi:.1e},".replace("e-0", "d-0"),
              f"   iverbosity = {iverbosity},"]
    if solo_contar:
        lineas.append("   determine_num_pert_only = .true.,")
    lineas.append(" /")
    return "\n".join(lineas) + "\n"


def elementos_hubbard(atoms) -> list:
    """Qué especies de la estructura son candidatas naturales a DFT+U."""
    return [s for s in dict.fromkeys(atoms.get_chemical_symbols())
            if s in ORBITAL_HUBBARD]


def prepare(atoms, outdir: str = "hubbard", especies=None,
            U_inicial: dict = None, qgrid=(2, 2, 2), pseudo_dir: str = None,
            ecutwfc: float = None, ecutrho: float = None,
            kspacing: float = None, insulator: bool = False,
            proyeccion: str = "ortho-atomic", nspin: int = 1,
            magnetization: dict = None, conv_thr: float = 1e-15,
            hubbard_style: str = "legacy") -> tuple:
    """Escribe el scf con U de arranque y el input de hp.x.

    `hubbard_style` es la sintaxis de DFT+U del scf: "legacy" (lda_plus_u +
    Hubbard_U(i) + U_projection_type en &SYSTEM, QE <= 7.0) o "card" (la
    tarjeta HUBBARD con la proyección en su cabecera, QE >= 7.1, donde la
    sintaxis vieja es un error). Por omisión "legacy", que es lo que se
    escribía siempre; es el mismo selector que `gen --hubbard-style`.
    """
    from qekit.modules import inputgen

    if hubbard_style not in ("legacy", "card"):
        raise ErrorDeUso(
            f"hubbard_style desconocido: {hubbard_style!r}. Usa 'legacy' "
            f"(QE <= 7.0) o 'card' (QE >= 7.1).")
    especies = list(especies or elementos_hubbard(atoms))
    if not especies:
        raise ErrorDeUso(
            "ninguna especie de la estructura es candidata habitual a DFT+U "
            f"({', '.join(dict.fromkeys(atoms.get_chemical_symbols()))}).\n"
            "DFT+U corrige orbitales LOCALIZADOS: d de metales de "
            "transición, f de tierras raras. Si de todas formas quieres "
            "perturbar otra especie, dila con --species.")

    hub = {s: float((U_inicial or {}).get(s, U_SEMILLA)) for s in especies}
    common = sweep.prepare_common(atoms, pseudo_dir, ecutwfc, ecutrho,
                                  insulator, tarea="hubbard")
    out = Path(outdir); out.mkdir(parents=True, exist_ok=True)
    grid = sweep.default_grid(atoms, kspacing)
    scf = inputgen.build_pw_input(
        atoms=atoms, pseudos=common["pseudos"], calculation="scf",
        prefix=common["prefix"], pseudo_dir=common["pseudo_dir"],
        ecutwfc=common["ecutwfc"], ecutrho=common["ecutrho"],
        kcard=f"K_POINTS automatic\n  {grid[0]} {grid[1]} {grid[2]} 0 0 0\n",
        insulator=insulator, degauss=common["degauss"],
        smearing=common["smearing"], hubbard=hub, hubbard_style=hubbard_style,
        nspin=nspin, magnetization=magnetization, conv_thr=conv_thr)
    if hubbard_style == "card":
        # en QE >= 7.1 la proyección va en la cabecera de la tarjeta, y
        # U_projection_type en &SYSTEM es un error
        scf = scf.replace("HUBBARD (ortho-atomic)", f"HUBBARD ({proyeccion})",
                          1)
    else:
        scf = _fijar_proyeccion(scf, proyeccion)
    sweep.write_input(out / "scf.in", scf)
    sweep.write_input(out / "hp.in",
                      build_hp_input(common["prefix"], qgrid))

    ncel = int(np.prod(qgrid))
    rep = ["--- U de Hubbard por respuesta lineal ---",
           f"Estructura: {atoms.get_chemical_formula()} ({len(atoms)} átomos)",
           f"Especies perturbadas: {', '.join(especies)}  "
           f"(orbital {', '.join(ORBITAL_HUBBARD.get(s, '?') for s in especies)})",
           "U de arranque: " + ", ".join(f"{k}={v:g} eV" for k, v in hub.items()),
           f"Proyección: {proyeccion}",
           f"Malla de q: {qgrid[0]}x{qgrid[1]}x{qgrid[2]}  "
           f"(equivale a una supercelda de {ncel} celdas)",
           "",
           f"Archivos en '{out.resolve()}': scf.in, hp.in",
           "Orden:  pw.x -in scf.in   ->   hp.x -in hp.in",
           ""]
    if ncel == 1:
        rep += ["AVISO: con nq = 1x1x1 la perturbación ve sus propias "
                "imágenes periódicas y\nel U sale mal. Usa al menos 2x2x2 y "
                "comprueba que el número no cambie al\nsubir la malla.", ""]
    rep += ["El U que salga de aquí es de PRIMERA ITERACIÓN: depende del U "
            "que se usó\nen el scf. El autoconsistente sale de repetir el "
            "ciclo (olla-dft hubbard --cycle).",
            "",
            f"El número solo vale con la MISMA proyección ('{proyeccion}'). "
            "Un U de la\nliteratura calculado con otra proyección no es "
            "comparable."]
    warn = sweep.missing_pseudo_warning(common)
    if warn:
        rep.append(warn)
    return common, "\n".join(rep)


def _fijar_proyeccion(scf: str, proyeccion: str) -> str:
    """Mete U_projection_type en &SYSTEM.

    Va aquí y no en el constructor general porque solo tiene sentido con
    DFT+U, y porque el valor por omisión de QE ('atomic') NO es el que
    usa hp.x en sus ejemplos: mezclarlos da un U que no corresponde.
    """
    if "U_projection_type" in scf:
        return scf
    return scf.replace("  lda_plus_u ",
                       f"  U_projection_type = '{proyeccion}'\n  lda_plus_u ",
                       1) if "lda_plus_u " in scf else re.sub(
        r"(&SYSTEM\n)", rf"\1  U_projection_type = '{proyeccion}'\n", scf,
        count=1)


# ----------------------------------------------------------------------
# Ejecución
# ----------------------------------------------------------------------
def run_hp(workdir, cmd: str = None, nproc: int = None,
           stem: str = "hp") -> str:
    work = Path(workdir)
    base = run_mod.build_command(cmd, nproc)
    exe = Path(base[-1]).parent / "hp.x" if "/" in base[-1] else Path("hp.x")
    if not shutil.which(str(exe)) and not Path(exe).exists():
        raise ErrorDeUso(
            f"no se encontró hp.x junto a pw.x ('{exe}'). Es parte de "
            "Quantum ESPRESSO pero no se compila por defecto:\n"
            "  cd <fuente de QE> && make hp")
    with open(work / f"{stem}.in") as fin, open(work / f"{stem}.out", "w") as fo:
        proc = subprocess.run(base[:-1] + [str(exe)], stdin=fin, stdout=fo,
                              stderr=subprocess.STDOUT, cwd=str(work))
    texto = (work / f"{stem}.out").read_text(errors="ignore")
    if proc.returncode != 0 or "JOB DONE" not in texto:
        raise RuntimeError(run_mod.failure_message(
            "hp.x", work / f"{stem}.out", texto))
    return str(work / f"{stem}.out")


def collect(path, qgrid=None, proyeccion: str = "ortho-atomic") -> HubbardRun:
    p = Path(path)
    dats = sorted(p.glob("*.Hubbard_parameters.dat"))
    if not dats:
        raise ErrorDeUso(
            f"no hay ningún *.Hubbard_parameters.dat en {p}. Corre primero:\n"
            "  pw.x -in scf.in  &&  hp.x -in hp.in")
    run = HubbardRun(sitios=leer_parametros(dats[0]), qgrid=qgrid,
                     proyeccion=proyeccion)
    run.iteraciones = [run.U]
    return run


def report(run: HubbardRun) -> str:
    lines = ["--- U de Hubbard por respuesta lineal ---"]
    if run.qgrid:
        lines.append(f"Malla de q: {run.qgrid[0]}x{run.qgrid[1]}x"
                     f"{run.qgrid[2]}")
    lines.append(f"Proyección: {run.proyeccion}")
    lines += ["", f"{'sitio':>6s} {'especie':>9s} {'U (eV)':>9s}"]
    for s in run.sitios:
        lines.append(f"{s.sitio:6d} {s.etiqueta:>9s} {s.U:9.4f}")

    if len(run.iteraciones) > 1:
        lines += ["", "Ciclo de autoconsistencia:"]
        etiquetas = sorted(run.iteraciones[-1])
        cab = "  iter  " + "  ".join(f"{e:>9s}" for e in etiquetas)
        lines.append(cab)
        for i, it in enumerate(run.iteraciones):
            fila = f"  {i:4d}  " + "  ".join(
                f"{it.get(e, float('nan')):9.4f}" for e in etiquetas)
            lines.append(fila)
        prim, ult = run.iteraciones[0], run.iteraciones[-1]
        mov = max(abs(ult[e] - prim.get(e, ult[e])) for e in etiquetas)
        lines += ["",
                  f"El U se movió {mov:.3f} eV entre la primera vuelta y la "
                  f"última."]
        if run.convergido:
            lines.append(f"Convergido: el cambio bajó de {run.tol:g} eV.")
        else:
            lines.append(
                "NO convergido: se agotaron las iteraciones. El número de "
                "abajo es\nprovisional; sube --max-iter o afloja --tol.")
        if mov > 0.5:
            lines.append(
                f"  Que se moviera {mov:.2f} eV es la razón por la que una "
                "sola vuelta no basta:\n  un U de primera iteración habría "
                "dado un número bastante distinto.")
    else:
        lines += ["",
                  "U de PRIMERA ITERACIÓN. Depende del U que llevaba el scf "
                  "de partida.\nEl autoconsistente sale de repetir el ciclo: "
                  "olla-dft hubbard ... --cycle"]

    for a in run.avisos:
        lines += ["", a]
    lines += ["",
              f"Este U solo vale con la proyección '{run.proyeccion}'. Un U "
              "de la literatura\ncalculado con otra proyección no es el mismo "
              "número, aunque el elemento\ny el compuesto coincidan.",
              "",
              "Cómo usarlo:  olla-dft gen estructura.cif --hubbard " +
              " ".join(f"{k}={v:.2f}" for k, v in run.U.items())]
    return "\n".join(lines)


def export(run: HubbardRun, outdir: str = ".") -> list:
    out = Path(outdir); out.mkdir(parents=True, exist_ok=True)
    f = out / "HUBBARD_U.dat"
    cab = provenance.header_plain(
        "U de Hubbard por respuesta lineal",
        {"proyeccion": run.proyeccion,
         "qgrid": "x".join(str(q) for q in (run.qgrid or ())) or None,
         "iteraciones": len(run.iteraciones),
         "convergido": run.convergido},
        titulo="Parametros de Hubbard (hp.x)")
    filas = ["# sitio  especie   U(eV)"]
    for s in run.sitios:
        filas.append(f"{s.sitio:7d} {s.etiqueta:>8s} {s.U:9.4f}")
    f.write_text(cab + "\n" + "\n".join(filas) + "\n")
    txt = out / "HUBBARD_U.txt"
    txt.write_text(report(run) + "\n")
    return [str(f), str(txt)]


# ----------------------------------------------------------------------
# Ciclo de autoconsistencia
# ----------------------------------------------------------------------
def ciclo(atoms, outdir: str = "hubbard", especies=None, qgrid=(2, 2, 2),
          max_iter: int = 6, tol: float = 0.05, mezcla: float = 1.0,
          pw_cmd: str = None, hp_cmd: str = None, nproc: int = None,
          **kwargs) -> HubbardRun:
    """Repite scf -> hp.x -> scf con el U nuevo, hasta que deje de moverse.

    Es el punto entero del módulo. Con una sola vuelta se obtiene el U que
    corresponde al scf de partida —normalmente uno sin U—, y ese número
    puede estar a más de 1 eV del autoconsistente.

    `mezcla` amortigua el paso: U_nuevo = (1-m)*U_viejo + m*U_calculado.
    Con m=1 se toma el valor tal cual, que converge rápido cuando
    converge; bajarlo a 0.5 ayuda cuando el ciclo oscila.

    Cada iteración vive en su propia subcarpeta para poder auditarla, y
    ninguna se pisa con la anterior.
    """
    base = Path(outdir)
    base.mkdir(parents=True, exist_ok=True)
    especies = list(especies or elementos_hubbard(atoms))
    U = {s: U_SEMILLA for s in especies}
    run = HubbardRun(qgrid=tuple(qgrid), tol=tol,
                     proyeccion=kwargs.get("proyeccion", "ortho-atomic"))
    run.U_entrada = dict(U)

    for it in range(max_iter):
        paso = base / f"iter{it:02d}"
        prepare(atoms, outdir=str(paso), especies=especies, U_inicial=U,
                qgrid=qgrid, **kwargs)
        run_mod.run_all(
            [run_mod.Job(name="scf", directory=paso, input_file="scf.in",
                          output_file="scf.out")],
            pw_cmd=pw_cmd, nproc=nproc, verbose=False)
        run_hp(paso, cmd=hp_cmd, nproc=nproc)
        parcial = collect(paso, qgrid=tuple(qgrid),
                          proyeccion=run.proyeccion)
        nuevo = parcial.U
        run.sitios = parcial.sitios
        run.iteraciones.append(dict(nuevo))

        cambio = max(abs(nuevo[e] - U[e]) for e in nuevo)
        # en la primera vuelta el U de partida es 1e-8, así que el "cambio"
        # es el propio valor: no cuenta como convergencia
        if it > 0 and cambio < tol:
            run.convergido = True
            U = nuevo
            break
        U = {e: (1 - mezcla) * U[e] + mezcla * nuevo[e] for e in nuevo}

    if not run.convergido:
        run.avisos.append(
            f"Se hicieron {len(run.iteraciones)} vueltas sin bajar de "
            f"{tol} eV. Mira la tabla: si el número oscila arriba y abajo, "
            "baja --mixing a 0.5; si baja despacio pero siempre en el mismo "
            "sentido, sube --max-iter.")
    return run
