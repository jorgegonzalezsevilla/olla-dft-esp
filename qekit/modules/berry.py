# Olla-DFT — command-line toolkit for Quantum ESPRESSO
# Copyright (C) 2026 Jorge Enrique González Sevilla
# SPDX-License-Identifier: AGPL-3.0-or-later
# This program is free software: you can redistribute it and/or modify it
# under the terms of the GNU Affero General Public License as published by the Free
# Software Foundation, either version 3 of the License, or (at your option)
# any later version. See the LICENSE file for details.

"""Polarización eléctrica por fase de Berry, y lo que se deduce de ella.

La polarización de un sólido periódico NO es la integral del dipolo en la
celda. Esa integral depende de dónde pongas los bordes de la celda, así que
no es una propiedad del material: el mismo cristal da valores distintos
según cómo lo cortes. Durante décadas eso pareció un problema del concepto;
resultó ser un problema de la fórmula.

King-Smith y Vanderbilt (1993) la arreglaron: lo que es observable es la
FASE DE BERRY de los estados ocupados a lo largo de una dirección,

    P_el = −(f·e/(2π·Ω)) · Im ln Π_j det S(k_j, k_{j+1}) · R

donde S son los solapes entre puntos k consecutivos de una cuerda. Esa fase
está definida módulo 2π, así que P está definida módulo el "cuanto de
polarización" e·R/Ω. No es un defecto de la teoría: es la afirmación de que
solo las DIFERENCIAS de polarización son medibles, que es exactamente lo
que dice el experimento —nadie mide P, se mide la carga que circula al
cambiar la estructura—.

De ahí salen las dos reglas que hacen falta para no decir tonterías:

  1. **Un valor de P, solo, no significa nada.** Hay que dar ΔP entre dos
     estructuras, y hay que llegar de una a otra por un camino con pasos
     pequeños, siguiendo la rama. Si te saltas la rama, el resultado sale
     desplazado un cuanto entero y sigue pareciendo razonable.
  2. **Solo vale para aislantes.** En un metal la fase de Berry no está
     definida: no hay una banda ocupada separada que seguir.

Lo que este módulo aporta sobre correr `lberry` a mano:

  - Escribe las cuerdas de puntos k bien. El espaciado dentro de la cuerda
    es G/(nppstr−1), no G/nppstr, y el último punto es el primero más G.
  - Fuerza `nosym` y `noinv`. Sin eso Quantum ESPRESSO EXPANDE la lista de
    puntos k con las operaciones de simetría —324 se convirtieron en 1438 en
    la prueba— y las cuerdas dejan de ser cuerdas. El error que da
    ("Wrong k-strings weights?") no menciona la simetría por ninguna parte.
  - Sigue la rama a lo largo de un camino adiabático y avisa si un paso se
    acerca al cuanto, que es cuando el seguimiento deja de ser fiable.
  - Comprueba la parte iónica contra su fórmula exacta, que es aritmética:
    φ_ion = Σ_a Z_a·f_a, con f_a la coordenada fraccionaria. Si esa no
    cuadra, la geometría o las valencias están mal y da igual lo demás.
  - Y la parte electrónica contra los centros de Wannier, que salen de un
    camino completamente distinto (`olla-dft wannier`). Que dos rutas
    independientes coincidan es la única validación que vale.
"""

from dataclasses import dataclass, field
from pathlib import Path
import re

import numpy as np

from qekit.core.errors import ErrorDeUso, FaltanDatos
from qekit.core import plataforma, style as qstyle

BOHR = 0.529177210903
# e/Å² a C/m²:  1.602176634e-19 C / (1e-10 m)² = 16.0217663
E_A2_A_C_M2 = 16.02176634
# fracción del cuanto por encima de la cual un paso del camino es sospechoso
FRACCION_SOSPECHOSA = 0.25


# ----------------------------------------------------------------------
# Aritmética de la fase
# ----------------------------------------------------------------------
def _nint(x):
    """El NINT de Fortran: el medio se va hacia afuera, no hacia el par.

    numpy redondea 0.5 al par y Fortran lo aleja del cero. La diferencia solo
    aparece justo en el borde —una fase de exactamente medio cuanto, que es
    precisamente lo que da un cristal centrosimétrico— y ahí hace que Olla-DFT
    diga +1 donde Quantum ESPRESSO dice −1. Son la misma fase, pero la
    comprobación contra QE fallaría.
    """
    x = np.asarray(x, float)
    return np.copysign(np.floor(np.abs(x) + 0.5), x)


def a_intervalo(fase, modulo=2.0):
    """Trae una fase al intervalo [−modulo/2, modulo/2), como hace QE.

    El borde es cerrado por abajo: una fase de exactamente medio cuanto sale
    NEGATIVA, que es lo que escribe pw.x para un cristal centrosimétrico
    (−1.00000 en el silicio).
    """
    f = np.asarray(fase, float)
    r = f - modulo * _nint(f / modulo)
    return float(r) if r.ndim == 0 else r


def modulo_de(valencias, atoms=None):
    """¿El cuanto vale 2 o 1? Depende de si alguna valencia es impar.

    Un ion con carga de valencia PAR, al moverse un vector de red, cambia la
    fase en un múltiplo de 2 y no se nota. Si la carga es impar, cambia en un
    número impar, y entonces la fase solo está definida módulo 1: **el cuanto
    de polarización es la mitad**. Quantum ESPRESSO lo pone en la salida como
    MOD_TOT, y pasa desapercibido con facilidad.

    Es un factor 2 en el cuanto, así que es también un factor 2 en cuántos
    puntos hace falta en un camino adiabático antes de que el seguimiento de
    la rama deje de ser fiable. En BN cúbico (Z = 3 y 5) el cuanto es la
    mitad del que tendría el silicio con la misma celda.
    """
    if atoms is not None:
        zs = [valencias[s] for s in atoms.get_chemical_symbols()
              if s in valencias]
    else:
        zs = list(valencias.values())
    if not zs:
        return 2.0
    return 1.0 if any(int(round(z)) % 2 == 1 for z in zs) else 2.0


def fase_ionica(atoms, valencias, gdir=3, plegar=True):
    """φ_ion = Σ_a Z_a · f_a, en las unidades de QE (el cuanto vale 2 o 1).

    Es aritmética pura: coordenadas fraccionarias por carga de valencia. No
    hace falta ningún cálculo para saberla, y por eso es la primera
    comprobación: si QE no da esto, el input no describe la estructura que
    crees.

    Con `plegar` se reproduce exactamente lo que hace QE, que pliega la
    contribución de CADA ion por separado —módulo 1 si su valencia es impar,
    módulo 2 si es par— y solo después suma. Sin ese plegado por ion el
    número sale desplazado un entero y la comprobación falla sin motivo.
    """
    frac = atoms.get_scaled_positions()
    simb = atoms.get_chemical_symbols()
    faltan = sorted({s for s in simb if s not in valencias})
    if faltan:
        raise FaltanDatos(
            f"no sé cuántos electrones de valencia tiene el pseudo de "
            f"{', '.join(faltan)}. Sale del propio UPF (z_valence) o de la "
            f"salida de pw.x.")
    j = int(gdir) - 1
    total = 0.0
    for i, s in enumerate(simb):
        z = valencias[s]
        x = z * frac[i, j]
        if plegar:
            m = 1.0 if int(round(z)) % 2 == 1 else 2.0
            x -= m * _nint(x / m)
        total += x
    if plegar:
        m = modulo_de(valencias, atoms)
        total -= m * _nint(total / m)
    return float(total)


def cuanto(atoms, gdir=3, modulo=2.0):
    """El cuanto de polarización: modulo · e·R/Ω, en e/Å² y en C/m².

    `modulo` es el MOD_TOT de Quantum ESPRESSO: 2 si todas las cargas de
    valencia son pares, 1 si alguna es impar (ver `modulo_de`). Darlo por
    hecho igual a 2 duplica el cuanto y hace creer que un salto de rama cabe
    dentro del margen cuando no cabe.
    """
    cell = np.array(atoms.cell.array, float)
    R = cell[int(gdir) - 1]
    vol = float(abs(np.linalg.det(cell)))
    q = float(modulo) * float(np.linalg.norm(R)) / vol      # e/Å²
    return q, q * E_A2_A_C_M2


def polarizacion(fase_total, atoms, gdir=3):
    """P a lo largo de R_gdir, a partir de la fase en unidades de QE.

    P·Ω/e = φ · R_gdir  (como VECTOR: es la proyección sobre R_gdir, no el
    módulo de P). De una sola dirección sale una sola componente; para el
    vector completo hacen falta las tres.
    """
    cell = np.array(atoms.cell.array, float)
    R = cell[int(gdir) - 1]
    vol = float(abs(np.linalg.det(cell)))
    P = float(fase_total) * float(np.linalg.norm(R)) / vol
    return P, P * E_A2_A_C_M2


def desenrollar(fases, modulo=2.0):
    """Sigue la rama a lo largo de un camino: el paso más pequeño gana.

    Cada fase llega plegada al intervalo (−1, 1]. Al recorrer un camino
    adiabático la fase física varía poco entre pasos consecutivos, así que
    de todas las imágenes φ + n·2 se elige la más cercana a la anterior.
    Esto es lo que convierte una lista de números sin sentido individual en
    una curva ΔP(λ) que sí lo tiene.

    Devuelve (fases_seguidas, saltos), donde `saltos` es el tamaño de cada
    paso ya seguido: si alguno se acerca al cuanto, el muestreo en λ es
    demasiado grueso y el seguimiento es una apuesta, no un cálculo.
    """
    f = np.asarray(fases, float).ravel()
    if len(f) == 0:
        return f, np.array([])
    fuera = [float(f[0])]
    for x in f[1:]:
        n = np.round((fuera[-1] - x) / modulo)
        fuera.append(float(x + n * modulo))
    fuera = np.array(fuera)
    return fuera, np.abs(np.diff(fuera))


def desde_wannier(centros, atoms, valencias, gdir=3, spin=2.0):
    """La misma fase, pero desde los centros de Wannier.

    φ_el = −f · Σ_n (r̄_n · b_gdir)/2π, con f = 2 por el espín. Es la MISMA
    fase de Berry: los centros de Wannier son literalmente
    −(1/2π)·Im ln⟨u|u⟩ integrado. Que este camino y el de `lberry` den lo
    mismo, calculados por rutinas que no comparten una línea de código, es
    la comprobación que vale.
    """
    cell = np.array(atoms.cell.array, float)
    frac = np.atleast_2d(np.asarray(centros, float)) @ np.linalg.inv(cell)
    el = -float(spin) * float(frac[:, int(gdir) - 1].sum())
    ion = fase_ionica(atoms, valencias, gdir)
    return a_intervalo(el), a_intervalo(ion), a_intervalo(el + ion)


# ----------------------------------------------------------------------
# Leer lo que escribe pw.x con lberry
# ----------------------------------------------------------------------
_RE_FASES = re.compile(
    r"Ionic Phase:\s*(-?[\d.]+)\s*.*?Electronic Phase:\s*(-?[\d.]+)"
    r"\s*.*?TOTAL PHASE:\s*(-?[\d.]+)\s*MOD_TOT:\s*(\d+)", re.S)
_RE_P = re.compile(r"P =\s*(-?[\d.]+)\s*\(mod\s*([\d.]+)\)\s*"
                   r"\(e/Omega\)\.bohr")
_RE_DIR = re.compile(r"direction of vector\s*(\d+)")
_RE_G = re.compile(r"G-vector along string \(2 pi/a\):\s*"
                   r"(-?[\d.]+)\s+(-?[\d.]+)\s+(-?[\d.]+)")
_RE_NPP = re.compile(r"Number of k-points per string:\s*(\d+)")
_RE_NSTR = re.compile(r"Number of different strings\s*:\s*(\d+)")


@dataclass
class Berry:
    """Una corrida de lberry: una dirección, una estructura."""
    gdir: int = 3
    fase_ion: float = float("nan")
    fase_el: float = float("nan")
    fase_total: float = float("nan")
    modulo: float = 2.0
    P_bohr: float = float("nan")       # (e/Ω)·bohr, como lo escribe QE
    cuanto_bohr: float = float("nan")
    nppstr: int = 0
    ncuerdas: int = 0
    etiqueta: str = ""
    fuente: str = ""


def leer_berry(salida):
    """Extrae las fases de una salida de pw.x con lberry = .true."""
    p = Path(salida)
    if p.is_dir():
        cand = sorted(p.glob("*.out"))
        cand = [c for c in cand if "POLARIZATION CALCULATION" in
                c.read_text(errors="replace")]
        if not cand:
            raise FaltanDatos(
                f"en {p} no hay ninguna salida de pw.x con la sección "
                f"POLARIZATION CALCULATION. ¿Corriste el paso con "
                f"lberry = .true.?")
        p = cand[0]
    if not p.exists():
        raise FaltanDatos(
            f"falta {p}. Es la salida del paso con lberry: córrelo con "
            f"`--run`, o a mano con\n  `bash correr.sh` dentro de la carpeta "
            f"del cálculo.")
    txt = p.read_text(errors="replace")
    m = _RE_FASES.search(txt)
    if not m:
        if "Wrong k-strings" in txt:
            raise FaltanDatos(
                "pw.x se paró con 'Wrong k-strings': la lista de puntos k no "
                "forma cuerdas.\nCasi siempre es que faltó nosym = .true. y "
                "noinv = .true., y QE expandió la\nlista con las operaciones "
                "de simetría. Vuelve a preparar el cálculo con Olla-DFT.")
        raise FaltanDatos(f"{p} no trae el resumen de fases de lberry.")
    b = Berry(fase_ion=float(m.group(1)), fase_el=float(m.group(2)),
              fase_total=float(m.group(3)), modulo=float(m.group(4)),
              fuente=str(p))
    mp = _RE_P.search(txt)
    if mp:
        b.P_bohr, b.cuanto_bohr = float(mp.group(1)), float(mp.group(2))
    md = _RE_DIR.search(txt)
    if md:
        b.gdir = int(md.group(1))
    mn = _RE_NPP.search(txt)
    if mn:
        b.nppstr = int(mn.group(1))
    ms = _RE_NSTR.search(txt)
    if ms:
        b.ncuerdas = int(ms.group(1))
    return b


def valencias_de(salida):
    """Carga de valencia de cada especie, leída de la salida de pw.x.

    Se saca de la tabla "atomic species / valence" que escribe pw.x, no del
    número atómico: lo que cuenta en la fase iónica es la carga que ve el
    pseudopotencial, y un pseudo con semicapa tiene otra.
    """
    p = Path(salida)
    if p.is_dir():
        cand = [c for c in sorted(p.glob("*.out"))
                if "atomic species" in c.read_text(errors="replace")]
        if not cand:
            raise FaltanDatos(f"no hay salidas de pw.x en {p}.")
        p = cand[0]
    txt = p.read_text(errors="replace")
    m = re.search(r"atomic species\s+valence\s+mass\s+pseudopotential\s*\n"
                  r"((?:\s+\S+\s+[\d.]+\s+[\d.]+\s+.*\n)+)", txt)
    fuera = {}
    if m:
        for linea in m.group(1).strip().split("\n"):
            t = linea.split()
            if len(t) >= 3:
                fuera[t[0]] = float(t[1])
    return fuera


# ----------------------------------------------------------------------
# Preparar: cuerdas de puntos k bien hechas
# ----------------------------------------------------------------------
def cuerdas(nppstr, gdir=3, kperp=(6, 6)):
    """Los puntos k que lberry necesita, en el orden que necesita.

    Una "cuerda" son nppstr puntos igualmente espaciados que recorren TODA
    la zona de Brillouin a lo largo de b_gdir: el último es el primero más
    un vector completo de la red recíproca. Por eso el espaciado es
    1/(nppstr−1) y no 1/nppstr, que es el error que hace que pw.x se pare
    con 'Wrong k-strings?' sin más explicación.

    El índice de la cuerda corre despacio y el de dentro de la cuerda
    deprisa; pw.x lo da por hecho y no lo comprueba más que por diferencias.
    """
    npp = int(nppstr)
    if npp < 3:
        raise ErrorDeUso(
            f"nppstr = {npp} es demasiado corto: con menos de tres puntos la "
            f"cuerda no muestrea la zona de Brillouin. Usa 7 o más, y "
            f"comprueba que la fase no cambia al subirlo.")
    perp = [int(x) for x in kperp]
    if len(perp) != 2 or min(perp) < 1:
        raise ErrorDeUso("--kperp son dos enteros positivos, por ejemplo 6x6.")
    otros = [d for d in (0, 1, 2) if d != int(gdir) - 1]
    ks = []
    for i in range(perp[0]):
        for j in range(perp[1]):
            for l in range(npp):
                k = [0.0, 0.0, 0.0]
                k[otros[0]] = i / perp[0]
                k[otros[1]] = j / perp[1]
                k[int(gdir) - 1] = l / (npp - 1)
                ks.append(k)
    return np.array(ks)


def _interpolar_estructuras(a, b, lam):
    """Estructura a fracción lam del camino de a a b, por imagen mínima.

    Las posiciones se interpolan en fraccionarias y con la imagen más
    cercana: si un átomo cruza el borde de la celda, interpolar en
    cartesianas lo mandaría de paseo por el centro del cristal.
    """
    if a.get_chemical_symbols() != b.get_chemical_symbols():
        raise ErrorDeUso(
            "las dos estructuras tienen especies distintas o en otro orden. "
            "El camino\nadiabático une la MISMA estructura en dos "
            "configuraciones; si no, no hay camino.")
    fa, fb = a.get_scaled_positions(), b.get_scaled_positions()
    d = fb - fa
    d -= np.round(d)                       # imagen mínima
    nueva = a.copy()
    nueva.set_cell(a.cell.array + lam * (b.cell.array - a.cell.array),
                   scale_atoms=False)
    nueva.set_scaled_positions(fa + lam * d)
    return nueva


@dataclass
class BerryRun:
    """Un camino adiabático de polarización, o un punto suelto de él."""
    formula: str = ""
    gdir: int = 3
    nppstr: int = 9
    kperp: tuple = (6, 6)
    lambdas: list = field(default_factory=list)
    estructuras: list = field(default_factory=list)
    jobs: list = field(default_factory=list)
    puntos: list = field(default_factory=list)      # objetos Berry
    valencias: dict = field(default_factory=dict)
    etiqueta_ref: str = "referencia"
    etiqueta_fin: str = "estructura polar"
    es_desplazamiento: bool = False
    desplazamiento: np.ndarray = None               # Å, si es un Z*
    avisos: list = field(default_factory=list)

    @property
    def cell(self):
        return np.array(self.estructuras[-1].cell.array, float) \
            if self.estructuras else None


def prepare(atoms, outdir="berry", gdir=3, nppstr=9, kperp=(6, 6),
            referencia=None, nlambda=5, desplazar=None,
            pseudo_dir=None, ecutwfc=None, ecutrho=None, kgrid_scf=None):
    """Escribe el scf y una corrida de lberry por cada punto del camino.

    Con `referencia` se construye el camino adiabático de la estructura de
    referencia (normalmente la centrosimétrica, donde P = 0 por simetría) a
    la polar, que es la ÚNICA forma correcta de dar una polarización
    espontánea. Con `desplazar` se construye un camino de desplazamiento de
    un átomo, y la pendiente de P frente al desplazamiento es la carga
    efectiva de Born.

    Sin ninguna de las dos se calcula un solo punto, que sirve para
    comprobaciones pero no para publicar: un valor de P aislado está
    definido módulo el cuanto y no significa nada por sí mismo.
    """
    from qekit.modules import inputgen, sweep

    gdir = int(gdir)
    if gdir not in (1, 2, 3):
        raise ErrorDeUso(f"--gdir es 1, 2 o 3 (el vector de la red "
                         f"recíproca); recibí {gdir}.")
    run = BerryRun(formula=atoms.get_chemical_formula(), gdir=gdir,
                   nppstr=int(nppstr), kperp=tuple(int(x) for x in kperp))

    if referencia is not None and desplazar is not None:
        raise ErrorDeUso(
            "elige un camino: o interpolas hacia una estructura de "
            "referencia, o desplazas\nun átomo. Los dos a la vez no definen "
            "un camino.")

    if referencia is not None:
        lams = np.linspace(0.0, 1.0, int(nlambda))
        run.estructuras = [_interpolar_estructuras(referencia, atoms, l)
                           for l in lams]
        run.lambdas = list(map(float, lams))
    elif desplazar is not None:
        idx, vec = desplazar
        if not 0 <= int(idx) < len(atoms):
            raise ErrorDeUso(
                f"el átomo {idx + 1} no existe: la estructura tiene "
                f"{len(atoms)}.")
        lams = np.linspace(0.0, 1.0, int(nlambda))
        run.estructuras = []
        for l in lams:
            a2 = atoms.copy()
            a2.positions[int(idx)] += l * np.asarray(vec, float)
            run.estructuras.append(a2)
        run.lambdas = list(map(float, lams))
        run.es_desplazamiento = True
        run.desplazamiento = np.asarray(vec, float)
        run.etiqueta_ref = "sin desplazar"
        run.etiqueta_fin = (f"átomo {int(idx) + 1} desplazado "
                            f"{np.linalg.norm(vec):.3f} Å")
    else:
        run.estructuras = [atoms]
        run.lambdas = [1.0]

    common = sweep.prepare_common(atoms, pseudo_dir, ecutwfc, ecutrho,
                                  insulator=True)
    out = Path(outdir); out.mkdir(parents=True, exist_ok=True)
    k_scf = tuple(kgrid_scf or sweep.default_grid(atoms, None))
    ks = cuerdas(nppstr, gdir, kperp)
    card = ("K_POINTS crystal\n" + f"{len(ks)}\n" +
            "".join(f" {a:.10f} {b:.10f} {c:.10f} {1.0 / len(ks):.10f}\n"
                    for a, b, c in ks))

    for i, (lam, at) in enumerate(zip(run.lambdas, run.estructuras)):
        d = out / f"p{i:02d}"
        d.mkdir(parents=True, exist_ok=True)
        for nombre, calc, kc, extra in (
                ("1_scf", "scf",
                 f"K_POINTS automatic\n  {k_scf[0]} {k_scf[1]} {k_scf[2]} "
                 "0 0 0\n", ""),
                ("2_berry", "nscf", card, "berry")):
            txt = inputgen.build_pw_input(
                atoms=at, pseudos=common["pseudos"], calculation=calc,
                prefix=common["prefix"], pseudo_dir=common["pseudo_dir"],
                ecutwfc=common["ecutwfc"], ecutrho=common["ecutrho"],
                kcard=kc, insulator=True, conv_thr=1e-10,
                nosym=(extra == "berry"))
            if extra == "berry":
                # lberry, gdir y nppstr van en &CONTROL, no en &SYSTEM: es
                # el primer sitio donde se equivoca todo el mundo
                txt = re.sub(
                    r"(&CONTROL\n)",
                    r"\1  lberry           = .true.\n"
                    f"  gdir             = {gdir}\n"
                    f"  nppstr           = {int(nppstr)}\n",
                    txt, count=1)
            sweep.write_input(d / f"{nombre}.in", txt)
        plataforma.escribir_par_de_guiones(d, [
            ("pw.x", "1_scf.in", "1_scf.out"),
            ("pw.x", "2_berry.in", "2_berry.out"),
        ])
        run.jobs.append(d)

    plataforma.escribir_script(
        out / "correr.sh",
        "#!/bin/bash\nset -e\n"
        "# En Windows, o sin bash:  python correr.py\n"
        "for d in p*/; do (cd \"$d\" && bash correr.sh) || exit 1; done\n")
    plataforma.escribir_script(
        out / "correr.py",
        "#!/usr/bin/env python3\n"
        "# Corre todos los puntos del camino, en orden. Generado por Olla-DFT.\n"
        "import subprocess, sys\n"
        "from pathlib import Path\n"
        "aqui = Path(__file__).resolve().parent\n"
        "for d in sorted(aqui.glob('p[0-9][0-9]')):\n"
        "    print('>>', d.name, flush=True)\n"
        "    r = subprocess.run([sys.executable, str(d / 'correr.py')])\n"
        "    if r.returncode:\n"
        "        sys.exit(r.returncode)\n"
        "print('Listo.')\n")

    q_eA2, q_cm2 = cuanto(atoms, gdir)
    rep = [f"--- Polarización por fase de Berry: {run.formula} ---",
           f"Dirección: vector {gdir} de la red recíproca",
           f"Cuerdas: {nppstr} puntos por cuerda × "
           f"{run.kperp[0]}×{run.kperp[1]} perpendiculares = {len(ks)} "
           f"puntos k",
           f"Cuanto de polarización: {q_eA2:.6f} e/Å² = {q_cm2:.4f} C/m²",
           ""]
    if len(run.lambdas) > 1:
        rep += [f"Camino adiabático de {len(run.lambdas)} puntos, de "
                f"«{run.etiqueta_ref}» a «{run.etiqueta_fin}».",
                "Solo la DIFERENCIA a lo largo del camino es física; cada "
                "valor suelto está",
                "definido módulo el cuanto.", ""]
    else:
        run.avisos.append(
            "Un solo punto. P está definida módulo el cuanto, así que este "
            "número por sí solo\n  no significa nada: sirve para comprobar la "
            "parte iónica y poco más. Para una\n  polarización espontánea "
            "hace falta --reference con la estructura centrosimétrica;\n  "
            "para una carga de Born, --displace.")
    rep += [f"Archivos en '{out.resolve()}':",
            f"  p00..p{len(run.jobs) - 1:02d}/   cada punto del camino "
            f"(1_scf.in y 2_berry.in)",
            "  bash correr.sh   los lanza todos en orden",
            "",
            "El nscf lleva nosym y noinv a la fuerza. Sin eso Quantum "
            "ESPRESSO expande la",
            "lista de puntos k con las operaciones de simetría, las cuerdas "
            "se rompen y pw.x",
            "se para con 'Wrong k-strings weights?', que no menciona la "
            "simetría por ningún lado."]
    warn = sweep.missing_pseudo_warning(common)
    if warn:
        rep.append(warn)
    return run, common, "\n".join(rep)


def correr(run, pw_cmd=None, nproc=None, timeout=None, verbose=True,
           rehacer=False):
    """Lanza scf + lberry en cada punto del camino, en orden."""
    import subprocess
    from qekit.core import runner as run_mod

    cmd = run_mod.build_command(pw_cmd, nproc)
    for d in run.jobs:
        d = Path(d)
        for nombre in ("1_scf", "2_berry"):
            ent, sal = d / f"{nombre}.in", d / f"{nombre}.out"
            if (sal.exists() and not rehacer
                    and "JOB DONE" in sal.read_text(errors="replace")):
                if verbose:
                    print(f"  {d.name}/{nombre}  (ya estaba)")
                continue
            if verbose:
                print(f"  {d.name}/{nombre} ...", end="", flush=True)
            with open(ent) as fi, open(sal, "w") as fo:
                r = subprocess.run(cmd, stdin=fi, stdout=fo,
                                   stderr=subprocess.STDOUT, cwd=str(d),
                                   timeout=timeout)
            txt = sal.read_text(errors="replace")
            ok = r.returncode == 0 and "JOB DONE" in txt
            if verbose:
                print("  ok" if ok else "  FALLÓ")
            if not ok:
                cola = "\n".join(txt.strip().split("\n")[-15:])
                raise FaltanDatos(
                    f"{d.name}/{nombre} falló. Últimas líneas:\n\n{cola}")
    return run


def collect(run, outdir="berry"):
    """Lee las fases de cada punto y sigue la rama a lo largo del camino."""
    out = Path(outdir)
    dirs = [Path(d) for d in run.jobs] if run.jobs else \
        sorted(p for p in out.glob("p[0-9][0-9]") if p.is_dir())
    if not dirs:
        raise FaltanDatos(
            f"en {out} no hay carpetas p00, p01... ¿Preparaste el cálculo "
            f"con `olla-dft berry`?")
    run.puntos = []
    for d in dirs:
        b = leer_berry(d / "2_berry.out")
        b.etiqueta = d.name
        run.puntos.append(b)
    if not run.valencias:
        try:
            run.valencias = valencias_de(dirs[0] / "1_scf.out")
        except Exception:                                   # noqa: BLE001
            run.valencias = {}
    if len(run.estructuras) != len(run.puntos):
        run.estructuras = run.estructuras[:len(run.puntos)] or []
        run.lambdas = run.lambdas[:len(run.puntos)] or \
            list(np.linspace(0, 1, len(run.puntos)))
    return run


def analizar(run):
    """Del camino de fases al número que se publica.

    Aquí está la aritmética que más fácil se hace mal:

      - El cuanto depende de MOD_TOT (ver `modulo_de`), no siempre es 2.
      - La carga efectiva de Born NO es (Ω/e)·ΔP/|u| con el módulo de P. La
        fase mide la PROYECCIÓN de P·Ω/e sobre R_gdir, así que

            Z* = 2π·Δφ / (u · B_gdir)

        con B_gdir el vector recíproco (2π·inv(celda)ᵀ). Usar |R| en vez de
        esa proyección da un factor √2 en la estructura zinc-blenda —salía
        2.46 en vez de 2.01 en BN cúbico— y el error no se nota porque el
        número sigue pareciendo razonable.
    """
    if not run.puntos:
        raise FaltanDatos("no hay puntos que analizar.")
    at = run.estructuras[-1] if run.estructuras else None
    if at is None:
        raise FaltanDatos("hace falta la estructura para convertir a C/m².")
    mod = float(run.puntos[0].modulo or 2.0)
    fases = np.array([p.fase_total for p in run.puntos])
    seguidas, saltos = desenrollar(fases, modulo=mod)
    q_eA2, q_cm2 = cuanto(at, run.gdir, mod)
    cell = np.array(at.cell.array, float)
    Rg = cell[run.gdir - 1]
    vol = float(abs(np.linalg.det(cell)))
    P = seguidas * float(np.linalg.norm(Rg)) / vol * E_A2_A_C_M2
    fuera = {"fases": fases, "seguidas": seguidas, "saltos": saltos,
             "P": P, "cuanto_cm2": q_cm2, "cuanto_eA2": q_eA2,
             "modulo": mod,
             "dP": float(P[-1] - P[0]) if len(P) > 1 else None}
    if len(saltos) and (saltos / mod).max() > FRACCION_SOSPECHOSA:
        run.avisos.append(
            "Un paso del camino mueve la fase "
            f"{(saltos / mod).max() * 100:.0f} % del cuanto. El seguimiento "
            "de la rama supone\n  que el paso es pequeño; con saltos así, "
            "elegir la imagen más cercana es una\n  apuesta. Sube --nlambda "
            "hasta que ΔP deje de cambiar.")
    if run.es_desplazamiento and run.desplazamiento is not None and len(P) > 1:
        Bg = (2 * np.pi * np.linalg.inv(cell).T)[run.gdir - 1]
        u = np.asarray(run.desplazamiento, float)
        proy = float(u @ Bg)
        if abs(proy) < 1e-8:
            run.avisos.append(
                "El desplazamiento es perpendicular a la dirección de la "
                "fase, así que esta\n  no ve nada y Z* no se puede sacar de "
                "aquí. Cambia --gdir o la dirección del\n  desplazamiento.")
        else:
            lam = np.array(run.lambdas[:len(seguidas)], float)
            if len(seguidas) > 2:
                coef = np.polyfit(lam, seguidas, 1)
                pend = float(coef[0])
                aj = np.polyval(coef, lam)
                ss = float(((seguidas - seguidas.mean()) ** 2).sum())
                fuera["zeff_r2"] = ((1.0 - float(((seguidas - aj) ** 2).sum())
                                     / ss) if ss > 0 else 1.0)
            else:
                pend = float((seguidas[-1] - seguidas[0]) /
                             (lam[-1] - lam[0]))
                fuera["zeff_r2"] = None
            fuera["zeff"] = 2 * np.pi * pend / proy
            fuera["u_proy"] = proy
    return fuera


def comprobar_ionica(run):
    """La parte iónica contra su fórmula exacta. O cuadra, o algo está mal."""
    if not run.valencias or not run.estructuras:
        return None
    fuera = []
    for at, b in zip(run.estructuras, run.puntos):
        try:
            teorica = a_intervalo(fase_ionica(at, run.valencias, run.gdir),
                                  b.modulo or 2.0)
        except FaltanDatos:
            return None
        fuera.append((b.fase_ion, float(teorica),
                      abs(a_intervalo(b.fase_ion - teorica,
                                      b.modulo or 2.0))))
    return fuera


def report(run) -> str:
    an = analizar(run)
    L = [f"--- Polarización por fase de Berry: {run.formula} ---",
         f"Dirección {run.gdir}   |   {run.nppstr} puntos por cuerda × "
         f"{run.puntos[0].ncuerdas or run.kperp[0] * run.kperp[1]} cuerdas",
         f"Cuanto de polarización: {an['cuanto_eA2']:.6f} e/Å² = "
         f"{an['cuanto_cm2']:.4f} C/m²   (MOD_TOT = {an['modulo']:g})",
         ""]
    if run.valencias:
        L.append("Valencias del pseudo: "
                 + ", ".join(f"{k} {v:g}" for k, v in
                             sorted(run.valencias.items())))
    L += ["", f"Fases (unidades de QE: el cuanto vale {an['modulo']:g}"
              + (", porque alguna valencia es impar)"
                 if an["modulo"] == 1 else ", todas las valencias son pares)"),
          "   λ      iónica   electrón.   total    seguida      P (C/m²)"]
    for i, b in enumerate(run.puntos):
        lam = run.lambdas[i] if i < len(run.lambdas) else float(i)
        L.append(f"  {lam:5.3f}  {b.fase_ion:9.5f} {b.fase_el:10.5f} "
                 f"{b.fase_total:9.5f} {an['seguidas'][i]:10.5f} "
                 f"{an['P'][i]:12.5f}")

    comp = comprobar_ionica(run)
    if comp:
        peor = max(c[2] for c in comp)
        L += ["", f"Comprobación de la parte iónica contra Σ Z_a·f_a: "
                  f"peor desviación {peor:.2e}"
                  + ("  ✓" if peor < 1e-4 else
                     "  ← MAL: la geometría o las valencias no son las que "
                     "crees")]

    if an["dP"] is not None:
        L += ["", f"ΔP a lo largo del camino ({run.etiqueta_ref} → "
                  f"{run.etiqueta_fin}):",
              f"  {an['dP']:+.5f} C/m²   "
              f"({an['dP'] / an['cuanto_cm2']:+.4f} cuantos)"]
        if abs(an["dP"]) > 0.9 * an["cuanto_cm2"]:
            L.append("  Ojo: ΔP es casi un cuanto entero. Comprueba con más "
                     "puntos que no es un salto de rama disfrazado.")
    if "zeff" in an:
        L += ["", f"Carga efectiva de Born a lo largo de la dirección "
                  f"{run.gdir}:",
              f"  Z* = {an['zeff']:+.4f} e"
              + (f"   (ajuste lineal, R² = {an['zeff_r2']:.6f})"
                 if an.get("zeff_r2") is not None else "")]
        L.append("  Z* es la carga que hay que mover para producir la "
                 "polarización observada;")
        L.append("  no es la carga iónica nominal, y en un cristal "
                 "homopolar vale cero exactamente.")
    if len(run.puntos) == 1:
        L += ["", "Con un solo punto no hay ΔP. El valor de P de arriba está "
                  "definido módulo el",
              "cuanto y no es publicable por sí mismo."]
    for a in run.avisos:
        L += ["", f"AVISO: {a}"]
    return "\n".join(L)


def export(run, outdir="berry") -> list:
    out = Path(outdir); out.mkdir(parents=True, exist_ok=True)
    an = analizar(run)
    f = out / "BERRY.dat"
    cols = np.column_stack([
        np.array(run.lambdas[:len(run.puntos)], float),
        [p.fase_ion for p in run.puntos], [p.fase_el for p in run.puntos],
        [p.fase_total for p in run.puntos], an["seguidas"], an["P"]])
    np.savetxt(f, cols, fmt="%14.7f",
               header="lambda    fase_ion    fase_el    fase_tot   "
                      "fase_seguida    P(C/m2)")
    g = out / "BERRY.txt"
    g.write_text(report(run) + "\n", encoding="utf-8")
    return [str(f), str(g)]


def polarizacion_plegada(fases_totales, an) -> np.ndarray:
    """P tal como la escribe pw.x (sin seguir la rama), en C/m².

    Una fase en unidades de QE vale `fase / MOD_TOT` cuantos, así que
    P = fase / MOD_TOT · cuanto. Es la misma aritmética que usa `analizar`
    para la rama seguida y BERRY.dat; la figura dividía por 2 a secas y con
    alguna valencia impar (MOD_TOT = 1) los marcadores de pw.x quedaban a la
    mitad de donde debían, sin coincidir con la rama seguida ni en λ = 0.
    """
    mod = float(an.get("modulo") or 2.0)
    return np.asarray(fases_totales, float) / mod * float(an["cuanto_cm2"])


def plot(run, outfile="berry", formats="pdf,png", theme=None, size=None,
         family=None, background=None, palette=None, usetex=None,
         width="single", journal="generic", aspect=0.72, mono=False,
         dpi=None) -> list:
    """P frente a λ, con el cuanto dibujado para que se vea la escala."""
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError as exc:                              # pragma: no cover
        raise RuntimeError("matplotlib no está instalado.") from exc
    an = analizar(run)
    if len(an["P"]) < 2:
        raise FaltanDatos("con un solo punto no hay curva que dibujar.")
    st = qstyle.apply(theme, size=size, family=family, background=background,
                      palette=palette, usetex=usetex, mono=mono)
    fig, ax = qstyle.new_figure(width, journal, aspect)
    cols = qstyle.palette(3, mono=mono)
    lam = np.array(run.lambdas[:len(an["P"])], float)
    ax.plot(lam, an["P"], marker="o", ms=4, lw=st["line"], color=cols[0],
            label="rama seguida")
    plegada = polarizacion_plegada(
        [p.fase_total for p in run.puntos], an)[:len(lam)]
    ax.plot(lam, plegada, ls="none", marker="x", ms=5, color=cols[1],
            label="lo que escribe pw.x (plegado)")
    ax.axhspan(an["P"][0] - an["cuanto_cm2"] / 2,
               an["P"][0] + an["cuanto_cm2"] / 2, color=cols[2], alpha=0.10,
               lw=0, zorder=0, label="un cuanto")
    ax.set_xlabel(r"$\lambda$")
    ax.set_ylabel("P (C/m$^2$)")
    ax.legend(frameon=False, fontsize=st["legend"])
    escritos = qstyle.save(fig, outfile, formats, dpi=dpi, modulo="berry")
    plt.close(fig)
    return escritos
