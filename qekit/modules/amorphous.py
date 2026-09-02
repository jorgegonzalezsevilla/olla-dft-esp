# Olla-DFT — command-line toolkit for Quantum ESPRESSO
# Copyright (C) 2026 Jorge Enrique González Sevilla
# SPDX-License-Identifier: GPL-3.0-or-later
# This program is free software: you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by the Free
# Software Foundation, either version 3 of the License, or (at your option)
# any later version. See the LICENSE file for details.

"""Sólidos amorfos por fundido y temple.

Un amorfo no se construye: se genera. El procedimiento estándar es fundir
el material muy por encima de su punto de fusión, dejar que pierda toda
memoria de la red cristalina, y enfriarlo lo bastante rápido para que no le
dé tiempo a recristalizar. Lo que queda depende de la VELOCIDAD DE TEMPLE,
y eso hay que decirlo: los temples que caben en un ordenador son de
10¹²–10¹⁴ K/s, entre seis y ocho órdenes de magnitud más rápidos que
cualquier experimento. Las estructuras que salen están más desordenadas y
menos relajadas que las reales, y sus densidades salen algo bajas.

El protocolo POR OMISIÓN (3000 → 300 K en 1000 pasos de 1 fs) es todavía
más rápido, 2.7×10¹⁵ K/s: es un protocolo de exploración, pensado para
tener un amorfo en minutos, y el reporte avisa de ello. Para bajar a
10¹⁴ K/s hacen falta ~27 000 pasos de temple (--quench-steps 27000), y a
10¹³ K/s diez veces más; el aviso desaparece por debajo de 10¹³ K/s.

Aquí el fundido y el temple se hacen con un potencial interatómico
aprendido (MACE y compañía), no con DFT. Es la única forma de que quepan
los miles de pasos que hace falta: con pw.x, mil pasos de una celda de
setenta átomos son semanas. Lo que se obtiene es un punto de partida
razonable que después se relaja con DFT, y el módulo lo dice en vez de
fingir que la estructura ya es de primeros principios.

Para el análisis (g(r), coordinación, difusión) se usa `olla-dft md`, que ya
sabe leer trayectorias.
"""

import time
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

from qekit.core import provenance
from qekit.core.errors import ErrorDeUso

# Factor sobre la suma de radios covalentes por debajo del cual dos átomos
# se consideran solapados al empaquetar. 0.75 deja sitio para que el
# potencial arregle el resto sin que la energía inicial se dispare.
FACTOR_MINIMO = 0.75


@dataclass
class Protocolo:
    T_fundido: float = 3000.0     # K
    T_final: float = 300.0
    pasos_fundido: int = 500
    pasos_temple: int = 1000
    pasos_recocido: int = 200
    dt_fs: float = 1.0

    @property
    def velocidad_temple(self) -> float:
        """K/s. Es el parámetro que decide cómo sale el amorfo."""
        t = self.pasos_temple * self.dt_fs * 1e-15
        return (self.T_fundido - self.T_final) / t if t > 0 else float("inf")

    @property
    def pasos(self) -> int:
        return self.pasos_fundido + self.pasos_temple + self.pasos_recocido

    @property
    def ps_totales(self) -> float:
        return self.pasos * self.dt_fs / 1000.0


@dataclass
class Amorfo:
    atoms: object = None
    protocolo: Protocolo = None
    modelo: str = "mace"
    densidad: float = None            # g/cm³
    energias: list = field(default_factory=list)
    temperaturas: list = field(default_factory=list)
    segundos: float = 0.0
    avisos: list = field(default_factory=list)


# ----------------------------------------------------------------------
# Empaquetado inicial
# ----------------------------------------------------------------------
def densidad_de(atoms) -> float:
    """g/cm³ de una celda."""
    masa_g = float(sum(atoms.get_masses())) * 1.66053906660e-24
    vol_cm3 = float(abs(np.linalg.det(atoms.cell.array))) * 1e-24
    return masa_g / vol_cm3 if vol_cm3 > 0 else 0.0


def celda_para_densidad(simbolos, densidad: float) -> float:
    """Arista de la celda cúbica que da esa densidad, en Å."""
    from ase.data import atomic_masses, atomic_numbers
    masa = sum(atomic_masses[atomic_numbers[s]] for s in simbolos)
    masa_g = masa * 1.66053906660e-24
    vol_cm3 = masa_g / densidad
    return float((vol_cm3 * 1e24) ** (1.0 / 3.0))


def empaquetar(simbolos, densidad: float, factor: float = FACTOR_MINIMO,
               intentos: int = 20000, semilla: int = 0):
    """Coloca los átomos al azar sin que se solapen, a la densidad pedida.

    Se rechaza cualquier posición más cercana que `factor` veces la suma de
    radios covalentes, con condiciones periódicas. Sin ese filtro la energía
    inicial se dispara y el primer paso de dinámica manda dos átomos al otro
    lado de la celda.
    """
    from ase import Atoms
    from ase.data import atomic_numbers, covalent_radii

    if densidad <= 0:
        raise ErrorDeUso(f"la densidad tiene que ser positiva; "
                         f"recibí {densidad}.")
    simbolos = list(simbolos)
    if not simbolos:
        raise ErrorDeUso("hace falta al menos un átomo.")
    L = celda_para_densidad(simbolos, densidad)
    rng = np.random.default_rng(semilla)
    radios = [covalent_radii[atomic_numbers[s]] for s in simbolos]

    pos, puestos = [], []
    fallos = 0
    for i, r_i in enumerate(radios):
        for _ in range(intentos):
            p = rng.random(3) * L
            ok = True
            for j, q in enumerate(pos):
                d = p - q
                d -= L * np.round(d / L)            # imagen mínima
                if np.linalg.norm(d) < factor * (r_i + radios[puestos[j]]):
                    ok = False
                    break
            if ok:
                pos.append(p)
                puestos.append(i)
                break
        else:
            fallos += 1
            pos.append(rng.random(3) * L)
            puestos.append(i)
    at = Atoms(symbols=simbolos, positions=np.array(pos),
               cell=np.eye(3) * L, pbc=True)
    if fallos:
        raise ErrorDeUso(
            f"no cupieron {fallos} de {len(simbolos)} átomos sin solaparse a "
            f"{densidad:g} g/cm³. O la densidad es demasiado alta para esa "
            f"composición, o hay que bajar el factor de distancia mínima "
            f"(--min-dist).")
    return at


def formula_a_simbolos(formula: str, unidades: int) -> list:
    """'SiO2' y 8 unidades -> ['Si']*8 + ['O']*16."""
    import re
    from ase.data import atomic_numbers

    piezas = re.findall(r"([A-Z][a-z]?)(\d*)", str(formula))
    piezas = [(el, int(n) if n else 1) for el, n in piezas if el]
    if not piezas:
        raise ErrorDeUso(
            f"no entiendo la fórmula '{formula}'. Se escribe como SiO2, "
            f"GeTe o Al2O3.")
    malos = [el for el, _ in piezas if el not in atomic_numbers]
    if malos:
        raise ErrorDeUso(f"elemento desconocido en la fórmula: "
                         f"{', '.join(malos)}.")
    fuera = []
    for el, n in piezas:
        fuera += [el] * (n * int(unidades))
    return fuera


# ----------------------------------------------------------------------
# Fundido y temple
# ----------------------------------------------------------------------
def fundir_y_templar(atoms, protocolo: Protocolo = None, modelo: str = "mace",
                     device: str = "cpu", semilla: int = 0,
                     traza: str = None, verbose: bool = True) -> Amorfo:
    """Funde, enfría a la velocidad pedida y relaja, con un potencial MLIP."""
    from ase import units
    from ase.md.langevin import Langevin
    from ase.md.velocitydistribution import MaxwellBoltzmannDistribution
    from qekit.modules import mlip

    p = protocolo or Protocolo()
    at = atoms.copy()
    at.calc = mlip.calculator(modelo, device=device)
    res = Amorfo(atoms=at, protocolo=p, modelo=modelo,
                 densidad=densidad_de(at))

    MaxwellBoltzmannDistribution(at, temperature_K=p.T_fundido,
                                 rng=np.random.default_rng(semilla))
    dyn = Langevin(at, p.dt_fs * units.fs, temperature_K=p.T_fundido,
                   friction=0.02, rng=np.random.default_rng(semilla))

    salida = open(traza, "w") if traza else None

    def _apuntar():
        e = at.get_potential_energy()
        T = at.get_temperature()
        res.energias.append(e)
        res.temperaturas.append(T)
        if salida:
            salida.write(f"{len(res.energias):6d} {T:10.2f} {e:16.6f}\n")

    t0 = time.time()
    dyn.attach(_apuntar, interval=10)
    if verbose:
        print(f"  fundiendo a {p.T_fundido:.0f} K "
              f"({p.pasos_fundido} pasos) ...", flush=True)
    dyn.run(p.pasos_fundido)

    if verbose:
        print(f"  templando hasta {p.T_final:.0f} K "
              f"({p.pasos_temple} pasos, "
              f"{p.velocidad_temple:.2e} K/s) ...", flush=True)
    n_tramos = max(1, p.pasos_temple // 20)
    for k in range(20):
        T = p.T_fundido + (p.T_final - p.T_fundido) * (k + 1) / 20.0
        dyn.set_temperature(temperature_K=T)
        dyn.run(n_tramos)

    if p.pasos_recocido:
        if verbose:
            print(f"  recociendo a {p.T_final:.0f} K "
                  f"({p.pasos_recocido} pasos) ...", flush=True)
        dyn.set_temperature(temperature_K=p.T_final)
        dyn.run(p.pasos_recocido)

    res.segundos = time.time() - t0
    if salida:
        salida.close()
    res.atoms = at
    res.densidad = densidad_de(at)

    # ¿siguió la temperatura a la rampa? Con fricción baja y un temple muy
    # rápido, el termostato se queda atrás y la estructura no llega a
    # solidificar: se comprobó en a-SiO2, donde un temple de 2.7e15 K/s con
    # fricción 0.02 acabó a 2257 K en vez de 300.
    if res.temperaturas:
        T_fin = float(np.mean(res.temperaturas[-3:]))
        if T_fin > p.T_final * 2.5 + 200:
            res.avisos.append(
                f"El sistema acabó a {T_fin:.0f} K, no a los {p.T_final:.0f} K "
                f"pedidos. El termostato\n  no consigue seguir una rampa tan "
                f"rápida: la estructura sigue siendo un\n  líquido, no un "
                f"vidrio. Sube --quench-steps, o alarga --anneal-steps para "
                f"que\n  termine de enfriar a temperatura fija.")

    if p.velocidad_temple > 1e13:
        res.avisos.append(
            f"Velocidad de temple {p.velocidad_temple:.1e} K/s. Un vidrio de "
            f"verdad se enfría\n  a 1-100 K/s: son diez órdenes de magnitud. "
            f"La estructura sale más\n  desordenada y menos densa que la real. "
            f"Para acercarse, alarga --quench-steps.")
    return res


# ----------------------------------------------------------------------
# Análisis
# ----------------------------------------------------------------------
def coordinaciones(atoms, corte: float = None, factor: float = 1.25) -> dict:
    """Número de coordinación medio de cada par de especies.

    El corte por defecto sale de los radios covalentes de cada par, no de un
    número único: en un óxido, el Si-O y el O-O no tienen la misma distancia
    de enlace ni de lejos, y un corte global cuenta como enlaces contactos
    que no lo son.
    """
    from ase.data import atomic_numbers, covalent_radii

    simbolos = atoms.get_chemical_symbols()
    especies = sorted(set(simbolos))
    d = atoms.get_all_distances(mic=True)
    fuera = {}
    for a in especies:
        ia = [i for i, s in enumerate(simbolos) if s == a]
        for b in especies:
            ib = [j for j, s in enumerate(simbolos) if s == b]
            if corte is not None:
                rc = corte
            else:
                rc = factor * (covalent_radii[atomic_numbers[a]]
                               + covalent_radii[atomic_numbers[b]])
            total = 0
            for i in ia:
                total += sum(1 for j in ib if j != i and d[i, j] < rc)
            fuera[(a, b)] = total / len(ia) if ia else 0.0
    return fuera


def distancia_media(atoms, a: str, b: str, factor: float = 1.25) -> float:
    """Distancia media del primer vecino entre dos especies."""
    from ase.data import atomic_numbers, covalent_radii

    simbolos = atoms.get_chemical_symbols()
    rc = factor * (covalent_radii[atomic_numbers[a]]
                   + covalent_radii[atomic_numbers[b]])
    d = atoms.get_all_distances(mic=True)
    vals = [d[i, j] for i, s in enumerate(simbolos) if s == a
            for j, t in enumerate(simbolos)
            if t == b and j != i and d[i, j] < rc]
    return float(np.mean(vals)) if vals else float("nan")


def report(res: Amorfo) -> str:
    at = res.atoms
    p = res.protocolo
    simbolos = at.get_chemical_symbols()
    L = ["--- Sólido amorfo por fundido y temple ---",
         f"Composición: {at.get_chemical_formula()} ({len(at)} átomos)",
         f"Densidad: {res.densidad:.4f} g/cm³   |   celda "
         f"{np.linalg.norm(at.cell.array[0]):.3f} Å",
         f"Potencial: {res.modelo}   |   {p.ps_totales:.2f} ps en "
         f"{res.segundos / 60:.1f} min",
         "",
         f"Protocolo: fundir a {p.T_fundido:.0f} K, templar a "
         f"{p.T_final:.0f} K, recocer",
         f"  Velocidad de temple: {p.velocidad_temple:.2e} K/s",
         ""]
    coord = coordinaciones(at)
    especies = sorted(set(simbolos))
    L.append("Coordinación media (primer vecino, corte por radios "
             "covalentes):")
    for a in especies:
        partes = [f"{b}: {coord[(a, b)]:.2f}" for b in especies]
        total = sum(coord[(a, b)] for b in especies)
        L.append(f"  {a:3s} -> {', '.join(partes)}   (total {total:.2f})")
    L.append("")
    L.append("Distancia media del primer vecino:")
    for i, a in enumerate(especies):
        for b in especies[i:]:
            d = distancia_media(at, a, b)
            if d == d:
                L.append(f"  {a}-{b}: {d:.3f} Å")

    if res.temperaturas:
        T_fin = float(np.mean(res.temperaturas[-3:]))
        L += ["", f"Temperatura final: {T_fin:.0f} K "
                  f"(objetivo {p.T_final:.0f} K)"]
    L += ["",
          "Esta estructura viene de un potencial aprendido, NO de DFT. Es un "
          "punto de\n  partida: relájala con 'olla-dft gen -p relax' antes de "
          "calcularle nada, y\n  compara varias realizaciones (--seed "
          "distintas), porque una sola no\n  representa a un amorfo."]
    for a in res.avisos:
        L.append(f"\nAVISO: {a}")
    return "\n".join(L)


def export(res: Amorfo, outdir: str = ".") -> list:
    from qekit.core import structure

    out = Path(outdir); out.mkdir(parents=True, exist_ok=True)
    cif = out / "amorfo.cif"
    structure.convert(res.atoms, str(cif))
    txt = out / "AMORFO.txt"
    txt.write_text(report(res) + "\n")
    dat = out / "AMORFO.dat"
    p = res.protocolo
    lineas = [provenance.header(
        "solido amorfo por fundido y temple",
        {"modelo": res.modelo, "densidad_g_cm3": f"{res.densidad:.4f}",
         "T_fundido_K": p.T_fundido, "T_final_K": p.T_final,
         "velocidad_temple_K_s": f"{p.velocidad_temple:.3e}",
         "ps": f"{p.ps_totales:.3f}"}),
        f"# {'muestra':>8s} {'T(K)':>10s} {'E(eV)':>16s}"]
    for i, (T, e) in enumerate(zip(res.temperaturas, res.energias), 1):
        lineas.append(f"{i:10d} {T:10.2f} {e:16.6f}")
    dat.write_text("\n".join(lineas) + "\n")
    return [str(cif), str(dat), str(txt)]
