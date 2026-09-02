# Olla-DFT — command-line toolkit for Quantum ESPRESSO
# Copyright (C) 2026 Jorge Enrique González Sevilla
# SPDX-License-Identifier: GPL-3.0-or-later
# This program is free software: you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by the Free
# Software Foundation, either version 3 of the License, or (at your option)
# any later version. See the LICENSE file for details.

"""Corrimientos de niveles de core (XPS) en aproximación de estado inicial.

`initial_state.x` de Quantum ESPRESSO calcula el corrimiento del nivel de
core de cada átomo respecto de una referencia, descompuesto en sus
contribuciones (Fermi, local, no local, iónica, core-correction, Hubbard).
Es la contraparte teórica de un espectro XPS: los corrimientos QUÍMICOS
entre átomos de la misma especie en entornos distintos.

QUÉ ES Y QUÉ NO ES
------------------
Esto es la aproximación de ESTADO INICIAL: solo tiene en cuenta el
potencial que ve el electrón de core ANTES de arrancarlo. El experimento
mide la energía de enlace, que incluye además la relajación del sistema
frente al hueco de core creado (el "estado final"). Por eso:

- los corrimientos RELATIVOS entre sitios equivalentes químicamente salen
  bien y son lo comparable con un XPS;
- las energías de enlace ABSOLUTAS no salen de aquí: para eso hace falta
  un cálculo ΔSCF con un pseudopotencial de hueco de core, y ese cálculo
  NO está implementado en este módulo (solo el estado inicial con
  initial_state.x; el UPF con hueco de core se usa aquí únicamente como
  la "especie excitada" que initial_state.x necesita para definir el
  corrimiento, no para relajar el sistema frente al hueco);
- la parte de relajación puede valer varias décimas de eV, así que un
  corrimiento calculado de 0.2 eV no es concluyente.

El módulo solo tiene sentido si hay VARIOS átomos de la misma especie en
entornos distintos: en un cristal donde todos son equivalentes por
simetría, el corrimiento es cero por construcción — y Olla-DFT lo dice en vez
de presentar una tabla de ceros como si fuera un resultado.

LO QUE HACE FALTA PARA QUE ESTO FUNCIONE
----------------------------------------
`initial_state.x` NO se activa con una bandera. Su variable `excite(nt)`
es el ÍNDICE DE OTRA ESPECIE ATÓMICA: la del mismo elemento pero con un
pseudopotencial de HUECO DE CORE (un electrón menos en el core, y por
tanto z_valence una unidad mayor). El código calcula el corrimiento a
partir de la diferencia entre las dos especies:

    delta_zv(nt) = zv(excite(nt)) - zv(nt)

Es decir, el input de pw.x tiene que declarar DOS especies para el mismo
elemento — la normal y la excitada — y `excite(1) = 2` significa "la
contraparte excitada del tipo 1 es el tipo 2".

Si se pone `excite(1) = 1`, delta_zv sale cero y el programa devuelve una
tabla entera de ceros SIN dar ningún error. Olla-DFT se niega a generar esa
configuración por eso: es un cero que parece un resultado.

Los pseudopotenciales de hueco de core no vienen en las tablas estándar;
hay que generarlos con `ld1.x` (por ejemplo, un Si con ocupación 1s¹ en
vez de 1s²) o conseguirlos ya hechos.
"""

import re
import shutil
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

from qekit.core import provenance
from qekit.core.errors import ErrorDeUso, FaltanDatos
from qekit.modules import sweep

_RE_SECCION = re.compile(r"The\s+(.+?)\s+contribution to shift", re.I)
_RE_ATOMO = re.compile(
    r"atom\s+(\d+)\s+type\s+(\d+)\s+shift\s*=\s*(-?[\d.]+)\s*Ry,\s*="
    r"\s*(-?[\d.]+)\s*eV")


@dataclass
class XPSResult:
    shifts: np.ndarray = None            # (nat,) eV, contribución total
    contributions: dict = field(default_factory=dict)   # nombre -> (nat,) eV
    types: np.ndarray = None
    symbols: list = field(default_factory=list)
    equivalentes: bool = False


def build_input(prefix: str, pares: dict) -> str:
    """Input de initial_state.x.

    `pares` mapea el índice del tipo en estado fundamental al índice del
    tipo con hueco de core, ambos base 1 y en el orden de ATOMIC_SPECIES:
    {1: 2} significa "la contraparte excitada del tipo 1 es el tipo 2".
    """
    if not pares:
        raise ErrorDeUso(
            "initial_state.x necesita al menos un par "
            "{tipo_fundamental: tipo_excitado}. Sin él, delta_zv vale cero "
            "y el programa devuelve una tabla de ceros sin avisar.")
    for gs, ex in pares.items():
        if int(gs) == int(ex):
            raise ErrorDeUso(
                f"excite({gs}) = {ex}: un tipo no puede ser su propia "
                "contraparte excitada. El tipo excitado es OTRA especie, "
                "la del mismo elemento con pseudopotencial de hueco de "
                "core.")
    lines = ["&inputpp", f"  prefix = '{prefix}'", "  outdir = './out'"]
    for gs, ex in sorted(pares.items()):
        lines.append(f"  excite({int(gs)}) = {int(ex)}")
    lines.append("/")
    return "\n".join(lines) + "\n"




def _verificar_par(elem: str, common: dict, upf_hueco, pseudo_dir) -> None:
    """El par normal/hueco tiene que diferir en exactamente un electrón.

    Es LA comprobación del módulo. `initial_state.x` calcula el
    corrimiento como zv(excitada) - zv(normal); si los dos pseudos son el
    mismo archivo, o difieren en otra cosa, la diferencia no es un hueco
    de core y el programa devuelve números sin quejarse.
    """
    from qekit.core import pseudo as ps

    base_nombre = common["pseudos"].get(elem, {}).get("filename")
    hueco = Path(upf_hueco)
    if base_nombre == hueco.name:
        raise ErrorDeUso(
            f"el pseudopotencial normal de {elem} y el de hueco de core son "
            f"el MISMO archivo ({hueco.name}). Los dos tipos quedarían "
            "idénticos y initial_state.x devolvería ceros. Usa --pseudo-dir "
            "con el pseudo normal, o genera el par con:\n"
            f"  olla-dft corehole {elem} --edge K")

    base_ruta = Path(pseudo_dir) / (base_nombre or "")
    zb = ps.z_valence(base_ruta) if base_ruta.is_file() else None
    zh = ps.z_valence(hueco) if hueco.is_file() else None
    if zb is None or zh is None:
        return
    if abs((zh - zb) - 1.0) > 1e-6:
        raise ErrorDeUso(
            f"z_valence de {elem}: {zb:g} (normal) y {zh:g} (hueco); la "
            f"diferencia es {zh - zb:+g} y tiene que ser exactamente +1.\n"
            "Un hueco de core es UN electrón menos en el core, así que el "
            "pseudo excitado lleva una carga de valencia más. Con otra "
            "diferencia el corrimiento que salga no es un corrimiento de "
            "nivel de core.\n"
            f"Genera el par consistente con:  olla-dft corehole {elem} --edge K")

def _copiar_pseudos(core_hole: dict, pseudo_dir) -> None:
    """Deja los UPF con hueco donde pw.x los va a buscar.

    Un pseudo recién generado suele estar en la carpeta de trabajo, no en
    la de pseudopotenciales. Copiarlo aquí evita el fallo más tonto de
    todos: 'cannot open file' después de esperar el scf.
    """
    destino = Path(pseudo_dir or ".")
    if not destino.is_dir():
        return
    for upf in core_hole.values():
        origen = Path(upf)
        if not origen.is_file():
            continue
        meta = destino / origen.name
        if meta.resolve() == origen.resolve():
            continue
        try:
            shutil.copy2(origen, meta)
        except OSError:
            pass

def prepare(atoms, outdir: str = "xps", pseudo_dir: str = None,
            ecutwfc: float = None, ecutrho: float = None,
            kspacing: float = None, insulator: bool = True,
            core_hole: dict = None, sites=None) -> tuple:
    """Escribe scf.in (con la especie excitada añadida) e initial_state.in.

    `core_hole` mapea elemento -> archivo UPF con hueco de core, por
    ejemplo {"Si": "Si.star1s.UPF"}. `sites` son los índices (base 0) de
    los átomos que se marcan como la especie excitada; si no se dan, se
    marca uno por cada sitio inequivalente.
    """
    from qekit.modules import inputgen

    common = sweep.prepare_common(
        atoms, pseudo_dir, ecutwfc, ecutrho, insulator,
        exclude_pseudos=list((core_hole or {}).values()))
    out = Path(outdir); out.mkdir(parents=True, exist_ok=True)
    grid = sweep.default_grid(atoms, kspacing)
    especies = list(dict.fromkeys(atoms.get_chemical_symbols()))

    # Las especies con hueco de core se declaran en ATOMIC_SPECIES aunque
    # ningún átomo las use: son el TIPO al que apunta excite(nt). Sin esto
    # el input decía excite(1)=2 y el tipo 2 no existía.
    extra, pares = [], {}
    for elem, upf in (core_hole or {}).items():
        if elem not in especies:
            raise ErrorDeUso(
                f"'{elem}' no está en la estructura ({', '.join(especies)}); "
                "el pseudo con hueco de core tiene que ser del mismo "
                "elemento que quieres excitar.")
        from qekit.modules.xanes import etiqueta_excitada
        etiqueta = etiqueta_excitada(elem)
        extra.append((etiqueta, elem, Path(upf).name))
        pares[especies.index(elem) + 1] = len(especies) + len(extra)
        _verificar_par(elem, common, upf, Path(common["pseudo_dir"]))

    scf = inputgen.build_pw_input(
        atoms=atoms, pseudos=common["pseudos"], calculation="scf",
        prefix=common["prefix"], pseudo_dir=common["pseudo_dir"],
        ecutwfc=common["ecutwfc"], ecutrho=common["ecutrho"],
        kcard=f"K_POINTS automatic\n  {grid[0]} {grid[1]} {grid[2]} 0 0 0\n",
        insulator=insulator, degauss=common["degauss"],
        smearing=common["smearing"], extra_species=extra)
    sweep.write_input(out / "scf.in", scf)
    if pares:
        sweep.write_input(out / "initial_state.in",
                          build_input(common["prefix"], pares))
        _copiar_pseudos(core_hole, common["pseudo_dir"])

    # ¿hay sitios inequivalentes? si no, el resultado será cero
    from qekit.core import structure as st
    try:
        ds = st.symmetry_dataset(atoms)
        n_orbitas = len(set(ds.equivalent_atoms))
    except Exception:
        n_orbitas = len(atoms)

    rep = ["--- Corrimientos de nivel de core (XPS) ---",
           f"Estructura: {atoms.get_chemical_formula()} "
           f"({len(atoms)} átomos)",
           f"Sitios inequivalentes por simetría: {n_orbitas}",
           "",
           f"Archivos en '{out.resolve()}': scf.in, initial_state.in",
           "Orden: pw.x -in scf.in  ->  initial_state.x -in "
           "initial_state.in"]
    if not core_hole:
        rep += ["",
                "NO se escribió initial_state.in: falta el pseudopotencial "
                "de HUECO DE CORE.",
                "initial_state.x compara la especie normal contra otra del "
                "mismo elemento con",
                "un electrón menos en el core; sin esa segunda especie "
                "devuelve una tabla de",
                "ceros sin dar error. Pásalo con --core-hole "
                "Si=Si.star1s.UPF (se genera con",
                "ld1.x, no viene en las tablas estándar)."]
    if n_orbitas <= 1:
        rep += ["",
                "AVISO: todos los átomos son equivalentes por simetría, así "
                "que todos los\ncorrimientos van a salir exactamente cero. "
                "Para un corrimiento químico\nhacen falta átomos de la misma "
                "especie en entornos distintos: una\nsuperficie, un defecto, "
                "un dopante o un compuesto con varios sitios."]
    rep += ["",
            "Recuerda: esto es la aproximación de ESTADO INICIAL. Sirve para "
            "los\ncorrimientos relativos, no para energías de enlace "
            "absolutas."]
    warn = sweep.missing_pseudo_warning(common)
    if warn:
        rep.append(warn)
    return common, "\n".join(rep)


def collect(path, symbols=None, tol: float = 1e-6) -> XPSResult:
    """Lee la salida de initial_state.x."""
    texto = Path(path).read_text(errors="ignore").splitlines()
    res = XPSResult(contributions={})
    seccion = "TOTAL"
    acumulado = {}
    for linea in texto:
        m = _RE_SECCION.search(linea)
        if m:
            seccion = m.group(1).strip().upper()
            continue
        a = _RE_ATOMO.search(linea)
        if a:
            iat, ityp, ev = (int(a.group(1)), int(a.group(2)),
                             float(a.group(4)))
            acumulado.setdefault(seccion, {})[iat] = (ityp, ev)
    if not acumulado:
        raise FaltanDatos(
            f"no hay corrimientos en '{path}'; ¿corrió initial_state.x?")

    for nombre, datos in acumulado.items():
        idx = sorted(datos)
        res.contributions[nombre] = np.array([datos[i][1] for i in idx])
        if res.types is None:
            res.types = np.array([datos[i][0] for i in idx])
    res.shifts = res.contributions.get(
        "TOTAL", list(res.contributions.values())[0])
    if symbols is not None:
        res.symbols = list(symbols)
    res.equivalentes = bool(np.all(np.abs(res.shifts) < tol))
    return res


def report(res: XPSResult) -> str:
    lines = ["--- Corrimientos de nivel de core (estado inicial) ---"]
    if res.equivalentes:
        lines += ["",
                  "Todos los corrimientos son cero: los átomos son "
                  "equivalentes por simetría.",
                  "No es un fallo del cálculo, es lo que la simetría exige. "
                  "Para ver un\ncorrimiento químico hace falta una "
                  "estructura con sitios inequivalentes\n(superficie, "
                  "defecto, dopante, o un compuesto con varios entornos)."]
        return "\n".join(lines)

    ref = float(np.min(res.shifts))
    lines += ["", f"{'átomo':>7s} {'especie':>8s} {'shift(eV)':>11s} "
                  f"{'rel. al mín.':>13s}"]
    for i, sh in enumerate(res.shifts):
        sym = res.symbols[i] if i < len(res.symbols) else \
            (f"tipo {res.types[i]}" if res.types is not None else "?")
        lines.append(f"{i+1:7d} {sym:>8s} {sh:11.4f} {sh - ref:+13.4f}")

    rango = float(np.ptp(res.shifts))
    lines += ["",
              f"Dispersión total: {rango:.3f} eV"]
    if rango < 0.1:
        lines.append(
            "  Por debajo de ~0.1 eV el corrimiento no es concluyente: la "
            "relajación de\n  estado final, que esta aproximación no "
            "incluye, es del mismo orden.")
    mayor = 0.0
    if len(res.contributions) > 1:
        lines += ["", "Descomposición (rango en eV por contribución):"]
        for nombre, v in res.contributions.items():
            if nombre == "TOTAL":
                continue
            r = float(np.ptp(v))
            mayor = max(mayor, r)
            lines.append(f"  {nombre.lower():18s} {r:8.4f}")

    # El corrimiento sale de la resta de dos números enormes que casi se
    # cancelan (local contra iónico). Cuando la razón es grande, la
    # precisión del resultado depende de que el scf esté bien convergido
    # mucho más de lo que sugiere el número final.
    if rango > 0 and mayor / rango > 20:
        lines += ["",
                  f"CUIDADO con la cancelacion: la contribucion mayor abarca "
                  f"{mayor:.1f} eV y el corrimiento final es de {rango:.2f} eV "
                  f"-- una cancelacion de 1 en {mayor / rango:.0f}.",
                  "El resultado hereda el error del scf amplificado por ese "
                  "factor: baja conv_thr (1e-10 o menos) y sube la malla k "
                  "antes de creerte la tercera cifra."]
    lines += ["",
              "Aproximación de estado inicial: los corrimientos RELATIVOS "
              "son lo comparable\ncon un XPS; las energías de enlace "
              "absolutas necesitan un ΔSCF con hueco\nde core."]
    return "\n".join(lines)


def export(res: XPSResult, outdir: str = ".") -> list:
    out = Path(outdir); out.mkdir(parents=True, exist_ok=True)
    f = out / "XPS_CORE.dat"
    lines = [provenance.header("corrimientos de core (XPS)",
                               titulo="Corrimientos de nivel de core"),
             f"# {'atomo':>6s} {'especie':>8s} " +
             " ".join(f"{k.lower():>14s}" for k in res.contributions)]
    for i in range(len(res.shifts)):
        sym = res.symbols[i] if i < len(res.symbols) else "?"
        vals = " ".join(f"{res.contributions[k][i]:14.6f}"
                        for k in res.contributions)
        lines.append(f"{i+1:8d} {sym:>8s} {vals}")
    f.write_text("\n".join(lines) + "\n")
    return [str(f)]
