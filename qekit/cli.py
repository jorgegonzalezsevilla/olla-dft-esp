# Olla-DFT — command-line toolkit for Quantum ESPRESSO
# Copyright (C) 2026 Jorge Enrique González Sevilla
# SPDX-License-Identifier: AGPL-3.0-or-later
# This program is free software: you can redistribute it and/or modify it
# under the terms of the GNU Affero General Public License as published by the Free
# Software Foundation, either version 3 of the License, or (at your option)
# any later version. See the LICENSE file for details.

"""CLI de Olla-DFT.

Dos modos de uso:
  olla-dft                  -> menú interactivo para principiantes
  olla-dft <comando> ...    -> uso directo (scriptable)

Comandos de pre-proceso:  gen, info, kpath, prim, conv, supercell, convert
Comandos de post-proceso: bands, dos, gap, plot
Módulos de cálculo:       converge, eos, elastic
Otros:                    config
"""

import argparse
import json
import os
import re
import sys
import tempfile
import textwrap

import numpy as np
from pathlib import Path

from qekit import __command_name__, __product_name__, __version__
from qekit import config as qcfg
from qekit.core import i18n, kpoints, structure
from qekit.core.errors import ErrorDeUso


def _prepare_matplotlib_cache() -> None:
    """Evita avisos ruidosos cuando la carpeta de usuario es de solo lectura."""
    if os.environ.get("MPLCONFIGDIR"):
        return
    default = Path.home() / ".config" / "matplotlib"
    # No crear una carpeta en el home solo por ejecutar ``olla-dft --version``.
    # Si ya existe y es escribible, o si su padre permite crearla, Matplotlib
    # puede usarla normalmente.
    if ((default.is_dir() and os.access(default, os.W_OK)) or
            (not default.exists() and default.parent.is_dir() and
             os.access(default.parent, os.W_OK))):
        return
    fallback = Path(tempfile.gettempdir()) / "olla-dft-matplotlib"
    try:
        fallback.mkdir(parents=True, exist_ok=True)
        os.environ["MPLCONFIGDIR"] = str(fallback)
    except OSError:
        # Matplotlib puede funcionar sin caché; no convertir una preferencia
        # de rendimiento en un fallo de CLI.
        pass


_prepare_matplotlib_cache()
from qekit.modules import inputgen


class _Perezoso:
    """Módulo que se importa en el primer uso.

    matplotlib (estilo y temas), strain y defects no hacen falta para la mayoría
    de los comandos; importarlos siempre costaba ~0.1 s por invocación.
    """

    def __init__(self, nombre):
        self._nombre = nombre
        self._mod = None

    def __getattr__(self, attr):
        if self._mod is None:
            import importlib
            self._mod = importlib.import_module(self._nombre)
        return getattr(self._mod, attr)


qstyle = _Perezoso("qekit.core.style")
qthemes = _Perezoso("qekit.core.themes")
strain_mod = _Perezoso("qekit.modules.strain")
defects_mod = _Perezoso("qekit.modules.defects")

WAVE_CHOICES = ("CuKa", "CuKa1", "CoKa", "MoKa", "FeKa", "CrKa", "AgKa")
_MENU_I18N_DIR = Path(__file__).resolve().parent / "data" / "i18n"


def _menu_labels(language="es") -> dict:
    """Carga las etiquetas del menú inicial desde un archivo independiente."""
    if language not in ("es", "en"):
        raise ErrorDeUso("language debe ser es o en")
    target = _MENU_I18N_DIR / f"menu_{language}.json"
    try:
        labels = json.loads(target.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise ErrorDeUso(f"no se pudo cargar el idioma {language}: {exc}") from None
    if not isinstance(labels, dict) or not isinstance(labels.get("items"), list):
        raise ErrorDeUso(f"el menú {language} no tiene un formato válido")
    return labels


def _menu_section(language, name):
    return _menu_labels(language)["submenus"][name]

# Títulos de los grupos de opciones que se repiten en decenas de comandos.
# Argparse los lista al final de la ayuda, después de las opciones propias
# de cada comando, que es donde interesa que estén.
GRUPO_FIGURA = "figura"
GRUPO_EJECUCION = "ejecución"
GRUPO_DFT = "parámetros DFT"

# Nombres alternativos en inglés de los subcomandos con nombre en español.
ALIASES = {"system": "sistema", "recipes": "recetas", "theory": "teoria",
           "actualizar": "update"}

PRESET_MENU = {
    "101": "scf",
    "102": "relax",
    "103": "vc-relax",
    "104": "nscf",
    "105": "bands",
    "106": "dos",
    "107": "all",
}


# El CLI conserva subcomandos planos (son fáciles de automatizar), pero la
# ayuda los presenta por tareas. Cada comando debe aparecer exactamente una
# vez: tests/test_cli_catalog.py protege esta tabla contra olvidos al crecer.
COMMAND_GROUPS = (
    ("Primeros pasos", ("start", "wizard", "recetas", "teoria", "docs",
                         "sistema", "selftest", "update")),
    ("Estructuras e inputs", ("gen", "info", "kpath", "prim", "conv",
                              "supercell", "convert")),
    ("Estructura electrónica", ("bands", "dos", "plot", "gap", "fermi",
                                "effmass", "wannier", "unfold", "topology",
                                "hubbard")),
    ("Espectros y respuesta", ("optics", "tddft", "xanes", "xps",
                               "corehole", "charge", "charges", "wf",
                               "berry")),
    ("Fonones, transporte y temperatura", ("phonons", "elph", "transport",
                                            "ballistic", "kappa", "qha",
                                            "thermochem", "md", "derived")),
    ("Mecánica y estabilidad", ("converge", "eos", "elastic", "strain",
                                "layers", "xrd", "exfoliate", "gamma")),
    ("Superficies, defectos y química", ("surface", "defect", "interface",
                                         "adsorb", "eform", "align", "esm",
                                         "echem", "neb", "amorphous")),
    ("Automatización y calidad", ("doctor", "audit", "crosscheck", "cost",
                                  "db", "hull", "mlip", "suggest",
                                  "datasheet", "report", "compare", "tune",
                                  "results", "campaign", "pseudos")),
    ("Proyecto", ("project", "resilient")),
    ("Apariencia y configuración", ("templates", "config")),
)


def _catalog_text(language="es") -> str:
    labels = _menu_labels(language)
    lines = [labels["catalog_title"]]
    for title, commands in COMMAND_GROUPS:
        display_title = labels["catalog_groups"].get(title, title)
        lines += [f"\n  {display_title}:", textwrap.fill(
            "  ".join(commands), width=78,
            initial_indent="    ", subsequent_indent="    ")]
    lines += ["", labels["catalog_instruction"], labels["catalog_example"]]
    return "\n".join(lines)


# ======================================================================
# Modo directo (subcomandos)
# ======================================================================
def _add_gen_parser(sub):
    p = sub.add_parser("gen", help="generar inputs de pw.x y post-proceso")
    p.add_argument("file", help="estructura (CIF, POSCAR, entrada de pw.x, ...)")
    p.add_argument(
        "-p", "--preset", default="scf", choices=inputgen.PRESETS,
        help="tipo de cálculo (default: scf)",
    )
    p.add_argument("-o", "--outdir", default=".", help="carpeta de salida")
    p.add_argument(
        "-k", "--klevel", choices=sorted(kpoints.KSPACING_LEVELS),
        help="densidad de la malla k (gamma/coarse/medium/fine/very-fine)",
    )
    p.add_argument("--kspacing", type=float, help="espaciado k en Å^-1 (anula --klevel)")
    p.add_argument("--kgrid", type=int, nargs=3, metavar="N",
                   help="malla k explícita para scf/relax (tres enteros; anula --kspacing y --klevel)")
    p.add_argument("--band-points", type=int, help="puntos por segmento del k-path")
    p.add_argument("--ecutwfc", type=float, help="cutoff de funciones de onda (Ry)")
    p.add_argument("--ecutrho", type=float, help="cutoff de densidad (Ry)")
    p.add_argument(
        "--insulator", action="store_true",
        help="occupations='fixed' (aislantes; default: smearing)",
    )
    p.add_argument(
        "--primitive", action="store_true",
        help="reducir a la celda primitiva estandarizada antes de generar",
    )
    p.add_argument("--pseudo-dir", help="carpeta de pseudopotenciales (anula config)")
    p.add_argument("--pseudo", action="append", metavar="EL=UPF",
                   help="forzar un pseudopotencial concreto, por ejemplo Fe=Fe.rel-pbe.UPF. Se puede repetir. Sin esto, Olla-DFT elige con 'olla-dft pseudos'")
    p.add_argument("--prefix", help="prefix del cálculo (default: fórmula)")
    p.add_argument("--nspin", type=int, default=1, choices=[1, 2],
                   help="2 activa la polarización de espín")
    p.add_argument("--mag",
                   help="magnetización inicial: un número (0.5) o por elemento "
                        "(Fe=0.7,O=0). Implica --nspin 2")
    p.add_argument("--vdw", default=None, choices=list(inputgen.VDW),
                   help="corrección de dispersión (van der Waals)")
    p.add_argument("--soc", action="store_true",
                   help="acoplamiento espín-órbita: cálculo no colineal con "
                        "lspinorb (exige pseudos totalmente relativistas)")
    p.add_argument("--hubbard", metavar="EL=U", action="append",
                   help="U de Hubbard en eV por elemento, por ejemplo Ni=4.1. "
                        "Se puede repetir. Para calcularlo en vez de "
                        "suponerlo:  olla-dft hubbard --cycle")
    p.add_argument("--hubbard-style", default="legacy",
                   choices=["legacy", "card"],
                   help="legacy = lda_plus_u (QE <= 7.0), card = tarjeta "
                        "HUBBARD (QE >= 7.1)")
    p.add_argument("--charge", type=float, default=None, metavar="Q",
                   help="carga total de la celda (tot_charge): +1 le quita un "
                        "electrón, -1 se lo añade")
    p.add_argument("--dipole", nargs="?", const=3, type=int, default=None,
                   choices=[1, 2, 3], metavar="EJE",
                   help="corrección dipolar para losas polares; sin valor usa "
                        "el eje c. Coloca el diente de sierra dentro del vacío")
    p.add_argument("--nosym", action="store_true",
                   help="desactivar la simetría (nosym y noinv)")
    p.add_argument("--functional", default=None,
                   choices=sorted(inputgen.HIBRIDOS),
                   help="funcional híbrido: hse, pbe0, b3lyp o gaupbe. "
                        "Cuesta entre uno y dos órdenes de magnitud más que "
                        "PBE, y el reporte lo dice con números")
    p.add_argument("--exx-grid", dest="exx_grid", metavar="NxNxN",
                   help="malla q del intercambio exacto (default 1x1x1). "
                        "Tiene que dividir la malla de k")
    p.add_argument("--exx-fraction", dest="exx_fraction", type=float,
                   default=None, help="fracción de intercambio exacto, si "
                                      "quieres cambiar la del funcional")
    # dinámica molecular (preset md)
    p.add_argument("--dt", type=float, default=1.0, metavar="FS",
                   help="paso de tiempo de la MD en fs (default: 1.0)")
    p.add_argument("--nstep", type=int, default=1000,
                   help="pasos de la MD (default: 1000)")
    p.add_argument("--thermostat", default="none",
                   choices=["none", "rescaling", "berendsen", "andersen",
                            "initial", "reduce-history"],
                   help="termostato de la MD; none = NVE (default)")
    p.add_argument("-T", "--temperature", type=float, default=300.0,
                   help="temperatura objetivo de la MD en K (default: 300)")


def _malla(texto: str, nombre: str = "--grid") -> tuple:
    """Convierte '8x8x8' en (8, 8, 8), con un error legible si no lo es.

    Sin esto, un '1x2' llega hasta el fondo del modulo y revienta con un
    IndexError sin decir que faltaba un numero. Lo detecto con el propio
    registro de incidencias de Olla-DFT.
    """
    if not texto:
        return None
    partes = str(texto).lower().replace(",", "x").replace(" ", "x").split("x")
    partes = [q for q in partes if q != ""]
    if len(partes) != 3:
        raise ErrorDeUso(
            f"{nombre} necesita TRES numeros separados por x, por ejemplo "
            f"8x8x8; recibi '{texto}' ({len(partes)} valor"
            f"{'es' if len(partes) != 1 else ''}).")
    try:
        vals = tuple(int(q) for q in partes)
    except ValueError:
        raise ErrorDeUso(
            f"{nombre} solo admite numeros enteros; recibi '{texto}'."
        ) from None
    if any(v < 1 for v in vals):
        raise ErrorDeUso(f"{nombre} debe ser positiva; recibi '{texto}'.")
    return vals


def _malla_2d(texto: str, nombre: str = "--grid") -> tuple:
    """Convierte ``40x40`` en una malla periódica bidimensional."""
    partes = str(texto).lower().replace(",", "x").replace(" ", "x").split("x")
    partes = [q for q in partes if q]
    if len(partes) != 2:
        raise ErrorDeUso(
            f"{nombre} necesita DOS enteros separados por x, por ejemplo "
            f"40x40; recibí '{texto}'.")
    try:
        vals = tuple(int(q) for q in partes)
    except ValueError:
        raise ErrorDeUso(
            f"{nombre} solo admite enteros; recibí '{texto}'.") from None
    if min(vals) < 3:
        raise ErrorDeUso(f"{nombre} debe ser de al menos 3x3.")
    return vals


# Opciones cuyo valor puede empezar por '-'. argparse las rompe: ve el guion
# y decide que '-4:4:5' es otra bandera, no el valor de --range. Se pegan
# antes de parsear para que '-r -5:5:11' funcione igual que '-r=-5:5:11',
# que es como cualquiera lo escribe la primera vez.
_OPC_NEGATIVAS = ("-r", "--range", "--span", "--delta", "--shift", "--emin",
                  "--emax", "--vmin", "--vmax", "-q", "--charges", "--dv",
                  "--charge", "--window", "--frozen", "--her", "--temps",
                  "--fixed", "--fermi")


def _pegar_negativos(argv: list) -> list:
    out, i = [], 0
    argv = list(argv)
    while i < len(argv):
        tok = argv[i]
        if (tok in _OPC_NEGATIVAS and i + 1 < len(argv)
                and re.match(r"^-[\d.]", argv[i + 1])):
            out.append(f"{tok}={argv[i + 1]}")
            i += 2
            continue
        out.append(tok)
        i += 1
    return out


def _duracion(texto) -> float:
    """'90m', '2h', '1h30m', '3600' -> segundos."""
    if texto is None:
        return None
    t = str(texto).strip().lower().replace(" ", "")
    if not t:
        return None
    m = re.fullmatch(r"(?:(\d+(?:\.\d+)?)h)?(?:(\d+(?:\.\d+)?)m)?"
                     r"(?:(\d+(?:\.\d+)?)s)?", t)
    if m and any(m.groups()):
        h, mi, se = (float(g) if g else 0.0 for g in m.groups())
        total = h * 3600 + mi * 60 + se
    else:
        try:
            total = float(t)
        except ValueError:
            raise ErrorDeUso(
                f"--max-time se escribe como 90m, 2h, 1h30m o un número de "
                f"segundos; recibí '{texto}'.") from None
    if total <= 0:
        raise ErrorDeUso(f"--max-time tiene que ser positivo; recibí '{texto}'.")
    return total


def _validar_ejecucion(args) -> None:
    """Revisa las banderas de ejecución ANTES de preparar nada.

    Si se validan dentro de _run_or_explain, un `--max-time manana` sin
    `--run` no llega nunca a comprobarse: el comando escribe los inputs,
    termina en 0, y la errata solo aparece media hora después cuando por fin
    se lanza el barrido. Se comprueban aquí, junto al resto del uso.
    """
    _duracion(getattr(args, "max_time", None))
    j = getattr(args, "jobs", None)
    if j is not None and j < 1:
        raise ErrorDeUso(f"--jobs es cuántos cálculos correr a la vez, así "
                         f"que al menos 1; recibí {j}.")
    n = getattr(args, "nproc", None)
    if n is not None and n < 1:
        raise ErrorDeUso(f"--nproc es cuántos procesos MPI por cálculo, así "
                         f"que al menos 1; recibí {n}.")
    t = getattr(args, "timeout", None)
    if t is not None and t <= 0:
        raise ErrorDeUso(f"--timeout es el límite en segundos por cálculo y "
                         f"tiene que ser positivo; recibí {t}.")
    if j and j > 1 and not getattr(args, "run", False):
        # no es un error: los inputs se escriben igual y ./run.sh N sirve
        pass


def _parse_hubbard(valores) -> dict:
    """Convierte ['Ni=4.1', 'O=0'] en {'Ni': 4.1, 'O': 0.0}."""
    if not valores:
        return None
    out = {}
    for item in valores:
        for trozo in str(item).replace(";", ",").split(","):
            trozo = trozo.strip()
            if not trozo:
                continue
            if "=" not in trozo:
                raise ErrorDeUso(
                    f"--hubbard se escribe ELEMENTO=U, por ejemplo Ni=4.1; "
                    f"recibí '{trozo}'.")
            sym, _, val = trozo.partition("=")
            sym = sym.strip().capitalize()
            try:
                out[sym] = float(val)
            except ValueError:
                raise ErrorDeUso(
                    f"el U de {sym} tiene que ser un número en eV; recibí "
                    f"'{val.strip()}'.") from None
    return out or None


def _print_prepare(rep: str) -> None:
    """Imprime el reporte de preparación, salvo en modo --collect.

    En --collect no se escribió nada: el reporte describiría un cálculo
    que no es el que corrió (el usuario pudo usar otros parámetros), así
    que se calla y se avisa de dónde se están leyendo los resultados.
    """
    from qekit.modules import sweep as _sweep
    if _sweep.writing_inputs():
        print(rep)
    else:
        print("Modo --collect: se leen los resultados ya existentes "
              "(los inputs no se tocan).")


def _cmd_gen(args) -> int:
    cfg = qcfg.load()
    atoms = structure.load(args.file)
    kspacing = float(cfg["kspacing"])
    if args.klevel:
        kspacing = kpoints.KSPACING_LEVELS[args.klevel]
        if kspacing is None:
            kspacing = 0.0  # solo Γ
    if args.kspacing is not None:
        kspacing = args.kspacing

    nspin = getattr(args, "nspin", 1)
    mag_arg = getattr(args, "mag", None)
    magnetization = {}
    if mag_arg:
        magnetization = inputgen.parse_magnetization(
            mag_arg, atoms.get_chemical_symbols()
        )
        nspin = 2  # indicar magnetización implica activar el espín

    hubbard = _parse_hubbard(getattr(args, "hubbard", None))
    if getattr(args, "soc", False) and nspin == 2:
        raise ErrorDeUso(
            "--soc y --nspin 2 no se combinan: el espín-órbita ya es un "
            "cálculo no colineal (noncolin), donde el espín no se separa en "
            "dos canales. Usa --soc solo, o --nspin 2 sin --soc.")
    md = None
    if args.preset == "md":
        md = dict(dt_fs=args.dt, nstep=args.nstep,
                  thermostat=args.thermostat, temperature=args.temperature)
    elif any(getattr(args, k, None) != d for k, d in
             (("dt", 1.0), ("nstep", 1000), ("thermostat", "none"))):
        raise ErrorDeUso(
            "--dt, --nstep y --thermostat solo tienen sentido con el preset de "
            "dinámica molecular. Añade  -p md.")

    opts = inputgen.GenOptions(
        preset=args.preset,
        outdir=args.outdir,
        kspacing=kspacing,
        kgrid=tuple(args.kgrid) if getattr(args, "kgrid", None) else None,
        kspacing_nscf=float(cfg["kspacing_nscf"]),
        band_points=args.band_points or int(cfg["band_points"]),
        ecutwfc=args.ecutwfc,
        ecutrho=args.ecutrho,
        insulator=args.insulator,
        use_primitive=True if args.primitive else None,
        pseudo_dir=args.pseudo_dir,
        prefix=args.prefix,
        nspin=nspin,
        magnetization=magnetization,
        vdw=args.vdw,
        soc=getattr(args, "soc", False),
        hubbard=hubbard,
        hubbard_style=getattr(args, "hubbard_style", "legacy"),
        tot_charge=args.charge,
        dipole=args.dipole if args.dipole else False,
        hibrido=getattr(args, "functional", None),
        exx_grid=_malla(getattr(args, "exx_grid", None), "--exx-grid"),
        exx_fraction=getattr(args, "exx_fraction", None),
        nosym=getattr(args, "nosym", False),
        md=md,
    )
    print(inputgen.generate(atoms, opts))
    return 0


def _cmd_info(args) -> int:
    atoms = structure.load(args.file)
    print(structure.info_text(atoms))
    return 0


def _cmd_kpath(args) -> int:
    atoms = structure.load(args.file)
    kp = kpoints.get_kpath(atoms)
    print(kpoints.kpath_text(kp))
    return 0


def _cmd_prim(args) -> int:
    atoms = structure.load(args.file)
    prim = structure.primitive(atoms)
    out = structure.convert(prim, args.output)
    print(f"Celda primitiva ({len(prim)} átomos) escrita en: {out}")
    return 0


def _cmd_conv(args) -> int:
    atoms = structure.load(args.file)
    conv = structure.conventional(atoms)
    out = structure.convert(conv, args.output)
    print(f"Celda convencional ({len(conv)} átomos) escrita en: {out}")
    return 0


def _cmd_supercell(args) -> int:
    atoms = structure.load(args.file)
    sc = structure.supercell(atoms, args.nx, args.ny, args.nz)
    out = structure.convert(sc, args.output)
    print(
        f"Supercelda {args.nx}x{args.ny}x{args.nz} ({len(sc)} átomos) escrita en: {out}"
    )
    return 0


def _cmd_convert(args) -> int:
    destino = args.output or args.output_flag
    if not destino:
        print("Error: falta el archivo de salida "
              "(olla-dft convert entrada.cif salida.vasp)", file=sys.stderr)
        return 1
    atoms = structure.load(args.file)
    out = structure.convert(atoms, destino)
    print(f"Estructura convertida: {args.file} -> {out}")
    return 0


# ----------------------------------------------------------------------
# Post-proceso
# ----------------------------------------------------------------------
def _cmd_gap(args) -> int:
    from qekit.modules import bands as bands_mod

    bs = bands_mod.load(args.path, prefix=args.prefix)
    print(bands_mod.gap_report(bs))
    if bs.result.converged is False:
        return 2
    return 0 if all(bands_mod.analyze_gap(bs, spin).is_metal or
                    bands_mod.analyze_gap(bs, spin).gap is not None
                    for spin in range(bs.result.nspin)) else 2


def _cmd_bands(args) -> int:
    from qekit.modules import bands as bands_mod

    bs = bands_mod.load(args.path, prefix=args.prefix)
    fat = fat_label = None
    if getattr(args, "fat", None):
        proy = bands_mod.leer_proyecciones(
            getattr(args, "projwfc", None) or args.path)
        bands_mod.comprobar_compatibilidad(bs, proy)
        fat = bands_mod.peso_de(proy, args.fat)
        fat_label = args.fat
        shift, _ = bands_mod.reference_energy(bs, args.ref)
        print(bands_mod.report_fat(proy, args.fat, bs, shift))
        print()
    print(bands_mod.gap_report(bs))
    print()
    written = bands_mod.export(bs, outdir=args.outdir, ref=args.ref)
    print("Datos exportados:")
    for f in written:
        print(f"  {f}")
    if not args.no_plot:
        imgs = bands_mod.plot(
            bs, outfile=str(Path(args.outdir) / "bandas"),
            ref=args.ref, emin=args.emin, emax=args.emax, dpi=args.dpi,
            formats=args.format, theme=args.template, size=args.size,
            family=args.font, background=args.background,
            palette=args.palette, usetex=args.usetex or None,
            width=args.width or "single", journal=args.journal,
            aspect=args.aspect or 0.88, mono=args.mono,
            title=args.title, gap_label=args.gap_label, panel=args.panel,
            fat=fat, fat_label=fat_label,
            fat_scale=getattr(args, "fat_scale", 55.0),
        )
        for f in imgs:
            print(f"  {f}")
    return 0


def _cmd_dos(args) -> int:
    from qekit.modules import dos as dos_mod

    dd = dos_mod.load(args.path, prefix=args.prefix)
    if getattr(args, "dband", None):
        pieza = str(args.dband).replace("_", "-").split("-")
        elemento = pieza[0].strip().capitalize()
        orbital = pieza[1].strip().lower() if len(pieza) > 1 else "d"
        m = dos_mod.momentos(dd, elemento, orbital,
                             emax=getattr(args, "dband_emax", None))
        print(dos_mod.report_momentos(m))
        return 0
    print(dos_mod.report(dd, ref=args.ref))
    print()
    written = dos_mod.export(dd, outdir=args.outdir, ref=args.ref)
    print("Datos exportados:")
    for f in written:
        print(f"  {f}")
    if not args.no_plot:
        imgs = dos_mod.plot(
            dd, outfile=str(Path(args.outdir) / "dos"),
            ref=args.ref, emin=args.emin, emax=args.emax, dpi=args.dpi,
            mode=args.mode, formats=args.format, theme=args.template,
            size=args.size, family=args.font, background=args.background,
            palette=args.palette, usetex=args.usetex or None,
            width=args.width or "single",
            journal=args.journal, aspect=args.aspect or 0.75,
            mono=args.mono, dash_mode=args.dashes,
            title=args.title, panel=args.panel,
        )
        for f in imgs:
            print(f"  {f}")
    return 0


def _cmd_plot(args) -> int:
    """Gráfica combinada bandas + DOS."""
    from qekit.modules import bands as bands_mod
    from qekit.modules import combined as combined_mod
    from qekit.modules import dos as dos_mod

    bs = bands_mod.load(args.path, prefix=args.prefix)
    dd = dos_mod.load(args.path, prefix=args.prefix)
    print(bands_mod.gap_report(bs))
    print()
    imgs = combined_mod.plot(
        bs, dd, outfile=str(Path(args.outdir) / "bandas_dos"),
        ref=args.ref, emin=args.emin, emax=args.emax, dpi=args.dpi,
        dos_mode=args.mode, formats=args.format, theme=args.template,
        size=args.size, family=args.font, background=args.background,
        palette=args.palette, usetex=args.usetex or None,
        width=args.width or "double", journal=args.journal,
        aspect=args.aspect or 0.46, mono=args.mono, dash_mode=args.dashes,
        title=args.title, gap_label=args.gap_label,
    )
    print("Figura combinada:")
    for f in imgs:
        print(f"  {f}")
    return 0


def _cmd_templates(args) -> int:
    action = args.action or "list"
    if action == "list":
        print("Plantillas disponibles:\n")
        user = qthemes.user_templates()
        for name in qthemes.names():
            origen = "usuario" if name in user else "incluida"
            desc = qthemes.load(name).get("description", "")
            print(f"  {name:12s} [{origen}]  {desc}")
        print(f"\nLas plantillas propias se leen de {qthemes.USER_DIR}")
        print("Para partir de una y modificarla:  olla-dft templates export <nombre>")
    elif action == "show":
        if not args.name:
            print("uso: olla-dft templates show <nombre>", file=sys.stderr)
            return 1
        print(qthemes.describe(args.name))
    elif action == "export":
        if not args.name:
            print("uso: olla-dft templates export <nombre> [archivo.json]",
                  file=sys.stderr)
            return 1
        out = qthemes.export(args.name, args.output)
        print(f"Plantilla escrita en: {out}")
        print("Edítala y úsala con:  olla-dft plot . --template " +
              Path(out).stem)
    return 0



# ----------------------------------------------------------------------
# Módulos de cálculo (barridos)
# ----------------------------------------------------------------------
def _figure_kwargs(args) -> dict:
    """Opciones de figura comunes, tomadas de los flags de la CLI.

    Con getattr y no con acceso directo: no todos los subcomandos declaran
    el juego completo de banderas de figura (los que usan _fig_opts_min no
    tienen --size), y un AttributeError aquí rompería el comando entero
    DESPUÉS de haber hecho todo el trabajo.
    """
    g = lambda n, d=None: getattr(args, n, d)                # noqa: E731
    return dict(
        formats=g("format", "pdf,png"), theme=g("template"), size=g("size"),
        family=g("font"), background=g("background"), palette=g("palette"),
        usetex=g("usetex") or None, journal=g("journal", "generic"),
        width=g("width") or "single", mono=g("mono", False), dpi=g("dpi"),
    )


def _run_or_explain(jobs, args, what: str):
    """Ejecuta el barrido si se pidió --run; si no, explica cómo correrlo."""
    from qekit.core import runner as run_mod

    par = getattr(args, "jobs", 1) or 1
    if getattr(args, "estimate", False) or (args.run and len(jobs) > 2):
        try:
            from qekit.modules import cost
            db = getattr(args, "db", None) or "olla-dft.db"
            modelo = cost.calibrar(db)
            est = cost.estimar_barrido(jobs, modelo, paralelo=par, db_path=db)
            print()
            print(cost.report(est, modelo))
        except Exception as exc:                            # noqa: BLE001
            if getattr(args, "estimate", False):
                print(f"\n(no se pudo estimar el coste: {exc})")
    if getattr(args, "estimate", False):
        print("\nNo se ha corrido nada. Quita --estimate para lanzarlo.")
        return None

    if getattr(args, "collect", False) and not args.run:
        # con --collect el usuario ya corrió los cálculos: explicarle cómo
        # lanzarlos, justo antes de analizarlos, solo confunde
        return None

    if not args.run:
        hilos = run_mod.nucleos()
        sugerido = max(2, hilos // 2) if hilos >= 4 and len(jobs) > 2 else None
        print()
        print("Los inputs están listos pero no se han corrido. Para ejecutarlos:")
        print(f"  olla-dft {what} ... --run          (Olla-DFT los lanza)")
        if sugerido:
            print(f"  olla-dft {what} ... --run -j {sugerido}     "
                  f"({sugerido} a la vez; son {len(jobs)} cálculos "
                  f"independientes)")
        print(f"  cd {args.outdir} && ./run.sh    (los lanzas tú"
              + (f"; ./run.sh {sugerido} para {sugerido} a la vez)"
                 if sugerido else ")"))
        print("Cuando terminen, vuelve a ejecutar el mismo comando con "
              "--collect para analizarlos.")
        return None
    print()
    return run_mod.run_all(jobs, pw_cmd=args.pw_cmd, nproc=args.nproc,
                           timeout=args.timeout,
                           paralelo=getattr(args, "jobs", 1) or 1,
                           rehacer=getattr(args, "redo", False),
                           presupuesto=_duracion(getattr(args, "max_time", None)))


def _cmd_converge(args) -> int:
    from qekit.modules import converge as conv_mod

    atoms = structure.load(args.file)
    values = None
    if args.values:
        values = [float(v) if args.kind != "kmesh" else v
                  for v in args.values.replace(";", ",").split(",") if v.strip()]
    run, rep = conv_mod.prepare(
        atoms, args.kind, outdir=args.outdir, values=values,
        threshold=args.threshold, pseudo_dir=args.pseudo_dir,
        insulator=args.insulator, ecutwfc=args.ecutwfc, ecutrho=args.ecutrho,
        kspacing=args.kspacing,
    )
    _print_prepare(rep)
    results = _run_or_explain(run.jobs, args, "converge")
    if results is None and not args.collect:
        return 0
    conv_mod.collect(run, results)
    print()
    print(conv_mod.report(run))
    print()
    written = conv_mod.export(run, args.outdir)
    for f in written:
        print(f"  {f}")
    if not args.no_plot:
        try:
            for f in conv_mod.plot(run, str(Path(args.outdir) / "convergencia"),
                                   aspect=args.aspect or 0.75, **_figure_kwargs(args)):
                print(f"  {f}")
        except Exception as exc:
            print(f"  (no se pudo graficar: {exc})")
    return 0


def _cmd_eos(args) -> int:
    from qekit.modules import eos as eos_mod

    atoms = structure.load(args.file)
    cell_a = atoms.cell.cellpar()[0]
    run, rep = eos_mod.prepare(
        atoms, outdir=args.outdir, npoints=args.npoints, span=args.span,
        pseudo_dir=args.pseudo_dir, insulator=args.insulator,
        ecutwfc=args.ecutwfc, ecutrho=args.ecutrho, kspacing=args.kspacing,
        relax_ions=args.relax_ions, center=args.scale,
    )
    _print_prepare(rep)
    results = _run_or_explain(run.jobs, args, "eos")
    if results is None and not args.collect:
        return 0
    eos_mod.collect(run, results)
    eos_mod.fit_all(run)
    print()
    print(eos_mod.report(run, cell_a))
    print()
    for f in eos_mod.export(run, args.outdir, cell_a):
        print(f"  {f}")
    if not args.no_plot:
        try:
            for f in eos_mod.plot(run, str(Path(args.outdir) / "eos"),
                                  equation=args.equation,
                                  aspect=args.aspect or 0.80,
                                  **_figure_kwargs(args)):
                print(f"  {f}")
        except Exception as exc:
            print(f"  (no se pudo graficar: {exc})")
    return 0


def _cmd_elastic(args) -> int:
    from qekit.modules import elastic as el_mod

    atoms = structure.load(args.file)
    run, rep = el_mod.prepare(
        atoms, outdir=args.outdir, delta=args.delta, npoints=args.npoints,
        pseudo_dir=args.pseudo_dir, insulator=args.insulator,
        ecutwfc=args.ecutwfc, ecutrho=args.ecutrho, kspacing=args.kspacing,
        ion_mode=args.ion_mode,
        dosd=args.dosd, espesor=args.thickness,
    )
    _print_prepare(rep)
    results = _run_or_explain(run.jobs, args, "elastic")
    if results is None and not args.collect:
        return 0
    el_mod.collect(run, results)
    el_mod.fit(run)
    print()
    print(el_mod.report(run))
    print()
    for f in el_mod.export(run, args.outdir):
        print(f"  {f}")
    if not args.no_plot:
        try:
            for f in el_mod.plot(run, str(Path(args.outdir) / "elastic"),
                                 aspect=args.aspect or 0.85,
                                 **_figure_kwargs(args)):
                print(f"  {f}")
        except Exception as exc:
            print(f"  (no se pudo graficar: {exc})")
    return 0



def _cmd_strain(args) -> int:
    from qekit.modules import strain as st_mod

    atoms = structure.load(args.file)
    nspin = getattr(args, "nspin", 1)
    mag = {}
    if getattr(args, "mag", None):
        mag = inputgen.parse_magnetization(args.mag, atoms.get_chemical_symbols())
        nspin = 2
    run, rep = st_mod.prepare(
        atoms, modo=args.mode, rangos=args.range, outdir=args.outdir,
        pseudo_dir=args.pseudo_dir, insulator=args.insulator,
        ecutwfc=args.ecutwfc, ecutrho=args.ecutrho, kspacing=args.kspacing,
        relax_ions=not args.fixed_ions, relax_perp=args.relax_perp,
        nspin=nspin, magnetization=mag,
        hubbard=_parse_hubbard(getattr(args, "hubbard", None)), vdw=args.vdw,
    )
    _print_prepare(rep)
    results = _run_or_explain(run.jobs, args, "strain")
    if results is None and not args.collect:
        return 0
    st_mod.collect(run, results)
    print()
    print(st_mod.report(run))
    print()
    for f in st_mod.export(run, args.outdir):
        print(f"  {f}")
    if not args.no_plot:
        try:
            for f in st_mod.plot(run, str(Path(args.outdir) / "strain"),
                                 aspect=args.aspect or 0.62,
                                 **_figure_kwargs(args)):
                print(f"  {f}")
        except Exception as exc:                            # noqa: BLE001
            print(f"  (no se pudo graficar: {exc})")
    return 0


def _cmd_adsorb(args) -> int:
    from qekit.modules import adsorb as ad_mod

    slab = structure.load(args.file)
    tipos = tuple(t.strip().lower() for t in args.sites.split(",") if t.strip())
    nspin = getattr(args, "nspin", 1)
    mag = {}
    if getattr(args, "mag", None):
        mag = inputgen.parse_magnetization(args.mag, slab.get_chemical_symbols())
        nspin = 2
    run, rep = ad_mod.prepare(
        slab, args.mol, outdir=args.outdir, altura=args.height, tipos=tipos,
        cara=args.face, rotaciones=args.rotations, ancla=args.anchor,
        pseudo_dir=args.pseudo_dir, insulator=args.insulator,
        ecutwfc=args.ecutwfc, ecutrho=args.ecutrho, kspacing=args.kspacing,
        relax_ions=not args.fixed_ions, vdw=args.vdw, dipolo=args.dipole,
        nspin=nspin, magnetization=mag,
    )
    _print_prepare(rep)
    results = _run_or_explain(run.jobs, args, "adsorb")
    if results is None and not args.collect:
        return 0
    ad_mod.collect(run, results)
    print()
    print(ad_mod.report(run))
    print()
    for f in ad_mod.export(run, args.outdir):
        print(f"  {f}")
    if not args.no_plot:
        try:
            for f in ad_mod.plot(run, str(Path(args.outdir) / "adsorcion"),
                                 aspect=args.aspect or 0.70,
                                 **_figure_kwargs(args)):
                print(f"  {f}")
        except Exception as exc:                            # noqa: BLE001
            print(f"  (no se pudo graficar: {exc})")
    return 0


def _cmd_eform(args) -> int:
    from qekit.modules import defects as df_mod

    atoms = structure.load(args.file)
    sc = _malla(args.supercell, "--supercell") or (2, 2, 2)
    try:
        cargas = [int(x) for x in args.charges.replace(";", ",").split(",")
                  if x.strip()]
    except ValueError:
        raise ErrorDeUso(
            "--charges son enteros separados por coma, por ejemplo "
            f"-2,-1,0,1,2; recibí '{args.charges}'.") from None
    if not cargas:
        raise ErrorDeUso("--charges necesita al menos un estado de carga.")
    pos = [float(x) for x in args.position.split(",")] if args.position else None
    nspin = getattr(args, "nspin", 1)
    mag = {}
    if getattr(args, "mag", None):
        mag = inputgen.parse_magnetization(args.mag, atoms.get_chemical_symbols())
        nspin = 2

    mu = {}
    for item in (args.mu or []):
        for trozo in str(item).replace(";", ",").split(","):
            if not trozo.strip():
                continue
            if "=" not in trozo:
                raise ErrorDeUso(
                    f"--mu se escribe ELEMENTO=eV, por ejemplo Si=-107.5; "
                    f"recibí '{trozo}'.")
            sym, _, val = trozo.partition("=")
            try:
                mu[sym.strip().capitalize()] = float(val)
            except ValueError:
                raise ErrorDeUso(
                    f"el μ de {sym.strip()} tiene que ser un número en eV; "
                    f"recibí '{val.strip()}'.") from None

    run, rep = df_mod.prepare(
        atoms, kind=args.kind, site=args.site, new_element=args.new_element,
        supercell=sc, position=pos, cargas=cargas, outdir=args.outdir,
        pseudo_dir=args.pseudo_dir, insulator=args.insulator,
        ecutwfc=args.ecutwfc, ecutrho=args.ecutrho, kspacing=args.kspacing,
        relax_ions=not args.fixed_ions, epsilon=args.epsilon,
        esquema=args.correction, vdw=args.vdw, nspin=nspin, magnetization=mag,
    )
    _print_prepare(rep)
    results = _run_or_explain(run.jobs, args, "eform")
    if results is None and not args.collect:
        return 0

    df_mod.collect(run, results, mu=mu)
    if run.aviso_mu:
        df_mod.asignar_mu_elemental(run, atoms.get_chemical_symbols())

    if args.align:
        al = df_mod.alineamiento(args.align[0], args.align[1])
        run.dV, run.dV_sigma = al["dV"], al["sigma"]
    if args.dv is not None:
        run.dV, run.dV_sigma = float(args.dv), None

    print()
    print(df_mod.report(run))
    print()
    for f in df_mod.export(run, args.outdir):
        print(f"  {f}")
    if not args.no_plot:
        try:
            for f in df_mod.plot(run, str(Path(args.outdir) / "formacion"),
                                 aspect=args.aspect or 0.78,
                                 **_figure_kwargs(args)):
                print(f"  {f}")
        except Exception as exc:                            # noqa: BLE001
            print(f"  (no se pudo graficar: {exc})")
    return 0


def _cmd_gamma(args) -> int:
    from qekit.modules import surfen

    atoms = structure.load(args.file)
    try:
        miller = [int(x) for x in args.miller.replace(",", " ").split()]
    except ValueError:
        raise ErrorDeUso(
            f"--miller son tres enteros, por ejemplo '1 1 1'; recibí "
            f"'{args.miller}'.") from None
    if len(miller) != 3:
        raise ErrorDeUso(
            f"--miller necesita TRES índices; recibí {len(miller)} "
            f"en '{args.miller}'.")
    try:
        capas = [int(x) for x in args.layers.replace(";", ",").split(",")
                 if x.strip()]
    except ValueError:
        raise ErrorDeUso(
            f"--layers son enteros separados por coma, por ejemplo 3,4,5,6; "
            f"recibí '{args.layers}'.") from None
    nspin = getattr(args, "nspin", 1)
    mag = {}
    if getattr(args, "mag", None):
        mag = inputgen.parse_magnetization(args.mag, atoms.get_chemical_symbols())
        nspin = 2

    run, rep = surfen.prepare(
        atoms, miller=miller, capas=capas, vacuum=args.vacuum,
        outdir=args.outdir, fijar=args.fix, relajar=args.relax,
        con_bulto=not args.no_bulk, reducir=not args.no_reduce,
        pseudo_dir=args.pseudo_dir,
        insulator=args.insulator, ecutwfc=args.ecutwfc, ecutrho=args.ecutrho,
        kspacing=args.kspacing, vdw=args.vdw, dipolo=args.dipole,
        nspin=nspin, magnetization=mag,
    )
    _print_prepare(rep)
    results = _run_or_explain(run.jobs, args, "gamma")
    if results is None and not args.collect:
        return 0
    surfen.collect(run, results)
    print()
    print(surfen.report(run))
    print()
    for f in surfen.export(run, args.outdir):
        print(f"  {f}")
    if not args.no_plot:
        try:
            for f in surfen.plot(run, str(Path(args.outdir) / "gamma"),
                                 aspect=args.aspect or 0.72,
                                 **_figure_kwargs(args)):
                print(f"  {f}")
        except Exception as exc:                            # noqa: BLE001
            print(f"  (no se pudo graficar: {exc})")
    return 0


# ----------------------------------------------------------------------
# Materiales laminares
# ----------------------------------------------------------------------
def _cmd_layers(args) -> int:
    from qekit.core import layers as layers_mod
    from qekit.modules import xrd as xrd_mod

    atoms = structure.load(args.file)
    res = layers_mod.analyze(atoms, tol=args.tol)
    lam = xrd_mod.wavelength_value(args.wavelength)
    print(layers_mod.report(atoms, res, wavelength=lam,
                            radiation=xrd_mod.wavelength_name(args.wavelength)))
    if args.slab and res.layers:
        slab = layers_mod.make_slab(atoms, res, vacuum=args.vacuum)
        out = structure.convert(slab, args.slab)
        print(f"\nMonocapa con {args.vacuum:g} Å de vacío escrita en: {out}")
    return 0


def _cmd_xrd(args) -> int:
    from qekit.modules import xrd as xrd_mod

    atoms = structure.load(args.file)
    rng = (args.tt_min, args.tt_max)
    pattern = xrd_mod.compute(atoms, wavelength=args.wavelength,
                              two_theta_range=rng, b_iso=args.biso,
                              basis=args.basis)
    xrd_mod.broaden(pattern, two_theta_range=rng, fwhm=args.fwhm,
                    size_nm=args.size)
    print(xrd_mod.report(pattern))
    exp = None
    if args.exp:
        exp = xrd_mod.read_experimental(args.exp)
        print(f"\nComparando con el difractograma experimental '{args.exp}'.")
    Path(args.outdir).mkdir(parents=True, exist_ok=True)
    for f in xrd_mod.export(pattern, args.outdir):
        print(f"  {f}")
    if args.suite:
        from qekit.modules import interop
        print("  " + interop.write(interop.from_xrd(pattern, atoms),
                                   Path(args.outdir) / "DRX_suite.json"))
    if not args.no_plot:
        for f in xrd_mod.plot(
            pattern, str(Path(args.outdir) / "xrd"), exp=exp,
            formats=args.format, theme=args.template, size=args.size_preset,
            family=args.font, background=args.background,
            palette=args.palette, usetex=args.usetex or None,
            width=args.width or "single", journal=args.journal,
            aspect=args.aspect or 0.62, mono=args.mono, dpi=args.dpi,
        ):
            print(f"  {f}")
    return 0


def _cmd_exfoliate(args) -> int:
    from qekit.modules import exfoliate as exf_mod

    atoms = structure.load(args.file)
    run, rep = exf_mod.prepare(
        atoms, outdir=args.outdir, vacuum=args.vacuum, vdw=args.vdw,
        tol=args.tol, pseudo_dir=args.pseudo_dir, ecutwfc=args.ecutwfc,
        ecutrho=args.ecutrho, kspacing=args.kspacing,
        insulator=args.insulator, relax_slab=args.relax_slab,
    )
    _print_prepare(rep)
    results = _run_or_explain(run.jobs, args, "exfoliate")
    if results is None and not args.collect:
        return 0
    exf_mod.collect(run, results)
    print()
    print(exf_mod.report_result(run))
    return 0



# ----------------------------------------------------------------------
# Campos, ópticas y fonones
# ----------------------------------------------------------------------
_AXES = {"a": 0, "b": 1, "c": 2, "x": 0, "y": 1, "z": 2}


def _cmd_wf(args) -> int:
    from qekit.core import qeout
    from qekit.modules import fields as f_mod

    axis = _AXES.get(args.axis.lower(), 2)
    cube_path = Path(args.path) / "potencial.cube"
    if not cube_path.exists() or args.rerun:
        print("Ejecutando pp.x para extraer el potencial electrostático...")
        cube_path = f_mod.run_pp(args.path, "potential", "potencial",
                                 pw_cmd=args.pw_cmd, nproc=args.nproc)
    cube = f_mod.read_cube(str(cube_path))
    qe = qeout.read_xml(args.path)
    if qe.fermi is None:
        print("Error: el XML no trae energía de Fermi (¿terminó el scf?)",
              file=sys.stderr)
        return 1
    # con las posiciones la meseta se busca en el hueco real sin átomos
    wf = f_mod.work_function(cube, qe.fermi, axis=axis,
                             positions=qe.positions)
    print(f_mod.report_wf(wf))
    Path(args.outdir).mkdir(parents=True, exist_ok=True)
    print()
    for f in f_mod.export_wf(wf, args.outdir):
        print(f"  {f}")
    if not args.no_plot:
        for f in f_mod.plot_profile(
            wf, str(Path(args.outdir) / "funcion_trabajo"),
            formats=args.format, theme=args.template, family=args.font,
            background=args.background, palette=args.palette,
            usetex=args.usetex or None, width=args.width or "single",
            journal=args.journal, mono=args.mono, dpi=args.dpi,
        ):
            print(f"  {f}")
    return 0


def _cmd_align(args) -> int:
    from qekit.modules import align as al_mod

    eje = _AXES.get(args.axis.lower(), 2)
    modo = "interfaz" if args.interface else "vacio"
    nombres = [n.strip() for n in (args.names or "").split(",") if n.strip()]
    na = nombres[0] if len(nombres) > 0 else None
    nb = nombres[1] if len(nombres) > 1 else None

    a = al_mod.leer_lado(args.a, na, modo=modo, eje=eje, rerun=args.rerun,
                         pw_cmd=args.pw_cmd, nproc=args.nproc)
    b = al_mod.leer_lado(args.b, nb, modo=modo, eje=eje, rerun=args.rerun,
                         pw_cmd=args.pw_cmd, nproc=args.nproc)
    puente = None
    if args.interface:
        pu = al_mod.puente_interfaz(args.interface, eje=eje, ancho=args.window,
                                    rerun=args.rerun, pw_cmd=args.pw_cmd,
                                    nproc=args.nproc)
        puente = pu["delta"]
    al = al_mod.alinear(a, b, modo=modo, puente=puente)
    print(al_mod.report(al))
    print()
    for f in al_mod.export(al, args.outdir):
        print(f"  {f}")
    if not args.no_plot:
        try:
            for f in al_mod.plot(al, str(Path(args.outdir) / "alineamiento"),
                                 **_figure_kwargs(args)):
                print(f"  {f}")
        except Exception as exc:                            # noqa: BLE001
            print(f"  (no se pudo graficar: {exc})")
    return 0


def _cmd_charge(args) -> int:
    from qekit.modules import fields as f_mod

    axis = _AXES.get(args.axis.lower(), 2)
    name = {"density": "densidad", "elf": "elf", "spin": "espin",
            "potential": "potencial", "vtotal": "vtotal"}[args.field]
    cube_path = Path(args.path) / f"{name}.cube"
    if not cube_path.exists() or args.rerun:
        print(f"Ejecutando pp.x ({args.field})...")
        cube_path = f_mod.run_pp(args.path, args.field, name,
                                 pw_cmd=args.pw_cmd, nproc=args.nproc)
    cube = f_mod.read_cube(str(cube_path))
    z, prof = f_mod.planar_average(cube, axis)
    desc = f_mod.PLOTS[args.field][1]
    print(f"Campo: {desc}")
    print(f"Malla: {cube.shape[0]}x{cube.shape[1]}x{cube.shape[2]}  |  "
          f"cube: {cube_path}")
    print("El .cube se abre directo en VESTA para las isosuperficies.")
    Path(args.outdir).mkdir(parents=True, exist_ok=True)
    import numpy as np
    np.savetxt(Path(args.outdir) / "PERFIL_PLANAR.dat",
               np.column_stack([z, prof]), fmt="%14.6e",
               header=f"# perfil planar de {args.field} a lo largo de "
                      f"{args.axis}\n# z(A)  campo", comments="")
    print(f"  {Path(args.outdir) / 'PERFIL_PLANAR.dat'}")
    if not args.no_plot:
        for f in f_mod.plot_density_profile(
            z, prof, str(Path(args.outdir) / f"perfil_{name}"), label=desc,
            formats=args.format, theme=args.template, family=args.font,
            background=args.background, palette=args.palette,
            usetex=args.usetex or None, width=args.width or "single",
            journal=args.journal, mono=args.mono, dpi=args.dpi,
        ):
            print(f"  {f}")
    return 0


def _cmd_optics(args) -> int:
    from qekit.modules import optics as opt_mod

    atoms = structure.load(args.file)
    run, rep = opt_mod.prepare(
        atoms, outdir=args.outdir, pseudo_dir=args.pseudo_dir,
        ecutwfc=args.ecutwfc, ecutrho=args.ecutrho,
        kspacing=args.kspacing or 0.12, insulator=not args.metal,
        wmax=args.wmax, intersmear=args.smear,
    )
    _print_prepare(rep)
    if args.run:
        from qekit.core import runner as run_mod
        print()
        results = run_mod.run_all(run.jobs, pw_cmd=args.pw_cmd,
                                  nproc=args.nproc, timeout=args.timeout)
        if not all(r.ok for r in results):
            return 1
        print("  epsilon.x ... ", end="", flush=True)
        opt_mod.run_epsilon(run, pw_cmd=args.pw_cmd, nproc=args.nproc)
        print("hecho")
    elif not args.collect:
        print("\nCorre con --run, o ejecuta scf, nscf y epsilon.x a mano y "
              "vuelve con --collect.")
        return 0
    run = opt_mod.collect(run)
    if args.scissor:
        run = opt_mod.scissor(run, args.scissor)
    print()
    print(opt_mod.report(run))
    print()
    for f in opt_mod.export(run, args.outdir):
        print(f"  {f}")
    if args.suite:
        from qekit.modules import interop
        print("  " + interop.write(
            interop.from_optics(run, atoms, args.tauc),
            Path(args.outdir) / "OPTICA_suite.json"))
    if not args.no_plot:
        for f in opt_mod.plot(
            run, str(Path(args.outdir) / "opticas"), formats=args.format,
            theme=args.template, family=args.font, background=args.background,
            palette=args.palette, usetex=args.usetex or None,
            width=args.width or "double", journal=args.journal,
            mono=args.mono, dpi=args.dpi, tauc_kind=args.tauc,
        ):
            print(f"  {f}")
    return 0


def _cmd_effmass(args) -> int:
    import glob as _glob

    from qekit.modules import bands as bands_mod
    from qekit.modules import effmass as em

    atoms = structure.load(args.file)

    if args.collect:
        meta = em.load_meta(args.outdir)
        xmls = _glob.glob(str(Path(args.outdir) / "out" / "*.xml"))
        if not xmls:
            print(f"Error: no hay XML en {args.outdir}/out; "
                  "¿corriste scf.in y masa.in?", file=sys.stderr)
            return 1
        run = em.collect_fine(xmls[0], meta)
        print(em.report(run))
        print()
        for f in em.export(run, args.outdir):
            print(f"  {f}")
        return 0

    if not args.bands_dir:
        print("Error: hace falta --bands-dir con un cálculo de bandas ya "
              "hecho\n(de ahí se localizan el VBM y el CBM).", file=sys.stderr)
        return 1

    bs = bands_mod.load(args.bands_dir)
    ventana = args.window if args.window is not None else em.WINDOW_DEFAULT
    quick = em.from_bands(bs, window=ventana, min_pts=args.min_points)
    print(em.report(quick))
    if quick.is_metal:
        return 0

    print()
    meta, rep = em.prepare(
        atoms, bs, outdir=args.outdir, pseudo_dir=args.pseudo_dir,
        ecutwfc=args.ecutwfc, ecutrho=args.ecutrho,
        half_width=args.half_width, npts=args.points,
    )
    _print_prepare(rep)

    if args.run:
        from qekit.core import runner as run_mod
        print()
        jobs = [run_mod.Job(name="scf", directory=Path(args.outdir),
                            input_file="scf.in", output_file="scf.out"),
                run_mod.Job(name="masa", directory=Path(args.outdir),
                            input_file="masa.in", output_file="masa.out")]
        results = run_mod.run_all(jobs, pw_cmd=args.pw_cmd,
                                  nproc=args.nproc, timeout=args.timeout)
        if not all(r.ok for r in results):
            return 1
        xmls = _glob.glob(str(Path(args.outdir) / "out" / "*.xml"))
        run = em.collect_fine(xmls[0], meta)
        print()
        print(em.report(run))
        print()
        for f in em.export(run, args.outdir):
            print(f"  {f}")
    return 0


def _cmd_surface(args) -> int:
    from qekit.modules import builder

    atoms = structure.load(args.file)
    miller = tuple(int(x) for x in args.miller.replace(",", " ").split())
    info = builder.surface(atoms, miller=miller, layers=args.layers,
                           vacuum=args.vacuum, fix_layers=args.fix)
    print(builder.report_slab(info))
    if args.output:
        structure.convert(info.atoms, args.output)
        print(f"\n  {args.output}")
        if info.fijados and not structure.conserva_fijos(args.output):
            print(
                f"\nAVISO: '{args.output}' no guarda qué átomos están "
                f"congelados (--fix {args.fix}): el CIF no tiene dónde "
                "ponerlo,\n  así que al volver a cargarlo se relajaría todo. "
                f"Escribe la losa en {builder.FORMATO_CON_FIJOS}, por "
                "ejemplo\n  '-o losa.vasp': ese formato conserva los fijos y "
                "'olla-dft inputgen' los traduce\n  a '0 0 0' en "
                "ATOMIC_POSITIONS. O usa directamente 'olla-dft gamma "
                "--fix'.", file=sys.stderr)
    return 0


def _cmd_defect(args) -> int:
    from qekit.modules import builder

    atoms = structure.load(args.file)
    sc = _malla(args.supercell, "--supercell") or (2, 2, 2)
    pos = None
    if args.position:
        pos = [float(x) for x in args.position.split(",")]
    perfecto, info = builder.defect(
        atoms, kind=args.kind, site=args.site,
        new_element=args.new_element, supercell=sc, position=pos)
    print(builder.report_defect(info))
    out = Path(args.outdir); out.mkdir(parents=True, exist_ok=True)
    f1 = out / "perfecto.cif"
    f2 = out / "defecto.cif"
    structure.convert(perfecto, str(f1))
    structure.convert(info.atoms, str(f2))
    print(f"\n  {f1}\n  {f2}")
    return 0


def _cmd_charges(args) -> int:


    from qekit.modules import charges as ch
    from qekit.modules import fields

    atoms = structure.load(args.file) if args.file else None
    hecho = False

    # Z_valencia por átomo desde los UPF: sin ella la columna 'neta' no
    # se puede calcular y queda en n/d.
    valence = None
    if atoms is not None and (args.lowdin or args.bader):
        pdir = getattr(args, "pseudo_dir", None) or qcfg.load()["pseudo_dir"]
        valence = ch.valence_from_pseudos(atoms.get_chemical_symbols(), pdir)
        if valence is None:
            print(f"Aviso: no pude leer z_valence de los UPF en '{pdir}'; la "
                  "columna 'neta' saldrá como n/d.\n  Pasa la carpeta con "
                  "--pseudo-dir para tener la carga neta por átomo.",
                  file=sys.stderr)

    if args.lowdin:
        res = ch.read_lowdin(
            args.lowdin,
            symbols=atoms.get_chemical_symbols() if atoms else None,
            valence=valence)
        print(ch.report_lowdin(res))
        hecho = True

    if args.bader:
        cube = fields.read_cube(args.bader)
        if atoms is None:
            print("Error: --bader necesita también la estructura",
                  file=sys.stderr)
            return 1
        res = ch.bader(cube, atoms.positions,
                       symbols=atoms.get_chemical_symbols(),
                       valence=valence)
        if hecho:
            print()
        print(ch.report_bader(res))
        hecho = True

    if args.difference:
        partes = [fields.read_cube(f) for f in args.difference[1:]]
        total = fields.read_cube(args.difference[0])
        d = ch.difference(total, partes)
        if hecho:
            print()
        print(ch.report_difference(d, axis=args.axis))
        if not args.no_plot:
            for f in ch.plot_difference(
                d, str(Path(args.outdir) / "diferencia_carga"),
                axis=args.axis, formats=args.format, theme=args.template,
                family=args.font, background=args.background,
                palette=args.palette, usetex=args.usetex or None,
                width=args.width or "single", journal=args.journal,
                mono=args.mono, dpi=args.dpi,
            ):
                print(f"  {f}")
        hecho = True

    if not hecho:
        print("Nada que hacer: usa --lowdin, --bader o --difference.",
              file=sys.stderr)
        return 1
    return 0


def _cmd_fermi(args) -> int:
    import glob as _glob

    from qekit.modules import transport as tr

    xmls = _glob.glob(str(Path(args.outdir) / "out" / "*.xml"))
    if not xmls:
        print(f"Error: no hay XML en {args.outdir}/out; corre primero "
              "'olla-dft transport ... --run'", file=sys.stderr)
        return 1
    run = tr.load(xmls[0])
    from qekit.core import qeout
    res = qeout.read_xml(xmls[0])
    cruzan = tr.crossing_bands(run)
    print("--- Superficie de Fermi ---")
    print(f"Malla: {run.grid[0]}x{run.grid[1]}x{run.grid[2]}  |  "
          f"E_F = {run.fermi:.4f} eV")
    if not cruzan:
        print("Ninguna banda cruza E_F: el sistema no es metálico y no "
              "tiene superficie de Fermi.")
        return 0
    print(f"Bandas que cruzan E_F: {[b + 1 for b in cruzan]}")
    destino = Path(args.outdir) / "superficie_fermi.bxsf"
    print("\n  " + tr.export_bxsf(run, res.cell, destino))
    print("\nÁbrelo con XCrySDen (xcrysden --bxsf archivo) o FermiSurfer.")
    return 0


def _cmd_xps(args) -> int:
    from qekit.modules import xps as xps_mod

    atoms = structure.load(args.file)
    if args.collect:
        f = Path(args.outdir) / "initial_state.out"
        if not f.exists():
            print(f"Error: falta {f}", file=sys.stderr)
            return 1
        res = xps_mod.collect(f, symbols=atoms.get_chemical_symbols())
        print(xps_mod.report(res))
        print()
        for x in xps_mod.export(res, args.outdir):
            print(f"  {x}")
        if args.suite:
            from qekit.modules import interop
            print("  " + interop.write(
                interop.from_xps(res, atoms),
                Path(args.outdir) / "XPS_suite.json"))
        return 0

    core_hole = {}
    for par in (args.core_hole or []):
        if "=" not in par:
            print(f"Error: --core-hole se escribe Elemento=archivo.UPF "
                  f"(recibí '{par}')", file=sys.stderr)
            return 1
        el, upf = par.split("=", 1)
        core_hole[el.strip()] = upf.strip()

    _c, rep = xps_mod.prepare(
        atoms, outdir=args.outdir, pseudo_dir=args.pseudo_dir,
        ecutwfc=args.ecutwfc, ecutrho=args.ecutrho,
        kspacing=args.kspacing, insulator=not args.metal,
        core_hole=core_hole or None)
    _print_prepare(rep)
    return 0


def _cmd_corehole(args) -> int:
    from qekit.modules import corehole as ch

    if args.core_wfc:
        destino = args.output or (Path(args.core_wfc).stem + ".wfc")
        f = ch.core_wfc(args.core_wfc, destino, orbital=args.orbital)
        print(f"Funcion de onda de core escrita en {f}")
        print("Se le pasa a xspectra.x en &pseudos como filecore='...'.")
        return 0

    g = ch.generar(args.element, borde=args.edge, outdir=args.outdir,
                   dft=args.functional, rcut=args.rcut, rel=args.rel,
                   semicore=args.semicore, pseudotype=args.pseudotype,
                   proyectores=args.projectors, solo_base=args.plain,
                   correr=not args.only_inputs, ld1_cmd=args.ld1_cmd)
    if args.only_inputs:
        print("Entradas de ld1.x escritas (no se ejecuto nada):")
        for e in g.entradas:
            print(f"  {e}")
        print("\nCorrelas con:  ld1.x < ld1_base.in > ld1_base.out")
        return 0
    print(ch.report(g))
    print()
    for x in ch.export(g, args.outdir):
        print(f"  {x}")
    return 0 if g.ok else 1


def _cmd_xanes(args) -> int:
    from qekit.modules import xanes as xa

    atoms = structure.load(args.file)
    if args.collect:
        run = xa.collect(args.outdir, elemento=args.element or "",
                         borde=args.edge)
        run.natoms = len(atoms)
        print(xa.report(run))
        print()
        salidas = xa.export(run, args.outdir)
        if not args.no_plot:
            salidas += xa.plot(
                run, str(Path(args.outdir) / "xanes"), formats=args.format,
                theme=args.template, family=args.font,
                background=args.background, palette=args.palette,
                usetex=args.usetex, width=args.width or "single",
                journal=args.journal, mono=args.mono, dpi=args.dpi)
        for x in salidas:
            print(f"  {x}")
        return 0

    if not args.element:
        raise ErrorDeUso(
            "falta --element: hay que decir QUE atomo absorbe.")
    borde = xa.validar_borde(args.edge)
    if not args.core_hole:
        raise ErrorDeUso(
            "falta --core-hole con el UPF de hueco de core. Sin el se "
            "calcularia el espectro del estado fundamental, que no es lo "
            "que mide el experimento.\n"
            f"  Genéralo con:  olla-dft corehole {args.element} "
            f"--edge {xa.BORDE_COREHOLE[borde]}")

    _c, rep = xa.prepare(
        atoms, args.element, args.core_hole, outdir=args.outdir,
        sitio=args.site, borde=args.edge, pseudo_dir=args.pseudo_dir,
        ecutwfc=args.ecutwfc, ecutrho=args.ecutrho, kspacing=args.kspacing,
        insulator=not args.metal, promedio=args.average,
        polarizacion=_vector3(args.polarization),
        xemin=args.emin, xemax=args.emax, xgamma=args.broadening,
        r_paw=args.r_paw)
    _print_prepare(rep)
    return 0


def _vector3(texto) -> tuple:
    if not texto:
        return (1.0, 0.0, 0.0)
    partes = [p for p in re.split(r"[,\s]+", str(texto).strip()) if p]
    if len(partes) != 3:
        raise ErrorDeUso(
            f"la polarizacion necesita TRES numeros, por ejemplo '0 0 1'; "
            f"recibi '{texto}'.")
    try:
        return tuple(float(p) for p in partes)
    except ValueError:
        raise ErrorDeUso(
            f"la polarizacion solo admite numeros; recibi '{texto}'.") from None


def _parse_mag(texto, simbolos=None):
    """Magnetizacion inicial desde la bandera --mag, o None."""
    if not texto:
        return None
    return inputgen.parse_magnetization(texto, simbolos or [])


def _buscar_hubbard_dat(carpeta):
    """El .Hubbard_parameters.dat que escribe hp.x, se llame como se llame."""
    base = Path(carpeta)
    for patron in ("*.Hubbard_parameters.dat", "**/*.Hubbard_parameters.dat"):
        encontrados = sorted(base.glob(patron))
        if encontrados:
            return encontrados[0]
    return None


def _cmd_hubbard(args) -> int:
    from qekit.modules import hubbard as hb

    atoms = structure.load(args.file)
    qgrid = _malla(args.qgrid, "--qgrid") or (2, 2, 2)
    if args.mag and args.nspin == 1:
        args.nspin = 2      # pedir magnetizacion implica activar el espin
    especies = [e.strip() for e in (args.species or "").split(",") if e.strip()]

    if args.collect:
        run = hb.collect(args.outdir, qgrid=qgrid, proyeccion=args.projection)
        print(hb.report(run))
        if getattr(args, "intersite", False):
            dat = _buscar_hubbard_dat(args.outdir)
            print()
            if dat is None:
                print("No encontré el <prefix>.Hubbard_parameters.dat de hp.x "
                      f"en '{args.outdir}'.")
            else:
                pares, sup = hb.leer_v(dat)
                run.v_pares, run.supercelda_v = pares, sup
                print(hb.report_v(pares, sup, umbral=args.v_threshold))
                if pares:
                    tarjeta = hb.tarjeta_hubbard(
                        run.sitios, pares, proyeccion=args.projection,
                        umbral_v=args.v_threshold)
                    destino = Path(args.outdir) / "HUBBARD.card"
                    destino.write_text(tarjeta)
                    print()
                    print("Tarjeta para el siguiente scf (QE >= 7.1):")
                    print("  " + str(destino))
                    print(tarjeta)
        print()
        for x in hb.export(run, args.outdir):
            print(f"  {x}")
        return 0

    if args.cycle:
        run = hb.ciclo(
            atoms, outdir=args.outdir, especies=especies or None,
            qgrid=qgrid, max_iter=args.max_iter, tol=args.tol,
            mezcla=args.mixing, pw_cmd=args.pw_cmd, nproc=args.nproc,
            pseudo_dir=args.pseudo_dir, ecutwfc=args.ecutwfc,
            ecutrho=args.ecutrho, kspacing=args.kspacing,
            insulator=not args.metal, proyeccion=args.projection,
            hubbard_style=getattr(args, "hubbard_style", "legacy"),
            nspin=args.nspin, magnetization=_parse_mag(args.mag,
                                 atoms.get_chemical_symbols()))
        print(hb.report(run))
        print()
        for x in hb.export(run, args.outdir):
            print(f"  {x}")
        return 0 if run.convergido else 1

    _c, rep = hb.prepare(
        atoms, outdir=args.outdir, especies=especies or None,
        qgrid=qgrid, pseudo_dir=args.pseudo_dir, ecutwfc=args.ecutwfc,
        ecutrho=args.ecutrho, kspacing=args.kspacing,
        insulator=not args.metal, proyeccion=args.projection,
        hubbard_style=getattr(args, "hubbard_style", "legacy"),
        nspin=args.nspin, magnetization=_parse_mag(args.mag,
                                 atoms.get_chemical_symbols()))
    _print_prepare(rep)
    return 0


def _cmd_interface(args) -> int:
    from qekit.modules import interface as itf

    a1 = structure.load(args.file1)
    a2 = structure.load(args.file2)
    if args.list:
        cands = itf.buscar(a1, a2, max_index=args.max_index, tol=args.tol,
                           max_atoms=args.max_atoms, n_mejores=args.top)
        if not cands:
            print("No hay ninguna coincidencia con esos limites.",
                  file=sys.stderr)
            return 1
        print(f"--- Superceldas comunes de {a1.get_chemical_formula()} y "
              f"{a2.get_chemical_formula()} ---")
        print(f"  {'#':>2s} {'atomos':>7s} {'deformacion':>12s} "
              f"{'n1':>4s} {'n2':>4s}  {'area(A2)':>9s}")
        for i, c in enumerate(cands):
            print(f"  {i:2d} {c.natoms:7d} {c.eps_pct:11.2f} % "
                  f"{c.n1:4d} {c.n2:4d}  {c.area:9.2f}")
        print("\nSe construye una con --index.")
        return 0

    het = itf.emparejar(
        a1, a2, max_index=args.max_index, tol=args.tol,
        max_atoms=args.max_atoms, indice=args.index,
        separacion=args.separation, vacio=args.vacuum,
        deformar=args.strain, desplazamiento=_par2(args.shift))
    print(itf.report(het))
    print()
    for x in itf.export(het, args.outdir, nombre=args.name):
        print(f"  {x}")
    return 0


def _par2(texto) -> tuple:
    if not texto:
        return (0.0, 0.0)
    partes = [p for p in re.split(r"[,\s]+", str(texto).strip()) if p]
    if len(partes) != 2:
        raise ErrorDeUso(
            f"--shift necesita DOS numeros (fracciones de la celda), por "
            f"ejemplo '0.33 0.33'; recibi '{texto}'.")
    try:
        return tuple(float(p) for p in partes)
    except ValueError:
        raise ErrorDeUso(
            f"--shift solo admite numeros; recibi '{texto}'.") from None


def _cmd_md(args) -> int:
    from qekit.modules import dynamics as dyn

    tray = dyn.leer_md(args.path, skip=args.skip)
    a = dyn.analizar(tray, rmax=args.rmax, nbins=args.bins,
                     equilibrado=args.skip)
    print(dyn.report(a))
    print()
    salidas = dyn.export(a, args.outdir)
    if not args.no_plot:
        salidas += dyn.plot(
            a, str(Path(args.outdir) / "md"), formats=args.format,
            theme=args.template, family=args.font,
            background=args.background, palette=args.palette,
            usetex=args.usetex, width=args.width or "double",
            journal=args.journal, mono=args.mono, dpi=args.dpi)
    for x in salidas:
        print(f"  {x}")
    return 0


def _cmd_neb(args) -> int:
    from qekit.modules import neb as nb

    if args.collect:
        run = nb.collect(args.outdir, prefix=args.prefix)
        print(nb.report(run))
        print()
        salidas = nb.export(run, args.outdir)
        if not args.no_plot:
            salidas += nb.plot(
                run, str(Path(args.outdir) / "neb"), formats=args.format,
                theme=args.template, family=args.font,
                background=args.background, palette=args.palette,
                usetex=args.usetex, width=args.width or "single",
                journal=args.journal, mono=args.mono, dpi=args.dpi)
        for x in salidas:
            print(f"  {x}")
        return 0 if run.convergido is not False else 1

    if not args.final:
        raise ErrorDeUso(
            "faltan las DOS estructuras: la inicial y la final.\n"
            "  olla-dft neb reactivo.cif producto.cif -o camino")
    ini = structure.load(args.file)
    fin = structure.load(args.final)
    fijos = [int(x) for x in re.split(r"[,\s]+", args.fix.strip())
             if x] if args.fix else None
    _c, rep = nb.prepare(
        ini, fin, outdir=args.outdir, n_imagenes=args.images,
        ci=not args.no_ci, pseudo_dir=args.pseudo_dir,
        ecutwfc=args.ecutwfc, ecutrho=args.ecutrho, kspacing=args.kspacing,
        insulator=not args.metal, path_thr=args.path_thr,
        nstep_path=args.nstep, nspin=args.nspin,
        magnetization=_parse_mag(args.mag, ini.get_chemical_symbols()),
        fijos=fijos)
    _print_prepare(rep)
    return 0


def _cmd_thermochem(args) -> int:
    from qekit.modules import thermochem as tc

    nu = _leer_frecuencias(args.freqs)
    # Una molécula en fase gas no tiene celda, y aquí no hace falta: solo
    # se usan masas y geometría. structure.load exige celda porque sirve
    # para cálculos periódicos, así que para esto se lee directo.
    atoms = None
    if args.structure:
        from ase.io import read as _read
        try:
            atoms = _read(args.structure)
        except Exception as exc:
            raise ErrorDeUso(
                f"no se pudo leer '{args.structure}': {exc}") from None
    tq = tc.corregir(nu, T=args.temp, fase=args.phase, atoms=atoms,
                     p=args.pressure * 1e5, simetria=args.symmetry,
                     multiplicidad=args.multiplicity, piso=args.floor)
    print(tc.report(tq, E_dft=args.energy))
    if args.outdir:
        out = Path(args.outdir); out.mkdir(parents=True, exist_ok=True)
        f = out / "TERMOQUIMICA.txt"
        f.write_text(tc.report(tq, E_dft=args.energy) + "\n")
        print(f"\n  {f}")
    return 0


def _leer_frecuencias(fuente):
    """Frecuencias desde un archivo (una por linea o en columnas) o inline."""
    import numpy as _np
    p = Path(fuente)
    if p.exists():
        datos = _np.loadtxt(p, comments="#")
        return _np.atleast_1d(datos if datos.ndim == 1 else datos[:, -1])
    partes = [x for x in re.split(r"[,\s]+", str(fuente).strip()) if x]
    try:
        return _np.array([float(x) for x in partes])
    except ValueError:
        raise ErrorDeUso(
            f"'{fuente}' no es un archivo ni una lista de numeros. Pasa el "
            "archivo de frecuencias (por ejemplo FONONES_GAMMA.dat) o los "
            "valores separados por coma.") from None


def _cmd_unfold(args) -> int:
    from qekit.modules import unfold as uf

    prim = structure.load(args.primitive)
    d = uf.desdoblar(args.path, prim.get_cell(), prefix=args.prefix,
                     bandas=range(args.bands) if args.bands else None,
                     spin=getattr(args, "spin", "up") or "up")
    print(uf.report(d))
    print()
    salidas = uf.export(d, args.outdir)
    if not args.no_plot:
        salidas += uf.plot(
            d, str(Path(args.outdir) / "unfold"), formats=args.format,
            emin=args.emin, emax=args.emax, theme=args.template,
            family=args.font, background=args.background,
            palette=args.palette, usetex=args.usetex,
            width=args.width or "single", journal=args.journal,
            mono=args.mono, dpi=args.dpi)
    for x in salidas:
        print(f"  {x}")
    return 0


def _cmd_elph(args) -> int:
    from qekit.modules import elph as ep

    if args.collect:
        run = ep.leer_elph_ph(args.outdir)
        lam_out = Path(args.outdir) / "lambda.out"
        if lam_out.exists():
            run2 = ep.leer_lambda_out(lam_out)
            for campo in ("lambdas", "dos_ef", "omega_log", "Tc", "sigmas"):
                v = getattr(run2, campo)
                if v is not None:
                    setattr(run, campo, v)
        a2f = sorted(Path(args.outdir).glob("a2F.dos*")) + \
            sorted(Path(args.outdir).glob("A2F.dat"))
        if a2f:
            run.freq, run.a2F = ep.leer_a2F(a2f[0])
            if run.lambdas is None:
                lam = ep.lambda_de_a2F(run.freq, run.a2F)
                run.lambdas = _np_array([lam])
                run.omega_log = _np_array(
                    [ep.omega_log_de_a2F(run.freq, run.a2F, lam)])
            # La media cuadrática solo sale del espectro; con ella la Tc
            # lleva el factor de forma f2 de Allen-Dynes, no solo f1.
            run.omega_2 = ep.omega_2(run.freq, run.a2F, run.lam)
        if run.lambdas is not None:
            run.i_plato = ep.plato(run.lambdas)
        print(ep.report(run, T_debye=args.debye))
        print()
        salidas = ep.export(run, args.outdir)
        if not args.no_plot and (run.a2F is not None or
                                 run.lambdas is not None):
            salidas += ep.plot(
                run, str(Path(args.outdir) / "elph"), formats=args.format,
                theme=args.template, family=args.font,
                background=args.background, palette=args.palette,
                usetex=args.usetex, width=args.width or "single",
                journal=args.journal, mono=args.mono, dpi=args.dpi)
        for x in salidas:
            print(f"  {x}")
        return 0

    atoms = structure.load(args.file)
    _c, rep = ep.prepare(
        atoms, outdir=args.outdir, qgrid=_malla(args.qgrid, "--qgrid"),
        kgrid_scf=_malla(args.kgrid, "--kgrid"),
        kgrid_nscf=_malla(args.kgrid_nscf, "--kgrid-nscf"),
        pseudo_dir=args.pseudo_dir, ecutwfc=args.ecutwfc,
        ecutrho=args.ecutrho, degauss=args.degauss,
        nsigma=args.nsigma, sigma_paso=args.sigma)
    _print_prepare(rep)
    return 0


def _np_array(x):
    import numpy as _np
    return _np.array(x, dtype=float)


def _cmd_wizard(args) -> int:
    from qekit.modules import wizard as wz

    lang = getattr(args, "language", None) or i18n.get_language()

    def T(s):
        return wz.t(s, lang)

    if args.list or (not args.goal and not args.ask and not args.file):
        print(wz.report_catalogo(lang))
        if args.list:
            return 0

    atoms = None
    if args.file:
        atoms = structure.load(args.file)
        d = wz.diagnosticar(atoms, pseudo_dir=args.pseudo_dir or
                            qcfg.load()["pseudo_dir"])
        print(wz.report_diagnostico(d, lang))
        print()

    archivo = args.file or "estructura.cif"

    if args.ask:
        cands = wz.buscar(args.ask, language=lang)
        if not cands:
            print(T('No encontre nada que encaje con "{texto}".').format(
                texto=args.ask), file=sys.stderr)
            print()
            print(wz.report_catalogo(lang))
            return 1
        if len(cands) > 1:
            print(T("Lo que has escrito puede ser varias cosas. La mas "
                    "probable primero:"))
            for i, m in enumerate(cands):
                print(f"  {i + 1}. {m.clave:18s} {m.pregunta}")
            print()
            print(T("Te cuento la primera ({clave}); las otras con "
                    "--goal <clave>.").format(clave=cands[0].clave))
            print()
        print(wz.report_meta(cands[0], archivo, glosario=not args.no_glossary,
                             language=lang))
        return 0

    if args.goal:
        m = wz.metas_por_clave(lang).get(args.goal)
        if m is None:
            raise ErrorDeUso(
                T("meta '{meta}' desconocida. Las disponibles salen con "
                  "'olla-dft wizard --list'.").format(meta=args.goal))
        print(wz.report_meta(m, archivo, glosario=not args.no_glossary,
                             language=lang))
        return 0

    if args.term:
        glos = wz.glosario(lang)
        t = args.term.lower()
        if t in glos:
            print(f"{t}: {glos[t]}")
            return 0
        parecidos = [k for k in glos if t in k or k in t]
        raise ErrorDeUso(
            T("'{term}' no esta en el glosario.").format(term=args.term)
            + (T(" Quiza: {lista}").format(lista=", ".join(parecidos))
               if parecidos else
               T(" Los terminos son: {lista}").format(
                   lista=", ".join(sorted(glos)))))

    print(wz.report_catalogo(lang))
    return 0


def _cmd_start(args) -> int:
    from qekit.modules import onboarding

    result = onboarding.guide(
        project_path=args.project, structure_path=args.structure,
        goal=args.goal, name=args.name,
        interactive=not args.non_interactive, validate=not args.no_validate,
        language=getattr(args, "language", "es"))
    print(onboarding.report(result))
    return 0


def _cmd_project(args) -> int:
    """Proyecto reproducible: manifiesto, workflow, calidad y dashboard."""
    from qekit.modules import dashboard, project, quality

    action = args.action
    if action == "init":
        root, data = project.init(args.target or args.project, args.name)
        print(f"Proyecto '{data['name']}' inicializado en:\n  {root}")
        print(f"Manifiesto:\n  {root / project.PROJECT_DIR / project.MANIFEST_NAME}")
        return 0

    root, data = project.load(args.project)
    if action == "cancel":
        marker = project.cancel(root, args.reason, args.cancel_file)
        print(f"Cancelación cooperativa solicitada en:\n  {marker.resolve()}")
        print("Las tareas en curso terminan su intento; las siguientes no arrancan.")
        return 0
    if action == "resume":
        changed = project.resume(root, data, args.cancel_file)
        print(f"Proyecto reanudable. Tareas devueltas a pendiente: {changed}.")
        return 0
    if action == "add":
        record = project.add_source(root, data, args.target)
        project.save(root, data)
        print(f"Fuente registrada: {record['path']}\n  SHA-256: {record['sha256']}")
        return 0
    if action == "plan":
        tasks = project.plan(root, data, args.target or "scf", args.task_commands)
        project.save(root, data)
        print(f"Plan '{args.target or 'scf'}' guardado con {len(tasks)} tareas.")
        print(project.report_status(root, data))
        return 0
    if action in ("show", "status"):
        if action == "show":
            print(json.dumps(data, ensure_ascii=False, indent=2, sort_keys=True))
        else:
            print(project.report_status(root, data))
        return 0
    if action == "validate":
        selftest_ok = None
        advanced_fail = False
        if args.selftest:
            from qekit.modules import selftest as st
            results = st.ejecutar(verbose=False)
            selftest_ok = bool(results) and all(result.bien for result in results)
            data["metadata"]["selftest"] = {
                "passed": selftest_ok,
                "at": project._now(),
                "tests": len(results),
                "failed": sum(not result.bien for result in results),
            }
            project.save(root, data)
            print(st.report(results))
            print()
        if args.advanced:
            from qekit.modules import validation
            checks = validation.check(root, data)
            advanced_fail = any(item["level"] == "fail" for item in checks)
            data["metadata"]["advanced_validation"] = {
                "passed": not advanced_fail, "at": project._now(),
                "fails": sum(item["level"] == "fail" for item in checks),
                "warnings": sum(item["level"] == "warn" for item in checks),
            }
            project.save(root, data)
            print(validation.report(checks))
            print()
        result = quality.evaluate(root, data)
        print(quality.report(result))
        return 1 if result["fails"] or selftest_ok is False or advanced_fail else 0
    if action in ("dashboard", "report"):
        if action == "report" and args.pdf:
            from qekit.modules import report as project_report
            target = project_report.generate_pdf(root, data, args.output)
            print(f"Informe PDF escrito en:\n  {target.resolve()}")
            return 0
        if getattr(args, "both", False):
            targets = dashboard.generate_pair(
                root, data, args.output, theme=getattr(args, "theme", "auto"))
            print("Dashboards escritos en:")
            for target in targets:
                print(f"  {target.resolve()}")
        else:
            target = dashboard.generate(root, data, args.output,
                                        theme=getattr(args, "theme", "auto"),
                                        language=getattr(args, "language", "es"))
            print(f"Dashboard escrito en:\n  {target.resolve()}")
        return 0
    if action == "export":
        target = project.export_snapshot(root, data, args.output)
        print(f"Snapshot reproducible escrito en:\n  {target.resolve()}")
        return 0
    if action == "environment":
        from qekit.modules import environment
        if args.verify_environment:
            result = environment.verify(root, args.target or None)
            print(environment.report(result))
            return 0 if result["ok"] else 1
        target = environment.write(root, args.target or None)
        print(f"Bloqueo de entorno escrito en:\n  {target.resolve()}")
        print("Verifícalo después con: olla-dft project environment --verify-environment")
        return 0
    if action == "diff":
        if not args.other:
            raise ErrorDeUso("project diff necesita --other SNAPSHOT_O_PROYECTO.")
        diff = project.diff(data, args.other)
        if args.json:
            print(json.dumps(diff, ensure_ascii=False, indent=2))
        else:
            print(project.diff_report(diff))
        return 0
    if action == "ingest":
        from qekit.modules import results
        paths = [args.target] if args.target else [
            root / project.PROJECT_DIR / "artifacts"]
        result = results.ingest_project(root, data, paths, tag="project")
        print(json.dumps(result, ensure_ascii=False, indent=2))
        return 0
    if action == "run":
        results = project.run(root, data, execute=args.execute, force=args.force,
                              parallel=getattr(args, "parallel", 1),
                              retries=getattr(args, "retries", 0),
                              timeout=getattr(args, "timeout", None),
                              cancel_file=getattr(args, "cancel_file", None))
        for task, code, detail in results:
            mark = task.get("status", "pending")
            print(f"[{mark:9s}] {task['id']}: {detail.splitlines()[-1]}")
        if not args.execute:
            print("\nSimulación solamente. Añade --execute para correr las tareas.")
        return 1 if any(task.get("status") == "failed" for task, _, _ in results) else 0
    raise ErrorDeUso(f"acción de project desconocida: {action}")


def _pseudos_forzados(args) -> dict:
    """Los --pseudo EL=archivo.UPF de la linea de comandos."""
    fuera = {}
    for par in (getattr(args, "pseudo", None) or []):
        if "=" not in par:
            raise ErrorDeUso(
                f"--pseudo se escribe Elemento=archivo.UPF; recibi '{par}'.")
        el, upf = par.split("=", 1)
        fuera[el.strip().capitalize()] = upf.strip()
    return fuera


def _cmd_pseudos(args) -> int:
    from qekit.modules import pseudos as pz

    pdir = args.pseudo_dir or qcfg.load()["pseudo_dir"]
    # La tarea se valida UNA vez, aqui: si se deja al bucle, un nombre mal
    # escrito sale como si fuera un problema del primer elemento.
    if args.task not in pz.TAREAS:
        raise ErrorDeUso(
            f"tarea '{args.task}' desconocida. Opciones: "
            + ", ".join(sorted(pz.TAREAS)))
    elementos = []
    if args.file:
        elementos = list(dict.fromkeys(
            structure.load(args.file).get_chemical_symbols()))
    if args.element:
        elementos = [e.strip().capitalize()
                     for e in re.split(r"[,\s]+", args.element) if e.strip()]
    if not elementos:
        raise ErrorDeUso(
            "dime de que elemento: pasa una estructura o --element Fe,O.")

    elegidos, rc = {}, 0
    for el in elementos:
        try:
            mejor, ev = pz.elegir(el, pdir, tarea=args.task,
                                  funcional=args.functional,
                                  prefiere_ligero=args.cheap)
            elegidos[el] = mejor
        except ErrorDeUso as exc:
            print(f"--- {el} ---\n{exc}\n", file=sys.stderr)
            rc = 1
            continue
        print(pz.report(el, ev, args.task))
        print()

    if len(elegidos) > 1:
        print(pz.report_coherencia(elegidos))
        print()
    if elegidos:
        print("Para usarlos tal cual:")
        print("  " + " ".join(f"--pseudo {k}={v.nombre}"
                              for k, v in elegidos.items()))
    return rc


def _cmd_tddft(args) -> int:
    from qekit.modules import tddft as td

    if args.collect:
        # el ensanchamiento fija el umbral de detección del excitón; si no
        # se repite en la línea de --collect, collect lo lee de spectrum.in
        run = td.collect(args.outdir, metodo=args.method, gap_ip=args.gap,
                         broadening=args.broadening)
        print(td.report(run))
        print()
        comparar = None
        if args.compare:
            # α por su nombre de columna: la última columna de OPTICS.dat
            # es la reflectividad, no la absorción.
            from qekit.modules import optics as _op
            cols = _op.read_optics_dat(args.compare)
            if "alpha(1/cm)" not in cols:
                raise ErrorDeUso(
                    f"'{args.compare}' no tiene la columna 'alpha(1/cm)'; "
                    "--compare espera el OPTICS.dat de 'olla-dft optics'.")
            comparar = (cols["E(eV)"], cols["alpha(1/cm)"])
        salidas = td.export(run, args.outdir)
        if not args.no_plot:
            salidas += td.plot(
                run, str(Path(args.outdir) / "tddft"), formats=args.format,
                comparar=comparar, theme=args.template, family=args.font,
                background=args.background, palette=args.palette,
                usetex=args.usetex, width=args.width or "single",
                journal=args.journal, mono=args.mono, dpi=args.dpi)
        for x in salidas:
            print(f"  {x}")
        return 0

    atoms = structure.load(args.file)
    _c, rep = td.prepare(
        atoms, outdir=args.outdir, metodo=args.method, itermax=args.iter,
        ipol=args.pol, n_estados=args.states, emin=args.emin,
        emax=args.emax,
        broadening=(args.broadening if args.broadening is not None
                    else td.BROADENING_DEFAULT),
        extrapolation=args.extrapolation, pseudo_dir=args.pseudo_dir,
        ecutwfc=args.ecutwfc, ecutrho=args.ecutrho,
        kspacing=args.kspacing, insulator=not args.metal, nbnd=args.nbnd,
        ltammd=args.tamm_dancoff, lrpa=args.rpa,
        gamma=None if not args.gamma else True, scissor=args.scissor)
    _print_prepare(rep)
    return 0


def _cmd_ballistic(args) -> int:
    from qekit.modules import ballistic as bl

    if args.collect:
        run = bl.collect(args.outdir)
        print(bl.report(run))
        print()
        salidas = bl.export(run, args.outdir)
        if not args.no_plot:
            salidas += bl.plot(
                run, str(Path(args.outdir) / "balistico"),
                formats=args.format, theme=args.template, family=args.font,
                background=args.background, palette=args.palette,
                usetex=args.usetex, width=args.width or "single",
                journal=args.journal, mono=args.mono, dpi=args.dpi)
        for x in salidas:
            print(f"  {x}")
        return 0

    if not args.file:
        raise ErrorDeUso(
            "falta la estructura del ELECTRODO (periodico en z).")
    electrodo = structure.load(args.file)
    dispersor = structure.load(args.scatterer) if args.scatterer else None
    _c, rep = bl.prepare(
        electrodo, outdir=args.outdir, dispersor=dispersor,
        ikind=args.ikind, emin=args.emin, emax=args.emax,
        npuntos=args.points, pseudo_dir=args.pseudo_dir,
        ecutwfc=args.ecutwfc, ecutrho=args.ecutrho,
        kspacing=args.kspacing, nz1=args.nz1)
    _print_prepare(rep)
    return 0


def _cmd_doctor(args) -> int:
    if args.system or not args.path:
        from qekit.modules import health

        result = health.check(args.path or ".", project_path=args.project)
        if args.json:
            print(json.dumps(result, ensure_ascii=False, indent=2))
        else:
            print(health.report(result))
        return 0 if result["ok"] else 1
    from qekit.modules import diagnose as dg

    d = dg.diagnose(args.path, prefix=args.prefix)
    print(dg.report(d))
    if not args.no_plot and (d.scf and d.scf.accuracy):
        print()
        for f in dg.plot(
            d, str(Path(args.outdir) / "diagnostico"), formats=args.format,
            theme=args.template, family=args.font,
            background=args.background, palette=args.palette,
            usetex=args.usetex or None, width=args.width or "double",
            journal=args.journal, mono=args.mono, dpi=args.dpi,
        ):
            print(f"  {f}")
    return 1 if d.problemas else 0


def _cmd_audit(args) -> int:
    from qekit.modules import audit as au

    runs = au.collect(args.paths)
    a = au.audit(runs)
    print(au.report(a))
    if args.index:
        nuevos, act = au.index(runs, args.db)
        print(f"\nBase '{args.db}': {nuevos} nuevos, {act} actualizados.")
    return 0 if a["comparables"] else 1


def _pares_ev(texto, nombre):
    if not texto:
        return None
    fuera = {}
    for trozo in str(texto).replace(";", ",").split(","):
        trozo = trozo.strip()
        if not trozo:
            continue
        if "=" not in trozo:
            raise ErrorDeUso(
                f"{nombre} se escribe CLAVE=valor, por ejemplo OH=0.77; "
                f"recibí '{trozo}'.")
        k, _, v = trozo.partition("=")
        try:
            fuera[k.strip()] = float(v)
        except ValueError:
            raise ErrorDeUso(
                f"el valor de {k.strip()} en {nombre} tiene que ser un número "
                f"en eV; recibí '{v.strip()}'.") from None
    return fuera


def _cmd_esm(args) -> int:
    from qekit.modules import esm as em

    atoms = structure.load(args.file)
    try:
        cargas = [float(x) for x in str(args.charge).replace(";", ",").split(",")
                  if x.strip()]
    except ValueError:
        raise ErrorDeUso(
            f"--charge son números separados por coma, por ejemplo "
            f"-0.2,0,0.2; recibí '{args.charge}'.") from None
    if not cargas:
        raise ErrorDeUso("--charge necesita al menos una carga.")

    run, _c, rep = em.prepare(
        atoms, outdir=args.outdir, bc=args.bc, cargas=cargas,
        campo=args.field, esm_w=args.width_esm, nfit=args.nfit,
        pseudo_dir=args.pseudo_dir, ecutwfc=args.ecutwfc,
        ecutrho=args.ecutrho, kspacing=args.kspacing)
    if not args.collect:
        _print_prepare(rep)
        for a in run.avisos:
            print()
            print(f"AVISO: {a}")
    if not args.collect:
        if _run_or_explain(run.jobs, args, "esm") is None:
            return 0
    em.collect(run, args.outdir)
    print()
    print(em.report(run))
    print()
    for f in em.export(run, args.outdir):
        print(f"  {f}")
    if not args.no_plot:
        try:
            for f in em.plot(run, str(Path(args.outdir) / "esm"),
                             **_figure_kwargs(args)):
                print(f"  {f}")
        except Exception as exc:                            # noqa: BLE001
            print(f"  (no se pudo graficar: {exc})")
    return 0


def _temperaturas(texto, nombre="--temps"):
    """'100:800:8' o '300,500,700' -> lista de temperaturas."""
    t = str(texto).strip()
    try:
        if ":" in t:
            partes = t.split(":")
            if len(partes) != 3:
                raise ValueError
            a, b, n = float(partes[0]), float(partes[1]), int(partes[2])
            if n < 1:
                raise ValueError
            return list(np.linspace(a, b, n))
        return [float(x) for x in t.replace(";", ",").split(",") if x.strip()]
    except ValueError:
        raise ErrorDeUso(
            f"{nombre} se escribe 100:800:8 (de, a, cuántas) o "
            f"300,500,700; recibí '{texto}'.") from None


def _cmd_kappa(args) -> int:
    from qekit.modules import kappa as kp, sweep

    atoms = structure.load(args.file)
    dim = _malla(args.dim, "--dim") or (2, 2, 2)
    dim2 = _malla(args.dim_fc2, "--dim-fc2")
    temps = _temperaturas(args.temps)
    out = Path(args.outdir)

    ph = kp.preparar(atoms, dim=dim, dim_fc2=dim2, distancia=args.distance)
    s3, s2 = kp.configuraciones(ph)
    run = kp.KappaRun(formula=atoms.get_chemical_formula(), dim=tuple(dim),
                      dim_fc2=tuple(dim2) if dim2 else None,
                      distancia=args.distance, n_config=len(s3),
                      n_atomos=len(s3[0]), isotopos=args.isotopes,
                      frontera=args.grain, directorio=str(out))

    if args.model:
        run.fuente = args.model.upper()
        print(f"{len(s3)} configuraciones de {len(s3[0])} átomos para la "
              f"fc3" + (f", {len(s2)} de {len(s2[0])} para la fc2"
                        if s2 else "") + ".")
        print(f"Fuerzas con {run.fuente} (esto NO es DFT):")
        F3 = kp.fuerzas_mlip(s3, args.model)
        F2 = kp.fuerzas_mlip(s2, args.model) if s2 else None
        tc, m = kp.resolver(ph, F3, F2, malla=args.mesh, temperaturas=temps,
                            isotopos=args.isotopes, frontera_um=args.grain)
        kp.recoger(run, ph, tc, m)
    elif args.collect:
        common = sweep.prepare_common(s3[0], args.pseudo_dir, args.ecutwfc,
                                      args.ecutrho, not args.metal)
        c3 = kp.escribir_inputs(s3, out / "fc3", common,
                                kspacing=args.kspacing)
        c2 = (kp.escribir_inputs(s2, out / "fc2", common,
                                 kspacing=args.kspacing) if s2 else [])
        F3 = kp.leer_fuerzas(c3, len(s3[0]))
        F2 = kp.leer_fuerzas(c2, len(s2[0])) if c2 else None
        run.fuente = "Quantum ESPRESSO"
        tc, m = kp.resolver(ph, F3, F2, malla=args.mesh, temperaturas=temps,
                            isotopos=args.isotopes, frontera_um=args.grain)
        kp.recoger(run, ph, tc, m)
    else:
        # aislante (occupations='fixed') salvo que se pida --metal: en un
        # metal las ocupaciones fijas hacen que los scf no converjan
        common = sweep.prepare_common(s3[0], args.pseudo_dir, args.ecutwfc,
                                      args.ecutrho, not args.metal)
        if len(s3) > kp.MUCHAS_CONFIGURACIONES and not args.force:
            raise ErrorDeUso(
                f"esta supercelda pide {len(s3)} cálculos de {len(s3[0])} "
                f"átomos cada uno.\nEn un portátil eso son días. Baja --dim, "
                f"prueba antes con --model mace para\nelegir el tamaño, o "
                f"insiste con --force si sabes lo que haces.")
        c3 = kp.escribir_inputs(s3, out / "fc3", common,
                                kspacing=args.kspacing)
        if s2:
            kp.escribir_inputs(s2, out / "fc2", common,
                               kspacing=args.kspacing)
        (out / "correr.sh").write_text(
            "#!/bin/bash\nfor d in fc3/d*/ fc2/d*/; do\n"
            "  [ -d \"$d\" ] || continue\n"
            "  if [ -f \"$d/pw.out\" ] && grep -q 'JOB DONE' "
            "\"$d/pw.out\"; then continue; fi\n"
            "  (cd \"$d\" && pw.x -in pw.in > pw.out 2>&1)\n"
            "done\n", encoding="utf-8")
        print(f"--- Conductividad térmica de red: {run.formula} ---")
        print(f"Supercelda fc3 {dim[0]}×{dim[1]}×{dim[2]}: {len(s3)} "
              f"configuraciones de {len(s3[0])} átomos")
        if s2:
            print(f"Supercelda fc2 {dim2[0]}×{dim2[1]}×{dim2[2]}: "
                  f"{len(s2)} configuraciones de {len(s2[0])} átomos")
        print(f"Desplazamiento: {args.distance} Å")
        print()
        print(f"Inputs en '{out.resolve()}'.")
        print("  bash correr.sh          los lanza todos, saltándose los "
              "que ya estén")
        print("Cuando terminen, el mismo comando con --collect.")
        print()
        print("Cada configuración es un scf independiente: se pueden lanzar "
              "en paralelo o en")
        print("otra máquina, y el orden da igual mientras estén TODAS.")
        return 0

    print()
    print(kp.report(run))
    print()
    for f in kp.export(run, args.outdir):
        print(f"  {f}")
    if not args.no_plot:
        try:
            for f in kp.plot(run, str(out / "kappa"), **_figure_kwargs(args)):
                print(f"  {f}")
        except Exception as exc:                            # noqa: BLE001
            print(f"  (no se pudo graficar: {exc})")
    return 0


def _cmd_berry(args) -> int:
    from qekit.modules import berry as bp

    atoms = structure.load(args.file)
    ref = structure.load(args.reference) if args.reference else None
    desp = None
    if args.displace:
        t = str(args.displace).split(":")
        if len(t) != 2:
            raise ErrorDeUso(
                "--displace se escribe ATOMO:dx,dy,dz, con el átomo en base 1 "
                "y el vector en Å. Por ejemplo 2:0,0,0.1")
        try:
            idx = int(t[0]) - 1
            vec = [float(x) for x in t[1].split(",")]
        except ValueError:
            raise ErrorDeUso(
                f"no entiendo '{args.displace}'; se escribe 2:0,0,0.1.") \
                from None
        if len(vec) != 3:
            raise ErrorDeUso("el desplazamiento necesita tres componentes.")
        desp = (idx, vec)
    kperp = args.kperp.lower().replace("x", " ").split()
    if len(kperp) != 2:
        raise ErrorDeUso(f"--kperp son dos enteros, por ejemplo 6x6; recibí "
                         f"'{args.kperp}'.")
    try:
        kperp = [int(x) for x in kperp]
    except ValueError:
        raise ErrorDeUso(f"--kperp son enteros; recibí '{args.kperp}'.") \
            from None

    run, _c, rep = bp.prepare(
        atoms, outdir=args.outdir, gdir=args.gdir, nppstr=args.nppstr,
        kperp=kperp, referencia=ref, nlambda=args.nlambda, desplazar=desp,
        pseudo_dir=args.pseudo_dir, ecutwfc=args.ecutwfc,
        ecutrho=args.ecutrho, kgrid_scf=_malla(args.kgrid, "--kgrid"))
    if not args.collect:
        _print_prepare(rep)
    if args.run:
        print()
        bp.correr(run, pw_cmd=args.pw_cmd, nproc=args.nproc,
                  timeout=args.timeout, rehacer=getattr(args, "redo", False))
    elif not args.collect:
        print()
        print("Los inputs están listos pero no se han corrido. Para "
              "ejecutarlos:")
        print("  olla-dft berry ... --run")
        print(f"  cd {args.outdir} && bash correr.sh")
        print("Cuando terminen, el mismo comando con --collect.")
        return 0
    bp.collect(run, args.outdir)
    print()
    print(bp.report(run))
    print()
    for f in bp.export(run, args.outdir):
        print(f"  {f}")
    if not args.no_plot and len(run.puntos) > 1:
        try:
            for f in bp.plot(run, str(Path(args.outdir) / "berry"),
                             **_figure_kwargs(args)):
                print(f"  {f}")
        except Exception as exc:                            # noqa: BLE001
            print(f"  (no se pudo graficar: {exc})")
    return 0


def _cmd_teoria(args) -> int:
    from qekit.modules import theory
    lang = getattr(args, "language", "es")
    if args.all:
        salida = theory.documento(lang)
    else:
        salida = theory.texto(args.comando, lang, crudo=bool(args.output))
    if args.output:
        Path(args.output).write_text(salida, encoding="utf-8")
        print(f"Escrito en: {Path(args.output).resolve()}")
        return 0
    print(salida)
    return 0


def _cmd_sistema(args) -> int:
    from qekit.core import plataforma

    print(plataforma.informe())
    return 0


def _cmd_update(args) -> int:
    from qekit.modules import update

    return update.run(check_only=args.check, yes=args.yes, target=args.version)


def _cmd_recetas(args) -> int:
    from qekit.modules import recipes as rec

    lang = getattr(args, "language", None) or i18n.get_language()

    def T(s):
        return rec.t(s, lang)

    if args.buscar:
        res = rec.buscar(args.buscar, language=lang)
        if not res:
            print(T("No encuentro ninguna receta para «{texto}».").format(
                texto=args.buscar))
            print()
            print(rec.listar(lang))
            return 0
        if len(res) > 1:
            print(T("Lo que has escrito puede ser varias cosas. La más "
                    "probable primero:"))
            for i, r in enumerate(res[:4]):
                print(f"  {i + 1}. {r.clave:16s} {r.pregunta}")
            print()
            print(T("Te cuento la primera ({clave}); las otras con "
                    "olla-dft recetas <clave>.").format(clave=res[0].clave))
            print()
        print(rec.report(res[0], lang))
        return 0

    if not args.receta:
        print(rec.listar(lang))
        return 0

    r = rec.obtener(args.receta, lang)
    if args.script is not None:
        destino = args.script or f"{r.clave}.sh"
        texto = rec.script(r, destino, lang)
        print(T("Guion escrito en {destino} ({n} líneas).").format(
            destino=destino, n=len(texto.splitlines())))
        print(T("Sale COMENTADO y con las rutas de ejemplo: ábrelo, cambia la "
                "estructura por"))
        print(T("la tuya y córrelo paso a paso. No lo lances a ciegas."))
        return 0
    print(rec.report(r, lang))
    return 0


def _ventana(texto, nombre="--window"):
    """'-10:20' -> (-10.0, 20.0). Es una ventana de energía en eV."""
    if not texto:
        return None
    t = str(texto).replace(" ", "")
    m = re.match(r"^(-?[\d.]+):(-?[\d.]+)$", t)
    if not m:
        raise ErrorDeUso(
            f"{nombre} se escribe MIN:MAX en eV, por ejemplo -10:20; "
            f"recibí '{texto}'.")
    lo, hi = float(m.group(1)), float(m.group(2))
    if hi <= lo:
        raise ErrorDeUso(
            f"{nombre}: el máximo ({hi:g}) tiene que ser mayor que el mínimo "
            f"({lo:g}).")
    return (lo, hi)


def _rango_bandas(texto):
    """'5-8,12' -> (5,6,7,8,12). Es como se citan bandas en todas partes."""
    if not texto:
        return ()
    fuera = []
    for trozo in str(texto).replace(" ", "").split(","):
        if not trozo:
            continue
        if "-" in trozo[1:]:
            a, b = trozo.split("-", 1)
            try:
                a, b = int(a), int(b)
            except ValueError:
                raise ErrorDeUso(
                    f"no entiendo el rango de bandas '{trozo}'; se escribe "
                    f"5-8.") from None
            if b < a:
                raise ErrorDeUso(f"el rango '{trozo}' va al revés.")
            fuera += list(range(a, b + 1))
        else:
            try:
                fuera.append(int(trozo))
            except ValueError:
                raise ErrorDeUso(
                    f"'{trozo}' no es un número de banda.") from None
    return tuple(sorted(set(fuera)))


def _cmd_wannier(args) -> int:
    from qekit.modules import wannier as wn

    atoms = None
    if getattr(args, "file", None):
        atoms = structure.load(args.file)

    if args.collect:
        if atoms is None:
            raise ErrorDeUso(
                "para analizar hace falta la estructura: "
                "olla-dft wannier <archivo> --collect -o <carpeta>.")
        bd = args.dft_bands
        if bd is None:
            for cand in ("out_bandas", "out_bands"):
                if (Path(args.outdir) / cand).exists():
                    bd = str(Path(args.outdir) / cand)
                    break
        run = wn.collect(args.outdir, minimizar_=not args.no_minimize,
                         pasos=args.iterations, bandas_dft=bd,
                         puntos_por_tramo=args.points, atoms=atoms,
                         exterior=_ventana(args.window, "--window"),
                         congelada=_ventana(args.frozen, "--frozen"))
        print(wn.report(run, atoms))
        print()
        salidas = wn.export(run, args.outdir)
        if args.dos:
            import numpy as _np
            e, d = wn.dos_interpolada(run, malla=args.dos, sigma=args.sigma)
            f = Path(args.outdir) / "WANNIER_dos.dat"
            _np.savetxt(f, _np.column_stack([e, d]),
                        header=f"E(eV)   DOS({wn.DOS_UNIDADES})")
            salidas.append(str(f))
            print(f"DOS interpolada en una malla {args.dos}³ = "
                  f"{args.dos ** 3} puntos k, sin volver a tocar pw.x.")
            print()
        for f in salidas:
            print(f"  {f}")
        if not args.no_plot:
            try:
                for f in wn.plot(run, str(Path(args.outdir) / "wannier"),
                                 **_figure_kwargs(args)):
                    print(f"  {f}")
            except Exception as exc:                        # noqa: BLE001
                print(f"  (no se pudo graficar: {exc})")
        return 0

    if atoms is None:
        raise ErrorDeUso("hace falta una estructura.")
    excluir = _rango_bandas(args.exclude)
    run, _common, rep = wn.prepare(
        atoms, outdir=args.outdir, malla=_malla(args.grid, "--grid"),
        proy=args.projections, nbnd=args.bands, excluir=excluir,
        pseudo_dir=args.pseudo_dir, ecutwfc=args.ecutwfc,
        ecutrho=args.ecutrho, kgrid_scf=_malla(args.kgrid, "--kgrid"),
        insulator=args.insulator, iteraciones=args.iterations)
    _print_prepare(rep)
    if not args.run:
        print()
        print("Los inputs están listos pero no se han corrido. Para "
              "ejecutarlos:")
        print("  olla-dft wannier ... --run        (Olla-DFT lanza los cuatro "
              "pasos en orden)")
        print(f"  cd {args.outdir} && bash correr.sh")
        print("Cuando terminen, el mismo comando con --collect.")
        return 0
    print()
    wn.correr(args.outdir, pw_cmd=args.pw_cmd, nproc=args.nproc,
              pw2wan_cmd=args.pw2wan_cmd, con_bandas=not args.no_dft_bands,
              timeout=args.timeout)
    args.collect = True
    return _cmd_wannier(args)


def _cmd_topology(args) -> int:
    from qekit.modules import topology as tp

    run = tp.analyze(
        args.model, occupied=args.occupied, fermi=args.fermi,
        grid=_malla_2d(args.grid), plane=args.plane, fixed=args.fixed,
        gap_tol=args.gap_tol)
    print(tp.report(run))
    print()
    written = tp.export(run, args.outdir)
    if not args.no_plot:
        try:
            figure_args = _figure_kwargs(args)
            figure_args["width"] = args.width or "double"
            written += tp.plot(
                run, str(Path(args.outdir) / "topology"), **figure_args)
        except Exception as exc:                            # noqa: BLE001
            print(f"  (no se pudo graficar: {exc})")
    for filename in written:
        print(f"  {filename}")
    return 0


def _cmd_amorphous(args) -> int:
    from qekit.core import structure as st
    from qekit.modules import amorphous as am

    simbolos = am.formula_a_simbolos(args.formula, args.units)
    at = am.empaquetar(simbolos, densidad=args.density,
                       factor=args.min_dist or am.FACTOR_MINIMO,
                       semilla=args.seed)
    out = Path(args.outdir); out.mkdir(parents=True, exist_ok=True)
    print(f"Empaquetados {len(at)} átomos de "
          f"{at.get_chemical_formula()} a {am.densidad_de(at):.4f} g/cm³")
    print(f"  celda cúbica de {np.linalg.norm(at.cell.array[0]):.3f} Å")

    if args.pack_only:
        f = out / "empaquetado.cif"
        st.convert(at, str(f))
        print(f"\n  {f}")
        print("\nSolo empaquetado: las posiciones son aleatorias y la energía "
              "es altísima.\nHace falta el fundido y temple para que esto sea "
              "un amorfo.")
        return 0

    p = am.Protocolo(T_fundido=args.melt, T_final=args.final,
                     pasos_fundido=args.melt_steps,
                     pasos_temple=args.quench_steps,
                     pasos_recocido=args.anneal_steps, dt_fs=args.dt)
    print(f"  {p.pasos} pasos de {p.dt_fs:g} fs = {p.ps_totales:.2f} ps  |  "
          f"temple a {p.velocidad_temple:.1e} K/s")
    print()
    res = am.fundir_y_templar(at, p, modelo=args.model, semilla=args.seed,
                              traza=str(out / "traza.dat"))
    print()
    print(am.report(res))
    print()
    for f in am.export(res, args.outdir):
        print(f"  {f}")
    return 0


def _cmd_docs(args) -> int:
    from qekit.modules import docs

    n = len(docs.extraer())
    if getattr(args, "both", False):
        base = Path(args.output)
        english = base.with_name(f"{base.stem}.en{base.suffix}")
        destinos = [docs.generar(str(base), language="es"),
                    docs.generar(str(english), language="en")]
        print(f"Referencias de {n} subcomandos escritas en:")
        for destino in destinos:
            print(f"  {Path(destino).resolve()}")
    else:
        destino = docs.generar(args.output,
                               language=getattr(args, "language", "es"))
        destinos = [destino]
        print(f"Referencia de {n} subcomandos escrita en:")
        print(f"  {Path(destino).resolve()}")
    print()
    print("Es una sola página, sin conexión ni dependencias: ábrela con doble "
          "clic.\nSe genera del árbol de argparse y de los docstrings, así que "
          "vuelve a\nejecutarla después de cada cambio y siempre estará al "
          "día.")
    if args.abrir:
        import webbrowser
        webbrowser.open(f"file://{Path(destinos[0]).resolve()}")
    return 0


def _cmd_echem(args) -> int:
    from qekit.modules import echem

    if (args.her is None) == (args.oer is None):
        raise ErrorDeUso(
            "elige una reacción: --her E_ads(H)  o  --oer OH=..,O=..,OOH=.. "
            "(una de las dos, no las dos ni ninguna).")
    corr = _pares_ev(args.corrections, "--corrections")
    if args.her is not None:
        e = echem.her(args.her, correccion=(corr or {}).get("H"),
                      T=args.temperature)
    else:
        e = echem.oer(_pares_ev(args.oer, "--oer"), correcciones=corr,
                      T=args.temperature)
    e.U, e.pH = args.potential, args.ph
    print(echem.report(e))
    print()
    for f in echem.export(e, args.outdir):
        print(f"  {f}")
    if not args.no_plot:
        try:
            for f in echem.plot(e, str(Path(args.outdir) / "echem"),
                                **_figure_kwargs(args)):
                print(f"  {f}")
        except Exception as exc:                            # noqa: BLE001
            print(f"  (no se pudo graficar: {exc})")
    return 0


def _cmd_selftest(args) -> int:
    from qekit.modules import selftest as st

    if args.list:
        print("Pruebas de validación:")
        for p in st.PRUEBAS:
            marca = ("pw.x" if p.necesita_qe
                     else "mlip" if getattr(p, "necesita_mlip", False)
                     else "----")
            print(f"  [{marca}] {p.clave:16s} {p.titulo}")
            print(f"           {p.magnitud} = {p.referencia:g} {p.unidad}"
                  f"  ±{p.tolerancia * 100:.0f} %   ({p.coste})")
            print(f"           {p.fuente}")
        return 0

    claves = None
    if args.only:
        claves = [c.strip() for c in args.only.split(",") if c.strip()]
        conocidas = {p.clave for p in st.PRUEBAS}
        malas = [c for c in claves if c not in conocidas]
        if malas:
            raise ErrorDeUso(
                f"no conozco la prueba {', '.join(malas)}. "
                f"Las que hay: {', '.join(sorted(conocidas))}.")
    if args.full and not args.pseudo_dir:
        raise ErrorDeUso(
            "las pruebas con --full corren pw.x de verdad y necesitan "
            "pseudopotenciales: pásalos con --pseudo-dir.")

    extras = []
    if args.full:
        extras.append("pw.x")
    if args.mlip:
        extras.append("MLIP")
    print("Corriendo las pruebas"
          + (f" (incluidas las de {' y '.join(extras)})" if extras
             else " rápidas") + " ...")
    res = st.ejecutar(claves=claves, con_qe=args.full, con_mlip=args.mlip,
                      pseudo_dir=args.pseudo_dir, pw_cmd=args.pw_cmd,
                      nproc=args.nproc, paralelo=args.jobs,
                      carpeta=args.keep)
    print()
    print(st.report(res))
    fallos = [r for r in res if not r.bien]
    return 1 if fallos else 0


def _cmd_cost(args) -> int:
    from qekit.modules import cost

    modelo = cost.calibrar(args.db)
    print(cost.report_modelo(modelo))
    return 0 if modelo.calibrado else 1


def _cmd_compare(args) -> int:
    from qekit.modules import compare

    result = compare.compare(args.paths, reference=args.reference)
    print(compare.report(result))
    if args.output:
        print(f"\nJSON escrito en: {compare.export(result, args.output).resolve()}")
    return 0


def _cmd_tune(args) -> int:
    from qekit.modules import tuning

    result = tuning.analyze(args.file, args.threshold)
    print(tuning.report(result))
    if args.output:
        print(f"\nJSON escrito en: {tuning.export(result, args.output).resolve()}")
    return 0


def _cmd_results(args) -> int:
    from qekit.modules import project, results

    if args.action == "explore" and args.db:
        root = Path(args.project).expanduser().resolve()
        data = {"name": root.name}
    else:
        root, data = project.load(args.project)
    db = Path(args.db).expanduser() if args.db else results.project_db(root)
    action = args.action
    if action == "ingest":
        paths = []
        if args.target:
            paths.append(args.target)
        paths.extend(args.extra_paths or [])
        if not paths:
            paths = [root / project.PROJECT_DIR / "artifacts"]
        result = results.ingest_project(root, data, paths, tag=args.tag, db_path=db)
        print(json.dumps(result, ensure_ascii=False, indent=2))
        return 0
    if action == "list":
        rows = results.list_results(
            db, formula=args.formula, calculation=args.calculation,
            status=args.status, limit=args.limit if args.limit is not None else 100)
        if args.json:
            print(json.dumps(rows, ensure_ascii=False, indent=2))
        else:
            print(results.report(rows, db))
        return 0
    if action == "show":
        if not args.target:
            raise ErrorDeUso("results show necesita un id de resultado.")
        print(json.dumps(results.get(db, args.target), ensure_ascii=False, indent=2))
        return 0
    if action == "review":
        if not args.target:
            raise ErrorDeUso("results review necesita un id de resultado.")
        if not args.review_status:
            raise ErrorDeUso("results review necesita --review-status.")
        row = results.review(db, args.target, args.review_status, args.note)
        print(json.dumps(row, ensure_ascii=False, indent=2))
        return 0
    if action == "explore":
        from qekit.modules import studio
        limit = args.limit if args.limit is not None else 10000
        if not 1 <= limit <= 10000:
            raise ErrorDeUso("--limit: 1–10000")
        filters = dict(formula=args.formula, calculation=args.calculation, status=args.status)
        rows = results.list_results(db, limit=limit, **filters)
        target = args.output or (root / project.PROJECT_DIR / "reports" / "results.html")
        target = studio.generate(rows, target, title=data.get("name", root.name),
                                 language=args.language, total_count=results.count_results(db, **filters),
                                 order="ingested_desc_path_asc")
        print(str(target.resolve()))
        return 0
    if action == "export":
        target = args.output or (root / project.PROJECT_DIR / "reports" /
                                 "results.json")
        print(f"Resultados exportados en:\n  {results.export(db, target).resolve()}")
        return 0
    raise ErrorDeUso(f"acción de results desconocida: {action}")


def _cmd_campaign(args) -> int:
    from qekit.modules import campaign, project

    root, data = project.load(args.project)
    action = args.action
    if action == "create":
        if not args.target:
            raise ErrorDeUso("campaign create necesita un nombre.")
        record = campaign.create(
            root, data, args.target, args.campaign_command, args.axis,
            goal=args.goal, convergence_file=args.convergence_file,
            adaptive=args.adaptive)
        project.save(root, data)
        print(f"Campaña '{record['id']}' guardada con {record['points']} puntos.")
        print(campaign.report(data, record["id"]))
        print("Nada se ejecutó. Revisa las tareas y usa 'olla-dft project run --execute'.")
        return 0
    if action == "list":
        print(campaign.report(data))
        return 0
    if action == "status":
        if not args.target:
            raise ErrorDeUso("campaign status necesita un id.")
        print(campaign.report(data, args.target))
        return 0
    if action == "run":
        if not args.target:
            raise ErrorDeUso("campaign run necesita un id.")
        runs = campaign.run(root, data, args.target, execute=args.execute,
                            force=args.force, parallel=args.parallel,
                            retries=args.retries, timeout=args.timeout,
                            cancel_file=args.cancel_file)
        for task, code, detail in runs:
            print(f"[{task.get('status', 'pending'):9s}] {task['id']}: "
                  f"{detail.splitlines()[-1]}")
        if not args.execute:
            print("\nSimulación solamente. Añade --execute para ejecutar la campaña.")
        return 1 if any(task.get("status") == "failed" for task, _, _ in runs) else 0
    if action == "extend":
        if not args.target or not args.convergence_file:
            raise ErrorDeUso("campaign extend necesita id y --convergence-file.")
        result = campaign.extend(root, data, args.target, args.convergence_file,
                                 threshold=args.threshold)
        project.save(root, data)
        print(f"Valor recomendado: {result['recommended_value']:g}")
        print(f"Puntos añadidos: {result.get('points_added', 0)}")
        if not result["extended"]:
            print(result["reason"])
        return 0
    if action == "export":
        target = args.output or (root / project.PROJECT_DIR / "reports" /
                                 "campaigns.json")
        print(f"Campañas exportadas en:\n  "
              f"{campaign.export(data, target, args.target).resolve()}")
        return 0
    raise ErrorDeUso(f"acción de campaign desconocida: {action}")


def _print_db_rows(rows):
    if not rows:
        print("(sin resultados)")
        return
    cols = list(rows[0])
    print("  ".join(f"{c:>20s}" for c in cols))
    for row in rows:
        print("  ".join(f"{str(row[c])[:20]:>20s}" for c in cols))


def _cmd_db(args) -> int:
    from qekit.modules import audit as au

    if args.paths:
        runs = au.collect(args.paths)
        nuevos, act = au.index(runs, args.db)
        print(f"Base '{args.db}': {nuevos} nuevos, {act} actualizados.")
        print()
    if args.query:
        filas = au.query(args.query, args.db)
        if not filas:
            print("(sin resultados)")
            return 0
        cols = list(filas[0])
        print("  ".join(f"{c:>16s}" for c in cols))
        for f in filas:
            print("  ".join(f"{str(f[c])[:16]:>16s}" for c in cols))
        return 0
    if any(value is not None for value in
           (args.formula, args.calculation, args.gap_min, args.gap_max)):
        _print_db_rows(au.search(
            args.db, formula=args.formula, calculation=args.calculation,
            gap_min=args.gap_min, gap_max=args.gap_max, limit=args.limit))
        return 0
    if args.export:
        print("  " + au.export_json(args.db, args.export))
        return 0
    print(au.summary(args.db))
    return 0


def _cmd_hull(args) -> int:
    from qekit.modules import audit as au
    from qekit.modules import thermo as th

    runs = au.collect(args.paths)
    a = au.audit(runs)
    if not a["comparables"]:
        print(au.report(a), file=sys.stderr)
        print("\nNo se construye el casco: con parámetros distintos las "
              "energías de formación\nno significan nada. Corrige lo de "
              "arriba, o usa --force si sabes lo que haces.",
              file=sys.stderr)
        if not args.force:
            return 1
        print()
    res = th.from_runs(a["buenos"], elementos=(args.elements.split(",")
                                               if args.elements else None))
    print(th.report(res, umbral=args.threshold))
    print()
    for f in th.export(res, args.outdir):
        print(f"  {f}")
    if not args.no_plot and len(res.elementos) == 2 and not res.faltan_ref:
        for f in th.plot(
            res, str(Path(args.outdir) / "casco"), formats=args.format,
            theme=args.template, family=args.font,
            background=args.background, palette=args.palette,
            usetex=args.usetex or None, width=args.width or "single",
            journal=args.journal, mono=args.mono, dpi=args.dpi,
        ):
            print(f"  {f}")
    return 0


def _cmd_report(args) -> int:
    from qekit.modules import feedback as fb

    if args.stats:
        print(fb.report_estadisticas(fb.estadisticas()))
        return 0
    if args.export:
        print("  " + fb.exportar(args.export,
                                 solo_abiertas=args.only_open))
        print("\nEse archivo lleva todo lo necesario para reproducir cada "
              "fallo:\ncomando, traza y versiones. Es lo que hay que "
              "entregar para que se arregle.")
        return 0
    if args.close:
        if fb.cerrar(args.close, nota=args.note or ""):
            print(f"Incidencia {args.close} cerrada.")
            return 0
        print(f"No existe la incidencia '{args.close}'.", file=sys.stderr)
        return 1
    if args.show:
        for i in fb.listar():
            if i.id == args.show:
                print(fb.report_detalle(i))
                return 0
        print(f"No existe la incidencia '{args.show}'.", file=sys.stderr)
        return 1
    if args.description:
        inc = fb.registrar(" ".join(args.description),
                           adjuntos=args.attach or [])
        print(f"Incidencia {inc.id} registrada en "
              f"{fb.DIR / inc.id}")
        if inc.adjuntos:
            print(f"  Adjuntos copiados: {', '.join(inc.adjuntos)}")
        print("\nSe guardó localmente. Olla-DFT no manda nada a ningún lado.")
        return 0
    print(fb.report_lista(fb.listar()))
    return 0


def _cmd_mlip(args) -> int:
    from qekit.modules import mlip as ml

    atoms = structure.load(args.file)
    if args.action == "relax":
        run = ml.relax(atoms, modelo=args.model, fmax=args.fmax,
                       steps=args.steps, cell=not args.fixed_cell,
                       device=args.device, model_size=args.size)
        print(ml.report_relax(run))
        destino = args.output or "relajado_mlip.cif"
        structure.convert(run.atoms_final, destino)
        marca = ml.write_provenance(run, destino)
        print(f"\n  {destino}\n  {marca}")
        return 0
    if args.action == "scan":
        d = ml.volume_scan(atoms, modelo=args.model, span=args.span,
                           npoints=args.npoints, device=args.device,
                           model_size=args.size)
        print(ml.report_scan(d, atoms))
        return 0
    d = ml.phonon_check(atoms, modelo=args.model,
                        supercell=_malla(args.supercell or "2x2x2",
                                         "--supercell"),
                        device=args.device, model_size=args.size)
    print(ml.report_phonon(d))
    return 0 if d["estable"] else 1


def _cmd_suggest(args) -> int:
    from qekit.modules import audit as au
    from qekit.modules import recommend as rc

    atoms = structure.load(args.file)
    els = list(dict.fromkeys(atoms.get_chemical_symbols()))
    try:
        filas = au.query("SELECT * FROM calculos", args.db)
    except FileNotFoundError as exc:
        print(f"{exc}", file=sys.stderr)
        return 1
    # ¿es una losa? mucho vacío en una dirección
    z = atoms.positions[:, 2]
    c = float(atoms.cell.array[2, 2])
    es_losa = bool(c > 0 and (c - (z.max() - z.min())) > 8.0)
    sug = rc.sugerir(filas, els, natoms=len(atoms), es_losa=es_losa)
    print(rc.report(sug, els, n_historial=len(filas)))
    return 0


def _cmd_crosscheck(args) -> int:
    from qekit.modules import crosscheck as cc

    kw = {}
    if args.file:
        atoms = structure.load(args.file)
        kw["masas"] = atoms.get_masses()
        kw["volumen"] = atoms.get_volume()
        kw["natoms"] = len(atoms)
        kw["cell"] = atoms.cell.array
    if args.gap_bandas is not None:
        kw["gap_bandas"] = args.gap_bandas
    if args.gap_tauc is not None:
        kw["gap_tauc"] = args.gap_tauc
    res = cc.run(project=args.project, **kw)
    print(cc.report(res))
    return 1 if any(c.ok is False for c in res.checks) else 0


def _cmd_derived(args) -> int:
    import numpy as np

    from qekit.modules import derived, elastic

    atoms = structure.load(args.file)
    C = np.loadtxt(args.cij, comments="#")
    if C.shape != (6, 6):
        print(f"'{args.cij}' no contiene una matriz 6x6", file=sys.stderr)
        return 1
    m = elastic.moduli(C)
    r = derived.analyze(m.B_hill, m.G_hill, atoms.get_masses(),
                        atoms.get_volume(), natoms=len(atoms), T=args.temp)
    print(derived.report(r))
    # v_L = sqrt(C11/rho) solo tiene sentido en un cristal cúbico: se mira la
    # simetría de la estructura y, por si el CIF viene sin ella, la forma
    # del propio tensor
    try:
        cubico = elastic.crystal_family(
            structure.symmetry_dataset(atoms).number) == "cúbico"
    except Exception:
        cubico = False
    d = (derived.cubic_directional(C, r.rho)
         if cubico or derived.is_cubic_tensor(C) else {})
    if d:
        print("\nEn un cristal cúbico, a lo largo de [100]:")
        print(f"  v_L = √(C₁₁/ρ) = {d['v_l_100']:.0f} m/s")
        print(f"  v_T = √(C₄₄/ρ) = {d['v_t_100']:.0f} m/s")
        print("  (son estas, no los promedios isótropos, las que se "
              "comparan\n   contra la pendiente de las ramas acústicas)")
    print()
    for f in derived.export(r, getattr(args, "outdir", ".") or "."):
        print(f"  {f}")
    return 0


def _cmd_datasheet(args) -> int:
    from qekit.modules import datasheet as ds

    f = ds.recoger(args.project)
    if args.methods:
        print(ds.metodos(f))
        print()
        print("Referencias:")
        for i, c in enumerate(ds.citas(f), start=1):
            print(f"  {i}. {c}")
        return 0
    if not f.resultados:
        print(f"No se encontró ningún resultado de Olla-DFT en "
              f"'{args.project}'.", file=sys.stderr)
        return 1
    print(f"Ficha de {f.formula or '?'}: "
          f"{len(f.resultados)} secciones con resultados")
    for s in f.resultados:
        print(f"  - {s}")
    for a in f.avisos:
        print(f"\nAVISO: {a}")
    print()
    for x in ds.escribir(f, args.outdir, args.name):
        print(f"  {x}")
    return 0


def _cmd_qha(args) -> int:
    import numpy as np

    from qekit.modules import qha as q

    datos = np.loadtxt(args.data, comments="#")
    if datos.ndim != 2 or datos.shape[1] < 3:
        print(f"'{args.data}' debe tener columnas: V(A^3) E(eV) w1 w2 ...",
              file=sys.stderr)
        return 1
    V, E, F = datos[:, 0], datos[:, 1], [f[f > -1e3] for f in datos[:, 2:]]
    # con la estructura se sabe cuántas celdas primitivas caben en la
    # convencional, y a(T) sale como parámetro de red de verdad; sin ella
    # solo es V_prim^(1/3) y el informe lo dice
    cubico, factor = args.cubic, None
    if getattr(args, "structure", None):
        atoms = structure.load(args.structure)
        factor = q.factor_convencional(atoms)
        cubico = cubico or q.es_cubico(atoms)
    res = q.run(V, E, F, T=np.arange(0.0, args.tmax + 1, args.dt),
                natoms=args.natoms, cubico=cubico,
                celdas_por_modo=args.cells, factor_conv=factor)
    print(q.report(res, T_ref=args.temp))
    print()
    for f in q.export(res, args.outdir):
        print(f"  {f}")
    if not args.no_plot:
        for f in q.plot(res, str(Path(args.outdir) / "qha"),
                        formats=args.format, theme=args.template,
                        family=args.font, background=args.background,
                        palette=args.palette, usetex=args.usetex or None,
                        width=args.width or "double", journal=args.journal,
                        mono=args.mono, dpi=args.dpi):
            print(f"  {f}")
    return 0


def _cmd_transport(args) -> int:
    import glob as _glob

    from qekit.modules import transport as tr

    atoms = structure.load(args.file)
    grid = _malla(args.grid, "--grid") or (16, 16, 16)

    if not args.collect:
        nspin = getattr(args, "nspin", 1)
        mag = _parse_mag(getattr(args, "mag", None),
                         atoms.get_chemical_symbols())
        if mag:
            nspin = 2       # pedir magnetización implica activar el espín
        _grid, rep = tr.prepare(
            atoms, outdir=args.outdir, pseudo_dir=args.pseudo_dir,
            ecutwfc=args.ecutwfc, ecutrho=args.ecutrho, grid=grid,
            insulator=not args.metal, kspacing=args.kspacing,
            nspin=nspin, magnetization=mag)
        _print_prepare(rep)
        if args.run:
            from qekit.core import runner as run_mod
            print()
            jobs = [run_mod.Job(name="scf", directory=Path(args.outdir),
                                input_file="scf.in", output_file="scf.out"),
                    run_mod.Job(name="nscf", directory=Path(args.outdir),
                                input_file="nscf.in", output_file="nscf.out")]
            res = run_mod.run_all(jobs, pw_cmd=args.pw_cmd,
                                  nproc=args.nproc, timeout=args.timeout)
            if not all(r.ok for r in res):
                return 1
        else:
            print("\nCuando termine, vuelve con --collect.")
            return 0

    xmls = _glob.glob(str(Path(args.outdir) / "out" / "*.xml"))
    if not xmls:
        print(f"Error: no hay XML en {args.outdir}/out", file=sys.stderr)
        return 1
    temps = [float(t) for t in args.temperatures.split(",")]
    run = tr.load(xmls[0])
    run = tr.compute(run, T=temps, mu_span=args.mu_span)
    print()
    print(tr.report(run, t=temps[0]))
    print()
    print(tr.report_lorenz(run, t=temps[0]))
    print()
    if getattr(args, "spin_resolved", False):
        from qekit.core import qeout as _qeout
        res = _qeout.read_xml(xmls[0])
        if res.nspin != 2:
            raise ErrorDeUso(
                "--spin-resolved necesita un cálculo con polarización de "
                "espín, y este tiene nspin = 1. Vuelve a preparar y correr "
                "con\n'olla-dft transport ESTRUCTURA --nspin 2 --mag "
                "EL=0.7 --run' (y luego\n--collect --spin-resolved).")
        dw = tr.compute(tr.load(xmls[0], spin=1), T=temps,
                        mu_span=args.mu_span)
        te = tr.TransporteEspin(up=run, dw=dw,
                                it=tr.at_temperature(run, temps[0]))
        print(tr.report_espin(te, t=temps[0]))
        print()
    for f in tr.export(run, args.outdir, t=temps[0]):
        print(f"  {f}")
    if not args.no_plot:
        for f in tr.plot(
            run, str(Path(args.outdir) / "transporte"), formats=args.format,
            theme=args.template, family=args.font, background=args.background,
            palette=args.palette, usetex=args.usetex or None,
            width=args.width or "double", journal=args.journal,
            mono=args.mono, dpi=args.dpi, t=temps[0],
        ):
            print(f"  {f}")
    return 0


def _cmd_phonons_tscan(args, atoms, qgrid) -> int:
    from qekit.core import runner as run_mod
    from qekit.modules import phonons as ph_mod
    from qekit.modules import tphonons as tp_mod

    try:
        temps = [float(x) for x in args.tscan.replace(";", ",").split(",")
                 if x.strip()]
    except ValueError:
        raise ErrorDeUso(
            f"--tscan son temperaturas en K separadas por coma, por ejemplo "
            f"300,1000,3000; recibí '{args.tscan}'.") from None

    run, rep = tp_mod.prepare(
        atoms, temps, outdir=args.outdir, gamma_only=args.gamma,
        pseudo_dir=args.pseudo_dir, ecutwfc=args.ecutwfc,
        ecutrho=args.ecutrho, kspacing=args.kspacing, qgrid=qgrid,
    )
    _print_prepare(rep)
    if args.run:
        for T, carpeta in zip(run.temperaturas, run.carpetas):
            print(f"\n--- {T:g} K ---")
            pr = ph_mod.PhononRun(prefix="", outdir=Path(carpeta),
                                  gamma_only=args.gamma)
            job = run_mod.Job(name=f"scf {T:g} K", directory=Path(carpeta),
                              input_file="scf.in", output_file="scf.out")
            res = run_mod.run_all([job], pw_cmd=args.pw_cmd,
                                  nproc=args.nproc, timeout=args.timeout,
                                  paralelo=getattr(args, "jobs", 1) or 1)
            if not all(r.ok for r in res):
                print(f"  el scf de {T:g} K falló; se salta esta temperatura")
                continue
            ph_mod.run_chain(pr, pw_cmd=args.pw_cmd, nproc=args.nproc)
    elif not args.collect:
        print("\nCorre con --run, o ejecuta las cadenas a mano y vuelve con "
              "--collect.")
        return 0
    tp_mod.collect(run)
    print()
    print(tp_mod.report(run))
    print()
    for f in tp_mod.export(run, args.outdir):
        print(f"  {f}")
    if not args.no_plot:
        try:
            for f in tp_mod.plot(run, str(Path(args.outdir) / "fonones_T"),
                                 **_figure_kwargs(args)):
                print(f"  {f}")
        except Exception as exc:                            # noqa: BLE001
            print(f"  (no se pudo graficar: {exc})")
    return 0


def _cmd_phonons(args) -> int:
    from qekit.core import runner as run_mod
    from qekit.modules import phonons as ph_mod

    atoms = structure.load(args.file)
    qgrid = _malla(args.qgrid, "--qgrid")
    if getattr(args, "tscan", None):
        return _cmd_phonons_tscan(args, atoms, qgrid)
    run, rep = ph_mod.prepare(
        atoms, outdir=args.outdir, pseudo_dir=args.pseudo_dir,
        ecutwfc=args.ecutwfc, ecutrho=args.ecutrho, kspacing=args.kspacing,
        insulator=args.insulator, qgrid=qgrid,
        gamma_only=args.gamma or args.raman, raman=args.raman,
    )
    _print_prepare(rep)
    natoms = len(structure.primitive(atoms))
    if args.run:
        print()
        results = run_mod.run_all(run.jobs, pw_cmd=args.pw_cmd,
                                  nproc=args.nproc, timeout=args.timeout)
        if not all(r.ok for r in results):
            return 1
        ph_mod.run_chain(run, pw_cmd=args.pw_cmd, nproc=args.nproc)
    elif not args.collect:
        print("\nCorre con --run (puede tardar), o ejecuta la cadena a mano "
              "y vuelve con --collect.")
        return 0
    ph_mod.collect(run)
    print()
    if run.gamma_only and run.modes:
        print(ph_mod.report_gamma_activities(run))
        if run.raman:
            w, inten, picos = ph_mod.raman_spectrum(run, laser_nm=args.laser)
            print()
            print(f"Espectro Raman simulado (láser {args.laser:g} nm, "
                  "300 K, Lorentz 5 cm⁻¹):")
            vistos = []
            for wp, ip in sorted(picos, key=lambda t: -t[1]):
                if any(abs(wp - v) < 1.0 for v in vistos):
                    continue
                vistos.append(wp)
                rel = ip / max(x[1] for x in picos) * 100.0
                print(f"  {wp:9.2f} cm⁻¹   I rel = {rel:6.1f}")
    else:
        print(ph_mod.report(run, natoms=natoms))
    print()
    for f in ph_mod.export(run, args.outdir, natoms=natoms):
        print(f"  {f}")
    if args.suite:
        if not run.gamma_only:
            print("  (--suite solo aplica con --gamma: el JSON de "
                  "intercambio lleva\n   frecuencias y actividad IR en Γ, "
                  "que es lo comparable con FTIR/Raman)")
        else:
            from qekit.modules import interop
            if run.raman and run.modes:
                doc = interop.from_raman(run, atoms, laser_nm=args.laser)
                nombre = "RAMAN_suite.json"
            else:
                doc = interop.from_phonons_gamma(run, atoms)
                nombre = "FONONES_suite.json"
            print("  " + interop.write(doc, Path(args.outdir) / nombre))
    # --raman fuerza gamma_only aunque no se pase --gamma: se decide por
    # el estado real de la corrida, no por la bandera
    if not args.no_plot and ph_mod.has_dispersion(run):
        for f in ph_mod.plot(
            run, str(Path(args.outdir) / "fonones"), formats=args.format,
            theme=args.template, family=args.font, background=args.background,
            palette=args.palette, usetex=args.usetex or None,
            width=args.width or "double", journal=args.journal,
            mono=args.mono, dpi=args.dpi,
        ):
            print(f"  {f}")
    return 0


def _cmd_config(args) -> int:
    if args.action == "show" or args.action is None:
        print(qcfg.show())
    elif args.action == "set":
        if not args.key or args.value is None:
            print("uso: olla-dft config set <clave> <valor>", file=sys.stderr)
            return 1
        qcfg.set_value(args.key, args.value)
        print(f"{args.key} = {args.value} guardado.")
    return 0


def build_parser(language=None) -> argparse.ArgumentParser:
    """Construye el árbol de argparse en el idioma pedido (es por defecto).

    Los textos de ayuda se escriben en español en el código; al final se
    completan las banderas que no traen ayuda propia y, si el idioma es
    inglés, se traducen todos con la tabla de ``qekit/data/i18n/cli_en.json``.
    """
    language = language or i18n.get_language()
    parser = argparse.ArgumentParser(
        prog=__command_name__,
        description=i18n.translate(
            f"{__product_name__} — toolkit para Quantum ESPRESSO", language),
        epilog=_catalog_text(language),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("--version", action="version",
                        version=f"{__command_name__} {__version__}",
                        help="mostrar la versión y salir")
    parser.add_argument(
        "--ascii", action="store_true",
        help="salida solo en ASCII (Å -> A, α -> alpha, → -> ->). Útil si tu "
             "terminal no admite UTF-8, o si vas a redirigir a un archivo")
    parser.add_argument(
        "--language", dest="language", choices=["es", "en"], default=language,
        help="idioma de la interfaz: es o en. También vale la variable "
             "OLLA_DFT_LANG o 'olla-dft config set language en'")
    sub = parser.add_subparsers(dest="command", metavar="COMANDO")

    _add_gen_parser(sub)

    p = sub.add_parser("info", help="información de estructura y simetría")
    p.add_argument("file")

    p = sub.add_parser("kpath", help="camino de alta simetría (seekpath)")
    p.add_argument("file")

    p = sub.add_parser("prim", help="celda primitiva estandarizada")
    p.add_argument("file")
    p.add_argument("-o", "--output", default="primitive.cif")

    p = sub.add_parser("conv", help="celda convencional estandarizada")
    p.add_argument("file")
    p.add_argument("-o", "--output", default="conventional.cif")

    p = sub.add_parser("supercell", help="construir supercelda")
    p.add_argument("file")
    p.add_argument("nx", type=int)
    p.add_argument("ny", type=int)
    p.add_argument("nz", type=int)
    p.add_argument("-o", "--output", default="supercell.cif")

    p = sub.add_parser("convert", help="convertir formato (CIF/POSCAR/XYZ)")
    p.add_argument("file")
    # el destino admite las dos formas: posicional (natural aquí) y -o, que
    # es como lo piden prim/conv/supercell — así ninguna de las dos falla
    p.add_argument("output", nargs="?")
    p.add_argument("-o", "--output-flag", dest="output_flag",
                   help="archivo de salida (alternativa a darlo posicional)")

    # ---------------- post-proceso ----------------
    def _post_opts(p, with_mode=False):
        p.add_argument("path", nargs="?", default=".",
                       help="carpeta del cálculo (o ruta al .xml)")
        if with_mode:
            p.add_argument("--mode", default="orbital",
                           choices=["orbital", "element", "total"],
                           help="cómo descomponer la PDOS")
        p.add_argument("-o", "--outdir", default=".", help="carpeta de salida")
        p.add_argument("--prefix", help="prefix del cálculo (se detecta solo)")
        p.add_argument("--ref", default="auto",
                       choices=["auto", "fermi", "vbm", "none"],
                       help="origen de energías (default: auto)")
        p.add_argument("--emin", type=float, default=-6.0,
                       help="límite inferior del eje de energía (eV)")
        p.add_argument("--emax", type=float, default=6.0,
                       help="límite superior del eje de energía (eV)")
        p.add_argument("--no-plot", action="store_true",
                       help="solo exportar datos, sin generar la gráfica")
        # --- apariencia de la figura: en su propio grupo de la ayuda ---
        p = p.add_argument_group(GRUPO_FIGURA)
        p.add_argument("--dpi", type=int, default=600,
                       help="resolución de los formatos de mapa de bits")
        p.add_argument("--format", default="pdf,png",
                       help="formatos separados por coma: pdf,png,svg,eps,tif")
        p.add_argument("-t", "--template", default=None,
                       help="plantilla visual: " + ", ".join(qthemes.names())
                            + " (o la ruta a un JSON propio)")
        p.add_argument("--size", default=None, choices=sorted(qstyle.STYLES),
                       help="escala tipográfica: paper / presentation / poster")
        p.add_argument("--font", default=None,
                       choices=["sans", "serif", "latex"],
                       help="familia tipográfica (latex = Computer Modern)")
        p.add_argument("--usetex", action="store_true",
                       help="renderizar el texto con LaTeX de verdad")
        p.add_argument("--palette", default=None,
                       help="paleta: " + ", ".join(sorted(qthemes.PALETTES))
                            + ", o colores hexadecimales separados por coma")
        p.add_argument("--background", default=None,
                       help="color de fondo, por ejemplo '#FFFFFF' o 'none'")
        p.add_argument("--journal", default="generic",
                       choices=sorted(qstyle.JOURNALS),
                       help="anchos de columna de la editorial")
        p.add_argument("--width", default=None,
                       help="ancho: single / onehalf / double, o un número en mm")
        p.add_argument("--aspect", type=float, default=None,
                       help="relación alto/ancho de la figura")
        p.add_argument("--mono", action="store_true",
                       help="monocromo: tinta negra y patrones de línea "
                            "(para revistas que cobran el color)")
        p.add_argument("--dashes", default="auto",
                       choices=["auto", "always", "never"],
                       help="patrones de línea como codificación secundaria")
        p.add_argument("--title", default=None,
                       help="título dentro de la figura (por defecto ninguno: "
                            "en un artículo el texto va en el pie de figura)")
        p.add_argument("--gap-label", action="store_true",
                       help="anotar el valor del gap dentro de la gráfica")
        p.add_argument("--panel", default=None,
                       help="etiqueta de panel, por ejemplo '(a)'")

    p = sub.add_parser("bands", help="analizar y graficar la estructura de bandas")
    _post_opts(p)

    p.add_argument("--fat", metavar="SELECTOR",
                   help="fatbands: peso de un orbital sobre cada banda. "
                        "Por ejemplo Ni-d, Si-p, O, d o atomo:3. Necesita la "
                        "salida de projwfc.x del MISMO cálculo de bandas")
    p.add_argument("--fat-scale", dest="fat_scale", type=float, default=55.0,
                   help="tamaño de los puntos de las fatbands")
    p.add_argument("--projwfc", metavar="ARCHIVO",
                   help="salida de projwfc.x (por omisión projwfc.out en la "
                        "misma carpeta)")

    p = sub.add_parser("dos", help="analizar y graficar DOS y PDOS")
    _post_opts(p, with_mode=True)

    p.add_argument("--dband", metavar="EL[-ORB]",
                   help="centro, anchura y llenado de una banda proyectada, "
                        "por ejemplo Pt (usa d) o Ni-p. Es el descriptor que "
                        "se correlaciona con la energía de adsorción")
    p.add_argument("--dband-emax", dest="dband_emax", type=float, default=None,
                   metavar="eV",
                   help="corte superior de la integral, respecto al Fermi")

    p = sub.add_parser("plot", help="gráfica combinada de bandas + DOS")
    _post_opts(p, with_mode=True)

    p = sub.add_parser("gap", help="solo el reporte de band gap (rápido)")
    p.add_argument("path", nargs="?", default=".",
                   help="carpeta del cálculo (o ruta al .xml)")
    p.add_argument("--prefix", help="prefix del cálculo (se detecta solo)")

    # ---------------- módulos de cálculo ----------------
    def _calc_opts(parser_, default_out):
        p = parser_
        p.add_argument("file", help="estructura (CIF, POSCAR, input de pw.x...)")
        p.add_argument("-o", "--outdir", default=default_out,
                       help="carpeta del barrido")
        p = parser_.add_argument_group(GRUPO_EJECUCION)
        p.add_argument("--run", action="store_true",
                       help="ejecutar los cálculos ahora, uno tras otro")
        p.add_argument("--collect", action="store_true",
                       help="solo analizar cálculos ya corridos")
        p.add_argument("--pw-cmd", dest="pw_cmd",
                       help="ejecutable de pw.x (anula la configuración)")
        p.add_argument("--nproc", type=int, help="procesos MPI por cálculo")
        p.add_argument("-j", "--jobs", type=int, default=1, metavar="N",
                       help="cálculos simultáneos (default: 1). Sin --nproc, "
                            "los hilos de la máquina se reparten entre ellos")
        p.add_argument("--redo", action="store_true",
                       help="rehacer también los cálculos que ya estaban "
                            "terminados")
        p.add_argument("--max-time", dest="max_time", metavar="T",
                       help="presupuesto TOTAL de tiempo: 90m, 2h, 3600. Al "
                            "agotarse no se lanzan más y el barrido queda "
                            "reanudable")
        p.add_argument("--estimate", action="store_true",
                       help="estimar cuánto va a tardar el barrido y salir, "
                            "usando el histórico de 'olla-dft db'")
        p.add_argument("--timeout", type=float,
                       help="límite en segundos por cálculo")
        p = parser_.add_argument_group(GRUPO_DFT)
        p.add_argument("--pseudo-dir", help="carpeta de pseudopotenciales")
        p.add_argument("--pseudo", action="append", metavar="EL=UPF",
                       help="forzar un pseudopotencial concreto, por ejemplo "
                            "Fe=Fe.rel-pbe.UPF. Se puede repetir. Sin esto, "
                            "Olla-DFT elige con 'olla-dft pseudos'")
        p.add_argument("--ecutwfc", type=float, help="cutoff de ondas (Ry)")
        p.add_argument("--ecutrho", type=float, help="cutoff de densidad (Ry)")
        p.add_argument("--kspacing", type=float, help="espaciado k en Å^-1")
        p.add_argument("--insulator", action="store_true",
                       help="occupations='fixed'")
        p = parser_.add_argument_group(GRUPO_FIGURA)
        p.add_argument("--dpi", type=int, default=600)
        p.add_argument("--format", default="pdf,png")
        p.add_argument("--no-plot", action="store_true")
        p.add_argument("-t", "--template", default=None,
                       help="plantilla visual de la figura")
        p.add_argument("--size", default=None, choices=sorted(qstyle.STYLES))
        p.add_argument("--font", default=None, choices=["sans", "serif", "latex"])
        p.add_argument("--usetex", action="store_true")
        p.add_argument("--palette", default=None)
        p.add_argument("--background", default=None)
        p.add_argument("--journal", default="generic",
                       choices=sorted(qstyle.JOURNALS))
        p.add_argument("--width", default=None)
        p.add_argument("--aspect", type=float, default=None)
        p.add_argument("--mono", action="store_true")

    p = sub.add_parser("converge",
                       help="pruebas de convergencia de cutoffs y malla k")
    _calc_opts(p, "convergencia")
    p.add_argument("-k", "--kind", default="ecutwfc",
                   choices=["ecutwfc", "ecutrho", "kmesh"],
                   help="qué parámetro se barre (default: ecutwfc)")
    p.add_argument("--values",
                   help="valores separados por coma; para kmesh admite 8x8x8")
    p.add_argument("--threshold", type=float, default=1.0,
                   help="umbral de convergencia en meV/átomo (default: 1)")

    p = sub.add_parser("eos", help="ecuación de estado E–V y módulo de bulk")
    _calc_opts(p, "eos")
    p.add_argument("--npoints", type=int, default=9,
                   help="número de volúmenes (default: 9)")
    p.add_argument("--scale", type=float, default=1.0,
                   help="factor lineal en que centrar el barrido "
                        "(lo devuelve 'olla-dft mlip scan')")
    p.add_argument("--span", type=float, default=0.10,
                   help="variación relativa de volumen a cada lado (default: 0.10)")
    p.add_argument("--equation", default="birch-murnaghan",
                   choices=["birch-murnaghan", "murnaghan", "vinet"],
                   help="ecuación que se grafica")
    p.add_argument("--relax-ions", action="store_true",
                   help="relajar posiciones internas en cada volumen")

    p = sub.add_parser("elastic",
                       help="constantes elásticas y propiedades mecánicas")
    _calc_opts(p, "elastic")
    p.add_argument("--delta", type=float, default=0.010,
                   help="deformación máxima aplicada (default: 0.010 = 1 %%)")
    p.add_argument("--npoints", type=int, default=4,
                   help="deformaciones no nulas por componente, par (default: 4)")
    p.add_argument("--2d", dest="dosd", action="store_true",
                   help="lámina: constantes en N/m (no en GPa), solo ε1, ε2 y "
                        "ε6, y criterios de Born en 2D")
    p.add_argument("--thickness", type=float, default=None, metavar="A",
                   help="espesor supuesto en Å para dar también el "
                        "equivalente en GPa (convenio, no medida)")
    p.add_argument("--ion-mode", default="auto",
                   choices=["auto", "relax", "fixed"],
                   help="posiciones internas: auto = fijas en deformaciones "
                        "normales y relajadas en cizallas (recomendado); "
                        "relax = relajar todas; fixed = clamped-ion")

    p = sub.add_parser("strain",
                       help="barrido de deformación: gap, energía y momento "
                            "en función de la deformación aplicada")
    _calc_opts(p, "strain")
    p.add_argument("-m", "--mode", default="biaxial",
                   choices=sorted(strain_mod.MODOS),
                   help="qué se deforma (default: biaxial)")
    p.add_argument("-r", "--range", dest="range", default="-5:5:11",
                   metavar="MIN:MAX:N",
                   help="rango en POR CIENTO, por ejemplo -5:5:11 "
                        "(de -5 %% a +5 %% en 11 puntos)")
    p.add_argument("--fixed-ions", action="store_true",
                   help="no relajar las posiciones internas en cada "
                        "deformación (más rápido y menos realista)")
    p.add_argument("--relax-perp", action="store_true",
                   help="dejar libre el eje perpendicular al plano deformado "
                        "(relajación de Poisson); imprescindible en láminas")
    p.add_argument("--nspin", type=int, default=1, choices=[1, 2],
                   help="2 activa la polarización de espín")
    p.add_argument("--mag", help="magnetización inicial (implica --nspin 2)")
    p.add_argument("--hubbard", metavar="EL=U", action="append",
                   help="U de Hubbard en eV por elemento")
    p.add_argument("--vdw", default=None, choices=list(inputgen.VDW),
                   help="corrección de dispersión")

    p = sub.add_parser("adsorb",
                       help="sitios de adsorción sobre una losa y su energía")
    _calc_opts(p, "adsorb")
    p.add_argument("--mol", required=True, metavar="MOLECULA",
                   help="adsorbato: nombre de la base de ASE (CO2, H2O, CO, "
                        "NH3, O2...) o un archivo con la molécula")
    p.add_argument("--sites", default="top,bridge,hollow",
                   help="tipos de sitio a probar (default: los tres)")
    p.add_argument("--height", type=float, default=2.0,
                   help="altura inicial del adsorbato sobre el sitio, en Å "
                        "(default: 2.0)")
    p.add_argument("--face", default="top", choices=["top", "bottom"],
                   help="cara de la losa donde adsorber")
    p.add_argument("--rotations", type=int, default=1,
                   help="orientaciones a probar girando alrededor de la "
                        "normal (default: 1)")
    p.add_argument("--anchor", type=int, default=0,
                   help="átomo de la molécula que se apoya en el sitio "
                        "(índice desde 0; default: 0)")
    p.add_argument("--fixed-ions", action="store_true",
                   help="no relajar: solo scf en la geometría inicial")
    p.add_argument("--vdw", default=None, choices=list(inputgen.VDW),
                   help="corrección de dispersión (casi obligatoria en "
                        "fisisorción)")
    p.add_argument("--dipole", action="store_true",
                   help="corrección dipolar: la losa con adsorbato en una "
                        "sola cara es polar")
    p.add_argument("--nspin", type=int, default=1, choices=[1, 2])
    p.add_argument("--mag", help="magnetización inicial (implica --nspin 2)")

    p = sub.add_parser("eform",
                       help="energía de formación de defectos cargados, "
                            "niveles de transición y diagrama E_f vs ε_F")
    _calc_opts(p, "formacion")
    p.add_argument("-k", "--kind", default="vacancy",
                   choices=["vacancy", "substitution", "interstitial"],
                   help="tipo de defecto")
    p.add_argument("--site", type=int, default=0,
                   help="índice del átomo afectado en la supercelda (base 0)")
    p.add_argument("--new-element", help="especie que entra")
    p.add_argument("--position", help="x,y,z fraccionarias (intersticial)")
    p.add_argument("--supercell", default="2x2x2",
                   help="tamaño de la supercelda (default: 2x2x2)")
    p.add_argument("-q", "--charges", default="0",
                   help="estados de carga separados por coma, por ejemplo "
                        "-2,-1,0,1,2")
    p.add_argument("--epsilon", type=float, default=None,
                   help="constante dieléctrica del material, para apantallar "
                        "la corrección de imagen")
    p.add_argument("--correction", default="lany-zunger",
                   choices=list(defects_mod.ESQUEMAS),
                   help="esquema de corrección de tamaño finito")
    p.add_argument("--mu", action="append", metavar="EL=eV",
                   help="potencial químico por elemento, en eV por átomo. "
                        "Se puede repetir")
    p.add_argument("--align", nargs=2, metavar=("POT_DEF", "POT_PERF"),
                   help="dos archivos cube de potencial electrostático "
                        "(defecto y perfecto) para el término ΔV")
    p.add_argument("--dv", type=float, default=None,
                   help="alineamiento ΔV en eV, si ya lo tienes calculado")
    p.add_argument("--fixed-ions", action="store_true",
                   help="no relajar el defecto en cada estado de carga")
    p.add_argument("--vdw", default=None, choices=list(inputgen.VDW))
    p.add_argument("--nspin", type=int, default=1, choices=[1, 2])
    p.add_argument("--mag", help="magnetización inicial (implica --nspin 2)")

    p = sub.add_parser("gamma",
                       help="energía de superficie y de escisión por el "
                            "ajuste lineal de Fiorentini–Methfessel, con la "
                            "convergencia contra el grosor de la losa")
    _calc_opts(p, "gamma")
    p.add_argument("-m", "--miller", default="1 0 0",
                   help="índices de Miller de la cara, por ejemplo '1 1 1'")
    p.add_argument("-l", "--layers", default="3,4,5,6",
                   help="grosores a calcular, separados por coma "
                        "(default: 3,4,5,6). Hacen falta al menos dos")
    p.add_argument("--vacuum", type=float, default=20.0,
                   help="vacío en Å (default: 20)")
    p.add_argument("--fix", type=int, default=0, metavar="N",
                   help="congelar N capas del fondo al relajar")
    p.add_argument("--relax", action="store_true",
                   help="relajar las posiciones (γ baja entre un 5 y un 20 %%)")
    p.add_argument("--no-bulk", action="store_true",
                   help="no calcular el bulto aparte; solo el ajuste lineal "
                        "E_losa(N) = 2γA + N·E_bulto")
    p.add_argument("--no-reduce", action="store_true",
                   help="no reducir la celda superficial a la mínima (por "
                        "omisión sí se reduce: mismo γ, mucho más barato)")
    p.add_argument("--vdw", default=None, choices=list(inputgen.VDW))
    p.add_argument("--dipole", action="store_true",
                   help="corrección dipolar, para losas polares")
    p.add_argument("--nspin", type=int, default=1, choices=[1, 2])
    p.add_argument("--mag", help="magnetización inicial (implica --nspin 2)")

    # ---------------- materiales laminares ----------------
    p = sub.add_parser("layers",
                       help="detectar capas, espaciado basal y hueco interlaminar")
    p.add_argument("file")
    p.add_argument("--tol", type=float, default=0.45,
                   help="tolerancia de enlace sobre radios covalentes (Å)")
    p.add_argument("--wavelength", default="CuKa",
                   help="radiación para las reflexiones basales (default CuKa)")
    p.add_argument("--slab", metavar="ARCHIVO",
                   help="además, escribir la monocapa con vacío a este archivo")
    p.add_argument("--vacuum", type=float, default=20.0,
                   help="vacío de la monocapa en Å (default 20)")

    p = sub.add_parser("xrd", help="difractograma de polvos simulado")
    p.add_argument("file")
    p.add_argument("-o", "--outdir", default=".")
    p.add_argument("--suite", action="store_true",
                   help="además, exportar JSON de intercambio para las "
                        "otras apps de la suite")
    p.add_argument("--basis", default="conventional",
                   choices=["conventional", "input"],
                   help="celda en que se indexan los hkl: 'conventional' "
                        "(por defecto, los índices de las fichas PDF) o "
                        "'input' (la celda del archivo tal cual)")
    p.add_argument("--wavelength", default="CuKa",
                   help=f"radiación: {', '.join(sorted(WAVE_CHOICES))} o λ en Å")
    p.add_argument("--tt-min", type=float, default=5.0, help="2θ mínimo (°)")
    p.add_argument("--tt-max", type=float, default=70.0, help="2θ máximo (°)")
    p.add_argument("--fwhm", type=float, default=0.15,
                   help="anchura instrumental (° 2θ, default 0.15)")
    p.add_argument("--size", type=float, default=None,
                   help="tamaño de cristalito en nm (anchura por Scherrer)")
    p.add_argument("--biso", type=float, default=0.0,
                   help="factor de temperatura global B (Å²)")
    p.add_argument("--exp", help="difractograma experimental (2θ, I) para superponer")
    p.add_argument("--dpi", type=int, default=600)
    p.add_argument("--format", default="pdf,png")
    p.add_argument("--no-plot", action="store_true")
    p.add_argument("-t", "--template", default=None)
    p.add_argument("--size-preset", default=None, choices=sorted(qstyle.STYLES),
                   help="escala tipográfica de la figura")
    p.add_argument("--font", default=None, choices=["sans", "serif", "latex"])
    p.add_argument("--usetex", action="store_true")
    p.add_argument("--palette", default=None)
    p.add_argument("--background", default=None)
    p.add_argument("--journal", default="generic", choices=sorted(qstyle.JOURNALS))
    p.add_argument("--width", default=None)
    p.add_argument("--aspect", type=float, default=None)
    p.add_argument("--mono", action="store_true")

    p = sub.add_parser("exfoliate",
                       help="energía de exfoliación (bulk vs monocapa)")
    _calc_opts(p, "exfoliacion")
    p.add_argument("--vacuum", type=float, default=20.0,
                   help="vacío de la monocapa en Å (default 20)")
    p.add_argument("--vdw", default=None, choices=list(inputgen.VDW),
                   help="corrección de dispersión para ambos cálculos")
    p.add_argument("--tol", type=float, default=0.45,
                   help="tolerancia de enlace para detectar las capas (Å)")
    p.add_argument("--relax-slab", action="store_true",
                   help="relajar las posiciones de la monocapa")

    # ---------------- campos, ópticas, fonones ----------------
    def _fig_opts_min(p):
        p = p.add_argument_group(GRUPO_FIGURA)
        p.add_argument("--dpi", type=int, default=600)
        p.add_argument("--format", default="pdf,png")
        p.add_argument("--no-plot", action="store_true")
        p.add_argument("-t", "--template", default=None)
        p.add_argument("--font", default=None, choices=["sans", "serif", "latex"])
        p.add_argument("--usetex", action="store_true")
        p.add_argument("--palette", default=None)
        p.add_argument("--background", default=None)
        p.add_argument("--journal", default="generic",
                       choices=sorted(qstyle.JOURNALS))
        p.add_argument("--width", default=None)
        p.add_argument("--mono", action="store_true")

    p = sub.add_parser("wf", help="función trabajo desde un cálculo con vacío")
    p.add_argument("path", nargs="?", default=".", help="carpeta del cálculo")
    p.add_argument("-o", "--outdir", default=".")
    p.add_argument("--axis", default="c", help="eje del vacío: a/b/c (default c)")
    p.add_argument("--rerun", action="store_true",
                   help="volver a correr pp.x aunque exista el cube")
    p.add_argument("--pw-cmd", dest="pw_cmd")
    p.add_argument("--nproc", type=int)
    _fig_opts_min(p)

    p = sub.add_parser("align",
                       help="alineamiento de bandas entre dos materiales: "
                            "offsets ΔE_v, ΔE_c y tipo I/II/III")
    p.add_argument("a", help="carpeta del cálculo del primer material")
    p.add_argument("b", help="carpeta del cálculo del segundo material")
    p.add_argument("--interface", metavar="CARPETA",
                   help="carpeta de la interfaz; activa el método riguroso "
                        "de Van de Walle-Martin")
    p.add_argument("--names", help="nombres para el reporte, separados por "
                                   "coma (por omisión, los de las carpetas)")
    p.add_argument("--axis", default="c", help="eje del perfil planar")
    p.add_argument("--window", type=float, default=None, metavar="A",
                   help="ventana del promedio macroscópico en Å "
                        "(por omisión, un octavo de la celda)")
    p.add_argument("--rerun", action="store_true",
                   help="volver a correr pp.x aunque ya exista el cube")
    p.add_argument("-o", "--outdir", default="alineamiento")
    p.add_argument("--pw-cmd", dest="pw_cmd")
    p.add_argument("--nproc", type=int)
    _fig_opts_min(p)

    p = sub.add_parser("charge",
                       help="densidad de carga / ELF / espín con pp.x")
    p.add_argument("path", nargs="?", default=".", help="carpeta del cálculo")
    p.add_argument("-o", "--outdir", default=".")
    p.add_argument("--field", default="density",
                   choices=["density", "elf", "spin", "potential", "vtotal"])
    p.add_argument("--axis", default="c", help="eje del perfil planar")
    p.add_argument("--rerun", action="store_true")
    p.add_argument("--pw-cmd", dest="pw_cmd")
    p.add_argument("--nproc", type=int)
    _fig_opts_min(p)

    p = sub.add_parser("optics",
                       help="ε(ω), absorción y Tauc con epsilon.x (pseudos NC)")
    _calc_opts(p, "opticas")
    p.add_argument("--wmax", type=float, default=20.0,
                   help="energía máxima del espectro (eV)")
    p.add_argument("--smear", type=float, default=0.10,
                   help="ensanchamiento interbanda (eV)")
    p.add_argument("--metal", action="store_true",
                   help="sistema metálico (ocupaciones con smearing)")
    p.add_argument("--suite", action="store_true",
                   help="además, exportar JSON de intercambio para las "
                        "otras apps de la suite")
    p.add_argument("--tauc", default="direct", choices=["direct", "indirect"],
                   help="tipo de transición para la gráfica de Tauc")
    p.add_argument("--scissor", type=float, default=0.0,
                   help="corrimiento rígido del gap en eV (gap experimental "
                        "o GW menos el gap del cálculo); desplaza ε2 y "
                        "rehace ε1 por Kramers-Kronig")

    p = sub.add_parser("effmass",
                       help="masa efectiva por ajuste parabólico de bandas")
    p.add_argument("file")
    p.add_argument("-o", "--outdir", default="masa_efectiva")
    p.add_argument("--bands-dir",
                   help="carpeta con un cálculo de bandas ya hecho "
                        "(de ahí salen VBM y CBM)")
    p.add_argument("--collect", action="store_true",
                   help="leer el cálculo fino ya corrido")
    p.add_argument("--run", action="store_true",
                   help="correr el cálculo fino en cuanto se prepare")
    p.add_argument("--half-width", type=float, default=0.06,
                   help="semiancho de cada línea en Å⁻¹")
    p.add_argument("--points", type=int, default=21,
                   help="puntos k por línea (impar)")
    # sin valor aquí: lo pone effmass.WINDOW_DEFAULT (la mitad del límite
    # parabólico del módulo), para que el ajuste por omisión no salga
    # siempre con el aviso de "fuera del régimen parabólico"
    p.add_argument("--window", type=float, default=None,
                   help="semiancho del ajuste rápido sobre el camino, en "
                        "Å⁻¹ a cada lado del extremo (por omisión, la "
                        "mitad del límite parabólico: ±0.06)")
    p.add_argument("--min-points", type=int, default=7,
                   help="puntos mínimos del ajuste rápido")
    p.add_argument("--pseudo-dir")
    p.add_argument("--pseudo", action="append", metavar="EL=UPF",
                   help="forzar un pseudopotencial concreto, por ejemplo Fe=Fe.rel-pbe.UPF. Se puede repetir. Sin esto, Olla-DFT elige con 'olla-dft pseudos'")
    p.add_argument("--ecutwfc", type=float)
    p.add_argument("--ecutrho", type=float)
    p.add_argument("--pw-cmd")
    p.add_argument("--nproc", type=int)
    p.add_argument("--timeout", type=int)

    p = sub.add_parser("surface", help="cortar una superficie (hkl) con vacío")
    p.add_argument("file")
    p.add_argument("-m", "--miller", default="1 0 0",
                   help="índices de Miller, por ejemplo '1 1 1' o 1,1,1")
    p.add_argument("-l", "--layers", type=int, default=6)
    p.add_argument("--vacuum", type=float, default=15.0,
                   help="vacío total en Å")
    p.add_argument("--fix", type=int, default=0,
                   help="planos atómicos del fondo a congelar")
    p.add_argument("-o", "--output", help="archivo de salida (CIF/POSCAR)")

    p = sub.add_parser("defect", help="crear un defecto puntual")
    p.add_argument("file")
    p.add_argument("-k", "--kind", default="vacancy",
                   choices=["vacancy", "substitution", "interstitial"])
    p.add_argument("--site", type=int, default=0,
                   help="índice del átomo afectado (base 0)")
    p.add_argument("--new-element", help="especie que entra")
    p.add_argument("--supercell", help="por ejemplo 3x3x3")
    p.add_argument("--position", help="x,y,z fraccionarias (intersticial)")
    p.add_argument("-o", "--outdir", default="defecto")

    p = sub.add_parser("charges",
                       help="cargas de Löwdin/Bader y diferencia de densidad")
    p.add_argument("file", nargs="?", help="estructura (para Bader)")
    p.add_argument("--lowdin", help="salida de projwfc.x")
    p.add_argument("--bader", help="cube de densidad (plot_num=0)")
    p.add_argument("--difference", nargs="+", metavar="CUBE",
                   help="total.cube parte1.cube parte2.cube ...")
    p.add_argument("--pseudo-dir",
                   help="carpeta con los UPF del cálculo: de ahí sale "
                        "Z_valencia para la columna 'neta' (anula config)")
    p.add_argument("--axis", type=int, default=2, choices=[0, 1, 2])
    p.add_argument("-o", "--outdir", default=".")
    p.add_argument("--dpi", type=int)
    p.add_argument("--format", default="pdf,png")
    p.add_argument("--no-plot", action="store_true")
    p.add_argument("-t", "--template")
    p.add_argument("--font", choices=["sans", "serif", "latex"])
    p.add_argument("--usetex", action="store_true")
    p.add_argument("--palette")
    p.add_argument("--background")
    p.add_argument("--journal", default="generic",
                   choices=sorted(["acs", "aps", "elsevier", "generic",
                                   "iop", "nature", "rsc", "wiley"]))
    p.add_argument("--width")
    p.add_argument("--mono", action="store_true")

    p = sub.add_parser("fermi",
                       help="exportar la superficie de Fermi en BXSF")
    p.add_argument("-o", "--outdir", default="transporte")

    p = sub.add_parser("xps",
                       help="corrimientos de nivel de core (estado inicial)")
    p.add_argument("file")
    p.add_argument("-o", "--outdir", default="xps")
    p.add_argument("--core-hole", action="append", metavar="EL=UPF",
                   help="pseudopotencial con hueco de core, por ejemplo "
                        "Si=Si.star1s.UPF. Se puede repetir. Sin esto, "
                        "initial_state.x devuelve una tabla de ceros")
    p.add_argument("--collect", action="store_true")
    p.add_argument("--suite", action="store_true")
    p.add_argument("--pseudo-dir")
    p.add_argument("--pseudo", action="append", metavar="EL=UPF",
                   help="forzar un pseudopotencial concreto, por ejemplo Fe=Fe.rel-pbe.UPF. Se puede repetir. Sin esto, Olla-DFT elige con 'olla-dft pseudos'")
    p.add_argument("--ecutwfc", type=float)
    p.add_argument("--ecutrho", type=float)
    p.add_argument("--kspacing", type=float)
    p.add_argument("--metal", action="store_true")

    p = sub.add_parser("corehole",
                       help="generar el par de pseudopotenciales normal + "
                            "hueco de core (ld1.x), para XPS y XANES")
    p.add_argument("element", nargs="?", default=None,
                   help="simbolo del elemento, por ejemplo Si")
    p.add_argument("-o", "--outdir", default="pseudos")
    p.add_argument("--edge", default="K",
                   help="borde/nivel del hueco: K (1s), L1 (2s), L23 (2p), "
                        "M1, M23, M45")
    p.add_argument("--functional", default="PBE",
                   help="funcional del pseudopotencial; tiene que ser el "
                        "mismo con el que vas a correr pw.x")
    p.add_argument("--rcut", type=float,
                   help="radio de corte en bohr (por omision, uno por fila "
                        "de la tabla periodica)")
    p.add_argument("--rel", type=int, default=0, choices=(0, 1, 2),
                   help="0 no relativista, 1 escalar, 2 completo")
    p.add_argument("--semicore", action="store_true",
                   help="mete la capa (n-1)s(n-1)p en la valencia")
    p.add_argument("--pseudotype", type=int, default=2, choices=(1, 2, 3),
                   help="1 y 2 conservan la norma, 3 es ultrasuave")
    p.add_argument("--plain", action="store_true",
                   help="generar SOLO el pseudopotencial normal, sin el de "
                        "hueco de core. Sirve para tener un pseudo "
                        "consistente de un elemento que no lo admita")
    p.add_argument("--only-inputs", action="store_true",
                   help="escribir los inputs de ld1.x sin ejecutarlos")
    p.add_argument("--projectors", type=int, default=1, choices=(1, 2),
                   help="proyectores GIPAW por canal. XSpectra recomienda 2, "
                        "pero con 2 el pseudo sale ultrasuave y suele haber "
                        "que ajustar --rcut a mano")
    p.add_argument("--ld1-cmd", help="ruta a ld1.x")
    p.add_argument("--core-wfc", metavar="UPF",
                   help="en vez de generar: extraer de un UPF la funcion de "
                        "onda de core en el formato que lee xspectra.x")
    p.add_argument("--orbital", help="orbital a verificar, por ejemplo 1S")
    p.add_argument("--output", help="archivo de salida para --core-wfc")

    p = sub.add_parser("xanes",
                       help="XANES/NEXAFS: absorcion de rayos X cerca del "
                            "borde (xspectra.x)")
    p.add_argument("file", help="estructura")
    p.add_argument("-o", "--outdir", default="xanes")
    p.add_argument("--element", help="elemento que absorbe")
    p.add_argument("--site", type=int, default=0,
                   help="cual atomo de ese elemento (desde 0)")
    p.add_argument("--edge", default="K",
                   help="borde: K, L1, L2, L3 o L23 (los que calcula "
                        "xspectra.x; los bordes M no)")
    p.add_argument("--core-hole", metavar="UPF",
                   help="pseudopotencial con hueco de core (olla-dft corehole)")
    p.add_argument("--polarization", default="1 0 0",
                   help="direccion del campo electrico, por ejemplo '0 0 1'")
    p.add_argument("--average", action="store_true",
                   help="tres direcciones ortogonales y promedio: es lo que "
                        "corresponde a una muestra en polvo")
    p.add_argument("--emin", type=float, default=-10.0)
    p.add_argument("--emax", type=float, default=30.0)
    p.add_argument("--broadening", type=float, default=0.8,
                   help="ensanchamiento en eV (xgamma)")
    p.add_argument("--r-paw", type=float, default=3.0)
    p.add_argument("--collect", action="store_true")
    p.add_argument("--pseudo-dir")
    p.add_argument("--pseudo", action="append", metavar="EL=UPF",
                   help="forzar un pseudopotencial concreto, por ejemplo Fe=Fe.rel-pbe.UPF. Se puede repetir. Sin esto, Olla-DFT elige con 'olla-dft pseudos'")
    p.add_argument("--ecutwfc", type=float)
    p.add_argument("--ecutrho", type=float)
    p.add_argument("--kspacing", type=float)
    p.add_argument("--metal", action="store_true")
    _fig_opts_min(p)

    p = sub.add_parser("hubbard",
                       help="U de Hubbard por respuesta lineal (hp.x), en vez "
                            "de copiarlo de un articulo")
    p.add_argument("file", help="estructura")
    p.add_argument("-o", "--outdir", default="hubbard")
    p.add_argument("--species",
                   help="especies a perturbar, separadas por coma. Por "
                        "omision, los metales de transicion y tierras raras "
                        "de la estructura")
    p.add_argument("--qgrid", default="2x2x2",
                   help="malla de q de la respuesta lineal; equivale a una "
                        "supercelda de nq1*nq2*nq3 celdas")
    p.add_argument("--projection", default="ortho-atomic",
                   choices=["atomic", "ortho-atomic", "norm-atomic",
                            "wannier", "pseudo"],
                   help="esquema de proyeccion. El U SOLO vale con el mismo "
                        "esquema con el que se calculo")
    p.add_argument("--hubbard-style", dest="hubbard_style", default="legacy",
                   choices=["legacy", "card"],
                   help="sintaxis de DFT+U del scf: legacy = lda_plus_u "
                        "(QE <= 7.0), card = tarjeta HUBBARD (QE >= 7.1, "
                        "donde la sintaxis vieja es un error)")
    p.add_argument("--cycle", action="store_true",
                   help="ciclo de autoconsistencia completo: scf -> hp.x -> "
                        "scf con el U nuevo, hasta que deje de moverse")
    p.add_argument("--max-iter", type=int, default=6)
    p.add_argument("--tol", type=float, default=0.05,
                   help="cambio en eV por debajo del cual se da por "
                        "convergido")
    p.add_argument("--mixing", type=float, default=1.0,
                   help="amortiguacion del paso; bajalo a 0.5 si oscila")
    p.add_argument("--collect", action="store_true")
    p.add_argument("--pw-cmd", dest="pw_cmd")
    p.add_argument("--nproc", type=int)
    p.add_argument("--pseudo-dir")
    p.add_argument("--pseudo", action="append", metavar="EL=UPF",
                   help="forzar un pseudopotencial concreto, por ejemplo Fe=Fe.rel-pbe.UPF. Se puede repetir. Sin esto, Olla-DFT elige con 'olla-dft pseudos'")
    p.add_argument("--ecutwfc", type=float)
    p.add_argument("--ecutrho", type=float)
    p.add_argument("--kspacing", type=float)
    p.add_argument("--metal", action="store_true")
    p.add_argument("--nspin", type=int, default=1, choices=(1, 2))
    p.add_argument("--mag")
    p.add_argument("--intersite", action="store_true",
                   help="además de las U, leer los V intersitio que hp.x ya "
                        "escribe y generar la tarjeta HUBBARD de QE >= 7.1")
    p.add_argument("--v-threshold", dest="v_threshold", type=float,
                   default=0.01, metavar="eV",
                   help="V por debajo de esto no se lista ni se escribe")

    p = sub.add_parser("interface",
                       help="heteroestructura: apilar dos materiales con la "
                            "menor deformacion de red posible")
    p.add_argument("file1", help="material de abajo (el sustrato)")
    p.add_argument("file2", help="material de arriba")
    p.add_argument("-o", "--outdir", default=".")
    p.add_argument("--name", default="heteroestructura")
    p.add_argument("--max-index", type=int, default=4,
                   help="mayor coeficiente entero de la supercelda; subirlo "
                        "encuentra celdas mas giradas pero tarda mucho mas")
    p.add_argument("--tol", type=float, default=0.05,
                   help="deformacion maxima aceptada (0.05 = 5 %%)")
    p.add_argument("--max-atoms", type=int, default=200)
    p.add_argument("--index", type=int, default=0,
                   help="cual de las candidatas construir")
    p.add_argument("--top", type=int, default=10)
    p.add_argument("--list", action="store_true",
                   help="solo listar las candidatas, sin construir nada")
    p.add_argument("--separation", type=float,
                   help="distancia inicial entre capas en A; por omision, de "
                        "los radios de van der Waals")
    p.add_argument("--vacuum", type=float, default=20.0)
    p.add_argument("--strain", default="second",
                   choices=["first", "second", "both"],
                   help="quien se deforma: el de abajo, el de arriba, o los "
                        "dos a medias")
    p.add_argument("--shift", help="desplazamiento lateral del material de "
                                   "arriba, en fracciones de la celda comun")

    p = sub.add_parser("md",
                       help="analizar una trayectoria de dinamica molecular: "
                            "g(r), difusion y espectro vibracional")
    p.add_argument("path", help="salida de pw.x con calculation='md', o su "
                                "carpeta")
    p.add_argument("-o", "--outdir", default=".")
    p.add_argument("--skip", type=int, default=0,
                   help="pasos iniciales a descartar (equilibrado)")
    p.add_argument("--rmax", type=float,
                   help="corte de g(r) en A; por omision, media arista de la "
                        "celda, que es hasta donde la normalizacion vale")
    p.add_argument("--bins", type=int, default=200)
    _fig_opts_min(p)

    p = sub.add_parser("neb",
                       help="camino de reaccion y barrera de activacion "
                            "(neb.x)")
    p.add_argument("file", help="estructura inicial (reactivo)")
    p.add_argument("final", nargs="?", help="estructura final (producto)")
    p.add_argument("-o", "--outdir", default="neb")
    p.add_argument("--images", type=int, default=7,
                   help="numero de imagenes de la cadena")
    p.add_argument("--no-ci", action="store_true",
                   help="sin imagen trepadora; la barrera saldra "
                        "SUBESTIMADA")
    p.add_argument("--path-thr", type=float, default=0.05,
                   help="umbral de fuerza del camino en eV/A")
    p.add_argument("--nstep", type=int, default=50)
    p.add_argument("--fix", help="indices de atomos a congelar (base 0)")
    p.add_argument("--prefix")
    p.add_argument("--collect", action="store_true")
    p.add_argument("--pseudo-dir")
    p.add_argument("--pseudo", action="append", metavar="EL=UPF",
                   help="forzar un pseudopotencial concreto, por ejemplo Fe=Fe.rel-pbe.UPF. Se puede repetir. Sin esto, Olla-DFT elige con 'olla-dft pseudos'")
    p.add_argument("--ecutwfc", type=float)
    p.add_argument("--ecutrho", type=float)
    p.add_argument("--kspacing", type=float)
    p.add_argument("--metal", action="store_true")
    p.add_argument("--nspin", type=int, default=1, choices=(1, 2))
    p.add_argument("--mag")
    _fig_opts_min(p)

    p = sub.add_parser("thermochem",
                       help="ZPE, entropia y energia libre: de una energia "
                            "DFT a una comparable con el experimento")
    p.add_argument("freqs",
                   help="archivo de frecuencias en cm-1, o la lista separada "
                        "por comas")
    p.add_argument("--phase", default="solido",
                   choices=["solido", "adsorbato", "gas", "transicion"],
                   help="gas anade traslaciones y rotaciones; transicion "
                        "exige exactamente una frecuencia imaginaria")
    p.add_argument("--structure", help="estructura (necesaria para la fase gas)")
    p.add_argument("--temp", type=float, default=298.15)
    p.add_argument("--pressure", type=float, default=1.0, help="en bar")
    p.add_argument("--symmetry", type=int, default=1,
                   help="numero de simetria del grupo puntual: 2 para H2O y "
                        "O2, 3 para NH3, 12 para CH4")
    p.add_argument("--multiplicity", type=int, default=1,
                   help="multiplicidad de espin del estado fundamental")
    p.add_argument("--floor", type=float,
                   help="sube los modos por debajo de este valor (cm-1); "
                        "100 es lo habitual")
    p.add_argument("--energy", type=float, help="E_DFT en eV, para dar G(T)")
    p.add_argument("-o", "--outdir")

    p = sub.add_parser("unfold",
                       help="desdoblar las bandas de una supercelda sobre la "
                            "zona de Brillouin primitiva")
    p.add_argument("path", help="carpeta del calculo de bandas de la "
                                "supercelda")
    p.add_argument("primitive", help="estructura de la celda PRIMITIVA")
    p.add_argument("-o", "--outdir", default=".")
    p.add_argument("--prefix")
    p.add_argument("--bands", type=int,
                   help="cuantas bandas desdoblar (desde la mas baja)")
    p.add_argument("--spin", choices=["up", "dw"], default="up",
                   help="canal de espin a desdoblar si el calculo es lsda "
                        "(se desdobla UN canal por corrida; por omision, up)")
    p.add_argument("--emin", type=float, default=-6.0)
    p.add_argument("--emax", type=float, default=6.0)
    _fig_opts_min(p)

    p = sub.add_parser("elph",
                       help="acoplamiento electron-fonon: lambda, Tc y un tau "
                            "de verdad para el transporte")
    p.add_argument("file", nargs="?", help="estructura")
    p.add_argument("-o", "--outdir", default="elph")
    p.add_argument("--qgrid", default="2x2x2")
    p.add_argument("--kgrid", help="malla de k del scf")
    p.add_argument("--kgrid-nscf", dest="kgrid_nscf",
                   help="malla de k del nscf denso; por omision, el doble de "
                        "la del scf redondeada a un multiplo de la de q")
    p.add_argument("--nsigma", type=int, default=10)
    p.add_argument("--sigma", type=float, default=0.005,
                   help="paso del barrido de ensanchamiento, en Ry")
    p.add_argument("--degauss", type=float, default=0.02)
    p.add_argument("--debye", type=float,
                   help="temperatura de Debye en K, para marcar el regimen "
                        "en el que vale la formula de tau")
    p.add_argument("--collect", action="store_true")
    p.add_argument("--pseudo-dir")
    p.add_argument("--pseudo", action="append", metavar="EL=UPF",
                   help="forzar un pseudopotencial concreto, por ejemplo Fe=Fe.rel-pbe.UPF. Se puede repetir. Sin esto, Olla-DFT elige con 'olla-dft pseudos'")
    p.add_argument("--ecutwfc", type=float)
    p.add_argument("--ecutrho", type=float)
    _fig_opts_min(p)

    p = sub.add_parser("teoria", aliases=["theory"],
                       help="el fundamento físico de un comando: qué responde, "
                            "las fórmulas que implementa, de qué módulo salen "
                            "y de dónde sale cada dato")
    p.add_argument("comando", nargs="?",
                   help="comando a explicar; sin nada, el índice")
    p.add_argument("--all", action="store_true",
                   help="el documento completo (todas las áreas)")
    p.add_argument("-o", "--output", metavar="ARCHIVO.md",
                   help="guardarlo en Markdown en vez de imprimirlo")

    p = sub.add_parser("sistema", aliases=["system"],
                       help="qué ve Olla-DFT de esta máquina: codificación, "
                            "dónde guarda la configuración, qué binarios de "
                            "QE encuentra y cómo lanzar los cálculos aquí")

    p = sub.add_parser("update", aliases=["actualizar"],
                       help="comprobar si hay una versión nueva de Olla-DFT y, "
                            "si la hay, instalarla con una confirmación; nunca "
                            "se ejecuta solo")
    p.add_argument("--check", action="store_true",
                   help="solo comprobar e informar, sin instalar nada")
    p.add_argument("--yes", action="store_true",
                   help="no preguntar; instalar directamente si hay versión nueva")
    p.add_argument("--version", metavar="TAG",
                   help="instalar una versión concreta (p. ej. v1.0.1) en vez de la última")

    p = sub.add_parser("start",
                       help="inicio guiado para crear un proyecto sin conocer la CLI")
    p.add_argument("--project", default=".", help="carpeta del proyecto")
    p.add_argument("--structure", help="CIF, POSCAR o input de pw.x")
    p.add_argument("--goal", help="relax, gap, dos, phonons, optics o scf")
    p.add_argument("--name", help="nombre visible del proyecto")
    p.add_argument("--non-interactive", action="store_true",
                   help="no preguntar; requiere --structure en un proyecto nuevo")
    p.add_argument("--no-validate", action="store_true",
                   help="no ejecutar la validación inicial")
    p.add_argument("--language", choices=["es", "en"], default=argparse.SUPPRESS,
                   help="idioma del inicio guiado (default: es)")

    p = sub.add_parser("recetas", aliases=["recipes"],
                       help="sesiones completas de principio a fin: qué "
                            "comando va después de cuál y qué archivo se "
                            "pasan entre ellos")
    p.add_argument("receta", nargs="?",
                   help="clave de la receta; sin nada, las lista todas")
    p.add_argument("--buscar", default=None, metavar="TEXTO",
                   help="buscarla con tus palabras, sin saber la clave")
    p.add_argument("--script", nargs="?", const="", default=None,
                   metavar="ARCHIVO",
                   help="escribir la receta como un guion de shell "
                        "comentado, listo para editar")

    p = sub.add_parser("wizard",
                       help="asistente: dime QUE quieres saber y te digo que "
                            "hay que correr, en orden y con los comandos")
    p.add_argument("file", nargs="?", help="tu estructura (opcional)")
    p.add_argument("--goal", help="clave de la meta; salen con --list")
    p.add_argument("--ask", metavar="TEXTO",
                   help="describelo con tus palabras, por ejemplo "
                        "'quiero saber si absorbe luz'")
    p.add_argument("--list", action="store_true",
                   help="listar todo lo que el asistente sabe hacer")
    p.add_argument("--term", help="que significa un termino")
    p.add_argument("--no-glossary", action="store_true")
    p.add_argument("--pseudo-dir")
    p.add_argument("--pseudo", action="append", metavar="EL=UPF",
                   help="forzar un pseudopotencial concreto, por ejemplo Fe=Fe.rel-pbe.UPF. Se puede repetir. Sin esto, Olla-DFT elige con 'olla-dft pseudos'")

    p = sub.add_parser("resilient", help="cálculos QE recuperables ante interrupciones del servidor")
    p.add_argument("action", choices=["init", "run", "status", "pause", "service"])
    p.add_argument("target", help="input para init; directorio persistente del trabajo para las demás acciones")
    p.add_argument("--state", help="directorio nuevo del trabajo en un disco persistente conservado")
    p.add_argument("--pw-cmd", default="pw.x", help="comando de QE o MPI con paralelismo fijo")
    p.add_argument("--runtime-id", help="identificador de la imagen inmutable del entorno")
    p.add_argument("--checkpoint-seconds", type=float, default=900)
    p.add_argument("--grace-seconds", type=float, default=300)
    p.add_argument("--max-failures", type=int, default=3)
    p.add_argument("--threads", type=int, default=1)
    p.add_argument("--keep", type=int, default=2, help="guardados íntegros que conservar (mínimo 2)")
    p.add_argument("--max-segments", type=int, default=0, help="detenerse tras este número de segmentos guardados; 0 significa sin límite")
    p.add_argument("--resume", action="store_true", help="retirar una pausa explícita antes de continuar")
    p.add_argument("--user", help="usuario sin privilegios del servicio systemd generado")
    p.add_argument("-o", "--output", help="archivo de servicio generado; se instala por separado")

    p = sub.add_parser(
        "project",
        help="gestionar un proyecto reproducible: fuentes, workflow, calidad y dashboard")
    p.add_argument(
        "action",
        choices=["init", "add", "plan", "show", "status", "validate",
                 "run", "dashboard", "report", "export", "ingest",
                 "environment", "diff", "cancel", "resume"],
        help="acción sobre el proyecto")
    p.add_argument("target", nargs="?",
                   help="directorio, archivo, objetivo, perfil o tarea según acción")
    p.add_argument("--project", default=".",
                   help="proyecto desde el que trabajar (default: .)")
    p.add_argument("--name", help="nombre al inicializar")
    # No usar el dest por defecto ``command``: ese nombre pertenece al
    # subparser global y argparse lo sobrescribiría con None al parsear este
    # subcomando, enviando ``olla-dft project ...`` al menú interactivo.
    p.add_argument("--command", dest="task_commands", action="append",
                   help="tarea olla-dft personalizada; se puede repetir con plan")
    p.add_argument("--execute", action="store_true",
                   help="ejecutar run/submit; por omisión solo simula o escribe")
    p.add_argument("--force", action="store_true",
                   help="en run, ignorar la caché y volver a preparar todas las tareas")
    p.add_argument("--parallel", type=int, default=1,
                   help="en run, tareas independientes simultáneas (default: 1)")
    p.add_argument("--retries", type=int, default=0,
                   help="reintentos por tarea fallida (default: 0)")
    p.add_argument("--timeout", type=float,
                   help="tiempo máximo por intento, en segundos")
    p.add_argument("--cancel-file", help="marker de cancelación cooperativa personalizado")
    p.add_argument("--reason", help="en cancel, motivo opcional")
    p.add_argument("--selftest", action="store_true",
                   help="en validate, ejecutar la validación rápida contra referencias físicas")
    p.add_argument("--advanced", action="store_true",
                   help="en validate, revisar estructura, comandos, unidades y colisiones")
    p.add_argument("-o", "--output", help="salida para dashboard, report o export")
    p.add_argument("--pdf", action="store_true",
                   help="en report, generar un informe PDF autocontenido")
    p.add_argument("--theme", choices=["auto", "light", "dark"], default="auto",
                   help="tema del dashboard (default: auto)")
    p.add_argument("--language", choices=["es", "en"], default=argparse.SUPPRESS,
                   help="idioma del dashboard (default: es)")
    p.add_argument("--both", action="store_true",
                   help="generar dashboard español e inglés en archivos separados")
    p.add_argument("--verify-environment", action="store_true",
                   help="en environment, comprobar el bloqueo guardado")
    p.add_argument("--other", help="en diff, snapshot o proyecto de comparación")
    p.add_argument("--json", action="store_true", help="en diff, imprimir JSON")

    p = sub.add_parser(
        "results",
        help="ingerir, consultar y exportar resultados normalizados del proyecto")
    p.add_argument("action", choices=["ingest", "list", "show", "review", "export", "explore"])
    p.add_argument("target", nargs="?",
                   help="ruta de entrada para ingest, o id para show")
    p.add_argument("extra_paths", nargs="*",
                   help="más carpetas/XML para ingest")
    p.add_argument("--project", default=".")
    p.add_argument("--db", help="SQLite alternativa; por defecto .qekit/results.sqlite3")
    p.add_argument("--tag", help="etiqueta de procedencia para ingest")
    p.add_argument("--formula", help="filtrar por fórmula")
    p.add_argument("--calculation", help="filtrar por tipo de cálculo")
    p.add_argument("--status", choices=["invalid", "not_converged", "parsed_no_energy",
                                         "parsed", "converged"])
    p.add_argument("--review-status", choices=["unreviewed", "accepted", "rejected"],
                   help="en review, estado de la revisión humana")
    p.add_argument("--note", help="en review, nota que acompaña la decisión")
    p.add_argument("--limit", type=int, default=None,
                   help="máximo de registros: list=100, explore=10000")
    p.add_argument("--json", action="store_true", help="en list, imprimir JSON")
    p.add_argument("-o", "--output", help="archivo de salida: export=JSON, explore=HTML interactivo")

    p = sub.add_parser(
        "campaign",
        help="crear matrices reproducibles de tareas parametrizadas")
    p.add_argument("action", choices=["create", "list", "status", "export", "run", "extend"])
    p.add_argument("target", nargs="?", help="nombre o id de campaña")
    p.add_argument("--project", default=".")
    p.add_argument("--command", dest="campaign_command",
                   help="plantilla del comando olla-dft; campos: {eje}, {index}, {id}, {structure}")
    p.add_argument("--axis", action="append", default=[],
                   help="eje nombre=v1,v2; se puede repetir")
    p.add_argument("--goal", help="objetivo científico de la campaña")
    p.add_argument("--convergence-file", help="CONVERGENCIA.dat para tomar una recomendación")
    p.add_argument("--adaptive", action="store_true",
                   help="añadir el siguiente valor recomendado al eje de convergencia")
    p.add_argument("--threshold", type=float, default=None,
                   help="umbral de convergencia al extender (meV/átomo)")
    p.add_argument("--execute", action="store_true",
                   help="en run, ejecutar los puntos seleccionados")
    p.add_argument("--force", action="store_true",
                   help="en run, ignorar la caché de tareas")
    p.add_argument("--parallel", type=int, default=1,
                   help="en run, puntos independientes simultáneos (default: 1)")
    p.add_argument("--retries", type=int, default=0,
                   help="reintentos por punto fallido (default: 0)")
    p.add_argument("--timeout", type=float,
                   help="tiempo máximo por intento, en segundos")
    p.add_argument("--cancel-file", help="marker de cancelación cooperativa personalizado")
    p.add_argument("-o", "--output", help="archivo JSON de export")

    p = sub.add_parser("pseudos",
                       help="comparar los pseudopotenciales disponibles y "
                            "elegir con criterio, no por orden alfabetico")
    p.add_argument("file", nargs="?", help="estructura")
    p.add_argument("--element", help="elementos separados por coma")
    p.add_argument("--task", default="general",
                   help="para que es: general, optics, soc, xanes, hubbard, "
                        "fonones. Cada tarea descarta los que no sirven")
    p.add_argument("--functional",
                   help="exigir un funcional concreto (PBE, PZ, PBEsol...)")
    p.add_argument("--cheap", action="store_true",
                   help="preferir ultrasuave/PAW, que necesitan menos ondas "
                        "planas")
    p.add_argument("--pseudo-dir")
    p.add_argument("--pseudo", action="append", metavar="EL=UPF",
                   help="forzar un pseudopotencial concreto, por ejemplo Fe=Fe.rel-pbe.UPF. Se puede repetir. Sin esto, Olla-DFT elige con 'olla-dft pseudos'")

    p = sub.add_parser("tddft",
                       help="absorcion optica con TDDFPT: deja que el "
                            "electron excitado y su hueco se vean")
    p.add_argument("file", nargs="?", help="estructura")
    p.add_argument("-o", "--outdir", default="tddft")
    p.add_argument("--method", default="lanczos",
                   choices=["lanczos", "davidson"],
                   help="lanczos da el espectro entero; davidson da las "
                        "primeras excitaciones una a una")
    p.add_argument("--iter", type=int, default=500,
                   help="iteraciones de Lanczos: manda la resolucion")
    p.add_argument("--pol", type=int, default=4, choices=(1, 2, 3, 4),
                   help="1/2/3 = xx/yy/zz, 4 = tensor completo")
    p.add_argument("--states", type=int, default=10,
                   help="excitaciones a buscar (davidson)")
    p.add_argument("--emin", type=float, default=0.0)
    p.add_argument("--emax", type=float, default=15.0)
    p.add_argument("--broadening", type=float, default=None,
                   help="ensanchamiento en eV (default 0.05). Con --collect "
                        "fija el umbral de deteccion del exciton; si se "
                        "omite se lee de spectrum.in")
    p.add_argument("--scissor", type=float, default=0.0,
                   help="corrimiento rigido de las bandas vacias en eV "
                        "(solo lanczos): compensa el gap subestimado")
    p.add_argument("--extrapolation", default="osc",
                   choices=["no", "constant", "osc"])
    p.add_argument("--tamm-dancoff", action="store_true",
                   help="aproximacion de Tamm-Dancoff: mas barata, no exacta")
    p.add_argument("--rpa", action="store_true",
                   help="apagar el kernel xc, para ver cuanto aporta")
    p.add_argument("--gamma", action="store_true",
                   help="forzar K_POINTS gamma. Se detecta solo cuando la "
                        "estructura es una molecula")
    p.add_argument("--gap", type=float,
                   help="gap de particulas independientes en eV, para "
                        "detectar si hay exciton ligado")
    p.add_argument("--compare", metavar="OPTICS.dat",
                   help="superponer el espectro de 'olla-dft optics'")
    p.add_argument("--nbnd", type=int)
    p.add_argument("--collect", action="store_true")
    p.add_argument("--pseudo-dir")
    p.add_argument("--pseudo", action="append", metavar="EL=UPF")
    p.add_argument("--ecutwfc", type=float)
    p.add_argument("--ecutrho", type=float)
    p.add_argument("--kspacing", type=float)
    p.add_argument("--metal", action="store_true")
    _fig_opts_min(p)

    p = sub.add_parser("ballistic",
                       help="conductancia balistica de Landauer (pwcond.x), "
                            "para nanocontactos y moleculas entre electrodos")
    p.add_argument("file", nargs="?",
                   help="electrodo: la celda periodica en z")
    p.add_argument("--scatterer",
                   help="region de dispersion (la molecula o el defecto). "
                        "Sin esto solo salen las bandas complejas")
    p.add_argument("-o", "--outdir", default="balistico")
    p.add_argument("--ikind", type=int, choices=(0, 1),
                   help="0 = solo bandas complejas, 1 = conductancia con el "
                        "mismo electrodo a los dos lados (default: 1 si hay "
                        "--scatterer, 0 si no). Electrodos distintos "
                        "(ikind=2 de pwcond.x) no están soportados")
    p.add_argument("--emin", type=float, default=-3.0)
    p.add_argument("--emax", type=float, default=3.0)
    p.add_argument("--points", type=int, default=61)
    p.add_argument("--nz1", type=int, default=3)
    p.add_argument("--collect", action="store_true")
    p.add_argument("--pseudo-dir")
    p.add_argument("--pseudo", action="append", metavar="EL=UPF")
    p.add_argument("--ecutwfc", type=float)
    p.add_argument("--ecutrho", type=float)
    p.add_argument("--kspacing", type=float)
    _fig_opts_min(p)

    p = sub.add_parser("doctor",
                       help="diagnosticar un cálculo: convergencia, fuerzas "
                            "y por qué no converge")
    p.add_argument("path", nargs="?", help="carpeta del cálculo o archivo de salida")
    p.add_argument("--system", action="store_true",
                   help="revisar instalación, recursos, QE y pseudopotenciales")
    p.add_argument("--project", help="además, revisar la puerta de calidad de este proyecto")
    p.add_argument("--json", action="store_true", help="imprimir el diagnóstico como JSON")
    p.add_argument("--prefix")
    p.add_argument("-o", "--outdir", default=".")
    p.add_argument("--no-plot", action="store_true")
    p = p.add_argument_group(GRUPO_FIGURA)
    p.add_argument("--dpi", type=int)
    p.add_argument("--format", default="pdf,png")
    p.add_argument("-t", "--template")
    p.add_argument("--font", choices=["sans", "serif", "latex"])
    p.add_argument("--usetex", action="store_true")
    p.add_argument("--palette")
    p.add_argument("--background")
    p.add_argument("--journal", default="generic",
                   choices=sorted(["acs", "aps", "elsevier", "generic",
                                   "iop", "nature", "rsc", "wiley"]))
    p.add_argument("--width")
    p.add_argument("--mono", action="store_true")

    p = sub.add_parser("audit",
                       help="verificar que un conjunto de cálculos sea "
                            "comparable antes de restar energías")
    p.add_argument("paths", nargs="+")
    p.add_argument("--index", action="store_true",
                   help="además, registrarlos en la base de datos")
    p.add_argument("--db", default="olla-dft.db")

    p = sub.add_parser("esm",
                       help="superficies cargadas con medio de "
                            "apantallamiento efectivo: función trabajo, "
                            "capacitancia y potencial de carga cero")
    _calc_opts(p, "esm")
    p.add_argument("--bc", default="bc1", choices=("bc1", "bc2", "bc3"),
                   help="bc1 vacío/vacío (losas neutras), bc2 metal/metal "
                        "(condensador), bc3 vacío/metal (electrodo, el "
                        "único junto con bc2 que admite carga neta)")
    p.add_argument("--charge", default="0",
                   help="cargas netas en e, separadas por coma: -0.2,0,0.2")
    p.add_argument("--field", type=float, default=0.0,
                   help="campo aplicado en Ry/u.a. (solo con bc2)")
    p.add_argument("--esm-w", dest="width_esm", type=float, default=0.0,
                   help="desplazamiento de la frontera de ESM en u.a.")
    p.add_argument("--nfit", type=int, default=4,
                   help="puntos de ajuste del potencial en la frontera")

    p = sub.add_parser("kappa",
                       help="conductividad térmica de red: fc3, ecuación de "
                            "Boltzmann de fonones y recorrido libre medio")
    p.add_argument("file", help="estructura (celda primitiva)")
    p.add_argument("-o", "--outdir", default="kappa")
    p.add_argument("--dim", default="2x2x2", metavar="NxNxN",
                   help="supercelda de la fc3 (default 2x2x2). Es lo que "
                        "más cuesta: el número de configuraciones crece "
                        "deprisa")
    p.add_argument("--dim-fc2", dest="dim_fc2", default=None, metavar="NxNxN",
                   help="supercelda MAYOR solo para la parte armónica, que "
                        "es barata y necesita más alcance")
    p.add_argument("--distance", type=float, default=0.03,
                   help="desplazamiento finito en Å (default 0.03)")
    p.add_argument("--mesh", type=int, default=13,
                   help="malla de q para la ecuación de Boltzmann "
                        "(default 13)")
    p.add_argument("--temps", default="100:800:8",
                   help="temperaturas: 100:800:8 o 300,500,700")
    p.add_argument("--isotopes", action="store_true",
                   help="añadir dispersión por isótopos con las abundancias "
                        "naturales (en Si son ~10 %%)")
    p.add_argument("--grain", type=float, default=None, metavar="UM",
                   help="tamaño de grano en µm: añade dispersión por "
                        "fronteras")
    p.add_argument("--model", default=None,
                   help="calcular las fuerzas con un potencial aprendido "
                        "(mace, chgnet, m3gnet) en vez de con pw.x: "
                        "segundos en vez de horas, pero el valor absoluto "
                        "puede estar lejos")
    p.add_argument("--collect", action="store_true",
                   help="leer las fuerzas ya calculadas y resolver")
    p.add_argument("--force", action="store_true",
                   help="escribir los inputs aunque sean muchísimos")
    p.add_argument("--metal", action="store_true",
                   help="sistema metálico (ocupaciones con smearing en los "
                        "scf de la fc2/fc3). Sin esto se usa "
                        "occupations='fixed', lo correcto para aislantes")
    p.add_argument("--pseudo-dir")
    p.add_argument("--ecutwfc", type=float)
    p.add_argument("--ecutrho", type=float)
    p.add_argument("--kspacing", type=float, default=0.35)
    _fig_opts_min(p)

    p = sub.add_parser("berry",
                       help="polarización por fase de Berry: ΔP espontánea, "
                            "cargas de Born y ferroelectricidad")
    p.add_argument("file", help="estructura (la polar, si hay camino)")
    p.add_argument("-o", "--outdir", default="berry")
    p.add_argument("--gdir", type=int, default=3, choices=(1, 2, 3),
                   help="dirección: vector de la red recíproca (default 3)")
    p.add_argument("--nppstr", type=int, default=9,
                   help="puntos por cuerda de k (default 9); sube hasta que "
                        "la fase deje de moverse")
    p.add_argument("--kperp", default="6x6", metavar="NxN",
                   help="malla perpendicular a la cuerda (default 6x6)")
    p.add_argument("-r", "--reference", default=None, metavar="ARCHIVO",
                   help="estructura de referencia, normalmente la "
                        "centrosimétrica: se interpola un camino adiabático "
                        "hasta la polar y ΔP es la polarización espontánea")
    p.add_argument("--displace", default=None, metavar="ATOMO:dx,dy,dz",
                   help="camino de desplazamiento de un átomo, en Å; la "
                        "pendiente de P da la carga efectiva de Born")
    p.add_argument("--nlambda", type=int, default=5,
                   help="puntos del camino (default 5)")
    p.add_argument("--run", action="store_true")
    p.add_argument("--collect", action="store_true")
    p.add_argument("--redo", action="store_true")
    p.add_argument("--pw-cmd", dest="pw_cmd")
    p.add_argument("--nproc", type=int)
    p.add_argument("--timeout", type=float)
    p.add_argument("--pseudo-dir")
    p.add_argument("--ecutwfc", type=float)
    p.add_argument("--ecutrho", type=float)
    p.add_argument("--kgrid", default=None, metavar="NxNxN")
    _fig_opts_min(p)

    p = sub.add_parser("wannier",
                       help="funciones de Wannier: interpolar bandas, "
                            "centros y dispersión, sin necesitar wannier90")
    p.add_argument("file", nargs="?", help="estructura")
    p.add_argument("-o", "--outdir", default="wannier")
    p.add_argument("-g", "--grid", default="4x4x4", metavar="NxNxN",
                   help="malla COMPLETA de puntos k (default 4x4x4). Es la "
                        "que fija la calidad de la interpolación")
    p.add_argument("-p", "--projections", default="auto",
                   metavar="SITIO:ORBITAL",
                   help="orbitales de prueba: 'Si:sp3', 'O:p;Ti:d', "
                        "'f=0.125,0.125,0.125:s'. Varias separadas por ';'. "
                        "Con 'auto' se ponen s y p en cada átomo")
    p.add_argument("--bands", type=int, default=None,
                   help="bandas del nscf (default: las que hagan falta)")
    p.add_argument("--exclude", default=None, metavar="5-8",
                   help="bandas que NO entran en la wannierización")
    p.add_argument("--window", default=None, metavar="MIN:MAX",
                   help="ventana exterior de desenredado en eV: de qué "
                        "bandas se puede elegir el subespacio. Hace falta "
                        "cuando las bandas están enredadas con otras "
                        "(conducción, metales)")
    p.add_argument("--frozen", default=None, metavar="MIN:MAX",
                   help="ventana congelada en eV: las bandas de dentro se "
                        "reproducen EXACTAS. Suele ser la valencia más el "
                        "trozo de conducción que te importe")
    p.add_argument("--no-minimize", dest="no_minimize", action="store_true",
                   help="quedarse en la gauge de proyección, sin minimizar "
                        "la dispersión")
    p.add_argument("--iterations", type=int, default=500,
                   help="pasos de minimización (default 500)")
    p.add_argument("--points", type=int, default=30,
                   help="puntos por tramo del camino interpolado")
    p.add_argument("--dft-bands", dest="dft_bands", default=None,
                   metavar="DIR",
                   help="carpeta con el cálculo de bandas de DFT con el que "
                        "comparar; sin esto no hay validación de verdad")
    p.add_argument("--no-dft-bands", dest="no_dft_bands", action="store_true",
                   help="con --run, saltarse el paso 4 de bandas")
    p.add_argument("--dos", type=int, default=None, metavar="N",
                   help="además, DOS interpolada en una malla NxNxN")
    p.add_argument("--sigma", type=float, default=0.05,
                   help="ensanchamiento de la DOS interpolada (eV)")
    p.add_argument("--run", action="store_true",
                   help="lanzar los cuatro pasos en orden")
    p.add_argument("--collect", action="store_true",
                   help="analizar lo que ya está corrido")
    p.add_argument("--pw-cmd", dest="pw_cmd")
    p.add_argument("--pw2wan-cmd", dest="pw2wan_cmd",
                   help="ejecutable de pw2wannier90.x (default: al lado de "
                        "pw.x)")
    p.add_argument("--nproc", type=int)
    p.add_argument("--timeout", type=float)
    p.add_argument("--pseudo-dir")
    p.add_argument("--ecutwfc", type=float)
    p.add_argument("--ecutrho", type=float)
    p.add_argument("--kgrid", default=None, metavar="NxNxN",
                   help="malla del scf inicial")
    p.add_argument("--insulator", action="store_true")
    _fig_opts_min(p)

    p = sub.add_parser(
        "topology",
        help="Chern y lazos de Wilson de un modelo Wannier")
    p.add_argument(
        "model", metavar="MODELO",
        help="archivo *_hr.dat o carpeta que contenga WANNIER_hr.dat")
    occupation = p.add_mutually_exclusive_group(required=True)
    occupation.add_argument(
        "--occupied", type=int, metavar="N",
        help="número de bandas ocupadas del subespacio aislado")
    occupation.add_argument(
        "--fermi", type=float, metavar="EV",
        help="nivel de Fermi; se rechaza si corta una banda")
    p.add_argument(
        "-g", "--grid", default="40x40", metavar="NxN",
        help="malla periódica de la sección 2D (default: 40x40)")
    p.add_argument(
        "--plane", choices=sorted(("xy", "xz", "yz")), default="xy",
        help="plano orientado de la sección del BZ (default: xy)")
    p.add_argument(
        "--fixed", type=float, default=0.0, metavar="K",
        help="coordenada fraccionaria perpendicular (default: 0)")
    p.add_argument(
        "--gap-tol", dest="gap_tol", type=float, default=1e-8, metavar="EV",
        help="gap directo mínimo para aceptar el invariante (default: 1e-8)")
    p.add_argument("-o", "--outdir", default="topology")
    _fig_opts_min(p)

    p = sub.add_parser("amorphous",
                       help="sólido amorfo por fundido y temple con un "
                            "potencial aprendido")
    p.add_argument("formula", help="fórmula de la unidad, por ejemplo SiO2")
    p.add_argument("-n", "--units", type=int, default=8,
                   help="unidades de fórmula en la celda (default: 8)")
    p.add_argument("-d", "--density", type=float, required=True,
                   metavar="G_CM3", help="densidad objetivo en g/cm³")
    p.add_argument("--melt", type=float, default=3000.0, metavar="K",
                   help="temperatura de fundido (default: 3000 K)")
    p.add_argument("--final", type=float, default=300.0, metavar="K",
                   help="temperatura final (default: 300 K)")
    p.add_argument("--melt-steps", dest="melt_steps", type=int, default=500)
    p.add_argument("--quench-steps", dest="quench_steps", type=int,
                   default=1000,
                   help="pasos del temple: son los que fijan la velocidad. "
                        "El default (1000) es un temple de exploración a "
                        "~3e15 K/s y el reporte lo avisa; 27000 baja a "
                        "1e14 K/s")
    p.add_argument("--anneal-steps", dest="anneal_steps", type=int, default=200)
    p.add_argument("--dt", type=float, default=1.0, metavar="FS")
    p.add_argument("--model", default="mace", help="potencial interatómico")
    p.add_argument("--min-dist", dest="min_dist", type=float,
                   default=None, metavar="F",
                   help="factor sobre la suma de radios covalentes al "
                        "empaquetar (default: 0.75)")
    p.add_argument("--seed", type=int, default=0,
                   help="semilla; cambia para generar otra realización")
    p.add_argument("--pack-only", dest="pack_only", action="store_true",
                   help="solo empaquetar, sin dinámica")
    p.add_argument("-o", "--outdir", default="amorfo")

    p = sub.add_parser("docs",
                       help="referencia navegable de todos los subcomandos, "
                            "generada del propio código")
    p.add_argument("-o", "--output", default="olla-dft-docs.html",
                   help="archivo HTML de salida")
    p.add_argument("--open", dest="abrir", action="store_true",
                   help="abrirla en el navegador al terminar")
    p.add_argument("--language", choices=["es", "en"], default=argparse.SUPPRESS,
                   help="idioma de la interfaz de referencia (default: es)")
    p.add_argument("--both", action="store_true",
                   help="generar referencias en español e inglés por separado")

    p = sub.add_parser("echem",
                       help="electrodo de hidrógeno computacional: HER, OER, "
                            "potencial limitante y sobrepotencial")
    p.add_argument("--her", type=float, metavar="E_ads",
                   help="energía de adsorción de H en eV (reacción HER)")
    p.add_argument("--oer", metavar="OH=..,O=..,OOH=..",
                   help="energías de adsorción de los tres intermedios de la "
                        "OER, en eV y referidas al agua")
    p.add_argument("--corrections", metavar="X=eV",
                   help="correcciones térmicas ZPE−TΔS por intermedio; sin "
                        "esto se usan las estándar de la literatura")
    p.add_argument("-U", "--potential", type=float, default=0.0,
                   help="potencial aplicado en V frente al SHE (a pH 0 es "
                        "el mismo que frente al RHE; el pH lo convierte)")
    p.add_argument("--ph", type=float, default=0.0, help="pH")
    p.add_argument("-T", "--temperature", type=float, default=298.15)
    p.add_argument("-o", "--outdir", default="echem")
    _fig_opts_min(p)

    p = sub.add_parser("selftest",
                       help="comprobar Olla-DFT contra valores publicados, no "
                            "contra sí mismo")
    p.add_argument("--full", action="store_true",
                   help="incluir las pruebas que corren pw.x de verdad "
                        "(unos diez minutos)")
    p.add_argument("--mlip", action="store_true",
                   help="incluir por separado la prueba con potencial "
                        "aprendido (requiere MACE)")
    p.add_argument("--only", help="solo estas pruebas, separadas por coma")
    p.add_argument("--list", action="store_true",
                   help="listar las pruebas y sus referencias, sin correr nada")
    p.add_argument("--pseudo-dir", help="pseudopotenciales para las de --full")
    p.add_argument("--pw-cmd", dest="pw_cmd")
    p.add_argument("--nproc", type=int)
    p.add_argument("-j", "--jobs", type=int, default=1)
    p.add_argument("--keep", metavar="CARPETA",
                   help="dejar los cálculos aquí en vez de borrarlos")

    p = sub.add_parser("cost",
                       help="qué sabe Olla-DFT de la velocidad de tu máquina")
    p.add_argument("--db", default="olla-dft.db", help="base de cálculos")

    p = sub.add_parser("db", help="índice local de cálculos")
    p.add_argument("paths", nargs="*", help="carpetas a registrar")
    p.add_argument("--db", default="olla-dft.db")
    p.add_argument("-q", "--query", help="consulta SQL (solo SELECT)")
    p.add_argument("--export", help="exportar todo a un JSON")
    p.add_argument("--formula", help="filtrar por fórmula, por ejemplo Si")
    p.add_argument("--calculation", help="filtrar por tipo: scf, relax, nscf...")
    p.add_argument("--gap-min", type=float, help="gap mínimo en eV")
    p.add_argument("--gap-max", type=float, help="gap máximo en eV")
    p.add_argument("--limit", type=int, default=100, help="máximo de filas filtradas")

    p = sub.add_parser("compare",
                       help="comparar corridas sin restar energías incompatibles")
    p.add_argument("paths", nargs="+", help="carpetas o XML de las corridas")
    p.add_argument("--reference", type=int, default=0,
                   help="índice de la corrida de referencia (default: 0)")
    p.add_argument("-o", "--output", help="guardar comparación en JSON")

    p = sub.add_parser("tune",
                       help="recomendar el siguiente punto de una convergencia")
    p.add_argument("file", help="CONVERGENCIA.dat")
    p.add_argument("--threshold", type=float, default=None,
                   help="umbral en meV/átomo (default: 1)")
    p.add_argument("-o", "--output", help="guardar recomendación en JSON")

    p = sub.add_parser("hull",
                       help="energías de formación y casco convexo")
    p.add_argument("paths", nargs="+")
    p.add_argument("-o", "--outdir", default=".")
    p.add_argument("--elements", help="orden de los elementos, por ejemplo Zn,Al")
    p.add_argument("--threshold", type=float, default=0.025,
                   help="umbral de metaestabilidad en eV/átomo")
    p.add_argument("--force", action="store_true",
                   help="construir el casco aunque la auditoría falle")
    p.add_argument("--dpi", type=int)
    p.add_argument("--format", default="pdf,png")
    p.add_argument("--no-plot", action="store_true")
    p.add_argument("-t", "--template")
    p.add_argument("--font", choices=["sans", "serif", "latex"])
    p.add_argument("--usetex", action="store_true")
    p.add_argument("--palette")
    p.add_argument("--background")
    p.add_argument("--journal", default="generic",
                   choices=sorted(["acs", "aps", "elsevier", "generic",
                                   "iop", "nature", "rsc", "wiley"]))
    p.add_argument("--width")
    p.add_argument("--mono", action="store_true")

    p = sub.add_parser("report",
                       help="registro local de fallas y confusiones")
    p.add_argument("description", nargs="*",
                   help="qué pasó (si se omite, lista las incidencias)")
    p.add_argument("--show", help="ver una incidencia por su id")
    p.add_argument("--close", help="marcar una incidencia como resuelta")
    p.add_argument("--note", help="nota al cerrar")
    p.add_argument("--stats", action="store_true",
                   help="qué subcomandos fallan más")
    p.add_argument("--export", help="empaquetar todo en un archivo JSON")
    p.add_argument("--only-open", action="store_true")
    p.add_argument("--attach", action="append",
                   help="adjuntar un archivo (se copia al registro local)")

    p = sub.add_parser("mlip",
                       help="potencial aprendido: pre-relajar y cribar "
                            "antes de gastar DFT")
    p.add_argument("action", choices=["relax", "scan", "phonons"])
    p.add_argument("file")
    p.add_argument("-o", "--output", help="estructura de salida (relax)")
    p.add_argument("--model", default="mace",
                   choices=["mace", "chgnet", "m3gnet"])
    p.add_argument("--size", default="small",
                   help="tamaño del modelo MACE (small/medium/large)")
    p.add_argument("--device", default="cpu")
    p.add_argument("--fmax", type=float, default=0.01,
                   help="fuerza objetivo en eV/Å")
    p.add_argument("--steps", type=int, default=300)
    p.add_argument("--fixed-cell", action="store_true",
                   help="no relajar la celda, solo las posiciones")
    p.add_argument("--span", type=float, default=0.10,
                   help="rango del barrido de volumen (scan)")
    p.add_argument("--npoints", type=int, default=15)
    p.add_argument("--supercell", help="supercelda del cribado, ej. 2x2x2")

    p = sub.add_parser("suggest",
                       help="sugerir parámetros a partir de tus cálculos "
                            "previos")
    p.add_argument("file")
    p.add_argument("--db", default="olla-dft.db")

    p = sub.add_parser("crosscheck",
                       help="cruzar la misma cantidad por rutas "
                            "independientes")
    p.add_argument("project", nargs="?", default=".",
                   help="carpeta del proyecto")
    p.add_argument("-f", "--file", help="estructura (para masas y volumen)")
    p.add_argument("--gap-bandas", type=float)
    p.add_argument("--gap-tauc", type=float)

    p = sub.add_parser("derived",
                       help="Debye, velocidades del sonido y Slack desde "
                            "las Cij")
    p.add_argument("file", help="estructura")
    p.add_argument("--cij", default="ELASTIC_C.dat",
                   help="archivo con la matriz elástica")
    p.add_argument("--temp", type=float, default=300.0)
    p.add_argument("-o", "--outdir", default=".",
                   help="carpeta donde dejar DERIVED.dat")

    p = sub.add_parser("datasheet",
                       help="ficha del material y párrafo de métodos")
    p.add_argument("project", nargs="?", default=".")
    p.add_argument("-o", "--outdir", default=".")
    p.add_argument("--name", help="nombre base de los archivos")
    p.add_argument("--methods", action="store_true",
                   help="solo el párrafo de metodología y las citas")

    p = sub.add_parser("qha",
                       help="cuasi-armónica: expansión térmica y a(T)")
    p.add_argument("data", help="tabla: V(A^3) E(eV) w1 w2 ... por volumen")
    p.add_argument("-o", "--outdir", default=".")
    p.add_argument("--natoms", type=int, default=1)
    p.add_argument("--cells", type=int, default=1,
                   help="celdas primitivas por supercelda de los modos")
    p.add_argument("--cubic", action="store_true",
                   help="además, a(T). Sin --structure es solo V_prim^(1/3)")
    p.add_argument("--structure", metavar="CIF",
                   help="estructura del material: con ella a(T) se convierte "
                        "al parámetro de red CONVENCIONAL (factor 4 en "
                        "fcc/diamante, 2 en bcc) y se detecta si es cúbica")
    p.add_argument("--tmax", type=float, default=1000.0)
    p.add_argument("--dt", type=float, default=5.0)
    p.add_argument("--temp", type=float, default=300.0)
    p.add_argument("--dpi", type=int)
    p.add_argument("--format", default="pdf,png")
    p.add_argument("--no-plot", action="store_true")
    p.add_argument("-t", "--template")
    p.add_argument("--font", choices=["sans", "serif", "latex"])
    p.add_argument("--usetex", action="store_true")
    p.add_argument("--palette")
    p.add_argument("--background")
    p.add_argument("--journal", default="generic",
                   choices=sorted(["acs", "aps", "elsevier", "generic",
                                   "iop", "nature", "rsc", "wiley"]))
    p.add_argument("--width")
    p.add_argument("--mono", action="store_true")

    p = sub.add_parser("transport",
                       help="Seebeck, sigma/tau y factor de potencia (CRTA)")
    _calc_opts(p, "transporte")
    p.add_argument("--grid", help="malla del nscf, por ejemplo 16x16x16")
    p.add_argument("--temperatures", default="300",
                   help="temperaturas en K separadas por comas")
    p.add_argument("--mu-span", type=float, default=1.0,
                   help="rango de potencial químico alrededor de E_F (eV)")
    p.add_argument("--metal", action="store_true")
    p.add_argument("--nspin", type=int, default=1, choices=[1, 2],
                   help="2 activa la polarización de espín en scf y nscf "
                        "(necesario para --spin-resolved)")
    p.add_argument("--mag", help="magnetización inicial, por ejemplo "
                                 "Fe=0.7 (implica --nspin 2)")

    p.add_argument("--spin-resolved", dest="spin_resolved",
                   action="store_true",
                   help="separar los dos canales de espín (modelo de dos "
                        "corrientes) y dar la polarización de la "
                        "conductividad y la termopotencia de espín")
    p = sub.add_parser("phonons",
                       help="fonones DFPT: dispersión, DOS, termodinámica, IR")
    _calc_opts(p, "fonones")
    p.add_argument("--qgrid", help="malla de q, por ejemplo 2x2x2")
    p.add_argument("--gamma", action="store_true",
                   help="solo Γ con dynmat.x: frecuencias y actividades IR")
    p.add_argument("--raman", action="store_true",
                   help="además, tensores e intensidades Raman en Γ "
                        "(lraman; solo pseudos de norma conservada, y es "
                        "bastante más caro)")
    p.add_argument("--laser", type=float, default=532.0,
                   help="longitud de onda del láser en nm para simular el "
                        "espectro Raman")
    p.add_argument("--suite", action="store_true",
                   help="además, exportar JSON de intercambio (solo con "
                        "--gamma) para las apps de FTIR y Raman")
    p.add_argument("--tscan", metavar="T1,T2,...",
                   help="barrido de temperatura ELECTRÓNICA en K: repite los "
                        "fonones con smearing fermi-dirac a cada una y mira "
                        "si un modo imaginario se estabiliza al calentar "
                        "(ondas de densidad de carga, transiciones "
                        "estructurales)")

    p = sub.add_parser("templates", help="listar, ver o exportar plantillas")
    p.add_argument("action", nargs="?", choices=["list", "show", "export"],
                   help="list (default), show o export")
    p.add_argument("name", nargs="?", help="nombre de la plantilla")
    p.add_argument("-o", "--output", help="archivo JSON de salida (export)")

    p = sub.add_parser("config", help="ver o cambiar la configuración")
    p.add_argument("action", nargs="?", choices=["show", "set"])
    p.add_argument("key", nargs="?")
    p.add_argument("value", nargs="?")

    # argparse imprime por omisión los 70+ comandos como una lista plana.
    # Se conserva _choices_actions porque la documentación extrae de ahí las
    # descripciones; solo se oculta esa lista al formatear la ayuda principal.
    sub._get_subactions = lambda: []
    parser._positionals.title = i18n.ui("arguments", language)
    parser._optionals.title = i18n.ui("global_options", language)
    for action in parser._actions:
        if action.dest == "help":
            action.help = i18n.ui("show_help", language)
        elif action.help:
            action.help = i18n.translate(action.help, language)
    for name, child in sub.choices.items():
        child._positionals.title = i18n.ui("arguments", language)
        child._optionals.title = i18n.ui("options", language)
        for group in child._action_groups[2:]:
            group.title = i18n.ui(group.title, language)
        if child.description:
            child.description = i18n.translate(child.description, language)
        for action in child._actions:
            if action.dest == "help":
                action.help = i18n.ui("show_help", language)
                continue
            if action.help is argparse.SUPPRESS:
                continue
            if not action.help:
                # Ayuda por defecto para las banderas que se repiten en
                # decenas de comandos (--dpi, --pseudo-dir, --outdir...).
                flag = action.option_strings[-1] if action.option_strings \
                    else f"<{action.dest}>"
                action.help = i18n.default_help(flag) or None
            if action.help and language != "es":
                action.help = i18n.translate(action.help, language)
    # El resumen de cada subcomando en la ayuda general y en `docs`.
    if language != "es":
        for choice in sub._choices_actions:
            choice.help = i18n.translate(choice.help, language)

    return parser


def _cmd_resilient(args):
    from qekit.modules import resilient
    return resilient.cli(args)


_DISPATCH = {
    "gen": _cmd_gen,
    "info": _cmd_info,
    "kpath": _cmd_kpath,
    "prim": _cmd_prim,
    "conv": _cmd_conv,
    "supercell": _cmd_supercell,
    "convert": _cmd_convert,
    "bands": _cmd_bands,
    "dos": _cmd_dos,
    "plot": _cmd_plot,
    "gap": _cmd_gap,
    "converge": _cmd_converge,
    "eos": _cmd_eos,
    "elastic": _cmd_elastic,
    "strain": _cmd_strain,
    "adsorb": _cmd_adsorb,
    "eform": _cmd_eform,
    "gamma": _cmd_gamma,
    "align": _cmd_align,
    "layers": _cmd_layers,
    "xrd": _cmd_xrd,
    "exfoliate": _cmd_exfoliate,
    "wf": _cmd_wf,
    "charge": _cmd_charge,
    "optics": _cmd_optics,
    "effmass": _cmd_effmass,
    "surface": _cmd_surface,
    "defect": _cmd_defect,
    "charges": _cmd_charges,
    "fermi": _cmd_fermi,
    "xps": _cmd_xps,
    "corehole": _cmd_corehole,
    "xanes": _cmd_xanes,
    "hubbard": _cmd_hubbard,
    "interface": _cmd_interface,
    "md": _cmd_md,
    "neb": _cmd_neb,
    "thermochem": _cmd_thermochem,
    "unfold": _cmd_unfold,
    "elph": _cmd_elph,
    "sistema": _cmd_sistema,
    "update": _cmd_update,
    "teoria": _cmd_teoria,
    "recetas": _cmd_recetas,
    "wizard": _cmd_wizard,
    "start": _cmd_start,
    "project": _cmd_project,
    "resilient": _cmd_resilient,
    "pseudos": _cmd_pseudos,
    "tddft": _cmd_tddft,
    "ballistic": _cmd_ballistic,
    "doctor": _cmd_doctor,
    "audit": _cmd_audit,
    "esm": _cmd_esm,
    "kappa": _cmd_kappa,
    "berry": _cmd_berry,
    "wannier": _cmd_wannier,
    "topology": _cmd_topology,
    "amorphous": _cmd_amorphous,
    "docs": _cmd_docs,
    "echem": _cmd_echem,
    "selftest": _cmd_selftest,
    "cost": _cmd_cost,
    "db": _cmd_db,
    "compare": _cmd_compare,
    "tune": _cmd_tune,
    "results": _cmd_results,
    "campaign": _cmd_campaign,
    "hull": _cmd_hull,
    "report": _cmd_report,
    "mlip": _cmd_mlip,
    "suggest": _cmd_suggest,
    "crosscheck": _cmd_crosscheck,
    "derived": _cmd_derived,
    "datasheet": _cmd_datasheet,
    "qha": _cmd_qha,
    "transport": _cmd_transport,
    "phonons": _cmd_phonons,
    "templates": _cmd_templates,
    "config": _cmd_config,
}


# ======================================================================
# Menú interactivo
# ======================================================================
def _ask(prompt: str, default: str = None) -> str:
    suffix = f" [{default}]" if default else ""
    try:
        answer = input(f"{prompt}{suffix}: ").strip()
    except (EOFError, KeyboardInterrupt):
        print()
        sys.exit(0)
    return answer or (default or "")


def _ask_structure(language="es"):
    labels = _menu_section(language, "structure")
    while True:
        fname = _ask(labels["file"])
        if not fname:
            return None, None
        try:
            return fname, structure.load(fname)
        except Exception as exc:
            print(f"{labels['error']} {exc}")


def _menu_recetas(language="es"):
    """Las recetas desde el menú, para quien no sabe ni que existen."""
    from qekit.modules import recipes as rec
    labels = _menu_section(language, "recipes")

    print()
    print(rec.listar(language))
    print()
    print(f"  {labels['intro']}")
    q = _ask(labels["prompt"])
    if not q.strip():
        return
    if q.strip().lower() in rec.RECETAS_POR_CLAVE:
        r = rec.obtener(q.strip().lower(), language)
    else:
        cand = rec.buscar(q, language=language)
        if not cand:
            print(f"  {labels['not_found']}")
            return
        r = cand[0]
        if len(cand) > 1:
            print(f"  {labels['show']} «{r.clave}»; {labels['also_lower']}: "
                  + ", ".join(x.clave for x in cand[1:3]))
    print()
    print(rec.report(r, language))


def _menu_asistente(language="es"):
    """El asistente desde el menú: de lo que quieres SABER a los comandos."""
    from qekit.modules import wizard as wz
    labels = _menu_section(language, "wizard")

    print()
    print(f"  {labels['intro']}")
    print(f"  {labels['examples']}")
    q = _ask(labels["prompt"])
    if not q.strip():
        return
    cands = wz.buscar(q, language=language)
    if not cands:
        print(f"  {labels['not_found']}")
        return
    m = cands[0]
    if len(cands) > 1:
        print(f"  {labels['show']} «{m.clave}»; {labels['also_lower']}: "
              + ", ".join(x.clave for x in cands[1:3]))
    fname, atoms = _ask_structure(language)
    print()
    print(wz.report_meta(m, fname or "estructura.cif", language=language))


def _menu_gen(language="es"):
    labels = _menu_section(language, "generate")
    fname, atoms = _ask_structure(language)
    if atoms is None:
        return
    print("\n" + "\n".join(labels["items"]))
    choice = _ask(labels["choice"], "101")
    preset = PRESET_MENU.get(choice)
    if preset is None:
        print(labels["invalid"])
        return
    outdir = _ask(labels["output"], ".")
    insulator = _ask(labels["insulator"], "n").lower().startswith(("s", "y"))
    mag = _ask(labels["magnetization"], "")
    args = argparse.Namespace(
        file=fname, preset=preset, outdir=outdir, klevel=None, kspacing=None,
        band_points=None, ecutwfc=None, ecutrho=None, insulator=insulator,
        primitive=False, pseudo_dir=None, prefix=None,
        nspin=2 if mag else 1, mag=mag or None,
    )
    print()
    try:
        _cmd_gen(args)
    except Exception as exc:
        print(f"{labels['error']} {exc}")


def _menu_structure_tools(language="es"):
    labels = _menu_section(language, "structure")
    print("\n" + "\n".join(labels["items"]))
    choice = _ask(labels["choice"], "401")
    fname, atoms = _ask_structure(language)
    if atoms is None:
        return
    try:
        if choice == "401":
            print(structure.info_text(atoms))
        elif choice == "402":
            out = _ask(labels["output"], "primitive.cif")
            prim = structure.primitive(atoms)
            print(f"  {labels['written']}: {structure.convert(prim, out)} "
                  f"({len(prim)} {labels['atoms']})")
        elif choice == "403":
            out = _ask(labels["output"], "conventional.cif")
            conv = structure.conventional(atoms)
            print(f"  {labels['written']}: {structure.convert(conv, out)} "
                  f"({len(conv)} {labels['atoms']})")
        elif choice == "404":
            reply = _ask(labels["factors"], "2 2 2").split()
            nx, ny, nz = (int(x) for x in reply[:3])
            out = _ask(labels["output"], "supercell.cif")
            sc = structure.supercell(atoms, nx, ny, nz)
            print(f"  {labels['written']}: {structure.convert(sc, out)} "
                  f"({len(sc)} {labels['atoms']})")
        elif choice == "405":
            out = _ask(labels["output_format"], "out.cif")
            print(f"  {labels['written']}: {structure.convert(atoms, out)}")
        else:
            print(labels["invalid"])
    except Exception as exc:
        print(f"{labels['error']} {exc}")


def _menu_postproc(language="es"):
    labels = _menu_section(language, "postprocess")
    print("\n" + "\n".join(labels["items"]))
    choice = _ask(labels["choice"], "201")
    path = _ask(labels["path"], ".")
    args = argparse.Namespace(
        path=path, outdir=".", prefix=None, ref="auto",
        emin=-6.0, emax=6.0, dpi=600, format="pdf,png",
        no_plot=False, mode="orbital",
        template=None, size=None, font=None, usetex=False, palette=None,
        background=None, journal="generic", width=None,
        aspect=None, mono=False, dashes="auto", title=None,
        gap_label=False, panel=None,
    )
    if choice not in ("204",):
        args.outdir = _ask(labels["output"], ".")
        rango = _ask(labels["range"], "-6 6").split()
        try:
            args.emin, args.emax = float(rango[0]), float(rango[1])
        except (IndexError, ValueError):
            print(labels["invalid_range"])
        tpl = _ask(labels["template"] + " (" + "/".join(qthemes.names()) + ")", "journal")
        if tpl:
            args.template = tpl
        rev = _ask(f"{labels['journal']} "
                   f"({'/'.join(sorted(qstyle.JOURNALS))})", "generic")
        if rev in qstyle.JOURNALS:
            args.journal = rev

        if _ask(labels["label_gap"], "n").lower().startswith(("s", "y")):
            args.gap_label = True
    print()
    try:
        if choice == "201":
            _cmd_bands(args)
        elif choice == "202":
            _cmd_dos(args)
        elif choice == "203":
            _cmd_plot(args)
        elif choice == "204":
            _cmd_gap(args)
        else:
            print(labels["invalid"])
    except Exception as exc:
        print(f"{labels['error']} {exc}")


def _menu_config(language="es"):
    labels = _menu_section(language, "config")
    print()
    print(qcfg.show())
    print()
    if _ask(labels["change"], "n").lower().startswith(("s", "y")):
        key = _ask(f"{labels['key']} ({', '.join(qcfg.VALID_KEYS)})")
        if key:
            value = _ask(f"{labels['value']} {key}")
            try:
                qcfg.set_value(key, value)
                print(labels["saved"])
            except KeyError as exc:
                print(f"{labels['error']} {exc}")


def _ask_yes(prompt: str, default: bool = False) -> bool:
    d = "s/N" if not default else "S/n"
    ans = _ask(f"{prompt} ({d})").lower()
    if not ans:
        return default
    return ans[0] in ("s", "y")


def _run_cli(argv: list, language="es") -> None:
    """Muestra el comando equivalente y lo ejecuta.

    Se enseña la línea de comandos para que el menú sirva también de guía:
    quien lo use unas cuantas veces ya puede escribir el comando directo o
    meterlo en un script.
    """
    labels = _menu_labels(language)
    shown = " ".join(a if " " not in a else f'"{a}"' for a in argv)
    print(f"\n  {labels['equivalent']}:\n    olla-dft {shown}\n")
    try:
        main(argv)
    except SystemExit:
        pass
    except Exception as exc:
        print(f"{labels['error']} ({type(exc).__name__}): {exc}")


def _menu_calc(language="es"):
    """Barridos: convergencia, ecuación de estado y elásticas."""
    labels = _menu_section(language, "calculation")
    print("\n" + "\n".join(labels["items"]))
    choice = _ask(labels["choice"])
    if choice not in ("1", "2", "3"):
        return
    fname, atoms = _ask_structure(language)
    if atoms is None:
        return
    outdir = _ask(labels["output"], {"1": "convergencia", "2": "eos",
                                     "3": "elastic"}[choice])
    if choice == "1":
        kind = _ask(labels["kind"], "ecutwfc")
        argv = ["converge", fname, "-k", kind, "-o", outdir]
    elif choice == "2":
        argv = ["eos", fname, "-o", outdir]
    else:
        argv = ["elastic", fname, "-o", outdir]

    print(f"\n  {labels['prepare']}")
    if _ask_yes(labels["run"]):
        argv.append("--run")
    elif _ask_yes(labels["collect"]):
        argv.append("--collect")
    _run_cli(argv, language)


def _menu_props(language="es"):
    """Ópticas, fonones, densidad de carga y función trabajo."""
    labels = _menu_section(language, "properties")
    print("\n" + "\n".join(labels["items"]))
    choice = _ask(labels["choice"])
    if choice not in ("1", "2", "3", "4"):
        return
    fname, atoms = _ask_structure(language)
    if atoms is None:
        return

    if choice == "1":
        argv = ["optics", fname, "-o", _ask(labels["optics_output"], "opticas")]
        sc = _ask(labels["scissor"])
        if sc:
            argv += ["--scissor", sc]
        print(f"\n  {labels['optics_note']}")
    elif choice == "2":
        argv = ["phonons", fname, "-o", _ask(labels["phonons_output"], "fonones")]
        if _ask_yes(labels["gamma"]):
            argv.append("--gamma")
        else:
            q = _ask(labels["qgrid"])
            if q:
                argv += ["--qgrid", q]
        print(f"\n  {labels['phonons_note']}")
    elif choice == "3":
        argv = ["charge", fname, "-o", _ask(labels["charge_output"], "campos")]
        plot = _ask(labels["plot"], "density")
        argv += ["--plot", plot]
    else:
        argv = ["wf", fname, "-o",
                _ask(labels["work_function_output"], "funcion_trabajo")]

    if choice in ("1", "2"):
        if _ask_yes(labels["run"]):
            argv.append("--run")
        elif _ask_yes(labels["collect"]):
            argv.append("--collect")
    _run_cli(argv, language)


def _menu_layered(language="es"):
    """Capas, difracción de polvos y exfoliación."""
    labels = _menu_section(language, "layered")
    print("\n" + "\n".join(labels["items"]))
    choice = _ask(labels["choice"])
    if choice not in ("1", "2", "3"):
        return
    fname, atoms = _ask_structure(language)
    if atoms is None:
        return
    if choice == "1":
        argv = ["layers", fname]
        slab = _ask(labels["slab"])
        if slab:
            argv += ["--slab", slab]
    elif choice == "2":
        argv = ["xrd", fname, "-o", _ask(labels["xrd_output"], ".")]
        lam = _ask(labels["radiation"], "CuKa")
        argv += ["--wavelength", lam]
        size = _ask(labels["size"])
        if size:
            argv += ["--size", size]
        exp = _ask(labels["experimental"])
        if exp:
            argv += ["--exp", exp]
    else:
        argv = ["exfoliate", fname, "-o",
                _ask(labels["exfoliation_output"], "exfoliacion")]
        vdw = _ask(labels["vdw"], "grimme-d2")
        if vdw and vdw != "none":
            argv += ["--vdw", vdw]
        if _ask_yes(labels["run"]):
            argv.append("--run")
        elif _ask_yes(labels["collect"]):
            argv.append("--collect")
    _run_cli(argv, language)


def _menu_catalog(language="es"):
    """Catálogo navegable: área primero y opciones completas al elegir."""
    labels = _menu_section(language, "catalog")
    print()
    print(_catalog_text(language))
    command = _ask(labels["prompt"]).strip()
    if not command:
        return
    parsers = build_parser()._subparsers._group_actions[0].choices
    if command not in parsers:
        print(f"{labels['not_found']} '{command}'.")
        return
    print()
    parsers[command].print_help()


def interactive_menu(language="es"):
    labels = _menu_labels(language)
    print(
        f"""
 ============================================================
    {labels['title'].format(version=__version__)}
    {labels['subtitle']}
 ============================================================"""
    )
    while True:
        print("\n" + "\n".join(labels["items"]))
        choice = _ask(labels["choice"])
        if choice in ("0", "q", "salir", "exit"):
            print(labels["goodbye"])
            return
        elif choice.lower() in ("r", "recetas"):
            _menu_recetas(language)
        elif choice.lower() in ("p", "proyecto", "inicio", "start"):
            _run_cli(["start", "--language", language])
        elif choice.lower() in ("a", "asistente", "wizard"):
            _menu_asistente(language)
        elif choice.lower() in ("c", "catalogo", "catálogo", "comandos"):
            _menu_catalog(language)
        elif choice.lower() in ("t", "teoria", "teoría", "theory"):
            cmd = _ask(_menu_section(language, "theory")["prompt"]).strip()
            _run_cli(["teoria"] + ([cmd] if cmd else []) +
                     ["--language", language], language)
        elif choice == "1":
            _menu_gen(language)
        elif choice == "2":
            _menu_postproc(language)
        elif choice == "3":
            _menu_calc(language)
        elif choice == "4":
            _menu_props(language)
        elif choice == "5":
            _menu_layered(language)
        elif choice == "6":
            fname, atoms = _ask_structure(language)
            if atoms is not None:
                print()
                print(structure.info_text(atoms))
        elif choice == "7":
            fname, atoms = _ask_structure(language)
            if atoms is not None:
                print()
                try:
                    print(kpoints.kpath_text(kpoints.get_kpath(atoms)))
                except Exception as exc:
                    print(f"  Error: {exc}")
        elif choice == "8":
            _menu_structure_tools(language)
        elif choice == "9":
            _menu_config(language)
        else:
            print(labels["invalid"])



def _validar_estilo(args) -> None:
    """Rechaza las banderas de figura MAL ESCRITAS antes de calcular nada.

    Antes, un '-t nature' se detectaba hasta el momento de dibujar: para
    entonces el difractograma ya estaba calculado y escrito, y el usuario
    se quedaba con media salida y un error. Como validar una plantilla no
    cuesta nada, se hace de entrada.
    """
    if not any(hasattr(args, k) for k in
               ("template", "journal", "palette", "font", "width")):
        return
    if getattr(args, "no_plot", False):
        return
    qthemes.load(getattr(args, "template", None) or None,
                 family=getattr(args, "font", None),
                 palette=getattr(args, "palette", None),
                 background=getattr(args, "background", None))
    journal = getattr(args, "journal", None)
    width = getattr(args, "width", None)
    if journal is not None or width is not None:
        qstyle.width_mm(width or "single", journal or "generic")


# ======================================================================
def main(argv=None) -> int:
    from qekit.core import consola, provenance
    crudo = argv if argv is not None else sys.argv[1:]
    # Antes que nada: dejar la salida en condiciones. Si no, el primer print
    # con una Å mata el comando en una consola de Windows heredada, y el
    # usuario cree que falló el cálculo.
    ascii_ = "--ascii" in list(crudo)
    consola.preparar(forzar_ascii=ascii_)
    provenance.record_argv(argv if argv is not None else sys.argv)
    # --ascii y --language son globales pero la gente las escribe donde le
    # sale: detrás del subcomando, delante, en medio. Se consumen aquí y se
    # quitan de la lista para que argparse no las vea fuera de sitio.
    limpio = [t for t in crudo if t != "--ascii"]
    try:
        limpio, idioma = i18n.extract_language(limpio)
    except ValueError as exc:
        print(f"Error: --language admite es o en, no '{exc}'.", file=sys.stderr)
        return 2
    idioma = i18n.set_language(idioma)
    parser = build_parser(idioma)
    args = parser.parse_args(_pegar_negativos(limpio))
    args.language = idioma
    args.command = ALIASES.get(args.command, args.command)
    if args.command is None:
        interactive_menu(idioma)
        return 0
    # en modo --collect el cálculo ya corrió: no se reescriben los inputs
    # (el usuario pudo haberlos editado a mano o usado otros parámetros).
    # Se fija en cada llamada, no solo cuando toca apagarlo: el menú
    # interactivo llama a main() varias veces en el mismo proceso y un
    # --collect previo dejaría mudos los prepare siguientes.
    from qekit.modules import sweep as _sweep
    _sweep.set_write_inputs(
        not (getattr(args, "collect", False)
             and not getattr(args, "run", False))
    )
    try:
        # Dentro del try a proposito: si esto lanza un ErrorDeUso fuera, se
        # escapa del manejador y sale una traza en vez del mensaje limpio.
        _sweep.set_pseudo_overrides(_pseudos_forzados(args))
        _validar_estilo(args)
        _validar_ejecucion(args)
        return _DISPATCH[args.command](args)
    except BrokenPipeError:
        # Ocurre al encadenar la salida con `head`, `less` y demás: no es un
        # fallo del programa, así que se cierra en silencio.
        try:
            sys.stdout.close()
        except Exception:
            pass
        return 0
    except KeyboardInterrupt:
        print("\nInterrumpido.", file=sys.stderr)
        return 130
    except FileNotFoundError as exc:
        print(f"Error: {exc}", file=sys.stderr)
        return 1
    except ErrorDeUso as exc:
        # El programa hizo lo correcto: el comando o el dato no encajaban y
        # el mensaje ya dice qué hacer. No es una falla, así que no se
        # archiva con traza ni se alarma al usuario; queda anotado como
        # "uso" nada más para poder ver luego qué banderas confunden.
        print(f"Error: {exc}", file=sys.stderr)
        try:
            from qekit.modules import feedback
            feedback.registrar(exc=exc, tipo="uso")
        except Exception:                              # noqa: BLE001
            pass
        return 2
    except Exception as exc:
        print(f"Error ({type(exc).__name__}): {exc}", file=sys.stderr)
        # Un fallo inesperado se registra solo, con el comando, la traza y
        # las versiones. Un reporte que dice "fallo optics" no sirve; uno
        # con esto se arregla en minutos.
        try:
            from qekit.modules import feedback
            inc = feedback.registrar(exc=exc)
            print(f"\nSe registró la incidencia {inc.id} con el comando, la "
                  "traza y las versiones.\n  Verla:     olla-dft report --show "
                  f"{inc.id}\n  Exportar:  olla-dft report --export "
                  "incidencias.json", file=sys.stderr)
        except Exception:                              # noqa: BLE001
            pass
        return 1


if __name__ == "__main__":
    sys.exit(main())
