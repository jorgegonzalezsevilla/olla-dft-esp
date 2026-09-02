# Olla-DFT — command-line toolkit for Quantum ESPRESSO
# Copyright (C) 2026 Jorge Enrique González Sevilla
# SPDX-License-Identifier: GPL-3.0-or-later
# This program is free software: you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by the Free
# Software Foundation, either version 3 of the License, or (at your option)
# any later version. See the LICENSE file for details.

"""Infraestructura común de los barridos.

Convergencia, ecuación de estado y constantes elásticas comparten la misma
mecánica: tomar una estructura, generar una familia de cálculos scf que
difieren en un parámetro, correrlos y leer los resultados. Aquí vive esa
parte compartida para no repetirla tres veces.
"""

from pathlib import Path

from ase import Atoms

from qekit import config as qcfg
from qekit.core import plataforma
from qekit.core import kpoints, pseudo
from qekit.core.runner import Job
from qekit.modules import inputgen
from qekit.core.errors import ErrorDeUso


_WRITE_INPUTS = True


def set_write_inputs(flag: bool):
    """Activa o desactiva la escritura de los inputs de QE.

    La CLI la apaga en modo --collect: ahí el usuario ya corrió el cálculo
    (posiblemente con otros parámetros, o editando los archivos a mano) y
    reescribir los inputs perdería ese trabajo y haría que el reporte
    describiera un cálculo distinto del que realmente se ejecutó.
    """
    global _WRITE_INPUTS
    _WRITE_INPUTS = bool(flag)


def writing_inputs() -> bool:
    return _WRITE_INPUTS


def write_input(path, text):
    """Escribe un input de QE respetando set_write_inputs()."""
    path = Path(path)
    if _WRITE_INPUTS:
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(text)
    return path


#: Pseudopotenciales que el usuario forzo con --pseudo EL=archivo.UPF.
#: Es un ajuste global, como _WRITE_INPUTS: la alternativa era pasarlo por
#: trece firmas de modulo y trece manejadores de la CLI.
_PSEUDOS_FORZADOS = {}


def set_pseudo_overrides(forzados: dict = None) -> None:
    """Fija (o limpia) los pseudopotenciales elegidos a mano."""
    global _PSEUDOS_FORZADOS
    _PSEUDOS_FORZADOS = dict(forzados or {})


def pseudo_overrides() -> dict:
    return dict(_PSEUDOS_FORZADOS)


def prepare_common(atoms: Atoms, pseudo_dir: str = None,
                   ecutwfc: float = None, ecutrho: float = None,
                   insulator: bool = False, prefix: str = None,
                   exclude_pseudos=None, tarea: str = None,
                   pseudos_forzados: dict = None) -> dict:
    """Resuelve pseudopotenciales y cutoffs una sola vez para todo el barrido."""
    cfg = qcfg.load()
    # Ruta ABSOLUTA a propósito. Los inputs se escriben en una subcarpeta y
    # pw.x corre desde ahí, así que un '--pseudo-dir ps' relativo se busca
    # dentro de la subcarpeta y falla con "file not found" DESPUÉS de que
    # el usuario ya se fue a hacer otra cosa.
    pdir = str(Path(pseudo_dir or cfg["pseudo_dir"]).expanduser().resolve())
    symbols = atoms.get_chemical_symbols()
    forzados = dict(_PSEUDOS_FORZADOS)
    forzados.update(pseudos_forzados or {})
    pseudos = pseudo.resolve(symbols, pdir, exclude=exclude_pseudos,
                             tarea=tarea, forzados=forzados)
    auto_wfc, auto_rho = pseudo.recommend_cutoffs(
        pseudos, float(cfg["ecutwfc"]), float(cfg["dual"])
    )
    return dict(
        pseudos=pseudos,
        pseudo_dir=pdir,
        ecutwfc=ecutwfc or auto_wfc,
        ecutrho=ecutrho or auto_rho,
        insulator=insulator,
        degauss=float(cfg["degauss"]),
        smearing=cfg["smearing"],
        prefix=prefix or atoms.get_chemical_formula(mode="hill", empirical=True),
        missing=[s for s, p in pseudos.items() if not p["found"]],
    )


def write_scf_job(atoms: Atoms, common: dict, directory: Path, name: str,
                  grid: tuple, ecutwfc: float = None, ecutrho: float = None,
                  meta: dict = None, calculation: str = "scf",
                  conv_thr: float = 1.0e-8,
                  forc_conv_thr: float = 1.0e-4,
                  vdw: str = None, cell_dofree: str = None,
                  nspin: int = 1, magnetization: dict = None,
                  hubbard: dict = None, nbnd: int = None,
                  tot_charge: float = None,
                  tot_magnetization: float = None,
                  dipole_correction=False) -> Job:
    """Escribe un scf en su propia carpeta y devuelve el Job correspondiente."""
    directory = Path(directory)
    directory.mkdir(parents=True, exist_ok=True)
    text = inputgen.build_pw_input(
        atoms=atoms,
        pseudos=common["pseudos"],
        calculation=calculation,
        prefix=common["prefix"],
        pseudo_dir=common["pseudo_dir"],
        ecutwfc=ecutwfc or common["ecutwfc"],
        ecutrho=ecutrho or common["ecutrho"],
        kcard=_kcard(grid),
        insulator=common["insulator"],
        degauss=common["degauss"],
        smearing=common["smearing"],
        conv_thr=conv_thr,
        forc_conv_thr=forc_conv_thr,
        vdw=vdw,
        cell_dofree=cell_dofree,
        nbnd=nbnd,
        tot_charge=tot_charge,
        tot_magnetization=tot_magnetization,
        dipole_correction=dipole_correction,
        nspin=nspin,
        magnetization=magnetization,
        hubbard=hubbard,
    )
    write_input(directory / "pw.in", text)
    return Job(name=name, directory=directory, meta=meta or {})


def check_soc_pseudos(common: dict) -> None:
    """Se niega a preparar un cálculo con SOC si los pseudos no lo permiten.

    lspinorb sobre pseudos escalar-relativistas no da error: da un
    desdoblamiento espín-órbita de CERO, indistinguible de un resultado
    legítimo. Es exactamente el tipo de fallo silencioso que hay que
    atajar en la preparación, no en el análisis.
    """
    malos = []
    for sym, p in common["pseudos"].items():
        if not p.get("found"):
            continue
        rel = (p.get("relativistic") or "").lower()
        if rel != "full":
            malos.append(f"{sym} ({p['filename']}: relativistic="
                         f"{rel or 'sin declarar'})")
    if malos:
        raise ErrorDeUso(
            "el acoplamiento espín-órbita necesita pseudopotenciales "
            "TOTALMENTE RELATIVISTAS (relativistic='full'), y estos no lo "
            "son:\n  " + "\n  ".join(malos) +
            "\n\nCon pseudos escalar-relativistas, lspinorb no falla: "
            "devuelve un desdoblamiento\nespín-órbita de cero que parece "
            "un resultado válido. Busca los pseudos 'rel-'\n"
            "correspondientes (por ejemplo rel-pbe en la tabla SSSP o en "
            "pslibrary)."
        )


def _kcard(grid) -> str:
    if grid is None or tuple(grid) == (1, 1, 1):
        return "K_POINTS gamma\n"
    n1, n2, n3 = grid
    return f"K_POINTS automatic\n  {n1} {n2} {n3} 0 0 0\n"


def default_grid(atoms: Atoms, kspacing: float = None) -> tuple:
    cfg = qcfg.load()
    return kpoints.kgrid_from_spacing(
        atoms, float(kspacing if kspacing is not None else cfg["kspacing"])
    )


def write_run_script(jobs: list, path: Path, nproc: int = None,
                     pw_cmd: str = None) -> str:
    """Escribe run.sh y, al lado, run.py, que corre en cualquier sistema.

    Se escriben LOS DOS a propósito:

      - `run.sh` es lo que espera quien trabaja en Linux o en macOS, y lo que
        se pega en un guion de Slurm sin pensar.
      - `run.py` hace exactamente lo mismo y corre igual en Windows, donde no
        hay bash. Es también el que sabe buscar `pw.exe` además de `pw.x`.

    Dos detalles del `.sh` que solo se notan fuera de Linux: `nproc` es de
    coreutils y en macOS no existe (allí es `sysctl -n hw.ncpu`), y el
    archivo se escribe con finales de línea POSIX a la fuerza, porque un
    `.sh` con CRLF falla en WSL con «bad interpreter: /bin/bash^M», que es un
    error famoso por lo ilegible que es.
    """
    cfg = qcfg.load()
    n = int(nproc if nproc is not None else cfg["nproc"])
    exe = (pw_cmd or cfg.get("pw_cmd") or "pw.x").strip()
    base = Path(path).parent
    carpetas = [str(Path(job.directory).relative_to(base)) for job in jobs]

    lines = [
        "#!/bin/bash",
        "# Generado por Olla-DFT — corre todos los puntos del barrido.",
        "#   ./run.sh      uno tras otro",
        "#   ./run.sh 4    cuatro a la vez",
        "# En Windows, o si prefieres no depender de bash:  python run.py",
        "PAR=${1:-1}",
        "HILOS=" + plataforma.cuenta_nucleos_shell(),
        '[ -z "$HILOS" ] && HILOS=%d' % n,
        'NP=${NPROC:-$(( HILOS / PAR ))}; [ "$NP" -lt 1 ] && NP=1',
        'PW=${PW:-%s}' % exe,
        'if [ "$PAR" -gt 1 ]; then',
        '  echo "Corriendo $PAR a la vez, $NP procesos cada uno."',
        "fi",
        "",
        "correr() {",
        '  cd "$1" || return 1',
        '  echo ">> $1"',
        '  if [ "$NP" -gt 1 ] && command -v mpirun >/dev/null 2>&1; then',
        '    mpirun -np "$NP" "$PW" -in pw.in > pw.out 2>&1',
        '  elif [ "$NP" -gt 1 ] && command -v mpiexec >/dev/null 2>&1; then',
        '    mpiexec -n "$NP" "$PW" -in pw.in > pw.out 2>&1',
        "  else",
        '    "$PW" -in pw.in > pw.out 2>&1',
        "  fi",
        "}",
        "export -f correr",
        "export NP PW",
        "",
        "CARPETAS=(",
    ]
    lines += [f'  "{c}"' for c in carpetas]
    lines += [
        ")",
        "",
        "estado=0",
        'if [ "$PAR" -le 1 ]; then',
        '  for d in "${CARPETAS[@]}"; do',
        '    ( correr "$d" ) || estado=1',
        '  done',
        "else",
        '  if ! printf "%s\\n" "${CARPETAS[@]}" | xargs -P "$PAR" -I{} '
        'bash -c \'correr "{}"\'; then',
        "    estado=1",
        "  fi",
        "fi",
        'echo "Barrido terminado."',
        'exit "$estado"',
        "",
    ]
    plataforma.escribir_script(path, "\n".join(lines) + "\n")
    write_run_py(carpetas, Path(path).with_name("run.py"), n, exe)
    return str(path)


_RUN_PY = r'''#!/usr/bin/env python3
# Corre todos los puntos de este barrido. Generado por Olla-DFT.
#
#     python run.py        uno tras otro
#     python run.py 4      cuatro a la vez
#
# Hace lo mismo que run.sh y funciona igual en Linux, macOS y Windows: no
# necesita bash, ni xargs, ni que los binarios se llamen .x. Si no encuentra
# MPI corre en serie, que en un portatil suele ser lo correcto de todos modos.
import os
import shutil
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

CARPETAS = __CARPETAS__
PW = os.environ.get("PW", "__PW__")
NPROC_DEF = __N__
AQUI = Path(__file__).resolve().parent


def binario(nombre):
    # el ejecutable, probando .x y .exe: en Windows es pw.exe
    base = nombre
    for suf in (".x", ".exe"):
        if base.endswith(suf):
            base = base[: -len(suf)]
    for cand in (nombre, base + ".exe", base + ".x"):
        hallado = shutil.which(cand)
        if hallado:
            return hallado
        if Path(cand).exists():
            return str(Path(cand).resolve())
    return None


def lanzador(np_):
    if np_ <= 1:
        return []
    for cand, bandera in (("mpirun", "-np"), ("mpiexec", "-n")):
        if shutil.which(cand):
            return [cand, bandera, str(np_)]
    return []


def nucleos():
    try:
        return len(os.sched_getaffinity(0))
    except AttributeError:
        return os.cpu_count() or NPROC_DEF


def correr(carpeta, np_, pw):
    d = AQUI / carpeta
    print(">> " + carpeta, flush=True)
    with open(d / "pw.in") as fi, open(d / "pw.out", "w") as fo:
        return subprocess.run(lanzador(np_) + [pw], stdin=fi, stdout=fo,
                              stderr=subprocess.STDOUT,
                              cwd=str(d)).returncode


def main():
    par = int(sys.argv[1]) if len(sys.argv) > 1 else 1
    pw = binario(PW)
    if pw is None:
        print("No encuentro el ejecutable '%s'." % PW)
        print("Ponlo en el PATH, o pasalo asi:  PW=/ruta/a/pw.x python run.py")
        return 2
    np_ = max(1, int(os.environ.get("NPROC", nucleos() // max(par, 1))))
    if par > 1:
        print("Corriendo %d a la vez, %d procesos cada uno." % (par, np_))
    if par <= 1:
        codigos = [correr(c, np_, pw) for c in CARPETAS]
    else:
        with ThreadPoolExecutor(max_workers=par) as ex:
            codigos = list(ex.map(lambda c: correr(c, np_, pw), CARPETAS))
    malos = sum(1 for c in codigos if c != 0)
    print("Barrido terminado. %d de %d fallaron." % (malos, len(CARPETAS)))
    return 1 if malos else 0


if __name__ == "__main__":
    sys.exit(main())
'''


def write_run_py(carpetas, path, nproc=1, exe="pw.x") -> str:
    """El mismo barrido, en Python: corre en Windows sin bash."""
    texto = (_RUN_PY.replace("__CARPETAS__", repr(list(carpetas)))
                    .replace("__PW__", str(exe))
                    .replace("__N__", str(int(nproc))))
    plataforma.escribir_script(path, texto)
    return str(path)


def missing_pseudo_warning(common: dict) -> str:
    if not common["missing"]:
        return ""
    return (
        "ATENCIÓN: faltan pseudopotenciales para "
        + ", ".join(common["missing"])
        + f"\nBúscalos en {common['pseudo_dir']} o cambia la ruta con:\n"
        "  olla-dft config set pseudo_dir /ruta/a/pseudos"
    )
