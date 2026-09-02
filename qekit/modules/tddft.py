# Olla-DFT — command-line toolkit for Quantum ESPRESSO
# Copyright (C) 2026 Jorge Enrique González Sevilla
# SPDX-License-Identifier: GPL-3.0-or-later
# This program is free software: you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by the Free
# Software Foundation, either version 3 of the License, or (at your option)
# any later version. See the LICENSE file for details.

"""Absorción óptica con TDDFPT: la parte que epsilon.x no ve.

QUÉ AÑADE SOBRE `olla-dft optics`
------------------------------
`epsilon.x` calcula la función dieléctrica en la aproximación de
PARTÍCULAS INDEPENDIENTES: suma transiciones verticales entre bandas
ocupadas y vacías, sin dejar que el electrón excitado y el hueco que deja
se vean entre sí.

Esa interacción existe y a veces domina. En un semiconductor de gap ancho
o en una molécula, el par electrón-hueco queda ligado —un excitón— y su
energía cae POR DEBAJO del gap. El resultado práctico:

- aparece un pico de absorción dentro del gap, que epsilon.x no tiene;
- el borde de absorción se corre a energías más bajas;
- las intensidades relativas cambian, a veces mucho.

TDDFPT resuelve la respuesta lineal dependiente del tiempo con el método
de Lanczos: no calcula estados vacíos explícitamente, así que es mucho
más barato de lo que uno esperaría para lo que da.

LOS DOS CAMINOS
---------------
- **Lanczos** (`turbo_lanczos.x` + `turbo_spectrum.x`): da el espectro
  entero de un tirón. Es lo que se quiere para una curva de absorción.
- **Davidson** (`turbo_davidson.x`): da las primeras N excitaciones una a
  una, con sus energías y fuerzas de oscilador. Es lo que se quiere para
  decir "el primer estado excitado está a 4.2 eV y es de tipo n->pi*".

LO QUE HAY QUE SABER
--------------------
1. **El número de iteraciones de Lanczos manda la resolución.** Con pocas,
   el espectro sale suave y falso. Se sube `--iter` hasta que los picos
   dejen de moverse; la extrapolación (`--extrapolation`) ayuda pero no
   sustituye.

2. **El funcional decide si hay excitón o no.** Con LDA o GGA, el kernel
   de intercambio-correlación adiabático NO liga excitones en un sólido:
   el espectro sale casi igual que el de partículas independientes. Para
   ver excitones de verdad en un cristal hace falta un funcional híbrido o
   ir a Bethe-Salpeter. En MOLÉCULAS, en cambio, TDDFT con GGA ya mejora
   bastante. Olla-DFT lo dice en el reporte en vez de dejar creer que el
   método ve algo que no ve.

3. **En un sólido periódico hay que fijar la polarización.** Con `--pol 4`
   se hacen las tres cadenas y sale el tensor completo.

4. **Una molécula necesita una caja grande.** Si las imágenes periódicas
   se ven, el espectro se ensucia. Olla-DFT mide el vacío y avisa.
"""

import re
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

from qekit.core import provenance
from qekit.core import style as qstyle
from qekit.core.errors import ErrorDeUso
from qekit.modules import sweep

RY_EV = 13.605693122994
#: Ensanchamiento por omisión (eV) de spectrum.in / davidson.in.
BROADENING_DEFAULT = 0.05
#: Nombres de la polarización por índice de ipol.
EJES = {1: "xx", 2: "yy", 3: "zz", 4: "tensor completo"}


@dataclass
class TddftRun:
    energias: np.ndarray = None      # eV
    total: np.ndarray = None         # coeficiente de absorción / S(omega)
    componentes: dict = field(default_factory=dict)   # 'x','y','z'
    picos: list = field(default_factory=list)         # (E_eV, altura rel.)
    excitaciones: list = field(default_factory=list)  # (E_eV, fuerza) Davidson
    polarizaciones: list = field(default_factory=list)  # (fx, fy, fz) por estado
    itermax: int = None
    broadening: float = None
    funcional: str = None
    gap_ip: float = None             # gap de partículas independientes, si se da
    metodo: str = "lanczos"
    avisos: list = field(default_factory=list)

    @property
    def onset(self) -> float:
        """Borde de absorción: el punto de INFLEXIÓN de la primera subida.

        Es la definición que se usa en el laboratorio, y la única robusta
        aquí. Tomar "donde supera el 5 % del máximo" parece razonable y no
        lo es: la cola gaussiana de un pico centrado justo en el gap cruza
        ese 5 % a 2.5 sigmas por debajo, así que un espectro SIN excitón
        parecería tenerlo solo por el ensanchamiento.
        """
        if self.total is None or self.energias is None:
            return float("nan")
        s = np.asarray(self.total, dtype=float)
        e = np.asarray(self.energias, dtype=float)
        if s.size < 5 or s.max() <= 0:
            return float("nan")
        d = np.gradient(s, e)
        # el primer maximo local de la derivada que sea apreciable
        umbral = 0.2 * d.max()
        for i in range(1, len(d) - 1):
            if d[i] >= d[i - 1] and d[i] > d[i + 1] and d[i] > umbral:
                return float(e[i])
        return float(e[int(np.argmax(d))])


# ----------------------------------------------------------------------
# Inputs
# ----------------------------------------------------------------------
def build_lanczos_input(prefix: str, itermax: int = 500, ipol: int = 4,
                        n_ipol: int = None, ltammd: bool = False,
                        lrpa: bool = False, no_hxc: bool = False,
                        scissor: float = 0.0) -> str:
    if n_ipol is None:
        n_ipol = 3 if ipol == 4 else 1
    lineas = [" &lr_input",
              f"   prefix = '{prefix}',", "   outdir = './out',",
              "   restart = .false.,", "   restart_step = 250,",
              " /", " &lr_control",
              f"   itermax = {itermax},",
              f"   ipol = {ipol},",
              f"   n_ipol = {n_ipol},"]
    if ltammd:
        lineas.append("   ltammd = .true.,")
    if lrpa:
        lineas.append("   lrpa = .true.,")
    if no_hxc:
        lineas.append("   no_hxc = .true.,")
    if scissor:
        lineas.append(f"   scissor = {scissor / RY_EV:.6f},")
    lineas.append(" /")
    return "\n".join(lineas) + "\n"


def build_spectrum_input(prefix: str, itermax: int = 500,
                         itermax0: int = None, emin: float = 0.0,
                         emax: float = 15.0, paso: float = 0.01,
                         broadening: float = 0.05, ipol: int = 4,
                         extrapolation: str = "osc") -> str:
    """Input de turbo_spectrum.x. Las energías en eV (units = 1)."""
    if extrapolation not in ("no", "constant", "osc"):
        raise ErrorDeUso(
            f"extrapolación '{extrapolation}' desconocida. Opciones: no, "
            "constant, osc.")
    lineas = [" &lr_input",
              f"   prefix = '{prefix}',", "   outdir = './out',",
              f"   itermax0 = {itermax0 or itermax},",
              f"   itermax = {max(itermax, itermax0 or 0) * 4},",
              f"   extrapolation = '{extrapolation}',",
              f"   epsil = {broadening / RY_EV:.6f},",
              "   units = 1,",
              f"   start = {emin},",
              f"   end = {emax},",
              f"   increment = {paso},",
              f"   ipol = {ipol},",
              " /"]
    return "\n".join(lineas) + "\n"


def build_davidson_input(prefix: str, n_estados: int = 10,
                         emin: float = 0.0, emax: float = 15.0,
                         paso: float = 0.005, broadening: float = 0.05,
                         nbnd_occ: int = None,
                         nbnd_virt: int = 15) -> str:
    lineas = [" &lr_input",
              f"   prefix = '{prefix}',", "   outdir = './out',",
              "   restart = .false.,", " /", " &lr_dav",
              f"   num_eign = {n_estados},",
              f"   num_init = {2 * n_estados},",
              f"   num_basis_max = {max(80, 8 * n_estados)},",
              "   residue_conv_thr = 1.0D-4,",
              f"   start = {emin / RY_EV:.6f},",
              f"   finish = {emax / RY_EV:.6f},",
              f"   step = {paso / RY_EV:.8f},",
              f"   broadening = {broadening / RY_EV:.6f},",
              f"   reference = {0.5 * (emin + emax) / RY_EV:.6f},",
              f"   p_nbnd_virt = {nbnd_virt},"]
    if nbnd_occ:
        lineas.append(f"   p_nbnd_occ = {nbnd_occ},")
    lineas.append(" /")
    return "\n".join(lineas) + "\n"


def _vacio_minimo(atoms) -> float:
    celda = np.array(atoms.get_cell(), dtype=float)
    pos = atoms.get_positions()
    huecos = []
    for eje in range(3):
        largo = float(np.linalg.norm(celda[eje]))
        if largo <= 0:
            continue
        proy = pos @ (celda[eje] / largo)
        huecos.append(largo - (proy.max() - proy.min()))
    return min(huecos) if huecos else 0.0


def prepare(atoms, outdir: str = "tddft", metodo: str = "lanczos",
            itermax: int = 500, ipol: int = 4, n_estados: int = 10,
            emin: float = 0.0, emax: float = 15.0, paso: float = 0.01,
            broadening: float = 0.05, extrapolation: str = "osc",
            pseudo_dir: str = None, ecutwfc: float = None,
            ecutrho: float = None, kspacing: float = None,
            insulator: bool = True, nbnd: int = None,
            ltammd: bool = False, lrpa: bool = False,
            gamma: bool = None, scissor: float = 0.0) -> tuple:
    """Escribe scf.in y los inputs de TDDFPT.

    `scissor` (eV) es el corrimiento rígido de las bandas vacías que
    turbo_lanczos.x aplica antes de la respuesta: sirve para compensar la
    subestimación del gap del funcional. Solo existe en el método lanczos.
    """
    from qekit.modules import inputgen

    if metodo not in ("lanczos", "davidson"):
        raise ErrorDeUso(
            f"método '{metodo}' desconocido. Opciones: lanczos (espectro "
            "completo), davidson (las primeras excitaciones una a una).")
    scissor = float(scissor or 0.0)
    if scissor < 0:
        raise ErrorDeUso(
            f"--scissor {scissor} no tiene sentido: el corrimiento abre el "
            "gap, así que es cero o positivo.")
    if scissor and metodo != "lanczos":
        raise ErrorDeUso(
            "--scissor solo existe en turbo_lanczos.x; turbo_davidson.x no "
            "lo admite. Usa --method lanczos o quita --scissor.")

    common = sweep.prepare_common(atoms, pseudo_dir, ecutwfc, ecutrho,
                                  insulator)
    out = Path(outdir); out.mkdir(parents=True, exist_ok=True)
    # TDDFPT solo tiene implementado el caso GAMMA: con una malla de k se
    # planta con "k-point algorithm is not tested yet". Y para una molecula
    # gamma es ademas lo correcto. Ojo: "K_POINTS gamma" NO es lo mismo que
    # una malla 1x1x1 — QE usa un algoritmo distinto, y es el unico que
    # TDDFPT acepta.
    if gamma is None:
        gamma = _vacio_minimo(atoms) > 5.0
    if gamma:
        kcard = "K_POINTS gamma\n"
        grid = (1, 1, 1)
    else:
        grid = sweep.default_grid(atoms, kspacing)
        kcard = (f"K_POINTS automatic\n  {grid[0]} {grid[1]} {grid[2]} "
                 "0 0 0\n")
    scf = inputgen.build_pw_input(
        atoms=atoms, pseudos=common["pseudos"], calculation="scf",
        prefix=common["prefix"], pseudo_dir=common["pseudo_dir"],
        ecutwfc=common["ecutwfc"], ecutrho=common["ecutrho"],
        kcard=kcard, insulator=insulator, degauss=common["degauss"],
        smearing=common["smearing"], nbnd=nbnd, nosym=True)
    # TDDFPT NO admite simetria: sin esto, turbo_*.x se planta con
    # "Linear response calculation is not implemented with symmetry"
    # despues de que el scf ya haya corrido entero.
    if "noinv" not in scf:
        scf = re.sub(r"(\n\s*nosym\s*=\s*\.true\.\n)",
                     r"\1  noinv            = .true.\n", scf, count=1)
    sweep.write_input(out / "scf.in", scf)

    if metodo == "lanczos":
        sweep.write_input(out / "lanczos.in", build_lanczos_input(
            common["prefix"], itermax=itermax, ipol=ipol, ltammd=ltammd,
            lrpa=lrpa, scissor=scissor))
        sweep.write_input(out / "spectrum.in", build_spectrum_input(
            common["prefix"], itermax=itermax, emin=emin, emax=emax,
            paso=paso, broadening=broadening, ipol=ipol,
            extrapolation=extrapolation))
        orden = ("pw.x -in scf.in  ->  turbo_lanczos.x -in lanczos.in  ->  "
                 "turbo_spectrum.x -in spectrum.in")
    else:
        sweep.write_input(out / "davidson.in", build_davidson_input(
            common["prefix"], n_estados=n_estados, emin=emin, emax=emax,
            broadening=broadening))
        orden = "pw.x -in scf.in  ->  turbo_davidson.x -in davidson.in"

    vac = _vacio_minimo(atoms)
    rep = ["--- Absorción óptica con TDDFPT ---",
           f"Estructura: {atoms.get_chemical_formula()} "
           f"({len(atoms)} átomos)",
           f"Método: {metodo}"
           + (f"   iteraciones de Lanczos: {itermax}"
              if metodo == "lanczos" else
              f"   excitaciones pedidas: {n_estados}"),
           f"Polarización: ipol = {ipol} ({EJES.get(ipol, '?')})",
           f"Ventana: {emin} a {emax} eV   ensanchamiento "
           f"{broadening * 1000:.0f} meV"
           + (f"   scissor {scissor:.3f} eV" if scissor else ""),
           ("Puntos k: solo GAMMA (es lo unico que TDDFPT implementa)"
            if gamma else
            f"Malla k: {grid[0]}x{grid[1]}x{grid[2]}  -- OJO: TDDFPT solo "
            "tiene implementado\nel caso gamma y se plantara al leer el "
            "input"),
           "",
           f"Archivos en '{out.resolve()}':",
           f"Orden:  {orden}",
           ""]
    if ltammd:
        rep += ["Aproximación de Tamm-Dancoff activada: se desprecian los "
                "términos de\ndesexcitación. Abarata el cálculo y suele "
                "corregir poco las energías, pero\nNO es exacta.", ""]
    if lrpa:
        rep += ["RPA: se apaga el kernel de intercambio-correlación. Sirve "
                "para VER cuánto\naporta ese kernel, comparando contra el "
                "cálculo completo.", ""]
    if vac < 6.0 and len(atoms) < 30:
        rep += [f"AVISO: solo hay {vac:.1f} Å de vacío. Si esto es una "
                "molécula, sus imágenes\nperiódicas se ven entre sí y el "
                "espectro sale contaminado: para una molécula\nhacen falta "
                "al menos 8-10 Å por todos lados.", ""]
    rep += ["Recuerda para qué sirve esto: `olla-dft optics` da el espectro de "
            "PARTÍCULAS\nINDEPENDIENTES; TDDFPT deja que el electrón "
            "excitado y su hueco se vean. La\ndiferencia entre los dos "
            "espectros ES el efecto de esa interacción.",
            "",
            "Y una advertencia honesta: con LDA o GGA el kernel adiabático "
            "NO liga\nexcitones en un SÓLIDO, así que el espectro se "
            "parecerá mucho al de\nepsilon.x. En MOLÉCULAS sí mejora. Para "
            "excitones en cristales hace falta un\nhíbrido o Bethe-Salpeter."]
    warn = sweep.missing_pseudo_warning(common)
    if warn:
        rep.append(warn)
    return common, "\n".join(rep)


# ----------------------------------------------------------------------
# Lectura
# ----------------------------------------------------------------------
def _broadening_de_inputs(p: Path):
    """Ensanchamiento (eV) escrito en spectrum.in o davidson.in, o None.

    Los dos inputs lo llevan en Ry ('epsil' en turbo_spectrum, 'broadening'
    en turbo_davidson); aquí se devuelve en eV.
    """
    for nombre, clave in (("spectrum.in", "epsil"),
                          ("davidson.in", "broadening")):
        f = p / nombre
        if not f.exists():
            continue
        m = re.search(clave + r"\s*=\s*([-+\d.eEdD]+)",
                      f.read_text(errors="ignore"))
        if m:
            try:
                return float(m.group(1).replace("d", "e").replace("D", "e")) \
                    * RY_EV
            except ValueError:
                return None
    return None


def collect(path, metodo: str = "lanczos", gap_ip: float = None,
            broadening: float = None) -> TddftRun:
    """Lee el espectro (lanczos) o las excitaciones (davidson) de `path`.

    `broadening` es el ensanchamiento en eV con que se generó el espectro
    (el --broadening de prepare): fija el umbral por debajo del cual un
    corrimiento del borde respecto al gap no se distingue de un excitón.
    Si no se da, se lee de spectrum.in / davidson.in en la misma carpeta.
    """
    p = Path(path)
    if broadening is None:
        broadening = _broadening_de_inputs(p)
    run = TddftRun(metodo=metodo, gap_ip=gap_ip, broadening=broadening)

    if metodo == "davidson":
        return _collect_davidson(p, run)

    dats = sorted(p.glob("*plot*.dat")) + sorted(p.glob("*.plot_chi.dat"))
    if not dats:
        raise ErrorDeUso(
            f"no hay ningún archivo de espectro en {p}. turbo_spectrum.x lo "
            "escribe como\n<prefix>.plot.dat; si no está, revisa "
            "spectrum.out.")
    d = np.loadtxt(dats[0], comments="#")
    if d.ndim == 1:
        d = d.reshape(1, -1)
    # turbo_spectrum con units=1 escribe la energía en eV
    run.energias = d[:, 0]
    run.total = d[:, 1]
    for i, nombre in enumerate("xyz"):
        if d.shape[1] > 2 + i:
            run.componentes[nombre] = d[:, 2 + i]

    salida = p / "lanczos.out"
    if salida.exists():
        texto = salida.read_text(errors="ignore")
        m = re.search(r"itermax\s*=?\s*(\d+)", texto)
        if m:
            run.itermax = int(m.group(1))
        m = re.search(r"Exchange-correlation\s*=\s*(\S+)", texto)
        if m:
            run.funcional = m.group(1)
    run.picos = _picos(run)
    _avisar(run)
    return run


def _collect_davidson(p: Path, run: TddftRun) -> TddftRun:
    """Lee las excitaciones de turbo_davidson.x.

    Salen en <prefix>.eigen, no en la salida de texto: una fila por
    excitacion con la energia en RYDBERG y la fuerza de oscilador total y
    por direccion. Buscar "Excitation No." en el .out no encuentra nada.
    """
    eig = sorted(p.glob("*.eigen"))
    if not eig:
        raise ErrorDeUso(
            f"no hay ningun archivo .eigen en {p}. turbo_davidson.x lo "
            "escribe al terminar\ncon las excitaciones que encontro; si no "
            "esta, revisa davidson.out.")
    d = np.loadtxt(eig[0], comments="#")
    if d.ndim == 1:
        d = d.reshape(1, -1)
    for fila in d:
        e_eV = float(fila[0]) * RY_EV
        f = float(fila[1]) if len(fila) > 1 else float("nan")
        run.excitaciones.append((e_eV, f))
        if len(fila) >= 5:
            # La polarizacion de CADA excitacion. No va en `componentes`,
            # que guarda espectros completos: mezclar seis tripletes con
            # arrays de miles de puntos rompe la exportacion.
            run.polarizaciones.append(tuple(float(x) for x in fila[2:5]))
    if not run.excitaciones:
        raise ErrorDeUso(
            f"{eig[0].name} esta vacio: turbo_davidson.x no convergio "
            "ninguna excitacion.")

    dats = sorted(p.glob("*plot*.dat"))
    if dats:
        pl = np.loadtxt(dats[0], comments="#")
        if pl.ndim == 2 and pl.shape[1] >= 2:
            run.energias = pl[:, 0] * RY_EV
            run.total = pl[:, 1]
    return run


def _picos(run: TddftRun, umbral: float = 0.05) -> list:
    s = np.asarray(run.total, dtype=float)
    e = np.asarray(run.energias, dtype=float)
    if s.size < 3 or s.max() <= 0:
        return []
    sn = s / s.max()
    idx = [i for i in range(1, len(sn) - 1)
           if sn[i] > sn[i - 1] and sn[i] >= sn[i + 1] and sn[i] > umbral]
    return sorted([(float(e[i]), float(sn[i])) for i in idx],
                  key=lambda t: -t[1])


#: Corrimiento minimo (eV) que se considera significativo frente al gap.
#: Por debajo del ensanchamiento no se puede distinguir nada.
UMBRAL_EXCITON = 0.10


def _avisar(run: TddftRun) -> None:
    if run.gap_ip is not None and np.isfinite(run.onset):
        d = run.onset - run.gap_ip
        umbral = max(UMBRAL_EXCITON, 2.0 * (run.broadening or 0.0))
        if d < -umbral:
            run.avisos.append(
                f"El borde de absorción cae {abs(d):.2f} eV POR DEBAJO del "
                f"gap de partículas\nindependientes ({run.gap_ip:.2f} eV). "
                "Eso es exactamente la firma de un excitón\nligado: el par "
                "electrón-hueco tiene menos energía que los dos por "
                "separado.")
        elif abs(d) < umbral:
            run.avisos.append(
                f"El borde cae a {umbral:.2f} eV o menos del gap de "
                "partículas independientes,\nque es el limite de lo que se "
                "puede distinguir con este ensanchamiento. Con\nLDA o GGA "
                "en un sólido es lo esperable: el kernel adiabático no liga "
                "excitones,\ny la diferencia con 'olla-dft optics' será "
                "pequeña.")


def report(run: TddftRun) -> str:
    lines = ["--- Absorción óptica con TDDFPT ---",
             f"Método: {run.metodo}"]
    if run.itermax:
        lines.append(f"Iteraciones de Lanczos: {run.itermax}")

    if run.metodo == "davidson" and run.excitaciones:
        hay_pol = len(run.polarizaciones) == len(run.excitaciones)
        lines += ["", f"{'n':>3s} {'E (eV)':>9s} {'f (osc.)':>11s} "
                  f"{'λ (nm)':>9s}" + ("  pol." if hay_pol else "")]
        for i, (e, f) in enumerate(run.excitaciones, start=1):
            lam = 1239.84 / e if e > 0 else float("nan")
            extra = ""
            if hay_pol:
                p3 = run.polarizaciones[i - 1]
                extra = ("  " + "xyz"[int(np.argmax(p3))]
                         if max(p3) > 0 else "  -")
            lines.append(f"{i:3d} {e:9.4f} {f:11.5f} {lam:9.1f}{extra}")
        brillantes = [(e, f) for e, f in run.excitaciones
                      if np.isfinite(f) and f > 0.01]
        if brillantes:
            e0, f0 = brillantes[0]
            lines += ["",
                      f"Primera excitación con fuerza apreciable: "
                      f"{e0:.3f} eV ({1239.84 / e0:.0f} nm)."]
        oscuras = len(run.excitaciones) - len(brillantes)
        if oscuras:
            lines.append(
                f"  {oscuras} de {len(run.excitaciones)} son OSCURAS "
                "(fuerza de oscilador casi cero):\n  existen, pero no se ven "
                "en un espectro de absorción.")

    if run.total is not None:
        lines += ["",
                  f"Rango: {run.energias[0]:.2f} a {run.energias[-1]:.2f} eV "
                  f"({len(run.energias)} puntos)",
                  f"Borde de absorción (punto de inflexión): "
                  f"{run.onset:.3f} eV "
                  f"({1239.84 / run.onset:.0f} nm)"
                  if np.isfinite(run.onset) and run.onset > 0 else ""]
        if run.picos:
            lines += ["", "Picos (energía y altura relativa):"]
            for e, h in run.picos[:6]:
                lines.append(f"  {e:8.3f} eV  ({1239.84 / e:6.0f} nm)   "
                             f"{h:.2f}")
    if run.componentes:
        lines += ["", f"Componentes: {', '.join(run.componentes)}"]
        anis = _anisotropia(run)
        if anis is not None:
            lines.append(f"Anisotropía entre direcciones: {anis * 100:.1f} % "
                         "del máximo")

    for a in run.avisos:
        lines += ["", a]
    lines += ["",
              "Para ver qué aporta la interacción electrón-hueco, compara "
              "este espectro con\nel de 'olla-dft optics' del MISMO cálculo: "
              "esa diferencia es el efecto.",
              "",
              "Y la resolución la manda el número de iteraciones de Lanczos. "
              "Si los picos se\nmueven al subirlo, todavía no está "
              "convergido."]
    return "\n".join(lines)


def _anisotropia(run: TddftRun):
    comps = [v for v in run.componentes.values()]
    if len(comps) < 2:
        return None
    a = np.array(comps, dtype=float)
    if a.max() <= 0:
        return 0.0
    return float(np.ptp(a, axis=0).max() / a.max())


def export(run: TddftRun, outdir: str = ".") -> list:
    out = Path(outdir); out.mkdir(parents=True, exist_ok=True)
    escritos = []
    cab = provenance.header_plain(
        "absorcion optica con TDDFPT",
        {"metodo": run.metodo, "itermax": run.itermax,
         "gap_particulas_indep_eV": run.gap_ip},
        titulo="TDDFPT")
    if run.total is not None:
        f = out / "TDDFT.dat"
        cols = [run.energias, run.total]
        nombres = ["E(eV)", "S(E)"]
        for k, v in run.componentes.items():
            cols.append(v); nombres.append(f"S_{k}")
        np.savetxt(f, np.column_stack(cols), fmt="%14.6e",
                   header=cab + "\n" + "  ".join(f"{n:>14s}" for n in nombres),
                   comments="# ")
        escritos.append(str(f))
    if run.excitaciones:
        f = out / "TDDFT_EXCITACIONES.dat"
        np.savetxt(f, np.array(run.excitaciones), fmt="%14.6f",
                   header=cab + "\n       E(eV)        fuerza_oscilador",
                   comments="# ")
        escritos.append(str(f))
    txt = out / "TDDFT.txt"
    txt.write_text(report(run) + "\n")
    escritos.append(str(txt))
    return escritos


def plot(run: TddftRun, outfile: str = "tddft", formats="pdf,png",
         comparar=None, theme: str = None, family: str = None,
         background: str = None, palette=None, usetex: bool = None,
         width="single", journal: str = "generic", mono: bool = False,
         dpi: int = None) -> list:
    """Espectro TDDFPT, opcionalmente sobre el de partículas independientes.

    `comparar` es (energias, alpha) de 'olla-dft optics'. Superponerlos es la
    forma de VER lo que aporta la interacción electrón-hueco, que es el
    motivo entero de correr esto.
    """
    try:
        import matplotlib
        matplotlib.use("Agg")
    except ImportError as exc:                          # pragma: no cover
        raise RuntimeError("matplotlib no está instalado.") from exc

    qstyle.apply(theme, family=family, background=background,
                 palette=palette, usetex=usetex, mono=mono)
    fig, ax = qstyle.new_figure(width, journal, 0.75)
    colores = qstyle.palette(4, mono=mono)

    if comparar is not None:
        e2, a2 = comparar
        a2 = np.asarray(a2, dtype=float)
        if a2.max() > 0 and run.total is not None and run.total.max() > 0:
            a2 = a2 / a2.max() * float(np.max(run.total))
        ax.plot(e2, a2, lw=1.1, color=qstyle.INK_FAINT, dashes=[4, 2],
                label="partículas independientes")

    if run.total is not None:
        ax.plot(run.energias, run.total, lw=1.5, color=colores[0],
                label="TDDFPT")
    for i, (k, v) in enumerate(run.componentes.items()):
        ax.plot(run.energias, v, lw=0.8, alpha=0.6,
                **qstyle.style_line(i + 1, 4, mono=mono),
                label=f"$S_{{{k}{k}}}$")
    for e, f in run.excitaciones[:12]:
        if np.isfinite(f) and f > 0:
            ax.vlines(e, 0, f, color=colores[2], lw=1.0, alpha=0.8)

    if run.gap_ip:
        ax.axvline(run.gap_ip, color=colores[3], lw=0.8, dashes=[2, 2])
        ax.text(run.gap_ip, ax.get_ylim()[1] * 0.95,
                " gap IP", fontsize="small", va="top")
    ax.set_xlabel("E (eV)")
    ax.set_ylabel("absorción (u. arb.)")
    ax.set_ylim(bottom=0)
    if run.energias is not None:
        ax.set_xlim(run.energias[0], run.energias[-1])
    ax.legend(frameon=False, fontsize="small")
    return qstyle.save(fig, outfile, formats, dpi=dpi)
