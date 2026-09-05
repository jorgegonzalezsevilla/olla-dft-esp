# Olla-DFT — command-line toolkit for Quantum ESPRESSO
# Copyright (C) 2026 Jorge Enrique González Sevilla
# SPDX-License-Identifier: AGPL-3.0-or-later
# This program is free software: you can redistribute it and/or modify it
# under the terms of the GNU Affero General Public License as published by the Free
# Software Foundation, either version 3 of the License, or (at your option)
# any later version. See the LICENSE file for details.

"""Conductividad térmica de red: el fonón que se dispersa contra otro fonón.

Los fonones armónicos no conducen calor de forma finita: si el cristal fuera
exactamente armónico, un paquete de fonones viajaría para siempre y κ sería
infinita. La conductividad sale del término CÚBICO de la energía, el que
permite que un fonón se parta en dos o que dos se fundan en uno. Por eso
hace falta la tercera derivada de la energía respecto de los
desplazamientos, la fc3, y por eso este cálculo es tan caro: la fc2 necesita
una derivada por par de átomos y la fc3, una por TRÍO.

En la práctica:

    Φ_ijk = ∂³E / ∂u_i ∂u_j ∂u_k

se obtiene por diferencias finitas, desplazando dos átomos a la vez en una
supercelda. El número de configuraciones crece deprisa con el tamaño de la
supercelda, y ahí está todo el coste. Con las fuerzas ya calculadas, la
ecuación de Boltzmann de fonones en aproximación de tiempo de relajación da

    κ_L = (1/NV) Σ_λ C_λ v_λ ⊗ v_λ τ_λ

que es lo que resuelve phono3py y lo que se reporta aquí.

Este módulo permite dos fuentes de fuerzas, y la diferencia importa:

  - **Quantum ESPRESSO.** Es el cálculo de verdad. En silicio, una supercelda
    2×2×2 de la primitiva son 57 configuraciones de 16 átomos.
  - **Un potencial aprendido (MACE y compañía).** Las mismas 57
    configuraciones en 8 segundos en vez de 40 minutos. Sirve para explorar
    y para decidir el tamaño de la supercelda antes de gastar el cálculo
    caro, pero **el valor absoluto puede estar muy lejos**: con MACE-MP
    pequeño el silicio da 51 W/mK a 300 K frente a los ~140 medidos. La
    dependencia con T sale bien (el 1/T de Umklapp), el número no. Se avisa
    en el informe, siempre.

Lo que hay que converger, y son tres cosas a la vez:

  1. El tamaño de la supercelda de la fc3. Es lo que fija hasta qué distancia
     se tienen en cuenta las interacciones anarmónicas.
  2. La malla de q sobre la que se resuelve la ecuación de Boltzmann.
  3. Y, si se compara con un experimento, los isótopos: el silicio natural
     conduce alrededor de un 10 % menos que el isotópicamente puro, y esa
     diferencia es mayor que muchas de las que se discuten en un artículo.

La aproximación de tiempo de relajación (RTA) subestima κ frente a la
solución exacta de la ecuación de Boltzmann, típicamente un 10-15 % en
silicio y mucho más en materiales con procesos normales dominantes, como el
grafeno o el diamante. No es un error: es la aproximación, y se dice.
"""

from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

from qekit.core.errors import FaltanDatos
from qekit.core import style as qstyle

# recorridos libres medios (Å) sobre los que se acumula κ
RECORRIDOS = np.logspace(0, 7, 141)
# por encima de esto se avisa antes de escribir cientos de inputs
MUCHAS_CONFIGURACIONES = 150


def _phonopy_atoms(atoms):
    from phonopy.structure.atoms import PhonopyAtoms
    return PhonopyAtoms(symbols=atoms.get_chemical_symbols(),
                        cell=np.array(atoms.cell.array, float),
                        scaled_positions=atoms.get_scaled_positions())


def _ase(sc):
    from ase import Atoms
    return Atoms(symbols=list(sc.symbols), positions=np.array(sc.positions),
                 cell=np.array(sc.cell), pbc=True)


def _importar():
    try:
        from phono3py import Phono3py
    except ImportError as exc:
        raise FaltanDatos(
            "hace falta phono3py, que es quien resuelve la ecuación de "
            "Boltzmann de fonones:\n  pip install phono3py\n"
            "Se instala con pip y no necesita compilar Quantum ESPRESSO otra "
            "vez.") from exc
    return Phono3py


@dataclass
class KappaRun:
    """Un cálculo de conductividad térmica de red."""
    formula: str = ""
    dim: tuple = (2, 2, 2)
    dim_fc2: tuple = None
    distancia: float = 0.03
    n_config: int = 0
    n_atomos: int = 0
    fuente: str = ""            # 'Quantum ESPRESSO' o el nombre del potencial
    malla: tuple = (11, 11, 11)
    temperaturas: np.ndarray = None
    kappa: np.ndarray = None    # (nT, 6) en notación de Voigt, W/mK
    isotopos: bool = False
    frontera: float = None      # tamaño de grano en µm, si se puso
    gamma: np.ndarray = None
    frecuencias: np.ndarray = None
    pesos: np.ndarray = None
    velocidades: np.ndarray = None
    cv: np.ndarray = None
    i300: int = None
    avisos: list = field(default_factory=list)
    directorio: str = ""

    @property
    def kappa_media(self):
        """(κ_xx + κ_yy + κ_zz)/3, que es lo que se compara con un policristal."""
        if self.kappa is None:
            return None
        return self.kappa[:, :3].mean(axis=1)


def preparar(atoms, dim=(2, 2, 2), dim_fc2=None, distancia=0.03,
             simetria=1e-5):
    """Genera las configuraciones desplazadas de la fc3 (y de la fc2).

    `dim_fc2` permite una supercelda MAYOR solo para la parte armónica, que
    es mucho más barata y a la vez la que más necesita alcance: las
    constantes de fuerza armónicas de un semiconductor llegan lejos, las
    anarmónicas no tanto. Es la forma estándar de no pagar la fc3 en una
    supercelda enorme.
    """
    Phono3py = _importar()
    ph = Phono3py(_phonopy_atoms(atoms),
                  supercell_matrix=list(map(int, dim)),
                  phonon_supercell_matrix=(list(map(int, dim_fc2))
                                           if dim_fc2 else None),
                  primitive_matrix="auto", symprec=simetria, log_level=0)
    ph.generate_displacements(distance=float(distancia))
    return ph


def configuraciones(ph):
    """Las superceldas desplazadas, como objetos de ASE."""
    sc3 = [_ase(s) for s in ph.supercells_with_displacements if s is not None]
    sc2 = []
    if ph.phonon_supercell_matrix is not None:
        sc2 = [_ase(s) for s in ph.phonon_supercells_with_displacements
               if s is not None]
    return sc3, sc2


def fuerzas_mlip(configs, modelo="mace", verbose=True):
    """Fuerzas con un potencial aprendido. Segundos en vez de horas."""
    from qekit.modules import mlip
    calc = mlip.calculator(modelo)
    fuera = []
    for i, a in enumerate(configs):
        a = a.copy()
        a.calc = calc
        fuera.append(a.get_forces())
        if verbose and (i % 25 == 0 or i == len(configs) - 1):
            print(f"    {i + 1}/{len(configs)}", end="\r", flush=True)
    if verbose:
        print(" " * 30, end="\r")
    return np.array(fuera)


def escribir_inputs(configs, destino, common, prefijo="d", kspacing=None,
                    ecutwfc=None, ecutrho=None, conv_thr=1e-10):
    """Un scf con cálculo de fuerzas por configuración, cada uno en su carpeta."""
    from qekit.modules import inputgen, sweep

    destino = Path(destino); destino.mkdir(parents=True, exist_ok=True)
    carpetas = []
    for i, a in enumerate(configs):
        d = destino / f"{prefijo}{i:04d}"
        d.mkdir(parents=True, exist_ok=True)
        grid = sweep.default_grid(a, kspacing)
        txt = inputgen.build_pw_input(
            atoms=a, pseudos=common["pseudos"], calculation="scf",
            prefix=common["prefix"], pseudo_dir=common["pseudo_dir"],
            ecutwfc=ecutwfc or common["ecutwfc"],
            ecutrho=ecutrho or common["ecutrho"],
            kcard=f"K_POINTS automatic\n  {grid[0]} {grid[1]} {grid[2]} "
                  "0 0 0\n",
            insulator=common["insulator"], degauss=common["degauss"],
            smearing=common["smearing"], conv_thr=conv_thr)
        sweep.write_input(d / "pw.in", txt)
        carpetas.append(d)
    return carpetas


def leer_fuerzas(carpetas, natomos):
    """Fuerzas de cada scf, en eV/Å y en el orden de los átomos del input."""
    from qekit.core import qeout
    fuera, faltan = [], []
    for d in carpetas:
        d = Path(d)
        try:
            res = qeout.read_xml(str(d))
            F = np.array(res.forces, float)
        except Exception:                                   # noqa: BLE001
            faltan.append(d.name)
            continue
        if F is None or F.shape != (natomos, 3):
            faltan.append(d.name)
            continue
        fuera.append(F)
    if faltan:
        raise FaltanDatos(
            f"faltan las fuerzas de {len(faltan)} configuraciones "
            f"({', '.join(faltan[:5])}{'...' if len(faltan) > 5 else ''}).\n"
            f"Sin TODAS no se puede construir la fc3: cada una aporta una "
            f"derivada distinta.")
    return np.array(fuera)


def resolver(ph, fuerzas, fuerzas_fc2=None, malla=11, temperaturas=None,
             isotopos=False, frontera_um=None, simetrizar=True):
    """Construye fc2 y fc3 y resuelve la ecuación de Boltzmann en RTA.

    `frontera_um` añade dispersión por el tamaño de grano, que es lo que
    convierte este cálculo en algo comparable con una película delgada o un
    nanohilo. `isotopos` añade la dispersión por masa con las abundancias
    naturales: en silicio son ~10 %, más de lo que suele discutirse en un
    artículo.
    """
    temperaturas = (np.asarray(temperaturas, float) if temperaturas is not None
                    else np.arange(100.0, 801.0, 100.0))
    ph.forces = np.asarray(fuerzas, float)
    if fuerzas_fc2 is not None and len(np.asarray(fuerzas_fc2)):
        ph.phonon_forces = np.asarray(fuerzas_fc2, float)
    ph.produce_fc3()
    ph.produce_fc2()
    if simetrizar:
        ph.symmetrize_fc3()
        ph.symmetrize_fc2()
    m = (int(malla),) * 3 if np.isscalar(malla) else tuple(int(x) for x in malla)
    ph.mesh_numbers = list(m)
    ph.init_phph_interaction()
    ph.run_thermal_conductivity(
        temperatures=temperaturas,
        is_isotope=bool(isotopos),
        boundary_mfp=(float(frontera_um) * 1e4 if frontera_um else 1e6),
        write_kappa=False)
    return ph.thermal_conductivity, m


def recoger(run, ph, tc, malla):
    """Pasa lo que devuelve phono3py a los campos del KappaRun."""
    run.malla = tuple(malla)
    run.temperaturas = np.asarray(tc.temperatures, float)
    run.kappa = np.asarray(tc.kappa[0], float)          # (nT, 6) Voigt
    run.frecuencias = np.asarray(tc.frequencies, float)
    run.pesos = np.asarray(tc.grid_weights, float)
    try:
        run.gamma = np.asarray(tc.gamma[0], float)
        run.velocidades = np.asarray(tc.group_velocities, float)
        run.cv = np.asarray(tc.mode_heat_capacities, float)
    except Exception:                                       # noqa: BLE001
        pass
    d = np.abs(run.temperaturas - 300.0)
    run.i300 = int(np.argmin(d)) if len(d) else None
    return run


def acumulada(run, iT=None):
    """κ acumulada frente al recorrido libre medio.

    Es la curva que dice si nanoestructurar sirve: si el 80 % de κ lo llevan
    fonones con Λ > 100 nm, un grano de 50 nm corta ese 80 %. Si lo llevan
    fonones de 5 nm, no hay nada que hacer con el tamaño de grano.
    """
    if run.gamma is None or run.velocidades is None or run.cv is None:
        return None, None
    iT = run.i300 if iT is None else int(iT)
    if iT is None:
        return None, None
    v2 = (run.velocidades ** 2).sum(axis=-1) / 3.0        # (grid, banda)
    vmod = np.sqrt((run.velocidades ** 2).sum(axis=-1))
    # Los modos acústicos en Γ tienen Γ = 0 exactamente: τ es infinito y el
    # producto τ·v² es 0·∞ = NaN. No es un fallo, es que esos modos no
    # transportan calor en un cristal infinito; se descartan más abajo. El
    # errestate tiene que cubrir TAMBIÉN los productos, no solo la división.
    with np.errstate(divide="ignore", invalid="ignore"):
        tau = 1.0 / (2.0 * run.gamma[iT])
        L = vmod * tau                                    # Å
        contrib = run.cv[iT] * v2 * tau                   # ∝ κ del modo
    w = run.pesos[:, None] * np.ones_like(contrib)
    ok = np.isfinite(L) & np.isfinite(contrib) & (contrib > 0)
    Lf, cf, wf = L[ok], contrib[ok], w[ok]
    orden = np.argsort(Lf)
    acum = np.cumsum(cf[orden] * wf[orden])
    if acum[-1] <= 0:
        return None, None
    return Lf[orden], acum / acum[-1]


def recorrido_representativo(run, fraccion=0.5, iT=None):
    """El Λ por debajo del cual se acumula `fraccion` de κ."""
    L, a = acumulada(run, iT)
    if L is None:
        return None
    return float(np.interp(float(fraccion), a, L))


def exponente_temperatura(run, T_min=200.0):
    """El n de κ ∝ T^(−n). Los procesos Umklapp puros dan n = 1."""
    if run.kappa is None or run.temperaturas is None:
        return None
    k = run.kappa_media
    m = (run.temperaturas >= T_min) & (k > 0)
    if m.sum() < 3:
        return None
    return float(-np.polyfit(np.log(run.temperaturas[m]), np.log(k[m]), 1)[0])


def report(run) -> str:
    L = [f"--- Conductividad térmica de red: {run.formula} ---",
         f"Fuerzas: {run.fuente}",
         f"Supercelda de la fc3: {run.dim[0]}×{run.dim[1]}×{run.dim[2]}"
         + (f"   |   fc2: {run.dim_fc2[0]}×{run.dim_fc2[1]}×{run.dim_fc2[2]}"
            if run.dim_fc2 else "")
         + f"   |   {run.n_config} configuraciones de {run.n_atomos} átomos",
         f"Malla de q: {run.malla[0]}×{run.malla[1]}×{run.malla[2]}"
         + ("   |   isótopos naturales" if run.isotopos else "")
         + (f"   |   granos de {run.frontera:g} µm" if run.frontera else ""),
         ""]
    if run.kappa is None:
        return "\n".join(L + ["Todavía no hay κ."])
    iso = np.allclose(run.kappa[:, 0], run.kappa[:, 1], rtol=0.02) and \
        np.allclose(run.kappa[:, 0], run.kappa[:, 2], rtol=0.02)
    L += ["   T (K)      κ_xx      κ_yy      κ_zz     media  (W/m·K)"]
    for i, T in enumerate(run.temperaturas):
        k = run.kappa[i]
        marca = "  ←" if i == run.i300 else ""
        L.append(f"  {T:6.0f}  {k[0]:9.2f} {k[1]:9.2f} {k[2]:9.2f} "
                 f"{run.kappa_media[i]:9.2f}{marca}")
    if iso:
        L.append("  (el tensor es isótropo, como corresponde a la simetría "
                 "cúbica)")
    n = exponente_temperatura(run)
    if n is not None:
        L += ["", f"Dependencia con la temperatura: κ ∝ T^−{n:.2f}"]
        if abs(n - 1.0) < 0.25:
            L.append("  Es el T⁻¹ de los procesos Umklapp: por encima de la "
                     "temperatura de Debye,")
            L.append("  cuantos más fonones hay, más se estorban entre sí.")
        else:
            L.append("  Se aleja del T⁻¹ de Umklapp puro; suele significar "
                     "que hay otro canal")
            L.append("  dominante (fronteras, isótopos) o que la malla de q "
                     "no está convergida.")
    L50 = recorrido_representativo(run, 0.5)
    L90 = recorrido_representativo(run, 0.9)
    if L50:
        L += ["", "Recorrido libre medio (a la T más cercana a 300 K):",
              f"  la mitad de κ la llevan fonones con Λ < {L50 / 10:.0f} nm",
              f"  el 90 %,                              Λ < {L90 / 10:.0f} nm",
              "  Es lo que dice si nanoestructurar sirve: un grano más "
              "pequeño que ese Λ corta",
              "  esa parte de κ; uno más grande no hace nada."]
    L += ["", "Lo que NO está incluido, y conviene tener presente:",
          "  · Es RTA, no la solución exacta de la ecuación de Boltzmann. La "
          "RTA subestima κ",
          "    (≈10-15 % en silicio, mucho más en grafeno o diamante).",
          "  · Solo hay dispersión de tres fonones. A alta temperatura los "
          "procesos de cuatro",
          "    fonones bajan κ, y en materiales muy anarmónicos no son un "
          "detalle."]
    if not run.isotopos:
        L.append("  · Sin dispersión por isótopos. El silicio natural conduce "
                 "~10 % menos que el")
        L.append("    isotópicamente puro: si comparas con un experimento, "
                 "pon --isotopes.")
    if run.fuente and "ESPRESSO" not in run.fuente.upper():
        run.avisos.append(
            f"Las fuerzas vienen de {run.fuente}, no de DFT. La forma de "
            f"κ(T) suele salir bien,\n  pero el valor absoluto puede estar "
            f"lejos: con MACE-MP pequeño el silicio da\n  ~51 W/mK a 300 K "
            f"donde el experimento son ~140. Úsalo para elegir la supercelda "
            f"y\n  la malla, y repite con Quantum ESPRESSO antes de "
            f"publicar nada.")
    if int(np.prod(run.dim)) <= 8:
        run.avisos.append(
            f"La supercelda de la fc3 es {run.dim[0]}×{run.dim[1]}×{run.dim[2]}"
            f", que es pequeña. κ tiene que converger\n  en el tamaño de la "
            f"supercelda Y en la malla de q a la vez: sube una, luego la "
            f"otra,\n  y no te fíes hasta que ninguna de las dos mueva el "
            f"resultado.")
    for a in run.avisos:
        L += ["", f"AVISO: {a}"]
    return "\n".join(L)


def export(run, outdir="kappa") -> list:
    out = Path(outdir); out.mkdir(parents=True, exist_ok=True)
    escritos = []
    if run.kappa is not None:
        f = out / "KAPPA.dat"
        np.savetxt(f, np.column_stack([run.temperaturas, run.kappa,
                                       run.kappa_media]), fmt="%12.5f",
                   header="T(K)  kxx  kyy  kzz  kyz  kxz  kxy  media "
                          "(W/m/K)")
        escritos.append(str(f))
    L, a = acumulada(run)
    if L is not None:
        f = out / "KAPPA_recorrido.dat"
        Lg = RECORRIDOS
        ac = np.interp(Lg, L, a)
        np.savetxt(f, np.column_stack([Lg / 10.0, ac]), fmt="%14.6e",
                   header="Lambda(nm)   fraccion acumulada de kappa")
        escritos.append(str(f))
    f = out / "KAPPA.txt"
    f.write_text(report(run) + "\n", encoding="utf-8")
    escritos.append(str(f))
    return escritos


def plot(run, outfile="kappa", formats="pdf,png", theme=None, size=None,
         family=None, background=None, palette=None, usetex=None,
         width="single", journal="generic", aspect=0.72, mono=False,
         dpi=None) -> list:
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError as exc:                              # pragma: no cover
        raise RuntimeError("matplotlib no está instalado.") from exc
    if run.kappa is None:
        raise FaltanDatos("no hay κ que dibujar.")
    st = qstyle.apply(theme, size=size, family=family, background=background,
                      palette=palette, usetex=usetex, mono=mono)
    fig, ax = qstyle.new_figure(width, journal, aspect)
    cols = qstyle.palette(3, mono=mono)
    ax.loglog(run.temperaturas, run.kappa_media, marker="o", ms=4,
              lw=st["line"], color=cols[0], label=r"$\kappa_L$ calculada")
    n = exponente_temperatura(run)
    if n is not None:
        i = run.i300 if run.i300 is not None else 0
        T0, k0 = run.temperaturas[i], run.kappa_media[i]
        ax.loglog(run.temperaturas, k0 * (run.temperaturas / T0) ** (-1.0),
                  lw=st["line"], color=cols[1], dashes=[4.0, 2.0],
                  label=r"$T^{-1}$ (Umklapp)")
    ax.set_xlabel("T (K)")
    ax.set_ylabel(r"$\kappa_L$ (W m$^{-1}$K$^{-1}$)")
    ax.legend(frameon=False, fontsize=st["legend"])
    escritos = qstyle.save(fig, outfile, formats, dpi=dpi, modulo="kappa")
    plt.close(fig)

    L, a = acumulada(run)
    if L is not None:
        fig2, ax2 = qstyle.new_figure(width, journal, aspect)
        ax2.semilogx(L / 10.0, a * 100.0, lw=st["line"], color=cols[0])
        for frac, et in ((0.5, "50 %"), (0.9, "90 %")):
            x = recorrido_representativo(run, frac) / 10.0
            ax2.axvline(x, color=cols[1], lw=0.7, dashes=[3.0, 2.0])
            ax2.annotate(f"{et}: {x:.0f} nm", (x, frac * 100),
                         textcoords="offset points", xytext=(5, -10),
                         fontsize=st["legend"], color=cols[1])
        ax2.set_xlabel(r"recorrido libre medio $\Lambda$ (nm)")
        ax2.set_ylabel(r"% de $\kappa_L$ acumulado")
        ax2.set_ylim(0, 102)
        escritos += qstyle.save(fig2, str(outfile) + "_recorrido", formats,
                                dpi=dpi, modulo="kappa")
        plt.close(fig2)
    return escritos
