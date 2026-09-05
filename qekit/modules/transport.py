# Olla-DFT — command-line toolkit for Quantum ESPRESSO
# Copyright (C) 2026 Jorge Enrique González Sevilla
# SPDX-License-Identifier: AGPL-3.0-or-later
# This program is free software: you can redistribute it and/or modify it
# under the terms of the GNU Affero General Public License as published by the Free
# Software Foundation, either version 3 of the License, or (at your option)
# any later version. See the LICENSE file for details.

"""Transporte electrónico en aproximación de tiempo de relajación constante.

De la estructura de bandas sobre una malla densa salen los coeficientes de
transporte, integrando la distribución de Fermi-Dirac sobre los estados:

    sigma/tau (mu,T) = e^2 * INT  v (x) v  * (-df/dE)  dE
    S (mu,T)         = (1/(e*T)) * INT v (x) v (E-mu) (-df/dE) dE  / (sigma/tau)
    kappa_e/tau      = (1/T) * INT v (x) v (E-mu)^2 (-df/dE) dE  -  S^2 sigma T
    PF               = S^2 * sigma

donde v = (1/hbar) dE/dk son las velocidades de banda, que Olla-DFT obtiene
derivando E(k) numéricamente sobre la malla.

LO QUE ESTO ES Y LO QUE NO ES
-----------------------------
La aproximación de tiempo de relajación constante (CRTA) supone tau
independiente de la energía y del punto k. Por eso:

- **S no depende de tau**: el tiempo de relajación se cancela en el
  cociente. El coeficiente Seebeck que sale de aquí es una predicción
  real, comparable con el experimento.
- **sigma y kappa_e SÍ dependen de tau**, y Olla-DFT no lo calcula: se
  reportan como sigma/tau y kappa_e/tau, en sus unidades (S/(m*s)). Para
  obtener sigma hay que multiplicar por un tau que venga de otro lado
  (ajuste a una medida, o un cálculo de electrón-fonón).
- El factor de potencia hereda esa dependencia: se reporta PF/tau.

Reportar "sigma = X S/m" sin decir de dónde salió tau es el error más
común al usar este tipo de cálculo, así que Olla-DFT no lo hace.

LA MALLA ES EL CUELLO DE BOTELLA
--------------------------------
Las derivadas de E(k) y las integrales necesitan una malla nscf uniforme y
MUY densa — bastante más que la de un scf. La razón: a 300 K la ventana de
-df/dE mide unos 25 meV, y si en ese rango caen pocos estados, sigma sale
como una serie de picos aislados en lugar de una curva. Se ve a simple
vista en la figura, y el reporte avisa por debajo de 24x24x24.

Códigos como BoltzTraP resuelven esto interpolando las bandas por
funciones estrella sobre una malla mucho más fina que la calculada. Olla-DFT
NO interpola: usa los autovalores tal cual, así que la única salida es
correr el nscf con más puntos. Es más lento, pero no introduce el error de
la interpolación — y sobre todo, no lo esconde.
"""

from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

from qekit.core import provenance, qeout
from qekit.core import style as qstyle
from qekit.core.errors import ErrorDeUso, FaltanDatos

# constantes en unidades SI salvo donde se indique
KB_EV = 8.617333262e-5          # eV/K
E_CHARGE = 1.602176634e-19      # C
HBAR_EVS = 6.582119569e-16      # eV*s
ANG_M = 1e-10


@dataclass
class TransportRun:
    energies: np.ndarray = None       # (nk, nbnd) eV
    weights: np.ndarray = None        # (nk,) pesos de la malla
    velocities: np.ndarray = None     # (nk, nbnd, 3) m/s
    volume: float = None              # A^3
    nelec: float = None
    fermi: float = None               # eV
    grid: tuple = None
    mu: np.ndarray = None             # eV, rejilla de potencial químico
    T: np.ndarray = None              # K
    sigma: np.ndarray = None          # (nT, nmu, 3, 3)  sigma/tau
    seebeck: np.ndarray = None        # (nT, nmu, 3, 3)  V/K
    kappa_e: np.ndarray = None        # (nT, nmu, 3, 3)  kappa/tau
    carriers: np.ndarray = None       # (nT, nmu) cm^-3 (positivo = huecos)
    warnings: list = field(default_factory=list)


def _fd_derivative(E_grid: np.ndarray, cell: np.ndarray, grid: tuple):
    """dE/dk por diferencias finitas sobre la malla uniforme de k.

    E_grid llega con forma (n1, n2, n3, nbnd) y coordenadas de k en
    fraccionarias; la derivada se pasa a cartesianas con la red recíproca.
    Se usa np.gradient con condiciones periódicas (la malla envuelve la
    zona de Brillouin), que es lo que hace correcta la derivada en el
    borde.
    """
    n1, n2, n3, nbnd = E_grid.shape
    recip = 2.0 * np.pi * np.linalg.inv(cell).T      # filas b_i, A^-1
    # envolver una capa para que la derivada en el borde sea periódica
    ext = np.pad(E_grid, ((1, 1), (1, 1), (1, 1), (0, 0)), mode="wrap")
    d1, d2, d3 = np.gradient(ext, 1.0 / n1, 1.0 / n2, 1.0 / n3,
                             axis=(0, 1, 2))
    d1, d2, d3 = d1[1:-1, 1:-1, 1:-1], d2[1:-1, 1:-1, 1:-1], d3[1:-1, 1:-1, 1:-1]
    # dE/dk_cart = sum_i (dE/dfrac_i) * (b_i / 2pi)... cadena completa:
    # k_cart = sum_i frac_i * b_i  =>  dE/dk_cart = B^{-T} dE/dfrac
    dfrac = np.stack([d1, d2, d3], axis=-1)          # (n1,n2,n3,nbnd,3)
    inv_bt = np.linalg.inv(recip.T)                  # (3,3)
    dcart = dfrac @ inv_bt.T                         # eV*A
    # v = (1/hbar) dE/dk : eV*A -> m/s
    return dcart * (ANG_M / HBAR_EVS)


def prepare(atoms, outdir: str = "transporte", pseudo_dir: str = None,
            ecutwfc: float = None, ecutrho: float = None,
            grid=(16, 16, 16), nbnd_factor: float = 2.0,
            insulator: bool = True, kspacing: float = None,
            nspin: int = 1, magnetization: dict = None) -> tuple:
    """Escribe scf.in y un nscf de malla uniforme densa SIN simetría.

    Las velocidades de banda se obtienen derivando E(k) sobre la malla, así
    que hace falta la malla COMPLETA: con la reducción por simetría que hace
    QE por omisión, los puntos no forman una rejilla y la derivada no se
    puede calcular. De ahí nosym/noinv, igual que en el módulo de ópticas.

    `nspin=2` (y `magnetization`, como la devuelve
    inputgen.parse_magnetization) escriben scf y nscf con polarización de
    espín, que es lo que necesita `--spin-resolved` para separar los dos
    canales al recoger.
    """
    from qekit.core import structure as struct_mod
    from qekit.modules import inputgen, sweep

    atoms = struct_mod.primitive(atoms)
    common = sweep.prepare_common(atoms, pseudo_dir, ecutwfc, ecutrho,
                                  insulator)
    out = Path(outdir); out.mkdir(parents=True, exist_ok=True)
    grid = tuple(int(g) for g in grid)
    if magnetization and nspin == 1:
        nspin = 2       # pedir magnetización implica activar el espín
    espin = dict(nspin=nspin, magnetization=magnetization or None)

    grid_scf = sweep.default_grid(atoms, kspacing)
    scf = inputgen.build_pw_input(
        atoms=atoms, pseudos=common["pseudos"], calculation="scf",
        prefix=common["prefix"], pseudo_dir=common["pseudo_dir"],
        ecutwfc=common["ecutwfc"], ecutrho=common["ecutrho"],
        kcard=f"K_POINTS automatic\n  {grid_scf[0]} {grid_scf[1]} "
              f"{grid_scf[2]} 0 0 0\n",
        insulator=insulator, degauss=common["degauss"],
        smearing=common["smearing"], **espin)
    sweep.write_input(out / "scf.in", scf)

    from qekit.modules.inputgen import _estimate_nbnd
    nb = _estimate_nbnd(atoms, common["pseudos"])
    # hacen falta bandas VACÍAS: sin ellas no hay rama tipo n y el Seebeck
    # sale truncado por encima del gap
    nbnd = int(nb * nbnd_factor) if nb else None
    nscf = inputgen.build_pw_input(
        atoms=atoms, pseudos=common["pseudos"], calculation="nscf",
        prefix=common["prefix"], pseudo_dir=common["pseudo_dir"],
        ecutwfc=common["ecutwfc"], ecutrho=common["ecutrho"],
        kcard=f"K_POINTS automatic\n  {grid[0]} {grid[1]} {grid[2]} "
              "0 0 0\n",
        insulator=insulator, degauss=common["degauss"],
        smearing=common["smearing"], nbnd=nbnd, nosym=True, **espin)
    # NO se pone disk_io='nowf' para ahorrar las funciones de onda: en
    # QE 6.x eso hace que punch() retorne ANTES de escribir el XML
    # (PW/src/punch.f90: "IF (io_level < 0) RETURN"), así que el cálculo
    # termina con JOB DONE y no guarda los autovalores. En versiones más
    # nuevas sí escribe el XML — o sea que el mismo input funciona o no
    # según la versión de QE, que es peor que gastar el disco.
    sweep.write_input(out / "nscf.in", nscf)

    rep = ["--- Transporte electrónico (CRTA) ---",
           f"Estructura: {atoms.get_chemical_formula()}  |  "
           f"malla del nscf: {grid[0]}x{grid[1]}x{grid[2]} "
           f"({grid[0]*grid[1]*grid[2]} puntos, sin simetría)",
           f"Bandas: {nbnd if nbnd else 'automáticas'}"
           + ("  |  con polarización de espín (nspin = 2)"
              if nspin == 2 else ""),
           "",
           f"Archivos en '{out.resolve()}': scf.in, nscf.in",
           "Orden: pw.x -in scf.in  ->  pw.x -in nscf.in",
           "",
           "La malla completa (nosym) es obligatoria: las velocidades salen",
           "de derivar E(k) sobre la rejilla, y una malla reducida por",
           "simetría no es una rejilla."]
    warn = sweep.missing_pseudo_warning(common)
    if warn:
        rep.append(warn)
    return grid, "\n".join(rep)


def load(xml_path, spin: int = 0) -> TransportRun:
    """Lee un nscf de malla uniforme y calcula las velocidades de banda."""
    res = qeout.read_xml(xml_path)
    if (res.calculation or "").lower() == "scf":
        raise FaltanDatos(
            f"'{xml_path}' es de un cálculo SCF, no del nscf de malla "
            "densa.\nEl nscf no llegó a escribir su XML y quedó el del "
            "scf. Suele pasar por dos\nmotivos: el nscf falló (mira "
            "nscf.out), o corrió con disk_io='nowf'/'none',\nque en QE 6.x "
            "impide escribir el XML aunque el cálculo termine bien.")
    run = TransportRun(volume=res.volume, nelec=res.nelec, fermi=res.fermi)
    E = res.eigenvalues[spin]                        # (nk, nbnd)
    # Llevar las fraccionarias a [0,1) CON TOLERANCIA: un cero que el XML
    # devuelve como 0.9999999 se contaria como un valor distinto de 0.0 y
    # la malla saldria de 9x9x9 en vez de 8x8x8.
    kfrac = np.asarray(res.kpoints_frac, dtype=float)
    kfrac = kfrac - np.floor(kfrac + 1e-6)
    kfrac = np.round(kfrac, 6)
    kfrac[np.abs(kfrac - 1.0) < 1e-5] = 0.0

    # reconstruir la malla uniforme a partir de los puntos k
    grid = []
    for d in range(3):
        vals = np.unique(kfrac[:, d])
        grid.append(len(vals))
    n1, n2, n3 = grid
    run.grid = (n1, n2, n3)
    if n1 * n2 * n3 != len(kfrac):
        raise FaltanDatos(
            f"los {len(kfrac)} puntos k no forman una malla uniforme "
            f"{n1}x{n2}x{n3}. El transporte necesita un nscf con "
            "K_POINTS automatic y nosym/noinv, no un camino de bandas ni "
            "una malla reducida por simetría."
        )
    idx = np.lexsort((kfrac[:, 2], kfrac[:, 1], kfrac[:, 0]))
    E_grid = E[idx].reshape(n1, n2, n3, -1)
    run.energies = E
    run.velocities = _fd_derivative(E_grid, res.cell, run.grid).reshape(
        -1, E.shape[1], 3)
    # reordenar las energías igual que las velocidades
    run.energies = E_grid.reshape(-1, E.shape[1])
    run.weights = np.full(len(run.energies), 1.0 / len(run.energies))
    ntot = n1 * n2 * n3
    if min(run.grid) < 24 or ntot < 12000:
        run.warnings.append(
            f"malla {n1}x{n2}x{n3} ({ntot} puntos): INSUFICIENTE. A 300 K "
            "la ventana de -df/dE mide unos 25 meV, y con una malla así "
            "solo caen dentro unos pocos estados sueltos: sigma sale en "
            "picos aislados en vez de una curva suave, y S hereda el "
            "ruido. Se nota a simple vista en la figura. Para transporte "
            "hacen falta mallas de 24x24x24 en adelante — bastante más "
            "densas que las de un scf.")
    return run


def compute(run: TransportRun, T=None, mu=None, mu_span: float = 1.0,
            nmu: int = 201) -> TransportRun:
    """Coeficientes de transporte en CRTA sobre rejillas de mu y T."""
    if T is None:
        T = np.array([300.0])
    T = np.atleast_1d(np.asarray(T, dtype=float))
    ef = run.fermi if run.fermi is not None else float(np.median(run.energies))
    if mu is None:
        mu = np.linspace(ef - mu_span, ef + mu_span, nmu)
    mu = np.atleast_1d(np.asarray(mu, dtype=float))

    E = run.energies                                  # (nk, nbnd) eV
    v = run.velocities                                # (nk, nbnd, 3) m/s
    w = run.weights[:, None]                          # (nk, 1)
    vol_m3 = run.volume * (ANG_M ** 3)

    # v (x) v para cada estado: (nk, nbnd, 3, 3)
    vv = v[..., :, None] * v[..., None, :]

    nT, nmu_ = len(T), len(mu)
    sigma = np.zeros((nT, nmu_, 3, 3))
    seebeck = np.zeros((nT, nmu_, 3, 3))
    kappa = np.zeros((nT, nmu_, 3, 3))
    carriers = np.zeros((nT, nmu_))

    for it, t in enumerate(T):
        kt = KB_EV * t
        for im, m in enumerate(mu):
            x = (E - m) / kt
            # -df/dE = 1/(4kT) sech^2(x/2); se usa la forma estable
            with np.errstate(over="ignore"):
                sech2 = 1.0 / (np.cosh(np.clip(x, -300, 300) / 2.0) ** 2)
            mdf = sech2 / (4.0 * kt)                   # 1/eV
            peso = (w * mdf)[..., None, None]          # (nk, nbnd, 1, 1)

            s0 = np.sum(peso * vv, axis=(0, 1))        # m^2/s^2 /eV
            s1 = np.sum(peso * vv * (E - m)[..., None, None], axis=(0, 1))
            s2 = np.sum(peso * vv * ((E - m) ** 2)[..., None, None],
                        axis=(0, 1))

            # sigma/tau  [S/(m*s)] = e^2/(V) * s0   (s0 en m^2/s^2/eV;
            # el eV^-1 y la carga se combinan: e^2/eV = e [C/V])
            sig = (E_CHARGE / vol_m3) * s0
            sigma[it, im] = sig
            # S = -(1/(eT)) * s1/s0   -> V/K   (s1/s0 sale en eV)
            # El signo MENOS viene de que la carga del portador es -e. Sin
            # él, un semiconductor tipo p daria un Seebeck negativo, que es
            # justo al reves de lo que se mide.
            s0_inv = np.linalg.pinv(s0)
            S_mat = -(s1 @ s0_inv) / t
            seebeck[it, im] = S_mat
            # kappa_e a corriente NULA: al integrando (E-mu)^2 hay que
            # restarle el termino S^2*sigma*T, que es el calor arrastrado
            # por la corriente termoelectrica. Sin restarlo se reporta
            # kappa_0, que es otra cantidad.
            kappa0 = (E_CHARGE / (vol_m3 * t)) * s2
            kappa[it, im] = kappa0 - (S_mat @ S_mat) @ sig * t

            # portadores: electrones por celda por encima de mu menos huecos
            f = 1.0 / (1.0 + np.exp(np.clip(x, -300, 300)))
            n_e = 2.0 * float(np.sum(run.weights[:, None] * f))
            carriers[it, im] = (run.nelec - n_e) / (run.volume * 1e-24)

    run.T, run.mu = T, mu
    run.sigma, run.seebeck, run.kappa_e = sigma, seebeck, kappa
    run.carriers = carriers
    return run


def crossing_bands(run: TransportRun, tol: float = 1e-6) -> list:
    """Indices de las bandas que cruzan el nivel de Fermi."""
    ef = run.fermi
    if ef is None:
        return []
    E = run.energies
    return [b for b in range(E.shape[1])
            if E[:, b].min() < ef - tol < ef + tol < E[:, b].max()]


def export_bxsf(run: TransportRun, cell, path, bands=None) -> str:
    """Escribe la superficie de Fermi en formato BXSF (XCrySDen/FermiSurfer).

    Detalle que estropea el archivo si se ignora: BXSF espera la rejilla
    con el punto final PERIODICO incluido, o sea (n+1) puntos por
    direccion, donde el ultimo plano repite el primero. Escribir solo los n
    puntos del calculo da un archivo que abre, pero con la superficie
    cortada por costuras en los bordes de la zona.
    """
    if run.fermi is None:
        raise ErrorDeUso("no hay nivel de Fermi: la superficie de Fermi solo "
                         "tiene sentido en un metal")
    idx = crossing_bands(run) if bands is None else list(bands)
    if not idx:
        raise ErrorDeUso(
            "ninguna banda cruza el nivel de Fermi: el sistema es aislante o "
            "semiconductor y no tiene superficie de Fermi")

    n1, n2, n3 = run.grid
    E = run.energies.reshape(n1, n2, n3, -1)
    recip = 2.0 * np.pi * np.linalg.inv(np.asarray(cell)).T

    path = Path(path)
    with open(path, "w") as f:
        f.write("BEGIN_INFO\n")
        f.write("  Fermi Energy: %.8f\n" % run.fermi)
        f.write("END_INFO\n\n")
        f.write("BEGIN_BLOCK_BANDGRID_3D\n")
        f.write(" superficie_de_fermi_QEkit\n")
        f.write(" BEGIN_BANDGRID_3D\n")
        f.write("  %d\n" % len(idx))
        f.write("  %d %d %d\n" % (n1 + 1, n2 + 1, n3 + 1))
        f.write("  0.0 0.0 0.0\n")
        for v in recip:
            f.write("  %.8f %.8f %.8f\n" % (v[0], v[1], v[2]))
        for b in idx:
            f.write("  BAND:  %d\n" % (b + 1))
            ext = np.pad(E[:, :, :, b], ((0, 1), (0, 1), (0, 1)), mode="wrap")
            vals = ext.ravel(order="C")
            for k in range(0, len(vals), 6):
                f.write("  " + " ".join("%14.8f" % v
                                        for v in vals[k:k + 6]) + "\n")
        f.write(" END_BANDGRID_3D\n")
        f.write("END_BLOCK_BANDGRID_3D\n")
    return str(path)


def power_factor(run: TransportRun) -> np.ndarray:
    """S^2 * sigma/tau, en W/(m*K^2*s)."""
    s = np.trace(run.seebeck, axis1=2, axis2=3) / 3.0
    sig = np.trace(run.sigma, axis1=2, axis2=3) / 3.0
    return (s ** 2) * sig


def at_temperature(run: TransportRun, t: float) -> int:
    return int(np.argmin(np.abs(run.T - t)))


def report(run: TransportRun, t: float = 300.0) -> str:
    it = at_temperature(run, t)
    s_iso = np.trace(run.seebeck[it], axis1=1, axis2=2) / 3.0 * 1e6   # uV/K
    np.trace(run.sigma[it], axis1=1, axis2=2) / 3.0
    pf = power_factor(run)[it]
    ef = run.fermi
    i_ef = int(np.argmin(np.abs(run.mu - ef)))

    lines = ["--- Transporte electrónico (CRTA) ---",
             f"Malla k: {run.grid[0]}x{run.grid[1]}x{run.grid[2]}  |  "
             f"volumen {run.volume:.2f} Å³  |  T = {run.T[it]:.0f} K"]
    if ef is not None:
        lines.append(f"En el nivel de Fermi ({ef:.3f} eV): "
                     f"S = {s_iso[i_ef]:+.1f} µV/K")
    j_p = int(np.argmax(np.where(run.carriers[it] > 0, s_iso, -np.inf)))
    j_n = int(np.argmin(np.where(run.carriers[it] < 0, s_iso, np.inf)))
    lines += ["",
              "Mejor coeficiente Seebeck en la ventana explorada:",
              f"  tipo p: S = {s_iso[j_p]:+7.1f} µV/K  en µ − E_F = "
              f"{run.mu[j_p] - (ef or 0):+.3f} eV  "
              f"(n = {run.carriers[it][j_p]:.2e} cm⁻³)",
              f"  tipo n: S = {s_iso[j_n]:+7.1f} µV/K  en µ − E_F = "
              f"{run.mu[j_n] - (ef or 0):+.3f} eV  "
              f"(n = {run.carriers[it][j_n]:.2e} cm⁻³)",
              "",
              f"Máximo del factor de potencia: {np.max(pf):.3e} "
              "W/(m·K²·s), es decir PF/τ",
              "",
              "IMPORTANTE: en CRTA el tiempo de relajación τ se cancela en S,",
              "que por tanto es una predicción real. En σ y κ_e NO se cancela:",
              "van como σ/τ y κ_e/τ. Para dar σ en S/m hace falta un τ que",
              "venga de un ajuste a una medida o de un cálculo de",
              "electrón-fonón — Olla-DFT no lo inventa."]
    for w in run.warnings:
        lines.append(f"\nAVISO: {w}")
    return "\n".join(lines)


def export(run: TransportRun, outdir: str = ".", t: float = 300.0) -> list:
    out = Path(outdir); out.mkdir(parents=True, exist_ok=True)
    it = at_temperature(run, t)
    s_iso = np.trace(run.seebeck[it], axis1=1, axis2=2) / 3.0 * 1e6
    sig_iso = np.trace(run.sigma[it], axis1=1, axis2=2) / 3.0
    k_iso = np.trace(run.kappa_e[it], axis1=1, axis2=2) / 3.0
    pf = power_factor(run)[it]
    f = out / "TRANSPORTE.dat"
    cab = provenance.header_plain(
        "transporte (CRTA)",
        {"T_K": float(run.T[it]),
         "malla": "x".join(map(str, run.grid))},
        titulo="Transporte electronico en CRTA")
    np.savetxt(f, np.column_stack([run.mu, s_iso, sig_iso, k_iso, pf,
                                   run.carriers[it]]),
               fmt="%16.6e", comments="# ",
               header=cab + "\n"
               f"{'mu(eV)':>14s} {'S(uV/K)':>16s} {'sigma/tau':>16s} "
               f"{'kappa_e/tau':>16s} {'PF/tau':>16s} {'n(cm^-3)':>16s}")
    return [str(f)]


def plot(run: TransportRun, outfile: str = "transporte", formats="pdf,png",
         theme: str = None, size: str = None, family: str = None,
         background: str = None, palette=None, usetex: bool = None,
         width="double", journal: str = "generic", aspect: float = 0.38,
         mono: bool = False, dpi: int = None, t: float = 300.0) -> list:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    st = qstyle.apply(theme, size=size, family=family, background=background,
                      palette=palette, usetex=usetex, mono=mono)
    fig = plt.figure(figsize=qstyle.figure_size(width, journal, aspect),
                     layout="constrained")
    axes = [qstyle.finish_axes(fig.add_subplot(1, 3, i + 1)) for i in range(3)]
    c = qstyle.palette(3, mono=mono)
    it = at_temperature(run, t)
    ef = run.fermi or 0.0
    x = run.mu - ef

    s_iso = np.trace(run.seebeck[it], axis1=1, axis2=2) / 3.0 * 1e6
    sig_iso = np.trace(run.sigma[it], axis1=1, axis2=2) / 3.0
    pf = power_factor(run)[it]

    for ax, y, lab, col in (
            (axes[0], s_iso, r"$S$ ($\mu$V/K)", c[0]),
            (axes[1], sig_iso, r"$\sigma/\tau$ (S m$^{-1}$s$^{-1}$)", c[1]),
            (axes[2], pf, r"$S^2\sigma/\tau$ (W m$^{-1}$K$^{-2}$s$^{-1}$)",
             c[2])):
        ax.plot(x, y, color=col, lw=st["line"])
        ax.axvline(0.0, color=qstyle.INK_FAINT, lw=st["axis_line"],
                   dashes=[3.5, 2.0])
        ax.set_xlabel(qstyle.tex_safe("mu - E_F (eV)").replace(
            "mu", r"$\mu$").replace("E_F", r"$E_\mathrm{F}$"))
        ax.set_ylabel(lab)
        ax.set_xlim(x.min(), x.max())
    axes[0].axhline(0.0, color=qstyle.INK_FAINT, lw=st["axis_line"])
    for ax, lab in zip(axes, ("(a)", "(b)", "(c)")):
        qstyle.panel_label(ax, lab)
    written = qstyle.save(fig, outfile, formats, dpi=dpi,
                          modulo="transporte (CRTA)")
    plt.close(fig)
    return written


# ----------------------------------------------------------------------
# Wiedemann-Franz y transporte por canal de espín
# ----------------------------------------------------------------------
# Número de Lorenz de Sommerfeld, (pi^2/3)(k_B/e)^2, en W·Ω/K²
L0_SOMMERFELD = 2.44e-8


def lorenz(run: TransportRun, it: int = 0) -> np.ndarray:
    """L = κ_e / (σ T), isótropo, en W·Ω/K².

    Es de las pocas magnitudes que la CRTA da en ABSOLUTO: κ_e y σ llevan
    los dos el mismo τ y se cancela en el cociente, igual que en el Seebeck.
    Se compara con L₀ = 2.44e-8: por debajo hay transporte bipolar o
    dispersión inelástica, y por encima suele haber contribución de fonones
    mal separada.
    """
    sig = np.trace(run.sigma[it], axis1=1, axis2=2) / 3.0
    kap = np.trace(run.kappa_e[it], axis1=1, axis2=2) / 3.0
    T = float(run.T[it])
    with np.errstate(divide="ignore", invalid="ignore"):
        L = np.where(np.abs(sig) > 0, kap / (sig * T), np.nan)
    return L


def cancelacion(run: TransportRun, it: int = 0) -> np.ndarray:
    """|κ_e| / κ⁰: cuánto sobrevive a la resta κ_e = κ⁰ − S²σT.

    Cuando el potencial químico cae dentro del gap, S se dispara y S²σT se
    acerca tanto a κ⁰ que la diferencia pierde todas sus cifras
    significativas. El número que sale entonces no es pequeño por física:
    es ruido de coma flotante. Sin esta comprobación el módulo devolvía
    L/L₀ = 0.001 en silicio con toda naturalidad.
    """
    sig = np.trace(run.sigma[it], axis1=1, axis2=2) / 3.0
    kap = np.trace(run.kappa_e[it], axis1=1, axis2=2) / 3.0
    See = np.trace(run.seebeck[it], axis1=1, axis2=2) / 3.0
    T = float(run.T[it])
    k0 = kap + See ** 2 * sig * T
    with np.errstate(divide="ignore", invalid="ignore"):
        return np.where(np.abs(k0) > 0, np.abs(kap) / np.abs(k0), np.nan)


def report_lorenz(run: TransportRun, t: float = 300.0) -> str:
    it = at_temperature(run, t)
    L = lorenz(run, it)
    canc = cancelacion(run, it)
    ef = run.fermi if run.fermi is not None else float(np.median(run.mu))
    i_ef = int(np.argmin(np.abs(run.mu - ef)))
    L_ef, c_ef = float(L[i_ef]), float(canc[i_ef])
    lineas = [f"--- Número de Lorenz (T = {run.T[it]:.0f} K) ---",
              f"L(E_F) = {L_ef:.3e} W·Ω/K²   "
              f"(L₀ de Sommerfeld = {L0_SOMMERFELD:.2e})",
              f"L/L₀   = {L_ef / L0_SOMMERFELD:.3f}"]

    if c_ef < 1e-3:
        lineas += [
            "",
            "NO TE FÍES DE ESTE NÚMERO. κ_e sale de restar S²σT a κ⁰, y "
            "aquí de esa resta\n"
            f"  sobrevive el {c_ef * 100:.4f} %: las dos cantidades coinciden "
            f"en casi todas sus" + "\n"
            "  cifras y lo que queda es ruido de coma flotante, no física. "
            "Pasa siempre que" + "\n"
            "  µ cae dentro del gap, donde S se dispara.",
            "  Mira L donde SÍ hay portadores: en la ventana de µ que sale "
            "abajo, o dopando."]
    else:
        r = L_ef / L0_SOMMERFELD
        if r > 1.15:
            lineas.append(
                "  Por encima de L₀: es la firma del transporte BIPOLAR. Con "
                "electrones y huecos" + "\n"
                "  a la vez, los dos llevan calor en el mismo sentido y carga "
                "en sentidos" + "\n"
                "  opuestos, así que κ_e crece sin que crezca σ.")
        elif r < 0.85:
            lineas.append(
                "  Por debajo de L₀: es lo normal cuando el gas de portadores "
                "NO es degenerado." + "\n"
                "  En el límite no degenerado y con τ constante, L tiende a "
                "2.5·(k_B/e)² = 1.86e-8," + "\n"
                "  o sea 0.76·L₀. Wiedemann-Franz solo vale en el límite "
                "metálico.")
        else:
            lineas.append(
                "  Se cumple Wiedemann-Franz dentro de un 15 %: el gas de "
                "portadores es" + "\n"
                "  degenerado y la dispersión elástica, que es justo lo que "
                "supone la CRTA.")

    # "fiable" pide que sobreviva al menos el 10 % de la resta. Con el 1 %
    # todavía entran puntos cuyo L es basura, y la ventana que se anuncia
    # sale más ancha de lo que es.
    buenos = np.where(np.isfinite(L) & (canc > 0.10))[0]
    if len(buenos) >= 3:
        Lb = L[buenos] / L0_SOMMERFELD
        p10, med, p90 = np.percentile(Lb, [10, 50, 90])
        lineas += ["",
                   f"Donde la resta conserva más del 10 % "
                   f"({len(buenos)} de {len(L)} puntos de µ):",
                   f"  L/L₀ mediana {med:.2f}, con el 80 % central entre "
                   f"{p10:.2f} y {p90:.2f}"]
    elif c_ef < 1e-3:
        lineas += ["",
                   "En NINGÚN punto de la ventana de µ sobrevive un 10 % de "
                   "la resta: este" + "\n"
                   "  cálculo no puede dar el número de Lorenz. Amplía "
                   "--mu-span para llegar a" + "\n"
                   "  zonas con portadores de verdad."]
    lineas.append("\nEste número NO depende de τ: κ_e y σ lo llevan los dos "
                  "y se cancela. Es de\n  las pocas cosas que la CRTA da en "
                  "absoluto, junto con el Seebeck.")
    return "\n".join(lineas)


@dataclass
class TransporteEspin:
    """Los dos canales de espín y su combinación en el modelo de dos corrientes."""

    up: TransportRun = None
    dw: TransportRun = None
    it: int = 0

    def _iso(self, run, campo):
        return np.trace(getattr(run, campo)[self.it], axis1=1, axis2=2) / 3.0

    @property
    def sigma_up(self):
        return self._iso(self.up, "sigma")

    @property
    def sigma_dw(self):
        return self._iso(self.dw, "sigma")

    @property
    def seebeck_up(self):
        return self._iso(self.up, "seebeck")

    @property
    def seebeck_dw(self):
        return self._iso(self.dw, "seebeck")

    @property
    def sigma_total(self):
        """Los dos canales conducen EN PARALELO: las conductancias se suman."""
        return self.sigma_up + self.sigma_dw

    @property
    def seebeck_total(self):
        """S de la mezcla: media de los dos, pesada por su conductancia.

        No es la media aritmética. Un canal que no conduce no aporta
        termopotencia por muy grande que sea su S, y promediar sin pesar es
        el error típico al juntar los dos canales a mano.
        """
        s = self.sigma_total
        with np.errstate(divide="ignore", invalid="ignore"):
            return np.where(np.abs(s) > 0,
                            (self.seebeck_up * self.sigma_up
                             + self.seebeck_dw * self.sigma_dw) / s, 0.0)

    @property
    def polarizacion(self):
        """(σ↑ − σ↓)/(σ↑ + σ↓): +1 es medio metal de espín arriba."""
        s = self.sigma_total
        with np.errstate(divide="ignore", invalid="ignore"):
            return np.where(np.abs(s) > 0,
                            (self.sigma_up - self.sigma_dw) / s, 0.0)

    @property
    def seebeck_de_espin(self):
        """S↑ − S↓, la termopotencia de espín."""
        return self.seebeck_up - self.seebeck_dw


def report_espin(te: TransporteEspin, t: float = 300.0) -> str:
    up = te.up
    ef = up.fermi if up.fermi is not None else float(np.median(up.mu))
    i = int(np.argmin(np.abs(up.mu - ef)))
    a_uv = 1e6
    L = [f"--- Transporte por canal de espín (T = {up.T[te.it]:.0f} K) ---",
         "",
         f"  {'':22s} {'espín ↑':>14s} {'espín ↓':>14s} {'mezcla':>14s}",
         "  " + "-" * 68,
         f"  {'S (µV/K)':22s} {te.seebeck_up[i] * a_uv:>14.2f} "
         f"{te.seebeck_dw[i] * a_uv:>14.2f} "
         f"{te.seebeck_total[i] * a_uv:>14.2f}",
         f"  {'σ/τ (S/(m·s))':22s} {te.sigma_up[i]:>14.3e} "
         f"{te.sigma_dw[i]:>14.3e} {te.sigma_total[i]:>14.3e}",
         "",
         f"Polarización de espín de la conductividad: "
         f"P = {te.polarizacion[i]:+.4f}",
         f"Termopotencia de espín  S↑ − S↓ = "
         f"{te.seebeck_de_espin[i] * a_uv:+.2f} µV/K"]
    p = abs(float(te.polarizacion[i]))
    if p > 0.95:
        L.append("  |P| > 0.95: prácticamente medio metal, solo conduce un "
                 "canal de espín.")
    elif p < 0.05:
        L.append("  P ≈ 0: los dos canales conducen igual; separarlos no "
                 "aporta nada aquí.")
    L += ["",
          "El S de la mezcla es la media de los dos canales PESADA POR SU "
          "CONDUCTANCIA,\n  no la media aritmética: un canal que no conduce "
          "no aporta termopotencia por\n  mucha que tenga.",
          "",
          "Aviso: esto supone que los dos canales son independientes (modelo "
          "de dos\n  corrientes). Vale mientras la dispersión con inversión "
          "de espín sea lenta\n  frente a la normal, que es lo habitual muy "
          "por debajo de la temperatura de\n  Curie y deja de valer cerca de "
          "ella."]
    return "\n".join(L)
