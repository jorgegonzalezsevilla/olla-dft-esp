# Olla-DFT — command-line toolkit for Quantum ESPRESSO
# Copyright (C) 2026 Jorge Enrique González Sevilla
# SPDX-License-Identifier: AGPL-3.0-or-later
# This program is free software: you can redistribute it and/or modify it
# under the terms of the GNU Affero General Public License as published by the Free
# Software Foundation, either version 3 of the License, or (at your option)
# any later version. See the LICENSE file for details.

"""Masa efectiva por ajuste parabólico de las bandas.

Cerca de un extremo, E(k) ≈ E₀ + ħ²(k−k₀)²/(2m*), así que la masa sale de
la curvatura:

    m*/mₑ = (ħ²/mₑ) / (d²E/dk²)     con ħ²/mₑ = 7.6199682 eV·Å²

El comando `olla-dft effmass ESTRUCTURA --bands-dir CARPETA` hace dos cosas
seguidas, porque hacen falta las dos:

1. **Ajuste rápido sobre el camino de bandas** (`from_bands`): usa el
   cálculo de bandas que ya tienes en `--bands-dir`. Es inmediato, pero
   solo da masas EN LAS DIRECCIONES DEL CAMINO, y con la densidad de
   puntos de un camino normal (~30 por tramo) el ajuste se hace con muy
   pocos puntos. `--window` (semiancho en Å⁻¹, por omisión
   WINDOW_DEFAULT) y `--min-points` controlan cuántos entran. Sirve para
   una primera estimación y para localizar los extremos.

2. **Preparar el cálculo dedicado** (`prepare`): escribe en `--outdir`
   un scf.in y un masa.in con líneas finas de puntos k que cruzan el
   extremo (`--half-width` Å⁻¹ a cada lado, `--points` por línea). Para un
   valle fuera de Γ (el Δ del silicio, por ejemplo) usa la dirección radial
   Γ→k₀ —la LONGITUDINAL— y dos perpendiculares a ella —las
   TRANSVERSALES—, que es la descomposición con la que se reportan m*_l y
   m*_t. Con `--run` se corre en el momento; si no, a mano:
   `pw.x -in scf.in` y después `pw.x -in masa.in`.

Después, `olla-dft effmass ESTRUCTURA --collect -o CARPETA` lee el cálculo
fino (`collect_fine`) y da los números publicables.

Lo que el reporte deja explícito, porque es donde se cuela el error:
- el número de puntos y el tramo en Å⁻¹ de cada ajuste, y su R²: una
  parábola ajustada a 3 puntos siempre da R² = 1 y no significa nada;
- si el extremo cae en un punto de alta simetría, cada lado es una
  dirección distinta y se reportan por separado, sin promediarlos;
- las bandas degeneradas con el extremo se ajustan por separado (en el
  silicio, hueco pesado y hueco ligero salen de bandas distintas).
"""

import json
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

from qekit.core import provenance, qeout
from qekit.modules import bands as bands_mod
from qekit.core.errors import ErrorDeUso, FaltanDatos

# ħ²/mₑ en eV·Å²  (= 2 × 3.80998212 eV·Å²)
HBAR2_OVER_ME = 7.6199682

DEGEN_TOL = 0.05        # eV: bandas "degeneradas" con el extremo
# más allá de ~0.1 Å⁻¹ la banda ya no es parabólica en un
# semiconductor típico y el ajuste deja de significar una masa.
# PARABOLIC_MAX es la EXTENSIÓN TOTAL del tramo ajustado (de k_min a k_max),
# que es lo que guarda MassFit.window y lo que se compara en el aviso.
PARABOLIC_MAX = 0.12    # Å⁻¹, extensión total del ajuste
# La ventana de `from_bands` y de --window es un SEMIANCHO: se toman los
# puntos con |k − k₀| ≤ window, así que el tramo ajustado mide hasta
# 2·window. Por eso la ventana por omisión es la mitad del límite
# parabólico: con ella un ajuste centrado en el extremo queda justo dentro
# del régimen y el aviso solo salta cuando de verdad hubo que ensanchar
# (min_pts) porque el camino no tiene puntos más finos. Antes valía 0.15
# (0.30 de tramo) contra un límite de 0.12, y el aviso saltaba siempre.
WINDOW_DEFAULT = PARABOLIC_MAX / 2    # Å⁻¹, semiancho: ±0.06
# holgura numérica al comparar el tramo con el límite (los k del camino
# se acumulan en coma flotante)
_TOL_VENTANA = 1e-6


@dataclass
class MassFit:
    carrier: str = ""            # "electrón" | "hueco"
    band: int = 0                # índice de banda (base 0)
    kindex: int = 0
    k_label: str = ""
    direction: str = ""
    mass: float = None           # m*/mₑ (con signo del ajuste)
    energy: float = None         # eV
    npts: int = 0
    window: float = 0.0          # Å⁻¹
    r2: float = None
    warning: str = ""


@dataclass
class EffMassRun:
    fits: list = field(default_factory=list)
    is_metal: bool = False
    vbm: float = None
    cbm: float = None
    source: str = ""


# ----------------------------------------------------------------------
def _r2(x, y, coef) -> float:
    pred = np.polyval(coef, x)
    ss_res = float(np.sum((y - pred) ** 2))
    ss_tot = float(np.sum((y - y.mean()) ** 2))
    return 1.0 - ss_res / ss_tot if ss_tot > 0 else float("nan")


def _mass_from_quadratic(a: float) -> float:
    """m*/mₑ a partir del coeficiente cuadrático a de E = a·k² + b·k + c.

    d²E/dk² = 2a, luego m* = ħ²/(2a). El signo se conserva: negativo
    significa curvatura hacia abajo, es decir un hueco.
    """
    if a == 0 or not np.isfinite(a):
        return float("nan")
    return HBAR2_OVER_ME / (2.0 * a)


def _segment_bounds(i: int, nk: int, stops: set) -> tuple:
    """Hasta dónde se puede extender el ajuste sin cambiar de dirección.

    Se frena en las discontinuidades del camino y en los puntos de alta
    simetría: pasado uno de ellos, la recta del camino gira y los puntos
    ya no pertenecen a la misma dirección.
    """
    lo = i
    while lo - 1 >= 0:
        lo -= 1
        if lo in stops:          # el punto de giro sí entra, lo de más allá no
            break
    hi = i
    while hi + 1 < nk:
        hi += 1
        if hi in stops:
            break
    return lo, hi


def _collect_window(k, i0, lo, hi, window, min_pts):
    """Índices a un lado o a ambos lados del extremo dentro de la ventana."""
    idx = [j for j in range(lo, hi + 1) if abs(k[j] - k[i0]) <= window]
    if len(idx) < min_pts:                    # ensanchar hasta juntar puntos
        order = sorted(range(lo, hi + 1), key=lambda j: abs(k[j] - k[i0]))
        idx = sorted(order[:min(min_pts, hi - lo + 1)])
    return idx


def from_bands(bs: bands_mod.BandStructure, spin: int = 0,
               window: float = WINDOW_DEFAULT, min_pts: int = 7,
               degen_tol: float = DEGEN_TOL) -> EffMassRun:
    """Ajusta masas efectivas sobre un cálculo de bandas ya hecho.

    `window` es el semiancho en Å⁻¹ alrededor del extremo (se ajustan los
    puntos con |k − k₀| ≤ window); si no juntan `min_pts` puntos se
    ensancha hasta juntarlos. El tramo ajustado (MassFit.window) se compara
    con PARABOLIC_MAX y, si lo supera, el ajuste lleva aviso.
    """
    if window is None:
        window = WINDOW_DEFAULT
    info = bands_mod.analyze_gap(bs, spin=spin)
    run = EffMassRun(is_metal=info.is_metal, vbm=info.vbm, cbm=info.cbm,
                     source="camino de bandas")
    if info.is_metal:
        return run

    res = bs.result
    E = res.eigenvalues[spin]            # (nk, nbnd)
    k = np.asarray(bs.kdist, dtype=float)
    nk = len(k)
    stops = set(int(b) for b in bs.breaks)
    stops |= set(int(i) for i, _ in bs.labels)

    tareas = []
    if info.vbm_band is not None:
        # huecos: la banda del VBM y las degeneradas con ella EN ESE k
        kv = info.vbm_kindex
        for b in range(info.vbm_band, -1, -1):
            if info.vbm - E[kv, b] > degen_tol:
                break
            tareas.append(("hueco", b, kv))
        cb = info.vbm_band + 1
        kc = info.cbm_kindex
        for b in range(cb, res.nbnd):
            if E[kc, b] - info.cbm > degen_tol:
                break
            tareas.append(("electrón", b, kc))

    for carrier, band, kext in tareas:
        y_all = E[:, band]
        # reubicar el extremo en esta banda (puede no ser el mismo k)
        if carrier == "hueco":
            i0 = int(np.argmax(y_all))
        else:
            i0 = int(np.argmin(y_all))
        if abs(y_all[i0] - y_all[kext]) > degen_tol:
            i0 = kext                     # banda degenerada con extremo lejano
        lo, hi = _segment_bounds(i0, nk, stops)
        interior = (i0 - lo >= 2) and (hi - i0 >= 2)

        lados = [("ambos", lo, hi)] if interior else []
        if not interior:
            if i0 - lo >= 2:
                lados.append(("izq", lo, i0))
            if hi - i0 >= 2:
                lados.append(("der", i0, hi))

        for etiqueta, a_lo, a_hi in lados:
            idx = _collect_window(k, i0, a_lo, a_hi, window, min_pts)
            if len(idx) < 3:
                continue
            x = k[idx] - k[i0]
            y = y_all[idx]
            coef = np.polyfit(x, y, 2)
            m = _mass_from_quadratic(float(coef[0]))
            fit = MassFit(
                carrier=carrier, band=band, kindex=i0,
                k_label=bands_mod._label_for(i0, bs.labels, res.kpoints_frac),
                direction=_direction_text(bs, i0, idx, etiqueta),
                mass=m, energy=float(y_all[i0]),
                npts=len(idx), window=float(x.max() - x.min()),
                r2=_r2(x, y, coef),
            )
            if len(idx) < 5:
                fit.warning = ("solo %d puntos: el ajuste no es confiable; "
                               "haz el cálculo dedicado (effmass sin "
                               "--collect y luego --collect)" % len(idx))
            elif fit.window > PARABOLIC_MAX + _TOL_VENTANA:
                # las dos son extensiones totales del tramo: se compara
                # lo mismo con lo mismo
                fit.warning = (
                    "tramo ajustado de %.3f Å⁻¹ (límite parabólico %.2f): "
                    "el camino no tiene puntos más finos; haz el cálculo "
                    "dedicado (effmass sin --collect y luego --collect)"
                    % (fit.window, PARABOLIC_MAX))
            run.fits.append(fit)
    return run


def _direction_text(bs, i0, idx, etiqueta) -> str:
    """Descripción legible de la dirección del ajuste en el espacio k."""
    res = bs.result
    kc = res.kpoints_cart
    if kc is None or len(idx) < 2:
        return etiqueta
    v = kc[idx[-1]] - kc[idx[0]]
    n = np.linalg.norm(v)
    if n < 1e-8:
        return etiqueta
    v = v / n
    comp = "  ".join(f"{c:+.3f}" for c in v)
    return f"[{comp}] ({etiqueta})"


# ----------------------------------------------------------------------
# Camino fino dedicado
# ----------------------------------------------------------------------
def valley_directions(k_ext_cart: np.ndarray) -> list:
    """(nombre, vector unitario) longitudinal y dos transversales.

    Longitudinal = dirección radial Γ→k₀, que es el eje del valle. Si el
    extremo está EN Γ no hay eje privilegiado y se devuelven los tres ejes
    cartesianos.
    """
    v = np.asarray(k_ext_cart, dtype=float)
    n = np.linalg.norm(v)
    if n < 1e-6:
        # Extremo EN Γ: x, y, z serían equivalentes por simetría en un
        # cristal cúbico y darían tres veces el mismo número. Las bandas de
        # valencia están fuertemente alabeadas ("warping"): en el silicio el
        # hueco pesado vale 0.28 mₑ a lo largo de [100] y 0.72 a lo largo de
        # [111]. Por eso se muestrean las tres direcciones cristalográficas
        # principales, que son las que capturan esa anisotropía.
        r2, r3 = np.sqrt(2.0), np.sqrt(3.0)
        return [("[100]", np.array([1.0, 0.0, 0.0])),
                ("[110]", np.array([1.0, 1.0, 0.0]) / r2),
                ("[111]", np.array([1.0, 1.0, 1.0]) / r3)]
    e1 = v / n
    # dos perpendiculares cualesquiera, ortonormales entre sí
    tmp = np.array([1.0, 0, 0])
    if abs(np.dot(tmp, e1)) > 0.9:
        tmp = np.array([0, 1.0, 0])
    e2 = np.cross(e1, tmp); e2 /= np.linalg.norm(e2)
    e3 = np.cross(e1, e2)
    return [("longitudinal", e1), ("transversal 1", e2), ("transversal 2", e3)]


def prepare(atoms, bs: bands_mod.BandStructure, outdir: str = "masa_efectiva",
            pseudo_dir: str = None, ecutwfc: float = None,
            ecutrho: float = None, half_width: float = 0.06,
            npts: int = 21, spin: int = 0, insulator: bool = True) -> tuple:
    """Escribe un nscf 'bands' con líneas finas que cruzan VBM y CBM.

    half_width en Å⁻¹ a cada lado del extremo; npts por línea (impar, para
    que el extremo caiga justo en un punto muestreado).
    """
    from qekit.core import structure as struct_mod
    from qekit.modules import inputgen, sweep

    if npts % 2 == 0:
        npts += 1
    info = bands_mod.analyze_gap(bs, spin=spin)
    if info.is_metal:
        raise ErrorDeUso("el sistema es metálico: la masa efectiva por "
                         "ajuste parabólico no aplica")

    atoms = struct_mod.primitive(atoms)
    common = sweep.prepare_common(atoms, pseudo_dir, ecutwfc, ecutrho,
                                  insulator)
    out = Path(outdir); out.mkdir(parents=True, exist_ok=True)
    res = bs.result
    recip = res.reciprocal                      # filas b_i en Å⁻¹

    lineas, meta = [], []
    for carrier, kidx in (("hueco", info.vbm_kindex),
                          ("electrón", info.cbm_kindex)):
        k0 = res.kpoints_cart[kidx]
        for nombre, e in valley_directions(k0):
            for t in np.linspace(-half_width, half_width, npts):
                kcart = k0 + t * e
                kfrac = kcart @ np.linalg.inv(recip)
                lineas.append(kfrac)
            meta.append((carrier, nombre, npts, kidx))

    card = ["K_POINTS crystal", str(len(lineas))]
    for kf in lineas:
        card.append(f"  {kf[0]:12.8f} {kf[1]:12.8f} {kf[2]:12.8f}  1.0")

    text = inputgen.build_pw_input(
        atoms=atoms, pseudos=common["pseudos"], calculation="bands",
        prefix=common["prefix"], pseudo_dir=common["pseudo_dir"],
        ecutwfc=common["ecutwfc"], ecutrho=common["ecutrho"],
        kcard="\n".join(card) + "\n", insulator=insulator,
        degauss=common["degauss"], smearing=common["smearing"],
        nbnd=res.nbnd,
    )
    sweep.write_input(out / "masa.in", text)

    grid = sweep.default_grid(atoms)
    scf = inputgen.build_pw_input(
        atoms=atoms, pseudos=common["pseudos"], calculation="scf",
        prefix=common["prefix"], pseudo_dir=common["pseudo_dir"],
        ecutwfc=common["ecutwfc"], ecutrho=common["ecutrho"],
        kcard=f"K_POINTS automatic\n  {grid[0]} {grid[1]} {grid[2]} 0 0 0\n",
        insulator=insulator, degauss=common["degauss"],
        smearing=common["smearing"],
    )
    sweep.write_input(out / "scf.in", scf)

    rep = ["--- Masa efectiva: camino fino ---",
           f"VBM en {info.vbm_label} ({info.vbm:.4f} eV)  |  "
           f"CBM en {info.cbm_label} ({info.cbm:.4f} eV)",
           f"{len(meta)} líneas ({npts} puntos cada una, "
           f"±{half_width} Å⁻¹):",
           "  " + "  |  ".join(
               f"{c}: {', '.join(d for cc, d, _, _ in meta if cc == c)}"
               for c in dict.fromkeys(x[0] for x in meta)),
           f"Total: {len(lineas)} puntos k",
           "",
           f"Archivos en '{out.resolve()}': scf.in, masa.in",
           "Orden: pw.x -in scf.in  ->  pw.x -in masa.in",
           "Después: olla-dft effmass estructura.cif --collect -o " + str(outdir)]
    if sweep.writing_inputs():
        (out / "masa_meta.json").write_text(json.dumps(
            {"lineas": [{"portador": c, "direccion": d, "npts": n,
                         "kindex": int(k)} for c, d, n, k in meta]},
            ensure_ascii=False, indent=2))
    warn = sweep.missing_pseudo_warning(common)
    if warn:
        rep.append(warn)
    return meta, "\n".join(rep)


def load_meta(outdir) -> list:
    """Recupera la descripción de las líneas escrita por prepare()."""
    f = Path(outdir) / "masa_meta.json"
    if not f.exists():
        raise FileNotFoundError(
            f"falta {f}: corre primero 'olla-dft effmass ... --bands-dir ...' "
            "para preparar el cálculo fino")
    d = json.loads(f.read_text())
    return [(x["portador"], x["direccion"], x["npts"], x["kindex"])
            for x in d["lineas"]]


def collect_fine(xml_path, meta, spin: int = 0,
                 degen_tol: float = DEGEN_TOL) -> EffMassRun:
    """Lee el cálculo fino y ajusta una masa por línea."""
    res = qeout.read_xml(xml_path)
    E = res.eigenvalues[spin]
    kc = res.kpoints_cart
    run = EffMassRun(source="camino fino dedicado")

    # Identificar valencia/conducción por el CONTEO DE ELECTRONES, no por el
    # nivel de Fermi: un cálculo 'bands' se corre sobre puntos k arbitrarios
    # y QE no calcula ahí un E_F utilizable. Con nelec la asignación es
    # exacta para un sistema de capa cerrada.
    vb = None
    if res.nspin == 1 and res.nelec:
        n_occ = int(round(res.nelec / 2))
        if 0 < n_occ <= E.shape[1]:
            vb = n_occ - 1
    if vb is None:                       # espín polarizado o nelec ausente
        ref = res.fermi if res.fermi is not None else res.homo
        if ref is None:
            ref = float(np.median(E))
        band_max = E.max(axis=0)
        candidatos = np.where(band_max <= ref + 1e-6)[0]
        vb = int(np.max(candidatos)) if len(candidatos) else 0
    cb = min(vb + 1, E.shape[1] - 1)
    if cb == vb:
        raise FaltanDatos("el cálculo no tiene bandas de conducción; "
                         "aumenta nbnd")

    pos = 0
    for carrier, nombre, npts, _ in meta:
        sl = slice(pos, pos + npts); pos += npts
        kk = kc[sl]
        c = len(kk) // 2
        t = np.linalg.norm(kk - kk[c], axis=1)
        t = np.where(np.arange(len(kk)) < c, -t, t)

        # Ajustar TODAS las bandas degeneradas con el extremo, no solo la
        # primera: en el silicio el hueco pesado y el ligero son bandas
        # distintas que coinciden en Γ y se separan al alejarse.
        if carrier == "hueco":
            candidatas = [b for b in range(vb, -1, -1)
                          if E[sl, vb][c] - E[sl, b][c] <= degen_tol]
        else:
            candidatas = [b for b in range(cb, E.shape[1])
                          if E[sl, b][c] - E[sl, cb][c] <= degen_tol]

        for band in candidatas:
            y = E[sl, band]
            coef = np.polyfit(t, y, 2)
            run.fits.append(MassFit(
                carrier=carrier, band=band, kindex=pos - npts + c,
                k_label="", direction=nombre,
                mass=_mass_from_quadratic(float(coef[0])),
                energy=float(y[c]), npts=npts,
                window=float(t.max() - t.min()), r2=_r2(t, y, coef),
            ))
    run.vbm = float(E[:, vb].max())
    run.cbm = float(E[:, cb].min())
    return run


# ----------------------------------------------------------------------
def report(run: EffMassRun) -> str:
    if run.is_metal:
        return ("--- Masa efectiva ---\nEl sistema es metálico: no hay un "
                "extremo de banda aislado que ajustar.\nPara un metal lo "
                "comparable es la masa de banda en la superficie de Fermi,\n"
                "que necesita otro tipo de cálculo.")
    lines = ["--- Masa efectiva (m*/mₑ) ---",
             f"Fuente: {run.source}"]
    if run.vbm is not None and run.cbm is not None:
        lines.append(f"VBM = {run.vbm:.4f} eV   CBM = {run.cbm:.4f} eV   "
                     f"gap = {run.cbm - run.vbm:.4f} eV")
    lines += ["",
              f"{'portador':10s} {'banda':>6s} {'m*/me':>9s} {'R²':>7s} "
              f"{'pts':>4s} {'Δk(Å⁻¹)':>9s}  dirección"]
    for f in run.fits:
        m = "   n/d" if f.mass is None or not np.isfinite(f.mass) \
            else f"{f.mass:9.3f}"
        r2 = "  n/d" if f.r2 is None or not np.isfinite(f.r2) \
            else f"{f.r2:7.4f}"
        etiqueta = f"{f.k_label} " if f.k_label else ""
        lines.append(f"{f.carrier:10s} {f.band + 1:6d} {m} {r2} "
                     f"{f.npts:4d} {f.window:9.4f}  {etiqueta}{f.direction}")
        if f.warning:
            lines.append(f"           ↳ {f.warning}")
    lines += ["",
              "El signo sale del ajuste: negativo = curvatura hacia abajo "
              "(hueco).",
              "Un R² de 1.0000 con 3 o 4 puntos no dice nada — una parábola "
              "pasa exacta\npor tres puntos cualesquiera."]
    if any(f.carrier == "hueco" for f in run.fits):
        lines += ["",
                  "OJO con los huecos: este cálculo NO incluye acoplamiento "
                  "espín-órbita, así\nque cerca de Γ hay un triplete "
                  "degenerado, no el par hueco pesado / hueco\nligero del "
                  "modelo de Luttinger. Los valores tabulados en la "
                  "literatura sí lo\nincluyen: coinciden bien en [100] y "
                  "[111] y pueden discrepar en [110]."]
    return "\n".join(lines)


def export(run: EffMassRun, outdir: str = ".") -> list:
    out = Path(outdir); out.mkdir(parents=True, exist_ok=True)
    f = out / "MASA_EFECTIVA.dat"
    lines = [provenance.header("masa efectiva", {"fuente": run.source},
                               titulo="Masa efectiva"),
             f"# {'portador':10s} {'banda':>6s} {'m*/me':>10s} {'R2':>8s} "
             f"{'pts':>4s} {'dk(A^-1)':>10s}  direccion"]
    for x in run.fits:
        lines.append(f"{x.carrier:12s} {x.band + 1:6d} {x.mass:10.4f} "
                     f"{x.r2:8.4f} {x.npts:4d} {x.window:10.4f}  "
                     f"{x.direction}")
    f.write_text("\n".join(lines) + "\n")
    return [str(f)]
