# Olla-DFT — command-line toolkit for Quantum ESPRESSO
# Copyright (C) 2026 Jorge Enrique González Sevilla
# SPDX-License-Identifier: GPL-3.0-or-later
# This program is free software: you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by the Free
# Software Foundation, either version 3 of the License, or (at your option)
# any later version. See the LICENSE file for details.

"""Constantes elásticas por el método de esfuerzo–deformación.

Se aplican deformaciones pequeñas a la celda, se calcula el tensor de
esfuerzos en cada una y se ajusta la relación lineal σ = C·ε. Con las seis
deformaciones independientes de la notación de Voigt queda determinada la
matriz completa 6×6.

Por qué esfuerzo–deformación y no energía–deformación: pw.x ya devuelve el
tensor de esfuerzos analítico, así que cada deformación aporta seis
ecuaciones en vez de una. Se necesitan muchos menos cálculos y no hay que
derivar numéricamente una energía.

Dos cuidados que determinan si el resultado sirve:

- **La deformación debe ser pequeña pero no diminuta.** Muy grande y se sale
  del régimen lineal; muy pequeña y el ruido numérico del esfuerzo domina.
  Se usan valores ±δ a ambos lados y se ajusta por mínimos cuadrados, lo que
  cancela la parte cuadrática.
- **La estructura de partida debe estar relajada.** Si la celda tiene
  esfuerzo residual, aparece como término independiente en el ajuste; se
  informa para que se vea.
"""

from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

from qekit.core import qeout
from qekit.core import provenance, structure
from qekit.core import style as qstyle
from qekit.modules import sweep
from qekit.core.errors import ErrorDeUso

# Deformaciones de Voigt: 1..3 normales, 4..6 de cizalla (γ = 2ε)
VOIGT = [(0, 0), (1, 1), (2, 2), (1, 2), (0, 2), (0, 1)]

# Familias cristalinas y sus constantes independientes
FAMILIES = {
    "cúbico": ["C11", "C12", "C44"],
    "hexagonal": ["C11", "C12", "C13", "C33", "C44"],
    "tetragonal": ["C11", "C12", "C13", "C33", "C44", "C66"],
    "ortorrómbico": ["C11", "C12", "C13", "C22", "C23", "C33",
                     "C44", "C55", "C66"],
    "trigonal": ["C11", "C12", "C13", "C14", "C33", "C44"],
    "monoclínico/triclínico": ["matriz completa"],
}


def crystal_family(spacegroup_number: int) -> str:
    n = spacegroup_number
    if n >= 195:
        return "cúbico"
    if n >= 168:
        return "hexagonal"
    if n >= 143:
        return "trigonal"
    if n >= 75:
        return "tetragonal"
    if n >= 16:
        return "ortorrómbico"
    return "monoclínico/triclínico"


# En una lámina solo tres componentes de la deformación tienen sentido:
# ε1 y ε2 en el plano y ε6 la cizalla del plano. Estirar ε3 deforma el VACÍO,
# no el material: da constantes que dependen de cuánto vacío se puso y no de
# la física. Restringirse a estas tres es además la mitad de cálculos.
COMPONENTES_2D = (0, 1, 5)

# 1 GPa · Å = 0.1 N/m
GPA_A_NM = 0.1


def constantes_2d(C: np.ndarray, altura: float) -> np.ndarray:
    """Pasa la matriz 3D (GPa) a constantes de lámina (N/m).

    Se multiplica por la ALTURA DE LA CELDA, no por un espesor supuesto. Es
    deliberado: pw.x calculó el esfuerzo dividiendo por el volumen A·c, así
    que C_3D es inversamente proporcional a c y el producto C_3D·c no depende
    del vacío que se haya puesto. Un "espesor de capa" (3.35 Å para el
    grafeno y demás) es un convenio, y usarlo aquí metería una arbitrariedad
    en el número principal.
    """
    idx = np.ix_(COMPONENTES_2D, COMPONENTES_2D)
    return np.array(C)[idx] * altura * GPA_A_NM


def born_2d(C2: np.ndarray) -> tuple:
    """Criterio de estabilidad mecánica de Born en 2D.

    Para un cristal bidimensional con C11, C22, C12 y C66:
        C11 > 0,  C66 > 0,  C11·C22 - C12² > 0
    """
    c11, c22, c12, c66 = C2[0, 0], C2[1, 1], C2[0, 1], C2[2, 2]
    det = c11 * c22 - c12 ** 2
    pruebas = [("C11 > 0", bool(c11 > 0)),
               ("C66 > 0", bool(c66 > 0)),
               ("C11·C22 − C12² > 0", bool(det > 0))]
    fallan = [n for n, ok in pruebas if not ok]
    return (not fallan), fallan


def modulos_2d(C2: np.ndarray) -> dict:
    """Young, Poisson, módulo de área y cizalla de una lámina, en N/m."""
    c11, c22, c12, c66 = (float(C2[0, 0]), float(C2[1, 1]),
                          float(C2[0, 1]), float(C2[2, 2]))
    out = {"C11": c11, "C22": c22, "C12": c12, "C66": c66}
    if c22 != 0:
        out["Y_x"] = (c11 * c22 - c12 ** 2) / c22
        out["nu_x"] = c12 / c22
    if c11 != 0:
        out["Y_y"] = (c11 * c22 - c12 ** 2) / c11
        out["nu_y"] = c12 / c11
    out["K"] = (c11 + c22 + 2 * c12) / 4.0     # módulo de área 2D
    out["G"] = c66
    return out


def strain_matrix(component: int, delta: float) -> np.ndarray:
    """Matriz de deformación para la componente de Voigt dada (0..5).

    Para las de cizalla se reparte γ/2 en las dos posiciones simétricas,
    que es la convención de la notación de Voigt: ε₄ = 2ε₂₃.
    """
    e = np.zeros((3, 3))
    i, j = VOIGT[component]
    if i == j:
        e[i, j] = delta
    else:
        e[i, j] = e[j, i] = delta / 2.0
    return e


@dataclass
class ElasticRun:
    deltas: list = field(default_factory=list)
    components: list = field(default_factory=list)
    jobs: list = field(default_factory=list)
    stresses: list = field(default_factory=list)   # (3,3) GPa o None
    natoms: int = 1
    volume: float = None
    family: str = ""
    spacegroup: str = ""
    reference_stress: np.ndarray = None            # esfuerzo de la celda sin deformar
    C: np.ndarray = None                           # matriz 6x6 en GPa
    dosd: bool = False                             # modo lámina (2D)
    altura: float = None                           # altura de la celda en c (Å)
    espesor: float = None                          # espesor supuesto para pasar a GPa

    @property
    def delta(self) -> float:
        """Deformación máxima aplicada, en fracción (0.010 = 1 %).

        No se guarda como campo aparte: se deduce de las deformaciones que de
        hecho se usaron, para que al releer un cálculo viejo con --collect el
        encabezado diga la deformación REAL del barrido y no la que traiga por
        omisión la versión actual del código.
        """
        vals = [abs(d) for d in self.deltas if d]
        return max(vals) if vals else 0.0


# ----------------------------------------------------------------------
# Preparación
# ----------------------------------------------------------------------
def prepare(atoms, outdir: str = "elastic", delta: float = 0.010,
            npoints: int = 4, pseudo_dir: str = None, insulator: bool = False,
            ecutwfc: float = None, ecutrho: float = None,
            kspacing: float = None, ion_mode: str = "auto",
            dosd: bool = False, espesor: float = None) -> tuple:
    """Genera la celda sin deformar más las deformadas.

    `npoints` es el número de deformaciones NO nulas por componente (pares,
    repartidas simétricamente: 4 → ±δ/2, ±δ).

    `ion_mode` decide qué pasa con las posiciones internas en cada celda
    deformada, y su valor por defecto viene de un problema real:

    - 'auto' (recomendado): deformaciones NORMALES (ε1–ε3) con iones fijos y
      CIZALLAS (ε4–ε6) con relajación. En muchas estructuras de alta
      simetría, bajo deformación normal las posiciones quedan fijadas por
      simetría y las fuerzas son puro ruido numérico; si se deja al BFGS
      perseguirlas, desplaza los átomos y contamina el esfuerzo (verificado
      en silicio: σyz espurio de ~0.3 GPa y C11 inflada de 159 a 190 GPa).
      Bajo cizalla, en cambio, el desplazamiento interno es física real
      (parámetro de Kleinman) y omitirlo sobreestima C44 seriamente.
    - 'relax': relajar en todas las componentes. Necesario si tu cristal
      tiene parámetros internos libres también bajo deformación normal
      (por ejemplo el parámetro u de la wurtzita bajo ε3).
    - 'fixed': iones fijos en todas: constantes "clamped-ion", útiles solo
      como referencia.
    """
    # ------------------------------------------------------------------
    # ORIENTACIÓN ESTÁNDAR — imprescindible, no cosmética.
    #
    # Las Cij se definen en el marco cristalofísico: para un cúbico, x∥[100];
    # para un hexagonal, c∥z. Pero un CIF puede traer la celda primitiva en
    # cualquier orientación (el convenio a₁∥x es común), y ahí "ε₁" deforma
    # a lo largo de una diagonal del cubo, no de [100]: la matriz medida es
    # la del tensor ROTADO, la simetrización por familia deja de ser válida
    # y aparecen σ fuera de la diagonal que parecen error sin serlo. Se
    # detectó comparando dos celdas del mismo silicio: la alineada daba
    # C11 = 158 GPa y la rotada 221, con B idéntico (94 GPa) porque B es
    # invariante ante rotaciones. Por eso aquí se lleva SIEMPRE la
    # estructura a la celda primitiva estandarizada de spglib, cuya
    # orientación es la del marco cristalofísico.
    # ------------------------------------------------------------------
    from qekit.core import kpoints as kp

    if dosd:
        # En una lámina NO se estandariza a la primitiva: spglib puede
        # reorientar la celda y dejar el vacío en otro eje, y entonces "el
        # plano ab" ya no es el plano del material. Se trabaja con la celda
        # tal como viene, que es la que tiene el vacío donde el usuario lo puso.
        if 2 not in kp.direcciones_con_vacio(atoms):
            raise ErrorDeUso(
                "--2d espera una lámina con vacío en c, y esta celda no lo "
                "tiene. Si es material en bulto, quita --2d; si es una "
                "monocapa, añádele vacío (por ejemplo 'olla-dft layers --slab "
                "mono.cif --vacuum 20').")
        atoms_in = atoms
        reoriented = False
    else:
        atoms_in = atoms
        atoms = structure.primitive(atoms)
        reoriented = (len(atoms) != len(atoms_in)
                      or not np.allclose(atoms.cell.array, atoms_in.cell.array,
                                         atol=1e-6))

    common = sweep.prepare_common(atoms, pseudo_dir, ecutwfc, ecutrho, insulator)
    out = Path(outdir); out.mkdir(parents=True, exist_ok=True)
    grid = sweep.default_grid(atoms, kspacing)

    ds = structure.symmetry_dataset(atoms)
    run = ElasticRun(
        natoms=len(atoms),
        volume=float(abs(np.linalg.det(atoms.cell.array))),
        family=crystal_family(ds.number),
        spacegroup=f"{ds.international} (N.º {ds.number})",
        dosd=dosd,
        altura=float(np.linalg.norm(atoms.cell.array[2])) if dosd else None,
        espesor=espesor,
    )

    if npoints % 2 or npoints < 2:
        raise ErrorDeUso("npoints debe ser un número par >= 2 "
                         "(las deformaciones van en pares ±)")
    if ion_mode not in ("auto", "relax", "fixed"):
        raise ErrorDeUso("ion_mode debe ser auto, relax o fixed")
    half = npoints // 2
    magnitudes = [delta * (k + 1) / half for k in range(half)]
    amounts = sorted([-m for m in magnitudes] + magnitudes)

    modo_txt = {
        "auto": "fijas en ε1–ε3, relajadas en ε4–ε6 (auto)",
        "relax": "relajadas en todas",
        "fixed": "fijas en todas (clamped-ion)",
    }[ion_mode]
    report = ["--- Constantes elásticas (esfuerzo–deformación) ---",
              f"Estructura: {atoms.get_chemical_formula()} ({len(atoms)} átomos)",
              f"Grupo espacial: {run.spacegroup}  ->  familia {run.family}",
              "Constantes independientes: "
              + ("C11, C22, C12, C66 (lámina)" if dosd
                 else ", ".join(FAMILIES.get(run.family, ['?']))),
              f"Deformaciones: {amounts} en cada una de las "
              + ("3 componentes del plano (ε1, ε2, ε6)" if dosd
                 else "6 componentes de Voigt"),
              f"Malla k: {grid[0]}x{grid[1]}x{grid[2]}  |  "
              f"posiciones internas: {modo_txt}"]
    if dosd:
        report.append(
            f"Modo lámina: las Cij se darán en N/m multiplicando por la altura "
            f"de la celda\n  (c = {run.altura:.3f} Å). No se toca ε3: estirar "
            "el vacío no mide nada.")
    if reoriented:
        report.append(
            "La estructura se llevó a la celda primitiva estandarizada para que\n"
            "los ejes cartesianos coincidan con los cristalofísicos (las Cij se\n"
            "definen en ese marco)."
        )
    warn = sweep.missing_pseudo_warning(common)
    if warn:
        report.append(warn)

    # punto de referencia sin deformar, para medir el esfuerzo residual
    job0 = sweep.write_scf_job(atoms, common, out / "ref", "sin deformar", grid,
                               meta={"component": None, "delta": 0.0},
                               calculation="scf")
    run.jobs.append(job0)
    run.components.append(None); run.deltas.append(0.0)

    componentes = COMPONENTES_2D if dosd else range(6)
    for comp in componentes:
        if ion_mode == "auto":
            calc = "scf" if comp < 3 else "relax"
        else:
            calc = "relax" if ion_mode == "relax" else "scf"
        for amt in amounts:
            e = strain_matrix(comp, amt)
            deformed = atoms.copy()
            new_cell = atoms.cell.array @ (np.eye(3) + e)
            deformed.set_cell(new_cell, scale_atoms=True)
            label = f"ε{comp + 1} = {amt:+.4f}"
            job = sweep.write_scf_job(
                deformed, common, out / f"e{comp + 1}_{amt:+.4f}".replace("+", "p").replace("-", "m"),
                label, grid, meta={"component": comp, "delta": amt},
                calculation=calc,
                # esfuerzo limpio: scf bien convergido en densidad
                conv_thr=1.0e-9,
            )
            run.jobs.append(job)
            run.components.append(comp); run.deltas.append(amt)

    sweep.write_run_script(run.jobs, out / "run.sh")
    report += ["", f"{len(run.jobs)} cálculos escritos en '{out.resolve()}'",
               "Córrelos con --run, o a mano con ./run.sh dentro de esa carpeta."]
    return run, "\n".join(report)


# ----------------------------------------------------------------------
# Recolección y ajuste
# ----------------------------------------------------------------------
def collect(run: ElasticRun, results: list = None) -> ElasticRun:
    run.stresses = []
    if results is not None:
        by_dir = {str(r.job.directory): r for r in results}
        for job in run.jobs:
            r = by_dir.get(str(job.directory))
            run.stresses.append(r.result.stress if (r and r.ok and r.result) else None)
    else:
        for job in run.jobs:
            try:
                run.stresses.append(qeout.read_xml(str(job.directory)).stress)
            except Exception:
                run.stresses.append(None)
    if run.stresses and run.stresses[0] is not None:
        run.reference_stress = run.stresses[0]
    return run


def _voigt(stress: np.ndarray) -> np.ndarray:
    """Tensor 3x3 -> vector de Voigt de 6 componentes."""
    return np.array([stress[0, 0], stress[1, 1], stress[2, 2],
                     stress[1, 2], stress[0, 2], stress[0, 1]])


def fit(run: ElasticRun) -> np.ndarray:
    """Ajusta σ = C·ε columna por columna. Devuelve C en GPa.

    **Signo.** pw.x no escribe el esfuerzo con el convenio de la elasticidad.
    Comprobado sobre silicio estirado un 1 % en x: la celda estirada está en
    tracción, que en la definición de Cij es σ > 0, y pw.x reporta
    σxx = −20.8 kbar. Además pw.x calcula su presión como P = +tr(σ)/3,
    mientras que el convenio estándar es P = −tr(σ)/3. Es decir, el tensor
    que escribe pw.x es el opuesto del de la elasticidad, y hay que
    invertirlo: C = −∂σ_pw/∂ε. Sin esta inversión toda la matriz sale
    negativa y el criterio de estabilidad de Born falla siempre.
    """
    C = np.full((6, 6), np.nan)
    ref = run.reference_stress
    for comp in (COMPONENTES_2D if run.dosd else range(6)):
        xs, ys = [], []
        for c, d, s in zip(run.components, run.deltas, run.stresses):
            if c != comp or s is None:
                continue
            sigma = _voigt(s)
            if ref is not None:
                sigma = sigma - _voigt(ref)     # descontar el esfuerzo residual
            xs.append(d); ys.append(sigma)
        if len(xs) < 2:
            continue
        xs = np.array(xs); ys = np.array(ys)      # (n, 6)
        # pendiente por mínimos cuadrados de cada componente de esfuerzo,
        # con el signo invertido por lo explicado en el docstring
        for row in range(6):
            slope = np.polyfit(xs, ys[:, row], 1)[0]
            C[row, comp] = -slope
    run.C = C
    return C


def symmetrize(C: np.ndarray, family: str) -> np.ndarray:
    """Impone la simetría de la familia cristalina promediando equivalentes."""
    C = np.array(C, dtype=float)
    S = 0.5 * (C + C.T)                 # la matriz elástica es simétrica

    def avg(pairs):
        vals = [S[i, j] for i, j in pairs if np.isfinite(S[i, j])]
        return np.mean(vals) if vals else np.nan

    out = np.array(S)
    if family == "cúbico":
        c11 = avg([(0, 0), (1, 1), (2, 2)])
        c12 = avg([(0, 1), (0, 2), (1, 2)])
        c44 = avg([(3, 3), (4, 4), (5, 5)])
        out = np.zeros((6, 6))
        out[0, 0] = out[1, 1] = out[2, 2] = c11
        out[0, 1] = out[1, 0] = out[0, 2] = out[2, 0] = out[1, 2] = out[2, 1] = c12
        out[3, 3] = out[4, 4] = out[5, 5] = c44
    elif family == "hexagonal":
        c11 = avg([(0, 0), (1, 1)])
        c12 = avg([(0, 1)])
        c13 = avg([(0, 2), (1, 2)])
        c33 = avg([(2, 2)])
        c44 = avg([(3, 3), (4, 4)])
        out = np.zeros((6, 6))
        out[0, 0] = out[1, 1] = c11
        out[0, 1] = out[1, 0] = c12
        out[0, 2] = out[2, 0] = out[1, 2] = out[2, 1] = c13
        out[2, 2] = c33
        out[3, 3] = out[4, 4] = c44
        out[5, 5] = (c11 - c12) / 2.0
    elif family == "tetragonal":
        c11 = avg([(0, 0), (1, 1)])
        out[0, 0] = out[1, 1] = c11
        out[0, 2] = out[2, 0] = out[1, 2] = out[2, 1] = avg([(0, 2), (1, 2)])
        out[3, 3] = out[4, 4] = avg([(3, 3), (4, 4)])
    return out


# ----------------------------------------------------------------------
# Propiedades mecánicas
# ----------------------------------------------------------------------
@dataclass
class Moduli:
    B_voigt: float = None
    B_reuss: float = None
    B_hill: float = None
    G_voigt: float = None
    G_reuss: float = None
    G_hill: float = None
    E: float = None            # módulo de Young (GPa)
    nu: float = None           # razón de Poisson
    pugh: float = None         # B/G
    anisotropy: float = None   # índice universal A^U
    stable: bool = None
    stability_note: str = ""


def moduli(C: np.ndarray) -> Moduli:
    """Promedios de Voigt–Reuss–Hill y propiedades derivadas."""
    m = Moduli()
    if not np.all(np.isfinite(C)):
        return m
    try:
        S = np.linalg.inv(C)                    # matriz de compliancia
    except np.linalg.LinAlgError:
        return m

    c = C
    m.B_voigt = ((c[0, 0] + c[1, 1] + c[2, 2])
                 + 2.0 * (c[0, 1] + c[1, 2] + c[0, 2])) / 9.0
    m.G_voigt = ((c[0, 0] + c[1, 1] + c[2, 2])
                 - (c[0, 1] + c[1, 2] + c[0, 2])
                 + 3.0 * (c[3, 3] + c[4, 4] + c[5, 5])) / 15.0
    m.B_reuss = 1.0 / ((S[0, 0] + S[1, 1] + S[2, 2])
                       + 2.0 * (S[0, 1] + S[1, 2] + S[0, 2]))
    m.G_reuss = 15.0 / (4.0 * (S[0, 0] + S[1, 1] + S[2, 2])
                        - 4.0 * (S[0, 1] + S[1, 2] + S[0, 2])
                        + 3.0 * (S[3, 3] + S[4, 4] + S[5, 5]))
    m.B_hill = 0.5 * (m.B_voigt + m.B_reuss)
    m.G_hill = 0.5 * (m.G_voigt + m.G_reuss)
    if m.B_hill and m.G_hill:
        m.E = 9.0 * m.B_hill * m.G_hill / (3.0 * m.B_hill + m.G_hill)
        m.nu = (3.0 * m.B_hill - 2.0 * m.G_hill) / (
            2.0 * (3.0 * m.B_hill + m.G_hill))
        m.pugh = m.B_hill / m.G_hill if m.G_hill else None
        if m.B_reuss and m.G_reuss:
            m.anisotropy = (5.0 * m.G_voigt / m.G_reuss
                            + m.B_voigt / m.B_reuss - 6.0)

    # Estabilidad mecánica: criterio general de Born — la matriz elástica
    # debe ser definida positiva.
    eig = np.linalg.eigvalsh(0.5 * (C + C.T))
    m.stable = bool(np.all(eig > 0))
    if m.stable:
        m.stability_note = ("todos los valores propios de C son positivos: "
                            "la estructura es mecánicamente estable.")
    else:
        neg = ", ".join(f"{e:.1f}" for e in eig if e <= 0)
        m.stability_note = (f"hay valores propios no positivos ({neg} GPa): la "
                            "estructura NO es\nmecánicamente estable, o el "
                            "cálculo no está bien convergido.")
    return m


# ----------------------------------------------------------------------
# Reporte y exportación
# ----------------------------------------------------------------------
def report(run: ElasticRun, symmetrized: bool = True) -> str:
    lines = ["--- Constantes elásticas ---",
             f"Estructura: {run.natoms} átomos, V = {run.volume:.4f} Å³",
             f"Grupo espacial: {run.spacegroup}  ->  familia {run.family}"]

    done = sum(1 for s in run.stresses if s is not None)
    lines.append(f"Cálculos con esfuerzo leído: {done} de {len(run.jobs)}")
    if run.reference_stress is not None:
        p = np.trace(run.reference_stress) / 3.0
        lines.append(f"Esfuerzo residual de la celda sin deformar: "
                     f"{p:.3f} GPa")
        if abs(p) > 0.5:
            lines.append("  AVISO: es alto. Relaja la celda con vc-relax antes "
                         "de calcular\n  las constantes elásticas, o los valores "
                         "saldrán sesgados.")
    if done < 3:
        lines.append("\nNo hay suficientes cálculos terminados para ajustar.")
        return "\n".join(lines)

    C = run.C if run.C is not None else fit(run)

    if run.dosd:
        return _report_2d(run, C, lines)

    Cs = symmetrize(C, run.family) if symmetrized else 0.5 * (C + C.T)

    lines += ["", "Matriz elástica C (GPa):"]
    header = "      " + "".join(f"{i + 1:>10d}" for i in range(6))
    lines.append(header)
    for i in range(6):
        row = "".join(f"{Cs[i, j]:10.2f}" if np.isfinite(Cs[i, j]) else f"{'—':>10s}"
                      for j in range(6))
        lines.append(f"  {i + 1:>2d}  {row}")

    indep = FAMILIES.get(run.family, [])
    if indep and indep != ["matriz completa"]:
        lines += ["", "Constantes independientes:"]
        for name in indep:
            i, j = int(name[1]) - 1, int(name[2]) - 1
            lines.append(f"  {name} = {Cs[i, j]:8.2f} GPa")

    m = moduli(Cs)
    if m.B_hill is not None:
        lines += ["", "Módulos elásticos (promedios de Voigt–Reuss–Hill):",
                  f"  {'':10s} {'Voigt':>10s} {'Reuss':>10s} {'Hill':>10s}",
                  f"  {'B (GPa)':10s} {m.B_voigt:10.2f} {m.B_reuss:10.2f} "
                  f"{m.B_hill:10.2f}",
                  f"  {'G (GPa)':10s} {m.G_voigt:10.2f} {m.G_reuss:10.2f} "
                  f"{m.G_hill:10.2f}",
                  "",
                  f"  Módulo de Young E = {m.E:.2f} GPa",
                  f"  Razón de Poisson ν = {m.nu:.4f}",
                  f"  Cociente de Pugh B/G = {m.pugh:.3f}  "
                  f"({'dúctil' if m.pugh > 1.75 else 'frágil'}, umbral 1.75)"]
        if m.anisotropy is not None:
            lines.append(f"  Anisotropía universal A^U = {m.anisotropy:.4f}  "
                         f"({'isótropo' if abs(m.anisotropy) < 0.01 else 'anisótropo'})")
        lines += ["", "Estabilidad mecánica (criterio de Born):",
                  f"  {m.stability_note}"]
    return "\n".join(lines)


def _report_2d(run: ElasticRun, C: np.ndarray, lines: list) -> str:
    """Reporte del modo lámina: todo en N/m, criterios de Born 2D."""
    C2 = constantes_2d(0.5 * (C + C.T), run.altura)
    m = modulos_2d(C2)
    etiquetas = ("1 (xx)", "2 (yy)", "6 (xy)")

    lines += ["", f"Altura de la celda: c = {run.altura:.4f} Å  "
                  f"(el vacío se cancela al multiplicar)",
              "", "Constantes elásticas de lámina (N/m):",
              "        " + "".join(f"{e:>12s}" for e in etiquetas)]
    for i, e in enumerate(etiquetas):
        fila = "".join(f"{C2[i, j]:12.2f}" if np.isfinite(C2[i, j])
                       else f"{'—':>12s}" for j in range(3))
        lines.append(f"  {e:>6s}{fila}")

    lines += ["", "Constantes independientes (N/m):",
              f"  C11 = {m['C11']:9.2f}",
              f"  C22 = {m['C22']:9.2f}",
              f"  C12 = {m['C12']:9.2f}",
              f"  C66 = {m['C66']:9.2f}"]
    iso = abs(m["C11"] - m["C22"]) < 0.02 * max(abs(m["C11"]), 1.0)
    if iso:
        c66_esp = (m["C11"] - m["C12"]) / 2.0
        lines.append(f"  C11 ≈ C22: lámina isótropa en el plano; para serlo del "
                     f"todo\n  debe cumplirse C66 = (C11−C12)/2 = "
                     f"{c66_esp:.2f} N/m  (calculado: {m['C66']:.2f})")
        if abs(c66_esp - m["C66"]) > 0.03 * max(abs(c66_esp), 1.0):
            lines.append(
                "  La diferencia es apreciable. Antes de leerla como "
                "anisotropía, mira que\n  con --ion-mode auto (el que viene "
                "por omisión) las cizallas se relajan y\n  las normales no: "
                "el relajamiento interno baja C66 y no toca C11 ni C12, así\n"
                "  que la identidad deja de cumplirse aunque la lámina sea "
                "isótropa. Para\n  compararlas en igualdad de condiciones:  "
                "--ion-mode fixed  (o  relax).")

    lines += ["", "Módulos de lámina (N/m):"]
    if "Y_x" in m:
        lines.append(f"  Módulo de Young 2D   Yx = {m['Y_x']:8.2f}   "
                     f"Yy = {m.get('Y_y', float('nan')):8.2f}")
        lines.append(f"  Razón de Poisson     νx = {m['nu_x']:8.4f}   "
                     f"νy = {m.get('nu_y', float('nan')):8.4f}")
    lines.append(f"  Módulo de área 2D     K = {m['K']:8.2f}")
    lines.append(f"  Módulo de cizalla     G = {m['G']:8.2f}")

    estable, fallan = born_2d(C2)
    lines += ["", "Estabilidad mecánica (criterio de Born en 2D):"]
    if estable:
        lines.append("  Estable: se cumplen C11 > 0, C66 > 0 y C11·C22 − C12² > 0.")
    else:
        lines.append("  INESTABLE: no se cumple " + "; ".join(fallan) + ".")
        lines.append("  Antes de concluir que la lámina no existe, comprueba "
                     "que la estructura\n  estaba relajada y que la malla k y "
                     "el cutoff están convergidos: una\n  celda con esfuerzo "
                     "residual da Cij sesgadas y puede fingir inestabilidad.")

    if run.espesor:
        f = 1.0 / (run.espesor * GPA_A_NM)
        lines += ["", f"Equivalente en 3D suponiendo un espesor de "
                      f"{run.espesor:g} Å:",
                  f"  C11 = {m['C11'] * f:8.2f} GPa    "
                  f"C12 = {m['C12'] * f:8.2f} GPa    "
                  f"C66 = {m['C66'] * f:8.2f} GPa"]
        if "Y_x" in m:
            lines.append(f"  Young Yx = {m['Y_x'] * f:8.2f} GPa")
        lines.append("  Este espesor es un CONVENIO, no una medida: los "
                     "números en GPa\n  cambian si se elige otro. Los N/m de "
                     "arriba no dependen de él.")
    return "\n".join(lines)


def export(run: ElasticRun, outdir: str = ".") -> list:
    out = Path(outdir); out.mkdir(parents=True, exist_ok=True)
    C = run.C if run.C is not None else fit(run)
    written = []
    fname = out / "ELASTIC_C.dat"
    if run.dosd:
        Cs = constantes_2d(0.5 * (C + C.T), run.altura)
        np.savetxt(fname, Cs, fmt="%14.6f",
                   header=provenance.header_plain(
                       "constantes de lamina", {"c_A": f"{run.altura:.4f}",
                                                "delta": run.delta},
                       titulo="Constantes elasticas 2D (N/m); orden 1(xx) 2(yy) 6(xy)"),
                   comments="# ")
    else:
        Cs = symmetrize(C, run.family)
        np.savetxt(fname, Cs, fmt="%14.6f",
                   header=provenance.header_plain(
                       "matriz elástica", {"familia": run.family,
                                           "delta": run.delta},
                       titulo="Matriz elastica C (GPa); Voigt 1..6"),
                   comments="# ")
    written.append(str(fname))
    txt = out / "ELASTIC.txt"
    txt.write_text(report(run) + "\n")
    written.append(str(txt))
    return written


def plot(run: ElasticRun, outfile: str = "elastic", formats="pdf,png",
         theme: str = None, size: str = None, family: str = None,
         background: str = None, palette=None, usetex: bool = None,
         width="single", journal: str = "generic", aspect: float = 0.85,
         mono: bool = False, dpi: int = None) -> list:
    """Rectas esfuerzo–deformación: se ve de un vistazo si el ajuste es lineal."""
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError as exc:  # pragma: no cover
        raise RuntimeError("matplotlib no está instalado.") from exc

    st = qstyle.apply(theme, size=size, family=family, background=background,
                      palette=palette, usetex=usetex, mono=mono)
    fig, ax = qstyle.new_figure(width, journal, aspect)
    qstyle.palette(6, mono=mono)
    ref = run.reference_stress

    for comp in (COMPONENTES_2D if run.dosd else range(6)):
        xs, ys = [], []
        for c, d, s in zip(run.components, run.deltas, run.stresses):
            if c != comp or s is None:
                continue
            sigma = _voigt(s)
            if ref is not None:
                sigma = sigma - _voigt(ref)
            xs.append(d * 100.0); ys.append(sigma[comp])
        if len(xs) < 2:
            continue
        order = np.argsort(xs)
        xs = np.array(xs)[order]; ys = np.array(ys)[order]
        kw = qstyle.style_line(comp, 6, mono=mono)
        ax.plot(xs, ys, marker="o", ms=3.2, lw=st["line"],
                label=rf"$\varepsilon_{comp + 1}$", **kw)

    ax.axhline(0.0, color=qstyle.INK_FAINT, lw=st["axis_line"])
    ax.axvline(0.0, color=qstyle.INK_FAINT, lw=st["axis_line"])
    ax.set_xlabel(r"deformación $\varepsilon$ (\%)" if qstyle.USETEX
                  else "deformación ε (%)")
    ax.set_ylabel(r"$\sigma$ (GPa)")
    ax.legend(ncol=2, loc="best")
    written = qstyle.save(fig, outfile, formats, dpi=dpi,
                          modulo="constantes elásticas")
    plt.close(fig)
    return written
