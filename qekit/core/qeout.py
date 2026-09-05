# Olla-DFT — command-line toolkit for Quantum ESPRESSO
# Copyright (C) 2026 Jorge Enrique González Sevilla
# SPDX-License-Identifier: AGPL-3.0-or-later
# This program is free software: you can redistribute it and/or modify it
# under the terms of the GNU Affero General Public License as published by the Free
# Software Foundation, either version 3 of the License, or (at your option)
# any later version. See the LICENSE file for details.

"""Lectura de las salidas de Quantum ESPRESSO.

La fuente principal es el archivo XML que escribe pw.x (`outdir/prefix.xml`,
o `outdir/prefix.save/data-file-schema.xml` en versiones antiguas), porque
contiene todo lo necesario en un solo lugar: celda, k-points, eigenvalores,
ocupaciones, energía de Fermi y niveles HOMO/LUMO.

Convenciones de unidades dentro del XML de QE:
  - energías en Hartree
  - longitudes en bohr
  - k_point en coordenadas cartesianas, en unidades de 2*pi/alat
Aquí todo se convierte a eV y Å⁻¹.
"""

import xml.etree.ElementTree as ET
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
from qekit.core.errors import FaltanDatos

HARTREE_EV = 27.211386245988
BOHR_ANG = 0.529177210903
RY_EV = HARTREE_EV / 2.0
# 1 Ha/bohr^3 en GPa: el XML de pw.x escribe el esfuerzo en esas unidades
HA_BOHR3_GPA = 29421.02648438959


def _tag(elem) -> str:
    """Nombre de la etiqueta sin el namespace."""
    return elem.tag.split("}")[-1]


def _children(elem, name: str) -> list:
    return [c for c in elem if _tag(c) == name]


def _child(elem, name: str):
    found = _children(elem, name)
    return found[0] if found else None


def _text(elem, name: str, default=None):
    c = _child(elem, name)
    return c.text.strip() if (c is not None and c.text) else default


def _floats(elem) -> np.ndarray:
    return np.array([float(x) for x in elem.text.split()])


@dataclass
class QEResult:
    """Resultado de un cálculo de pw.x leído desde el XML."""

    prefix: str = ""
    xml_path: str = ""
    calculation: str = ""
    # estructura
    alat: float = 0.0                       # bohr
    cell: np.ndarray = None                 # (3,3) en Å
    reciprocal: np.ndarray = None           # (3,3) en Å⁻¹, con factor 2*pi
    symbols: list = field(default_factory=list)
    positions: np.ndarray = None            # (nat,3) en Å
    # bandas
    nbnd: int = 0
    nelec: float = 0.0
    nspin: int = 1                          # 1 o 2 (lsda)
    noncolin: bool = False
    kpoints_frac: np.ndarray = None         # (nk,3) fraccionarias
    kpoints_cart: np.ndarray = None         # (nk,3) en Å⁻¹
    eigenvalues: np.ndarray = None          # (nspin, nk, nbnd) en eV
    occupations: np.ndarray = None          # (nspin, nk, nbnd)
    weights: np.ndarray = None              # (nk,)
    fermi: float = None                     # eV
    homo: float = None                      # eV
    lumo: float = None                      # eV
    occupations_kind: str = ""
    total_energy: float = None              # eV
    stress: np.ndarray = None               # (3,3) en GPa
    pressure: float = None                  # GPa (traza/3, signo de QE)
    volume: float = None                    # Å³
    forces: np.ndarray = None               # (nat,3) en eV/Å
    # --- diagnóstico y procedencia (todo esto ya estaba en el XML) -----
    converged: bool = None                  # convergence_achieved
    n_scf_steps: int = None
    scf_error: float = None                 # Ry (comparable con conv_thr)
    n_bfgs_steps: int = None                # pasos iónicos de un relax
    functional: str = ""
    pseudo_files: dict = field(default_factory=dict)   # símbolo -> archivo
    ecutwfc: float = None                   # Ry
    ecutrho: float = None                   # Ry
    kgrid: tuple = None                     # (n1,n2,n3) si fue automática
    kshift: tuple = None
    smearing: str = ""
    degauss: float = None                   # Ry
    energy_terms: dict = field(default_factory=dict)   # eV, por término
    total_magnetization: float = None
    absolute_magnetization: float = None
    n_sym: int = None                       # operaciones que USÓ QE
    equivalent_atoms: np.ndarray = None
    cpu_time: float = None                  # s
    wall_time: float = None                 # s

    @property
    def max_force(self) -> float:
        """Fuerza residual máxima en eV/Å (None si no hay fuerzas)."""
        if self.forces is None or len(self.forces) == 0:
            return None
        return float(np.max(np.linalg.norm(self.forces, axis=1)))

    @property
    def fingerprint(self) -> tuple:
        """Los parámetros que TIENEN que coincidir para restar energías.

        Comparar energías totales de cálculos que difieren en cualquiera de
        estos es inválido, y QE no avisa: devuelve dos números perfectamente
        formados cuya diferencia no significa nada.
        """
        return (self.functional,
                tuple(sorted(self.pseudo_files.items())),
                self.ecutwfc, self.ecutrho,
                self.smearing, self.degauss,
                self.occupations_kind, self.nspin)

    @property
    def nk(self) -> int:
        return 0 if self.kpoints_frac is None else len(self.kpoints_frac)

    @property
    def is_spin_polarized(self) -> bool:
        return self.nspin == 2


def find_xml(path: str = ".", prefix: str = None) -> Path:
    """Localiza el XML de salida de pw.x.

    `path` puede ser el XML mismo, la carpeta del cálculo o el outdir.
    Busca en este orden: el archivo dado, ./out/*.xml, ./*.xml y
    */*.save/data-file-schema.xml.
    """
    p = Path(path)
    if p.is_file() and p.suffix == ".xml":
        return p
    if not p.exists():
        raise FileNotFoundError(f"no existe la ruta '{path}'")

    candidates = []
    search_dirs = [p, p / "out"]
    for d in search_dirs:
        if not d.is_dir():
            continue
        if prefix:
            direct = d / f"{prefix}.xml"
            if direct.is_file():
                return direct
            schema = d / f"{prefix}.save" / "data-file-schema.xml"
            if schema.is_file():
                return schema
        candidates += sorted(d.glob("*.xml"))
        candidates += sorted(d.glob("*.save/data-file-schema.xml"))

    # descartar XML que no sean de QE
    for cand in candidates:
        try:
            with open(cand, "r", errors="ignore") as fh:
                head = fh.read(2000)
            if "quantum-espresso" in head or "<qes:espresso" in head or "espresso" in head:
                return cand
        except OSError:
            continue
    raise FileNotFoundError(
        f"no se encontró el XML de salida de pw.x en '{path}'.\n"
        "Ejecuta primero el cálculo, o indica la ruta al archivo .xml "
        "(suele estar en la carpeta 'out/')."
    )


def read_xml(path: str = ".", prefix: str = None) -> QEResult:
    """Lee el XML de pw.x y devuelve un QEResult con todo en eV y Å."""
    xml_path = find_xml(path, prefix)
    root = ET.parse(xml_path).getroot()

    out = _child(root, "output")
    if out is None:
        raise FaltanDatos(
            f"'{xml_path}' no contiene una sección <output>. "
            "El cálculo probablemente no terminó."
        )

    res = QEResult(xml_path=str(xml_path))

    # --- prefijo y tipo de cálculo (de la sección de entrada) ---
    inp = _child(root, "input")
    if inp is not None:
        ctrl = _child(inp, "control_variables")
        if ctrl is not None:
            res.prefix = _text(ctrl, "prefix", "") or ""
            res.calculation = _text(ctrl, "calculation", "") or ""

    # --- estructura ---
    ast = _child(out, "atomic_structure")
    res.alat = float(ast.attrib.get("alat", 0.0))
    cell_el = _child(ast, "cell")
    cell_bohr = np.array(
        [_floats(_child(cell_el, v)) for v in ("a1", "a2", "a3")]
    )
    res.cell = cell_bohr * BOHR_ANG
    # recíproca en Å⁻¹ con factor 2*pi (filas = b1, b2, b3)
    res.reciprocal = 2.0 * np.pi * np.linalg.inv(res.cell).T

    apos = _child(ast, "atomic_positions")
    if apos is not None:
        atoms = _children(apos, "atom")
        res.symbols = [a.attrib.get("name", "X") for a in atoms]
        res.positions = np.array([_floats(a) for a in atoms]) * BOHR_ANG

    # --- energía total ---
    tot = _child(out, "total_energy")
    if tot is not None:
        etot = _text(tot, "etot")
        if etot:
            res.total_energy = float(etot) * HARTREE_EV

    # --- volumen ---
    res.volume = float(abs(np.linalg.det(res.cell)))

    # --- tensor de esfuerzos ---
    # pw.x solo lo escribe si el input pidió tstress=.true.; viene en
    # Ha/bohr^3 y en orden Fortran (por columnas).
    stress_el = _child(out, "stress")
    if stress_el is not None and stress_el.text:
        vals = _floats(stress_el)
        if vals.size == 9:
            order = stress_el.attrib.get("order", "F").upper()
            mat = vals.reshape((3, 3), order="F" if order.startswith("F") else "C")
            res.stress = mat * HA_BOHR3_GPA
            res.pressure = float(np.trace(res.stress) / 3.0)

    # --- estructura de bandas ---
    bs = _child(out, "band_structure")
    if bs is None:
        raise FaltanDatos(f"'{xml_path}' no contiene <band_structure>.")

    res.nbnd = int(float(_text(bs, "nbnd", "0")))
    res.nelec = float(_text(bs, "nelec", "0"))
    lsda = (_text(bs, "lsda", "false") or "false").lower() == "true"
    res.noncolin = (_text(bs, "noncolin", "false") or "false").lower() == "true"
    res.nspin = 2 if lsda else 1
    res.occupations_kind = _text(bs, "occupations_kind", "") or ""

    for attr, tag in (
        ("fermi", "fermi_energy"),
        ("homo", "highestOccupiedLevel"),
        ("lumo", "lowestUnoccupiedLevel"),
    ):
        val = _text(bs, tag)
        if val:
            setattr(res, attr, float(val.split()[0]) * HARTREE_EV)

    # con lsda, nbnd puede venir como nbnd_up/nbnd_dw
    if lsda and res.nbnd == 0:
        nup = _text(bs, "nbnd_up")
        if nup:
            res.nbnd = int(float(nup))

    ks_list = _children(bs, "ks_energies")
    if not ks_list:
        raise FaltanDatos(f"'{xml_path}' no contiene puntos k con eigenvalores.")

    kcart_2pi_alat = []
    weights = []
    eig_raw = []
    occ_raw = []
    for ks in ks_list:
        kp = _child(ks, "k_point")
        kcart_2pi_alat.append(_floats(kp))
        weights.append(float(kp.attrib.get("weight", 1.0)))
        eig_raw.append(_floats(_child(ks, "eigenvalues")))
        occ_el = _child(ks, "occupations")
        occ_raw.append(_floats(occ_el) if occ_el is not None else None)

    kcart_2pi_alat = np.array(kcart_2pi_alat)
    res.weights = np.array(weights)

    # cartesianas: de 2*pi/alat (bohr⁻¹) a Å⁻¹
    res.kpoints_cart = kcart_2pi_alat * (2.0 * np.pi / (res.alat * BOHR_ANG))
    # fraccionarias respecto a la recíproca
    res.kpoints_frac = res.kpoints_cart @ np.linalg.inv(res.reciprocal)

    nk = len(ks_list)
    total_per_k = len(eig_raw[0])
    if res.nspin == 2:
        nbnd = total_per_k // 2
    else:
        nbnd = total_per_k
    res.nbnd = nbnd

    eig = np.zeros((res.nspin, nk, nbnd))
    occ = np.zeros((res.nspin, nk, nbnd))
    for ik in range(nk):
        vals = eig_raw[ik]
        ovals = occ_raw[ik] if occ_raw[ik] is not None else np.zeros_like(vals)
        if res.nspin == 2:
            eig[0, ik] = vals[:nbnd]
            eig[1, ik] = vals[nbnd: 2 * nbnd]
            occ[0, ik] = ovals[:nbnd]
            occ[1, ik] = ovals[nbnd: 2 * nbnd]
        else:
            eig[0, ik] = vals[:nbnd]
            occ[0, ik] = ovals[:nbnd]
    res.eigenvalues = eig * HARTREE_EV
    res.occupations = occ

    _read_diagnostics(root, res)
    return res


def _read_diagnostics(root, res) -> None:
    """Lee del XML lo que sirve para diagnosticar y auditar el cálculo.

    Todo esto ya estaba ahí desde siempre: convergencia, fuerzas, la
    descomposición de la energía, los pseudos, los cutoffs y las
    operaciones de simetría que QE realmente usó. Sin ello no se puede
    decir si un cálculo sirve, ni si dos cálculos son comparables entre sí.
    """
    out = _child(root, "output")
    inp = _child(root, "input")

    # --- convergencia ---
    for parent in (root, out):
        if parent is None:
            continue
        conv = _child(parent, "convergence_info")
        if conv is None:
            continue
        scf = _child(conv, "scf_conv")
        if scf is not None:
            val = _text(scf, "convergence_achieved")
            if val is not None:
                res.converged = str(val).strip().lower() in ("true", "t", "1")
            n = _text(scf, "n_scf_steps")
            if n:
                res.n_scf_steps = int(float(n))
            e = _text(scf, "scf_error")
            if e:
                # el XML lo escribe en Hartree; se guarda en Ry, que es la
                # unidad en que el usuario escribe conv_thr en el input
                res.scf_error = float(e) * 2.0
        opt = _child(conv, "opt_conv")
        if opt is not None:
            n = _text(opt, "n_opt_steps")
            if n:
                res.n_bfgs_steps = int(float(n))

    # --- fuerzas (Ha/bohr en el XML) ---
    if out is not None:
        f = _child(out, "forces")
        if f is not None:
            arr = _floats(f)
            if arr.size and arr.size % 3 == 0:
                res.forces = arr.reshape(-1, 3) * (HARTREE_EV / BOHR_ANG)

        # --- descomposición de la energía ---
        te = _child(out, "total_energy")
        if te is not None:
            for tag in ("etot", "eband", "ehart", "vtxc", "etxc", "ewald",
                        "demet", "efieldcorr", "potentiostat_contr"):
                v = _text(te, tag)
                if v:
                    try:
                        res.energy_terms[tag] = float(v) * HARTREE_EV
                    except ValueError:
                        pass

        # --- magnetización ---
        mag = _child(out, "magnetization")
        if mag is not None:
            for tag, attr in (("total", "total_magnetization"),
                              ("absolute", "absolute_magnetization")):
                v = _text(mag, tag)
                if v:
                    try:
                        setattr(res, attr, float(v.split()[0]))
                    except ValueError:
                        pass

    # --- funcional, pseudos y cutoffs ---
    for parent in (out, inp):
        if parent is None:
            continue
        dft = _child(parent, "dft")
        if dft is not None and not res.functional:
            fx = _text(dft, "functional")
            if fx:
                res.functional = fx.strip()
        basis = _child(parent, "basis_set")
        if basis is not None:
            if res.ecutwfc is None:
                v = _text(basis, "ecutwfc")
                if v:
                    res.ecutwfc = float(v) * 2.0        # Ha -> Ry
            if res.ecutrho is None:
                v = _text(basis, "ecutrho")
                if v:
                    res.ecutrho = float(v) * 2.0
        species = _child(parent, "atomic_species")
        if species is not None and not res.pseudo_files:
            for sp in _children(species, "species"):
                nombre = sp.attrib.get("name", "?")
                upf = _text(sp, "pseudo_file")
                if upf:
                    res.pseudo_files[nombre] = upf.strip()

    # --- malla k y smearing (del bloque de entrada) ---
    if inp is not None:
        kp = _child(inp, "k_points_IBZ")
        if kp is not None:
            mesh = _child(kp, "monkhorst_pack")
            if mesh is not None:
                a = mesh.attrib
                try:
                    res.kgrid = (int(a.get("nk1", 0)), int(a.get("nk2", 0)),
                                 int(a.get("nk3", 0)))
                    res.kshift = (int(a.get("k1", 0)), int(a.get("k2", 0)),
                                  int(a.get("k3", 0)))
                except ValueError:
                    pass
        bs_in = _child(inp, "bands")
        if bs_in is not None:
            sm = _child(bs_in, "smearing")
            if sm is not None:
                res.smearing = (sm.text or "").strip()
                try:
                    res.degauss = float(sm.attrib.get("degauss", 0.0)) * 2.0
                except ValueError:
                    pass

    # --- simetría que QE realmente usó ---
    for parent in (out, root):
        if parent is None:
            continue
        sym = _child(parent, "symmetries")
        if sym is None:
            continue
        n = _text(sym, "nsym")
        if n:
            res.n_sym = int(float(n))
        ops = _children(sym, "symmetry")
        if ops and res.equivalent_atoms is None:
            eq = _child(ops[0], "equivalent_atoms")
            if eq is not None:
                res.equivalent_atoms = _floats(eq).astype(int)
        break

    # --- tiempos ---
    ti = _child(root, "timing_info")
    if ti is not None:
        tot = _child(ti, "total")
        if tot is not None:
            for tag, attr in (("cpu", "cpu_time"), ("wall", "wall_time")):
                v = _text(tot, tag)
                if v:
                    try:
                        setattr(res, attr, float(v))
                    except ValueError:
                        pass


# ----------------------------------------------------------------------
# Etiquetas de puntos de alta simetría
# ----------------------------------------------------------------------
def read_kpath_labels(path: str) -> list:
    """Lee KPATH.txt (generado por `olla-dft gen`) -> [(etiqueta, frac), ...]."""
    labels = []
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split()
            if len(parts) < 5:
                continue
            try:
                coords = np.array([float(x) for x in parts[-3:]])
            except ValueError:
                continue
            labels.append((parts[1], coords))
    return labels


def read_crystal_b_card(path: str) -> list:
    """Extrae los puntos de una tarjeta K_POINTS crystal_b de un input de pw.x.

    Devuelve [(etiqueta_o_None, frac), ...]. Sirve como respaldo cuando no
    existe KPATH.txt: acepta las etiquetas en comentario al final de la
    línea (`0.0 0.0 0.0 20 ! G`), como suelen escribirse a mano.
    """
    points = []
    lines = Path(path).read_text().splitlines()
    start = None
    for i, line in enumerate(lines):
        if line.strip().upper().startswith("K_POINTS") and "crystal_b" in line.lower():
            start = i
            break
    if start is None:
        return []
    try:
        npts = int(lines[start + 1].split()[0])
    except (IndexError, ValueError):
        return []
    for line in lines[start + 2: start + 2 + npts]:
        label = None
        if "!" in line:
            line, comment = line.split("!", 1)
            label = comment.strip().split()[0] if comment.strip() else None
        parts = line.split()
        if len(parts) < 3:
            continue
        try:
            coords = np.array([float(x) for x in parts[:3]])
        except ValueError:
            continue
        points.append((label, coords))
    return points


def match_labels_to_kpoints(
    kpoints_frac: np.ndarray, labels: list, tol: float = 1e-3
) -> list:
    """Ubica cada punto especial dentro de la lista de k-points calculados.

    Recorre el camino en orden y busca hacia adelante, de modo que las
    etiquetas repetidas (Γ aparece varias veces) caigan en el punto correcto
    y no siempre en el primero.

    Devuelve [(índice, etiqueta), ...] ordenado por índice.
    """
    matched = []
    cursor = 0
    nk = len(kpoints_frac)
    for label, frac in labels:
        # tolerar traslaciones por vectores de la red recíproca
        diff = kpoints_frac[cursor:] - frac
        diff -= np.round(diff)
        dist = np.linalg.norm(diff, axis=1)
        if len(dist) == 0:
            break
        idx_rel = int(np.argmin(dist))
        if dist[idx_rel] > tol:
            continue
        idx = cursor + idx_rel
        matched.append((idx, label))
        cursor = min(idx, nk - 1)
    return matched
