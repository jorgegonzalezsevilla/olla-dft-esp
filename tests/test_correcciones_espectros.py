"""Regresiones de la ronda de correcciones de unidades, banderas y textos.

Cada prueba fija un comportamiento que estaba mal o era engañoso:
unidades de Bader y del alineamiento de defectos, banderas que no hacían
nada (--dipole, --fix, --scissor), columnas equivocadas al comparar
espectros, y textos de ayuda que prometían otra cosa.
"""

import re
from pathlib import Path

import numpy as np
import pytest
from ase import Atoms
from ase.build import bulk

from qekit.cli import build_parser, main
from qekit.core.errors import ErrorDeUso


@pytest.fixture(autouse=True)
def _inputs_se_escriben():
    """main() con --collect apaga la escritura de inputs para todo el
    proceso; cada prueba de aquí arranca (y termina) con ella encendida."""
    from qekit.modules import sweep
    sweep.set_write_inputs(True)
    yield
    sweep.set_write_inputs(True)


# ----------------------------------------------------------------------
# Bader en unidades de pp.x (e/bohr³) y Z_valencia desde los UPF
# ----------------------------------------------------------------------
def _upf_minimo(carpeta, elemento, zval):
    f = Path(carpeta) / f"{elemento}.pbe.UPF"
    f.write_text(
        '<UPF version="2.0.1">\n'
        f'  <PP_HEADER element="{elemento}" pseudo_type="NC" '
        f'relativistic="scalar" functional="PBE" z_valence="{zval:.6e}" '
        'wfc_cutoff="4.000000e+01" rho_cutoff="1.600000e+02" '
        'mesh_size="1000"/>\n</UPF>\n')
    return f


def _cube_dos_gaussianas(path, L=12.0, n=36, cargas=(3.0, 5.0)):
    """Un .cube como el de pp.x: rejilla en bohr, densidad en e/bohr³."""
    from qekit.modules import fields
    g = (np.arange(n) + 0.5) * L / n
    X, Y, Z = np.meshgrid(g, g, g, indexing="ij")
    dv_bohr3 = (L / n) ** 3 / fields.BOHR ** 3
    centros = [(3.0, 6.0, 6.0), (9.0, 6.0, 6.0)]
    rho = np.zeros((n, n, n))
    for c, q in zip(centros, cargas):
        d2 = sum(((A - ci + L / 2) % L - L / 2) ** 2
                 for A, ci in zip((X, Y, Z), c))
        r = np.exp(-d2 / (2 * 0.9 ** 2))
        rho += q * r / (r.sum() * dv_bohr3)
    paso = L / n / fields.BOHR
    lineas = ["densidad sintetica", "e/bohr3",
              f"{2:5d} {0.0:12.6f} {0.0:12.6f} {0.0:12.6f}"]
    for k in range(3):
        v = [0.0, 0.0, 0.0]
        v[k] = paso
        lineas.append(f"{n:5d} {v[0]:12.6f} {v[1]:12.6f} {v[2]:12.6f}")
    for c, zn in zip(centros, (3, 5)):
        cb = [x / fields.BOHR for x in c]
        lineas.append(f"{zn:5d} {0.0:12.6f} {cb[0]:12.6f} {cb[1]:12.6f} "
                      f"{cb[2]:12.6f}")
    vals = rho.ravel()
    for i in range(0, vals.size, 6):
        lineas.append(" ".join(f"{x:13.5e}" for x in vals[i:i + 6]))
    Path(path).write_text("\n".join(lineas) + "\n")
    return centros


def test_read_cube_y_bader_conservan_los_electrones(tmp_path):
    """Leer el cube de pp.x y sumar las cuencas da el número de electrones."""
    from qekit.modules import charges, fields
    cube_f = tmp_path / "rho.cube"
    centros = _cube_dos_gaussianas(cube_f)
    cube = fields.read_cube(cube_f)
    assert np.linalg.norm(cube.axes[0]) == pytest.approx(12.0 / 36, rel=1e-5)
    res = charges.bader(cube, np.array(centros), symbols=["Li", "B"],
                        valence=[3.0, 5.0])
    assert res.total_grid == pytest.approx(8.0, abs=1e-3)
    assert res.charges[0] == pytest.approx(3.0, abs=0.03)
    assert res.charges[1] == pytest.approx(5.0, abs=0.03)


def test_valencia_desde_los_upf(tmp_path):
    from qekit.modules import charges
    _upf_minimo(tmp_path, "Li", 3.0)
    _upf_minimo(tmp_path, "B", 3.0)
    v = charges.valence_from_pseudos(["Li", "B", "Li"], tmp_path)
    assert list(v) == [3.0, 3.0, 3.0]
    # sin UPF para un elemento, no se inventa nada
    assert charges.valence_from_pseudos(["Li", "Xx"], tmp_path) is None
    assert charges.valence_from_pseudos(["Li"], tmp_path / "no_existe") is None


def test_cli_charges_rellena_la_columna_neta_con_pseudo_dir(tmp_path, capsys):
    from ase.io import write as ase_write
    pdir = tmp_path / "ps"
    pdir.mkdir()
    _upf_minimo(pdir, "Li", 3.0)
    _upf_minimo(pdir, "B", 5.0)
    cube_f = tmp_path / "rho.cube"
    centros = _cube_dos_gaussianas(cube_f)
    at = Atoms("LiB", positions=centros, cell=np.eye(3) * 12.0, pbc=True)
    cif = tmp_path / "lib.cif"
    ase_write(str(cif), at)

    rc = main(["charges", str(cif), "--bader", str(cube_f),
               "--pseudo-dir", str(pdir)])
    out = capsys.readouterr().out
    assert rc == 0
    assert "n/d" not in out
    assert "Electrones de valencia según los UPF: 8.0000 e" in out
    assert "no coincide" not in out
    # una carga neta de ~0 para cada átomo (recuperan su Z_valencia)
    filas = [l for l in out.splitlines() if l.strip().startswith(("1 ", "2 "))]
    assert len(filas) == 2
    for fila in filas:
        assert abs(float(fila.split()[3])) < 0.05


# ----------------------------------------------------------------------
# Alineamiento de potencial: pp.x escribe Ry, E_f suma eV
# ----------------------------------------------------------------------
def _cube_potencial(path, perfil_ry, n_xy=4):
    """Un .cube de potencial (rejilla en bohr) con el perfil planar dado."""
    from qekit.modules import fields
    nz = len(perfil_ry)
    paso = 0.5 / fields.BOHR
    lineas = ["potencial sintetico", "Ry",
              f"{1:5d} {0.0:12.6f} {0.0:12.6f} {0.0:12.6f}",
              f"{n_xy:5d} {paso:12.6f} {0.0:12.6f} {0.0:12.6f}",
              f"{n_xy:5d} {0.0:12.6f} {paso:12.6f} {0.0:12.6f}",
              f"{nz:5d} {0.0:12.6f} {0.0:12.6f} {paso:12.6f}",
              f"{14:5d} {0.0:12.6f} {0.0:12.6f} {0.0:12.6f} {0.0:12.6f}"]
    datos = np.tile(np.asarray(perfil_ry, dtype=float), (n_xy, n_xy, 1))
    vals = datos.ravel()
    for i in range(0, vals.size, 6):
        lineas.append(" ".join(f"{x:13.5e}" for x in vals[i:i + 6]))
    Path(path).write_text("\n".join(lineas) + "\n")
    return path


def test_alineamiento_devuelve_dv_en_ev(tmp_path):
    from qekit.core import qeout
    from qekit.modules import defects
    nz = 40
    z = np.arange(nz)
    perfecto = 0.02 * np.sin(2 * np.pi * z / nz)
    # el defecto perturba el potencial cerca de z=0 y desplaza TODO 0.1 Ry
    defecto = perfecto + 0.1 + 0.5 * np.exp(-((z - 0.0) ** 2) / 8.0)
    fd = _cube_potencial(tmp_path / "def.cube", defecto)
    fp = _cube_potencial(tmp_path / "perf.cube", perfecto)

    al = defects.alineamiento(str(fd), str(fp))
    assert al["unidades"] == "eV"
    assert al["dV"] == pytest.approx(0.1 * qeout.RY_EV, abs=1e-3)
    assert al["sigma"] < 1e-3
    # el mismo cube declarado en eV no se convierte
    al_ev = defects.alineamiento(str(fd), str(fp), unidades_cube="eV")
    assert al_ev["dV"] == pytest.approx(0.1, abs=1e-4)
    with pytest.raises(ErrorDeUso, match="unidades"):
        defects.alineamiento(str(fd), str(fp), unidades_cube="Ha")


def test_reporte_de_defectos_etiqueta_dv_en_ev():
    from qekit.modules import defects
    run = defects.DefectRun(kind="vacancy", cargas=[0, 1],
                            cell=np.eye(3) * 10.0, epsilon=11.7,
                            natoms_perf=8, supercell=(2, 2, 2))
    run.energies = {0: -100.0, 1: -95.0}
    run.converged = {0: True, 1: True}
    run.E_perfecto, run.vbm, run.gap = -101.0, 5.0, 1.1
    run.n_especies, run.mu = {"Si": -1}, {"Si": -5.0}
    run.dV, run.dV_sigma = 1.3606, 0.01
    rep = defects.report(run)
    assert "ΔV = +1.3606 eV" in rep
    assert "Ry" in rep and "q·ΔV" in rep


# ----------------------------------------------------------------------
# adsorb --dipole escribe de verdad la corrección dipolar
# ----------------------------------------------------------------------
def test_adsorb_dipole_escribe_dipfield_en_los_tres_calculos(tmp_path):
    from ase.build import fcc111
    from qekit.modules import adsorb
    sl = fcc111("Al", size=(2, 2, 3), vacuum=10.0)
    sl.pbc = (True, True, True)
    run, rep = adsorb.prepare(sl, "CO", outdir=str(tmp_path / "con"),
                              tipos=("top",), pseudo_dir=str(tmp_path),
                              dipolo=True)
    assert run.jobs
    for job in run.jobs:
        txt = (Path(job.directory) / "pw.in").read_text()
        assert re.search(r"dipfield\s*=\s*\.true\.", txt)
        assert re.search(r"tefield\s*=\s*\.true\.", txt)
        assert re.search(r"edir\s*=\s*3\b", txt)
        assert "emaxpos" in txt and "eopreg" in txt
    assert "--dipole" not in rep

    run2, rep2 = adsorb.prepare(sl, "CO", outdir=str(tmp_path / "sin"),
                                tipos=("top",), pseudo_dir=str(tmp_path))
    for job in run2.jobs:
        assert "dipfield" not in (Path(job.directory) / "pw.in").read_text()
    assert "--dipole" in rep2


# ----------------------------------------------------------------------
# surface --fix llega a ATOMIC_POSITIONS (qekit_fijo y FixAtoms)
# ----------------------------------------------------------------------
def _pw_in_de(atoms):
    from qekit.modules import inputgen
    pseudos = {s: {"filename": f"{s}.UPF", "found": True}
               for s in dict.fromkeys(atoms.get_chemical_symbols())}
    return inputgen.build_pw_input(
        atoms=atoms, pseudos=pseudos, calculation="relax", prefix="x",
        pseudo_dir=".", ecutwfc=40, ecutrho=320, kcard="K_POINTS gamma\n",
        insulator=False)


def _filas_posiciones(txt):
    bloque = txt.split("ATOMIC_POSITIONS crystal")[1].split("CELL_PARAMETERS")[0]
    return [l for l in bloque.splitlines() if l.strip()]


def test_fixatoms_pone_ceros_en_atomic_positions():
    from ase.build import fcc111
    from ase.constraints import FixAtoms
    sl = fcc111("Al", size=(1, 1, 4), vacuum=8.0)
    sl.pbc = (True, True, True)
    sl.set_constraint(FixAtoms(indices=[0, 1]))
    filas = _filas_posiciones(_pw_in_de(sl))
    assert len(filas) == 4
    assert filas[0].rstrip().endswith("0 0 0")
    assert filas[1].rstrip().endswith("0 0 0")
    assert not filas[2].rstrip().endswith("0 0 0")
    assert not filas[3].rstrip().endswith("0 0 0")
    # una máscara booleana también vale
    sl.set_constraint(FixAtoms(mask=[False, False, False, True]))
    filas = _filas_posiciones(_pw_in_de(sl))
    assert [f.rstrip().endswith("0 0 0") for f in filas] == [False] * 3 + [True]


def test_surface_fix_marca_los_atomos_y_llega_al_input(tmp_path):
    from qekit.core import structure
    from qekit.modules import builder, inputgen
    info = builder.surface(bulk("Si", "diamond", a=5.43), miller=(1, 0, 0),
                           layers=6, vacuum=15.0, fix_layers=2)
    fijos = inputgen.fixed_atoms(info.atoms)
    assert len(fijos) == info.fijados > 0
    filas = _filas_posiciones(_pw_in_de(info.atoms))
    assert sum(f.rstrip().endswith("0 0 0") for f in filas) == info.fijados
    # sin marca no aparece la columna
    assert not any(f.rstrip().endswith("0 0 0")
                   for f in _filas_posiciones(_pw_in_de(bulk("Si", "diamond",
                                                              a=5.43))))
    # POSCAR conserva los fijos al recargar; el CIF no
    poscar = tmp_path / "losa.vasp"
    structure.convert(info.atoms, str(poscar))
    assert structure.conserva_fijos(str(poscar))
    assert len(inputgen.fixed_atoms(structure.load(str(poscar)))) == info.fijados
    cif = tmp_path / "losa.cif"
    structure.convert(info.atoms, str(cif))
    assert not structure.conserva_fijos(str(cif))
    assert inputgen.fixed_atoms(structure.load(str(cif))) == set()


def test_cli_surface_fix_avisa_si_el_formato_pierde_los_fijos(tmp_path, capsys):
    from ase.io import write as ase_write
    cif_in = tmp_path / "si.cif"
    ase_write(str(cif_in), bulk("Si", "diamond", a=5.43))
    rc = main(["surface", str(cif_in), "-m", "1 0 0", "-l", "6", "--fix", "2",
               "-o", str(tmp_path / "losa.cif")])
    err = capsys.readouterr().err
    assert rc == 0
    assert "no guarda qué átomos están congelados" in err
    assert "losa.vasp" in err
    rc = main(["surface", str(cif_in), "-m", "1 0 0", "-l", "6", "--fix", "2",
               "-o", str(tmp_path / "losa.vasp")])
    err = capsys.readouterr().err
    assert rc == 0 and "congelados" not in err


# ----------------------------------------------------------------------
# tddft: --compare lee α por nombre, --broadening y --scissor llegan
# ----------------------------------------------------------------------
def _optics_dat(path, E, alpha):
    """Un OPTICS.dat con el orden de columnas real de optics.export."""
    from qekit.modules import optics
    cab = ("# funciones ópticas de prueba\n# " +
           " ".join(f"{c:>10s}" for c in optics.OPTICS_COLUMNS))
    eps1 = 10 + 0 * E
    eps2 = 0.1 * E
    R = 0.3 + 0.01 * E             # la ÚLTIMA columna, que no es alpha
    datos = np.column_stack([E, eps1, eps2, eps1 ** 0.5, 0 * E + 0.01,
                             alpha, R])
    np.savetxt(path, datos, header=cab, comments="")
    return path


def test_read_optics_dat_devuelve_alpha_por_nombre(tmp_path):
    from qekit.modules import optics
    E = np.linspace(0, 10, 11)
    alpha = 1e4 * E ** 2
    f = _optics_dat(tmp_path / "OPTICS.dat", E, alpha)
    cols = optics.read_optics_dat(f)
    assert list(cols) == list(optics.OPTICS_COLUMNS)
    assert np.allclose(cols["alpha(1/cm)"], alpha)
    assert not np.allclose(cols["R"], alpha)
    # export real -> read: ida y vuelta
    run = optics.OpticsRun(energies=E, eps1=np.full_like(E, 10.0),
                           eps2=0.1 * E)
    optics.export(run, str(tmp_path / "real"))
    cols2 = optics.read_optics_dat(tmp_path / "real" / "OPTICS.dat")
    assert np.allclose(cols2["alpha(1/cm)"], optics.derived(run)["alpha"])
    # sin encabezado de nombres se asume el orden estándar
    np.savetxt(tmp_path / "sin.dat", np.column_stack([E, E, E, E, E, alpha, E]))
    assert np.allclose(optics.read_optics_dat(tmp_path / "sin.dat")
                       ["alpha(1/cm)"], alpha)


def test_cli_tddft_compare_usa_alpha_y_no_la_ultima_columna(tmp_path,
                                                             monkeypatch):
    from qekit.modules import tddft
    E = np.linspace(0, 10, 11)
    alpha = 1e4 * E ** 2
    _optics_dat(tmp_path / "OPTICS.dat", E, alpha)
    # espectro TDDFPT sintético (turbo_spectrum con units=1 -> eV)
    e = np.linspace(0, 10, 201)
    s = np.exp(-((e - 5.0) ** 2) / 0.1)
    np.savetxt(tmp_path / "x.plot.dat", np.column_stack([e, s, s, 0 * s, 0 * s]))

    visto = {}

    def plot_falso(run, outfile, comparar=None, **kw):
        visto["comparar"] = comparar
        return []

    monkeypatch.setattr(tddft, "plot", plot_falso)
    rc = main(["tddft", "-o", str(tmp_path), "--collect",
               "--compare", str(tmp_path / "OPTICS.dat")])
    assert rc == 0
    e2, a2 = visto["comparar"]
    assert np.allclose(e2, E) and np.allclose(a2, alpha)


def test_collect_fija_el_ensanchamiento_y_el_umbral_del_exciton(tmp_path):
    from qekit.modules import tddft
    e = np.linspace(0, 10, 2001)
    s = np.exp(-((e - 5.0) ** 2) / (2 * 0.05 ** 2))
    np.savetxt(tmp_path / "x.plot.dat", np.column_stack([e, s]))
    # borde ~ 4.95 eV, gap 5.1: 0.15 eV por debajo
    run = tddft.collect(tmp_path, gap_ip=5.1, broadening=0.05)
    assert run.broadening == pytest.approx(0.05)
    assert any("excitón" in a for a in run.avisos)
    # con ensanchamiento grande la misma diferencia NO se distingue
    run2 = tddft.collect(tmp_path, gap_ip=5.1, broadening=0.5)
    assert run2.broadening == pytest.approx(0.5)
    assert not any("excitón" in a for a in run2.avisos)
    assert any("no se puede distinguir" in a or "limite" in a
               for a in run2.avisos)
    # sin valor explícito se lee de spectrum.in (epsil en Ry)
    (tmp_path / "spectrum.in").write_text(
        tddft.build_spectrum_input("x", broadening=0.3))
    run3 = tddft.collect(tmp_path, gap_ip=5.1)
    assert run3.broadening == pytest.approx(0.3, abs=2e-5)


def test_prepare_tddft_pasa_scissor_a_lanczos(tmp_path):
    from ase.build import molecule
    from qekit.modules import tddft
    mol = molecule("H2O")
    mol.set_cell(np.eye(3) * 12.0)
    mol.center()
    mol.pbc = True
    tddft.prepare(mol, outdir=str(tmp_path / "a"), pseudo_dir=str(tmp_path),
                  scissor=1.0)
    txt = (tmp_path / "a" / "lanczos.in").read_text()
    assert f"scissor = {1.0 / tddft.RY_EV:.6f}" in txt
    tddft.prepare(mol, outdir=str(tmp_path / "b"), pseudo_dir=str(tmp_path))
    assert "scissor" not in (tmp_path / "b" / "lanczos.in").read_text()
    with pytest.raises(ErrorDeUso, match="scissor"):
        tddft.prepare(mol, outdir=str(tmp_path / "c"),
                      pseudo_dir=str(tmp_path), metodo="davidson", scissor=1.0)
    with pytest.raises(ErrorDeUso, match="scissor"):
        tddft.prepare(mol, outdir=str(tmp_path / "d"),
                      pseudo_dir=str(tmp_path), scissor=-1.0)


def test_parser_tddft_tiene_scissor():
    args = build_parser().parse_args(["tddft", "x.cif", "--scissor", "0.8"])
    assert args.scissor == pytest.approx(0.8)
    assert args.broadening is None


# ----------------------------------------------------------------------
# diagnose: en un relax cada paso iónico es un ciclo SCF aparte
# ----------------------------------------------------------------------
def _relax_falso(tmp_path, ciclos, beta=0.4):
    """stdout de pw.x con varios ciclos SCF; cada uno arranca en
    'iteration #  1'. `ciclos` = [(accuracies, convergio), ...]."""
    txt = ["     Program PWSCF v.7.2 starts"]
    for accs, ok in ciclos:
        txt.append("     Self-consistent Calculation")
        for i, a in enumerate(accs, start=1):
            txt += [f"     iteration # {i:3d}     ecut=    60.00 Ry     "
                    f"beta= {beta:.2f}",
                    "     total energy              =     -15.8000000 Ry",
                    f"     estimated scf accuracy    <    {a:.8f} Ry"]
        txt.append("     convergence has been achieved in %d iterations"
                   % len(accs) if ok else
                   "     convergence NOT achieved after 20 iterations")
        if ok:
            txt.append("!    total energy              =     -15.8100000 Ry")
    f = tmp_path / "relax.out"
    f.write_text("\n".join(txt))
    return f


def test_relax_clasifica_solo_el_ultimo_ciclo_scf(tmp_path):
    from qekit.modules import diagnose
    lento = [1e-1 * (0.45 ** i) for i in range(10)]     # monótono, no llega
    bien = [1e-2, 1e-4, 1e-6, 1e-8]                      # convergió
    # dos ciclos que convergen: concatenados parecerían una oscilación
    h = diagnose.read_scf_history(_relax_falso(tmp_path, [(bien, True),
                                                          (bien, True)]))
    assert h.n_ciclos == 2
    assert h.n_iter == 4 and h.accuracy == pytest.approx(bien)
    assert h.converged is True and h.patologia == ""
    # el último ciclo es el que manda
    h2 = diagnose.read_scf_history(_relax_falso(tmp_path, [(bien, True),
                                                           (lento, False)],
                                                beta=0.2))
    assert h2.n_ciclos == 2 and h2.n_iter == 10
    assert h2.converged is False and h2.patologia == "lenta"
    # y un ciclo malo al principio no condena a un final bueno
    h3 = diagnose.read_scf_history(_relax_falso(tmp_path, [(lento, False),
                                                           (bien, True)]))
    assert h3.converged is True and h3.patologia == ""


def test_reporte_de_diagnostico_dice_cuantos_ciclos_vio(tmp_path):
    from qekit.modules import diagnose
    bien = [1e-2, 1e-4, 1e-6, 1e-8]
    f = _relax_falso(tmp_path, [(bien, True)] * 3)
    d = diagnose.Diagnosis(scf=diagnose.read_scf_history(f),
                           traj=diagnose.read_trajectory(f))
    rep = diagnose.report(d)
    assert "3 ciclos SCF" in rep and "solo el último" in rep
    # un scf normal no menciona ciclos
    d1 = diagnose.Diagnosis(scf=diagnose.read_scf_history(
        _relax_falso(tmp_path, [(bien, True)])))
    assert "ciclos SCF" not in diagnose.report(d1)


# ----------------------------------------------------------------------
# función trabajo: la planitud se mide en la región de vacío
# ----------------------------------------------------------------------
def test_planitud_de_la_funcion_trabajo_se_mide_en_el_vacio():
    """Losa entre z=5 y z=15 de una celda de 40 Å. El potencial sube
    despacio en la cola de la losa y solo es plano lejos de ella: con
    posiciones, la meseta cae en el centro del hueco (z ~ 30) y es plana;
    la ventana ciega alrededor del máximo pisa la cola si el máximo está
    cerca de la losa."""
    from qekit.core import qeout
    from qekit.modules import fields
    nz, L = 400, 40.0
    z = np.linspace(0, L, nz, endpoint=False)
    # potencial (Ry): pozo en la losa, cola que decae hacia el vacío, meseta
    # exacta en 0.5 Ry lejos de la losa, y un máximo pequeño pegado a la losa
    dist = np.minimum(np.abs(z - 5.0), np.abs(z - 15.0))
    dentro = (z > 5) & (z < 15)
    V = np.where(dentro, -1.0, 0.5 - 0.3 * np.exp(-dist / 2.0))
    V = V + np.where(dentro, 0.0, 0.15 * np.exp(-((z - 17.0) ** 2) / 0.5))
    cube = fields.CubeData(origin=np.zeros(3), axes=np.diag([1, 1, L / nz]),
                           shape=(2, 2, nz), data=np.tile(V, (2, 2, 1)),
                           natoms=2)
    pos = np.array([[0.0, 0.0, 5.0], [0.0, 0.0, 15.0]])
    wf = fields.work_function(cube, fermi_ev=-3.0, axis=2, positions=pos)
    # la meseta es el 20 % central del hueco de 30 Å: ±3 Å en torno a z=30
    assert wf.vacuum_z[0] == pytest.approx(27.0, abs=0.2)
    assert wf.vacuum_z[1] == pytest.approx(33.0, abs=0.2)
    assert wf.v_vacuum == pytest.approx(0.5 * qeout.RY_EV, abs=0.05)
    assert wf.flatness < 0.05
    assert wf.phi == pytest.approx(wf.v_vacuum + 3.0)
    assert "evaluada en z" in fields.report_wf(wf)
    # sin posiciones se usa la ventana alrededor del máximo (documentado)
    wf0 = fields.work_function(cube, fermi_ev=-3.0, axis=2)
    assert wf0.vacuum_z[0] < 20.0 < wf0.vacuum_z[1]     # centrada en z~17
    assert wf0.flatness > wf.flatness


# ----------------------------------------------------------------------
# echem: U entra frente al SHE; a pH 0 coincide con el RHE
# ----------------------------------------------------------------------
def test_echem_ph_cero_iguala_she_y_rhe_y_el_ph_convierte():
    from qekit.modules import echem
    e = echem.her(-0.33)
    # a pH 0 las dos escalas son la misma
    assert e.U_rhe(0.4, 0.0) == pytest.approx(0.4)
    assert e.dG(0.4, 0.0) == e.dG(echem.u_rhe(0.4, 0.0), 0.0)
    # a pH 7, U vs SHE + 0.0592·7 es el U vs RHE que da los mismos ΔG
    u_she, ph = 0.1, 7.0
    g_she = [g for _, g in e.dG(u_she, ph)]
    g_rhe = [g for _, g in e.dG(echem.u_rhe(u_she, ph, e.T), 0.0)]
    assert g_she == pytest.approx(g_rhe)
    assert echem.u_rhe(0.0, 7.0) == pytest.approx(0.0592 * 7, abs=2e-3)
    # el sobrepotencial lleva signo y aquí es positivo (cuesta arriba)
    assert e.sobrepotencial > 0
    e.U, e.pH = u_she, ph
    rep = echem.report(e)
    assert "V vs SHE" in rep and "V vs RHE" in rep
    assert "positivo = en U_eq" in rep


# ----------------------------------------------------------------------
# selftest de la relación de escala llama al módulo echem de verdad
# ----------------------------------------------------------------------
def test_selftest_escala_oer_usa_echem(monkeypatch):
    from qekit.modules import echem, selftest as st
    pruebas = {p.clave: p for p in st.PRUEBAS}
    # el valor real del RuO2 (3.87 - 0.77 = 3.10 eV) contra la referencia
    # 3.2 con la tolerancia RELATIVA del selftest (10 %)
    v = pruebas["escala_oer"].fn(None)
    assert v == pytest.approx(3.10, abs=1e-9)
    assert abs(v - 3.2) / 3.2 <= pruebas["escala_oer"].tolerancia
    assert pruebas["escala_eta_min"].fn(None) == pytest.approx(0.37, abs=0.02)
    # si la función de echem cambia, el selftest lo nota
    monkeypatch.setattr(echem, "escala_ooh_oh", lambda e: 99.0)
    assert pruebas["escala_oer"].fn(None) == 99.0
    e = echem.oer({"OH": 0.77, "O": 2.16, "OOH": 3.87},
                  correcciones={"OH": 0.1, "O": 0, "OOH": 0.4})
    monkeypatch.undo()
    assert echem.escala_ooh_oh(e) == pytest.approx(3.87 + 0.4 - 0.77 - 0.1)
    with pytest.raises(Exception, match="OER"):
        echem.escala_ooh_oh(echem.her(-0.1))


# ----------------------------------------------------------------------
# ayudas y mensajes: --edge, bordes M, neb sin banderas de Hubbard
# ----------------------------------------------------------------------
def test_xanes_rechaza_bordes_m_y_acepta_los_de_xspectra():
    from qekit.modules import xanes
    for b in ("K", "L1", "L2", "L3", "L23", "k"):
        assert xanes.validar_borde(b) == b.upper()
    with pytest.raises(ErrorDeUso, match="bordes M"):
        xanes.validar_borde("M45")
    with pytest.raises(ErrorDeUso, match="desconocido"):
        xanes.validar_borde("N7")
    with pytest.raises(ErrorDeUso, match="bordes M"):
        xanes.build_xspectra_input("x", 1, (1, 0, 0), "x.wfc", borde="M23")
    ayuda = build_parser().parse_args(["xanes", "x.cif"])
    assert ayuda.edge == "K"


def test_cli_xanes_borde_m_da_error_de_uso(tmp_path, capsys):
    from ase.io import write as ase_write
    cif = tmp_path / "si.cif"
    ase_write(str(cif), bulk("Si", "diamond", a=5.43))
    rc = main(["xanes", str(cif), "--element", "Si", "--edge", "M45",
               "-o", str(tmp_path / "x")])
    err = capsys.readouterr().err
    assert rc != 0 and "bordes M" in err


def test_mensajes_sugieren_edge_y_no_borde():
    import inspect
    from qekit.modules import corehole, xps
    from qekit import cli
    for mod in (corehole, xps, cli):
        assert "--borde" not in inspect.getsource(mod)
    # y el hueco que se sugiere para L2/L3 es el 2p (L23)
    from qekit.modules import xanes
    assert xanes.BORDE_COREHOLE["L3"] == "L23"


def test_neb_no_lleva_banderas_de_hubbard():
    p = build_parser()
    with pytest.raises(SystemExit):
        p.parse_args(["neb", "a.cif", "b.cif", "--intersite"])
    assert p.parse_args(["hubbard", "x.cif", "--intersite"]).intersite


def test_recomendador_distingue_uno_de_dos_casos():
    from qekit.modules import recommend as rc
    s1 = rc.Sugerencia(campo="ecutwfc", valor=50.0, n_casos=1, confianza="baja",
                       razon="x")
    s2 = rc.Sugerencia(campo="ecutwfc", valor=50.0, n_casos=2, confianza="baja",
                       razon="x")
    assert "UN SOLO CASO" in rc.report([s1], ["Si"], 1)
    t2 = rc.report([s2], ["Si"], 2)
    assert "UN SOLO CASO" not in t2 and "SOLO 2 CASOS" in t2
