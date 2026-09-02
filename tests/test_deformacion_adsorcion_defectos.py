"""Barrido de deformación, adsorción, constantes elásticas 2D, defectos
cargados y las banderas de física de 'gen'.

Barrido de deformación, sitios de adsorción, constantes de lámina,
correcciones de Madelung y energías de formación de defectos, y las
banderas de física que 'gen' expone (dipolo, MD, malla k de losas, U).
"""

import numpy as np
import pytest
from ase.build import bulk, fcc111, fcc100, bcc110, hcp0001, molecule

from qekit.core.errors import ErrorDeUso
from qekit.modules import adsorb, elastic, inputgen, strain


# ----------------------------------------------------------------------
# strain
# ----------------------------------------------------------------------
def test_rango_en_porciento():
    v = strain.rango("-5:5:11")
    assert len(v) == 11
    assert v[0] == pytest.approx(-0.05)
    assert v[-1] == pytest.approx(0.05)


def test_rango_mete_el_cero_si_falta():
    v = strain.rango("-4:4:4")
    assert any(abs(x) < 1e-12 for x in v), "el punto de referencia es obligatorio"


@pytest.mark.parametrize("malo", ["-5:5:2", "5:-5:11", "-50:50:11", "hola", "-5:5"])
def test_rango_rechaza_lo_que_no_es(malo):
    with pytest.raises(ErrorDeUso):
        strain.rango(malo)


def test_matriz_biaxial_no_toca_c():
    e = strain.matriz("biaxial", 0.03)
    assert e[0, 0] == pytest.approx(0.03)
    assert e[1, 1] == pytest.approx(0.03)
    assert e[2, 2] == pytest.approx(0.0)


def test_matriz_cizalla_reparte_la_mitad():
    e = strain.matriz("cizalla", 0.02)
    assert e[0, 1] == pytest.approx(0.01)
    assert e[1, 0] == pytest.approx(0.01)


def test_modo_desconocido():
    with pytest.raises(ErrorDeUso):
        strain.matriz("diagonal", 0.01)


def test_potencial_de_deformacion_recupera_la_pendiente():
    run = strain.StrainRun(natoms=2)
    run.strains = [-0.02, -0.01, 0.0, 0.01, 0.02]
    run.energies = [0.0] * 5
    run.gaps = [1.0 + 3.5 * e for e in run.strains]
    m, r2 = strain.potencial_deformacion(run)
    assert m == pytest.approx(3.5, abs=1e-6)
    assert r2 > 0.999


def test_minimo_de_la_parabola():
    run = strain.StrainRun(natoms=1)
    run.strains = [-0.02, -0.01, 0.0, 0.01, 0.02]
    # minimo desplazado a +0.005
    run.energies = [100.0 * (e - 0.005) ** 2 for e in run.strains]
    run.gaps = [None] * 5
    eps, _ = strain.minimo(run)
    assert eps == pytest.approx(0.005, abs=1e-6)


def test_cierre_de_gap_se_interpola():
    run = strain.StrainRun(natoms=1)
    run.strains = [0.0, 0.01, 0.02]
    run.energies = [0.0, 0.0, 0.0]
    run.gaps = [0.4, 0.2, 0.0]
    assert strain.cierre_de_gap(run) is not None


def test_modulo_biaxial_recupera_la_curvatura():
    """E = E0 + (1/2)·A0·M·eps^2  =>  el modulo devuelto tiene que ser M."""
    run = strain.StrainRun(modo="biaxial", natoms=1, area0=5.0, laminar=True)
    M_eV = 12.0                       # eV/A^2
    run.strains = [-0.02, -0.01, 0.0, 0.01, 0.02]
    run.energies = [0.5 * 5.0 * M_eV * e ** 2 for e in run.strains]
    run.gaps = [None] * 5
    assert strain.modulo_biaxial(run) == pytest.approx(M_eV * 16.021766, rel=1e-6)


# ----------------------------------------------------------------------
# adsorb: los sitios de cada superficie son los del libro
# ----------------------------------------------------------------------
@pytest.mark.parametrize("constructor,kwargs,esperado", [
    (fcc111, dict(symbol="Al", size=(2, 2, 4)), {"top": 1, "bridge": 1, "hollow": 2}),
    (fcc111, dict(symbol="Al", size=(3, 3, 4)), {"top": 1, "bridge": 1, "hollow": 2}),
    (fcc100, dict(symbol="Al", size=(2, 2, 4)), {"top": 1, "bridge": 1, "hollow": 1}),
    (bcc110, dict(symbol="Fe", size=(2, 2, 4)), {"top": 1, "bridge": 2, "hollow": 1}),
    (hcp0001, dict(symbol="Co", size=(2, 2, 4)), {"top": 1, "bridge": 1, "hollow": 2}),
])
def test_sitios_por_superficie(constructor, kwargs, esperado):
    sl = constructor(vacuum=10.0, **kwargs)
    sl.pbc = (True, True, True)
    cuenta = {}
    for s in adsorb.sitios(sl):
        cuenta[s.tipo] = cuenta.get(s.tipo, 0) + 1
    assert cuenta == esperado


def test_los_sitios_no_dependen_del_tamano_de_la_supercelda():
    """2x2 y 3x3 de la misma superficie tienen los mismos sitios distintos."""
    a = fcc111("Al", size=(2, 2, 4), vacuum=10.0); a.pbc = (True, True, True)
    b = fcc111("Al", size=(3, 3, 4), vacuum=10.0); b.pbc = (True, True, True)
    tipos = lambda sl: sorted(s.tipo for s in adsorb.sitios(sl))
    assert tipos(a) == tipos(b)


def test_colocar_pone_el_ancla_a_la_altura_pedida():
    sl = fcc111("Al", size=(2, 2, 3), vacuum=8.0); sl.pbc = (True, True, True)
    s = adsorb.sitios(sl)[0]
    mol = molecule("CO")
    fuera = adsorb.colocar(sl, mol, s, altura=1.8, ancla=0)
    assert len(fuera) == len(sl) + len(mol)
    z_ancla = fuera.get_positions()[len(sl)][2]
    assert z_ancla == pytest.approx(s.z + 1.8, abs=1e-6)


def test_molecula_inexistente_da_error_de_uso():
    with pytest.raises(ErrorDeUso):
        adsorb.cargar_molecula("unobtanio")


def test_adsorcion_sin_vacio_se_niega():
    with pytest.raises(ErrorDeUso):
        adsorb.prepare(bulk("Al", "fcc", 4.05), "CO", outdir="/tmp/no_deberia")


def test_energia_de_adsorcion_es_la_resta():
    run = adsorb.AdsorbRun(natoms_slab=4, molecula="CO")
    run.sitios = [adsorb.Sitio("top", (0.0, 0.0), 0.0, etiqueta="top1")]
    run.E_slab, run.E_mol = -100.0, -20.0
    run.energies = [-121.5]
    assert run.energias_ads[0] == pytest.approx(-1.5)


# ----------------------------------------------------------------------
# elastic 2D
# ----------------------------------------------------------------------
def test_constantes_2d_no_dependen_del_vacio():
    """C_3D va como 1/c, asi que C_3D * c tiene que salir igual con dos vacios."""
    C_a = np.zeros((6, 6)); C_a[0, 0] = 100.0        # GPa con c = 20 A
    C_b = np.zeros((6, 6)); C_b[0, 0] = 50.0         # el mismo material con c = 40 A
    assert (elastic.constantes_2d(C_a, 20.0)[0, 0]
            == pytest.approx(elastic.constantes_2d(C_b, 40.0)[0, 0]))


def test_conversion_gpa_a_nm():
    C = np.zeros((6, 6)); C[0, 0] = 1.0              # 1 GPa
    assert elastic.constantes_2d(C, 10.0)[0, 0] == pytest.approx(1.0)  # 1 GPa*10A = 1 N/m


def test_born_2d_acepta_una_lamina_estable():
    C2 = np.array([[350.0, 60.0, 0.0], [60.0, 350.0, 0.0], [0.0, 0.0, 145.0]])
    estable, fallan = elastic.born_2d(C2)
    assert estable and not fallan


def test_born_2d_detecta_c12_demasiado_grande():
    C2 = np.array([[100.0, 200.0, 0.0], [200.0, 100.0, 0.0], [0.0, 0.0, 50.0]])
    estable, fallan = elastic.born_2d(C2)
    assert not estable and any("C12" in f for f in fallan)


def test_modulos_2d_de_una_lamina_isotropa():
    c11, c12 = 352.0, 60.0
    C2 = np.array([[c11, c12, 0.0], [c12, c11, 0.0],
                   [0.0, 0.0, (c11 - c12) / 2]])
    m = elastic.modulos_2d(C2)
    assert m["Y_x"] == pytest.approx(c11 - c12 ** 2 / c11, rel=1e-9)
    assert m["nu_x"] == pytest.approx(c12 / c11, rel=1e-9)
    assert m["K"] == pytest.approx((c11 + c12) / 2.0, rel=1e-9)


def test_2d_sobre_bulto_se_niega():
    with pytest.raises(ErrorDeUso):
        elastic.prepare(bulk("Si", "diamond", 5.43), outdir="/tmp/no_deberia",
                        dosd=True)


# ----------------------------------------------------------------------
# banderas de gen
# ----------------------------------------------------------------------
def test_dipolo_necesita_vacio():
    with pytest.raises(ErrorDeUso):
        inputgen._region_vacio(bulk("Si", "diamond", 5.43), 3)


def test_dipolo_pone_la_sierra_en_el_vacio():
    sl = fcc111("Al", size=(1, 1, 4), vacuum=10.0); sl.pbc = (True, True, True)
    emaxpos, eopreg = inputgen._region_vacio(sl, 3)
    assert 0.0 <= emaxpos <= 1.0
    assert 0.02 <= eopreg <= 0.2
    # la sierra no puede caer donde hay atomos
    z = sorted(f % 1.0 for f in sl.get_scaled_positions()[:, 2])
    assert all(abs(((emaxpos - f + 0.5) % 1.0) - 0.5) > eopreg / 2 for f in z)


def test_md_convierte_el_paso_a_unidades_atomicas():
    sl = fcc111("Al", size=(1, 1, 4), vacuum=10.0); sl.pbc = (True, True, True)
    txt = inputgen.build_pw_input(
        atoms=sl, pseudos={"Al": {"filename": "Al.UPF", "found": True}},
        calculation="md", prefix="al", pseudo_dir=".", ecutwfc=30, ecutrho=240,
        kcard="K_POINTS gamma\n", md=dict(dt_fs=1.0, nstep=500,
                                          thermostat="berendsen",
                                          temperature=300.0),
        nosym=True)
    assert "dt               = 20.6706" in txt or "dt " in txt
    dt = float([l for l in txt.splitlines() if l.strip().startswith("dt ")][0]
               .split("=")[1])
    assert dt == pytest.approx(1.0 / 4.8378e-2, rel=1e-3)
    assert "ion_dynamics     = 'verlet'" in txt
    assert "nosym            = .true." in txt


def test_md_fuerza_nosym_por_el_preset():
    """La MD rompe la simetria al primer paso; pw.x aborta si arranco con ella."""
    sl = fcc111("Al", size=(2, 2, 3), vacuum=8.0); sl.pbc = (True, True, True)
    opts = inputgen.GenOptions(preset="md", outdir="/tmp/qekit_md_test")
    assert opts.nosym is False          # no viene forzado en las opciones...
    # ...lo fuerza generate() por ser preset md
    import inspect
    src = inspect.getsource(inputgen.generate)
    assert 'nosym=opts.nosym or preset in ("md",)' in src


def test_kmesh_no_gasta_puntos_en_el_vacio():
    from qekit.core import kpoints as kp
    sl = fcc111("Al", size=(2, 2, 4), vacuum=12.0); sl.pbc = (True, True, True)
    n = kp.kgrid_from_spacing(sl, 0.20)
    assert n[2] == 1, "una direccion que solo tiene vacio no necesita mas de un k"
    assert n[0] > 1 and n[1] > 1


def test_kmesh_de_bulto_no_se_toca():
    from qekit.core import kpoints as kp
    n = kp.kgrid_from_spacing(bulk("Si", "diamond", 5.43), 0.20)
    assert all(v > 1 for v in n)


def test_hubbard_se_lee_bien():
    from qekit.cli import _parse_hubbard
    assert _parse_hubbard(["Ni=4.6", "Fe=3.0"]) == {"Ni": 4.6, "Fe": 3.0}
    assert _parse_hubbard(["Ni=4.6,Fe=3"]) == {"Ni": 4.6, "Fe": 3.0}
    assert _parse_hubbard(None) is None
    with pytest.raises(ErrorDeUso):
        _parse_hubbard(["Ni:4.6"])
    with pytest.raises(ErrorDeUso):
        _parse_hubbard(["Ni=mucho"])


def test_valores_negativos_se_pegan_a_su_bandera():
    from qekit.cli import _pegar_negativos
    assert _pegar_negativos(["strain", "a.cif", "-r", "-5:5:11"]) == \
        ["strain", "a.cif", "-r=-5:5:11"]
    # no debe tocar una bandera normal seguida de otra bandera
    assert _pegar_negativos(["gen", "a.cif", "-r", "--run"]) == \
        ["gen", "a.cif", "-r", "--run"]


# ----------------------------------------------------------------------
# defectos cargados
# ----------------------------------------------------------------------
from qekit.modules import defects


def test_madelung_de_la_cubica_simple():
    """Valor de libro: alpha = 2.8372974 para la red cubica simple."""
    for L in (1.0, 5.0, 13.7):
        assert defects.constante_madelung(np.eye(3) * L) == pytest.approx(
            2.8372974, abs=1e-5)


def test_madelung_no_depende_de_la_escala():
    a = defects.constante_madelung(np.eye(3) * 3.0)
    b = defects.constante_madelung(np.eye(3) * 30.0)
    assert a == pytest.approx(b, rel=1e-7)


def test_madelung_fcc_y_bcc():
    a = 5.0
    fcc = np.array([[0, a / 2, a / 2], [a / 2, 0, a / 2], [a / 2, a / 2, 0]])
    bcc = np.array([[-a / 2, a / 2, a / 2], [a / 2, -a / 2, a / 2],
                    [a / 2, a / 2, -a / 2]])
    assert defects.constante_madelung(fcc) == pytest.approx(2.8883, abs=2e-3)
    assert defects.constante_madelung(bcc) == pytest.approx(2.8883, abs=2e-3)


def test_xi_es_negativo():
    """La carga se estabiliza con su propio fondo: la energia es negativa."""
    assert defects.madelung_xi(np.eye(3) * 10.0) < 0


def test_correccion_va_como_q_cuadrado():
    cell = np.eye(3) * 10.0
    c1 = defects.correccion_imagen(1, cell, 5.0, "makov-payne")["E_corr"]
    c2 = defects.correccion_imagen(2, cell, 5.0, "makov-payne")["E_corr"]
    assert c2 == pytest.approx(4.0 * c1, rel=1e-9)


def test_correccion_va_como_uno_sobre_epsilon():
    cell = np.eye(3) * 10.0
    a = defects.correccion_imagen(1, cell, 2.0, "makov-payne")["E_corr"]
    b = defects.correccion_imagen(1, cell, 8.0, "makov-payne")["E_corr"]
    assert a == pytest.approx(4.0 * b, rel=1e-9)


def test_correccion_va_como_uno_sobre_L():
    a = defects.correccion_imagen(1, np.eye(3) * 10.0, 5.0, "makov-payne")["E_corr"]
    b = defects.correccion_imagen(1, np.eye(3) * 20.0, 5.0, "makov-payne")["E_corr"]
    assert a == pytest.approx(2.0 * b, rel=1e-9)


def test_correccion_contra_la_formula_a_mano():
    """E = ke*q^2*alpha/(2*eps*L), con ke = 14.399645 eV*A."""
    L, eps, q = 10.86, 11.7, 1
    esperado = 14.399645 * q ** 2 * 2.8372974 / (2 * eps * L)
    got = defects.correccion_imagen(q, np.eye(3) * L, eps, "makov-payne")
    assert got["E_corr"] == pytest.approx(esperado, rel=1e-4)


def test_carga_cero_no_lleva_correccion():
    c = defects.correccion_imagen(0, np.eye(3) * 10.0, 11.7)
    assert c["E_corr"] == 0.0


def test_correccion_sin_epsilon_se_niega():
    with pytest.raises(ErrorDeUso):
        defects.correccion_imagen(1, np.eye(3) * 10.0, None)


def test_esquema_desconocido():
    with pytest.raises(ErrorDeUso):
        defects.correccion_imagen(1, np.eye(3) * 10.0, 5.0, "magia")


def _run_de_prueba():
    """Un DefectRun sintetico con numeros redondos, para revisar la formula."""
    run = defects.DefectRun(kind="vacancy", cargas=[-1, 0, 1],
                            cell=np.eye(3) * 10.0, epsilon=1e9,
                            esquema="makov-payne", natoms_perf=8)
    run.n_especies = {"Si": -1}
    run.E_perfecto = -100.0
    run.vbm = 5.0
    run.gap = 1.0
    run.mu = {"Si": -12.0}
    run.energies = {-1: -83.0, 0: -88.0, 1: -94.0}
    run.converged = {q: True for q in run.cargas}
    return run


def test_formula_de_la_energia_de_formacion():
    """E_f = E_def - E_perf - sum(n_i mu_i) + q(vbm + eF) + E_corr."""
    run = _run_de_prueba()
    # con epsilon enorme la correccion es despreciable
    assert run.E_f(0, 0.0) == pytest.approx(-88.0 + 100.0 - (-1 * -12.0), abs=1e-6)
    # q=+1 anade vbm + eF
    assert run.E_f(1, 0.0) == pytest.approx(
        -94.0 + 100.0 - 12.0 + 1 * 5.0, abs=1e-4)
    assert run.E_f(1, 0.4) == pytest.approx(run.E_f(1, 0.0) + 0.4, abs=1e-9)
    assert run.E_f(-1, 0.4) == pytest.approx(run.E_f(-1, 0.0) - 0.4, abs=1e-9)


def test_la_pendiente_de_E_f_es_la_carga():
    run = _run_de_prueba()
    for q in run.cargas:
        pendiente = (run.E_f(q, 1.0) - run.E_f(q, 0.0)) / 1.0
        assert pendiente == pytest.approx(q, abs=1e-9)


def test_nivel_de_transicion_es_donde_se_cruzan():
    run = _run_de_prueba()
    for t in defects.niveles_transicion(run):
        e1 = run.E_f(t["q1"], t["eps"])
        e2 = run.E_f(t["q2"], t["eps"])
        assert e1 == pytest.approx(e2, abs=1e-8), "en el nivel, las dos E_f coinciden"


def test_envolvente_es_el_minimo():
    run = _run_de_prueba()
    ef = np.linspace(0, 1.0, 21)
    env, qs = defects.envolvente(run, ef)
    for k, e in enumerate(ef):
        todas = [run.E_f(q, e) for q in run.cargas]
        assert env[k] == pytest.approx(min(todas), abs=1e-9)
        assert run.E_f(int(qs[k]), e) == pytest.approx(env[k], abs=1e-9)


def test_mu_elemental_solo_en_cristal_de_una_especie():
    run = _run_de_prueba()
    run.mu = {}
    assert defects.asignar_mu_elemental(run, ["Si"] * 8)
    assert run.mu["Si"] == pytest.approx(-100.0 / 8)
    # en un compuesto NO se inventa
    run2 = _run_de_prueba(); run2.mu = {}
    assert not defects.asignar_mu_elemental(run2, ["Si"] * 4 + ["C"] * 4)


def test_termino_mu_de_cada_tipo_de_defecto():
    """Vacancia suma mu, intersticial lo resta, sustitucion hace las dos."""
    base = _run_de_prueba()
    base.mu = {"Si": -12.0, "P": -7.0}
    base.n_especies = {"Si": -1}                       # vacancia
    vac = base.E_f(0, 0.0)
    base.n_especies = {"P": +1}                        # intersticial
    inter = base.E_f(0, 0.0)
    base.n_especies = {"Si": -1, "P": +1}              # sustitucion
    sust = base.E_f(0, 0.0)
    assert vac == pytest.approx(-88.0 + 100.0 + (-12.0) * 1, abs=1e-6)
    assert inter == pytest.approx(-88.0 + 100.0 + 7.0, abs=1e-6)
    assert sust == pytest.approx(vac + (inter - (-88.0 + 100.0)), abs=1e-6)


def test_alineamiento_entra_multiplicado_por_q():
    run = _run_de_prueba()
    sin_dv = {q: run.E_f(q, 0.0) for q in run.cargas}
    run.dV = 0.25
    for q in run.cargas:
        assert run.E_f(q, 0.0) == pytest.approx(sin_dv[q] + q * 0.25, abs=1e-9)
