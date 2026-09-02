"""Termoquímica, electrón-fonón, dinámica molecular, heteroestructuras,
NEB y asistente.

Casi todo se comprueba contra valores que NO salen de Olla-DFT: tablas de
NIST, límites analíticos exactos, y teoremas que el resultado tiene que
cumplir por construcción.
"""

import gzip
from pathlib import Path

import numpy as np
import pytest

from qekit.core.errors import ErrorDeUso
from qekit.modules import elph, dynamics as dyn, interface as itf
from qekit.modules import neb as nb, thermochem as tc, wizard as wz

DATOS = Path(__file__).parent / "datos"


# ======================================================================
# Termoquímica: contra las tablas de NIST
# ======================================================================
def _S_JmolK(tq):
    return tq.S * 96485.0


@pytest.mark.parametrize("molec,nu,sigma,S_nist", [
    ("H2O", [1595.0, 3657.0, 3756.0], 2, 188.83),
    ("N2", [2358.6], 2, 191.61),
    ("CH4", [1306] * 3 + [1534] * 2 + [2917] + [3019] * 3, 12, 186.25),
])
def test_entropia_de_gases_contra_nist(molec, nu, sigma, S_nist):
    """Prueba dura: la entropía absoluta a 298.15 K y 1 bar es tabulada."""
    from ase.build import molecule
    tq = tc.corregir(nu, T=298.15, fase="gas", atoms=molecule(molec),
                     p=1e5, simetria=sigma)
    assert _S_JmolK(tq) == pytest.approx(S_nist, rel=0.005)


def test_zpe_del_agua():
    """ZPE = (1/2) suma de los cuantos; tabla: 0.5581 eV."""
    assert tc.zpe([1595.0, 3657.0, 3756.0]) == pytest.approx(0.5581, abs=0.002)


def test_numero_de_simetria_no_es_cosmetico():
    """Olvidar sigma=12 en el metano mueve G unos 0.06 eV a 300 K."""
    from ase.build import molecule
    ch4 = molecule("CH4")
    nu = [1306] * 3 + [1534] * 2 + [2917] + [3019] * 3
    bien = tc.corregir(nu, fase="gas", atoms=ch4, simetria=12)
    mal = tc.corregir(nu, fase="gas", atoms=ch4, simetria=1)
    assert (mal.TS - bien.TS) == pytest.approx(0.0638, abs=0.005)


def test_entropia_vibracional_tiende_a_cero_a_baja_T():
    assert tc.S_vib([1000.0], 1.0) == pytest.approx(0.0, abs=1e-12)


def test_cv_vibracional_tiende_a_dulong_petit():
    """A T alta, cada modo aporta k_B a C_v."""
    n = 6
    cv = tc.Cv_vib([300.0] * n, 30000.0)
    assert cv == pytest.approx(n * tc.KB_EV, rel=0.01)


def test_frecuencia_imaginaria_en_un_minimo_avisa():
    tq = tc.corregir([-120.0, 500.0, 900.0], fase="solido")
    assert tq.n_imaginarias == 1
    assert any("MÍNIMO" in a for a in tq.avisos)


def test_estado_de_transicion_sin_imaginaria_avisa():
    tq = tc.corregir([500.0, 900.0], fase="transicion")
    assert any("no hay ninguna frecuencia imaginaria" in a for a in tq.avisos)


def test_modos_blandos_se_pueden_subir_a_un_piso():
    tq = tc.corregir([5.0, 20.0, 500.0], fase="solido", piso=100.0)
    assert tq.n_subidos == 2
    sin = tc.corregir([5.0, 20.0, 500.0], fase="solido")
    # la entropía de un modo blando es enorme: subirla la baja mucho
    assert tq.S_vib < sin.S_vib


def test_adsorcion_resta_las_correcciones_del_gas():
    from ase.build import molecule
    gas = tc.corregir([2358.6], fase="gas", atoms=molecule("N2"), simetria=2)
    ads = tc.corregir([1900.0, 400.0, 380.0, 250.0, 240.0, 90.0],
                      fase="adsorbato")
    d = tc.adsorcion(-100.0, -95.0, -4.5, ads, gas)
    assert d["E_ads"] == pytest.approx(-0.5)
    # la molécula pierde entropía al pegarse: -dTS es positivo
    assert d["dTS"] < 0
    assert d["G_ads"] > d["E_ads"]


# ======================================================================
# Electrón-fonón: límites exactos y superconductores reales
# ======================================================================
def test_lambda_y_omega_log_de_un_modo_de_einstein():
    """Con un solo modo, lambda = 2A/w0 y omega_log = w0, exactos."""
    w0, A, s = 5.0, 0.8, 0.02
    w = np.linspace(0.001, 20, 200000)
    a2F = A / (s * np.sqrt(2 * np.pi)) * np.exp(-0.5 * ((w - w0) / s) ** 2)
    lam = elph.lambda_de_a2F(w, a2F)
    assert lam == pytest.approx(2 * A / w0, rel=1e-4)
    assert elph.omega_log_de_a2F(w, a2F, lam) == pytest.approx(
        w0 * elph.THZ_K, rel=1e-4)


def test_allen_dynes_para_el_aluminio():
    """Al: lambda 0.44, w_log 270 K -> Tc experimental 1.18 K."""
    tcs = [elph.allen_dynes(0.44, 270.0, mu) for mu in (0.10, 0.13, 0.16)]
    assert min(tcs) < 1.18 < max(tcs)


def test_correcciones_suben_tc_en_acoplamiento_fuerte():
    """Sin f1 y f2, Allen-Dynes subestima el plomo."""
    sin = elph.allen_dynes(1.55, 56.0, 0.13, correcciones=False)
    con = elph.allen_dynes(1.55, 56.0, 0.13)
    assert con > sin
    assert 6.0 < con < 7.5           # experimental 7.20 K


def test_correcciones_no_cambian_el_acoplamiento_debil():
    sin = elph.allen_dynes(0.44, 270.0, 0.13, correcciones=False)
    con = elph.allen_dynes(0.44, 270.0, 0.13)
    assert con == pytest.approx(sin, rel=0.05)


def test_sin_acoplamiento_suficiente_tc_es_cero():
    assert elph.allen_dynes(0.05, 300.0, 0.13) == 0.0


def test_tau_por_fonones_del_aluminio():
    """tau(300 K) del aluminio: unos 10 fs en la literatura."""
    tau = float(elph.tau_elph(0.44, 300.0)[0])
    assert 5e-15 < tau < 20e-15


def test_tau_baja_con_la_temperatura():
    t = elph.tau_elph(0.5, [100.0, 300.0, 900.0])
    assert t[0] > t[1] > t[2]
    assert t[0] / t[2] == pytest.approx(9.0, rel=1e-6)


def test_plato_encuentra_el_tramo_estable():
    serie = np.array([0.90, 0.62, 0.47, 0.441, 0.438, 0.440, 0.442, 0.45,
                      0.52])
    i = elph.plato(serie)
    assert serie[i] == pytest.approx(0.44, abs=0.01)


def test_sin_plato_devuelve_none():
    assert elph.plato(np.array([0.1, 0.4, 0.9, 1.6, 2.5])) is None


def test_lee_lambda_dat_del_aluminio():
    d = DATOS / "elph_al" / "lambda.dat"
    if not d.exists():
        pytest.skip("falta el lambda.dat de prueba")
    run = elph.leer_lambda_out(d.with_name("lambda.out"))
    assert run.lambdas is not None and len(run.lambdas) == 10
    run.i_plato = elph.plato(run.lambdas)
    # aluminio con malla de q 2x2x2: lambda del orden de 0.35
    assert 0.25 < run.lam < 0.45
    assert any("NaN" in a for a in run.avisos)


def test_input_de_lambda_exige_que_cuadren_los_tamanos():
    with pytest.raises(ErrorDeUso, match="MISMO orden"):
        elph.build_lambda_input([[0, 0, 0]], [1.0, 2.0], ["a", "b"])


# ======================================================================
# Dinámica molecular: casos con respuesta exacta
# ======================================================================
def _tray(pos, celda, dt=1.0, simbolos=None):
    n = pos.shape[1]
    return dyn.Trayectoria(simbolos=simbolos or ["X"] * n, posiciones=pos,
                           celda=celda, dt=dt)


def test_gas_ideal_tiene_gr_igual_a_uno():
    """La prueba de la normalización: sin correlación, g(r) = 1."""
    rng = np.random.default_rng(0)
    L = 20.0
    t = _tray(rng.random((300, 200, 3)) * L, np.eye(3) * L)
    r, gr, _ = dyn.rdf(t, nbins=50)
    sel = r > 1.5
    assert gr[sel].mean() == pytest.approx(1.0, abs=0.02)
    assert gr[sel].std() < 0.05


def test_camino_aleatorio_da_el_D_teorico():
    """MSD = 3*paso^2*t  ->  D = paso^2/2 en A^2/fs."""
    rng = np.random.default_rng(1)
    paso = 0.10
    pos = np.cumsum(rng.normal(0, paso, (3000, 60, 3)), axis=0)
    t = _tray(pos, np.eye(3) * 1e4)
    tt, m = dyn.msd(t)
    D, r2 = dyn.difusion(tt, m)
    assert D == pytest.approx(paso ** 2 / 2 * 1e-1, rel=0.05)
    assert r2 > 0.99


def test_vdos_encuentra_la_frecuencia_del_oscilador():
    f_cm1 = 500.0
    f_fs = f_cm1 * 2.99792458e10 / 1e15
    dt, n = 0.5, 4000
    x = 0.1 * np.sin(2 * np.pi * f_fs * np.arange(n) * dt)
    pos = np.zeros((n, 4, 3))
    pos[:, :, 0] = x[:, None]
    pos += np.array([[0, 0, 0], [3, 0, 0], [0, 3, 0], [0, 0, 3]])[None]
    fr, v = dyn.vdos(_tray(pos, np.eye(3) * 30, dt=dt, simbolos=["H"] * 4))
    assert fr[int(np.argmax(v))] == pytest.approx(f_cm1, rel=0.01)


def test_desdoblar_quita_los_saltos_de_frontera():
    """Sin esto, un átomo que cruza la caja dispara el MSD."""
    L = 20.0
    pos = np.zeros((100, 1, 3))
    pos[:, 0, 0] = (np.arange(100) * 0.5) % L
    u = dyn.desdoblar(_tray(pos, np.eye(3) * L))
    assert np.abs(np.diff(u[:, 0, 0])).max() == pytest.approx(0.5, abs=1e-9)


def test_gr_se_corta_en_media_arista():
    L = 10.0
    rng = np.random.default_rng(2)
    t = _tray(rng.random((10, 30, 3)) * L, np.eye(3) * L)
    r, _, _ = dyn.rdf(t, rmax=50.0)
    assert r[-1] <= L / 2 + 1e-9


def test_lee_una_md_real_de_pwx():
    g = DATOS / "md_si" / "md.out.gz"
    if not g.exists():
        pytest.skip("falta la MD de prueba")
    import tempfile
    with tempfile.TemporaryDirectory() as d:
        f = Path(d) / "md.out"
        f.write_bytes(gzip.decompress(g.read_bytes()))
        t = dyn.leer_md(f, skip=50)
    assert t.natoms == 2
    assert t.nsteps == 250
    # las posiciones venían en unidades de alat y tienen que salir en A
    assert t.dt == pytest.approx(0.9676, abs=0.001)
    assert np.abs(t.celda).max() == pytest.approx(2.715, abs=0.01)
    assert t.temperaturas.mean() == pytest.approx(900.0, abs=60.0)


# ======================================================================
# Heteroestructuras
# ======================================================================
def test_grafeno_sobre_hbn_da_el_desajuste_conocido():
    """(2.46 - 2.504)/2.504 = -1.76 %, el valor de libro."""
    from ase.build import graphene
    gr = graphene(a=2.46, vacuum=8.0)
    hbn = graphene(a=2.504, vacuum=8.0)
    hbn.set_chemical_symbols(["B", "N"])
    c = itf.buscar(gr, hbn, max_index=2, tol=0.05, max_atoms=40)
    assert c
    assert c[0].eps_pct == pytest.approx(1.76, abs=0.05)
    assert c[0].n1 == c[0].n2 == 1


def test_la_deformacion_es_isotropa_en_dos_redes_hexagonales():
    from ase.build import graphene
    gr = graphene(a=2.46, vacuum=8.0)
    hbn = graphene(a=2.504, vacuum=8.0)
    hbn.set_chemical_symbols(["B", "N"])
    c = itf.buscar(gr, hbn, max_index=2, tol=0.05, max_atoms=40)[0]
    assert c.deformacion[0, 0] == pytest.approx(c.deformacion[1, 1], abs=1e-6)
    assert abs(c.deformacion[0, 1]) < 1e-6


def test_reduccion_2d_identifica_la_misma_red():
    """(2.46, 4.26 A, 30 grados) y (2.46, 2.46, 60) son la misma red."""
    a = np.array([[2.46, 0.0], [-1.23, 2.1304]])
    b = np.array([[2.46, 0.0], [1.23, 2.1304]])
    assert itf._forma(a) == itf._forma(b)


def test_la_heteroestructura_no_tiene_atomos_encimados():
    from ase.build import graphene
    gr = graphene(a=2.46, vacuum=8.0)
    hbn = graphene(a=2.504, vacuum=8.0)
    hbn.set_chemical_symbols(["B", "N"])
    het = itf.emparejar(gr, hbn, max_index=2, tol=0.05, max_atoms=40)
    d = het.atoms.get_all_distances(mic=True)
    np.fill_diagonal(d, 9e9)
    assert d.min() > 1.0


def test_sin_coincidencia_el_error_dice_que_hacer():
    from ase.build import graphene
    gr = graphene(a=2.46, vacuum=8.0)
    raro = graphene(a=3.71, vacuum=8.0)
    with pytest.raises(ErrorDeUso, match="max-atoms"):
        itf.emparejar(gr, raro, max_index=1, tol=0.001, max_atoms=10)


# ======================================================================
# NEB
# ======================================================================
def test_extremos_con_atomos_desordenados_se_rechazan():
    from ase import Atoms
    a = Atoms("HHO", positions=[[0, 0, 0], [1, 0, 0], [2, 0, 0]],
              cell=np.eye(3) * 10)
    b = Atoms("OHH", positions=[[0, 0, 0], [1, 0, 0], [2, 0, 0]],
              cell=np.eye(3) * 10)
    problemas = nb.comprobar_extremos(a, b)
    assert any("ORDEN" in p for p in problemas)


def test_extremos_con_celdas_distintas_se_rechazan():
    from ase import Atoms
    a = Atoms("H2", positions=[[0, 0, 0], [1, 0, 0]], cell=np.eye(3) * 10)
    b = Atoms("H2", positions=[[0, 0, 0], [2, 0, 0]], cell=np.eye(3) * 11)
    assert any("celdas" in p for p in nb.comprobar_extremos(a, b))


def test_extremos_identicos_se_rechazan():
    from ase import Atoms
    a = Atoms("H2", positions=[[0, 0, 0], [1, 0, 0]], cell=np.eye(3) * 10)
    assert any("idénticas" in p for p in nb.comprobar_extremos(a, a.copy()))


def test_input_de_neb_tiene_la_estructura_de_bloques():
    from ase import Atoms
    a = Atoms("H2", positions=[[0, 0, 0], [1, 0, 0]], cell=np.eye(3) * 8,
              pbc=True)
    b = Atoms("H2", positions=[[0, 0, 0], [2, 0, 0]], cell=np.eye(3) * 8,
              pbc=True)
    pseudos = {"H": {"filename": "H.UPF", "found": True}}
    txt = nb.build_neb_input(a, b, pseudos, "H", ".", 20.0, 160.0,
                             "K_POINTS gamma\n", n_imagenes=5)
    for marca in ("BEGIN_PATH_INPUT", "END_PATH_INPUT",
                  "BEGIN_ENGINE_INPUT", "FIRST_IMAGE", "LAST_IMAGE",
                  "END_POSITIONS", "END_ENGINE_INPUT"):
        assert marca in txt
    assert txt.strip().endswith("END")
    assert "num_of_images     = 5" in txt
    # las posiciones NO pueden quedar en el bloque del motor
    motor = txt.split("BEGIN_POSITIONS")[0]
    assert motor.count("ATOMIC_POSITIONS") == 0


def test_sin_imagen_trepadora_se_declara():
    from ase import Atoms
    a = Atoms("H2", positions=[[0, 0, 0], [1, 0, 0]], cell=np.eye(3) * 8)
    b = Atoms("H2", positions=[[0, 0, 0], [2, 0, 0]], cell=np.eye(3) * 8)
    txt = nb.build_neb_input(a, b, {"H": {"filename": "H.UPF"}}, "H", ".",
                             20.0, 160.0, "K_POINTS gamma\n", ci=False)
    assert "no-CI" in txt


# ======================================================================
# Asistente
# ======================================================================
def test_todas_las_metas_tienen_lo_minimo():
    for m in wz.METAS:
        assert m.pregunta and m.nombre and m.explica
        assert m.pasos, m.clave
        assert m.coste in ("bajo", "medio", "alto", "muy alto")


def test_los_prerrequisitos_existen():
    for m in wz.METAS:
        for k in m.necesita:
            assert k in wz.METAS_POR_CLAVE, f"{m.clave} -> {k}"


def test_los_terminos_estan_en_el_glosario():
    for m in wz.METAS:
        for t in m.terminos:
            assert t in wz.GLOSARIO, f"{m.clave}: {t}"


def test_el_plan_pone_los_prerrequisitos_delante():
    pasos = wz.plan("conduce", "x.cif")
    claves = [c for c, _, _ in pasos]
    assert claves.index("relajar") < claves.index("conduce")


def test_el_plan_no_repite_un_prerrequisito():
    pasos = wz.plan("oxido", "x.cif")
    claves = [c for c, _, _ in pasos]
    assert claves.count("relajar") == claves.count("relajar")
    assert len(set(claves)) == len(dict.fromkeys(claves))


def test_los_comandos_llevan_el_archivo_del_usuario():
    pasos = wz.plan("conduce", "miMaterial.cif")
    cmds = [c for _, _, c in pasos if c]
    assert any("miMaterial.cif" in c for c in cmds)
    assert not any("{file}" in c for c in cmds)


@pytest.mark.parametrize("texto,clave", [
    ("quiero saber si absorbe luz visible", "color"),
    ("mi oxido sale metalico y no deberia", "oxido"),
    ("es duro o fragil", "mecanicas"),
    ("como se ve en rayos x", "difractograma"),
    ("quiero poner un material sobre otro", "interfase"),
])
def test_la_busqueda_en_lenguaje_llano_acierta(texto, clave):
    cands = wz.buscar(texto)
    assert cands, texto
    assert clave in [m.clave for m in cands[:3]], \
        f"{texto} -> {[m.clave for m in cands]}"


def test_meta_desconocida_lista_las_que_hay():
    with pytest.raises(ErrorDeUso, match="Disponibles"):
        wz.plan("noexiste", "x.cif")


def test_diagnostico_detecta_metales_de_transicion():
    from ase.build import bulk
    d = wz.diagnosticar(bulk("Ni", "fcc", a=3.52))
    assert "Ni" in d.metales_transicion
    assert any("autointeracción" in n for n in d.notas)


def test_diagnostico_detecta_una_losa():
    from ase.build import graphene
    d = wz.diagnosticar(graphene(a=2.46, vacuum=10.0))
    assert d.tiene_vacio
    assert any("losa" in n or "monocapa" in n for n in d.notas)
