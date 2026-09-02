"""Regresiones de mecánica de la CLI y de los módulos de propiedades.

Cada prueba fija un comportamiento que en su día falló: se documenta en la
docstring qué pasaba, para que si alguien vuelve a tocarlo entienda por qué
la prueba existe.
"""

from pathlib import Path

import numpy as np
import pytest
from ase.build import bulk
from ase.io import write

PSEUDO_DIR = "/usr/share/espresso/pseudo"


# ----------------------------------------------------------------------
# 1. phonons: --raman fuerza Γ aunque no se pase --gamma
# ----------------------------------------------------------------------
def test_una_corrida_en_gamma_no_intenta_dibujar_la_dispersion():
    """Con --raman y sin --gamma, gamma_only queda True pero la CLI decidía
    dibujar con la bandera --gamma: llamaba a plot() sin band_freqs y moría
    con AttributeError después de toda la cadena."""
    from qekit.modules import phonons
    run = phonons.PhononRun(gamma_only=True, raman=True)
    assert run.band_freqs is None
    assert phonons.has_dispersion(run) is False


def test_una_corrida_con_dispersion_si_se_dibuja():
    from qekit.modules import phonons
    run = phonons.PhononRun(gamma_only=False, qgrid=(2, 2, 2))
    assert phonons.has_dispersion(run) is False       # aún sin collect
    run.band_freqs = np.zeros((5, 6))
    run.qdist = np.linspace(0, 1, 5)
    assert phonons.has_dispersion(run) is True


# ----------------------------------------------------------------------
# 2. kappa: --metal
# ----------------------------------------------------------------------
def _kappa_falso(monkeypatch):
    """phono3py es opcional: se sustituye la parte que lo necesita por dos
    configuraciones fijas, que es todo lo que hace falta para ver qué
    ocupaciones llevan los inputs."""
    from qekit.modules import kappa
    si = bulk("Si", "diamond", 5.43)
    sc = si.repeat((2, 2, 2))
    monkeypatch.setattr(kappa, "preparar", lambda *a, **k: object())
    monkeypatch.setattr(kappa, "configuraciones",
                        lambda ph: ([sc.copy(), sc.copy()], []))
    return si


def _ocupaciones(pw_in: Path) -> str:
    for linea in pw_in.read_text().splitlines():
        if "occupations" in linea:
            return linea
    raise AssertionError(f"sin occupations en {pw_in}")


def test_kappa_sin_metal_escribe_ocupaciones_fijas(tmp_path, monkeypatch):
    from qekit.cli import main
    si = _kappa_falso(monkeypatch)
    write(tmp_path / "si.cif", si)
    out = tmp_path / "k"
    assert main(["kappa", str(tmp_path / "si.cif"), "-o", str(out),
                 "--pseudo-dir", PSEUDO_DIR]) == 0
    assert "fixed" in _ocupaciones(out / "fc3" / "d0000" / "pw.in")


def test_kappa_con_metal_escribe_smearing(tmp_path, monkeypatch):
    """Antes el aislante iba fijo a True: un metal recibía
    occupations='fixed' en los 57 scf y ninguno convergía."""
    from qekit.cli import main
    si = _kappa_falso(monkeypatch)
    write(tmp_path / "si.cif", si)
    out = tmp_path / "k"
    assert main(["kappa", str(tmp_path / "si.cif"), "--metal", "-o", str(out),
                 "--pseudo-dir", PSEUDO_DIR]) == 0
    txt = (out / "fc3" / "d0000" / "pw.in").read_text()
    assert "smearing" in _ocupaciones(out / "fc3" / "d0000" / "pw.in")
    assert "degauss" in txt


# ----------------------------------------------------------------------
# 3. qha: a(T) convencional
# ----------------------------------------------------------------------
def _modos(V, w0=300.0, gamma=1.5, V0=40.0, n=6):
    return np.full(n, w0 * (V0 / V) ** gamma)


def test_factor_convencional_de_al_y_si():
    """Antes el factor 4 se aplicaba solo con natoms == 2: el Al (fcc, un
    átomo en la primitiva) devolvía V^(1/3) y lo llamaba parámetro de red."""
    from qekit.modules import qha
    assert qha.factor_convencional(bulk("Al", "fcc", 4.05)) == pytest.approx(4.0)
    assert qha.factor_convencional(bulk("Si", "diamond", 5.43)) == pytest.approx(4.0)
    assert qha.factor_convencional(bulk("Fe", "bcc", 2.87)) == pytest.approx(2.0)
    assert qha.es_cubico(bulk("Al", "fcc", 4.05))
    assert not qha.es_cubico(bulk("Mg", "hcp", a=3.21, c=5.21))


def test_a_de_t_del_aluminio_es_el_parametro_convencional():
    from qekit.modules import qha
    al = bulk("Al", "fcc", 4.05)
    V0 = al.get_volume()                          # 4.05^3 / 4
    V = np.linspace(0.92, 1.08, 9) * V0
    E = 0.02 * (V - V0) ** 2
    F = [_modos(v, V0=V0) for v in V]
    r = qha.run(V, E, F, T=np.array([0.0, 10.0]), natoms=1, cubico=True,
                factor_conv=qha.factor_convencional(al))
    assert r.a_convencional
    assert r.a_T[0] == pytest.approx(4.05, abs=0.02)
    assert "convencional" in qha.report(r, T_ref=0.0)


def test_a_de_t_del_silicio_con_estructura():
    from qekit.modules import qha
    si = bulk("Si", "diamond", 5.43)
    V0 = si.get_volume()
    V = np.linspace(0.92, 1.08, 9) * V0
    F = [_modos(v, V0=V0) for v in V]
    r = qha.run(V, 0.02 * (V - V0) ** 2, F, T=np.array([0.0, 10.0]),
                natoms=2, cubico=True, factor_conv=qha.factor_convencional(si))
    assert r.a_T[0] == pytest.approx(5.43, abs=0.03)


def test_sin_estructura_a_de_t_se_etiqueta_como_primitiva(tmp_path):
    """Sin estructura no se puede saber el factor: se devuelve V_prim^(1/3),
    se dice en el informe, en la columna del .dat y en un aviso."""
    from qekit.modules import qha
    V = np.linspace(36.0, 44.0, 9)
    F = [_modos(v) for v in V]
    r = qha.run(V, 0.02 * (V - 40) ** 2, F, T=np.array([0.0, 10.0]),
                natoms=2, cubico=True)
    assert not r.a_convencional
    assert r.a_T[0] == pytest.approx(r.V_T[0] ** (1 / 3))
    assert any("PRIMITIVA" in a for a in r.avisos)
    assert "NO es el parámetro de red" in qha.report(r, T_ref=0.0)
    (f,) = qha.export(r, tmp_path)
    assert "Vprim^1/3" in Path(f).read_text()


def test_qha_cli_acepta_structure(tmp_path):
    from qekit.cli import main
    al = bulk("Al", "fcc", 4.05)
    write(tmp_path / "al.cif", al)
    V0 = al.get_volume()
    V = np.linspace(0.92, 1.08, 9) * V0
    filas = [[v, 0.02 * (v - V0) ** 2, *_modos(v, V0=V0)] for v in V]
    np.savetxt(tmp_path / "qha.dat", np.array(filas))
    assert main(["qha", str(tmp_path / "qha.dat"), "-o", str(tmp_path),
                 "--structure", str(tmp_path / "al.cif"), "--tmax", "20",
                 "--dt", "10", "--temp", "0", "--no-plot"]) == 0
    txt = (tmp_path / "QHA.dat").read_text()
    assert "a_conv(A)" in txt
    ultima = txt.strip().splitlines()[-1].split()
    assert float(ultima[-1]) == pytest.approx(4.05, abs=0.03)


# ----------------------------------------------------------------------
# 4 y 8. layers: monocapa que cruza la frontera periódica; radiación
# ----------------------------------------------------------------------
def _bicapa_mos2_partida():
    """Bicapa de MoS2 desplazada de modo que una de las capas queda con un S
    en z ≈ 0.6 Å y el Mo y el otro S en z ≈ 10–11 Å (cruza z = 0 / z = c)."""
    from ase.build import mx2
    m = mx2("MoS2", a=3.16, thickness=3.19, vacuum=None)
    bl = m.copy()
    bl.set_cell([m.cell[0], m.cell[1], [0, 0, 12.3]], scale_atoms=False)
    top = m.copy()
    top.set_cell(bl.cell)
    top.translate([0, 0, 6.15])
    bl += top
    bl.pbc = True
    bl.translate([0, 0, -1.0])
    bl.wrap()
    return bl


def test_la_monocapa_que_cruza_la_frontera_se_desenrolla_y_centra():
    """Antes make_slab usaba las posiciones envueltas: el grosor salía ~c
    entero y el centrado partía la capa en dos, y exfoliate heredaba una
    monocapa rota."""
    from qekit.core import layers
    bl = _bicapa_mos2_partida()
    res = layers.analyze(bl)
    assert len(res.layers) == 2
    zs = bl.positions[:, 2]
    partida = [k for k, L in enumerate(res.layers)
               if np.ptp(zs[L.indices]) > 6.0]
    assert partida, "la estructura de prueba debe tener una capa partida"
    for k in range(2):
        slab = layers.make_slab(bl, res, layer_index=k, vacuum=20.0)
        z = slab.positions[:, 2]
        assert z.max() - z.min() == pytest.approx(3.19, abs=1e-6)
        assert (z.max() + z.min()) / 2 == pytest.approx(slab.cell[2, 2] / 2,
                                                        abs=1e-6)
        assert slab.cell[2, 2] == pytest.approx(3.19 + 20.0, abs=1e-6)


def test_el_rotulo_de_la_reflexion_basal_lleva_la_radiacion_pedida():
    """Decía 'Cu Kα' fuera cual fuera --wavelength."""
    from qekit.core import layers
    from qekit.modules import xrd
    bl = _bicapa_mos2_partida()
    res = layers.analyze(bl)
    lam = xrd.wavelength_value("MoKa")
    txt = layers.report(bl, res, wavelength=lam,
                        radiation=xrd.wavelength_name("MoKa"))
    assert "Mo Kα" in txt and "Cu Kα" not in txt
    assert f"{lam:.4f}" in txt
    txt = layers.report(bl, res, wavelength=0.9, radiation=xrd.wavelength_name(0.9))
    assert "Cu Kα" not in txt and "λ dada" in txt


def test_nombre_de_la_radiacion():
    from qekit.modules import xrd
    assert xrd.wavelength_name("CuKa") == "Cu Kα"
    assert xrd.wavelength_name("CuKa1") == "Cu Kα1"
    assert xrd.wavelength_name("AgKa") == "Ag Kα"
    assert xrd.wavelength_name(1.54) == "λ dada"


# ----------------------------------------------------------------------
# 5. ballistic: bdl/bds son longitudes de celda, no alturas atómicas
# ----------------------------------------------------------------------
def _electrodo_y_dispersor():
    from ase import Atoms
    # electrodo: dos átomos en z = 0 y z = 2.0 en una celda de 4.0 Å
    e = Atoms("Al2", positions=[[0, 0, 0.0], [0, 0, 2.0]],
              cell=[[4, 0, 0], [0, 4, 0], [0, 0, 4.0]], pbc=True)
    # dispersor: tres átomos, celda de 10 Å en z
    d = Atoms("Al3", positions=[[0, 0, 0.0], [0, 0, 2.0], [0, 0, 4.0]],
              cell=[[4, 0, 0], [0, 4, 0], [0, 0, 10.0]], pbc=True)
    return e, d


def test_longitud_z_es_la_celda_en_unidades_de_alat():
    from qekit.modules import ballistic as bl
    e, d = _electrodo_y_dispersor()
    assert bl.longitud_z(e) == pytest.approx(1.0)      # 4.0 / alat 4.0
    assert bl.longitud_z(d) == pytest.approx(2.5)      # 10.0 / 4.0
    # y NO es la altura del último átomo, que era lo que se escribía antes
    assert bl.limites_z(e)[1] == pytest.approx(0.5)


def test_bdl_y_bds_salen_de_la_longitud_de_celda(tmp_path):
    """pwcond.x sitúa cada región de z = 0 a z = bd*: con la altura del
    último átomo (0.5 en vez de 1.0) el empalme quedaba a media celda."""
    from qekit.modules import ballistic as bl
    e, d = _electrodo_y_dispersor()
    bl.prepare(e, outdir=str(tmp_path), dispersor=d, pseudo_dir=PSEUDO_DIR)
    txt = (tmp_path / "cond.in").read_text()
    assert "bdl=1.000000" in txt
    assert "bds=2.500000" in txt
    assert "prefixr" not in txt and "bdr" not in txt


def test_ikind_2_se_rechaza_con_explicacion(tmp_path):
    from qekit.core.errors import ErrorDeUso
    from qekit.modules import ballistic as bl
    e, d = _electrodo_y_dispersor()
    with pytest.raises(ErrorDeUso, match="ikind=2"):
        bl.prepare(e, outdir=str(tmp_path), dispersor=d, ikind=2,
                   pseudo_dir=PSEUDO_DIR)


def test_ikind_1_sin_dispersor_se_rechaza(tmp_path):
    from qekit.core.errors import ErrorDeUso
    from qekit.modules import ballistic as bl
    e, _ = _electrodo_y_dispersor()
    with pytest.raises(ErrorDeUso, match="scatterer"):
        bl.prepare(e, outdir=str(tmp_path), ikind=1, pseudo_dir=PSEUDO_DIR)


def test_el_parser_de_ballistic_ya_no_ofrece_ikind_2():
    from qekit.cli import build_parser
    p = build_parser()
    with pytest.raises(SystemExit):
        p.parse_args(["ballistic", "e.cif", "--ikind", "2"])
    args = p.parse_args(["ballistic", "e.cif", "--ikind", "1"])
    assert args.ikind == 1


# ----------------------------------------------------------------------
# 6. transport: --nspin/--mag y --kspacing llegan a los inputs
# ----------------------------------------------------------------------
def test_transport_prepare_escribe_nspin_2_con_magnetizacion(tmp_path):
    from qekit.modules import transport as tr
    fe = bulk("Fe", "bcc", 2.87)
    tr.prepare(fe, outdir=str(tmp_path), pseudo_dir=PSEUDO_DIR, ecutwfc=30,
               ecutrho=240, grid=(4, 4, 4), insulator=False,
               magnetization={"Fe": 0.7})
    for nombre in ("scf.in", "nscf.in"):
        txt = (tmp_path / nombre).read_text()
        assert "nspin            = 2" in txt
        assert "starting_magnetization(1) = 0.7" in txt


def test_transport_cli_acepta_nspin_mag_y_kspacing(tmp_path):
    """El mensaje de --spin-resolved mandaba usar '--nspin 2 --mag', pero el
    parser de transport no tenía esas banderas; y --kspacing se aceptaba y
    se ignoraba."""
    from qekit.cli import main
    fe = bulk("Fe", "bcc", 2.87)
    write(tmp_path / "fe.cif", fe)
    out = tmp_path / "tr"
    assert main(["transport", str(tmp_path / "fe.cif"), "-o", str(out),
                 "--pseudo-dir", PSEUDO_DIR, "--metal", "--mag", "Fe=0.7",
                 "--grid", "4x4x4", "--kspacing", "0.6"]) == 0
    scf = (out / "scf.in").read_text()
    assert "nspin" in scf and "starting_magnetization(1)" in scf
    grueso = [l for l in scf.splitlines() if l.strip()
              and l.strip()[0].isdigit()][-1]
    out2 = tmp_path / "tr2"
    assert main(["transport", str(tmp_path / "fe.cif"), "-o", str(out2),
                 "--pseudo-dir", PSEUDO_DIR, "--metal", "--grid", "4x4x4",
                 "--kspacing", "0.15"]) == 0
    scf2 = (out2 / "scf.in").read_text()
    fino = [l for l in scf2.splitlines() if l.strip()
            and l.strip()[0].isdigit()][-1]
    assert "nspin" not in scf2
    assert grueso != fino          # el kspacing cambia la malla del scf


# ----------------------------------------------------------------------
# 7. elph: la columna Tc(K) de la tabla por ensanchamiento
# ----------------------------------------------------------------------
DATOS = Path(__file__).parent / "datos"

_LAMBDA_OUT_CON_TC = """\
     lambda = 1.500000 (   1.480000 )  <log w>=   56.000 K  N(Ef)=  2.0 at degauss= 0.010
     lambda = 1.550000 (   1.530000 )  <log w>=   55.000 K  N(Ef)=  2.1 at degauss= 0.020
lambda        omega_log          T_c
   1.50000      56.00000           7.12345
   1.55000      55.00000           7.54321
"""

_LAMBDA_OUT_SIN_TC = """\
     lambda = 1.500000 (   1.480000 )  <log w>=   56.000 K  N(Ef)=  2.0 at degauss= 0.010
     lambda = 1.550000 (   1.530000 )  <log w>=   55.000 K  N(Ef)=  2.1 at degauss= 0.020
"""


def test_la_muestra_del_aluminio_trae_la_tabla_de_tc():
    """run.Tc no se asignaba nunca: la columna Tc(K) salía NaN aunque
    lambda.out tuviera la tabla final de lambda.x."""
    from qekit.modules import elph
    f = DATOS / "elph_al" / "lambda.out"
    if not f.exists():
        pytest.skip("falta el lambda.out de prueba")
    run = elph.leer_lambda_out(f)
    assert run.Tc is not None and len(run.Tc) == len(run.lambdas) == 10
    assert "lambda.x" in run.Tc_fuente
    # en esta muestra lambda.x dejó omega_log y T_c en NaN, y así se refleja
    assert not np.any(np.isfinite(run.Tc))
    assert "Tc(K) de la tabla: lambda.x" in elph.report(run)


def test_tc_se_lee_de_la_tabla_de_lambda_x(tmp_path):
    from qekit.modules import elph
    (tmp_path / "lambda.out").write_text(_LAMBDA_OUT_CON_TC)
    (tmp_path / "lambda.in").write_text("20 0.12 1\n1\n 0 0 0 1.0\nelph.1\n0.13\n")
    run = elph.leer_lambda_out(tmp_path / "lambda.out")
    assert run.mustar == pytest.approx(0.13)          # el mu* de lambda.in
    assert run.Tc == pytest.approx([7.12345, 7.54321])
    assert run.omega_log == pytest.approx([56.0, 55.0])
    assert "lambda.x" in run.Tc_fuente and "0.13" in run.Tc_fuente


def test_sin_tabla_tc_se_calcula_con_allen_dynes_y_se_dice(tmp_path):
    from qekit.modules import elph
    (tmp_path / "lambda.out").write_text(_LAMBDA_OUT_SIN_TC)
    run = elph.leer_lambda_out(tmp_path / "lambda.out")
    esperado = [elph.allen_dynes(1.5, 56.0, 0.10, correcciones=False),
                elph.allen_dynes(1.55, 55.0, 0.10, correcciones=False)]
    assert run.Tc == pytest.approx(esperado)
    assert run.Tc[0] > 0
    assert "calculada por Olla-DFT" in run.Tc_fuente
    txt = elph.report(run)
    assert "calculada por Olla-DFT" in txt
    assert f"{esperado[0]:8.3f}" in txt


# ----------------------------------------------------------------------
# 9. derived: la temperatura de Slack y el bloque cúbico
# ----------------------------------------------------------------------
def _C_cubica(c11=165.6, c12=63.9, c44=79.5):
    C = np.zeros((6, 6))
    for i in range(3):
        C[i, i] = c11
        for j in range(3):
            if i != j:
                C[i, j] = c12
    for i in range(3, 6):
        C[i, i] = c44
    return C


def _C_hexagonal():
    # MoS2 aproximado: claramente no cúbico
    C = np.zeros((6, 6))
    C[0, 0] = C[1, 1] = 238.0
    C[2, 2] = 52.0
    C[0, 1] = C[1, 0] = 55.0
    C[0, 2] = C[2, 0] = C[1, 2] = C[2, 1] = 23.0
    C[3, 3] = C[4, 4] = 19.0
    C[5, 5] = (238.0 - 55.0) / 2
    return C


def test_kappa_de_slack_se_etiqueta_con_la_temperatura_pedida(tmp_path):
    """Decía '300 K' y 'kappa_Slack_300K' aunque se pasara --temp 600."""
    from qekit.modules import derived, elastic
    si = bulk("Si", "diamond", 5.43)
    m = elastic.moduli(_C_cubica())
    r = derived.analyze(m.B_hill, m.G_hill, si.get_masses(), si.get_volume(),
                        natoms=2, T=600.0)
    r300 = derived.analyze(m.B_hill, m.G_hill, si.get_masses(),
                           si.get_volume(), natoms=2, T=300.0)
    assert r.kappa_slack == pytest.approx(r300.kappa_slack / 2)
    txt = derived.report(r)
    assert "Slack, 600 K" in txt and "300 K" not in txt
    (f,) = derived.export(r, tmp_path)
    dat = Path(f).read_text()
    assert "kappa_Slack_600K" in dat and "300K" not in dat


def test_is_cubic_tensor_distingue_cubico_de_hexagonal():
    from qekit.modules import derived
    assert derived.is_cubic_tensor(_C_cubica())
    ruidoso = _C_cubica() + np.random.default_rng(0).normal(0, 0.5, (6, 6))
    assert derived.is_cubic_tensor(ruidoso)
    assert not derived.is_cubic_tensor(_C_hexagonal())


def test_el_bloque_cubico_no_se_imprime_para_un_hexagonal(tmp_path, capsys):
    from qekit.cli import main
    from ase.build import mx2
    m = mx2("MoS2", a=3.16, thickness=3.19, vacuum=6.0)
    m.pbc = True
    write(tmp_path / "mos2.cif", m)
    np.savetxt(tmp_path / "C.dat", _C_hexagonal())
    assert main(["derived", str(tmp_path / "mos2.cif"), "--cij",
                 str(tmp_path / "C.dat"), "-o", str(tmp_path)]) == 0
    assert "En un cristal cúbico" not in capsys.readouterr().out


def test_el_bloque_cubico_si_se_imprime_para_el_silicio(tmp_path, capsys):
    from qekit.cli import main
    write(tmp_path / "si.cif", bulk("Si", "diamond", 5.43))
    np.savetxt(tmp_path / "C.dat", _C_cubica())
    assert main(["derived", str(tmp_path / "si.cif"), "--cij",
                 str(tmp_path / "C.dat"), "-o", str(tmp_path),
                 "--temp", "500"]) == 0
    out = capsys.readouterr().out
    assert "En un cristal cúbico" in out
    assert "Slack, 500 K" in out


# ----------------------------------------------------------------------
# 10. eos: a0 se fija en el ajuste y el informe usa ese mismo valor
# ----------------------------------------------------------------------
def test_eos_fit_fija_a0_y_el_informe_lo_usa():
    """EOSFit.a0 quedaba siempre en None y a₀ se recalculaba en report();
    había además una expresión muerta en fit()."""
    from qekit.modules import eos
    V0, B0, Bp, E0 = 40.05, 94.2 / 160.21766208, 4.2, -20.0
    V = np.linspace(0.88 * V0, 1.12 * V0, 11)
    eta = (V0 / V) ** (2.0 / 3.0)
    E = E0 + 9 * V0 * B0 / 16 * (
        (eta - 1) ** 3 * Bp + (eta - 1) ** 2 * (6 - 4 * eta))
    run = eos.EOSRun(volumes=list(V), energies=list(E), natoms=2)
    run.cubic, run.conv_ratio = True, 4.0            # diamante: 4 primitivas
    f = eos.fit(run, equation="birch-murnaghan")
    assert f.a0 == pytest.approx((V0 * 4.0) ** (1 / 3), rel=1e-4)
    eos.fit_all(run)
    txt = eos.report(run)
    assert f"a₀ = {run.fits['birch-murnaghan'].a0:.5f} Å" in txt
    # sin celda cúbica no hay a0
    run2 = eos.EOSRun(volumes=list(V), energies=list(E), natoms=2)
    assert eos.fit(run2).a0 is None
