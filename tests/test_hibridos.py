"""Funcionales híbridos: el input, la malla de EXX y el aviso de coste."""

import pytest
from ase.build import bulk

from qekit.core.errors import ErrorDeUso
from qekit.modules import inputgen


PSEUDOS = {"Si": {"filename": "Si.UPF", "found": True, "z_valence": 4.0}}


def _input(**kw):
    return inputgen.build_pw_input(
        atoms=bulk("Si", "diamond", 5.43), pseudos=PSEUDOS,
        calculation="scf", prefix="si", pseudo_dir=".", ecutwfc=30,
        ecutrho=120, kcard="K_POINTS automatic\n  4 4 4 0 0 0\n", **kw)


def test_hse_escribe_sus_parametros():
    txt = _input(hibrido="hse")
    assert "input_dft        = 'hse'" in txt
    assert "exx_fraction     = 0.25" in txt
    assert "screening_parameter = 0.106" in txt


def test_pbe0_no_lleva_apantallamiento():
    txt = _input(hibrido="pbe0")
    assert "input_dft        = 'pbe0'" in txt
    assert "screening_parameter" not in txt


def test_la_fraccion_se_puede_cambiar():
    assert "exx_fraction     = 0.4" in _input(hibrido="pbe0", exx_fraction=0.4)


def test_la_malla_de_exx_va_al_input():
    txt = _input(hibrido="hse", exx_grid=(2, 2, 2))
    for i in (1, 2, 3):
        assert f"nqx{i}             = 2" in txt


def test_sin_malla_de_exx_se_usa_gamma():
    txt = _input(hibrido="hse")
    assert "nqx1             = 1" in txt


def test_funcional_inventado():
    with pytest.raises(ErrorDeUso):
        _input(hibrido="magico")


# ----------------------------------------------------------------------
# a través de generate: divisibilidad y avisos
# ----------------------------------------------------------------------
def _generar(tmp_path, **kw):
    opts = inputgen.GenOptions(preset="scf", outdir=str(tmp_path),
                               pseudo_dir="/no/existe", kspacing=0.55, **kw)
    return inputgen.generate(bulk("Si", "diamond", 5.4073), opts)


def test_la_malla_de_exx_tiene_que_dividir_la_de_k(tmp_path):
    """pw.x se detiene con 'nqx must be a divisor of nk'."""
    with pytest.raises(ErrorDeUso) as exc:
        _generar(tmp_path, hibrido="hse", exx_grid=(3, 3, 3))
    assert "DIVIDIR" in str(exc.value)


def test_una_malla_valida_pasa(tmp_path):
    rep = _generar(tmp_path, hibrido="hse", exx_grid=(2, 2, 2))
    assert "HSE06" in rep


def test_gen_escribe_los_dos_guiones_portables(tmp_path):
    _generar(tmp_path)
    sh = (tmp_path / "run.sh").read_text()
    py = (tmp_path / "run.py").read_text()
    assert "pipefail" in sh
    assert "scf.in" in py and "pw.x" in py
    compile(py, str(tmp_path / "run.py"), "exec")


def test_magnetizacion_por_elemento_invalida_es_error_de_uso():
    with pytest.raises(ErrorDeUso) as exc:
        inputgen.parse_magnetization("Fe=mucho,O=0", ["Fe", "O"])
    assert "Fe" in str(exc.value) and "mucho" in str(exc.value)


def test_el_coste_se_da_con_un_numero(tmp_path):
    """factor ~ 3 + 2.6 n_q; medido: nq=8 dio x23.5, nq=64 dio x168."""
    r1 = _generar(tmp_path / "a", hibrido="hse", exx_grid=(1, 1, 1))
    r8 = _generar(tmp_path / "b", hibrido="hse", exx_grid=(2, 2, 2))
    r64 = _generar(tmp_path / "c", hibrido="hse", exx_grid=(4, 4, 4))
    assert "unas 6 veces" in r1
    assert "unas 24 veces" in r8
    assert "unas 169 veces" in r64


def test_avisa_de_que_la_malla_q_cambia_el_resultado(tmp_path):
    rep = _generar(tmp_path, hibrido="hse", exx_grid=(2, 2, 2))
    assert "CONVERGENCIA" in rep
    assert "2.68" in rep and "1.41" in rep, "los numeros medidos, no adjetivos"


def test_con_una_sola_q_avisa_de_que_no_sirve_para_citar(tmp_path):
    rep = _generar(tmp_path, hibrido="hse", exx_grid=(1, 1, 1))
    assert "sobrestimado" in rep


def test_avisa_de_que_no_hay_bandas_con_hibrido(tmp_path):
    opts = inputgen.GenOptions(preset="bands", outdir=str(tmp_path),
                               pseudo_dir="/no/existe", kspacing=0.55,
                               hibrido="hse")
    rep2 = inputgen.generate(bulk("Si", "diamond", 5.4073), opts)
    assert "NO puede hacer un cálculo 'bands' con EXX" in rep2


def test_el_divisor_sugerido_es_valido():
    for n in (4, 6, 8, 11, 12):
        d = inputgen._divisor(n)
        assert n % d == 0 and d <= n
