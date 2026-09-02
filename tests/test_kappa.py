"""Conductividad térmica de red.

Las pruebas que tocan phono3py se saltan si no está instalado: es una
dependencia opcional, y el resto del módulo (informe, acumulada, ajustes)
tiene que funcionar y probarse sin ella.
"""

import numpy as np
import pytest

from qekit.core.errors import ErrorDeUso, FaltanDatos
from qekit.modules import kappa as K

p3 = pytest.importorskip("phono3py", reason="phono3py es opcional")


def _si():
    from ase.build import bulk
    return bulk("Si", "diamond", 5.43)


# ----------------------------------------------------------------------
# Configuraciones desplazadas
# ----------------------------------------------------------------------
def test_el_silicio_2x2x2_pide_57_configuraciones():
    """Es un número fijo: lo fija la simetría, no una elección.

    Si cambia, o cambió phono3py o cambió la detección de simetría, y en
    cualquiera de los dos casos hay que enterarse.
    """
    ph = K.preparar(_si(), (2, 2, 2))
    s3, s2 = K.configuraciones(ph)
    assert len(s3) == 57
    assert len(s3[0]) == 16
    assert s2 == []


def test_una_supercelda_aparte_para_la_parte_armonica():
    ph = K.preparar(_si(), (2, 2, 2), dim_fc2=(3, 3, 3))
    s3, s2 = K.configuraciones(ph)
    assert len(s3) == 57 and len(s3[0]) == 16
    assert len(s2) >= 1 and len(s2[0]) == 54


def test_las_configuraciones_son_estructuras_de_ase_periodicas():
    ph = K.preparar(_si(), (2, 2, 2))
    s3, _ = K.configuraciones(ph)
    a = s3[0]
    assert all(a.pbc)
    assert a.get_chemical_symbols().count("Si") == 16
    # y están DESPLAZADAS respecto de la supercelda perfecta
    b = s3[1]
    assert not np.allclose(a.get_positions(), b.get_positions())


def test_el_desplazamiento_es_el_pedido():
    """Se mide sobre las propias estructuras, no sobre la API de phono3py.

    Los nombres internos de phono3py cambian entre versiones; la distancia
    que separa la supercelda desplazada de la perfecta, no.
    """
    ph = K.preparar(_si(), (2, 2, 2), distancia=0.05)
    perfecta = K._ase(ph.supercell)
    d = (K.configuraciones(ph)[0][0].get_positions()
         - perfecta.get_positions())
    movidos = np.linalg.norm(d, axis=1)
    assert movidos.max() == pytest.approx(0.05, abs=1e-9)
    assert (movidos > 1e-9).sum() == 1        # solo uno se mueve en la fc3


def test_escribir_inputs_hace_una_carpeta_por_configuracion(tmp_path):
    from qekit.modules import sweep
    ph = K.preparar(_si(), (2, 2, 2))
    s3, _ = K.configuraciones(ph)
    s3 = s3[:4]
    common = sweep.prepare_common(s3[0], "/usr/share/espresso/pseudo",
                                  25, 100, True)
    carp = K.escribir_inputs(s3, tmp_path, common, kspacing=0.4)
    assert len(carp) == 4
    for d in carp:
        assert (d / "pw.in").exists()
        txt = (d / "pw.in").read_text()
        assert "calculation" in txt and "K_POINTS" in txt
        assert txt.count("Si ") >= 16 or txt.count("Si") >= 16


def test_leer_fuerzas_dice_cuales_faltan(tmp_path):
    d1, d2 = tmp_path / "d0000", tmp_path / "d0001"
    d1.mkdir(); d2.mkdir()
    with pytest.raises(FaltanDatos) as e:
        K.leer_fuerzas([d1, d2], 16)
    assert "d0000" in str(e.value) and "d0001" in str(e.value)


# ----------------------------------------------------------------------
# Post-proceso, sin necesidad de phono3py
# ----------------------------------------------------------------------
def _run_sintetico(n=1.0, nT=8):
    T = np.linspace(100, 800, nT)
    k = 100.0 * (T / 300.0) ** (-n)
    kap = np.zeros((nT, 6))
    kap[:, 0] = kap[:, 1] = kap[:, 2] = k
    return K.KappaRun(formula="Si2", dim=(3, 3, 3), temperaturas=T,
                      kappa=kap, i300=int(np.argmin(abs(T - 300))),
                      fuente="Quantum ESPRESSO", n_config=57, n_atomos=16)


def test_la_media_del_tensor_es_la_traza_entre_tres():
    run = _run_sintetico()
    run.kappa[:, 1] *= 2.0
    esperado = (run.kappa[:, 0] + run.kappa[:, 1] + run.kappa[:, 2]) / 3
    assert np.allclose(run.kappa_media, esperado)


@pytest.mark.parametrize("n", [0.5, 1.0, 1.5, 2.0])
def test_el_exponente_de_la_temperatura_se_recupera(n):
    run = _run_sintetico(n)
    assert K.exponente_temperatura(run) == pytest.approx(n, abs=1e-6)


def test_sin_datos_no_hay_exponente():
    run = K.KappaRun()
    assert K.exponente_temperatura(run) is None


def test_el_informe_dice_que_es_rta_y_que_faltan_los_isotopos():
    r = K.report(_run_sintetico())
    assert "RTA" in r
    assert "isótopos" in r
    assert "cuatro" in r          # los procesos de cuatro fonones


def test_el_informe_avisa_de_que_el_potencial_no_es_dft():
    run = _run_sintetico()
    run.fuente = "MACE"
    r = K.report(run)
    assert any("no de DFT" in a for a in run.avisos)
    assert "MACE" in r


def test_el_informe_avisa_de_una_supercelda_pequena():
    run = _run_sintetico()
    run.dim = (2, 2, 2)
    K.report(run)
    assert any("pequeña" in a for a in run.avisos)


def test_una_supercelda_grande_no_dispara_el_aviso():
    run = _run_sintetico()
    run.dim = (3, 3, 3)
    K.report(run)
    assert not any("pequeña" in a for a in run.avisos)


def test_el_informe_reconoce_el_uno_partido_por_t_de_umklapp():
    assert "Umklapp" in K.report(_run_sintetico(1.0))


def test_export_escribe_kappa_y_el_informe(tmp_path):
    f = K.export(_run_sintetico(), str(tmp_path))
    assert any(x.endswith("KAPPA.dat") for x in f)
    assert any(x.endswith("KAPPA.txt") for x in f)
    dat = np.loadtxt(tmp_path / "KAPPA.dat")
    assert dat.shape == (8, 8)


# ----------------------------------------------------------------------
# Acumulada frente al recorrido libre medio
# ----------------------------------------------------------------------
def _run_con_modos(nq=40, nb=6, semilla=0):
    rng = np.random.default_rng(semilla)
    run = _run_sintetico()
    nT = len(run.temperaturas)
    run.pesos = np.ones(nq)
    run.velocidades = rng.random((nq, nb, 3)) * 50.0        # Å/ps
    run.gamma = rng.random((nT, nq, nb)) * 0.5 + 0.01       # THz
    run.cv = rng.random((nT, nq, nb)) + 0.1
    return run


def test_la_acumulada_va_de_cero_a_uno_y_no_baja():
    L, a = K.acumulada(_run_con_modos())
    assert a[-1] == pytest.approx(1.0)
    assert np.all(np.diff(a) >= -1e-12)
    assert np.all(np.diff(L) >= -1e-12)


def test_el_recorrido_al_noventa_por_ciento_es_mayor_que_al_cincuenta():
    run = _run_con_modos()
    assert (K.recorrido_representativo(run, 0.9)
            > K.recorrido_representativo(run, 0.5))


def test_sin_gammas_no_hay_acumulada():
    L, a = K.acumulada(_run_sintetico())
    assert L is None and a is None


def test_export_incluye_la_acumulada_cuando_la_hay(tmp_path):
    f = K.export(_run_con_modos(), str(tmp_path))
    assert any("recorrido" in x for x in f)
    d = np.loadtxt(tmp_path / "KAPPA_recorrido.dat")
    assert d.shape[1] == 2 and 0.0 <= d[:, 1].min() and d[:, 1].max() <= 1.0


# ----------------------------------------------------------------------
# CLI
# ----------------------------------------------------------------------
@pytest.mark.parametrize("txt,esp", [
    ("300", [300.0]),
    ("300,500", [300.0, 500.0]),
    ("100:300:3", [100.0, 200.0, 300.0]),
    ("300;500", [300.0, 500.0]),
])
def test_lista_de_temperaturas(txt, esp):
    from qekit.cli import _temperaturas
    assert _temperaturas(txt) == pytest.approx(esp)


@pytest.mark.parametrize("mal", ["100:800", "a:b:c", "100:800:0", "abc"])
def test_lista_de_temperaturas_mal_escrita(mal):
    from qekit.cli import _temperaturas
    with pytest.raises(ErrorDeUso):
        _temperaturas(mal)
