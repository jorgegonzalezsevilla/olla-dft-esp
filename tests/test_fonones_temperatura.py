"""Fonones a temperatura electrónica."""

import numpy as np
import pytest

from qekit.core.errors import ErrorDeUso, FaltanDatos
from qekit.modules import tphonons as tp


def test_degauss_es_kb_por_T():
    """k_B = 6.3336e-6 Ry/K; 300 K son 1.90 mRy."""
    assert tp.degauss_de_T(300.0) == pytest.approx(0.0019001, abs=1e-6)
    assert tp.degauss_de_T(1.0 / tp.KB_RY) == pytest.approx(1.0, rel=1e-9)


def test_la_conversion_va_en_los_dos_sentidos():
    for T in (100.0, 300.0, 6000.0):
        assert tp.T_de_degauss(tp.degauss_de_T(T)) == pytest.approx(T, rel=1e-9)


def test_temperatura_negativa_se_rechaza():
    with pytest.raises(ErrorDeUso):
        tp.degauss_de_T(-5.0)


def _barrido(datos):
    """datos: {T: [frecuencias en cm-1]}."""
    r = tp.BarridoT(temperaturas=sorted(datos), gamma_only=True)
    for T, f in datos.items():
        r.frecuencias[float(T)] = np.array(f, dtype=float)
    return r


def test_cuenta_los_modos_imaginarios():
    r = _barrido({300: [-120.0, -50.0, 0.0, 480.0]})
    assert len(r.imaginarias(300.0)) == 2


def test_el_ruido_numerico_no_cuenta_como_imaginario():
    """Un -3 cm-1 en un modo acustico es ruido, no una inestabilidad."""
    r = _barrido({300: [-3.0, -1.0, 0.0, 500.0]})
    assert len(r.imaginarias(300.0)) == 0


def test_encuentra_la_temperatura_de_estabilizacion():
    r = _barrido({300: [-100.0, 400.0], 600: [-50.0, 410.0],
                  900: [50.0, 420.0]})
    T = tp.temperatura_de_estabilizacion(r)
    assert T is not None
    assert 600 < T < 900, f"salio {T}"


def test_si_nunca_hay_imaginarios_no_hay_estabilizacion():
    r = _barrido({300: [100.0, 400.0], 900: [110.0, 410.0]})
    assert tp.temperatura_de_estabilizacion(r) is None
    assert "No hay modos imaginarios" in tp.report(r)


def test_si_siempre_hay_imaginarios_se_dice():
    r = _barrido({300: [-100.0, 400.0], 3000: [-90.0, 410.0]})
    assert tp.temperatura_de_estabilizacion(r) is None
    txt = tp.report(r)
    assert "no es de origen electrónico" in txt


def test_avisa_si_el_numero_de_imaginarios_no_es_monotono():
    r = _barrido({300: [-100.0, 400.0], 600: [50.0, 400.0],
                  900: [-80.0, 400.0]})
    assert not tp.monotono(r)
    assert "NO baja de forma monótona" in tp.report(r)
    assert "sería inventada" in tp.report(r)


def test_el_modo_blando_ignora_los_acusticos():
    """En Gamma hay tres modos a cero que no dicen nada del reblandecimiento."""
    r = _barrido({300: [0.0, 0.0, 0.0, -80.0, 500.0, 500.0]})
    assert r.modo_blando(300.0) == pytest.approx(-80.0)


def test_la_estabilizacion_cae_ENTRE_dos_puntos():
    """Con -50 a 600 K y +50 a 900 K, el cruce esta en 750, no en 900."""
    r = _barrido({300: [-100.0, 400.0], 600: [-50.0, 410.0],
                  900: [50.0, 420.0]})
    assert tp.temperatura_de_estabilizacion(r) == pytest.approx(750.0, abs=1.0)


def test_el_reporte_recuerda_que_los_iones_no_se_mueven():
    r = _barrido({300: [100.0, 400.0], 900: [110.0, 410.0]})
    txt = tp.report(r)
    assert "temperatura ELECTRÓNICA" in txt
    assert "iones siguen estando quietos" in txt


def test_sin_datos_se_queja():
    with pytest.raises(FaltanDatos):
        tp.report(tp.BarridoT(temperaturas=[300.0]))


def test_una_sola_temperatura_no_es_un_barrido(tmp_path):
    from ase.build import bulk
    with pytest.raises(ErrorDeUso):
        tp.prepare(bulk("Al", "fcc", 4.05), [300.0], outdir=str(tmp_path))


def test_una_temperatura_absurda_se_rechaza(tmp_path):
    from ase.build import bulk
    with pytest.raises(ErrorDeUso):
        tp.prepare(bulk("Al", "fcc", 4.05), [300.0, 99000.0],
                   outdir=str(tmp_path))


def test_prepare_impone_fermi_dirac(tmp_path):
    from ase.build import bulk
    run, rep = tp.prepare(bulk("Al", "fcc", 4.05), [300.0, 3000.0],
                          outdir=str(tmp_path), pseudo_dir="/no/existe")
    assert "fermi-dirac" in rep
    scf = (tmp_path / "T00300" / "scf.in").read_text()
    assert "'fermi-dirac'" in scf
    assert "degauss" in scf
    # y la de 3000 K tiene un degauss diez veces mayor
    scf2 = (tmp_path / "T03000" / "scf.in").read_text()
    d1 = float([l for l in scf.splitlines() if "degauss" in l][0].split("=")[1])
    d2 = float([l for l in scf2.splitlines() if "degauss" in l][0].split("=")[1])
    assert d2 / d1 == pytest.approx(10.0, rel=1e-3)
