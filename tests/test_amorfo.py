"""Generación de sólidos amorfos."""

import numpy as np
import pytest
from ase.data import atomic_numbers, covalent_radii

from qekit.core.errors import ErrorDeUso
from qekit.modules import amorphous as am


# ----------------------------------------------------------------------
# fórmula y empaquetado
# ----------------------------------------------------------------------
def test_la_formula_se_expande():
    assert am.formula_a_simbolos("SiO2", 8) == ["Si"] * 8 + ["O"] * 16
    assert am.formula_a_simbolos("GeTe", 4) == ["Ge"] * 4 + ["Te"] * 4
    assert len(am.formula_a_simbolos("Al2O3", 2)) == 10


def test_una_formula_absurda_se_queja():
    with pytest.raises(ErrorDeUso):
        am.formula_a_simbolos("Xx2", 4)
    with pytest.raises(ErrorDeUso):
        am.formula_a_simbolos("", 4)


def test_la_densidad_sale_exacta():
    at = am.empaquetar(am.formula_a_simbolos("SiO2", 8), densidad=2.2)
    assert am.densidad_de(at) == pytest.approx(2.2, rel=1e-6)


def test_la_densidad_fija_la_celda():
    """Doblar la densidad encoge la arista un factor 2^(1/3)."""
    s = am.formula_a_simbolos("SiO2", 8)
    a = am.celda_para_densidad(s, 2.2)
    b = am.celda_para_densidad(s, 4.4)
    assert a / b == pytest.approx(2 ** (1 / 3), rel=1e-9)


def test_ningun_par_queda_solapado():
    """La razon de ser del empaquetado: sin esto el primer paso de MD revienta."""
    simbolos = am.formula_a_simbolos("SiO2", 8)
    at = am.empaquetar(simbolos, densidad=2.2, semilla=3)
    d = at.get_all_distances(mic=True)
    for i in range(len(at)):
        for j in range(i + 1, len(at)):
            rmin = am.FACTOR_MINIMO * (
                covalent_radii[atomic_numbers[simbolos[i]]]
                + covalent_radii[atomic_numbers[simbolos[j]]])
            assert d[i, j] >= rmin - 1e-9, (i, j, d[i, j], rmin)


def test_semillas_distintas_dan_estructuras_distintas():
    s = am.formula_a_simbolos("SiO2", 8)
    a = am.empaquetar(s, 2.2, semilla=1)
    b = am.empaquetar(s, 2.2, semilla=2)
    assert not np.allclose(a.get_positions(), b.get_positions())


def test_la_misma_semilla_da_la_misma_estructura():
    s = am.formula_a_simbolos("SiO2", 8)
    a = am.empaquetar(s, 2.2, semilla=7)
    b = am.empaquetar(s, 2.2, semilla=7)
    assert np.allclose(a.get_positions(), b.get_positions())


def test_una_densidad_imposible_se_queja():
    """A 30 g/cm3 el SiO2 no cabe sin solaparse."""
    with pytest.raises(ErrorDeUso):
        am.empaquetar(am.formula_a_simbolos("SiO2", 8), densidad=30.0,
                      intentos=300)


def test_densidad_negativa():
    with pytest.raises(ErrorDeUso):
        am.empaquetar(["Si"], densidad=-1.0)


# ----------------------------------------------------------------------
# protocolo
# ----------------------------------------------------------------------
def test_la_velocidad_de_temple_sale_de_los_pasos():
    p = am.Protocolo(T_fundido=3000, T_final=300, pasos_temple=1000, dt_fs=1.0)
    # 2700 K en 1 ps = 2.7e15 K/s
    assert p.velocidad_temple == pytest.approx(2.7e15, rel=1e-6)


def test_mas_pasos_de_temple_es_mas_lento():
    a = am.Protocolo(pasos_temple=1000).velocidad_temple
    b = am.Protocolo(pasos_temple=10000).velocidad_temple
    assert b == pytest.approx(a / 10, rel=1e-9)


def test_el_tiempo_total_cuadra():
    p = am.Protocolo(pasos_fundido=400, pasos_temple=1200, pasos_recocido=200,
                     dt_fs=1.0)
    assert p.pasos == 1800
    assert p.ps_totales == pytest.approx(1.8)


# ----------------------------------------------------------------------
# análisis estructural
# ----------------------------------------------------------------------
def _cristobalita():
    """SiO2 cristalina: cada Si con 4 O y cada O con 2 Si, exacto."""
    from ase.spacegroup import crystal
    return crystal(("Si", "O"), basis=[(0, 0, 0), (0.125, 0.125, 0.125)],
                   spacegroup=227, cellpar=[7.16, 7.16, 7.16, 90, 90, 90])


def test_la_coordinacion_de_la_silice_es_4_y_2():
    """En cualquier SiO2 de red tetraedrica: Si-O = 4, O-Si = 2."""
    at = _cristobalita()
    c = am.coordinaciones(at)
    assert c[("Si", "O")] == pytest.approx(4.0, abs=0.01)
    assert c[("O", "Si")] == pytest.approx(2.0, abs=0.01)
    assert c[("O", "O")] == pytest.approx(0.0, abs=0.01)


def test_la_distancia_Si_O_de_la_silice():
    at = _cristobalita()
    d = am.distancia_media(at, "Si", "O")
    assert 1.5 < d < 1.8, f"salio {d:.3f} A; el enlace Si-O mide 1.61 A"


def test_el_corte_por_radios_no_es_un_corte_unico():
    """Con un corte global, el O-O de la silice contaria como enlace."""
    at = _cristobalita()
    por_radios = am.coordinaciones(at)
    global_3A = am.coordinaciones(at, corte=3.0)
    assert por_radios[("O", "O")] == pytest.approx(0.0, abs=0.01)
    assert global_3A[("O", "O")] > 3, "con corte unico, los O se 'enlazan'"
