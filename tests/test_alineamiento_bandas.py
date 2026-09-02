"""Alineamiento de bandas: quitar el cero arbitrario y clasificar el tipo."""

import pytest

from qekit.core.errors import FaltanDatos
from qekit.modules import align


def _lado(nombre, vbm, gap, ref):
    """Un lado con VBM y gap dados, en una escala con el cero en `ref`."""
    return align.Lado(nombre=nombre, vbm=ref + vbm, cbm=ref + vbm + gap,
                      gap=gap, referencia=ref, ref_tipo="nivel de vacío")


def test_el_cero_arbitrario_se_va():
    """Dos calculos del MISMO material con distinto cero dan offset cero.

    Es la razon de ser del modulo: en hBN monocapa, cambiar el vacio de 16 a
    22 A movio el VBM crudo 0.60 eV sin que cambiara nada fisico.
    """
    a = _lado("A", -5.81, 4.63, ref=0.0)
    b = _lado("B", -5.81, 4.63, ref=137.4)      # otro cero cualquiera
    al = align.alinear(a, b)
    assert al.delta_v == pytest.approx(0.0, abs=1e-9)
    assert al.delta_c == pytest.approx(0.0, abs=1e-9)


def test_el_offset_es_la_diferencia_de_potenciales_de_ionizacion():
    a = _lado("A", -5.0, 2.0, ref=10.0)
    b = _lado("B", -6.0, 3.0, ref=-40.0)
    al = align.alinear(a, b)
    assert al.delta_v == pytest.approx(1.0)     # A esta 1 eV mas arriba
    assert al.delta_c == pytest.approx(0.0)     # -5+2 = -3 ; -6+3 = -3


def test_tipo_I_anidado():
    """El gap de A cae dentro del de B: los dos portadores van a A."""
    a = _lado("A", -5.0, 1.0, 0.0)
    b = _lado("B", -6.0, 3.0, 0.0)
    assert align.alinear(a, b).tipo == "I"


def test_tipo_II_escalonado():
    a = _lado("A", -5.0, 2.0, 0.0)      # VBM -5, CBM -3
    b = _lado("B", -6.0, 2.0, 0.0)      # VBM -6, CBM -4
    al = align.alinear(a, b)
    assert al.tipo == "II"
    txt = align.report(al)
    assert "electrón se va a B" in txt, "el CBM mas bajo es el de B"
    assert "hueco a A" in txt, "el VBM mas alto es el de A"


def test_tipo_III_roto():
    """El CBM de A queda por DEBAJO del VBM de B."""
    a = _lado("A", -9.0, 1.0, 0.0)      # VBM -9, CBM -8
    b = _lado("B", -6.0, 2.0, 0.0)      # VBM -6, CBM -4
    assert align.alinear(a, b).tipo == "III"


def test_sin_offset_apreciable_no_se_inventa_un_tipo():
    """20 meV es ruido de DFT, no un alineamiento escalonado."""
    a = _lado("A", -5.000, 2.000, 0.0)
    b = _lado("B", -5.010, 2.005, 0.0)
    al = align.alinear(a, b)
    assert al.tipo == "="
    assert "sería leer ruido" in align.report(al)


def test_el_cbm_se_situa_bien_en_la_escala_del_otro():
    """El CBM de A esta en gap_B + dEc, no en dEc.

    Con la confusion, dos materiales identicos salian clasificados como
    tipo III, porque dEc ~ 0 quedaba 'por debajo' del VBM de B.
    """
    a = _lado("A", -5.0, 4.63, 0.0)
    b = _lado("B", -5.0, 4.65, 0.0)
    al = align.alinear(a, b)
    assert al.tipo != "III"


def test_un_metal_no_tiene_cbm():
    a = align.Lado(nombre="metal", vbm=-5.0, cbm=None, referencia=0.0,
                   es_metal=True)
    b = _lado("B", -6.0, 2.0, 0.0)
    al = align.alinear(a, b)
    assert al.delta_v == pytest.approx(1.0)
    assert al.delta_c is None
    assert any("no tiene CBM" in x for x in al.avisos)


def test_avisa_si_la_meseta_de_vacio_no_es_plana():
    a = _lado("A", -5.0, 2.0, 0.0); a.planitud = 0.4
    b = _lado("B", -6.0, 2.0, 0.0); b.planitud = 0.01
    al = align.alinear(a, b)
    assert any("meseta de vacío de A" in x for x in al.avisos)
    assert not any("meseta de vacío de B" in x for x in al.avisos)


def test_el_modo_vacio_avisa_de_lo_que_ignora():
    al = align.alinear(_lado("A", -5.0, 2.0, 0.0), _lado("B", -6.0, 2.0, 0.0))
    assert any("dipolo de interfaz" in x for x in al.avisos)


def test_el_puente_de_la_interfaz_desplaza_el_offset():
    a = _lado("A", -5.0, 2.0, 0.0)
    b = _lado("B", -6.0, 2.0, 0.0)
    sin = align.alinear(a, b, modo="interfaz")
    con = align.alinear(a, b, modo="interfaz", puente=0.35)
    assert con.delta_v - sin.delta_v == pytest.approx(0.35)
    assert con.delta_c - sin.delta_c == pytest.approx(0.35)


def test_sin_referencia_no_hay_alineamiento():
    a = align.Lado(nombre="A", vbm=-5.0, referencia=None)
    b = _lado("B", -6.0, 2.0, 0.0)
    with pytest.raises(FaltanDatos):
        align.alinear(a, b)


def test_export_escribe_los_dos_archivos(tmp_path):
    al = align.alinear(_lado("A", -5.0, 2.0, 0.0), _lado("B", -6.0, 3.0, 0.0))
    f = align.export(al, str(tmp_path))
    assert len(f) == 2
    assert (tmp_path / "ALINEAMIENTO.txt").exists()
    assert "dEv_eV" in (tmp_path / "ALINEAMIENTO.dat").read_text()


def test_la_figura_pone_el_cbm_de_a_donde_lo_pone_el_reporte():
    """La figura dibujaba el CBM de A en dEc a secas, no en gap_B + dEc.

    Con A = B (mismo material, dEc = 0) el CBM de A salia en 0 eV: pegado a
    su propio VBM, cuando el reporte decia gap de 2 eV. La posicion que usa
    la figura tiene que ser exactamente la que usa `alinear` para clasificar.
    """
    a = _lado("A", -5.0, 2.0, 0.0)
    b = _lado("B", -5.5, 3.0, 0.0)
    al = align.alinear(a, b)
    pos = align.posiciones_en_escala_de_b(al)
    assert pos["v_b"] == 0.0
    assert pos["c_b"] == pytest.approx(b.gap)
    assert pos["v_a"] == pytest.approx(al.delta_v)
    assert pos["c_a"] == pytest.approx(b.gap + al.delta_c)
    # y, por tanto, el gap dibujado de A es el gap real de A
    assert pos["c_a"] - pos["v_a"] == pytest.approx(a.gap)


def test_la_figura_se_genera_con_el_cbm_correcto(tmp_path):
    pytest.importorskip("matplotlib")
    al = align.alinear(_lado("A", -5.0, 2.0, 0.0), _lado("B", -5.0, 2.0, 0.0))
    escritos = align.plot(al, outfile=str(tmp_path / "ali"), formats="png")
    assert escritos
    assert (tmp_path / "ali.png").exists()
