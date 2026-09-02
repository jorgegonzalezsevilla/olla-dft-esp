"""Parámetros V intersitio: leerlos de hp.x y escribir la tarjeta HUBBARD."""

import pytest

from qekit.core.errors import ErrorDeUso
from qekit.modules import hubbard as hb


# Formato del <prefix>.Hubbard_parameters.dat de hp.x cuando el scf pidio
# DFT+U+V. Los indices del segundo atomo estan en la numeracion de la
# SUPERCELDA de vecinos que monta hp.x, no en la de la celda de entrada.
CON_V = """
  =-------------------------------------------------------------------=

                           Hubbard U parameters:

       site n.  type  label  spin  new_type  new_label  Hubbard U (eV)
         1        1    Ni      1      1         Ni         5.0104
         2        2    O       1      2         O          0.0000

  =-------------------------------------------------------------------=

       Hubbard V parameters:
       (adapted for a supercell 3x3x3)

         Atom 1     Atom 2     Distance (Bohr)   Hubbard V (eV)

           1 Ni       1 Ni        0.000000          5.0104
           1 Ni      19 O         3.947603          0.7521
           1 Ni      23 O         3.947603          0.7519
           1 Ni      55 O        11.842810          0.0031
           2 O        2 O         0.000000          0.0000

  =-------------------------------------------------------------------=
"""

SIN_V = """
  =-------------------------------------------------------------------=

                           Hubbard U parameters:

       site n.  type  label  spin  new_type  new_label  Hubbard U (eV)
         1        1    Ni      1      1         Ni         5.4429

  =-------------------------------------------------------------------=

          chi0 matrix :
   -0.384115    0.057138
"""


@pytest.fixture
def con_v(tmp_path):
    f = tmp_path / "NiO.Hubbard_parameters.dat"
    f.write_text(CON_V)
    return f


def test_lee_los_pares_v(con_v):
    pares, sup = hb.leer_v(con_v)
    assert len(pares) == 5
    assert sup == (3, 3, 3)


def test_distingue_el_sitio_del_vecino(con_v):
    pares, _ = hb.leer_v(con_v)
    sitio = [p for p in pares if p.es_sitio]
    vecino = [p for p in pares if not p.es_sitio]
    assert len(sitio) == 2 and len(vecino) == 3
    assert sitio[0].V == pytest.approx(5.0104), "V(i,i) es la U del sitio"


def test_los_indices_de_la_supercelda_se_conservan(con_v):
    """El 19 y el 23 son atomos de la supercelda de vecinos de hp.x."""
    pares, _ = hb.leer_v(con_v)
    js = sorted(p.j for p in pares if not p.es_sitio)
    assert js == [19, 23, 55]


def test_sin_seccion_v_devuelve_lista_vacia(tmp_path):
    f = tmp_path / "x.Hubbard_parameters.dat"
    f.write_text(SIN_V)
    pares, sup = hb.leer_v(f)
    assert pares == [] and sup is None


def test_archivo_inexistente(tmp_path):
    with pytest.raises(ErrorDeUso):
        hb.leer_v(tmp_path / "no_existe.dat")


def test_el_reporte_sin_v_explica_por_que(tmp_path):
    f = tmp_path / "x.Hubbard_parameters.dat"; f.write_text(SIN_V)
    pares, sup = hb.leer_v(f)
    txt = hb.report_v(pares, sup)
    assert "lda_plus_u_kind=2" in txt


def test_el_reporte_ordena_por_magnitud(con_v):
    pares, sup = hb.leer_v(con_v)
    txt = hb.report_v(pares, sup, umbral=0.01)
    assert "0.7521" in txt
    assert "0.0031" not in txt, "por debajo del umbral no se lista"
    assert "1 pares por debajo" in txt


# ----------------------------------------------------------------------
# tarjeta HUBBARD
# ----------------------------------------------------------------------
def test_la_tarjeta_lleva_U_y_V(con_v):
    pares, _ = hb.leer_v(con_v)
    sitios = hb.leer_parametros(con_v)
    card = hb.tarjeta_hubbard(sitios, pares)
    assert card.startswith("HUBBARD (ortho-atomic)")
    assert "U Ni-3d 5.0104" in card
    assert "V Ni-3d O-2p 1 19 0.7521" in card


def test_la_tarjeta_no_repite_la_U_de_cada_sitio(con_v):
    sitios = hb.leer_parametros(con_v)
    card = hb.tarjeta_hubbard(sitios, [])
    assert card.count("U Ni-") == 1


def test_la_tarjeta_descarta_los_V_pequenos(con_v):
    pares, _ = hb.leer_v(con_v)
    sitios = hb.leer_parametros(con_v)
    assert "0.0031" not in hb.tarjeta_hubbard(sitios, pares, umbral_v=0.01)
    assert "0.0031" in hb.tarjeta_hubbard(sitios, pares, umbral_v=0.001)


def test_la_tarjeta_no_mete_V_de_un_sitio_consigo_mismo(con_v):
    """V(i,i) ya esta puesto como U; repetirlo lo contaria dos veces."""
    pares, _ = hb.leer_v(con_v)
    sitios = hb.leer_parametros(con_v)
    card = hb.tarjeta_hubbard(sitios, pares)
    assert "V Ni-3d Ni-3d 1 1" not in card


def test_la_proyeccion_va_en_la_cabecera(con_v):
    sitios = hb.leer_parametros(con_v)
    card = hb.tarjeta_hubbard(sitios, [], proyeccion="atomic")
    assert card.startswith("HUBBARD (atomic)")
