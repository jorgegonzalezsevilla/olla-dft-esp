"""Fatbands: leer las proyecciones de projwfc.x y sumarlas bien."""

import numpy as np
import pytest

from qekit.core.errors import ErrorDeUso
from qekit.modules import bands


# Salida de projwfc.x reducida: silicio en Gamma y en un punto cualquiera.
# Los numeros son los de una corrida real (Si diamante, LDA): en Gamma la
# banda 1 es s pura y las bandas 2-4, degeneradas, son p puras.
PROJ = """     Atomic states used for projection
     state #   1: atom   1 (Si ), wfc  1 (l=0 m= 1)
     state #   2: atom   1 (Si ), wfc  2 (l=1 m= 1)
     state #   3: atom   1 (Si ), wfc  2 (l=1 m= 2)
     state #   4: atom   1 (Si ), wfc  2 (l=1 m= 3)
     state #   5: atom   2 (Si ), wfc  1 (l=0 m= 1)
     state #   6: atom   2 (Si ), wfc  2 (l=1 m= 1)
     state #   7: atom   2 (Si ), wfc  2 (l=1 m= 2)
     state #   8: atom   2 (Si ), wfc  2 (l=1 m= 3)

 k =   0.0000000000  0.0000000000  0.0000000000
==== e(   1) =    -5.88600 eV ====
     psi = 0.498*[#   1]+0.498*[#   5]
    |psi|^2 = 0.996
==== e(   2) =     6.04700 eV ====
     psi = 0.160*[#   2]+0.160*[#   3]+0.160*[#   4]+0.160*[#   6]+0.160*[#   7]
          +0.160*[#   8]
    |psi|^2 = 0.960
==== e(   3) =     6.04700 eV ====
     psi = 0.160*[#   2]+0.160*[#   3]+0.160*[#   4]+0.160*[#   6]+0.160*[#   7]
          +0.160*[#   8]
    |psi|^2 = 0.960

 k =   0.1000000000  0.0000000000  0.0000000000
==== e(   1) =    -5.70000 eV ====
     psi = 0.400*[#   1]+0.400*[#   5]+0.050*[#   2]+0.050*[#   6]
    |psi|^2 = 0.900
==== e(   2) =     5.90000 eV ====
     psi = 0.300*[#   3]+0.300*[#   7]
    |psi|^2 = 0.600
==== e(   3) =     5.95000 eV ====
     psi = 0.250*[#   4]+0.250*[#   8]
    |psi|^2 = 0.500
"""


@pytest.fixture
def proy(tmp_path):
    f = tmp_path / "projwfc.out"
    f.write_text(PROJ)
    return bands.leer_proyecciones(f)


def test_lee_los_estados_atomicos(proy):
    assert len(proy.estados) == 8
    assert proy.estados[0] == {"n": 1, "atomo": 1, "elemento": "Si",
                               "l": 0, "orb": "s"}
    assert proy.etiquetas == ["Si-s", "Si-p"]


def test_la_forma_del_arreglo(proy):
    assert proy.pesos.shape == (1, 2, 3, 8)     # spin, k, banda, estado


def test_en_gamma_la_banda_1_es_s_pura(proy):
    """Gamma_1 del silicio: enlace s, sin nada de p."""
    s = bands.peso_de(proy, "Si-s")[0]
    p = bands.peso_de(proy, "Si-p")[0]
    assert s[0, 0] == pytest.approx(0.996, abs=1e-3)
    assert p[0, 0] == pytest.approx(0.0, abs=1e-6)


def test_en_gamma_las_bandas_2_y_3_son_p_puras(proy):
    """Gamma_25': triplemente degenerada y de caracter p."""
    s = bands.peso_de(proy, "Si-s")[0]
    p = bands.peso_de(proy, "Si-p")[0]
    for ib in (1, 2):
        assert p[0, ib] == pytest.approx(0.960, abs=1e-3)
        assert s[0, ib] == pytest.approx(0.0, abs=1e-6)


def test_los_orbitales_suman_el_total(proy):
    s = bands.peso_de(proy, "Si-s")
    p = bands.peso_de(proy, "Si-p")
    todo = bands.peso_de(proy, "Si")
    assert np.allclose(s + p, todo, atol=1e-6)


def test_seleccionar_por_atomo(proy):
    """Los dos silicios son equivalentes: cada uno lleva la mitad."""
    a1 = bands.peso_de(proy, "atomo:1")
    a2 = bands.peso_de(proy, "atomo:2")
    assert np.allclose(a1[0, 0], a2[0, 0], atol=1e-6)


def test_seleccionar_solo_por_orbital(proy):
    assert np.allclose(bands.peso_de(proy, "s"), bands.peso_de(proy, "Si-s"))


def test_selector_inexistente(proy):
    with pytest.raises(ErrorDeUso) as exc:
        bands.peso_de(proy, "Ni-d")
    assert "Si-s" in str(exc.value), "el error tiene que decir qué SÍ hay"


def test_atomo_inexistente(proy):
    with pytest.raises(ErrorDeUso):
        bands.peso_de(proy, "atomo:9")


def test_un_archivo_que_no_es_de_projwfc(tmp_path):
    f = tmp_path / "projwfc.out"
    f.write_text("     JOB DONE.\n")
    with pytest.raises(ErrorDeUso):
        bands.leer_proyecciones(f)


def test_se_niega_a_mezclar_dos_calculos(proy):
    """Proyecciones del nscf sobre bandas del camino: el fallo silencioso."""
    class FalsoRes:
        eigenvalues = np.zeros((1, 40, 3))      # 40 puntos k, no 2
    class FalsoBS:
        result = FalsoRes()
        energies = FalsoRes.eigenvalues
    with pytest.raises(ErrorDeUso) as exc:
        bands.comprobar_compatibilidad(FalsoBS(), proy)
    assert "40" in str(exc.value) and "2" in str(exc.value)


def test_el_reporte_avisa_del_peso_que_falta(proy):
    """En la banda 3 del segundo k solo se proyecta el 50 %: hay que decirlo."""
    txt = bands.report_fat(proy, "Si")
    assert "NO cae dentro de ninguna" in txt


def test_lee_una_carpeta(tmp_path):
    (tmp_path / "projwfc.out").write_text(PROJ)
    p = bands.leer_proyecciones(tmp_path)
    assert len(p.estados) == 8


def test_carpeta_sin_projwfc(tmp_path):
    with pytest.raises(FileNotFoundError):
        bands.leer_proyecciones(tmp_path)
