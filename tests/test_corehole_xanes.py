"""Pseudopotenciales con hueco de core, XPS y XANES.

Lo que se comprueba aquí sin necesitar Quantum ESPRESSO: las
configuraciones electrónicas, el formato de los inputs de ld1.x, la
extracción de la función de onda de core de un UPF real, y la lectura y
el análisis de un espectro XANES ya calculado.
"""

import gzip
from pathlib import Path

import numpy as np
import pytest

from qekit.core import atomconf
from qekit.core.compat import trapezoid
from qekit.core.errors import ErrorDeUso
from qekit.modules import corehole, xanes

DATOS = Path(__file__).parent / "datos"


@pytest.fixture(scope="module")
def upf_hueco(tmp_path_factory):
    """El UPF con hueco de core del silicio, generado con ld1.x."""
    origen = DATOS / "Si.hueco1s.UPF.gz"
    if not origen.exists():
        pytest.skip("falta el UPF de prueba")
    destino = tmp_path_factory.mktemp("ps") / "Si.hueco1s.UPF"
    destino.write_bytes(gzip.decompress(origen.read_bytes()))
    return destino


# ----------------------------------------------------------------------
# Configuraciones atómicas
# ----------------------------------------------------------------------
def test_aufbau_cuenta_todos_los_electrones():
    for z in range(1, 87):
        conf = atomconf.aufbau(z)
        assert sum(o for _, _, o in conf) == pytest.approx(z), f"Z={z}"


def test_configuracion_del_silicio():
    conf = atomconf.configuracion("Si")
    assert conf == [(1, 0, 2.0), (2, 0, 2.0), (2, 1, 6.0),
                    (3, 0, 2.0), (3, 1, 2.0)]


def test_excepciones_de_aufbau():
    """El cobre es 3d10 4s1, no 3d9 4s2."""
    conf = dict(((n, l), o) for n, l, o in atomconf.configuracion("Cu"))
    assert conf[(3, 2)] == 10.0
    assert conf[(4, 0)] == 1.0
    assert sum(conf.values()) == pytest.approx(29)


@pytest.mark.parametrize("el,core,val", [
    ("Si", ["1S", "2S", "2P"], ["3S", "3P"]),
    ("O", ["1S"], ["2S", "2P"]),
    ("Fe", ["1S", "2S", "2P", "3S", "3P"], ["3D", "4S"]),
])
def test_particion_core_valencia(el, core, val):
    c, v = atomconf.particion(el)
    assert [atomconf.etiqueta(n, l) for n, l, _ in c] == core
    assert [atomconf.etiqueta(n, l) for n, l, _ in v] == val


def test_siempre_hay_canal_s_p_y_d():
    """Sin el canal d vacío salen estados fantasma; sin el p, el sodio falla."""
    for el in ("Si", "O", "Na", "Fe", "Ti", "H"):
        ls = {l for _, _, l, _ in atomconf.canales_pseudo(el)}
        assert ls == {0, 1, 2}, el


def test_el_canal_d_vacio_del_oxigeno_es_3d_no_2d():
    """No existe un orbital 2d."""
    canales = atomconf.canales_pseudo("O")
    d = [c for c in canales if c[2] == 2][0]
    assert d[0] == "3D"


def test_hueco_quita_exactamente_un_electron():
    base = atomconf.configuracion("Si")
    hueco, nivel = atomconf.config_hueco("Si", "K")
    assert nivel == "1S"
    assert sum(o for _, _, o in base) - sum(o for _, _, o in hueco) == \
        pytest.approx(1.0)


def test_hueco_en_nivel_de_valencia_se_rechaza():
    """Un hueco en 3s del silicio no es un hueco de CORE."""
    with pytest.raises(ValueError, match="core"):
        atomconf.config_hueco("Si", "M1")


def test_borde_inexistente_para_el_elemento():
    with pytest.raises(ValueError):
        atomconf.config_hueco("H", "L23")


# ----------------------------------------------------------------------
# Inputs de ld1.x
# ----------------------------------------------------------------------
def test_input_ld1_tiene_la_forma_que_ld1_espera():
    conf = atomconf.configuracion("Si") + [(3, 2, -1.0)]
    canales = atomconf.canales_pseudo("Si")
    txt = corehole.input_ld1("Si", sorted(conf), canales, "Si.UPF", 1.6,
                             dft="PZ")
    lineas = txt.strip().split("\n")
    assert "iswitch=3" in txt
    assert "lgipaw_reconstruction=.true." in txt
    assert "zed=14.0" in txt
    # el número de estados va en su propia línea antes de la lista
    i = lineas.index("6")
    assert lineas[i + 1].startswith("1S")
    # el canal vacío lleva energía distinta de cero
    d = [l for l in lineas if l.startswith("3D")][-1]
    assert float(d.split()[4]) > 0


def test_dos_proyectores_comparten_el_numero_principal():
    """Numerarlos por línea hace que ld1.x falle con 'Zero norm'."""
    canales = atomconf.canales_pseudo("Si", proyectores=2)
    conf = sorted(atomconf.configuracion("Si") +
                  [(4, 0, -1.0), (4, 1, -1.0), (3, 2, -1.0)])
    txt = corehole.input_ld1("Si", conf, canales, "Si.UPF", 1.6,
                             pseudotype=3)
    tarjeta = [l for l in txt.strip().split("\n")
               if l[:2] in ("3S", "4S", "3P", "4P", "3D")]
    ns = [int(l.split()[1]) for l in tarjeta]
    # 3S y 4S son el mismo canal l=0 -> mismo índice
    assert ns[0] == ns[1]


# ----------------------------------------------------------------------
# Lectura del UPF y función de onda de core
# ----------------------------------------------------------------------
def test_lee_upf_con_cabecera_autocerrada(upf_hueco):
    """ld1.x de QE 6.6 escribe <PP_HEADER .../>, sin etiqueta de cierre."""
    p = corehole.leer_upf(upf_hueco)
    assert p.elemento == "Si"
    assert p.z_valence == pytest.approx(5.0)
    assert p.mesh == 1141
    assert "1S" in p.orbitales_core


def test_core_wfc_escribe_bloques_separados(upf_hueco, tmp_path):
    destino = tmp_path / "Si.wfc"
    corehole.core_wfc(upf_hueco, destino, orbital="1S")
    texto = destino.read_text()
    bloques = [b for b in texto.split("\n\n") if b.strip()]
    assert len(bloques) >= 3
    assert texto.startswith("#")
    for b in bloques:
        filas = [l for l in b.split("\n") if l.strip() and not l.startswith("#")]
        assert len(filas) == 1141


def test_las_funciones_de_core_estan_normalizadas(upf_hueco, tmp_path):
    destino = tmp_path / "Si.wfc"
    corehole.core_wfc(upf_hueco, destino)
    bloques = [b for b in destino.read_text().split("\n\n") if b.strip()]
    for b in bloques:
        filas = [l for l in b.split("\n")
                 if l.strip() and not l.startswith("#")]
        d = np.array([[float(x) for x in l.split()] for l in filas])
        norma = trapezoid(d[:, 1] ** 2, d[:, 0])
        assert norma == pytest.approx(1.0, abs=0.01)


def test_los_orbitales_de_core_son_distintos(upf_hueco, tmp_path):
    """1S, 2S y 2P tienen que ser funciones diferentes, con nodos distintos."""
    destino = tmp_path / "Si.wfc"
    corehole.core_wfc(upf_hueco, destino)
    bloques = [b for b in destino.read_text().split("\n\n") if b.strip()]
    picos = []
    for b in bloques:
        filas = [l for l in b.split("\n")
                 if l.strip() and not l.startswith("#")]
        d = np.array([[float(x) for x in l.split()] for l in filas])
        picos.append(d[np.argmax(np.abs(d[:, 1])), 0])
    # el 1s está más cerca del núcleo que el 2s, y ese que el 3s
    assert picos[0] < picos[1] < picos[-1]


def test_core_wfc_sin_gipaw_avisa(tmp_path):
    falso = tmp_path / "malo.UPF"
    falso.write_text('<UPF version="2.0.1">\n<PP_HEADER mesh_size="10"/>\n'
                     "<PP_R>1 2 3</PP_R>\n</UPF>\n")
    with pytest.raises(ErrorDeUso, match="lgipaw_reconstruction"):
        corehole.core_wfc(falso, tmp_path / "x.wfc")


def test_core_wfc_orbital_inexistente(upf_hueco, tmp_path):
    with pytest.raises(ErrorDeUso, match="4F"):
        corehole.core_wfc(upf_hueco, tmp_path / "x.wfc", orbital="4F")


# ----------------------------------------------------------------------
# Etiquetas de especie
# ----------------------------------------------------------------------
def test_etiqueta_excitada_cabe_en_tres_caracteres():
    """QE usa CHARACTER(LEN=3): 'Si_h' se trunca en silencio y luego falla."""
    for el in ("H", "C", "Si", "Fe", "Zn"):
        assert len(xanes.etiqueta_excitada(el)) <= 3
        assert xanes.etiqueta_excitada(el).startswith(el)


# ----------------------------------------------------------------------
# XANES
# ----------------------------------------------------------------------
def test_input_de_xspectra_pone_xanes_file_en_plot():
    """En &input_xspectra hace fallar la lectura del namelist entero."""
    txt = xanes.build_xspectra_input("Si", 1, (1, 0, 0), "Si.wfc")
    antes, despues = txt.split("&plot", 1)
    assert "xanes_file" not in antes
    assert "xanes_file" in despues
    assert "calculation='xanes_dipole'" in antes


def test_input_de_xspectra_lleva_la_malla_al_final():
    txt = xanes.build_xspectra_input("Si", 1, (0, 0, 1), "Si.wfc",
                                     kgrid=(6, 6, 6))
    assert txt.strip().split("\n")[-1].split() == ["6", "6", "6", "1", "1", "1"]


@pytest.fixture
def espectro_si(tmp_path):
    src = DATOS / "xanes_si"
    if not src.is_dir():
        pytest.skip("faltan los espectros de prueba")
    for f in src.glob("xanes_*.dat"):
        (tmp_path / f.name).write_bytes(f.read_bytes())
    return xanes.collect(tmp_path, elemento="Si", borde="K")


def test_silicio_cubico_no_tiene_anisotropia(espectro_si):
    """Prueba fuerte: en un cristal cúbico las tres polarizaciones DEBEN
    dar el mismo espectro. Si el manejo de xepsilon estuviera mal, no.

    El umbral no es cero sino 1e-3: xspectra.x corta la fracción continua
    de Lanczos cuando el error baja de `xerror` (1e-3 aquí), así que las
    tres direcciones coinciden hasta ese ruido y no más. Medido: 3.7e-4.
    """
    assert len(espectro_si.componentes) == 3
    assert xanes._anisotropia(espectro_si) < 1e-3


def test_el_borde_cae_cerca_del_nivel_de_fermi(espectro_si):
    """El silicio es semiconductor: el borde arranca justo sobre E_F."""
    assert 0.0 < xanes.onset(espectro_si) < 3.0


def test_estructura_del_espectro_del_silicio(espectro_si):
    """Las tres primeras estructuras del borde K del silicio cristalino."""
    picos = sorted(e for e, _ in xanes._picos(espectro_si))
    assert 2.5 < picos[0] < 5.0          # máximo principal
    assert any(8.0 < p < 14.0 for p in picos)
    assert any(15.0 < p < 20.0 for p in picos)


def test_el_promedio_es_el_promedio(espectro_si):
    comps = np.array(list(espectro_si.componentes.values()))
    assert np.allclose(espectro_si.sigma, comps.mean(axis=0))


def test_espectro_positivo_y_creciente_en_el_borde(espectro_si):
    s = espectro_si.sigma
    e = espectro_si.energias
    assert (s >= 0).all()
    antes = s[e < -5].mean()
    despues = s[(e > 2) & (e < 6)].mean()
    assert despues > 20 * antes


def test_collect_sin_datos_avisa(tmp_path):
    with pytest.raises(ErrorDeUso, match="xanes_"):
        xanes.collect(tmp_path)


def test_export_escribe_las_componentes(espectro_si, tmp_path):
    salidas = xanes.export(espectro_si, str(tmp_path))
    d = np.loadtxt(tmp_path / "XANES.dat")
    assert d.shape[1] == 5          # E, promedio, x, y, z
    assert len(salidas) == 2
