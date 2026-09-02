"""Errores de uso, pistas de fallo y el arreglo de ElasticRun.delta.

Todo esto salió del barrido de regresión de la 0.12.0: son casos que el
programa ya vivía, no hipótesis.
"""

import sqlite3

import numpy as np
import pytest

from qekit.core import runner, themes
from qekit.core.errors import ErrorDeUso, FaltanDatos
from qekit.modules import audit, elastic, feedback


# ----------------------------------------------------------------------
# ElasticRun.delta — el fallo real: export() lo pedía y no existía
# ----------------------------------------------------------------------
def test_elastic_delta_sale_de_las_deformaciones_usadas():
    run = elastic.ElasticRun(deltas=[0.0, 0.005, -0.005, 0.010, -0.010])
    assert run.delta == pytest.approx(0.010)


def test_elastic_delta_sin_deformaciones_es_cero():
    assert elastic.ElasticRun().delta == 0.0


def test_elastic_export_escribe_el_encabezado(tmp_path):
    """La regresión que se coló: export() reventaba con AttributeError."""
    C = np.diag([160.0, 160.0, 160.0, 76.6, 76.6, 76.6])
    C[0, 1] = C[1, 0] = C[0, 2] = C[2, 0] = C[1, 2] = C[2, 1] = 61.7
    run = elastic.ElasticRun(C=C, family="cúbico", natoms=2, volume=40.05,
                             deltas=[0.0, 0.01, -0.01],
                             components=[None, 0, 0],
                             stresses=[None, None, None])
    escritos = elastic.export(run, outdir=str(tmp_path))
    assert len(escritos) == 2
    texto = (tmp_path / "ELASTIC_C.dat").read_text()
    assert "delta" in texto
    leido = np.loadtxt(tmp_path / "ELASTIC_C.dat")
    assert leido.shape == (6, 6)
    assert leido[0, 0] == pytest.approx(160.0, abs=0.1)


# ----------------------------------------------------------------------
# Errores de uso
# ----------------------------------------------------------------------
def test_error_de_uso_sigue_siendo_valueerror():
    """Para no romper el código (ni las pruebas) que atrapaban ValueError."""
    assert issubclass(ErrorDeUso, ValueError)
    assert issubclass(FaltanDatos, ErrorDeUso)


def test_malla_mal_escrita_es_error_de_uso():
    from qekit.cli import _malla
    assert _malla("8x8x8") == (8, 8, 8)
    assert _malla("4,4,2") == (4, 4, 2)
    with pytest.raises(ErrorDeUso):
        _malla("1x2")
    with pytest.raises(ErrorDeUso):
        _malla("8x8xocho")
    with pytest.raises(ErrorDeUso):
        _malla("8x8x0")


def test_plantilla_desconocida_delata_la_confusion_con_revista():
    """'-t nature' es el error típico: nature es revista, no plantilla."""
    with pytest.raises(ErrorDeUso) as exc:
        themes.load("nature")
    assert "--journal nature" in str(exc.value)


def test_plantilla_de_verdad_desconocida_lista_las_disponibles():
    with pytest.raises(ErrorDeUso) as exc:
        themes.load("noexiste123")
    assert "Disponibles" in str(exc.value)
    assert "--journal" not in str(exc.value)


def test_incidencia_de_uso_no_guarda_traza(tmp_path):
    try:
        raise ErrorDeUso("--grid necesita tres numeros")
    except ErrorDeUso as exc:
        inc = feedback.registrar(exc=exc, tipo="uso", dir_=str(tmp_path))
    assert inc.tipo == "uso"
    assert inc.traceback == ""
    assert "ErrorDeUso" in inc.excepcion


def test_incidencia_de_error_si_guarda_traza(tmp_path):
    try:
        raise RuntimeError("reventó de verdad")
    except RuntimeError as exc:
        inc = feedback.registrar(exc=exc, dir_=str(tmp_path))
    assert inc.tipo == "error"
    assert "Traceback" in inc.traceback


def test_estadisticas_separan_uso_de_error(tmp_path):
    try:
        raise ErrorDeUso("mala bandera")
    except ErrorDeUso as e:
        feedback.registrar(exc=e, tipo="uso", dir_=str(tmp_path))
    try:
        raise RuntimeError("falla real")
    except RuntimeError as e:
        feedback.registrar(exc=e, dir_=str(tmp_path))
    st = feedback.estadisticas(dir_=str(tmp_path))
    assert st["total"] == 2
    assert st["errores"] == 1
    assert st["uso"] == 1
    # el de uso no contamina el conteo de "dónde falla el programa"
    assert sum(st["por_comando"].values()) == 1
    assert sum(st["uso_por_comando"].values()) == 1
    texto = feedback.report_estadisticas(st)
    assert "de uso: 1" in texto


# ----------------------------------------------------------------------
# Pistas de fallo de los binarios de QE
# ----------------------------------------------------------------------
@pytest.mark.parametrize("salida,clave", [
    ("mpirun has detected an attempt to run as root", "root"),
    ("bash: pw.x: command not found", "PATH"),
    ("Error in routine readpp (1):\n cannot open file Si.upf", "pseudopotencial"),
    ("Error in routine electrons (1):\n charge is wrong", "carga"),
    ("out of memory", "memoria"),
    ("*** buffer overflow detected ***: terminated", "compilación"),
])
def test_pistas_reconocen_las_causas_frecuentes(salida, clave):
    assert clave in runner.failure_hint(salida)


def test_sin_causa_reconocida_no_inventa_pista():
    assert runner.failure_hint("todo salió perfecto") == ""


def test_mensaje_de_fallo_prefiere_el_bloque_de_error_de_qe(tmp_path):
    log = tmp_path / "scf.out"
    log.write_text(
        "Program PWSCF starts\n"
        "     lots of banner\n" * 5 +
        "     Error in routine cdiaghg (1):\n"
        "     S matrix not positive definite\n"
        "     stopping ...\n"
        "     cola irrelevante\n")
    msg = runner.failure_message("scf", log)
    assert "Error in routine cdiaghg" in msg
    assert "Causa probable" in msg
    assert "ecutwfc" in msg
    assert "banner" not in msg


def test_mensaje_de_fallo_sin_bloque_de_error_muestra_la_cola(tmp_path):
    log = tmp_path / "pp.out"
    log.write_text("\n".join(f"linea {i}" for i in range(20)))
    msg = runner.failure_message("pp.x", log, lineas=3)
    assert "linea 19" in msg
    assert "linea 5" not in msg


def test_mensaje_de_fallo_con_archivo_inexistente_no_revienta(tmp_path):
    msg = runner.failure_message("scf", tmp_path / "no_existe.out")
    assert "no_existe.out" in msg


# ----------------------------------------------------------------------
# La base de datos: una columna mal escrita debe decir cuáles hay
# ----------------------------------------------------------------------
def _base_vacia(tmp_path):
    db = tmp_path / "q.db"
    con = sqlite3.connect(str(db))
    con.executescript(audit.ESQUEMA)
    con.commit(); con.close()
    return db


def test_columna_inexistente_lista_las_reales(tmp_path):
    db = _base_vacia(tmp_path)
    with pytest.raises(ErrorDeUso) as exc:
        audit.query("SELECT prefix FROM calculos", db_path=str(db))
    msg = str(exc.value)
    assert "no such column" in msg
    assert "formula" in msg and "ecutwfc" in msg


def test_columnas_devuelve_el_esquema(tmp_path):
    db = _base_vacia(tmp_path)
    cols = audit.columnas(str(db))
    assert "ruta" in cols and "energia_eV" in cols


def test_solo_select(tmp_path):
    db = _base_vacia(tmp_path)
    with pytest.raises(ErrorDeUso):
        audit.query("DELETE FROM calculos", db_path=str(db))


# ----------------------------------------------------------------------
# Las banderas de figura se validan ANTES de calcular
# ----------------------------------------------------------------------
def test_plantilla_mala_no_deja_salida_a_medias(tmp_path):
    """El fallo real: el difractograma se escribía y DESPUÉS reventaba."""
    from qekit.cli import main
    cif = tmp_path / "Si.cif"
    cif.write_text(_CIF_SI)
    salida = tmp_path / "fuera"
    rc = main(["xrd", str(cif), "-o", str(salida), "-t", "nature",
               "--format", "png"])
    assert rc == 2
    assert not salida.exists()


def test_plantilla_buena_si_corre(tmp_path):
    from qekit.cli import main
    cif = tmp_path / "Si.cif"
    cif.write_text(_CIF_SI)
    salida = tmp_path / "fuera"
    rc = main(["xrd", str(cif), "-o", str(salida), "-t", "journal",
               "--no-plot"])
    assert rc == 0
    assert (salida / "XRD.dat").exists()


def test_no_plot_no_valida_el_estilo(tmp_path):
    """Sin figura, la plantilla da igual: no debe estorbar."""
    from qekit.cli import main
    cif = tmp_path / "Si.cif"
    cif.write_text(_CIF_SI)
    rc = main(["xrd", str(cif), "-o", str(tmp_path / "f2"),
               "-t", "nature", "--no-plot"])
    assert rc == 0


_CIF_SI = """data_Si
_cell_length_a 5.43070
_cell_length_b 5.43070
_cell_length_c 5.43070
_cell_angle_alpha 90.0
_cell_angle_beta 90.0
_cell_angle_gamma 90.0
_symmetry_space_group_name_H-M 'P 1'
loop_
_atom_site_label
_atom_site_fract_x
_atom_site_fract_y
_atom_site_fract_z
Si1 0.00000 0.00000 0.00000
Si2 0.00000 0.50000 0.50000
Si3 0.50000 0.00000 0.50000
Si4 0.50000 0.50000 0.00000
Si5 0.25000 0.25000 0.25000
Si6 0.25000 0.75000 0.75000
Si7 0.75000 0.25000 0.75000
Si8 0.75000 0.75000 0.25000
"""
