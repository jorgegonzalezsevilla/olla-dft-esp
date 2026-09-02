"""La licencia pública tiene que seguir siendo la que se decidió.

Olla-DFT se publica bajo GPL-3.0. Cambiarla por descuido (copiando un
pyproject de otro proyecto, o aceptando el valor por defecto de una
plantilla) no produce ningún error visible, así que se comprueba.
"""

from pathlib import Path

import pytest

RAIZ = Path(__file__).resolve().parent.parent


def _leer(nombre: str) -> str:
    f = RAIZ / nombre
    if not f.exists():
        pytest.skip(f"{nombre} no está en el árbol instalado")
    return f.read_text(encoding="utf-8")


def test_existe_el_archivo_de_licencia_gpl3():
    texto = _leer("LICENSE")
    assert "GNU GENERAL PUBLIC LICENSE" in texto
    assert "Version 3, 29 June 2007" in texto
    assert "END OF TERMS AND CONDITIONS" in texto


def test_pyproject_declara_gpl3_y_al_autor():
    texto = _leer("pyproject.toml")
    assert "GNU General Public License v3" in texto
    assert 'license = { file = "LICENSE" }' in texto
    assert "Jorge Enrique González Sevilla" in texto
    assert "Private :: Do Not Upload" not in texto


def test_se_conserva_la_atribucion_de_los_datos_de_terceros():
    """Los datos de dispersión atómica vienen de pymatgen (MIT).

    La MIT exige conservar su aviso de copyright al redistribuir.
    """
    from qekit import data as _data
    aviso = Path(_data.__file__).parent / "LEEME_datos.txt"
    assert aviso.exists(), "falta la atribución de los datos incluidos"
    texto = aviso.read_text(encoding="utf-8")
    assert "pymatgen" in texto
    assert "MIT" in texto
    assert "Permission is hereby granted" in texto


def test_inventario_de_terceros_cubre_las_dependencias():
    texto = _leer("THIRD_PARTY_NOTICES.md")
    for dep in ("numpy", "scipy", "matplotlib", "spglib", "seekpath",
                "ASE", "pymatgen", "Quantum ESPRESSO", "phono3py", "mace"):
        assert dep in texto, f"{dep} no está en el inventario de terceros"
    assert "LGPL" in texto
