"""Configuración de pytest.

Dos niveles de prueba, deliberadamente separados:

- las PURAS no necesitan Quantum ESPRESSO ni datos de cálculo: corren en
  segundos y son las que hay que poder correr en cada cambio;
- las marcadas @pytest.mark.qe necesitan binarios de QE o salidas ya
  calculadas, y se saltan solas si no están.

    pytest                     -> todo lo que se pueda correr aquí
    pytest -m "not qe"         -> solo las puras (segundos)
    OLLA_DFT_TEST_DATA=/ruta pytest -m qe   -> las que leen salidas de QE
"""

import os
import shutil
import tempfile
from pathlib import Path

import pytest

# Las pruebas nunca deben leer ni escribir la configuración real del usuario
# (idioma, pseudo_dir...). Se fija ANTES de importar qekit.config, que calcula
# la carpeta al importarse.
if "OLLA_DFT_CONFIG_DIR" not in os.environ and "QEKIT_CONFIG_DIR" not in os.environ:
    _CFG_TMP = tempfile.mkdtemp(prefix="olla-dft-tests-config-")
    os.environ["OLLA_DFT_CONFIG_DIR"] = _CFG_TMP
    os.environ.setdefault("OLLA_DFT_DATA_DIR", tempfile.mkdtemp(prefix="olla-dft-tests-data-"))


def pytest_configure(config):
    config.addinivalue_line(
        "markers", "qe: necesita binarios de Quantum ESPRESSO o salidas ya "
                   "calculadas")


@pytest.fixture(autouse=True)
def _idioma_limpio():
    """Cada prueba empieza en español, sea cual sea el idioma por defecto.

    Las cadenas de origen del código están en español y las pruebas las
    comparan literalmente; el inglés se pide explícitamente con
    ``language="en"`` en las pruebas de paridad. Además, ``main(["...",
    "--language", "en"])`` fija el idioma del proceso y, sin esto, se
    colaría en la siguiente prueba.
    """
    from qekit.core import i18n
    i18n.set_language("es")
    yield
    i18n.set_language(None)


@pytest.fixture(scope="session")
def datadir():
    """Carpeta con salidas de QE ya calculadas, o None."""
    d = os.environ.get("OLLA_DFT_TEST_DATA")
    return Path(d) if d and Path(d).is_dir() else None


@pytest.fixture(scope="session")
def pw_disponible():
    return shutil.which("pw.x") is not None


def necesita(path):
    """Salta la prueba si el dato calculado no está disponible."""
    if path is None or not Path(path).exists():
        pytest.skip(f"falta {path}; define OLLA_DFT_TEST_DATA")
    return Path(path)
