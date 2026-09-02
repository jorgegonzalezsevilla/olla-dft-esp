"""Pruebas sobre salidas REALES de Quantum ESPRESSO.

Los archivos de tests/datos/ son salidas de verdad de QE 6.6 sobre silicio
LDA, recortadas a lo mínimo (216 KB en total). Así estas pruebas corren sin
tener QE instalado y, aun así, ejercitan los parsers y los análisis con los
datos que realmente van a llegarles — que es donde aparecieron casi todos
los errores de este proyecto: unidades, orden de ejes, coordenadas
cartesianas contra fraccionarias, identificación de bandas.

Los valores esperados están en referencias.py, cada uno con la fuente
externa contra la que se validó.
"""

import gzip
import shutil
from pathlib import Path

import numpy as np
import pytest

from tests import referencias as ref

DATOS = Path(__file__).parent / "datos"


@pytest.fixture(scope="session")
def xml_bandas(tmp_path_factory):
    return _descomprimir(DATOS / "bandas_si" / "out" / "Si.xml.gz",
                         tmp_path_factory.mktemp("bandas") / "Si.xml")


@pytest.fixture(scope="session")
def dir_bandas(tmp_path_factory, xml_bandas):
    d = tmp_path_factory.mktemp("bandas_dir")
    (d / "out").mkdir()
    shutil.copy(xml_bandas, d / "out" / "Si.xml")
    kpath = DATOS / "bandas_si" / "KPATH.txt"
    if kpath.exists():          # sin esto no hay etiquetas de alta simetría
        shutil.copy(kpath, d / "KPATH.txt")
    return d


@pytest.fixture(scope="session")
def xml_masa(tmp_path_factory):
    return _descomprimir(DATOS / "mef_si" / "out" / "Si.xml.gz",
                         tmp_path_factory.mktemp("masa") / "Si.xml")


def _descomprimir(origen, destino):
    if not origen.exists():
        pytest.skip(f"falta el dato {origen}")
    with gzip.open(origen, "rb") as f, open(destino, "wb") as g:
        shutil.copyfileobj(f, g)
    return destino


# ----------------------------------------------------------------------
# Lectura del XML y análisis de gap
# ----------------------------------------------------------------------
def test_gap_de_silicio(dir_bandas):
    from qekit.modules import bands as bm
    bs = bm.load(str(dir_bandas))
    info = bm.analyze_gap(bs)
    assert not info.is_metal
    assert info.gap == pytest.approx(ref.SI_GAP_INDIRECTO, abs=0.005)
    assert info.direct_gap == pytest.approx(ref.SI_GAP_DIRECTO, abs=0.005)
    assert not info.is_direct
    assert info.vbm_label == "Γ"          # el VBM del Si está en Gamma


def test_kdist_en_angstrom_inversos(dir_bandas):
    """kdist debe estar en Å⁻¹, no en 2*pi/alat ni en fraccionarias."""
    from qekit.modules import bands as bm
    bs = bm.load(str(dir_bandas))
    # el camino completo de la FCC del Si mide unos 5 Å⁻¹
    assert 3.0 < bs.kdist[-1] < 8.0
    assert np.all(np.diff(bs.kdist) >= -1e-12)   # monótono


# ----------------------------------------------------------------------
# Masa efectiva sobre el cálculo fino
# ----------------------------------------------------------------------
def test_masas_efectivas_de_silicio(xml_masa):
    """El electrón contra experimento; el hueco contra Luttinger."""
    import json

    from qekit.modules import effmass as em
    meta_f = DATOS / "mef_si" / "masa_meta.json"
    if not meta_f.exists():
        pytest.skip("falta masa_meta.json")
    meta = [(x["portador"], x["direccion"], x["npts"], x["kindex"])
            for x in json.loads(meta_f.read_text())["lineas"]]
    run = em.collect_fine(str(xml_masa), meta)

    def masa(portador, direccion, banda=None):
        for f in run.fits:
            if (f.carrier == portador and f.direction == direccion
                    and (banda is None or f.band == banda)):
                return f.mass
        raise AssertionError(f"no se ajustó {portador} en {direccion}")

    assert masa("electrón", "longitudinal") == pytest.approx(
        ref.SI_ME_LONGITUDINAL, abs=0.01)
    assert masa("electrón", "transversal 1") == pytest.approx(
        ref.SI_ME_TRANSVERSAL, abs=0.005)
    # las dos transversales tienen que salir iguales por simetría
    assert masa("electrón", "transversal 1") == pytest.approx(
        masa("electrón", "transversal 2"), abs=1e-6)
    assert masa("hueco", "[100]") == pytest.approx(
        ref.SI_MH_100_PESADO, abs=0.01)
    assert masa("hueco", "[111]") == pytest.approx(
        ref.SI_MH_111_PESADO, abs=0.01)
    # el gap del cálculo fino debe coincidir con el del camino de bandas
    assert run.cbm - run.vbm == pytest.approx(ref.SI_GAP_INDIRECTO, abs=0.005)


def test_bandas_se_identifican_por_conteo_de_electrones(xml_masa):
    """En un cálculo 'bands' no hay E_F utilizable; con nelec sí se acierta.

    Si esto se rompe, las masas salen de bandas equivocadas y el gap del
    reporte sale absurdo — que es exactamente como se detectó el error.
    """
    import json

    from qekit.modules import effmass as em
    meta_f = DATOS / "mef_si" / "masa_meta.json"
    if not meta_f.exists():
        pytest.skip("falta masa_meta.json")
    meta = [(x["portador"], x["direccion"], x["npts"], x["kindex"])
            for x in json.loads(meta_f.read_text())["lineas"]]
    run = em.collect_fine(str(xml_masa), meta)
    huecos = {f.band for f in run.fits if f.carrier == "hueco"}
    electrones = {f.band for f in run.fits if f.carrier == "electrón"}
    assert max(huecos) == 3           # Si: 8 electrones -> 4 bandas llenas
    assert min(electrones) == 4
    assert run.vbm == pytest.approx(6.205, abs=0.01)


# ----------------------------------------------------------------------
# Ópticas sobre la salida real de epsilon.x
# ----------------------------------------------------------------------
@pytest.fixture(scope="session")
def optica():
    from qekit.modules import optics
    d = DATOS / "opt_si"
    if not (d / "epsi_Si.dat").exists():
        pytest.skip("faltan los datos de epsilon.x")
    return optics.collect(optics.OpticsRun(prefix="Si", outdir=d))


def test_regla_de_suma_f(optica):
    """int E*eps2 dE = (pi/2)*(hbar*wp)^2 con la frecuencia de plasmón que
    reporta el propio epsilon.x. Es una verificación independiente de que
    eps2 está bien normalizado."""
    from qekit.core.compat import trapezoid
    E, e2 = optica.energies, optica.eps2
    integral = trapezoid(E * e2, E)
    esperado = np.pi / 2 * ref.SI_PLASMON ** 2
    assert integral / esperado == pytest.approx(1.0, abs=0.01)


def test_kramers_kronig_reproduce_eps1_de_epsilon_x(optica):
    """La KK propia de QEkit contra el eps1 que escribió QE: dos caminos
    independientes al mismo número."""
    from qekit.modules.optics import kramers_kronig
    e1_kk = kramers_kronig(optica.energies, optica.eps2)
    assert e1_kk[0] == pytest.approx(float(optica.eps1[0]), rel=0.01)


def test_eps1_estatico_y_scissor(optica):
    from qekit.modules import optics
    assert optica.eps1[0] == pytest.approx(ref.SI_EPS1_0_SIN_SCISSOR, abs=0.05)
    con = optics.scissor(optica, 0.65)
    assert con.eps1[0] == pytest.approx(ref.SI_EPS1_0_SCISSOR_065, abs=0.05)
    pico = con.energies[int(np.argmax(con.eps2))]
    assert pico == pytest.approx(ref.SI_PICO_EPS2_SCISSOR, abs=0.05)


def test_tauc_sobre_datos_reales_da_el_gap_directo(optica):
    """La extrapolación de Tauc sobre el espectro real del Si debe caer en
    el gap DIRECTO (2.569 eV en LDA), porque epsilon.x no incluye
    transiciones asistidas por fonones."""
    from qekit.modules import optics
    gap, _, _, _ = optics.tauc_gap(optica, "direct")
    assert gap == pytest.approx(ref.SI_GAP_DIRECTO, abs=0.05)
    con = optics.tauc_gap(optics.scissor(optica, 0.65), "direct")[0]
    assert con == pytest.approx(ref.SI_GAP_DIRECTO + 0.65, abs=0.05)


# ----------------------------------------------------------------------
# Fonones sobre la salida real de matdyn.x
# ----------------------------------------------------------------------
@pytest.fixture(scope="session")
def fonones():
    from qekit.modules import phonons
    d = DATOS / "fon_si"
    if not (d / "bandas.freq").exists():
        pytest.skip("faltan los datos de matdyn.x")
    run = phonons.PhononRun(prefix="Si", outdir=d, qgrid=(2, 2, 2))
    q, fr = phonons._read_flfrq(d / "bandas.freq")
    run.band_q, run.band_freqs = q, fr
    datos = np.loadtxt(d / "fonones.dos", comments="#")
    run.dos_w, run.dos = datos[:, 0], datos[:, 1]
    run.qdist = np.arange(len(q), dtype=float)
    return run


def test_frecuencias_de_silicio(fonones):
    """Contra dispersión inelástica de neutrones (Nilsson & Nelin 1972)."""
    fr = fonones.band_freqs
    assert fr.shape[1] == 6                     # 2 átomos -> 6 ramas
    # Gamma: tres acústicas en cero y tres ópticas degeneradas
    assert np.allclose(fr[0, :3], 0.0, atol=0.5)
    assert np.allclose(fr[0, 3:], ref.SI_FONON_GAMMA_TO, atol=1.0)
    assert fr[30, 0] == pytest.approx(ref.SI_FONON_X_TA, abs=1.0)   # X
    assert fr[30, 2] == pytest.approx(ref.SI_FONON_X_LA, abs=1.0)
    assert fr[121, 0] == pytest.approx(ref.SI_FONON_L_TA, abs=1.0)  # L
    assert fr[121, 5] == pytest.approx(ref.SI_FONON_L_TO, abs=1.0)


def test_sin_frecuencias_imaginarias(fonones):
    """Una imaginaria significa estructura mal relajada o inestable."""
    assert fonones.band_freqs.min() > -1.0


def test_degeneraciones_por_simetria(fonones):
    """En X y en L las ramas transversales son dobles: si el parser
    desordenara las ramas, esto se rompería."""
    fr = fonones.band_freqs
    assert fr[30, 0] == pytest.approx(fr[30, 1], abs=0.1)     # X: TA doble
    assert fr[121, 0] == pytest.approx(fr[121, 1], abs=0.1)   # L: TA doble
    assert fr[121, 4] == pytest.approx(fr[121, 5], abs=0.1)   # L: TO doble


def test_termodinamica_armonica(fonones):
    from qekit.modules import phonons
    th = phonons.thermodynamics(fonones, natoms=2)
    assert th["ZPE"] * 1000 == pytest.approx(ref.SI_ZPE_MEV, abs=0.5)
    i300 = int(np.argmin(np.abs(th["T"] - 300.0)))
    assert th["Cv"][i300] * 1000 == pytest.approx(ref.SI_CV_300K, abs=0.005)
    # Cv -> 0 cuando T -> 0, y crece monótonamente
    assert th["Cv"][0] == pytest.approx(0.0, abs=1e-9)
    assert np.all(np.diff(th["Cv"]) >= -1e-12)


def test_dos_de_fonones_integra_a_3N(fonones):
    """La integral de la DOS tiene que dar 3N = 6 modos."""
    from qekit.core.compat import trapezoid
    total = trapezoid(fonones.dos, fonones.dos_w)
    assert total == pytest.approx(6.0, rel=0.05)
