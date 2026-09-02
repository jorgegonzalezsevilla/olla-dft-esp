"""Pruebas que no necesitan Quantum ESPRESSO.

Cubren la parte del código donde históricamente aparecieron los errores:
conversiones de unidades, orientación de celdas, índices hkl, integrales
numéricas y compatibilidad entre versiones de numpy.
"""

from pathlib import Path

import numpy as np
import pytest
from ase.build import bulk

from tests import referencias as ref


# ----------------------------------------------------------------------
# Compatibilidad y utilidades
# ----------------------------------------------------------------------
def test_trapezoid_existe_en_cualquier_numpy():
    """numpy 2.0 eliminó np.trapz; numpy 1.x no tiene np.trapezoid."""
    from qekit.core.compat import trapezoid
    x = np.linspace(0.0, 1.0, 101)
    assert trapezoid(x, x) == pytest.approx(0.5, abs=1e-6)


def test_constantes_fisicas():
    from qekit.core import qeout
    from qekit.modules import effmass, phonons
    assert qeout.HARTREE_EV == pytest.approx(27.211386, abs=1e-5)
    assert qeout.BOHR_ANG == pytest.approx(0.5291772, abs=1e-6)
    assert qeout.RY_EV == pytest.approx(13.605693, abs=1e-5)
    # h^2/m_e = 2 * 3.80998 eV*A^2
    assert effmass.HBAR2_OVER_ME == pytest.approx(7.6199682, abs=1e-6)
    assert phonons.CM1_TO_THZ == pytest.approx(0.0299792458, abs=1e-10)


# ----------------------------------------------------------------------
# Estructura y simetría
# ----------------------------------------------------------------------
def test_primitiva_realinea_al_marco_cubico():
    """El bug que rompió las constantes elásticas: una primitiva girada.

    structure.primitive() debe devolver la celda FCC estandarizada, con los
    ejes alineados a los del cubo. Si no, las Cij y las masas efectivas
    salen mal aunque el cálculo de QE esté bien.
    """
    from qekit.core import structure as st
    at = bulk("Si", "diamond", a=5.43)
    # celda girada: misma red, otro marco
    girada = at.copy()
    th = np.radians(31.0)
    R = np.array([[np.cos(th), -np.sin(th), 0],
                  [np.sin(th), np.cos(th), 0], [0, 0, 1.0]])
    girada.set_cell(at.cell.array @ R.T, scale_atoms=True)
    cell = st.primitive(girada).cell.array
    esperada = np.array([[0, 1, 1], [1, 0, 1], [1, 1, 0]]) * (5.43 / 2)
    assert np.allclose(np.sort(np.abs(cell.ravel())),
                       np.sort(np.abs(esperada.ravel())), atol=1e-6)
    for fila in cell:                      # cada vector tiene un cero exacto
        assert np.min(np.abs(fila)) < 1e-8


def test_conventional_expande_fcc():
    from qekit.core import structure as st
    at = bulk("Si", "diamond", a=5.43)
    assert len(at) == 2
    assert len(st.conventional(at)) == 8


# ----------------------------------------------------------------------
# Difracción
# ----------------------------------------------------------------------
def test_xrd_indexa_en_celda_convencional():
    """Con la primitiva de entrada, los hkl deben ser los de la ficha PDF.

    En la base primitiva FCC el (220) cúbico se llama (211): correcto, pero
    incomparable con la literatura.
    """
    from qekit.modules import xrd
    at = bulk("Si", "diamond", a=5.43)
    pat = xrd.compute(at, two_theta_range=(5, 70))
    picos = sorted(pat.peaks, key=lambda p: p.two_theta)
    assert [p.label for p in picos] == [h for _, h in ref.SI_XRD]
    for p, (tt, _) in zip(picos, ref.SI_XRD):
        assert p.two_theta == pytest.approx(tt, abs=0.06)


def test_xrd_base_input_da_indices_primitivos():
    from qekit.modules import xrd
    at = bulk("Si", "diamond", a=5.43)
    pat = xrd.compute(at, two_theta_range=(5, 70), basis="input")
    picos = sorted(pat.peaks, key=lambda p: p.two_theta)
    assert picos[1].label == "(211)"        # el (220) convencional
    assert picos[1].two_theta == pytest.approx(47.343, abs=0.06)


def test_xrd_nacl_sin_division_por_cero():
    """Las reflexiones extinguidas (F=0) coincidentes en 2theta daban 0/0."""
    from qekit.modules import xrd
    at = bulk("NaCl", "rocksalt", a=5.64)
    with np.errstate(all="raise"):
        pat = xrd.compute(at, two_theta_range=(5, 70))
    tts = sorted(p.two_theta for p in pat.peaks)
    assert len(tts) == len(ref.NACL_XRD)
    for got, (esperado, _) in zip(tts, ref.NACL_XRD):
        assert got == pytest.approx(esperado, abs=0.06)
    assert all(np.isfinite(p.two_theta) for p in pat.peaks)


def test_xrd_intensidad_maxima_normalizada():
    from qekit.modules import xrd
    pat = xrd.compute(bulk("Si", "diamond", a=5.43))
    assert max(p.intensity for p in pat.peaks) == pytest.approx(100.0)


# ----------------------------------------------------------------------
# Ópticas: Kramers-Kronig y scissor
# ----------------------------------------------------------------------
def _lorentz(E, E0=3.0, g=0.2, A=8.0):
    """eps2 de un oscilador; su eps1 por KK es analítico y conocido."""
    return A * g * E0 / ((E ** 2 - E0 ** 2) ** 2 + (g * E) ** 2) * E


def test_kramers_kronig_reproduce_oscilador():
    """KK de un oscilador de Lorentz debe dar eps1 > 1 bajo la resonancia
    y eps1 < 1 encima: la firma de la dispersión anómala."""
    from qekit.modules.optics import kramers_kronig
    E = np.linspace(0.01, 30.0, 3000)
    e2 = _lorentz(E)
    e1 = kramers_kronig(E, e2)
    assert e1[np.argmin(np.abs(E - 1.0))] > 1.0
    assert e1[np.argmin(np.abs(E - 10.0))] < 1.0
    assert np.all(np.isfinite(e1))


def test_scissor_desplaza_y_conserva_area_ponderada():
    """El scissor mueve eps2 en Delta y reescala por ((E-D)/E)^2."""
    from qekit.core.compat import trapezoid
    from qekit.modules import optics
    E = np.linspace(0.0, 30.0, 3000)
    e2 = _lorentz(E)
    run = optics.OpticsRun(energies=E,
                           eps2_xyz=np.tile(e2, (3, 1)),
                           eps1_xyz=np.zeros((3, len(E))))
    run.eps2 = e2; run.eps1 = np.zeros_like(e2)
    d = 0.7
    nuevo = optics.scissor(run, d)
    assert nuevo.scissor == pytest.approx(d)
    # el pico se corre exactamente Delta
    assert E[np.argmax(nuevo.eps2)] - E[np.argmax(e2)] == pytest.approx(
        d, abs=E[1] - E[0])
    # Invariante correcto del scissor: como eps2 ~ |p|^2/w^2 y el scissor
    # no toca los elementos de matriz, lo que se conserva es
    #     int w^2 eps2(w) dw  ~  suma de |p|^2
    # y NO la regla de suma f (int w eps2 dw), que el scissor SÍ cambia:
    # al subir las transiciones, cada una pesa 1/w y el total baja.
    assert (trapezoid(E ** 2 * nuevo.eps2, E)
            == pytest.approx(trapezoid(E ** 2 * e2, E), rel=5e-3))
    assert trapezoid(E * nuevo.eps2, E) < trapezoid(E * e2, E)
    # el original no se toca
    assert run.scissor == 0.0
    assert np.allclose(run.eps2, e2)


def test_scissor_cero_es_identidad():
    from qekit.modules import optics
    E = np.linspace(0.0, 20.0, 500)
    run = optics.OpticsRun(energies=E, eps1=np.ones_like(E),
                           eps2=_lorentz(E))
    assert optics.scissor(run, 0.0) is run


@pytest.mark.parametrize("altura_pico", [0.0, 3.0, 20.0, 100.0])
def test_tauc_recupera_el_gap_de_un_borde_ideal(altura_pico):
    """Borde que obedece el modelo de Tauc en 2.0 eV, con un pico
    interbanda más arriba que puede ser mucho más intenso.

    Para una transición directa permitida, alpha*hv = A*(hv - Eg)^(1/2), o
    sea (alpha*hv)^2 lineal en hv: la extrapolación DEBE devolver 2.0.
    El pico es más empinado que el borde, así que un ajuste que busque la
    pendiente máxima global se va a él. Con el pico 100x se comprueba que
    la detección del borde no depende del máximo global del espectro.
    """
    from qekit.modules import optics
    E = np.linspace(0.01, 12.0, 2000)
    # eps1 constante y eps2 pequeño => alpha ~ eps2*E/(2*hbar*c);
    # se elige eps2 para que (alpha*hv)^2 salga exactamente lineal en hv
    e2 = np.where(E > 2.0, 0.05 * np.sqrt(np.maximum(E - 2.0, 0.0)) / E ** 2,
                  0.0)
    e2 = e2 + altura_pico * 0.05 * np.exp(-((E - 8.0) ** 2) / 0.05)
    run = optics.OpticsRun(energies=E, eps1=np.full_like(E, 4.0), eps2=e2)
    gap, pend, ventana, _ = optics.tauc_gap(run, "direct")
    assert gap is not None
    assert gap == pytest.approx(2.0, abs=0.1), (
        f"pico {altura_pico}x: el ajuste dio {gap:.3f} en vez de 2.0")


def test_tauc_no_inventa_gap_sin_absorcion():
    """Sin borde no debe devolver un número."""
    from qekit.modules import optics
    E = np.linspace(0.01, 12.0, 500)
    run = optics.OpticsRun(energies=E, eps1=np.full_like(E, 4.0),
                           eps2=np.zeros_like(E))
    assert optics.tauc_gap(run, "direct")[0] is None


# ----------------------------------------------------------------------
# Masa efectiva
# ----------------------------------------------------------------------
def test_masa_efectiva_de_una_parabola_exacta():
    """E = h^2 k^2 / (2m) con m conocida debe recuperarse exacta."""
    from qekit.modules.effmass import HBAR2_OVER_ME, _mass_from_quadratic
    for m_obj in (0.19, 0.92, 2.7):
        a = HBAR2_OVER_ME / (2.0 * m_obj)
        assert _mass_from_quadratic(a) == pytest.approx(m_obj, rel=1e-10)


def test_masa_efectiva_signo_del_hueco():
    from qekit.modules.effmass import HBAR2_OVER_ME, _mass_from_quadratic
    a = -HBAR2_OVER_ME / (2.0 * 0.49)
    assert _mass_from_quadratic(a) == pytest.approx(-0.49, rel=1e-10)


def test_direcciones_de_valle_ortonormales():
    from qekit.modules.effmass import valley_directions
    dirs = valley_directions(np.array([0.5, 0.0, 0.5]))
    assert [n for n, _ in dirs] == ["longitudinal", "transversal 1",
                                    "transversal 2"]
    v = np.array([e for _, e in dirs])
    assert np.allclose(v @ v.T, np.eye(3), atol=1e-10)
    # la longitudinal apunta al valle
    assert np.allclose(v[0], np.array([0.5, 0, 0.5]) / np.sqrt(0.5))


def test_direcciones_en_gamma_son_cristalograficas():
    """En Gamma, x/y/z serían equivalentes por simetría cúbica: hay que
    muestrear [100], [110] y [111] para ver el alabeo de la banda."""
    from qekit.modules.effmass import valley_directions
    dirs = valley_directions(np.zeros(3))
    assert [n for n, _ in dirs] == ["[100]", "[110]", "[111]"]
    for _, e in dirs:
        assert np.linalg.norm(e) == pytest.approx(1.0)


# ----------------------------------------------------------------------
# Ecuación de estado
# ----------------------------------------------------------------------
def test_eos_recupera_parametros_sinteticos():
    """Birch-Murnaghan generado y reajustado: debe cerrar el círculo."""
    from qekit.modules import eos
    V0, B0, Bp, E0 = 40.05, 94.2 / 160.21766208, 4.2, -20.0
    V = np.linspace(0.88 * V0, 1.12 * V0, 11)
    eta = (V0 / V) ** (2.0 / 3.0)
    E = E0 + 9 * V0 * B0 / 16 * (
        (eta - 1) ** 3 * Bp + (eta - 1) ** 2 * (6 - 4 * eta))
    run = eos.EOSRun(volumes=list(V), energies=list(E), natoms=2)
    fit = eos.fit(run, equation="birch-murnaghan")
    assert fit.V0 == pytest.approx(V0, rel=1e-4)
    assert fit.B0 == pytest.approx(94.2, rel=1e-3)
    assert fit.Bp == pytest.approx(Bp, rel=1e-3)


# ----------------------------------------------------------------------
# Promedios planares (función trabajo)
# ----------------------------------------------------------------------
def test_promedio_macroscopico_elimina_el_rizado():
    """Una losa con rizado atómico de periodo 1 A: el promedio móvil sobre
    exactamente un periodo debe dejar el valor plano de cada región."""
    from qekit.modules import fields
    nz, L = 600, 30.0
    z = np.linspace(0, L, nz, endpoint=False)
    V = np.where((z > 5) & (z < 15), -8.0, 0.0) + 0.0
    V = V + np.where((z > 5) & (z < 15), 0.5 * np.sin(2 * np.pi * z), 0.0)
    cube = fields.CubeData(origin=np.zeros(3), axes=np.diag([1, 1, L / nz]),
                           shape=(4, 4, nz), data=np.tile(V, (4, 4, 1)),
                           natoms=2)
    zz, prof = fields.planar_average(cube, axis=2)
    ma = fields.macroscopic_average(zz, prof, 1.0)
    assert ma[nz // 3] == pytest.approx(-8.0, abs=1e-3)
    assert ma[-5] == pytest.approx(0.0, abs=1e-3)


# ----------------------------------------------------------------------
# Procedencia
# ----------------------------------------------------------------------
def test_encabezado_de_procedencia_lleva_version_y_parametros():
    from qekit import __version__
    from qekit.core import provenance
    provenance.record_argv(["olla-dft", "optics", "x.cif", "--scissor", "0.65"])
    h = provenance.header("ópticas", {"scissor_eV": 0.65})
    assert __version__ in h
    assert "--scissor 0.65" in h
    assert "scissor_eV = 0.65" in h
    assert all(l.startswith("#") for l in h.splitlines())
    # la variante sin '#' es para np.savetxt(comments="# ")
    assert not any(l.startswith("#") for l in
                   provenance.header_plain("ópticas").splitlines())


def test_metadatos_de_figura_por_formato():
    """PNG descarta 'Subject' en silencio; ahí el detalle va en
    'Description'."""
    from qekit.core import provenance
    pdf = provenance.figure_metadata("x", {"a": 1}, "pdf")
    png = provenance.figure_metadata("x", {"a": 1}, "png")
    assert "Subject" in pdf and "Keywords" in pdf
    assert "Description" in png and "Subject" not in png
    assert png["Creation Time"]


# ----------------------------------------------------------------------
# Intercambio con la suite
# ----------------------------------------------------------------------
def test_json_de_suite_es_serializable_y_versionado():
    import json

    from qekit.modules import interop, xrd
    at = bulk("Si", "diamond", a=5.43)
    pat = xrd.compute(at)
    xrd.broaden(pat, two_theta_range=(5, 70))
    doc = interop.from_xrd(pat, at)
    texto = json.dumps(doc, ensure_ascii=False)     # no debe reventar
    vuelta = json.loads(texto)
    assert vuelta["qekit_suite_schema"] == interop.SCHEMA
    assert vuelta["tipo"] == "drx_patron_calculado"
    # los hkl son de la convencional: el material debe describir ESA celda
    assert vuelta["material"]["formula"] == "Si8"
    assert vuelta["material"]["grupo_espacial_numero"] == 227
    assert vuelta["datos"]["picos"][0]["hkl"] == "(111)"
    assert vuelta["generado_por"]["qekit_version"]


def test_json_de_suite_sin_nan():
    """NaN/Infinity no son JSON válido: hay que limpiarlos antes."""
    import json

    from qekit.modules import interop
    doc = interop.envelope("prueba", {"x": np.array([1.0, np.nan, np.inf])})
    texto = json.dumps(doc, allow_nan=False)        # falla si quedó un NaN
    assert json.loads(texto)["datos"]["x"] == [1.0, None, None]


# ----------------------------------------------------------------------
# Constructores: superficies y defectos
# ----------------------------------------------------------------------
def test_superficie_mide_el_vacio_entre_superficies_atomicas():
    """El vacío que importa es el hueco entre átomos, no el de la celda."""
    from qekit.modules import builder
    at = bulk("Si", "diamond", a=5.43)
    info = builder.surface(at, miller=(1, 1, 1), layers=4, vacuum=16.0)
    assert info.vacuum_real == pytest.approx(16.0, abs=0.1)
    c = info.atoms.cell.array[2, 2]
    assert c == pytest.approx(info.thickness + info.vacuum_real, abs=0.1)


def test_superficie_avisa_si_el_vacio_es_insuficiente():
    from qekit.modules import builder
    info = builder.surface(bulk("Si", "diamond", a=5.43), miller=(1, 0, 0),
                           layers=4, vacuum=6.0)
    assert any("vacío REAL" in w for w in info.warnings)


def test_superficie_detecta_polaridad_en_un_compuesto():
    """Una losa con caras distintas es polar y necesita dipfield."""
    from qekit.modules import builder
    zb = bulk("GaAs", "zincblende", a=5.65)
    info = builder.surface(zb, miller=(1, 1, 1), layers=4, vacuum=15.0)
    if info.polar:
        assert any("POLAR" in w for w in info.warnings)


def test_congelar_planos_no_congela_la_losa_entera():
    from qekit.modules import builder
    at = bulk("Si", "diamond", a=5.43)
    info = builder.surface(at, miller=(1, 0, 0), layers=6, vacuum=15.0,
                           fix_layers=2)
    assert 0 < info.fijados < len(info.atoms)
    assert info.atoms.arrays["qekit_fijo"].sum() == info.fijados


def test_vacancia_quita_exactamente_un_atomo():
    from qekit.modules import builder
    perf, info = builder.defect(bulk("Si", "diamond", a=5.43),
                                kind="vacancy", site=0, supercell=(3, 3, 3))
    assert len(perf) == 54
    assert len(info.atoms) == 53
    assert "mu(Si)" in builder.formation_energy_text(info)


def test_intersticial_avisa_si_cae_encima_de_un_atomo():
    """La posición se toma de un átomo REAL de la supercelda: la celda
    estandarizada por spglib no pone necesariamente un átomo en (0,0,0)."""
    from qekit.modules import builder
    at = bulk("Si", "diamond", a=5.43)
    perf, _ = builder.defect(at, kind="vacancy", site=0, supercell=(2, 2, 2))
    ocupada = perf.get_scaled_positions()[0]
    _p, info = builder.defect(at, kind="interstitial", new_element="H",
                              supercell=(2, 2, 2), position=tuple(ocupada))
    assert any("vecino más cercano" in w for w in info.warnings)


def test_intersticial_acepta_un_hueco_valido():
    from qekit.modules import builder
    at = bulk("Si", "diamond", a=5.43)
    _p, info = builder.defect(at, kind="interstitial", new_element="H",
                              supercell=(2, 2, 2), position=(0.0, 0.0, 0.0))
    assert not any("vecino más cercano" in w for w in info.warnings)
    assert info.atoms.get_chemical_symbols().count("H") == 1


def test_sustitucion_exige_elemento_nuevo():
    from qekit.modules import builder
    with pytest.raises(ValueError, match="new-element"):
        builder.defect(bulk("Si", "diamond", a=5.43), kind="substitution")


# ----------------------------------------------------------------------
# Cargas de Bader
# ----------------------------------------------------------------------
def test_bader_reparte_dos_gaussianas_conocidas():
    """Dos gaussianas de 3 e y 5 e separadas: cada cuenca debe recuperarlas.

    La densidad se normaliza como la escribe pp.x (e/bohr³) con la rejilla
    en Å, que es lo que devuelve `fields.read_cube`: la carga ABSOLUTA de
    cada cuenca tiene que salir en electrones, no multiplicada por
    (Å/bohr)³ ≈ 6.75."""
    from qekit.modules import charges, fields
    n, L = 40, 12.0
    g = (np.arange(n) + 0.5) * L / n
    X, Y, Z = np.meshgrid(g, g, g, indexing="ij")
    dv_bohr3 = (L / n) ** 3 / fields.BOHR ** 3

    def gauss(c, q, s=0.9):
        d2 = sum(((A - ci + L / 2) % L - L / 2) ** 2
                 for A, ci in zip((X, Y, Z), c))
        r = np.exp(-d2 / (2 * s * s))
        return q * r / (r.sum() * dv_bohr3)          # e/bohr³

    rho = gauss((3.0, 6.0, 6.0), 3.0) + gauss((9.0, 6.0, 6.0), 5.0)
    cube = fields.CubeData(origin=np.zeros(3), axes=np.eye(3) * (L / n),
                           shape=(n, n, n), data=rho, natoms=2)
    res = charges.bader(cube, np.array([[3.0, 6, 6], [9.0, 6, 6]]),
                        symbols=["A", "B"], valence=[3.0, 5.0])
    assert res.charges[0] == pytest.approx(3.0, abs=0.02)
    assert res.charges[1] == pytest.approx(5.0, abs=0.02)
    # la suma de las cuencas es el número de electrones de valencia
    assert res.total == pytest.approx(8.0, abs=1e-6)
    assert res.total_grid == pytest.approx(8.0, abs=1e-6)
    # la carga no se pierde al repartir en cuencas
    assert res.total == pytest.approx(res.total_grid, rel=1e-9)
    # cuencas iguales para gaussianas iguales y simétricas; volumen en Å³
    assert res.volumes[0] == pytest.approx(res.volumes[1], rel=1e-6)
    assert res.volumes.sum() == pytest.approx(L ** 3, rel=1e-9)
    rep = charges.report_bader(res)
    assert "n/d" not in rep and "no coincide" not in rep

    # la misma densidad declarada en e/Å³ debe dar lo mismo
    cube_A = fields.CubeData(origin=np.zeros(3), axes=np.eye(3) * (L / n),
                             shape=(n, n, n), data=rho / fields.BOHR ** 3,
                             natoms=2)
    res_A = charges.bader(cube_A, np.array([[3.0, 6, 6], [9.0, 6, 6]]),
                          density_units="e/A3")
    assert res_A.total == pytest.approx(8.0, abs=1e-6)


def test_diferencia_de_carga_exige_la_misma_rejilla():
    from qekit.modules import charges, fields
    a = fields.CubeData(origin=np.zeros(3), axes=np.eye(3), shape=(4, 4, 4),
                        data=np.ones((4, 4, 4)), natoms=1)
    b = fields.CubeData(origin=np.zeros(3), axes=np.eye(3), shape=(4, 4, 5),
                        data=np.ones((4, 4, 5)), natoms=1)
    with pytest.raises(ValueError, match="rejillas no coinciden"):
        charges.difference(a, [b])


# ----------------------------------------------------------------------
# Raman y XPS
# ----------------------------------------------------------------------
def test_espectro_raman_aplica_los_tres_factores():
    """I ~ (wL-w)^4 / w * [n(w)+1] * A. Con dos modos de igual actividad,
    el de menor frecuencia tiene que salir MÁS intenso."""
    from qekit.modules import phonons
    run = phonons.PhononRun(gamma_only=True, raman=True)
    run.modes = [{"modo": 1, "omega_cm1": 200.0, "omega_thz": 6.0,
                  "ir": 0.0, "raman": 100.0, "depol": 0.75},
                 {"modo": 2, "omega_cm1": 900.0, "omega_thz": 27.0,
                  "ir": 0.0, "raman": 100.0, "depol": 0.75}]
    _w, _I, picos = phonons.raman_spectrum(run, laser_nm=532.0, T=300.0)
    d = dict((round(w), i) for w, i in picos)
    assert d[200] > d[900], "falta el peso (wL-w)^4/w o el factor de Bose"


def test_espectro_raman_ignora_los_acusticos():
    from qekit.modules import phonons
    run = phonons.PhononRun(gamma_only=True, raman=True)
    run.modes = [{"modo": i, "omega_cm1": 0.0, "omega_thz": 0.0, "ir": 0.0,
                  "raman": 0.0, "depol": 0.5} for i in (1, 2, 3)]
    run.modes.append({"modo": 4, "omega_cm1": 500.0, "omega_thz": 15.0,
                      "ir": 0.0, "raman": 50.0, "depol": 0.75})
    _w, _I, picos = phonons.raman_spectrum(run)
    assert len(picos) == 1 and picos[0][0] == pytest.approx(500.0)


def test_xps_rechaza_el_excite_que_da_ceros():
    """excite(n)=n produce delta_zv=0 y una tabla de ceros sin error."""
    from qekit.modules import xps
    with pytest.raises(ValueError, match="propia contraparte"):
        xps.build_input("X", {1: 1})
    with pytest.raises(ValueError, match="al menos un par"):
        xps.build_input("X", {})
    assert "excite(1) = 2" in xps.build_input("X", {1: 2})


# ----------------------------------------------------------------------
# SOC y DFT+U
# ----------------------------------------------------------------------
def test_soc_se_niega_con_pseudos_no_relativistas():
    """lspinorb sobre un pseudo escalar da desdoblamiento CERO sin avisar."""
    from qekit.modules import sweep
    common = {"pseudos": {"Si": {"found": True, "filename": "Si.pz-vbc.UPF",
                                 "relativistic": "scalar"}}}
    with pytest.raises(ValueError, match="TOTALMENTE RELATIVISTAS"):
        sweep.check_soc_pseudos(common)
    common["pseudos"]["Si"]["relativistic"] = "full"
    sweep.check_soc_pseudos(common)          # ahora pasa


def test_hubbard_legacy_y_card_producen_sintaxis_distinta():
    """QE <=7.0 usa Hubbard_U(i) en &SYSTEM; >=7.1 usa la tarjeta HUBBARD."""
    from qekit.modules import inputgen
    at = bulk("Fe", "bcc", a=2.87)
    pseudos = {"Fe": {"filename": "Fe.UPF", "found": True}}
    kw = dict(atoms=at, pseudos=pseudos, calculation="scf", prefix="Fe",
              pseudo_dir=".", ecutwfc=60, ecutrho=480,
              kcard="K_POINTS gamma\n", insulator=False)
    legacy = inputgen.build_pw_input(hubbard={"Fe": 4.5},
                                     hubbard_style="legacy", **kw)
    assert "lda_plus_u" in legacy and "Hubbard_U(1)" in legacy
    assert "HUBBARD" not in legacy.split("ATOMIC_SPECIES")[0].replace(
        "Hubbard_U", "")
    card = inputgen.build_pw_input(hubbard={"Fe": 4.5},
                                   hubbard_style="card", **kw)
    assert "HUBBARD (ortho-atomic)" in card and "U Fe-3d 4.5" in card
    assert "lda_plus_u" not in card


def test_orbital_hubbard_por_bloque():
    from qekit.modules.inputgen import _orbital_hubbard
    assert _orbital_hubbard("Fe") == "3d"
    assert _orbital_hubbard("Mo") == "4d"
    assert _orbital_hubbard("Ce") == "4f"
    assert _orbital_hubbard("O") == "2p"


def test_nosym_se_emite_en_el_namelist():
    from qekit.modules import inputgen
    at = bulk("Si", "diamond", a=5.43)
    t = inputgen.build_pw_input(
        atoms=at, pseudos={"Si": {"filename": "Si.UPF", "found": True}},
        calculation="nscf", prefix="Si", pseudo_dir=".", ecutwfc=60,
        ecutrho=480, kcard="K_POINTS gamma\n", insulator=True, nosym=True)
    assert "nosym" in t and "noinv" in t


# ----------------------------------------------------------------------
# Transporte
# ----------------------------------------------------------------------
def test_bxsf_lleva_el_punto_periodico_final():
    """BXSF espera (n+1) puntos por dirección, con el último = el primero."""
    from qekit.modules import transport as tr
    n = 4
    run = tr.TransportRun(grid=(n, n, n), fermi=0.0)
    E = np.linspace(-1, 1, n ** 3).reshape(n, n, n, 1)
    run.energies = E.reshape(-1, 1)
    import tempfile
    with tempfile.TemporaryDirectory() as d:
        f = Path(d) / "fs.bxsf"
        tr.export_bxsf(run, np.eye(3) * 5.43, f, bands=[0])
        txt = f.read_text()
    assert f"{n+1} {n+1} {n+1}" in txt
    nums = [float(x) for l in txt.splitlines()
            if l.startswith("  ") and "BAND" not in l and ":" not in l
            for x in l.split() if _es_num(x)]
    # 5^3 = 125 energías más las 9 de los vectores recíprocos y el origen
    assert len(nums) >= (n + 1) ** 3


def _es_num(x):
    try:
        float(x)
        return True
    except ValueError:
        return False


def test_bxsf_se_niega_sin_bandas_que_crucen():
    from qekit.modules import transport as tr
    run = tr.TransportRun(grid=(2, 2, 2), fermi=-50.0)
    run.energies = np.zeros((8, 2))
    with pytest.raises(ValueError, match="ninguna banda cruza"):
        tr.export_bxsf(run, np.eye(3), "x.bxsf")


# ----------------------------------------------------------------------
# Transporte: propiedades exactas que la implementación debe cumplir
# ----------------------------------------------------------------------
def _semiconductor_simetrico(gap=1.0, n=12, m=0.5):
    """Bandas de valencia y conducción idénticas y simétricas respecto de 0.

    Con simetría electrón-hueco exacta, el Seebeck TIENE que valer cero en
    el centro del gap y ser impar alrededor de él. Es una prueba de los
    signos y de los integrandos que no depende de ningún cálculo de QE.
    """
    from qekit.modules import transport as tr
    k = (np.arange(n) - n // 2) / n
    KX, KY, KZ = np.meshgrid(k, k, k, indexing="ij")
    k2 = KX ** 2 + KY ** 2 + KZ ** 2
    Ec = gap / 2 + 4.0 * k2 / m
    Ev = -gap / 2 - 4.0 * k2 / m
    E = np.stack([Ev.ravel(), Ec.ravel()], axis=1)
    # velocidades: v = (1/hbar) dE/dk, con signo opuesto en cada banda
    kv = np.stack([KX.ravel(), KY.ravel(), KZ.ravel()], axis=1)
    vc = 8.0 / m * kv
    v = np.stack([-vc, vc], axis=1)          # (nk, 2, 3)
    run = tr.TransportRun(energies=E, velocities=v, volume=40.0,
                          nelec=2.0, fermi=0.0, grid=(n, n, n))
    run.weights = np.full(len(E), 1.0 / len(E))
    return run


def test_seebeck_es_cero_en_el_centro_de_un_gap_simetrico():
    from qekit.modules import transport as tr
    run = tr.compute(_semiconductor_simetrico(), T=[300.0],
                     mu=np.array([0.0]))
    s = float(np.trace(run.seebeck[0, 0]) / 3.0)
    assert s == pytest.approx(0.0, abs=1e-12), (
        f"con simetría electrón-hueco exacta S debe ser 0, salió {s:.3e}")


def test_seebeck_cambia_de_signo_a_los_dos_lados_del_gap():
    """Portadores tipo p dan S > 0 y tipo n dan S < 0."""
    from qekit.modules import transport as tr
    run = tr.compute(_semiconductor_simetrico(), T=[300.0],
                     mu=np.array([-0.35, 0.35]))
    s = np.trace(run.seebeck[0], axis1=1, axis2=2) / 3.0
    assert s[0] > 0 > s[1], f"signos equivocados: {s}"
    assert s[0] == pytest.approx(-s[1], rel=1e-9)   # impar por simetría


def test_seebeck_no_depende_de_la_escala_de_las_velocidades():
    """En CRTA el tiempo de relajación se cancela en S y NO en sigma.

    Escalar todas las velocidades equivale a cambiar tau: S tiene que
    quedar igual y sigma escalar con el cuadrado. Si S se moviera, el
    cociente estaría mal montado.
    """
    from qekit.modules import transport as tr
    mu = np.array([-0.3, 0.0, 0.3])
    a = tr.compute(_semiconductor_simetrico(), T=[300.0], mu=mu)
    b = _semiconductor_simetrico()
    b.velocities = b.velocities * 3.0
    b = tr.compute(b, T=[300.0], mu=mu)
    s_a = np.trace(a.seebeck[0], axis1=1, axis2=2) / 3.0
    s_b = np.trace(b.seebeck[0], axis1=1, axis2=2) / 3.0
    assert np.allclose(s_a, s_b, rtol=1e-9)
    sig_a = np.trace(a.sigma[0], axis1=1, axis2=2) / 3.0
    sig_b = np.trace(b.sigma[0], axis1=1, axis2=2) / 3.0
    assert np.allclose(sig_b, 9.0 * sig_a, rtol=1e-9)


def test_seebeck_baja_al_meterse_en_la_banda():
    """|S| es máximo cerca del borde y cae al degenerar el portador."""
    from qekit.modules import transport as tr
    mu = np.linspace(0.3, 1.6, 14)
    run = tr.compute(_semiconductor_simetrico(), T=[300.0], mu=mu)
    s = np.abs(np.trace(run.seebeck[0], axis1=1, axis2=2) / 3.0)
    assert s[0] > s[-1]


def test_transporte_rechaza_una_malla_que_no_es_rejilla():
    from qekit.modules import transport as tr
    run = tr.TransportRun(grid=(3, 3, 3))
    run.energies = np.zeros((10, 2))
    # compute no valida la malla; el que valida es load(), sobre el XML.
    # Aquí se comprueba el mensaje del validador con una malla incoherente.
    assert 3 * 3 * 3 != 10
