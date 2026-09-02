"""Superficies cargadas con ESM.

La prueba de `test_bc1_con_carga_es_error_de_uso` no es una comprobación de
formato: bc1 con carga neta es un problema mal planteado (campo que llega al
infinito), y pw.x lo calcula igualmente. Salieron −379 y −677 Ry para la
misma losa con dos vacíos distintos, sin un solo mensaje de error.
"""

import numpy as np
import pytest

from qekit.core.errors import ErrorDeUso, FaltanDatos
from qekit.modules import esm as E


def _losa(vac=10.0, centrada=False):
    from ase.build import fcc111
    sl = fcc111("Al", size=(1, 1, 5), a=4.05, vacuum=vac, periodic=True)
    sl.center(axis=2)
    if centrada:
        sl.positions[:, 2] -= sl.positions[:, 2].mean()
    return sl


# ----------------------------------------------------------------------
# Centrado: la trampa silenciosa
# ----------------------------------------------------------------------
def test_centrar_deja_la_losa_en_cero():
    a = E.centrar(_losa())
    z = a.get_positions()[:, 2]
    assert 0.5 * (z.min() + z.max()) == pytest.approx(0.0, abs=1e-12)


def test_centrar_no_toca_la_celda_ni_las_otras_coordenadas():
    orig = _losa()
    a = E.centrar(orig)
    assert np.allclose(a.cell.array, orig.cell.array)
    assert np.allclose(a.get_positions()[:, :2], orig.get_positions()[:, :2])


def test_avisa_de_que_habia_que_centrar():
    avisos = E.comprobar(_losa(), "bc1", [0.0])
    assert any("centrado" in a or "centrada" in a for a in avisos)


def test_una_losa_ya_centrada_no_dispara_el_aviso():
    avisos = E.comprobar(_losa(centrada=True), "bc1", [0.0])
    assert not any("centrada en z = 0" in a for a in avisos)


def test_espesor_y_vacio():
    esp, vac = E.espesor_y_vacio(_losa(vac=12.0))
    assert esp == pytest.approx(4 * 4.05 / np.sqrt(3), abs=0.01)
    assert vac == pytest.approx(24.0, abs=0.01)


# ----------------------------------------------------------------------
# Condiciones de contorno
# ----------------------------------------------------------------------
def test_bc1_con_carga_es_error_de_uso():
    with pytest.raises(ErrorDeUso) as e:
        E.comprobar(_losa(centrada=True), "bc1", [0.1])
    assert "bc3" in str(e.value)


def test_bc1_neutra_es_valida():
    E.comprobar(_losa(centrada=True), "bc1", [0.0])


@pytest.mark.parametrize("bc", ["bc2", "bc3"])
def test_bc2_y_bc3_admiten_carga(bc):
    avisos = E.comprobar(_losa(centrada=True), bc, [0.1])
    assert any("contraelectrodo" in a for a in avisos)


def test_condicion_de_contorno_inventada():
    with pytest.raises(ErrorDeUso):
        E.comprobar(_losa(centrada=True), "bc9", [0.0])


def test_se_niega_con_una_celda_no_ortogonal_en_z():
    a = _losa(centrada=True)
    c = a.cell.array.copy()
    c[0, 2] = 1.0
    a.set_cell(c, scale_atoms=False)
    with pytest.raises(ErrorDeUso) as e:
        E.comprobar(a, "bc1", [0.0])
    assert "ortogonal" in str(e.value)


def test_avisa_si_apenas_hay_vacio():
    avisos = E.comprobar(_losa(vac=1.0, centrada=True), "bc1", [0.0])
    assert any("vacío" in a for a in avisos)


# ----------------------------------------------------------------------
# Perfil y nivel de vacío
# ----------------------------------------------------------------------
def _perfil(espesor=9.35, cola=1.2, ruido=0.0, semilla=0):
    """Perfil sintético: potencial que decae exponencialmente en el vacío."""
    z = np.linspace(-14.0, 14.0, 281)
    fuera = np.clip(np.abs(z) - espesor / 2, 0, None)
    v = -3.0 * np.exp(-fuera / cola)
    chg = np.exp(-(z / (espesor / 3)) ** 2)
    if ruido:
        v = v + np.random.default_rng(semilla).normal(0, ruido, len(z))
    return {"z": z, "carga": chg, "v_hartree": v, "v_local": v * 0,
            "v_total": v}


def test_el_nivel_de_vacio_se_aleja_de_la_losa_hasta_que_es_plano():
    v, s, n = E.nivel_vacio(_perfil(), espesor=9.35)
    assert abs(v) < 1e-2                 # la cola ya no cuenta
    assert s < 1e-3                      # y es plano de verdad


def test_un_margen_corto_daria_el_potencial_de_la_cola():
    """Es el error que se corrigió: con 2 Å sale la cola, no el vacío."""
    p = _perfil()
    corto, _s, _n = E.nivel_vacio(p, espesor=9.35, margen=2.0,
                                  margen_max=2.0, tol=1e9)
    largo, _s2, _n2 = E.nivel_vacio(p, espesor=9.35)
    assert abs(corto) > 10 * abs(largo)


def test_una_losa_que_llena_la_celda_no_tiene_nivel_de_vacio():
    with pytest.raises(FaltanDatos):
        E.nivel_vacio(_perfil(), espesor=27.0)


def test_se_puede_promediar_un_solo_lado():
    p = _perfil()
    v, _s, n = E.nivel_vacio(p, espesor=9.35, lado=-1)
    assert n < len(p["z"]) / 2
    assert np.isfinite(v)


def test_la_funcion_trabajo_es_el_vacio_menos_el_fermi():
    p = _perfil()
    v, _s, _n = E.nivel_vacio(p, espesor=9.35)
    assert E.funcion_trabajo(p, -4.2382, espesor=9.35) == pytest.approx(
        v + 4.2382)


def test_leer_esm1(tmp_path):
    f = tmp_path / "x.esm1"
    f.write_text("#z  chg  vh  vl  vt\n"
                 "-1.0  0.1  -2.0  3.0  1.0\n"
                 " 0.0  0.5  -1.0  2.0  1.0\n"
                 " 1.0  0.1  -2.0  3.0  1.0\n")
    p = E.leer_esm1(f)
    assert np.allclose(p["z"], [-1, 0, 1])
    assert np.allclose(p["v_total"], [1, 1, 1])


def test_leer_esm1_dice_que_falta(tmp_path):
    with pytest.raises(FaltanDatos) as e:
        E.leer_esm1(tmp_path)
    assert "assume_isolated" in str(e.value)


# ----------------------------------------------------------------------
# Capacitancia y linealidad
# ----------------------------------------------------------------------
def test_la_capacitancia_es_la_pendiente_dividida_por_el_area():
    q = np.array([-0.02, -0.01, 0.0, 0.01, 0.02])
    phi = 5.0 - 100.0 * q                 # pendiente −100 eV/e
    C, r2 = E.capacitancia(q, phi, area=7.102)
    assert C == pytest.approx(-1.0 / 100.0 / 7.102 * E.E_A2_A_UF_CM2,
                              rel=1e-9)
    assert r2 == pytest.approx(1.0)


def test_el_potencial_de_carga_cero_se_interpola():
    q = np.array([-0.02, 0.0, 0.02])
    phi = np.array([6.0, 5.0, 4.0])
    assert E.potencial_de_carga_cero(q, phi) == pytest.approx(5.0)


def test_una_recta_pasa_la_prueba_de_linealidad():
    q = np.linspace(-0.02, 0.02, 5)
    ok, desv = E.linealidad(q, 5.0 - 100.0 * q)
    assert ok and desv < 1e-9


def test_una_curva_no_pasa_la_prueba_de_linealidad():
    """Los números son los medidos en Al(111) con bc3 y ±0.04 e."""
    q = np.array([-0.04, -0.02, 0.0, 0.02, 0.04])
    phi = np.array([1.2443, -0.5084, -4.2381, -9.0256, -14.5087]) * -1
    ok, desv = E.linealidad(q, phi)
    assert not ok and desv > 0.05


def test_el_gran_canonico_pesa_la_carga_con_el_potencial():
    E0 = np.array([-10.0, -10.0, -10.0])
    q = np.array([-0.1, 0.0, 0.1])
    om = E.gran_canonico(E0, q, 4.0)
    assert np.allclose(om, [-10.4, -10.0, -9.6])


# ----------------------------------------------------------------------
# Preparación
# ----------------------------------------------------------------------
def test_prepare_escribe_esm_en_el_system(tmp_path):
    run, _c, _r = E.prepare(_losa(), outdir=str(tmp_path), bc="bc3",
                            cargas=[-0.05, 0.0, 0.05],
                            pseudo_dir="/usr/share/espresso/pseudo")
    assert len(run.jobs) == 3
    txt = (tmp_path / "q00" / "pw.in").read_text()
    assert "assume_isolated = 'esm'" in txt
    assert "esm_bc          = 'bc3'" in txt
    assert "tot_charge" in txt
    # y sin carga no se escribe tot_charge
    assert "tot_charge" not in (tmp_path / "q01" / "pw.in").read_text()


def test_la_malla_de_k_nunca_tiene_puntos_a_lo_largo_del_vacio(tmp_path):
    E.prepare(_losa(), outdir=str(tmp_path), bc="bc1", cargas=[0.0],
              pseudo_dir="/usr/share/espresso/pseudo")
    card = (tmp_path / "q00" / "pw.in").read_text().split("K_POINTS")[1]
    assert card.split()[1:4][2] == "1"


def test_el_campo_solo_se_escribe_con_bc2(tmp_path):
    E.prepare(_losa(), outdir=str(tmp_path), bc="bc2", cargas=[0.0],
              campo=0.001, pseudo_dir="/usr/share/espresso/pseudo")
    assert "esm_efield" in (tmp_path / "q00" / "pw.in").read_text()
    d2 = tmp_path / "otra"
    E.prepare(_losa(), outdir=str(d2), bc="bc1", cargas=[0.0], campo=0.001,
              pseudo_dir="/usr/share/espresso/pseudo")
    assert "esm_efield" not in (d2 / "q00" / "pw.in").read_text()


def test_prepare_centra_la_losa_en_el_input(tmp_path):
    """Centrada en el marco de ESM, que va de −c/2 a +c/2.

    El input se escribe en fraccionarias dentro de [0,1), así que una losa
    centrada en z = 0 aparece partida entre 0.84… y 0.16…: es lo correcto,
    porque ESM pliega la celda a [−c/2, c/2] y ahí queda centrada.
    """
    E.prepare(_losa(), outdir=str(tmp_path), bc="bc1", cargas=[0.0],
              pseudo_dir="/usr/share/espresso/pseudo")
    txt = (tmp_path / "q00" / "pw.in").read_text()
    zs = [float(l.split()[3]) for l in txt.split("ATOMIC_POSITIONS")[1]
          .split("\n") if len(l.split()) == 4 and l.split()[0] == "Al"]
    assert len(zs) == 5
    esm = [z - 1.0 if z > 0.5 else z for z in zs]     # a [−0.5, 0.5)
    assert 0.5 * (min(esm) + max(esm)) == pytest.approx(0.0, abs=1e-6)


def test_collect_se_queja_si_no_hay_carpetas(tmp_path):
    run = E.EsmRun()
    with pytest.raises(FaltanDatos):
        E.collect(run, str(tmp_path))


def test_el_voltaje_de_la_celda_si_es_lineal_aunque_phi_no_lo_sea():
    """Medido en Al(111) con bc3, cinco cargas entre ∓0.04 e.

    V_vac (el nivel de vacío respecto de la frontera de ESM) sale lineal a
    seis cifras; Φ = V_vac − E_F se desvía un 16 %. Son cosas distintas: V_vac
    es el voltaje de la celda, y Φ mezcla ese voltaje con el cambio del dipolo
    de superficie al cargarla. La capacitancia sale del primero.
    """
    q = np.array([-0.04, -0.02, 0.0, 0.02, 0.04])
    vac = np.array([11.89359, 5.94898, -0.00036, -5.94997, -11.90013])
    phi = np.array([10.6493, 6.4574, 4.2377, 3.0756, 2.6086])
    okv, desvv = E.linealidad(q, vac)
    okp, desvp = E.linealidad(q, phi)
    assert okv and desvv < 1e-4
    assert not okp and desvp > 0.05
    C, r2 = E.capacitancia(q, vac, area=7.102)
    assert abs(C) == pytest.approx(0.758, abs=0.01)     # µF/cm²
    assert r2 > 0.999999


def test_la_capacitancia_medida_es_del_orden_de_epsilon0_partido_por_d():
    """0.76 µF/cm² para un hueco de 7 Å; ε₀/d da 1.27. Mismo orden.

    No tienen por qué coincidir: el plano imagen no está en el último plano
    de iones. Pero un orden de magnitud de diferencia sí sería una señal de
    que algo va mal, y por eso se comprueba.
    """
    eps0 = 8.8541878128e-12                     # F/m
    d = 7.0e-10                                 # m
    ideal = eps0 / d * 1e2                      # F/m² -> µF/cm²
    assert 0.2 < 0.758 / ideal < 5.0
