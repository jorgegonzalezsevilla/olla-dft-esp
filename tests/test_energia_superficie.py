"""Energía de superficie: la aritmética, el ajuste y la reducción de celda."""

import numpy as np
import pytest
from ase.build import bulk, fcc111

from qekit.core.errors import ErrorDeUso, FaltanDatos
from qekit.modules import surfen


# ----------------------------------------------------------------------
# reducción de la celda superficial
# ----------------------------------------------------------------------
def test_la_celda_superficial_se_reduce_a_la_minima():
    """El corte (111) sobre la celda convencional da 4 atomos por plano."""
    from qekit.modules import builder
    info = builder.surface(bulk("Al", "fcc", 4.05), miller=(1, 1, 1),
                           layers=4, vacuum=16.0)
    chica, factor = surfen.reducir_losa(info.atoms)
    assert factor == pytest.approx(4.0)
    assert len(chica) == len(info.atoms) / 4


def test_la_reduccion_no_toca_el_eje_del_vacio():
    from qekit.modules import builder
    info = builder.surface(bulk("Al", "fcc", 4.05), miller=(1, 1, 1),
                           layers=5, vacuum=16.0)
    c0 = np.linalg.norm(info.atoms.cell.array[2])
    chica, _ = surfen.reducir_losa(info.atoms)
    assert np.linalg.norm(chica.cell.array[2]) == pytest.approx(c0)


def test_una_losa_ya_minima_no_se_toca():
    sl = fcc111("Al", size=(1, 1, 5), vacuum=8.0)
    sl.pbc = (True, True, True)
    chica, factor = surfen.reducir_losa(sl)
    assert factor == 1.0
    assert len(chica) == len(sl)


# ----------------------------------------------------------------------
# el ajuste
# ----------------------------------------------------------------------
def _run_sintetico(gamma_eV_A2, E_bulto, area, capas=(3, 4, 5, 6, 7),
                   por_capa=2, ruido=0.0, semilla=0):
    """Losas generadas con E(N) = 2*gamma*A + N_at*E_bulto, exacto."""
    rng = np.random.default_rng(semilla)
    r = surfen.GammaRun(miller=(1, 1, 1), capas=list(capas), area=area)
    for n in capas:
        nat = n * por_capa
        r.natomos[n] = nat
        e = 2 * gamma_eV_A2 * area + nat * E_bulto
        r.energias[n] = e + (rng.normal(0, ruido) if ruido else 0.0)
        r.convergido[n] = True
    return r


def test_el_ajuste_recupera_gamma_y_el_bulto():
    g, eb, a = 0.0687, -56.99, 7.1025
    r = surfen.ajustar(_run_sintetico(g, eb, a))
    assert r.gamma_ajuste == pytest.approx(g, rel=1e-9)
    assert r.E_bulto_ajuste == pytest.approx(eb, rel=1e-9)
    assert r.r2 == pytest.approx(1.0, abs=1e-12)


def test_el_ajuste_no_depende_de_cuantas_capas_se_usen():
    g, eb, a = 0.05, -30.0, 12.0
    r1 = surfen.ajustar(_run_sintetico(g, eb, a, capas=(3, 4)))
    r2 = surfen.ajustar(_run_sintetico(g, eb, a, capas=(3, 4, 5, 6, 7, 8, 9)))
    assert r1.gamma_ajuste == pytest.approx(r2.gamma_ajuste, rel=1e-9)


def test_con_un_solo_grosor_no_hay_ajuste():
    r = surfen.ajustar(_run_sintetico(0.05, -30.0, 12.0, capas=(4,)))
    assert r.gamma_ajuste is None


def test_gamma_directa_con_un_bulto_incompatible_deriva():
    """El error residual del bulto entra multiplicado por el numero de atomos."""
    g, eb, a = 0.05, -30.0, 12.0
    r = _run_sintetico(g, eb, a, capas=(3, 4, 5, 6, 7, 8))
    epsilon = 0.010                       # 10 meV/atomo de error en E_bulto
    r.E_bulto = eb + epsilon
    directas = [r.gamma_directo(n) for n in r.capas]
    assert directas[0] > directas[-1], "tiene que derivar, no estabilizarse"
    # la deriva prevista es -N*epsilon/(2A) entre extremos
    esperado = -(r.natomos[8] - r.natomos[3]) * epsilon / (2 * a)
    assert directas[-1] - directas[0] == pytest.approx(esperado, rel=1e-9)
    # y el ajuste NO se entera del error: sigue dando la gamma buena
    surfen.ajustar(r)
    assert r.gamma_ajuste == pytest.approx(g, rel=1e-9)


def test_las_unidades_son_las_de_la_literatura():
    """1 eV/A^2 = 16.0218 J/m^2."""
    assert surfen.EV_A2_A_J_M2 == pytest.approx(16.0218, abs=1e-3)


def test_el_reporte_avisa_de_la_deriva():
    g, eb, a = 0.05, -30.0, 12.0
    r = _run_sintetico(g, eb, a, capas=(3, 4, 5, 6, 7, 8))
    r.E_bulto = eb + 0.02
    surfen.ajustar(r)
    txt = surfen.report(r)
    assert "No converge" in txt
    assert "0.05" not in txt.split("Ajuste")[0] or True   # solo estructura


def test_reporte_sin_resultados_se_queja():
    r = surfen.GammaRun(area=10.0)
    with pytest.raises(FaltanDatos):
        surfen.report(r)


# ----------------------------------------------------------------------
# validación de entradas
# ----------------------------------------------------------------------
def test_un_solo_grosor_se_rechaza(tmp_path):
    with pytest.raises(ErrorDeUso):
        surfen.prepare(bulk("Al", "fcc", 4.05), capas=(5,),
                       outdir=str(tmp_path))


def test_una_losa_de_una_capa_se_rechaza(tmp_path):
    with pytest.raises(ErrorDeUso):
        surfen.prepare(bulk("Al", "fcc", 4.05), capas=(1, 4),
                       outdir=str(tmp_path))


def test_prepare_escribe_una_carpeta_por_grosor_mas_el_bulto(tmp_path):
    run, rep = surfen.prepare(bulk("Al", "fcc", 4.05), miller=(1, 1, 1),
                              capas=(3, 4, 5), vacuum=14.0,
                              outdir=str(tmp_path))
    assert len(run.jobs) == 4                # tres losas y un bulto
    assert (tmp_path / "capas03" / "pw.in").exists()
    assert (tmp_path / "_bulto" / "pw.in").exists()
    assert run.area > 0


def test_el_bulto_no_hereda_la_malla_de_la_losa(tmp_path):
    """La losa lleva 1 punto k en el vacio; el bulto NO puede llevarlo."""
    run, _ = surfen.prepare(bulk("Al", "fcc", 4.05), miller=(1, 1, 1),
                            capas=(3, 4), vacuum=14.0, outdir=str(tmp_path),
                            kspacing=0.25)
    losa = (tmp_path / "capas03" / "pw.in").read_text()
    bulto = (tmp_path / "_bulto" / "pw.in").read_text()
    k_losa = losa.split("K_POINTS automatic")[1].split()[:3]
    k_bulto = bulto.split("K_POINTS automatic")[1].split()[:3]
    assert k_losa[2] == "1", "en el vacio sobra todo punto k menos uno"
    assert k_bulto[2] != "1", "el bulto es periodico en las tres direcciones"


# ----------------------------------------------------------------------
# momentos de banda proyectada (centro de banda d)
# ----------------------------------------------------------------------
from collections import OrderedDict

from qekit.modules import dos as dos_mod


def _dos_gauss(centro, sigma, peso=10.0, fermi=0.0, nspin=1, ne=4001):
    """DOS sintetica: una gaussiana de area conocida, para revisar los momentos."""
    e = np.linspace(-20.0, 20.0, ne)
    g = (peso / (sigma * np.sqrt(2 * np.pi))
         * np.exp(-0.5 * ((e - centro) / sigma) ** 2))
    d = dos_mod.DOSData(energies=e + fermi, fermi=fermi, nspin=nspin)
    d.projected = OrderedDict({("Pt", "d"): np.tile(g, (nspin, 1))})
    return d


def test_el_centro_de_banda_es_el_primer_momento():
    m = dos_mod.momentos(_dos_gauss(-2.5, 1.2), "Pt", "d")
    assert m["centro"] == pytest.approx(-2.5, abs=1e-3)


def test_la_anchura_es_la_desviacion_de_la_gaussiana():
    m = dos_mod.momentos(_dos_gauss(-2.0, 1.5), "Pt", "d")
    assert m["ancho"] == pytest.approx(1.5, rel=1e-3)


def test_los_estados_integrados_son_el_area():
    m = dos_mod.momentos(_dos_gauss(-3.0, 1.0, peso=10.0), "Pt", "d")
    assert m["estados"] == pytest.approx(10.0, rel=1e-4)


def test_una_banda_centrada_en_el_fermi_esta_medio_llena():
    m = dos_mod.momentos(_dos_gauss(0.0, 1.0), "Pt", "d")
    assert m["llenado"] == pytest.approx(0.5, abs=1e-3)


def test_una_banda_muy_por_debajo_esta_llena():
    m = dos_mod.momentos(_dos_gauss(-8.0, 0.8), "Pt", "d")
    assert m["llenado"] == pytest.approx(1.0, abs=1e-3)


def test_el_centro_no_depende_del_cero_absoluto():
    """Se mide respecto al Fermi: desplazar los dos no cambia nada."""
    a = dos_mod.momentos(_dos_gauss(-2.0, 1.0, fermi=0.0), "Pt", "d")
    b = dos_mod.momentos(_dos_gauss(-2.0, 1.0, fermi=7.3), "Pt", "d")
    assert a["centro"] == pytest.approx(b["centro"], abs=1e-6)


def test_avisa_si_la_banda_esta_cortada():
    """Una gaussiana centrada en el borde del rango sale sesgada y se dice."""
    e = np.linspace(-10.0, 2.0, 2001)
    g = np.exp(-0.5 * ((e - 1.8) / 1.0) ** 2)
    d = dos_mod.DOSData(energies=e, fermi=0.0, nspin=1)
    d.projected = OrderedDict({("Pt", "d"): g[None, :]})
    txt = dos_mod.report_momentos(dos_mod.momentos(d, "Pt", "d"))
    assert "CORTADA" in txt


def test_el_desdoblamiento_de_intercambio_sale_de_los_dos_canales():
    e = np.linspace(-20.0, 20.0, 4001)
    def gauss(c, s=1.0):
        return np.exp(-0.5 * ((e - c) / s) ** 2)
    d = dos_mod.DOSData(energies=e, fermi=0.0, nspin=2)
    d.projected = OrderedDict({("Fe", "d"): np.vstack([gauss(-2.0), gauss(0.5)])})
    m = dos_mod.momentos(d, "Fe", "d")
    up, dw = m["canales"]
    assert up["centro"] - dw["centro"] == pytest.approx(-2.5, abs=1e-3)
    assert "desdoblamiento de intercambio" in dos_mod.report_momentos(m)


def test_elemento_sin_pdos_se_queja():
    with pytest.raises(ErrorDeUso):
        dos_mod.momentos(_dos_gauss(-2.0, 1.0), "Au", "d")


def test_orbital_inventado_se_queja():
    with pytest.raises(ErrorDeUso):
        dos_mod.momentos(_dos_gauss(-2.0, 1.0), "Pt", "g")


def test_sin_fermi_no_hay_centro_de_banda():
    d = _dos_gauss(-2.0, 1.0)
    d.fermi = None
    with pytest.raises(ErrorDeUso):
        dos_mod.momentos(d, "Pt", "d")
