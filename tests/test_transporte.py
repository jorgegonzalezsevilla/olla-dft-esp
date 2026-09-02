"""Número de Lorenz y transporte por canal de espín."""

import numpy as np
import pytest

from qekit.modules import transport as tr

HBAR = 6.582119569e-16                       # eV s
ME = 0.51099895e6 / (2.99792458e8) ** 2      # eV s^2 / m^2


def gas_de_electrones(n=20, a=4.0, T=300.0):
    """Gas de electrones libres: E = ħ²k²/2m, v = ħk/m.

    Para un metal degenerado el numero de Lorenz vale exactamente
    L0 = (pi^2/3)(k_B/e)^2 = 2.44e-8. Es el unico caso con respuesta
    analitica, asi que es EL test del modulo.
    """
    cell = np.eye(3) * a
    recip = 2 * np.pi * np.linalg.inv(cell).T
    ks = (np.arange(n) + 0.5) / n - 0.5
    K = np.array(np.meshgrid(ks, ks, ks, indexing="ij")).reshape(3, -1).T
    kSI = (K @ recip) * 1e10
    E = HBAR ** 2 * np.sum(kSI ** 2, axis=1) / (2 * ME)
    v = HBAR * kSI / ME
    run = tr.TransportRun(volume=a ** 3, nelec=2.0,
                          fermi=float(np.median(E)), grid=(n, n, n))
    run.energies = E[:, None]
    run.velocities = v[:, None, :]
    run.weights = np.full(len(E), 2.0 / len(E))
    return tr.compute(run, T=[T], mu=np.array([run.fermi]))


def test_el_gas_de_electrones_libres_da_L0():
    """La comprobación decisiva: L -> 2.44e-8 en un metal degenerado."""
    run = gas_de_electrones()
    L = float(tr.lorenz(run, 0)[0])
    assert L / tr.L0_SOMMERFELD == pytest.approx(1.0, abs=0.10), \
        f"L/L0 = {L / tr.L0_SOMMERFELD:.3f}"


def test_el_seebeck_de_un_metal_es_pequeno():
    run = gas_de_electrones()
    S = float(np.trace(run.seebeck[0][0]) / 3.0)
    assert abs(S) < 60e-6, "un metal degenerado da decenas de uV/K, no cientos"


def test_L_no_depende_de_tau():
    """Multiplicar sigma y kappa por el mismo factor no cambia L."""
    run = gas_de_electrones()
    L1 = float(tr.lorenz(run, 0)[0])
    run.sigma = run.sigma * 7.3
    run.kappa_e = run.kappa_e * 7.3
    assert float(tr.lorenz(run, 0)[0]) == pytest.approx(L1, rel=1e-12)


def test_L_no_depende_de_la_temperatura_en_un_metal():
    a = float(tr.lorenz(gas_de_electrones(T=200.0), 0)[0])
    b = float(tr.lorenz(gas_de_electrones(T=500.0), 0)[0])
    assert a == pytest.approx(b, rel=0.10)


def test_el_metal_no_dispara_la_alarma_de_cancelacion():
    run = gas_de_electrones()
    assert float(tr.cancelacion(run, 0)[0]) > 0.1
    assert "NO TE FÍES" not in tr.report_lorenz(run)


# ----------------------------------------------------------------------
# cancelación catastrófica
# ----------------------------------------------------------------------
def _run_falso(sig, See, kap, T=300.0):
    r = tr.TransportRun(volume=40.0, fermi=0.0, grid=(4, 4, 4))
    n = len(sig)
    r.mu = np.linspace(-1, 1, n)
    r.T = np.array([T])
    eye = np.eye(3)
    r.sigma = np.array([[s * eye for s in sig]])
    r.seebeck = np.array([[s * eye for s in See]])
    r.kappa_e = np.array([[k * eye for k in kap]])
    r.carriers = np.zeros((1, n))
    return r


def test_detecta_la_cancelacion_catastrofica():
    """kappa_e = kappa0 - S^2 sigma T; si de la resta queda el 0.003 %, es ruido."""
    sig, S, T = 6.14e15, 976.5e-6, 300.0
    k0 = S ** 2 * sig * T
    kap = k0 * 3.4e-5                       # lo que sobrevive de verdad
    r = _run_falso([sig] * 3, [S] * 3, [kap] * 3)
    c = tr.cancelacion(r, 0)
    assert c[1] < 1e-3
    assert "NO TE FÍES" in tr.report_lorenz(r)


def test_sin_cancelacion_no_avisa():
    sig, S, T = 1e20, 20e-6, 300.0
    kap = tr.L0_SOMMERFELD * sig * T
    r = _run_falso([sig] * 3, [S] * 3, [kap] * 3)
    txt = tr.report_lorenz(r)
    assert "NO TE FÍES" not in txt
    assert "Wiedemann-Franz dentro de un 15" in txt


def test_por_encima_de_L0_habla_de_bipolar():
    sig, S, T = 1e20, 20e-6, 300.0
    r = _run_falso([sig] * 3, [S] * 3, [3.0 * tr.L0_SOMMERFELD * sig * T] * 3)
    assert "BIPOLAR" in tr.report_lorenz(r)


def test_por_debajo_de_L0_habla_de_no_degenerado():
    sig, S, T = 1e20, 20e-6, 300.0
    r = _run_falso([sig] * 3, [S] * 3, [0.6 * tr.L0_SOMMERFELD * sig * T] * 3)
    txt = tr.report_lorenz(r)
    assert "NO es degenerado" in txt and "0.76" in txt


# ----------------------------------------------------------------------
# dos canales de espín
# ----------------------------------------------------------------------
def _espin(sig_up, sig_dw, S_up, S_dw):
    up = _run_falso([sig_up] * 3, [S_up] * 3, [1.0] * 3)
    dw = _run_falso([sig_dw] * 3, [S_dw] * 3, [1.0] * 3)
    up.fermi = dw.fermi = 0.0
    return tr.TransporteEspin(up=up, dw=dw, it=0)


def test_las_conductancias_se_suman():
    te = _espin(3.0, 1.0, 10e-6, 10e-6)
    assert te.sigma_total[1] == pytest.approx(4.0)


def test_el_seebeck_se_pesa_por_la_conductancia():
    """Media pesada, no aritmetica: 3:1 de conductancia da (3*100+1*(-100))/4."""
    te = _espin(3.0, 1.0, 100e-6, -100e-6)
    assert te.seebeck_total[1] * 1e6 == pytest.approx(50.0)
    media_mala = 0.0
    assert te.seebeck_total[1] * 1e6 != pytest.approx(media_mala)


def test_un_canal_que_no_conduce_no_aporta_termopotencia():
    te = _espin(1.0, 0.0, 10e-6, 9999e-6)
    assert te.seebeck_total[1] * 1e6 == pytest.approx(10.0)


def test_la_polarizacion_va_de_menos_uno_a_uno():
    assert _espin(1.0, 0.0, 0, 0).polarizacion[1] == pytest.approx(1.0)
    assert _espin(0.0, 1.0, 0, 0).polarizacion[1] == pytest.approx(-1.0)
    assert _espin(1.0, 1.0, 0, 0).polarizacion[1] == pytest.approx(0.0)


def test_un_medio_metal_se_reconoce():
    te = _espin(1.0, 0.001, 20e-6, 20e-6)
    assert "medio metal" in tr.report_espin(te)


def test_sin_polarizacion_lo_dice():
    te = _espin(1.0, 1.0, 20e-6, 20e-6)
    assert "no aporta nada" in tr.report_espin(te)


def test_la_termopotencia_de_espin_es_la_diferencia():
    te = _espin(1.0, 1.0, 80e-6, 20e-6)
    assert te.seebeck_de_espin[1] * 1e6 == pytest.approx(60.0)
