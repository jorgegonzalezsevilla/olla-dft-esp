"""Funciones de Wannier: geometría, gauge, dispersión e interpolación.

La prueba que de verdad importa aquí es la del gradiente por diferencias
finitas: el factor 1/N_k y el signo del gradiente de Marzari-Vanderbilt no
se pueden adivinar leyendo el artículo, y con el factor equivocado el
descenso no baja —da tumbos— sin que nada avise. Se comprobó así, y así se
queda comprobado.
"""

import numpy as np
import pytest

from qekit.core.errors import ErrorDeUso, FaltanDatos
from qekit.modules import wannier as W


def _bg(cell):
    return 2 * np.pi * np.linalg.inv(np.asarray(cell, float)).T


# ----------------------------------------------------------------------
# Malla y capas de vecinos
# ----------------------------------------------------------------------
def test_la_malla_va_en_el_orden_de_quantum_espresso():
    k = W.malla_completa((2, 2, 3))
    assert len(k) == 12
    # el último índice corre más rápido, como en K_POINTS crystal
    assert np.allclose(k[0], [0, 0, 0])
    assert np.allclose(k[1], [0, 0, 1 / 3])
    assert np.allclose(k[3], [0, 0.5, 0])


def test_malla_de_cero_es_error_de_uso():
    with pytest.raises(ErrorDeUso):
        W.malla_completa((0, 2, 2))


@pytest.mark.parametrize("nombre,cell,n", [
    ("sc", np.eye(3) * 4.0, (4, 4, 4)),
    ("fcc", np.array([[0, .5, .5], [.5, 0, .5], [.5, .5, 0]]) * 5.43, (4, 4, 4)),
    ("bcc", np.array([[-.5, .5, .5], [.5, -.5, .5], [.5, .5, -.5]]) * 3.0,
     (4, 4, 4)),
    ("hex", np.array([[2.46, 0, 0], [-1.23, 2.13, 0], [0, 0, 6.7]]),
     (6, 6, 2)),
    ("tetragonal", np.diag([3.9, 3.9, 4.1]), (4, 4, 4)),
])
def test_las_capas_cumplen_la_condicion_de_completitud(nombre, cell, n):
    capas, pesos = W.capas_b(_bg(cell), n)
    assert W.residuo_completitud(capas, pesos) < 1e-10, nombre


def test_el_peso_de_la_capa_unica_de_una_malla_cubica_es_el_analitico():
    # con seis vecinos a distancia b0, Σ_b w b_a b_b = 2 w b0² δ_ab = δ_ab
    a = 4.0
    n = (4, 4, 4)
    capas, pesos = W.capas_b(_bg(np.eye(3) * a), n)
    b0 = 2 * np.pi / (a * n[0])
    assert len(capas) == 1
    assert sum(len(l) for _, l in capas) == 6
    assert pesos[0] == pytest.approx(1.0 / (2 * b0 ** 2), rel=1e-12)


def test_no_se_guardan_capas_de_peso_nulo():
    # celda muy alargada: la primera capa (a lo largo de c) no aporta
    capas, pesos = W.capas_b(_bg(np.diag([3.0, 3.0, 12.0])), (4, 4, 4))
    assert all(abs(w) > W.TOL_PESO for w in pesos)
    assert W.residuo_completitud(capas, pesos) < 1e-10


def test_los_vecinos_cierran_exactamente():
    cell = np.array([[0, .5, .5], [.5, 0, .5], [.5, .5, 0]]) * 5.43
    n = (4, 4, 4)
    k = W.malla_completa(n)
    capas, _ = W.capas_b(_bg(cell), n)
    idx, G, _ = W.vecinos(k, capas)
    hs = [h for _, lst in capas for h, _ in lst]
    for ik in range(len(k)):
        for ib, h in enumerate(hs):
            destino = k[idx[ik, ib] - 1] + G[ik, ib]
            assert np.allclose(destino, k[ik] + np.array(h) / np.array(n))


def test_cada_vector_b_tiene_su_opuesto():
    cell = np.diag([4.0, 4.0, 4.0])
    capas, _ = W.capas_b(_bg(cell), (4, 4, 4))
    bs = [tuple(np.round(b, 8)) for _, lst in capas for _, b in lst]
    for b in bs:
        assert tuple(-np.array(b) + 0.0) in [tuple(np.round(-np.array(x), 8))
                                             for x in bs]


# ----------------------------------------------------------------------
# Celda de Wigner-Seitz
# ----------------------------------------------------------------------
@pytest.mark.parametrize("cell,n", [
    (np.eye(3) * 4.0, (3, 3, 3)),
    (np.array([[0, .5, .5], [.5, 0, .5], [.5, .5, 0]]) * 5.43, (4, 4, 4)),
    (np.diag([3.0, 3.0, 8.0]), (4, 4, 2)),
])
def test_la_suma_de_pesos_de_wigner_seitz_da_el_numero_de_puntos_k(cell, n):
    R, deg = W.celda_wigner_seitz(cell, n)
    assert float((1.0 / deg).sum()) == pytest.approx(np.prod(n), abs=1e-9)


def test_wigner_seitz_es_simetrica_bajo_inversion():
    R, deg = W.celda_wigner_seitz(np.eye(3) * 4.0, (3, 3, 3))
    tabla = {tuple(r): d for r, d in zip(R, deg)}
    for r, d in tabla.items():
        assert tabla[tuple(-np.array(r))] == d


# ----------------------------------------------------------------------
# Interpolación: contra una fórmula cerrada
# ----------------------------------------------------------------------
def test_interpola_exactamente_una_cadena_de_enlace_fuerte():
    # E(k) = e0 + 2t cos(2πk): la transformada tiene solo R = 0, ±1, así que
    # una malla de 8 puntos la reproduce EXACTA en cualquier k, no aproximada
    e0, t, n = -1.5, 0.7, (8, 1, 1)
    cell = np.diag([2.5, 20.0, 20.0])
    k = W.malla_completa(n)
    Hk = (e0 + 2 * t * np.cos(2 * np.pi * k[:, 0]))[:, None, None] + 0j
    R, deg = W.celda_wigner_seitz(cell, n)
    HR = W.a_reales(Hk, k, R)
    kx = np.linspace(0, 1, 37)
    kk = np.column_stack([kx, np.zeros_like(kx), np.zeros_like(kx)])
    E = W.interpolar(HR, R, deg, kk)[:, 0]
    assert np.abs(E - (e0 + 2 * t * np.cos(2 * np.pi * kx))).max() < 1e-10


def test_interpola_exactamente_una_red_cuadrada_a_segundos_vecinos():
    t1, t2, n = 1.0, 0.25, (10, 10, 1)
    cell = np.diag([3.0, 3.0, 20.0])
    k = W.malla_completa(n)
    ex = lambda kk: (-2 * t1 * (np.cos(2 * np.pi * kk[:, 0]) +
                                np.cos(2 * np.pi * kk[:, 1]))
                     - 4 * t2 * np.cos(2 * np.pi * kk[:, 0]) *
                     np.cos(2 * np.pi * kk[:, 1]))
    Hk = ex(k)[:, None, None] + 0j
    R, deg = W.celda_wigner_seitz(cell, n)
    HR = W.a_reales(Hk, k, R)
    rng = np.random.default_rng(1)
    kk = np.column_stack([rng.random(40), rng.random(40), np.zeros(40)])
    assert np.abs(W.interpolar(HR, R, deg, kk)[:, 0] - ex(kk)).max() < 1e-10


def test_la_interpolacion_reproduce_siempre_los_puntos_de_partida():
    # sea cual sea H(k), volver a los puntos de la malla tiene que ser exacto
    rng = np.random.default_rng(7)
    n, nw = (4, 4, 2), 3
    cell = np.diag([4.0, 4.0, 6.0])
    k = W.malla_completa(n)
    A = rng.normal(size=(len(k), nw, nw)) + 1j * rng.normal(size=(len(k), nw, nw))
    Hk = 0.5 * (A + np.conj(np.transpose(A, (0, 2, 1))))
    R, deg = W.celda_wigner_seitz(cell, n)
    HR = W.a_reales(Hk, k, R)
    assert np.abs(np.sort(W.interpolar(HR, R, deg, k), axis=1)
                  - np.sort(np.linalg.eigvalsh(Hk), axis=1)).max() < 1e-10


def test_hr_es_hermitico_bajo_R_a_menos_R():
    rng = np.random.default_rng(11)
    n, nw = (4, 4, 4), 2
    k = W.malla_completa(n)
    A = rng.normal(size=(len(k), nw, nw)) + 1j * rng.normal(size=(len(k), nw, nw))
    Hk = 0.5 * (A + np.conj(np.transpose(A, (0, 2, 1))))
    R, deg = W.celda_wigner_seitz(np.eye(3) * 4.0, n)
    HR = W.a_reales(Hk, k, R)
    tabla = {tuple(r): i for i, r in enumerate(R)}
    for r, i in tabla.items():
        j = tabla[tuple(-np.array(r))]
        assert np.abs(HR[i] - HR[j].conj().T).max() < 1e-10


# ----------------------------------------------------------------------
# Un juego de solapes M sintético pero MATEMÁTICAMENTE consistente
# ----------------------------------------------------------------------
def _sistema_sintetico(n=(4, 4, 4), D=6, nw=3, semilla=0):
    """Estados de Bloch inventados, pero con la estructura correcta.

    Se toma una familia suave de matrices unitarias V(k) de dimensión D y se
    usan sus primeras nw columnas como |u_nk⟩. Entonces M^{k,b} = V(k)† V(k+b)
    ES el solape exacto de unos estados ortonormales de verdad, así que todo
    lo que se calcule con él (centros, dispersión, gradiente) es consistente
    aunque la física sea inventada. Es la única forma de probar el gradiente
    sin necesitar Quantum ESPRESSO.
    """
    rng = np.random.default_rng(semilla)
    cell = np.eye(3) * 4.0
    k = W.malla_completa(n)
    capas, pes = W.capas_b(_bg(cell), n)
    idx, _, b = W.vecinos(k, capas)
    wb = W.pesos_por_b(capas, pes)
    G = [rng.normal(size=(D, D)) + 1j * rng.normal(size=(D, D))
         for _ in range(3)]
    G = [0.5 * (g + g.conj().T) for g in G]
    V = np.empty((len(k), D, D), complex)
    for ik, kk in enumerate(k):
        Hk = sum(np.sin(2 * np.pi * kk[i]) * G[i] for i in range(3)) \
            + sum(np.cos(2 * np.pi * kk[i]) * G[i].T.conj() for i in range(3))
        w_, v_ = np.linalg.eigh(Hk)
        V[ik] = v_
    M = np.empty((len(k), len(b), nw, nw), complex)
    for ik in range(len(k)):
        for ib in range(len(b)):
            M[ik, ib] = V[ik][:, :nw].conj().T @ V[idx[ik, ib] - 1][:, :nw]
    U = np.tile(np.eye(nw, dtype=complex), (len(k), 1, 1))
    return M, U, idx, b, wb


def test_omega_se_parte_en_I_mas_D_mas_OD():
    M, U, idx, b, wb = _sistema_sintetico()
    d = W.dispersion(W._m_gauge(M, U, idx), b, wb)
    assert d.omega == pytest.approx(d.omega_I + d.omega_D + d.omega_OD,
                                    abs=1e-10)


def test_las_tres_piezas_de_omega_son_no_negativas():
    M, U, idx, b, wb = _sistema_sintetico(semilla=4)
    d = W.dispersion(W._m_gauge(M, U, idx), b, wb)
    assert d.omega_I > -1e-12 and d.omega_D > -1e-12 and d.omega_OD > -1e-12


def test_el_gradiente_coincide_con_la_derivada_numerica():
    """El factor 1/N_k y el signo, comprobados contra la propia Ω.

    Sin esta prueba el gradiente puede estar 64 veces mal en una malla
    4×4×4 y parecer correcto: sigue siendo antihermítico, sigue dejando Ω_I
    quieto, y solo se nota en que el descenso no baja.
    """
    rng = np.random.default_rng(2)
    M, U, idx, b, wb = _sistema_sintetico(n=(3, 3, 3), semilla=5)
    nk, nw = U.shape[0], U.shape[2]
    # apartarse un poco de la identidad para que el gradiente no sea trivial
    dW0 = np.array([(lambda z: (z - z.conj().T) / 2)(
        0.2 * (rng.normal(size=(nw, nw)) + 1j * rng.normal(size=(nw, nw))))
        for _ in range(nk)])
    U = W._rotar(U, dW0)
    G = W._gradiente(W._m_gauge(M, U, idx), b, wb)
    dW = np.array([(lambda z: (z - z.conj().T) / 2)(
        rng.normal(size=(nw, nw)) + 1j * rng.normal(size=(nw, nw)))
        for _ in range(nk)])
    Om = lambda Ux: W.dispersion(W._m_gauge(M, Ux, idx), b, wb).omega
    eps = 1e-5
    num = (Om(W._rotar(U, eps * dW)) - Om(W._rotar(U, -eps * dW))) / (2 * eps)
    ana = sum(float(np.real(np.trace(G[i].conj().T @ dW[i])))
              for i in range(nk))
    assert num == pytest.approx(ana, rel=1e-5, abs=1e-9)


def test_el_gradiente_es_antihermitico():
    M, U, idx, b, wb = _sistema_sintetico(semilla=6)
    G = W._gradiente(W._m_gauge(M, U, idx), b, wb)
    assert np.abs(G + np.conj(np.transpose(G, (0, 2, 1)))).max() < 1e-10


def test_rotar_conserva_la_unitariedad():
    rng = np.random.default_rng(3)
    U = np.tile(np.eye(4, dtype=complex), (5, 1, 1))
    dW = np.array([(lambda z: (z - z.conj().T) / 2)(
        rng.normal(size=(4, 4)) + 1j * rng.normal(size=(4, 4)))
        for _ in range(5)])
    Un = W._rotar(U, dW)
    for u in Un:
        assert np.abs(u.conj().T @ u - np.eye(4)).max() < 1e-12


def test_minimizar_baja_omega_y_deja_quieto_omega_I():
    M, U, idx, b, wb = _sistema_sintetico(semilla=8)
    d0 = W.dispersion(W._m_gauge(M, U, idx), b, wb)
    Um, d1, hist = W.minimizar(M, U, idx, b, wb, pasos=200)
    assert d1.omega < d0.omega
    assert np.all(np.diff(hist) <= 1e-12)          # monótona, nunca sube
    assert abs(d1.omega_I - d0.omega_I) < 1e-9     # Ω_I es invariante
    assert d1.omega_D + d1.omega_OD < d0.omega_D + d0.omega_OD


def test_minimizar_no_rompe_la_unitariedad():
    M, U, idx, b, wb = _sistema_sintetico(semilla=9)
    Um, _, _ = W.minimizar(M, U, idx, b, wb, pasos=60)
    for u in Um:
        assert np.abs(u.conj().T @ u - np.eye(u.shape[1])).max() < 1e-10


def test_la_gauge_de_proyeccion_sale_unitaria():
    rng = np.random.default_rng(12)
    A = rng.normal(size=(9, 6, 4)) + 1j * rng.normal(size=(9, 6, 4))
    U, sv = W.gauge_proyeccion(A)
    assert U.shape == (9, 6, 4)
    for u in U:
        assert np.abs(u.conj().T @ u - np.eye(4)).max() < 1e-12
    assert sv > 0


def test_gauge_de_proyeccion_deja_intacta_una_A_ya_unitaria():
    rng = np.random.default_rng(13)
    z = rng.normal(size=(4, 4)) + 1j * rng.normal(size=(4, 4))
    q, r = np.linalg.qr(z)
    q = q * (np.diag(r) / np.abs(np.diag(r)))
    U, sv = W.gauge_proyeccion(q[None, :, :])
    assert np.abs(U[0] - q).max() < 1e-10
    assert sv == pytest.approx(1.0, abs=1e-10)


# ----------------------------------------------------------------------
# Orbitales de prueba y el archivo .nnkp
# ----------------------------------------------------------------------
def _si():
    from ase.build import bulk
    return bulk("Si", "diamond", 5.43)


def test_una_proyeccion_por_atomo_y_por_orbital():
    p = W.proyecciones("Si:sp3", _si())
    assert len(p) == 8                       # 2 átomos × 4 híbridos
    assert {x.l for x in p} == {-3}
    assert sorted(x.mr for x in p) == [1, 1, 2, 2, 3, 3, 4, 4]


def test_se_pueden_encadenar_varias_proyecciones():
    p = W.proyecciones("Si:s;Si:p", _si())
    assert len(p) == 2 * (1 + 3)


def test_proyeccion_en_coordenadas_explicitas():
    p = W.proyecciones("f=0.125,0.125,0.125:s", _si())
    assert len(p) == 1
    assert np.allclose(p[0].centro, (0.125, 0.125, 0.125))


def test_auto_pone_s_y_p_en_cada_atomo():
    p = W.proyecciones("auto", _si())
    assert len(p) == 2 * 4


@pytest.mark.parametrize("mala", [
    "Si", "Si:sp9", "Ge:sp3", "f=0.1,0.2:s", "f=a,b,c:s",
])
def test_proyecciones_mal_escritas_son_error_de_uso(mala):
    with pytest.raises(ErrorDeUso):
        W.proyecciones(mala, _si())


def test_el_nnkp_se_puede_volver_a_leer(tmp_path):
    si = _si()
    p = W.proyecciones("Si:sp3", si)
    info = W.escribir_nnkp(tmp_path / "x.nnkp", si, (4, 4, 2), p,
                           excluir=(5, 6))
    cell, n, kpts, ex = W._leer_nnkp(tmp_path / "x.nnkp")
    assert np.allclose(cell, si.cell.array, atol=1e-6)
    assert n == (4, 4, 2)
    assert ex == (5, 6)
    assert len(kpts) == 32 == info["nk"]


def test_el_nnkp_lleva_los_bloques_en_el_orden_que_los_lee_pw2wannier90():
    # pw2wannier90 NO rebobina el archivo entre bloque y bloque
    si = _si()
    from pathlib import Path as _P
    p = W.proyecciones("Si:s", si)
    W.escribir_nnkp("/tmp/_orden.nnkp", si, (2, 2, 2), p)
    txt = _P("/tmp/_orden.nnkp").read_text()
    orden = ["real_lattice", "recip_lattice", "kpoints", "projections",
             "nnkpts", "exclude_bands"]
    pos = [txt.index("begin " + b) for b in orden]
    assert pos == sorted(pos)


def test_el_nnkp_declara_tantos_vecinos_como_escribe():
    si = _si()
    p = W.proyecciones("Si:s", si)
    info = W.escribir_nnkp("/tmp/_cuenta.nnkp", si, (3, 3, 3), p)
    from pathlib import Path as _P
    lin = _P("/tmp/_cuenta.nnkp").read_text().split("\n")
    i = lin.index("begin nnkpts")
    j = lin.index("end nnkpts")
    assert int(lin[i + 1]) == info["nnb"]
    assert j - i - 2 == info["nk"] * info["nnb"]


# ----------------------------------------------------------------------
# Ida y vuelta por el formato _hr.dat
# ----------------------------------------------------------------------
def test_el_hr_que_escribimos_lo_volvemos_a_leer(tmp_path):
    rng = np.random.default_rng(21)
    n, nw = (3, 3, 3), 2
    cell = np.eye(3) * 4.0
    k = W.malla_completa(n)
    A = rng.normal(size=(len(k), nw, nw)) + 1j * rng.normal(size=(len(k), nw, nw))
    Hk = 0.5 * (A + np.conj(np.transpose(A, (0, 2, 1))))
    R, deg = W.celda_wigner_seitz(cell, n)
    run = W.WannierRun(HR=W.a_reales(Hk, k, R), R=R, deg=deg, cell=cell,
                       nwann=nw, malla=n, nk=len(k))
    W.export(run, str(tmp_path))
    H2, R2, d2 = W.leer_hr(tmp_path / "WANNIER_hr.dat")
    assert np.array_equal(R2, R)
    assert np.array_equal(d2, deg)
    assert np.abs(H2 - run.HR).max() < 1e-5      # el archivo lleva 6 decimales


def test_collect_se_queja_si_no_hay_nada(tmp_path):
    with pytest.raises(FaltanDatos):
        W.collect(str(tmp_path))


def test_collect_avisa_si_falta_lo_que_escribe_pw2wannier90(tmp_path):
    si = _si()
    W.escribir_nnkp(tmp_path / "si.nnkp", si, (2, 2, 2),
                    W.proyecciones("Si:s", si))
    with pytest.raises(FaltanDatos) as e:
        W.collect(str(tmp_path))
    assert "pw2wannier90" in str(e.value)


# ----------------------------------------------------------------------
# Lo que sale del modelo
# ----------------------------------------------------------------------
def test_la_dos_interpolada_integra_el_numero_de_funciones():
    from qekit.core.compat import trapezoid
    e0, t, n = 0.0, 1.0, (6, 6, 6)
    cell = np.eye(3) * 3.0
    k = W.malla_completa(n)
    Hk = (e0 - 2 * t * np.cos(2 * np.pi * k).sum(axis=1))[:, None, None] + 0j
    R, deg = W.celda_wigner_seitz(cell, n)
    run = W.WannierRun(HR=W.a_reales(Hk, k, R), R=R, deg=deg, cell=cell,
                       nwann=1, malla=n, nk=len(k))
    e, d = W.dos_interpolada(run, malla=12, sigma=0.1)
    assert float(trapezoid(d, e)) == pytest.approx(1.0, abs=1e-3)


def test_asignar_pone_la_funcion_del_enlace_en_el_enlace():
    si = _si()
    centro = np.array([0.125, 0.125, 0.125]) @ si.cell.array
    (sitio, dist), = W.asignar(centro[None, :], si)
    assert "–" in sitio                     # es un enlace, no un átomo
    assert dist < 1e-6


def test_asignar_pone_la_funcion_atomica_en_el_atomo():
    si = _si()
    (sitio, dist), = W.asignar(si.get_positions()[1][None, :], si)
    assert sitio == "Si2" and dist < 1e-8


def test_el_decaimiento_sale_ordenado_por_distancia():
    rng = np.random.default_rng(31)
    R, deg = W.celda_wigner_seitz(np.eye(3) * 4.0, (3, 3, 3))
    HR = rng.normal(size=(len(R), 2, 2)) + 0j
    d, a = W.decaimiento(HR, R, np.eye(3) * 4.0)
    assert np.all(np.diff(d) >= -1e-12)
    assert len(d) == len(R)


# ----------------------------------------------------------------------
# Preparación
# ----------------------------------------------------------------------
def test_prepare_escribe_los_cuatro_pasos(tmp_path):
    run, common, rep = W.prepare(
        _si(), outdir=str(tmp_path), malla=(2, 2, 2), proy="Si:sp3",
        pseudo_dir="/usr/share/espresso/pseudo")
    for f in ("1_scf.in", "2_nscf.in", "3_pw2wan.in", "4_bands.in",
              "correr.sh"):
        assert (tmp_path / f).exists(), f
    assert list(tmp_path.glob("*.nnkp")) and list(tmp_path.glob("*.win"))
    assert run.nwann == 8


def test_el_nscf_lleva_nosym_y_la_malla_completa(tmp_path):
    W.prepare(_si(), outdir=str(tmp_path), malla=(2, 2, 2), proy="Si:s",
              pseudo_dir="/usr/share/espresso/pseudo")
    txt = (tmp_path / "2_nscf.in").read_text()
    assert "nosym" in txt and "noinv" in txt
    assert "K_POINTS crystal" in txt
    assert "\n8\n" in txt                    # 2×2×2 = 8 puntos, sin reducir


def test_no_se_pueden_pedir_mas_wannier_que_bandas(tmp_path):
    with pytest.raises(ErrorDeUso):
        W.prepare(_si(), outdir=str(tmp_path), malla=(2, 2, 2),
                  proy="Si:sp3", nbnd=4,
                  pseudo_dir="/usr/share/espresso/pseudo")


def test_proyecciones_automaticas_avisan(tmp_path):
    run, _, _ = W.prepare(_si(), outdir=str(tmp_path), malla=(2, 2, 2),
                          proy="auto",
                          pseudo_dir="/usr/share/espresso/pseudo")
    assert any("automáticas" in a for a in run.avisos)


def test_el_win_declara_las_mismas_funciones_que_el_nnkp(tmp_path):
    W.prepare(_si(), outdir=str(tmp_path), malla=(2, 2, 2), proy="Si:sp3",
              excluir=(9, 10), nbnd=10,
              pseudo_dir="/usr/share/espresso/pseudo")
    win = next(tmp_path.glob("*.win")).read_text()
    assert "num_wann        = 8" in win
    assert "num_bands       = 10" in win
    assert "exclude_bands   = 9,10" in win


# ----------------------------------------------------------------------
# CLI
# ----------------------------------------------------------------------
@pytest.mark.parametrize("txt,esp", [
    ("5-8", (5, 6, 7, 8)),
    ("3", (3,)),
    ("1,2,5-7", (1, 2, 5, 6, 7)),
    ("", ()),
    (None, ()),
    ("5-8,5-6", (5, 6, 7, 8)),
])
def test_rangos_de_bandas(txt, esp):
    from qekit.cli import _rango_bandas
    assert _rango_bandas(txt) == esp


@pytest.mark.parametrize("mala", ["8-5", "a-b", "x"])
def test_rangos_de_bandas_malos(mala):
    from qekit.cli import _rango_bandas
    with pytest.raises(ErrorDeUso):
        _rango_bandas(mala)


# ----------------------------------------------------------------------
# Desenredado (Souza-Marzari-Vanderbilt)
# ----------------------------------------------------------------------
def _sistema_enredado(n=(3, 3, 3), D=8, nb=6, nw=3, semilla=0):
    """Como el sintético de arriba, pero con MÁS bandas que funciones.

    Devuelve además energías crecientes por banda, para poder poner
    ventanas con sentido.
    """
    rng = np.random.default_rng(semilla)
    cell = np.eye(3) * 4.0
    k = W.malla_completa(n)
    capas, pes = W.capas_b(_bg(cell), n)
    idx, _, b = W.vecinos(k, capas)
    wb = W.pesos_por_b(capas, pes)
    G = [rng.normal(size=(D, D)) + 1j * rng.normal(size=(D, D))
         for _ in range(3)]
    G = [0.5 * (g + g.conj().T) for g in G]
    V = np.empty((len(k), D, D), complex)
    E = np.empty((len(k), nb))
    for ik, kk in enumerate(k):
        Hk = sum(np.sin(2 * np.pi * kk[i]) * G[i] for i in range(3))
        Hk = 0.5 * (Hk + Hk.conj().T)
        w_, v_ = np.linalg.eigh(Hk)
        V[ik] = v_
        E[ik] = w_[:nb]
    M = np.empty((len(k), len(b), nb, nb), complex)
    for ik in range(len(k)):
        for ib in range(len(b)):
            M[ik, ib] = V[ik][:, :nb].conj().T @ V[idx[ik, ib] - 1][:, :nb]
    A = rng.normal(size=(len(k), nb, nw)) + 1j * rng.normal(size=(len(k), nb, nw))
    return M, A, E, idx, b, wb


def test_la_ventana_selecciona_por_energia():
    E = np.array([[-5.0, -1.0, 2.0, 8.0]])
    assert list(W.ventana(E, (-6, 3))[0]) == [True, True, True, False]
    assert W.ventana(E, None).all()


def test_una_ventana_al_reves_es_error_de_uso():
    with pytest.raises(ErrorDeUso):
        W.ventana(np.zeros((1, 3)), (5.0, 1.0))


def test_omega_I_no_depende_del_gauge():
    """Ω_I es invariante: rotar U en cada k no puede cambiarla.

    Es lo que separa el problema del desenredado (elegir el subespacio) del
    de la localización (elegir la base dentro de él). Si esta prueba falla,
    los dos problemas están mezclados.
    """
    rng = np.random.default_rng(5)
    M, A, E, idx, b, wb = _sistema_enredado()
    U, _h, _m, _f = W.gauge_desenredo(M, A, idx, wb, pasos=20)
    o1 = W.omega_I(M, U, idx, wb)
    Ur = np.empty_like(U)
    for ik in range(len(U)):
        z = rng.normal(size=(U.shape[2],) * 2) + 1j * rng.normal(
            size=(U.shape[2],) * 2)
        q, r = np.linalg.qr(z)
        Ur[ik] = U[ik] @ (q * (np.diag(r) / np.abs(np.diag(r))))
    assert W.omega_I(M, Ur, idx, wb) == pytest.approx(o1, abs=1e-10)


def test_el_desenredado_baja_omega_I_y_no_sube_nunca():
    M, A, E, idx, b, wb = _sistema_enredado(semilla=2)
    U, hist, _m, _f = W.gauge_desenredo(M, A, idx, wb, pasos=100)
    assert hist[-1] < hist[0]
    assert np.all(np.diff(hist) <= 1e-12)


def test_el_subespacio_desenredado_es_ortonormal():
    M, A, E, idx, b, wb = _sistema_enredado(semilla=3)
    U, _h, _m, _f = W.gauge_desenredo(M, A, idx, wb, pasos=30)
    for u in U:
        assert np.abs(u.conj().T @ u - np.eye(u.shape[1])).max() < 1e-10


def test_las_bandas_congeladas_quedan_dentro_del_subespacio():
    """Lo que promete la ventana congelada, comprobado sobre el proyector.

    Si la banda i está congelada en el punto k, el proyector del subespacio
    tiene que dejarla intacta: P|i⟩ = |i⟩. Es lo que garantiza que esas
    bandas se reproduzcan EXACTAS, y en el silicio se midió 3·10⁻¹³ eV.
    """
    M, A, E, idx, b, wb = _sistema_enredado(semilla=4)
    lo = float(E.min()) - 1.0
    # el corte se elige por DEBAJO de la tercera banda en todos los k, así
    # que congela como mucho dos: caben en las tres funciones de Wannier
    corte = float(E[:, 2].min()) - 1e-9
    U, _h, _m, fro = W.gauge_desenredo(M, A, idx, wb, E=E,
                                       exterior=(lo, float(E.max()) + 1),
                                       congelada=(lo, corte), pasos=40)
    assert fro.sum() > 0
    for ik in range(len(U)):
        P = U[ik] @ U[ik].conj().T
        for i in np.where(fro[ik])[0]:
            e = np.zeros(U.shape[1], complex); e[i] = 1.0
            assert np.abs(P @ e - e).max() < 1e-9


def test_una_ventana_exterior_demasiado_estrecha_es_error_de_uso():
    M, A, E, idx, b, wb = _sistema_enredado(nw=4, semilla=6)
    with pytest.raises(ErrorDeUso) as e:
        W.gauge_desenredo(M, A, idx, wb, E=E,
                          exterior=(float(E.min()) - 1, float(E.min()) - 0.5))
    assert "ventana" in str(e.value)


def test_congelar_mas_bandas_que_funciones_de_wannier_es_error_de_uso():
    M, A, E, idx, b, wb = _sistema_enredado(nw=2, semilla=7)
    with pytest.raises(ErrorDeUso) as e:
        W.gauge_desenredo(M, A, idx, wb, E=E,
                          exterior=(float(E.min()) - 1, float(E.max()) + 1),
                          congelada=(float(E.min()) - 1, float(E.max()) + 1))
    assert "congelada" in str(e.value)


def test_el_desenredado_deja_un_gauge_de_partida_suave():
    """Tras elegir el subespacio hay que volver a proyectar.

    El gauge que sale de diagonalizar Z es arbitrario. Sin la reproyección,
    Ω arrancaba en 51.8 Å² en el silicio y el descenso posterior no lo
    arreglaba; con ella arranca en 14.3.
    """
    M, A, E, idx, b, wb = _sistema_enredado(semilla=8)
    U, _h, _m, _f = W.gauge_desenredo(M, A, idx, wb, pasos=40)
    d = W.dispersion(W._m_gauge(M, U, idx), b, wb)
    # el gauge de partida ya tiene Ω_D pequeño frente a Ω_OD si es suave
    assert d.omega == pytest.approx(d.omega_I + d.omega_D + d.omega_OD,
                                    abs=1e-9)
    assert d.omega_D < d.omega


@pytest.mark.parametrize("txt,esp", [
    ("-10:20", (-10.0, 20.0)),
    ("0:5", (0.0, 5.0)),
    ("-3.5:-1.25", (-3.5, -1.25)),
    (None, None),
    ("", None),
])
def test_ventanas_del_cli(txt, esp):
    from qekit.cli import _ventana
    assert _ventana(txt) == esp


@pytest.mark.parametrize("mal", ["10:2", "abc", "10", "1:2:3", "1;2"])
def test_ventanas_del_cli_mal_escritas(mal):
    from qekit.cli import _ventana
    with pytest.raises(ErrorDeUso):
        _ventana(mal)
