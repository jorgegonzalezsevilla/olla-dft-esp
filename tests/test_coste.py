"""Estimador de coste: la física del modelo y su honestidad."""

import sqlite3

import numpy as np
import pytest

from qekit.modules import cost


# ----------------------------------------------------------------------
# la parte física, que no necesita historial
# ----------------------------------------------------------------------
def test_ondas_planas_contra_quantum_espresso():
    """Celda primitiva de Si (V=39.5 A^3) a 30 Ry: pw.x reporta ~725 PWs."""
    npw = cost.n_ondas_planas(39.53, 30.0)
    assert npw == pytest.approx(725, rel=0.06), f"salieron {npw:.0f}"


def test_ondas_planas_escalan_con_ecut_a_la_tres_medios():
    a = cost.n_ondas_planas(100.0, 20.0)
    b = cost.n_ondas_planas(100.0, 80.0)
    assert b / a == pytest.approx(8.0, rel=1e-9)   # (80/20)^1.5 = 8


def test_ondas_planas_escalan_con_el_volumen():
    a = cost.n_ondas_planas(50.0, 40.0)
    b = cost.n_ondas_planas(150.0, 40.0)
    assert b / a == pytest.approx(3.0, rel=1e-9)


def test_ondas_planas_sin_datos_es_cero():
    assert cost.n_ondas_planas(0, 30) == 0.0
    assert cost.n_ondas_planas(40, 0) == 0.0


def test_k_irreducibles_reduce_una_red_cubica():
    from ase.build import bulk
    a = bulk("Si", "diamond", 5.43)
    n = cost.k_irreducibles(a, (8, 8, 8))
    assert 1 < n < 512, "la simetria tiene que reducir, pero no a un punto"


def test_los_dos_terminos_del_trabajo():
    d = {"nk": 10, "nspin": 1, "npw": 1000.0, "nbnd": 20.0}
    w1, w2 = cost.trabajo(d)
    assert w1 == pytest.approx(10 * 1000.0 * 20.0)
    assert w2 == pytest.approx(w1 * 20.0), "el segundo va con bandas al cuadrado"


def test_el_espin_duplica_el_trabajo():
    base = {"nk": 4, "npw": 500.0, "nbnd": 8.0}
    w1_sin = cost.trabajo({**base, "nspin": 1})[0]
    w1_con = cost.trabajo({**base, "nspin": 2})[0]
    assert w1_con == pytest.approx(2 * w1_sin)


def test_no_colineal_cuesta_cuatro_veces():
    base = {"nk": 4, "npw": 500.0, "nbnd": 8.0}
    assert (cost.trabajo({**base, "nspin": 4})[0]
            == pytest.approx(4 * cost.trabajo({**base, "nspin": 1})[0]))


# ----------------------------------------------------------------------
# leer un input
# ----------------------------------------------------------------------
INPUT = """&CONTROL
  calculation      = 'relax'
  prefix           = 'Si'
/
&SYSTEM
  ibrav            = 0
  nat              = 2
  ntyp             = 1
  ecutwfc          = 30
  ecutrho          = 120
  nbnd             = 9
  nspin            = 2
/
&ELECTRONS
/
&IONS
/

ATOMIC_SPECIES
  Si     28.0850  Si.UPF

ATOMIC_POSITIONS crystal
  Si    0.0000000000   0.0000000000   0.0000000000
  Si    0.2500000000   0.2500000000   0.2500000000

CELL_PARAMETERS angstrom
    0.0000000000   2.7036500000   2.7036500000
    2.7036500000   0.0000000000   2.7036500000
    2.7036500000   2.7036500000   0.0000000000

K_POINTS automatic
  8 8 8 0 0 0
"""


def test_descriptores_de_un_input(tmp_path):
    f = tmp_path / "pw.in"
    f.write_text(INPUT)
    d = cost.descriptores_de_input(f)
    assert d["calculation"] == "relax"
    assert d["nat"] == 2
    assert d["ecutwfc"] == 30.0
    assert d["nbnd"] == 9.0
    assert d["nspin"] == 2
    assert d["volumen"] == pytest.approx(39.53, rel=0.01)
    assert 1 < d["nk"] < 512
    assert d["npw"] > 0


def test_el_nk_de_la_salida_gana_al_estimado(tmp_path):
    """pw.x y spglib no siempre ven la misma simetria; manda pw.x."""
    f = tmp_path / "pw.in"; f.write_text(INPUT)
    d1 = cost.descriptores_de_input(f)
    assert d1["nk_fuente"].startswith("simetría")
    (tmp_path / "pw.out").write_text(
        "     number of k points=    85  Marzari-Vanderbilt smearing\n")
    d2 = cost.descriptores_de_input(f)
    assert d2["nk"] == 85
    assert d2["nk_fuente"] == "el que usó pw.x"


def test_k_points_gamma(tmp_path):
    f = tmp_path / "pw.in"
    f.write_text(INPUT.replace("K_POINTS automatic\n  8 8 8 0 0 0",
                               "K_POINTS gamma"))
    assert cost.descriptores_de_input(f)["nk"] == 1


# ----------------------------------------------------------------------
# calibración
# ----------------------------------------------------------------------
def _base(tmp_path, filas):
    p = tmp_path / "q.db"
    con = sqlite3.connect(str(p))
    con.execute("""CREATE TABLE calculos (
        ruta TEXT PRIMARY KEY, natoms INTEGER, ecutwfc REAL, kgrid TEXT,
        nspin INTEGER, volumen_A3 REAL, n_scf INTEGER, nk INTEGER,
        nbnd INTEGER, n_bfgs INTEGER, wall_s REAL, calculation TEXT)""")
    for i, f in enumerate(filas):
        con.execute("INSERT INTO calculos VALUES (?,?,?,?,?,?,?,?,?,?,?,?)",
                    (f"c{i}", f["nat"], f["ecut"], "4x4x4", f.get("nspin", 1),
                     f["vol"], f.get("n_scf", 10), f["nk"], f["nbnd"],
                     f.get("n_bfgs", 1), f["t"], f.get("calc", "scf")))
    con.commit(); con.close()
    return str(p)


def _sintetico(t0, C1, C2, n=40, semilla=0, n_scf_fijo=None):
    """Cálculos falsos generados con un modelo conocido, para ver si lo recupera.

    Con `n_scf_fijo` todas las corridas usan las mismas iteraciones. Importa:
    el modelo no puede adivinar cuántas iteraciones necesitará un cálculo
    futuro, así que usa la mediana del historial; si en los datos varían,
    aparece una dispersión que es REAL y no un fallo del ajuste.
    """
    rng = np.random.default_rng(semilla)
    filas = []
    for _ in range(n):
        nat = int(rng.integers(2, 60))
        ecut = float(rng.choice([20, 30, 40, 60]))
        vol = 20.0 * nat
        nk = int(rng.integers(1, 60))
        nbnd = max(4, nat * 2)
        n_scf = n_scf_fijo or int(rng.integers(6, 20))
        d = {"nk": nk, "nspin": 1, "npw": cost.n_ondas_planas(vol, ecut),
             "nbnd": float(nbnd)}
        w1, w2 = cost.trabajo(d)
        filas.append(dict(nat=nat, ecut=ecut, vol=vol, nk=nk, nbnd=nbnd,
                          n_scf=n_scf, t=t0 + (C1 * w1 + C2 * w2) * n_scf))
    return filas


def _con_relax(filas, cuantos, pasos, t0, C1, C2):
    """Convierte las primeras filas en relajaciones, con su coste coherente."""
    for f in filas[:cuantos]:
        f["calc"] = "relax"
        f["n_bfgs"] = pasos
        f["t"] = t0 + (f["t"] - t0) * pasos
    return filas


def test_base_vacia_no_calibra(tmp_path):
    m = cost.calibrar(_base(tmp_path, []))
    assert not m.calibrado
    assert "Sin calibrar" in cost.report_modelo(m)


def test_recupera_un_modelo_conocido(tmp_path):
    """Con datos sin ruido generados por t = t0 + C1 w1 + C2 w2, los recupera."""
    t0, C1, C2 = 3.0, 2.0e-6, 5.0e-9
    m = cost.calibrar(_base(tmp_path, _sintetico(t0, C1, C2)))
    assert m.calibrado
    assert m.C1 == pytest.approx(C1, rel=0.02)
    assert m.C2 == pytest.approx(C2, rel=0.05)
    assert m.t0 == pytest.approx(t0, abs=1.0)


def test_con_datos_perfectos_la_dispersion_es_casi_uno(tmp_path):
    m = cost.calibrar(_base(tmp_path,
                            _sintetico(3.0, 2e-6, 5e-9, n_scf_fijo=10)))
    assert m.dispersion < 1.2, f"salio x{m.dispersion:.2f}"


def test_las_iteraciones_variables_meten_dispersion_de_verdad(tmp_path):
    """El modelo no puede saber cuantas iteraciones hara falta: eso se ve."""
    a = tmp_path / "a"; a.mkdir()
    b = tmp_path / "b"; b.mkdir()
    fija = cost.calibrar(_base(a, _sintetico(3.0, 2e-6, 5e-9, n_scf_fijo=10)))
    variable = cost.calibrar(_base(b, _sintetico(3.0, 2e-6, 5e-9)))
    assert variable.dispersion > fija.dispersion


def test_predice_bien_lo_que_no_ha_visto(tmp_path):
    t0, C1, C2 = 3.0, 2.0e-6, 5.0e-9
    m = cost.calibrar(_base(tmp_path, _sintetico(t0, C1, C2, n=40)))
    d = {"nk": 20, "nspin": 1, "npw": cost.n_ondas_planas(600.0, 45.0),
         "nbnd": 60.0, "calculation": "scf"}
    w1, w2 = cost.trabajo(d)
    it = m.n_scf_mediana
    esperado = t0 + C1 * w1 * it + C2 * w2 * it
    assert cost.estimar(d, m)["segundos"] == pytest.approx(esperado, rel=0.05)


def test_historial_poco_variado_se_declara_flojo(tmp_path):
    """Veinte calculos identicos calibran la escala pero no la pendiente."""
    filas = [dict(nat=2, ecut=30, vol=40, nk=10, nbnd=8, n_scf=10, t=5.0)
             for _ in range(20)]
    m = cost.calibrar(_base(tmp_path, filas))
    assert m.calibrado
    assert not m.extrapola_bien
    assert "poco variado" in cost.report_modelo(m)


def test_los_pasos_ionicos_se_aprenden(tmp_path):
    filas = _con_relax(_sintetico(3.0, 2e-6, 5e-9, n=30, n_scf_fijo=10),
                       10, 9, 3.0, 2e-6, 5e-9)
    m = cost.calibrar(_base(tmp_path, filas))
    assert m.pasos_ionicos["relax"] == pytest.approx(9.0)


def test_una_relajacion_cuesta_mas_que_un_scf(tmp_path):
    filas = _con_relax(_sintetico(3.0, 2e-6, 5e-9, n=30, n_scf_fijo=10),
                       12, 8, 3.0, 2e-6, 5e-9)
    m = cost.calibrar(_base(tmp_path, filas))
    d = {"nk": 40, "nspin": 1, "npw": 20000.0, "nbnd": 60.0}
    t_scf = cost.estimar({**d, "calculation": "scf"}, m)["segundos"]
    t_rel = cost.estimar({**d, "calculation": "relax"}, m)["segundos"]
    assert t_rel > 6 * t_scf, f"scf {t_scf:.1f}s, relax {t_rel:.1f}s"


def test_los_coeficientes_nunca_son_negativos(tmp_path):
    """Un C negativo diria que anadir bandas acelera el calculo."""
    rng = np.random.default_rng(7)
    filas = _sintetico(3.0, 2e-6, 5e-9, n=30)
    for f in filas:
        f["t"] *= float(rng.uniform(0.3, 3.0))      # mucho ruido
    m = cost.calibrar(_base(tmp_path, filas))
    assert m.C1 >= 0 and (m.C2 or 0) >= 0 and m.t0 >= 0


def test_humano():
    assert cost.humano(45) == "45 s"
    assert cost.humano(600) == "10 min"
    assert cost.humano(7200) == "2.0 h"
    assert cost.humano(None) == "?"
