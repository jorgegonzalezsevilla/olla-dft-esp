"""Polarización por fase de Berry.

Dos de estas pruebas existen porque el código estuvo mal:

  - `test_el_cuanto_se_parte_por_la_mitad_...`: con alguna valencia impar el
    cuanto es e·R/Ω y no 2e·R/Ω. Suponerlo siempre 2 duplica el margen y
    esconde los saltos de rama.
  - `test_la_carga_de_born_usa_la_proyeccion_...`: Z* sale de la proyección
    de ΔP sobre el vector recíproco, no del módulo de R. Con el módulo salía
    2.46 en BN cúbico donde la respuesta es 2.01, y el número seguía
    pareciendo razonable.
"""

import numpy as np
import pytest

from qekit.core.errors import ErrorDeUso, FaltanDatos
from qekit.modules import berry as B


def _si(d=0.0):
    from ase.build import bulk
    a = bulk("Si", "diamond", 5.43)
    if d:
        a.positions[1] += np.array([0.0, 0.0, float(d)])
    return a


def _bn():
    """BN cúbico en la celda que sale al pasar por un CIF (triangular).

    Es la misma estructura que `bulk("BN","zincblende")` pero con los ejes
    girados, y es la que se midió con pw.x. La orientación importa: la fase
    de Berry se mide a lo largo de b_gdir, así que Δφ y el desplazamiento
    tienen que venir de la MISMA celda.
    """
    from ase import Atoms
    cell = [[2.5561910140, 0.0, 0.0],
            [1.2780955070, 2.2137263550, 0.0],
            [1.2780955070, 0.7379087850, 2.0871212231]]
    return Atoms("BN", cell=cell, pbc=True,
                 scaled_positions=[(0, 0, 0), (0.25, 0.25, 0.25)])


# ----------------------------------------------------------------------
# Aritmética de fases
# ----------------------------------------------------------------------
@pytest.mark.parametrize("x,esp", [
    (0.0, 0.0), (0.5, 0.5), (1.5, -0.5), (2.0, 0.0),
    (-1.5, 0.5), (3.25, -0.75),
])
def test_plegado_al_intervalo_de_qe(x, esp):
    assert B.a_intervalo(x) == pytest.approx(esp, abs=1e-12)


def test_medio_cuanto_sale_negativo_como_en_quantum_espresso():
    """El borde: pw.x escribe −1.00000 para el silicio, no +1.00000."""
    assert B.a_intervalo(1.0) == pytest.approx(-1.0, abs=1e-12)
    assert B.a_intervalo(-1.0) == pytest.approx(1.0, abs=1e-12)
    # el silicio cae JUSTO en el borde: +1 y −1 son la misma fase, y cuál de
    # los dos sale depende del último bit del producto. Lo que sí es fijo es
    # que vale medio cuanto, y que la diferencia con QE es cero módulo 2
    f = B.fase_ionica(_si(), {"Si": 4.0}, 3)
    assert abs(f) == pytest.approx(1.0, abs=1e-9)
    assert abs(B.a_intervalo(f - (-1.0))) < 1e-9


def test_desenrollar_recupera_una_fase_que_crece():
    fisica = np.linspace(0, 3.6, 19)
    plegada = B.a_intervalo(fisica)
    seguida, saltos = B.desenrollar(plegada)
    assert np.allclose(seguida - seguida[0], fisica - fisica[0], atol=1e-12)
    assert saltos.max() < 0.3


def test_desenrollar_no_puede_con_pasos_de_medio_cuanto():
    """Si el paso es medio cuanto el seguimiento es ambiguo, y hay que verlo.

    No es un fallo del código: es el límite del método. Por eso se devuelven
    los saltos, para poder avisar.
    """
    fisica = np.array([0.0, 1.1, 2.2, 3.3])
    seguida, saltos = B.desenrollar(B.a_intervalo(fisica))
    assert not np.allclose(seguida - seguida[0], fisica - fisica[0])


def test_desenrollar_con_un_solo_punto_no_revienta():
    s, saltos = B.desenrollar([0.3])
    assert len(s) == 1 and len(saltos) == 0


def test_desenrollar_vacio():
    s, saltos = B.desenrollar([])
    assert len(s) == 0 and len(saltos) == 0


# ----------------------------------------------------------------------
# Parte iónica: se conoce exactamente
# ----------------------------------------------------------------------
def test_la_fase_ionica_del_silicio_es_la_de_quantum_espresso():
    # QE, con el Si a 0.227901 en fraccionarias: 0.91160
    assert B.fase_ionica(_si(0.12), {"Si": 4.0}, 3) == pytest.approx(
        0.911602, abs=1e-5)


def test_la_fase_ionica_del_bn_pliega_cada_ion_por_separado():
    """Con Z impar QE pliega la contribución de CADA ion módulo 1.

    Sin ese plegado por ion sale 1.25 donde QE da 0.25, y la comprobación
    contra QE falla aunque los dos números describan la misma estructura.
    """
    assert B.fase_ionica(_bn(), {"B": 3.0, "N": 5.0}, 3) == pytest.approx(
        0.25, abs=1e-9)
    crudo = B.fase_ionica(_bn(), {"B": 3.0, "N": 5.0}, 3, plegar=False)
    assert crudo == pytest.approx(1.25, abs=1e-9)


def test_la_fase_ionica_de_una_estructura_centrosimetrica_es_media_o_cero():
    f = B.fase_ionica(_si(), {"Si": 4.0}, 3)
    assert abs(B.a_intervalo(2 * f)) < 1e-9      # 2φ ≡ 0: φ es 0 o medio cuanto


def test_sin_valencia_no_hay_fase_ionica():
    with pytest.raises(FaltanDatos):
        B.fase_ionica(_si(), {"Ge": 4.0}, 3)


# ----------------------------------------------------------------------
# El cuanto
# ----------------------------------------------------------------------
def test_el_cuanto_se_parte_por_la_mitad_con_valencias_impares():
    assert B.modulo_de({"Si": 4.0}) == 2.0
    assert B.modulo_de({"B": 3.0, "N": 5.0}) == 1.0
    assert B.modulo_de({"Ga": 3.0, "As": 5.0}) == 1.0
    assert B.modulo_de({"Zn": 12.0, "O": 6.0}) == 2.0


def test_el_cuanto_del_silicio_coincide_con_el_de_quantum_espresso():
    # QE escribió: P = -7.2557731 (mod 14.5115464) (e/Omega).bohr,
    # o sea 0.0537251 e/bohr² = 0.1918558 e/Å²
    q, _ = B.cuanto(_si(0.12), 3, modulo=2.0)
    assert q == pytest.approx(0.1918558, rel=1e-5)


def test_el_cuanto_del_bn_coincide_con_el_de_quantum_espresso():
    # QE escribió: (mod 4.8305009) (e/Omega).bohr con MOD_TOT = 1
    q, _ = B.cuanto(_bn(), 3, modulo=1.0)
    esperado = 4.8305009 * 0.529177210903 / abs(np.linalg.det(_bn().cell.array))
    assert q == pytest.approx(esperado, rel=1e-4)


def test_el_cuanto_impar_es_la_mitad_del_par():
    a, _ = B.cuanto(_bn(), 3, modulo=1.0)
    b, _ = B.cuanto(_bn(), 3, modulo=2.0)
    assert b == pytest.approx(2 * a)


# ----------------------------------------------------------------------
# Cuerdas de puntos k
# ----------------------------------------------------------------------
def test_la_cuerda_recorre_toda_la_zona_de_brillouin():
    ks = B.cuerdas(9, gdir=3, kperp=(2, 2))
    assert len(ks) == 2 * 2 * 9
    # el último de cada cuerda es el primero más un vector completo
    for i in range(4):
        c = ks[i * 9:(i + 1) * 9]
        assert c[-1][2] - c[0][2] == pytest.approx(1.0)
        d = np.diff(c[:, 2])
        assert np.allclose(d, 1.0 / 8)          # 1/(nppstr−1), no 1/nppstr


def test_el_indice_de_la_cuerda_corre_despacio():
    ks = B.cuerdas(5, gdir=1, kperp=(2, 3))
    assert len(ks) == 30
    assert np.allclose(ks[:5, 0], np.linspace(0, 1, 5))
    assert np.allclose(ks[:5, 1], 0.0) and np.allclose(ks[:5, 2], 0.0)


def test_la_direccion_de_la_cuerda_es_la_pedida():
    for g in (1, 2, 3):
        ks = B.cuerdas(5, gdir=g, kperp=(2, 2))
        var = [len(np.unique(np.round(ks[:, i], 8))) for i in range(3)]
        assert var[g - 1] == 5


@pytest.mark.parametrize("mal", [2, 1, 0, -3])
def test_una_cuerda_demasiado_corta_es_error_de_uso(mal):
    with pytest.raises(ErrorDeUso):
        B.cuerdas(mal)


def test_kperp_mal_escrito_es_error_de_uso():
    with pytest.raises(ErrorDeUso):
        B.cuerdas(7, kperp=(6, 6, 6))
    with pytest.raises(ErrorDeUso):
        B.cuerdas(7, kperp=(0, 6))


# ----------------------------------------------------------------------
# Lectura de la salida de pw.x
# ----------------------------------------------------------------------
_SALIDA = """
                            POLARIZATION CALCULATION
       G-vector along string (2 pi/a):  0.70711  0.70711 -0.70711
       Modulus of the vector (1/bohr):  1.06058
       Number of k-points per string:   9
       Number of different strings  :  36

                               SUMMARY OF PHASES

                          Ionic Phase:  0.91160
                     Electronic Phase:  0.08840
                          TOTAL PHASE:  1.00000 MOD_TOT: 2

        The calculation of phases done along the direction of vector 3

           P =   7.2557729  (mod  14.5115464)  (e/Omega).bohr
"""


def test_lee_las_fases_de_una_salida_de_pw(tmp_path):
    f = tmp_path / "2_berry.out"
    f.write_text(_SALIDA)
    b = B.leer_berry(f)
    assert b.fase_ion == pytest.approx(0.91160)
    assert b.fase_el == pytest.approx(0.08840)
    assert b.fase_total == pytest.approx(1.0)
    assert b.modulo == 2.0
    assert b.gdir == 3 and b.nppstr == 9 and b.ncuerdas == 36
    assert b.P_bohr == pytest.approx(7.2557729)


def test_traduce_el_fallo_de_las_cuerdas(tmp_path):
    f = tmp_path / "2_berry.out"
    f.write_text("POLARIZATION CALCULATION\n Error in routine c_phase (1):\n"
                 "     Wrong k-strings weights?\n")
    with pytest.raises(FaltanDatos) as e:
        B.leer_berry(f)
    assert "nosym" in str(e.value)


def test_se_queja_si_no_hay_seccion_de_polarizacion(tmp_path):
    (tmp_path / "x.out").write_text("JOB DONE.\n")
    with pytest.raises(FaltanDatos):
        B.leer_berry(tmp_path)


def test_lee_las_valencias_de_la_tabla_de_pw(tmp_path):
    f = tmp_path / "1_scf.out"
    f.write_text("     atomic species   valence    mass     pseudopotential\n"
                 "        B              3.00    10.81100     B( 1.00)\n"
                 "        N              5.00    14.00700     N( 1.00)\n\n")
    assert B.valencias_de(f) == {"B": 3.0, "N": 5.0}


# ----------------------------------------------------------------------
# Carga efectiva de Born: la proyección correcta
# ----------------------------------------------------------------------
def _run_falso(atoms, fases, u, gdir=3, modulo=2.0):
    run = B.BerryRun(formula=atoms.get_chemical_formula(), gdir=gdir)
    run.estructuras = [atoms] * len(fases)
    run.lambdas = list(np.linspace(0, 1, len(fases)))
    run.puntos = [B.Berry(gdir=gdir, fase_total=f, fase_ion=0.0, fase_el=f,
                          modulo=modulo) for f in fases]
    run.es_desplazamiento = True
    run.desplazamiento = np.asarray(u, float)
    return run


def test_la_carga_de_born_usa_la_proyeccion_sobre_el_vector_reciproco():
    """Z* = 2π·Δφ/(u·B_gdir), no (Ω/e)·ΔP/|u|.

    Los números son los del BN cúbico medido con pw.x: la fase total va de
    0.25212 a 0.34825 al desplazar el boro 0.10 Å. Con la proyección salen
    2.01 e, que es el valor de la literatura (1.9); con el módulo de R salían
    2.46, un factor √2 de más que en zinc-blenda pasa desapercibido.
    """
    at = _bn()
    run = _run_falso(at, [0.25212, 0.29911, 0.34825], [0.0, 0.0, 0.10],
                     modulo=1.0)
    an = B.analizar(run)
    Bg = (2 * np.pi * np.linalg.inv(at.cell.array).T)[2]
    esperado = 2 * np.pi * (0.34825 - 0.25212) / (0.10 * Bg[2])
    assert an["zeff"] == pytest.approx(esperado, rel=1e-6)
    assert an["zeff"] == pytest.approx(2.006, abs=0.01)   # literatura ~1.9


def test_una_fase_que_no_se_mueve_da_carga_de_born_cero():
    at = _si(0.0)
    run = _run_falso(at, [-1.0] * 5, [0.0, 0.0, 0.16])
    an = B.analizar(run)
    assert an["zeff"] == pytest.approx(0.0, abs=1e-12)
    assert an["dP"] == pytest.approx(0.0, abs=1e-12)


def test_un_desplazamiento_perpendicular_avisa_en_vez_de_inventar():
    at = _si(0.0)
    cell = at.cell.array
    Bg = (2 * np.pi * np.linalg.inv(cell).T)[2]
    perp = np.cross(Bg, [1.0, 0.0, 0.0])
    perp = perp / np.linalg.norm(perp) * 0.1
    run = _run_falso(at, [-1.0, -1.0, -1.0], perp)
    an = B.analizar(run)
    assert "zeff" not in an
    assert any("perpendicular" in a for a in run.avisos)


def test_avisa_cuando_un_paso_del_camino_es_demasiado_grande():
    at = _si(0.0)
    run = _run_falso(at, [0.0, 0.9, 1.8], [0.0, 0.0, 0.1])
    B.analizar(run)
    assert any("del cuanto" in a for a in run.avisos)


def test_el_cuanto_del_analisis_sale_del_mod_tot_de_la_salida():
    at = _bn()
    run = _run_falso(at, [0.1, 0.2], [0.0, 0.0, 0.1], modulo=1.0)
    an = B.analizar(run)
    assert an["modulo"] == 1.0
    assert an["cuanto_eA2"] == pytest.approx(B.cuanto(at, 3, 1.0)[0])


def test_los_marcadores_de_pw_de_la_figura_usan_el_mod_tot_real():
    """La figura dividía la fase por 2 aunque MOD_TOT fuera 1.

    En BN cúbico (valencias 3 y 5) el cuanto es la mitad y una fase de QE
    vale fase/1 cuantos, no fase/2: los marcadores "lo que escribe pw.x"
    quedaban a la mitad de la rama seguida. Con fases pequeñas (sin salto de
    rama) el plegado tiene que coincidir con la rama seguida exactamente.
    """
    at = _bn()
    fases = [0.10, 0.15, 0.20]
    run = _run_falso(at, fases, [0.0, 0.0, 0.1], modulo=1.0)
    an = B.analizar(run)
    plegada = B.polarizacion_plegada(fases, an)
    assert plegada == pytest.approx(an["P"], rel=1e-12)
    assert plegada == pytest.approx(np.array(fases) / 1.0 * an["cuanto_cm2"])
    # y con valencias pares sigue valiendo la mitad de la fase
    run2 = _run_falso(_si(0.0), fases, [0.0, 0.0, 0.1], modulo=2.0)
    an2 = B.analizar(run2)
    assert B.polarizacion_plegada(fases, an2) == pytest.approx(an2["P"])


# ----------------------------------------------------------------------
# Caminos y preparación
# ----------------------------------------------------------------------
def test_el_camino_interpola_por_la_imagen_mas_cercana():
    a = _si()
    b = a.copy()
    b.positions[1] += np.array([0.0, 0.0, 0.1])
    b.wrap()
    med = B._interpolar_estructuras(a, b, 0.5)
    d = med.get_positions()[1] - a.get_positions()[1]
    assert np.linalg.norm(d) < 0.06        # medio paso, no media celda


def test_no_se_interpola_entre_estructuras_distintas():
    from ase.build import bulk
    with pytest.raises(ErrorDeUso):
        B._interpolar_estructuras(_si(), bulk("Ge", "diamond", 5.65), 0.5)


def test_prepare_pone_lberry_en_control_y_nosym(tmp_path):
    run, _c, rep = B.prepare(_si(0.1), outdir=str(tmp_path), gdir=3,
                             nppstr=7, kperp=(2, 2),
                             pseudo_dir="/usr/share/espresso/pseudo")
    txt = (tmp_path / "p00" / "2_berry.in").read_text()
    cab = txt.split("&SYSTEM")[0]
    assert "lberry" in cab and "gdir" in cab and "nppstr" in cab
    assert "nosym" in txt and "noinv" in txt
    assert "K_POINTS crystal" in txt


def test_prepare_hace_un_punto_por_lambda(tmp_path):
    run, _c, _r = B.prepare(_si(0.1), outdir=str(tmp_path),
                            referencia=_si(0.0), nlambda=4, nppstr=5,
                            kperp=(2, 2),
                            pseudo_dir="/usr/share/espresso/pseudo")
    assert len(run.jobs) == 4
    for i in range(4):
        assert (tmp_path / f"p{i:02d}" / "2_berry.in").exists()


def test_un_punto_suelto_avisa_de_que_no_significa_nada(tmp_path):
    run, _c, _r = B.prepare(_si(0.1), outdir=str(tmp_path), nppstr=5,
                            kperp=(2, 2),
                            pseudo_dir="/usr/share/espresso/pseudo")
    assert len(run.jobs) == 1
    assert any("no significa nada" in a for a in run.avisos)


def test_no_se_pueden_pedir_los_dos_caminos_a_la_vez(tmp_path):
    with pytest.raises(ErrorDeUso):
        B.prepare(_si(0.1), outdir=str(tmp_path), referencia=_si(0.0),
                  desplazar=(1, [0, 0, 0.1]),
                  pseudo_dir="/usr/share/espresso/pseudo")


def test_desplazar_un_atomo_que_no_existe(tmp_path):
    with pytest.raises(ErrorDeUso):
        B.prepare(_si(), outdir=str(tmp_path), desplazar=(9, [0, 0, 0.1]),
                  nppstr=5, kperp=(2, 2),
                  pseudo_dir="/usr/share/espresso/pseudo")


@pytest.mark.parametrize("g", [0, 4, -1])
def test_gdir_fuera_de_rango(tmp_path, g):
    with pytest.raises(ErrorDeUso):
        B.prepare(_si(), outdir=str(tmp_path), gdir=g, nppstr=5, kperp=(2, 2),
                  pseudo_dir="/usr/share/espresso/pseudo")


# ----------------------------------------------------------------------
# El contraste con los centros de Wannier
# ----------------------------------------------------------------------
def test_la_fase_electronica_de_los_centros_de_wannier_coincide_con_lberry():
    """Dos rutas que no comparten una línea de código.

    Los números son medidos: silicio con el segundo átomo desplazado 0.12 Å
    en z. `lberry` (36 cuerdas de 9 puntos) da fase electrónica 0.08840 y
    iónica 0.91160. Los cuatro centros de Wannier del mismo sistema, sacados
    de una malla 6×6×6 con proyección + minimización, están en las
    fraccionarias de abajo. Que las dos den lo mismo a 10⁻³ es la validación
    del método, no una casualidad: son la MISMA fase de Berry.
    """
    at = _si(0.12)
    frac = np.array([[0.136138, 0.136138, 0.113261],
                     [-0.364424, 0.136738, 0.113843],
                     [0.136738, -0.364424, 0.113843],
                     [0.136138, 0.136138, -0.385538]])
    centros = frac @ at.cell.array
    el, ion, tot = B.desde_wannier(centros, at, {"Si": 4.0}, gdir=3)
    assert el == pytest.approx(0.08840, abs=2e-3)
    assert ion == pytest.approx(0.91160, abs=1e-4)
    assert abs(B.a_intervalo(tot - 1.0)) < 3e-3     # total ≡ medio cuanto


def test_los_centros_de_un_cristal_centrosimetrico_dan_fase_de_medio_cuanto():
    at = _si()
    cell = at.cell.array
    frac = np.array([[0.125, 0.125, 0.125], [-0.375, 0.125, 0.125],
                     [0.125, -0.375, 0.125], [0.125, 0.125, -0.375]])
    el, ion, tot = B.desde_wannier(frac @ cell, at, {"Si": 4.0}, gdir=3)
    assert el == pytest.approx(0.0, abs=1e-9)
    assert abs(abs(tot) - 1.0) < 1e-9


def test_el_espin_multiplica_la_fase_electronica():
    at = _si()
    c = np.array([[0.1, 0.1, 0.1]]) @ at.cell.array
    e1, _, _ = B.desde_wannier(c, at, {"Si": 4.0}, gdir=3, spin=1.0)
    e2, _, _ = B.desde_wannier(c, at, {"Si": 4.0}, gdir=3, spin=2.0)
    assert e2 == pytest.approx(2 * e1)
