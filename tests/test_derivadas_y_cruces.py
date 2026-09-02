"""Pruebas de las derivadas termoelásticas, la validación cruzada, la QHA
y la ficha del material.

Casi todo se comprueba contra valores experimentales de silicio o contra
respuestas exactas de casos sintéticos: son los únicos contrastes que
sirven para algo en un módulo que solo hace aritmética sobre resultados de
otros.
"""

from pathlib import Path

import numpy as np
import pytest

# Constantes elásticas LDA del silicio validadas en este proyecto
SI_C11, SI_C12, SI_C44 = 159.9, 61.7, 76.6
SI_VOL, SI_NAT = 40.05, 2
SI_MASAS = [28.0855, 28.0855]


def _C_cubica(c11=SI_C11, c12=SI_C12, c44=SI_C44):
    C = np.zeros((6, 6))
    for i in range(3):
        C[i, i] = c11
        for j in range(3):
            if i != j:
                C[i, j] = c12
    for i in range(3, 6):
        C[i, i] = c44
    return C


# ----------------------------------------------------------------------
# Derivadas termoelásticas
# ----------------------------------------------------------------------
def test_densidad_del_silicio():
    from qekit.modules import derived
    rho = derived.density(SI_MASAS, SI_VOL)
    assert rho == pytest.approx(2329.0, rel=0.01)      # exp. 2329 kg/m³


def test_velocidades_direccionales_contra_experimento():
    """v_L = √(C₁₁/ρ) y v_T = √(C₄₄/ρ) a lo largo de [100]."""
    from qekit.modules import derived
    rho = derived.density(SI_MASAS, SI_VOL)
    d = derived.cubic_directional(_C_cubica(), rho)
    assert d["v_l_100"] == pytest.approx(8433, rel=0.04)   # exp. 8433 m/s
    assert d["v_t_100"] == pytest.approx(5843, rel=0.04)   # exp. 5843 m/s


def test_temperatura_de_debye_del_silicio():
    from qekit.modules import derived, elastic
    m = elastic.moduli(_C_cubica())
    r = derived.analyze(m.B_hill, m.G_hill, SI_MASAS, SI_VOL, natoms=SI_NAT)
    assert r.theta_D == pytest.approx(645, rel=0.05)   # exp. 645 K
    assert r.poisson == pytest.approx(0.22, abs=0.02)  # exp. 0.22


def test_slack_da_un_orden_de_magnitud_razonable():
    """El prefactor está ajustado para delta en Å; en metros da cero."""
    from qekit.modules import derived, elastic
    m = elastic.moduli(_C_cubica())
    r = derived.analyze(m.B_hill, m.G_hill, SI_MASAS, SI_VOL, natoms=SI_NAT)
    # Si real: 148 W/(m·K). Slack subestima, pero no puede dar 0 ni 10^4.
    assert 20.0 < r.kappa_slack < 400.0


def test_conversion_de_pendiente_a_velocidad():
    """v[m/s] = 2*pi*c*(dnu/dq) = 18.836 * (dnu/dq) con nu en cm⁻¹, q en Å⁻¹.

    Un factor 100 de más aquí infla las velocidades del sonido a cientos
    de km/s sin que nada más falle.
    """
    from qekit.modules.derived import CM1_A_A_MS, acoustic_velocities
    assert CM1_A_A_MS == pytest.approx(18.8365, rel=1e-4)
    # rama perfectamente lineal de pendiente conocida
    q = np.linspace(0.0, 0.1, 6)
    pend = 448.0                       # cm⁻¹/Å⁻¹
    f = np.column_stack([q * pend] * 3)
    a = acoustic_velocities(q, f, n_puntos=5)
    assert a["v_l"] == pytest.approx(pend * 18.8365, rel=1e-6)


def test_debye_de_la_dos_usa_todo_el_espectro():
    """La Debye de la DOS es mayor que la elástica porque cuenta ópticas."""
    from qekit.modules import derived
    w = np.linspace(1, 500, 500)
    dos = np.ones_like(w)              # DOS plana
    td = derived.debye_from_dos(w, dos, natoms=2)
    assert td > 0 and np.isfinite(td)


# ----------------------------------------------------------------------
# Validación cruzada
# ----------------------------------------------------------------------
def test_cruce_detecta_acuerdo_y_desacuerdo():
    from qekit.modules import crosscheck as cc
    ok = cc.Check(nombre="x", valor_a=100.0, valor_b=102.0, tolerancia=0.05)
    mal = cc.Check(nombre="y", valor_a=100.0, valor_b=150.0, tolerancia=0.05)
    assert ok.ok is True and ok.desvio == pytest.approx(0.02)
    assert mal.ok is False and mal.desvio == pytest.approx(0.50)


def test_cruce_b0_eos_contra_elasticas():
    from qekit.modules import crosscheck as cc
    r = cc.run(project=None, C=_C_cubica(), b0_eos=94.2)
    b = [c for c in r.checks if "volumétrico" in c.nombre]
    assert len(b) == 1 and b[0].ok is True


def test_cruce_marca_la_transversal_con_malla_pobre():
    """Una rama TA con pendiente muy baja tiene que disparar el aviso."""
    from qekit.modules import crosscheck as cc
    q = np.linspace(0.0, 0.15, 6)
    # LA correcta (~440) y TA a la mitad de lo que debería
    f = np.column_stack([q * 170.0, q * 170.0, q * 440.0])
    r = cc.run(project=None, C=_C_cubica(), masas=SI_MASAS, volumen=SI_VOL,
               natoms=SI_NAT, qdist=q, band_freqs=f)
    tv = [c for c in r.checks if "transversal" in c.nombre][0]
    assert tv.ok is False
    assert "malla de q" in tv.diagnostico


def test_cruce_numero_de_modos_y_dulong_petit():
    from qekit.modules import crosscheck as cc
    # DOS plana normalizada a 3N = 6 modos entre 1 y 500 cm⁻¹
    w = np.linspace(1.0, 500.0, 400)
    dos = np.full_like(w, 6.0 / (500.0 - 1.0))
    r = cc.run(project=None, dos_w=w, dos=dos, natoms=2)
    modos = [c for c in r.checks if "modos" in c.nombre][0]
    assert modos.ok is True
    dp = [c for c in r.checks if "Dulong" in c.ruta_a][0]
    assert dp.ok is True          # a 1500 K con w<500 cm⁻¹ ya es clásico


def test_cruce_sin_datos_no_inventa():
    from qekit.modules import crosscheck as cc
    r = cc.run(project=None)
    assert r.checks == []
    assert "dos rutas independientes" in cc.report(r)


# ----------------------------------------------------------------------
# Cuasi-armónica
# ----------------------------------------------------------------------
def _modos_con_gruneisen(V, w0=300.0, gamma=1.5, V0=40.0, n=6):
    """Modos cuya frecuencia sigue w = w0 (V0/V)^gamma exactamente."""
    return np.full(n, w0 * (V0 / V) ** gamma)


def test_qha_recupera_el_gruneisen_impuesto():
    """Si se construyen modos con gamma conocido, la QHA debe medirlo."""
    from qekit.modules import qha
    V = np.linspace(36.0, 44.0, 7)
    E = 0.02 * (V - 40.0) ** 2                     # parábola en el mínimo
    F = [_modos_con_gruneisen(v, gamma=1.5) for v in V]
    r = qha.run(V, E, F, T=np.arange(0, 601, 20), natoms=2)
    assert r.gruneisen == pytest.approx(1.5, rel=0.02)


def test_qha_da_expansion_positiva_con_gruneisen_positivo():
    from qekit.modules import qha
    V = np.linspace(36.0, 44.0, 9)
    E = 0.02 * (V - 40.0) ** 2
    F = [_modos_con_gruneisen(v, gamma=2.0) for v in V]
    r = qha.run(V, E, F, T=np.arange(0, 601, 20), natoms=2)
    i = int(np.argmin(np.abs(r.T - 300.0)))
    assert r.alpha[i] > 0
    assert r.V_T[i] > r.V_T[0]           # se dilata al calentar


def test_qha_da_expansion_negativa_con_gruneisen_negativo():
    """El caso del silicio a baja T: gamma negativo => contrae al calentar.

    Si la implementación no reprodujera este signo, estaría mal.
    """
    from qekit.modules import qha
    V = np.linspace(36.0, 44.0, 9)
    E = 0.02 * (V - 40.0) ** 2
    F = [_modos_con_gruneisen(v, gamma=-1.0) for v in V]
    r = qha.run(V, E, F, T=np.arange(0, 601, 20), natoms=2)
    i = int(np.argmin(np.abs(r.T - 300.0)))
    assert r.alpha[i] < 0


def test_qha_cv_tiende_a_dulong_petit():
    from qekit.modules import qha
    cv = qha.cv_modos(np.full(6, 200.0), T=3000.0)
    assert cv == pytest.approx(3 * 2 * 8.617333262e-5 * 1000, rel=0.02)


def test_energia_libre_a_cero_kelvin_es_el_punto_cero():
    from qekit.modules import qha
    w = np.array([100.0, 200.0, 300.0])
    zpe = 0.5 * np.sum(w * 1.239841984e-4)
    assert qha.f_vib(w, 0.0) == pytest.approx(zpe, rel=1e-12)
    # y a T alta la energía libre baja
    assert qha.f_vib(w, 800.0) < zpe


def test_qha_avisa_con_pocos_volumenes():
    from qekit.modules import qha
    V = np.array([39.0, 40.0, 41.0])
    F = [_modos_con_gruneisen(v) for v in V]
    r = qha.run(V, 0.02 * (V - 40) ** 2, F, T=np.array([300.0]), natoms=2)
    assert any("volúmenes" in a for a in r.avisos)


# ----------------------------------------------------------------------
# Ficha del material
# ----------------------------------------------------------------------
def test_ficha_recoge_lo_que_encuentra(tmp_path):
    from qekit.modules import datasheet as ds
    np.savetxt(tmp_path / "ELASTIC_C.dat", _C_cubica(), header="C (GPa)")
    (tmp_path / "EOS.txt").write_text("a0 = 5.402 A\nB0 = 94.2 GPa\n")
    f = ds.recoger(tmp_path)
    assert "Elásticas" in f.resultados
    assert "Ecuación de estado" in f.resultados
    magnitudes = [r["magnitud"] for r in f.resultados["Elásticas"]]
    assert any("C₁₁" in m for m in magnitudes)


def test_metodos_usa_los_parametros_reales():
    from qekit.modules import datasheet as ds
    f = ds.Ficha(formula="Si2", parametros={
        "funcional": "PBE", "ecutwfc_Ry": 60.0, "ecutrho_Ry": 480.0,
        "malla_k": "8x8x8", "pseudos": {"Si": "Si.upf"},
        "ocupaciones": "fixed", "nspin": 1})
    texto = ds.metodos(f)
    assert "PBE" in texto and "60.0" in texto and "8x8x8" in texto
    assert "Si.upf" in texto and "ocupaciones fijas" in texto


def test_metodos_declara_el_uso_de_mlip():
    from qekit.modules import datasheet as ds
    f = ds.Ficha(parametros={"funcional": "PBE"}, codigos=["qe", "mace"])
    texto = ds.metodos(f)
    assert "MACE" in texto
    assert "resultados reportados provienen" in texto   # deja claro el rol


def test_ficha_escribe_markdown_y_html(tmp_path):
    from qekit.modules import datasheet as ds
    f = ds.Ficha(formula="Si2",
                 resultados={"Elásticas": [
                     {"magnitud": "C₁₁", "valor": 159.9, "unidad": "GPa",
                      "nota": ""}]},
                 parametros={"funcional": "PZ"})
    archivos = ds.escribir(f, tmp_path)
    assert len(archivos) == 2
    html = Path(archivos[1]).read_text(encoding="utf-8")
    assert "<!doctype html>" in html.lower()
    assert "159.9" in html
    assert "borrador" in html          # el aviso sobre el párrafo de métodos
