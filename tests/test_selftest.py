"""La suite de validación física, en su parte que no necesita pw.x."""

import pytest

from qekit.core.errors import ErrorDeUso
from qekit.modules import selftest as st


def test_hay_pruebas_registradas():
    assert len(st.PRUEBAS) >= 10
    claves = [p.clave for p in st.PRUEBAS]
    assert len(claves) == len(set(claves)), "las claves tienen que ser únicas"


def test_toda_prueba_declara_su_fuente():
    """Un valor 'de la literatura' sin decir de dónde sale no vale nada."""
    for p in st.PRUEBAS:
        assert p.fuente and len(p.fuente) > 20, p.clave
        assert p.unidad is not None
        assert 0 <= p.tolerancia <= 1.0, p.clave


def test_las_rapidas_pasan_todas():
    """Es la prueba de las pruebas: si esto falla, algo de física cambió."""
    res = st.ejecutar(con_qe=False, verbose=False)
    malas = [r for r in res if not r.bien]
    detalle = "; ".join(
        f"{r.prueba.clave}: {r.error or f'{r.valor:.5g} vs {r.prueba.referencia:g}'}"
        for r in malas)
    assert not malas, detalle


def test_ninguna_rapida_necesita_pw_x():
    res = st.ejecutar(con_qe=False, verbose=False)
    assert all(not r.prueba.necesita_qe for r in res)


def test_full_y_mlip_son_opciones_independientes(monkeypatch):
    comunes = dict(titulo="t", magnitud="m", referencia=1.0, unidad="",
                   tolerancia=0.0, fuente="identidad sintética para el test",
                   fn=lambda _ctx: 1.0)
    pruebas = [
        st.Prueba(clave="quick", **comunes),
        st.Prueba(clave="qe", necesita_qe=True, **comunes),
        st.Prueba(clave="mlip", necesita_mlip=True, **comunes),
    ]
    monkeypatch.setattr(st, "PRUEBAS", pruebas)
    solo_qe = st.ejecutar(con_qe=True, verbose=False)
    assert [r.prueba.clave for r in solo_qe] == ["quick", "qe"]
    con_ambas = st.ejecutar(con_qe=True, con_mlip=True, verbose=False)
    assert [r.prueba.clave for r in con_ambas] == ["quick", "qe", "mlip"]


def test_se_puede_pedir_una_sola():
    res = st.ejecutar(claves=["madelung"], verbose=False)
    assert len(res) == 1 and res[0].prueba.clave == "madelung"


def test_una_clave_inventada_se_queja():
    with pytest.raises(ErrorDeUso):
        st.ejecutar(claves=["no_existe"], verbose=False)


def test_la_referencia_cero_usa_desviacion_absoluta():
    """Con referencia 0 no hay desviacion relativa: dividir daba un crash."""
    p = st.Prueba(clave="x", titulo="t", magnitud="m", referencia=0.0,
                  unidad="", tolerancia=1e-6, fuente="identidad exacta")
    r = st.Resultado(prueba=p, valor=1e-9)
    assert r.relativa is False
    assert r.desviacion == pytest.approx(1e-9)
    assert r.bien


def test_una_prueba_con_error_no_cuenta_como_buena():
    p = st.Prueba(clave="x", titulo="t", magnitud="m", referencia=1.0,
                  unidad="", tolerancia=0.1, fuente="lo que sea")
    r = st.Resultado(prueba=p, valor=1.0, error="revento")
    assert not r.bien


def test_el_reporte_lleva_las_fuentes():
    res = st.ejecutar(claves=["madelung", "lorenz"], verbose=False)
    txt = st.report(res)
    assert "fuente:" in txt
    assert "Sommerfeld" in txt


def test_el_reporte_cuenta_bien():
    res = st.ejecutar(con_qe=False, verbose=False)
    txt = st.report(res)
    assert f"{len(res)} pruebas" in txt
