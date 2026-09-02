"""La referencia navegable que se genera del propio código."""

from qekit.modules import docs


def test_extrae_todos_los_subcomandos():
    cs = docs.extraer()
    assert len(cs) > 50
    nombres = {c.nombre for c in cs}
    for esperado in ("gen", "bands", "strain", "adsorb", "eform", "echem"):
        assert esperado in nombres


def test_ningun_comando_se_queda_sin_grupo():
    """Si se añade un subcomando y se olvida agruparlo, esto lo pilla."""
    cs = docs.extraer()
    agrupados = {n for _, _, ns in docs.GRUPOS for n in ns}
    huerfanos = sorted(c.nombre for c in cs if c.nombre not in agrupados)
    assert not huerfanos, f"sin grupo en docs.GRUPOS: {huerfanos}"


def test_los_grupos_no_citan_comandos_que_no_existen():
    cs = {c.nombre for c in docs.extraer()}
    fantasmas = sorted(n for _, _, ns in docs.GRUPOS for n in ns if n not in cs)
    assert not fantasmas, f"en GRUPOS pero no en la CLI: {fantasmas}"


def test_cada_comando_trae_sus_banderas():
    c = {x.nombre: x for x in docs.extraer()}["strain"]
    nombres = " ".join(b.nombres for b in c.banderas)
    assert "--range" in nombres and "--mode" in nombres
    assert any(b.nombres == "file" for b in c.posicionales)


def test_la_ayuda_agrupada_no_borra_las_descripciones_cortas():
    """La ayuda principal puede ocultar la lista plana, pero docs la usa."""
    comandos = {c.nombre for c in docs.extraer()}
    ayudas = docs.ayudas_cortas()
    assert set(ayudas) == comandos
    assert all(ayudas.values())
    assert "Chern" in ayudas["topology"]


def test_la_fisica_sale_del_docstring_del_modulo():
    c = {x.nombre: x for x in docs.extraer()}["strain"]
    assert "deformación" in c.fisica.lower()


# ----------------------------------------------------------------------
# formato del docstring a HTML
# ----------------------------------------------------------------------
def test_un_parrafo_normal():
    h = docs._parrafos("Una frase\ny su continuación.")
    assert h == "<p>Una frase y su continuación.</p>"


def test_las_vinetas_salen_como_lista():
    """El fallo original: las continuaciones sangradas salían como código."""
    txt = ("Lo que sale de aquí:\n"
           "  - primero, que sigue\n"
           "    en la línea de abajo.\n"
           "  - segundo.\n")
    h = docs._parrafos(txt)
    assert "<ul>" in h and h.count("<li>") == 2
    assert "<pre>" not in h
    assert "que sigue en la línea de abajo" in h


def test_un_bloque_sangrado_si_es_codigo():
    txt = "La fórmula es:\n\n    E = m c^2\n    con m en kg\n"
    h = docs._parrafos(txt)
    assert "<pre>" in h and "m c^2" in h


def test_el_texto_de_entrada_de_una_lista_no_se_come():
    h = docs._parrafos("Tres cosas:\n  - una\n  - dos\n")
    assert "<p>Tres cosas:</p>" in h


def test_las_comillas_invertidas_son_codigo():
    assert "<code>olla-dft gen</code>" in docs._parrafos("Usa `olla-dft gen` así.")


def test_los_asteriscos_dobles_son_negrita():
    assert "<b>importante</b>" in docs._parrafos("Esto es **importante**.")


def test_se_escapa_el_html():
    assert "&lt;script&gt;" in docs._parrafos("cuidado con <script>")


def test_docstring_vacio():
    assert docs._parrafos("") == ""
    assert docs._parrafos(None) == ""


# ----------------------------------------------------------------------
# la página
# ----------------------------------------------------------------------
def test_genera_una_pagina_autocontenida(tmp_path):
    f = docs.generar(str(tmp_path / "d.html"))
    txt = open(f, encoding="utf-8").read()
    assert txt.startswith("<!doctype html>")
    assert "http://" not in txt and "https://" not in txt, \
        "tiene que abrirse sin conexión"
    assert "olla-dft strain" in txt
    assert 'id="q"' in txt, "el buscador"


def test_la_pagina_lleva_la_version(tmp_path):
    from qekit import __version__
    f = docs.generar(str(tmp_path / "d.html"))
    assert __version__ in open(f, encoding="utf-8").read()


def test_cada_comando_tiene_su_ancla(tmp_path):
    f = docs.generar(str(tmp_path / "d.html"))
    txt = open(f, encoding="utf-8").read()
    for c in docs.extraer():
        assert f'id="c-{c.nombre}"' in txt, c.nombre


def test_genera_referencia_ingles_con_los_mismos_comandos(tmp_path):
    es = docs.generar(str(tmp_path / "docs.html"), language="es")
    en = docs.generar(str(tmp_path / "docs.en.html"), language="en")
    spanish = open(es, encoding="utf-8").read()
    english = open(en, encoding="utf-8").read()
    assert '<html lang="es">' in spanish
    assert '<html lang="en">' in english
    assert "Generate calculations" in english
    assert "generate pw.x inputs and post-processing" in english
    assert 'role="search"' in spanish and 'aria-live="polite"' in spanish
    assert 'aria-label="Reference navigation"' in english
    assert "No commands match the search." in english
    for command in ("gen", "bands", "topology", "project", "phonons"):
        assert f'id="c-{command}"' in spanish
        assert f'id="c-{command}"' in english


def test_contrato_de_traduccion_de_referencia_es_en_es_completo():
    assert set(docs._labels("es")) == set(docs._labels("en"))
