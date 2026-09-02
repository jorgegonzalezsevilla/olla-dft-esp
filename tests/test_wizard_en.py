"""El asistente en inglés: el mismo catálogo, traducido al vuelo.

No hay una copia inglesa de METAS ni de GLOSARIO: hay una tabla
``qekit/data/i18n/wizard_en.json`` y ``i18n.translate_data`` recorre la
estructura original. La prueba que mantiene la tabla al día es
``test_la_tabla_inglesa_cubre_todo_el_catalogo``: añadir una meta o un
término sin su traducción falla aquí antes de que lo vea nadie.
"""

import re

import pytest

from qekit.core import i18n
from qekit.core.errors import ErrorDeUso
from qekit.modules import wizard as wz


def _tabla():
    return i18n.load_table("wizard_en")


def _cadenas_traducibles(m):
    """Lo que un usuario lee de una meta y no es comando ni clave."""
    fuera = [m.pregunta, m.nombre, m.explica, m.coste, m.error_tipico]
    fuera += [desc for desc, _cmd in m.pasos]
    return [s for s in fuera if s]


# ----------------------------------------------------------------------
# translate_data, el mecanismo
# ----------------------------------------------------------------------
def test_translate_data_recorre_dataclasses_listas_tuplas_y_dicts():
    tabla = {"hola": "hello", "adiós": "bye"}
    m = wz.Meta("k", "hola", "adiós", "hola", pasos=[("hola", "cmd hola")],
                necesita=["hola"], terminos=["adiós"])
    t = i18n.translate_data(m, tabla)
    assert t.pregunta == "hello" and t.nombre == "bye" and t.clave == "k"
    assert t.pasos == [("hello", "cmd hola")]        # tupla sigue siendo tupla
    assert isinstance(t.pasos[0], tuple)
    assert t.necesita == ["hello"] and t.terminos == ["bye"]
    assert i18n.translate_data({"hola": ["adiós", 3, None]}, tabla) \
        == {"hello": ["bye", 3, None]}
    assert m.pregunta == "hola"                       # el original no se toca
    assert i18n.load_table("no_existe_esta_tabla") == {}


# ----------------------------------------------------------------------
# La tabla
# ----------------------------------------------------------------------
def test_la_tabla_inglesa_cubre_todo_el_catalogo():
    strings = _tabla()["strings"]
    faltan = [f"{m.clave}: {s[:60]!r}" for m in wz.METAS
              for s in _cadenas_traducibles(m) if s not in strings]
    faltan += [f"glosario: {s[:60]!r}" for k, v in wz.GLOSARIO.items()
               for s in (k, v) if s not in strings]
    assert not faltan, "faltan en wizard_en.json:\n  " + "\n  ".join(faltan)
    assert all(v.strip() for v in strings.values())


def test_la_tabla_tiene_palabras_de_busqueda_para_cada_meta():
    kw = _tabla()["keywords"]
    assert set(kw) == {m.clave for m in wz.METAS}
    assert all(kw[m.clave] for m in wz.METAS)


def test_traducir_dos_veces_no_cambia_nada():
    strings = _tabla()["strings"]
    valores = set(strings.values())
    assert not [k for k in strings if k in valores and strings[k] != k]
    for m in wz.METAS:
        una = i18n.translate_data(m, strings)
        assert i18n.translate_data(una, strings) == una


# ----------------------------------------------------------------------
# Lo que se imprime
# ----------------------------------------------------------------------
_CASTELLANO = re.compile(
    r"\b(del|de la|cada|archivo|carpeta|paso|pasos|después|cuando|hasta)\b",
    re.IGNORECASE)


@pytest.mark.parametrize("clave", [m.clave for m in wz.METAS])
def test_el_informe_en_ingles_no_deja_castellano(clave):
    m = wz.metas_por_clave("en")[clave]
    txt = wz.report_meta(m, "x.cif", language="en")
    restos = _CASTELLANO.findall(txt)
    assert not restos, f"{clave}: {sorted(set(restos))}"
    assert "The question:" in txt and "Cost:" in txt and "Steps:" in txt
    if m.terminos:
        assert "Terms that appear here:" in txt
    if m.necesita:
        assert "Needed first:" in txt and "[from " in txt


def test_los_terminos_traducidos_estan_en_el_glosario_traducido():
    glos = wz.glosario("en")
    assert "k-mesh" in glos and "band gap" in glos and "malla k" not in glos
    for m in wz.metas("en"):
        for t in m.terminos:
            assert t in glos, f"{m.clave}: {t}"
    txt = wz.report_meta(wz.metas_por_clave("en")["conduce"], "x.cif",
                         language="en")
    assert "  band gap: The energy distance" in txt


def test_en_ingles_los_comandos_no_cambian():
    for es, en in zip(wz.METAS, wz.metas("en")):
        assert es.clave == en.clave and es.necesita == en.necesita
        assert [c for _d, c in es.pasos] == [c for _d, c in en.pasos]
    p_es = wz.plan("conduce", "mi.cif")
    p_en = wz.plan("conduce", "mi.cif", language="en")
    assert [c for _k, _d, c in p_es] == [c for _k, _d, c in p_en]
    assert [d for _k, d, _c in p_es] != [d for _k, d, _c in p_en]


def test_en_espanol_nada_cambia():
    assert wz.metas("es") is wz.METAS
    assert wz.metas() is wz.METAS
    assert wz.glosario("es") is wz.GLOSARIO
    m = wz.METAS_POR_CLAVE["conduce"]
    txt = wz.report_meta(m, "x.cif")
    assert txt == wz.report_meta(m, "x.cif", language="es")
    assert 'La pregunta: "' in txt and "Coste: medio" in txt
    assert "Antes hace falta: relajación de la estructura" in txt
    assert "Términos que salen aquí:" in txt and "  gap: " in txt
    assert wz.report_catalogo().startswith("--- ¿Qué quieres saber? ---")
    assert m.coste == "medio"


def test_el_catalogo_en_ingles():
    txt = wz.report_catalogo("en")
    assert txt.startswith("--- What do you want to know? ---")
    assert "  conduce " in txt and "--goal <key>" in txt
    assert not _CASTELLANO.findall(txt)


def test_la_meta_desconocida_se_explica_en_ingles():
    with pytest.raises(ErrorDeUso, match="unknown goal 'nada'. Available"):
        wz.plan("nada", "x.cif", language="en")
    with pytest.raises(ErrorDeUso, match="Disponibles"):
        wz.plan("nada", "x.cif")


def test_el_diagnostico_se_traduce_con_sus_numeros():
    from ase.build import bulk, graphene
    d = wz.diagnosticar(bulk("Ni", "fcc", a=3.52))
    en = wz.report_diagnostico(d, "en")
    assert en.startswith("--- What I see in your structure ---")
    assert "self-interaction" in en and "(Ni)" in en
    assert "autointeracción" in wz.report_diagnostico(d)
    d2 = wz.diagnosticar(graphene(a=2.46, vacuum=10.0))
    en2 = wz.report_diagnostico(d2, "en")
    assert "of vacuum along the c axis" in en2
    assert re.search(r"\d+\.\d Å", en2)
    # un Diagnostico armado a mano solo con `notas` sale tal cual
    d3 = wz.Diagnostico(formula="X", natoms=1, grupo="no determinado",
                        elementos=["X"], notas=["una nota"])
    txt = wz.report_diagnostico(d3, "en")
    assert "una nota" in txt and "not determined" in txt


# ----------------------------------------------------------------------
# Búsqueda en los dos idiomas
# ----------------------------------------------------------------------
@pytest.mark.parametrize("texto,clave", [
    ("does it absorb light", "color"),
    ("what colour will it be", "color"),
    ("is it a metal", "conduce"),
    ("band gap", "conduce"),
    ("is it hard or brittle", "mecanicas"),
    ("x-ray diffraction pattern", "difractograma"),
    ("I want to put one material on top of another", "interfase"),
    ("my oxide comes out metallic", "oxido"),
    ("does it conduct heat", "conduce_calor"),
    ("how do I know my results are right", "fiarme"),
])
def test_la_busqueda_tambien_entiende_ingles(texto, clave):
    cands = wz.buscar(texto)
    assert cands, texto
    assert clave in [m.clave for m in cands[:3]], [m.clave for m in cands]


def test_la_busqueda_devuelve_las_metas_en_el_idioma_pedido():
    es = wz.buscar("absorbe luz")[0]
    en = wz.buscar("absorbe luz", language="en")[0]
    assert es.clave == en.clave == "color"
    assert es.pregunta.startswith("¿de qué color")
    assert en.pregunta.startswith("what colour")


def test_la_busqueda_en_espanol_sigue_igual():
    """Las palabras inglesas no desordenan lo que ya funcionaba."""
    for texto, clave in [("quiero saber si absorbe luz visible", "color"),
                         ("mi oxido sale metalico y no deberia", "oxido"),
                         ("es duro o fragil", "mecanicas"),
                         ("como se ve en rayos x", "difractograma"),
                         ("quiero poner un material sobre otro", "interfase")]:
        assert wz.buscar(texto)[0].clave == clave, texto


# ----------------------------------------------------------------------
# Desde el CLI
# ----------------------------------------------------------------------
def test_el_comando_wizard_habla_ingles_con_language_en(capsys):
    from qekit.cli import main
    assert main(["wizard", "--list", "--language", "en"]) == 0
    assert "What do you want to know?" in capsys.readouterr().out
    assert main(["wizard", "--ask", "does it absorb light",
                 "--language", "en"]) == 0
    out = capsys.readouterr().out
    assert "--- optical properties ---" in out and "norm-conserving" in out
    assert main(["wizard", "--term", "k-mesh", "--language", "en"]) == 0
    assert "k-mesh: How many points" in capsys.readouterr().out
    assert main(["wizard", "--goal", "nada", "--language", "en"]) == 2
    assert "unknown goal 'nada'" in capsys.readouterr().err
    assert main(["wizard", "--goal", "color", "--language", "es"]) == 0
    assert "--- propiedades ópticas ---" in capsys.readouterr().out
