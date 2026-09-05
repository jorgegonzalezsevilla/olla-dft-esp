"""El fundamento físico escrito tiene que cubrir los comandos y estar al día."""

from pathlib import Path

import pytest

from qekit.cli import build_parser, main
from qekit.core.i18n import DEFAULT_LANGUAGE
from qekit.modules import theory

RAIZ = Path(__file__).resolve().parent.parent

# Comandos de infraestructura: no tienen física propia que explicar.
SIN_FISICA = {
    "start", "wizard", "recetas", "teoria", "docs", "sistema", "templates",
    "config", "project", "results", "campaign", "compare", "report",
    "datasheet", "convert", "update", "resilient",
}


def _comandos_cli():
    sub = build_parser("es")._subparsers._group_actions[0]
    return set(sub.choices) - set(("system", "recipes", "theory", "actualizar"))


def test_todos_los_comandos_cientificos_tienen_fundamento():
    documentados = {c for s in theory.secciones("es") for c in s.comandos}
    faltan = sorted(_comandos_cli() - documentados - SIN_FISICA)
    assert not faltan, f"comandos sin sección de teoría: {faltan}"


def test_las_secciones_no_citan_comandos_inexistentes():
    reales = _comandos_cli()
    for lang in ("es", "en"):
        for sec in theory.secciones(lang):
            for c in sec.comandos:
                assert c in reales, f"[{lang}] la teoría documenta '{c}', que no existe"


def test_espanol_e_ingles_en_paridad():
    es = [s.comandos for s in theory.secciones("es")]
    en = [s.comandos for s in theory.secciones("en")]
    assert es == en


def test_cada_seccion_trae_los_apartados_obligatorios():
    # Fórmulas y la tabla de datos pueden faltar en los comandos de control
    # de calidad, que listan reglas en vez de ecuaciones; el resto, siempre.
    apartados = {
        "es": ["Qué responde", "Fundamento para no expertos",
               "Cómo lo calcula Olla-DFT", "Límites y trampas", "Referencias"],
        "en": ["What it answers", "Background for non-experts",
               "How Olla-DFT computes it", "Limits and pitfalls", "References"],
    }
    for lang, nombres in apartados.items():
        for sec in theory.secciones(lang):
            for nombre in nombres:
                assert f"**{nombre}" in sec.texto, \
                    f"[{lang}] {sec.comandos}: falta el apartado '{nombre}'"


def test_las_formulas_estan_balanceadas():
    for lang in ("es", "en"):
        for sec in theory.secciones(lang):
            assert sec.texto.count("$$") % 2 == 0, f"[{lang}] {sec.comandos}"


def test_los_documentos_publicados_estan_al_dia():
    nombre = {"es": "TEORIA.md", "en": "THEORY.md"}[DEFAULT_LANGUAGE]
    publicado = RAIZ / "docs" / nombre
    if not publicado.exists():
        pytest.skip("docs/ no está en el árbol instalado")
    assert publicado.read_text(encoding="utf-8") == theory.documento(DEFAULT_LANGUAGE), \
        f"docs/{nombre} está desactualizado: corre python tools/build_docs.py"


def test_el_comando_imprime_una_seccion_y_el_indice(capsys):
    assert main(["teoria", "eos"]) == 0
    out = capsys.readouterr().out
    assert "Birch" in out and "**" not in out
    assert main(["theory", "--language", "en"]) == 0
    assert "Electronic structure" in capsys.readouterr().out


def test_un_comando_sin_fundamento_da_error_de_uso(capsys):
    assert main(["teoria", "inventado"]) == 2
    assert "no hay fundamento" in capsys.readouterr().err


def test_la_referencia_de_comandos_publicada_esta_al_dia():
    import importlib.util
    guion = RAIZ / "tools" / "build_docs.py"
    nombre = {"es": "COMANDOS.md", "en": "COMMANDS.md"}[DEFAULT_LANGUAGE]
    if not guion.exists() or not (RAIZ / "docs" / nombre).exists():
        pytest.skip("tools/ o docs/ no están en el árbol instalado")
    spec = importlib.util.spec_from_file_location("build_docs", guion)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    publicado = (RAIZ / "docs" / nombre).read_text(encoding="utf-8")
    assert publicado == mod.comandos_md(DEFAULT_LANGUAGE), \
        f"docs/{nombre} está desactualizado: corre python tools/build_docs.py"
