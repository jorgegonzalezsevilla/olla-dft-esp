"""Recetas: sesiones completas que enseñan cómo encajan los módulos.

La prueba que justifica todo el archivo es
`test_todos_los_comandos_de_las_recetas_existen_de_verdad`: recorre cada
comando de cada receta y lo valida contra el PROPIO árbol de argparse. Una
sección de ejemplos escrita a mano se queda obsoleta a la tercera versión y
nadie se entera; así, renombrar una bandera rompe la prueba antes de que
llegue a un usuario.
"""

import argparse
import re
import shlex

import pytest

from qekit.core.errors import ErrorDeUso
from qekit.modules import recipes as R


def _subparsers(parser):
    for acc in parser._actions:
        if isinstance(acc, argparse._SubParsersAction):
            return acc.choices
    return {}


def _opciones(p):
    fuera = set()
    for acc in p._actions:
        fuera.update(acc.option_strings)
    return fuera


def _comandos_de_receta(r):
    """Los comandos de la receta que invocan a olla-dft."""
    for i, paso in enumerate(r.pasos):
        if paso.comando.startswith("olla-dft "):
            yield i, paso.comando
    for _texto, cmd in r.despues:
        if cmd.startswith("olla-dft "):
            yield None, cmd


def test_hay_recetas_y_todas_tienen_lo_minimo():
    assert len(R.RECETAS) >= 8
    for r in R.RECETAS:
        assert r.clave and r.titulo and r.pregunta and r.para_que
        assert r.pasos, r.clave
        for p in r.pasos:
            assert p.comando and p.hace, f"{r.clave}: {p.comando}"


def test_las_claves_no_se_repiten():
    claves = [r.clave for r in R.RECETAS]
    assert len(claves) == len(set(claves))


def test_todos_los_comandos_de_las_recetas_existen_de_verdad():
    """Contra el árbol de argparse, no contra una lista escrita aparte."""
    from qekit.cli import build_parser
    subs = _subparsers(build_parser())
    malos = []
    for r in R.RECETAS:
        for i, cmd in _comandos_de_receta(r):
            trozos = shlex.split(cmd)[1:]
            if not trozos:
                continue
            nombre = trozos[0]
            if nombre not in subs:
                malos.append(f"{r.clave}: subcomando '{nombre}' no existe")
                continue
            validas = _opciones(subs[nombre])
            for t in trozos[1:]:
                if not t.startswith("-") or t == "-":
                    continue
                base = t.split("=")[0]
                if base in validas:
                    continue
                # puede ser un valor negativo de la opción anterior
                if len(base) > 1 and base[1].isdigit():
                    continue
                malos.append(f"{r.clave}: 'olla-dft {nombre}' no acepta {base}")
    assert not malos, "\n".join(malos)


def test_los_pasos_no_leen_de_pasos_que_no_existen():
    for r in R.RECETAS:
        etiquetas = set(R.CIRCULOS[:len(r.pasos)])
        for p in r.pasos:
            for _texto, de in p.lee:
                # o cita un paso de esta receta, o dice de dónde viene
                cita_paso = any(c in de for c in etiquetas)
                assert cita_paso or "receta" in de or "tu" in de or "todos" in de \
                    or "cualquier" in de, f"{r.clave}: referencia '{de}'"


def test_las_recetas_relacionadas_existen():
    for r in R.RECETAS:
        for c in r.ver_tambien:
            assert c in R.RECETAS_POR_CLAVE, f"{r.clave} -> {c}"


def test_una_receta_verificada_lleva_salidas_reales():
    """Si dice ✓, tiene que enseñar lo que sale de verdad en algún paso."""
    for r in R.RECETAS:
        if r.verificada:
            assert any(p.salida for p in r.pasos), r.clave


def test_obtener_una_receta_que_no_existe_es_error_de_uso():
    with pytest.raises(ErrorDeUso) as e:
        R.obtener("no_existe")
    assert "bandas" in str(e.value)          # lista las que sí hay


@pytest.mark.parametrize("texto,esperada", [
    ("quiero saber si conduce el calor", "termoelectrico"),
    ("acabo de instalarlo", "primero"),
    ("cuanto vale su gap", "bandas"),
    ("como se que no son basura", "fiarme"),
])
def test_la_busqueda_en_lenguaje_llano_encuentra_la_receta(texto, esperada):
    res = R.buscar(texto)
    assert res, texto
    assert esperada in [x.clave for x in res[:3]]


def test_el_informe_ensena_el_flujo_de_archivos():
    """Lo que distingue una receta de una lista de comandos."""
    txt = R.report(R.obtener("mecanicas"))
    assert "→ deja" in txt and "← lee" in txt
    assert "ELASTIC_C.dat" in txt
    assert "OJO:" in txt


def test_el_informe_avisa_si_la_receta_no_se_ha_corrido_entera():
    txt = R.report(R.obtener("termoelectrico"))
    assert "no se ha corrido entera" in txt


def test_el_guion_es_shell_valido_y_lleva_los_comandos(tmp_path):
    r = R.obtener("bandas")
    f = tmp_path / "s.sh"
    txt = R.script(r, str(f))
    assert txt.startswith("#!/bin/bash")
    assert f.exists()
    for p in r.pasos:
        assert p.comando in txt
    # y los avisos van comentados, no sueltos
    for linea in txt.split("\n"):
        if linea.strip() and not linea.startswith("#"):
            assert linea.startswith("olla-dft ") or linea.startswith("bash ") \
                or linea in ("set -e", "#!/bin/bash"), linea


def test_el_guion_marca_los_pasos_que_tardan():
    txt = R.script(R.obtener("mecanicas"))
    assert "CORRE pw.x" in txt


def test_cada_receta_cita_al_menos_un_archivo_que_produce():
    for r in R.RECETAS:
        escritos = [f for p in r.pasos for f in p.escribe
                    if "nada" not in f]
        assert escritos, r.clave


# ======================================================================
# En inglés: la misma estructura, traducida al vuelo por una tabla
# ======================================================================
def _cadenas_traducibles(r):
    """Todo lo que un usuario lee de una receta y NO es comando, archivo
    ni salida real. Es la definición de "lo que tiene que estar en
    recipes_en.json"; si se añade un campo de texto a Receta o Paso, hay
    que añadirlo aquí."""
    fuera = [r.titulo, r.pregunta, r.para_que, r.coste]
    for p in r.pasos:
        fuera.append(p.hace)
        fuera += [f for f in p.escribe if " " in f]     # "nada, solo imprime"
        for texto, de in p.lee:
            fuera.append(texto)
            if re.search(r"[A-Za-zÀ-ÿ]", de):         # no es solo un ①
                fuera.append(de)
        if p.ojo:
            fuera.append(p.ojo)
    fuera += [texto for texto, _cmd in r.despues]
    return [s for s in fuera if s]


def _tabla_en():
    from qekit.core import i18n
    return i18n.load_table("recipes_en")


def test_la_tabla_inglesa_cubre_todas_las_cadenas_de_todas_las_recetas():
    """Lo que mantiene recipes_en.json al día: añadir una receta sin su
    traducción falla aquí, no delante de un usuario."""
    strings = _tabla_en()["strings"]
    faltan = [f"{r.clave}: {s[:60]!r}" for r in R.RECETAS
              for s in _cadenas_traducibles(r) if s not in strings]
    assert not faltan, "faltan en recipes_en.json:\n  " + "\n  ".join(faltan)
    assert all(v.strip() for v in strings.values())


def test_la_tabla_inglesa_tiene_palabras_de_busqueda_para_cada_receta():
    kw = _tabla_en()["keywords"]
    assert set(kw) == {r.clave for r in R.RECETAS}
    assert all(kw[r.clave] for r in R.RECETAS)


def test_traducir_dos_veces_no_cambia_nada():
    """report() traduce lo que recibe, venga en español o ya en inglés."""
    from qekit.core import i18n
    strings = _tabla_en()["strings"]
    valores = set(strings.values())
    assert not [k for k in strings if k in valores and strings[k] != k]
    for r in R.RECETAS:
        una = i18n.translate_data(r, strings)
        assert i18n.translate_data(una, strings) == una


_CASTELLANO = re.compile(
    r"\b(del|de la|cada|archivo|carpeta|paso|pasos|después|cuando|hasta)\b",
    re.IGNORECASE)


def _sin_extractos(txt):
    """El informe sin las líneas de salida real, que siguen en español."""
    return "\n".join(l for l in txt.split("\n") if "│" not in l)


@pytest.mark.parametrize("clave", [r.clave for r in R.RECETAS])
def test_el_informe_en_ingles_no_deja_castellano(clave):
    txt = R.report(R.obtener(clave, "en"), "en")
    restos = _CASTELLANO.findall(_sin_extractos(txt))
    assert not restos, f"{clave}: {sorted(set(restos))}"
    assert "The question:" in txt and "Cost:" in txt
    assert "→ writes" in txt and "WATCH OUT:" in txt


def test_en_ingles_los_comandos_y_la_salida_real_quedan_igual():
    es = R.obtener("bandas")
    en = R.obtener("bandas", "en")
    assert [p.comando for p in es.pasos] == [p.comando for p in en.pasos]
    assert [p.salida for p in es.pasos] == [p.salida for p in en.pasos]
    assert [c for _t, c in es.despues] == [c for _t, c in en.despues]
    assert es.clave == en.clave and es.ver_tambien == en.ver_tambien
    txt = R.report(en, "en")
    # el extracto sigue en español, pero avisado
    assert "Gap indirecto: 0.524 eV" in txt
    assert "(real output of the run, in Spanish)" in txt
    assert R.report(es).count("│") == txt.count("│")


def test_en_espanol_nada_cambia():
    """La traducción no toca los datos: en español se devuelve el original."""
    assert R.recetas("es") is R.RECETAS
    assert R.obtener("mecanicas") is R.RECETAS_POR_CLAVE["mecanicas"]
    assert R.obtener("mecanicas", "es") is R.obtener("mecanicas")
    txt = R.report(R.obtener("mecanicas"))
    assert txt == R.report(R.obtener("mecanicas"), "es")
    assert "La pregunta: «" in txt and "Coste:" in txt and "OJO:" in txt
    assert "→ deja" in txt and "← lee" in txt and "Y después:" in txt
    assert "(real output" not in txt
    assert R.listar().startswith("--- Recetas: sesiones completas")
    assert "Guionizar:" in R.listar()


def test_el_listado_y_el_guion_en_ingles(tmp_path):
    lst = R.listar("en")
    assert lst.startswith("--- Recipes: complete sessions")
    assert "steps ·" in lst and "As a script:" in lst
    assert not _CASTELLANO.findall(lst)
    r = R.obtener("mecanicas", "en")
    f = tmp_path / "s.sh"
    txt = R.script(r, str(f), "en")
    assert f.exists() and txt.startswith("#!/bin/bash")
    assert "THIS STEP RUNS pw.x" in txt and "# --- step 1:" in txt
    assert not _CASTELLANO.findall(txt)
    for p in R.obtener("mecanicas").pasos:
        assert p.comando in txt


def test_una_receta_desconocida_se_explica_en_ingles():
    with pytest.raises(ErrorDeUso) as e:
        R.obtener("no_existe", "en")
    assert "there is no recipe called 'no_existe'" in str(e.value)
    assert "bandas" in str(e.value)


@pytest.mark.parametrize("texto,esperada", [
    ("band gap", "bandas"),
    ("is it a metal or a semiconductor", "bandas"),
    ("does it conduct heat", "termoelectrico"),
    ("how do I know my results are not garbage", "fiarme"),
    ("I just installed it, where do I start", "primero"),
    ("how hard is it", "mecanicas"),
])
def test_la_busqueda_tambien_entiende_ingles(texto, esperada):
    res = R.buscar(texto)
    assert res, texto
    assert esperada in [x.clave for x in res[:3]], [x.clave for x in res]
    # y devuelve la receta en el idioma pedido
    assert R.buscar(texto, language="en")[0].pregunta \
        == R.obtener(res[0].clave, "en").pregunta


def test_el_comando_recetas_habla_ingles_con_language_en(capsys):
    from qekit.cli import main
    assert main(["recetas", "primero", "--language", "en"]) == 0
    out = capsys.readouterr().out
    assert "The first session" in out and "WATCH OUT:" in out
    assert main(["recipes", "--buscar", "band gap", "--language", "en"]) == 0
    assert "Band structure, DOS and band gap" in capsys.readouterr().out
    # Sin bandera manda el idioma por defecto del repositorio; en español
    # se pide explícito para que la prueba valga en los dos repositorios.
    assert main(["recetas", "primero", "--language", "es"]) == 0
    assert "La primera sesión" in capsys.readouterr().out
