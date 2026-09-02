"""Los README de ``examples/`` no pueden quedarse obsoletos.

Cada carpeta de ejemplo enseña los comandos ``olla-dft`` exactos que
produjeron sus datos. Esta prueba recorre todos los ``examples/**/README.md``,
saca cada línea que empieza por ``olla-dft `` (en bloques de código o
indentada con cuatro espacios) y la valida contra el PROPIO árbol de
argparse, igual que ``test_recetas.py`` hace con las recetas: un
subcomando que ya no existe o una bandera renombrada rompen la prueba
antes de que lleguen a un usuario.
"""

import argparse
import re
import shlex
from pathlib import Path

import pytest

RAIZ = Path(__file__).resolve().parent.parent
EJEMPLOS = RAIZ / "examples"


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


def comandos_de_readme(texto):
    """Las líneas ``olla-dft ...`` de un README, con su número de línea.

    Se aceptan dentro de bloques ``` ``` y las indentadas con cuatro
    espacios; un comentario ``# ...`` al final de la línea se descarta.
    """
    dentro = False
    for n, linea in enumerate(texto.splitlines(), 1):
        if linea.lstrip().startswith("```"):
            dentro = not dentro
            continue
        if not (dentro or linea.startswith("    ")):
            continue
        cmd = linea.strip()
        if not cmd.startswith("olla-dft "):
            continue
        cmd = re.sub(r"\s+#.*$", "", cmd)
        yield n, cmd


def _readmes():
    return sorted(EJEMPLOS.glob("**/README.md"))


def _validar(cmd, subs):
    """Devuelve la lista de problemas de un comando (vacía si está bien)."""
    trozos = shlex.split(cmd)[1:]
    if not trozos:
        return ["línea vacía"]
    nombre = trozos[0]
    if nombre not in subs:
        return [f"subcomando '{nombre}' no existe"]
    validas = _opciones(subs[nombre])
    malos = []
    for t in trozos[1:]:
        if not t.startswith("-") or t == "-":
            continue
        base = t.split("=")[0]
        if base in validas:
            continue
        # puede ser un valor negativo de la opción anterior
        if len(base) > 1 and base[1].isdigit():
            continue
        malos.append(f"'olla-dft {nombre}' no acepta {base}")
    return malos


def test_hay_readmes_en_todas_las_carpetas_de_ejemplo():
    carpetas = sorted(p for p in EJEMPLOS.iterdir() if p.is_dir())
    assert carpetas
    sin = [c.name for c in carpetas if not (c / "README.md").exists()]
    assert not sin, f"carpetas sin README.md: {sin}"
    assert (EJEMPLOS / "README.md").exists()
    assert not list(EJEMPLOS.glob("**/LEEME.txt")), "quedan LEEME.txt viejos"


def test_los_readmes_no_usan_el_nombre_viejo():
    for rd in _readmes():
        texto = rd.read_text(encoding="utf-8")
        assert not re.search(r"\bqekit\b", texto, re.IGNORECASE), rd
        assert not re.search(r"\bv?0\.\d+\.\d+\b", texto), \
            f"{rd}: no anclar el ejemplo a una versión"


def test_cada_readme_enseña_al_menos_un_comando():
    for rd in _readmes():
        cmds = list(comandos_de_readme(rd.read_text(encoding="utf-8")))
        assert cmds, f"{rd} no tiene ningún comando olla-dft"


def test_cada_readme_esta_en_el_idioma_del_repositorio():
    """Cada repositorio (español o inglés) lleva sus README en un solo idioma."""
    from qekit.core.i18n import DEFAULT_LANGUAGE
    marcas = {"es": ("Español", "English"), "en": ("English", "Español")}
    propia, ajena = marcas[DEFAULT_LANGUAGE]
    for rd in _readmes():
        texto = rd.read_text(encoding="utf-8")
        assert f"## {ajena}" not in texto, f"{rd} mezcla idiomas"


@pytest.mark.parametrize("readme", _readmes(),
                         ids=lambda p: str(p.relative_to(EJEMPLOS)))
def test_todos_los_comandos_de_los_ejemplos_existen_de_verdad(readme):
    """Contra el árbol de argparse, no contra una lista escrita aparte."""
    from qekit.cli import build_parser
    subs = _subparsers(build_parser())
    malos = []
    for n, cmd in comandos_de_readme(readme.read_text(encoding="utf-8")):
        for problema in _validar(cmd, subs):
            malos.append(f"{readme.name}:{n}: {problema}  ({cmd})")
    assert not malos, "\n".join(malos)


def test_los_archivos_citados_en_cada_readme_existen():
    """Todo `archivo` en la tabla de archivos tiene que estar en la carpeta."""
    for rd in _readmes():
        texto = rd.read_text(encoding="utf-8")
        faltan = []
        for fila in texto.splitlines():
            if not fila.startswith("| `"):
                continue
            celda = fila.split("|")[1]
            for nombre in re.findall(r"`([^`]+)`", celda):
                if nombre.startswith("."):          # `.png` abreviado
                    continue
                if not (rd.parent / nombre).exists():
                    faltan.append(nombre)
        assert not faltan, f"{rd}: cita archivos que no existen: {faltan}"


def test_el_readme_general_lista_todas_las_carpetas_y_estructuras():
    texto = (EJEMPLOS / "README.md").read_text(encoding="utf-8")
    for p in sorted(EJEMPLOS.iterdir()):
        if p.name == "README.md":
            continue
        nombre = p.name + "/" if p.is_dir() else p.name
        assert f"`{nombre}`" in texto, f"examples/README.md no lista {nombre}"
