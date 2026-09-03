# Olla-DFT — command-line toolkit for Quantum ESPRESSO
# Copyright (C) 2026 Jorge Enrique González Sevilla
# SPDX-License-Identifier: GPL-3.0-or-later
# This program is free software: you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by the Free
# Software Foundation, either version 3 of the License, or (at your option)
# any later version. See the LICENSE file for details.

"""Referencia navegable, generada del propio código.

Un README de mil ochocientas líneas con sesenta subcomandos no es
documentación: es un archivo. Y una documentación escrita aparte se
desincroniza del código a la tercera versión.

Esto la saca de donde ya está la verdad: el árbol de argparse da los
comandos y sus banderas, y los docstrings de los módulos dan la física —
que en esta suite no son cuatro líneas de cortesía, sino los avisos que
explican por qué un cálculo sale mal. Nada que haya que mantener a mano,
así que nada que se pueda quedar viejo.

Sale una sola página HTML, sin dependencias externas ni conexión, con
buscador. Se abre con doble clic.
"""

import html
import inspect
import json
import re
from dataclasses import dataclass, field
from pathlib import Path

from qekit import __command_name__, __product_name__, __version__

_I18N_DIR = Path(__file__).resolve().parent.parent / "data" / "i18n"

# Los subcomandos agrupados por para qué sirven. El orden es el del flujo de
# trabajo real, no alfabético: quien busca algo suele saber en qué fase está,
# no cómo se llama el comando.
GRUPOS = [
    ("Estructura", "Leer, convertir y construir celdas.",
     ["info", "kpath", "prim", "conv", "supercell", "convert", "surface",
      "defect", "interface", "layers", "exfoliate"]),
    ("Generar cálculos", "Escribir los inputs de pw.x y su post-proceso.",
     ["gen", "templates", "pseudos", "config"]),
    ("Barridos", "Familias de cálculos que responden una pregunta.",
     ["converge", "eos", "elastic", "strain", "gamma", "adsorb", "eform"]),
    ("Estructura electrónica", "Bandas, densidad de estados y gap.",
     ["bands", "dos", "plot", "gap", "effmass", "fermi", "unfold", "align",
      "wannier", "topology"]),
    ("Espectroscopía", "Lo que se compara con un espectro medido.",
     ["optics", "tddft", "xrd", "xps", "xanes", "corehole"]),
    ("Vibraciones y térmica", "Fonones y todo lo que sale de ellos.",
     ["phonons", "qha", "thermochem", "derived", "md", "kappa"]),
    ("Transporte y electroquímica", "Portadores, conductancia y catálisis.",
     ["transport", "elph", "ballistic", "echem", "wf", "charge", "charges",
      "berry", "esm"]),
    ("Correlación", "DFT+U y sus parámetros.",
     ["hubbard"]),
    ("Caminos de reacción", "Barreras y estados de transición.",
     ["neb"]),
    ("Control de calidad", "Comprobar que lo calculado se sostiene.",
     ["doctor", "audit", "crosscheck", "selftest", "cost", "db", "hull",
      "report", "compare", "tune", "results", "campaign"]),
    ("Aprendizaje automático", "Modelos interatómicos para cribar barato.",
     ["mlip", "suggest", "amorphous"]),
    ("Salida", "Lo que se entrega al final.",
     ["datasheet", "docs"]),
    ("Ayuda", "Para cuando no sabes qué comando buscar.",
     ["start", "wizard", "recetas", "teoria", "sistema", "update"]),
    ("Proyecto", "Organizar y validar un proyecto reproducible.",
     ["project"]),
]

# De qué módulo sale la física de cada comando, para poder citar su docstring.
MODULO_DE = {
    "strain": "strain", "adsorb": "adsorb", "eform": "defects",
    "gamma": "surfen", "align": "align", "echem": "echem",
    "elastic": "elastic", "eos": "eos", "converge": "converge",
    "bands": "bands", "dos": "dos", "transport": "transport",
    "phonons": "phonons", "elph": "elph", "hubbard": "hubbard",
    "xanes": "xanes", "corehole": "corehole", "xps": "xps",
    "tddft": "tddft", "ballistic": "ballistic", "md": "dynamics",
    "thermochem": "thermochem", "qha": "qha", "neb": "neb",
    "unfold": "unfold", "interface": "interface", "mlip": "mlip",
    "selftest": "selftest", "cost": "cost", "optics": "optics",
    "doctor": "diagnose", "audit": "audit", "wizard": "wizard",
    "exfoliate": "exfoliate", "xrd": "xrd", "docs": "docs",
    "amorphous": "amorphous", "wannier": "wannier", "topology": "topology",
    "berry": "berry",
    "kappa": "kappa", "esm": "esm",
    "recetas": "recipes", "sistema": "plataforma", "teoria": "theory",
    "update": "update",
    "start": "onboarding",
    "project": "project",
    "compare": "compare", "tune": "tuning",
    "results": "results", "campaign": "campaign",
}

_GROUP_KEYS = (
    "structure", "generation", "sweeps", "electronic", "spectroscopy",
    "vibrations", "transport", "correlation", "reaction", "quality",
    "machine_learning", "output", "help", "project",
)


def _labels(language="es") -> dict:
    """Carga la capa de presentación de la referencia sin traducir la física."""
    if language not in ("es", "en"):
        raise ValueError("language debe ser es o en")
    if language == "en":
        target = _I18N_DIR / "docs_en.json"
        try:
            value = json.loads(target.read_text(encoding="utf-8"))
        except (OSError, json.JSONDecodeError) as exc:
            raise RuntimeError(f"no se pudo cargar la traducción {target}: {exc}") from None
        if not isinstance(value, dict):
            raise RuntimeError(f"la traducción {target} no es un objeto JSON")
        return value
    return {
        "lang": "es",
        "page_title": "referencia",
        "subtitle": "{total} subcomandos. Esta página se genera del propio código con olla-dft docs, así que no puede quedarse vieja.",
        "search": "Buscar comando u opción…",
        "search_label": "Buscar en la referencia",
        "no_matches": "No hay comandos que coincidan con la búsqueda.",
        "navigation": "Navegación de la referencia",
        "start_here": "Empezar aquí",
        "recipes_title": "Empezar aquí: recetas",
        "recipes_lede": "Sesiones completas, de la estructura al resultado. Lo que enseñan no son los comandos sueltos, sino qué archivo deja cada paso y qué paso posterior lo lee.",
        "full_run": "corrida entera",
        "not_run": "sin correr entera",
        "steps": "pasos",
        "reads": "lee",
        "from": "de",
        "leaves": "deja",
        "watch": "Ojo:",
        "after": "Y después:",
        "arguments": "Argumentos",
        "options": "opciones",
        "physics": "La física detrás",
        "one_of": "una de",
        "default": "por omisión",
        "footer": "Generado por olla-dft docs · suite privada, sin telemetría",
        "other_group": ["Otros", "Comandos que no encajan en los grupos de arriba."],
        "groups": {key: [title, lema] for key, (title, lema, _names)
                   in zip(_GROUP_KEYS, GRUPOS)},
        "command_summaries": {},
        "option_summaries": {},
    }


@dataclass
class Bandera:
    nombres: str
    ayuda: str
    metavar: str = ""
    opciones: list = field(default_factory=list)
    defecto: object = None
    requerida: bool = False


@dataclass
class Comando:
    nombre: str
    ayuda: str = ""
    posicionales: list = field(default_factory=list)
    banderas: list = field(default_factory=list)
    fisica: str = ""


def _limpio(t):
    return " ".join(str(t or "").split())


def extraer() -> list:
    """Recorre el árbol de argparse y devuelve los comandos con sus banderas."""
    from qekit.cli import build_parser

    parser = build_parser()
    sub = None
    for a in parser._actions:
        choices = getattr(a, "choices", None)
        # argparse usa una lista para opciones como --language y un dict para
        # subcomandos. Solo el segundo puede contener parsers hijos.
        if isinstance(choices, dict) and choices and all(
                hasattr(v, "_actions") for v in choices.values()):
            sub = choices
            break
    if not sub:
        return []

    fuera = []
    vistos = set()
    for nombre, p in sub.items():
        # Los alias en inglés (system, recipes, theory) comparten el parser
        # con el nombre en español: se documentan una sola vez.
        if id(p) in vistos:
            continue
        vistos.add(id(p))
        c = Comando(nombre=nombre, ayuda=_limpio(p.description or ""))
        for a in p._actions:
            if a.dest == "help":
                continue
            b = Bandera(
                nombres=", ".join(a.option_strings) if a.option_strings
                        else a.dest,
                ayuda=_limpio(a.help),
                metavar=str(a.metavar or ""),
                opciones=[str(x) for x in (a.choices or [])],
                defecto=a.default if a.default not in (None, False) else None,
                requerida=bool(getattr(a, "required", False))
                          or not a.option_strings)
            (c.posicionales if not a.option_strings else c.banderas).append(b)
        mod = MODULO_DE.get(nombre)
        if mod:
            try:
                m = __import__(f"qekit.modules.{mod}", fromlist=["x"])
                c.fisica = inspect.getdoc(m) or ""
            except Exception:                               # noqa: BLE001
                c.fisica = ""
        fuera.append(c)
    return fuera


def ayudas_cortas() -> dict:
    """El help de una línea que cada subcomando declara al registrarse."""
    from qekit.cli import build_parser

    parser = build_parser()
    for a in parser._actions:
        if getattr(a, "choices", None) and hasattr(a, "_choices_actions"):
            # ca.dest es "sistema (system)" cuando hay alias: se queda el nombre.
            return {ca.dest.split(" ")[0]: _limpio(ca.help)
                    for ca in a._choices_actions}
    return {}


# ----------------------------------------------------------------------
# HTML
# ----------------------------------------------------------------------
E = lambda s: html.escape(str(s), quote=True)          # noqa: E731


_RE_VINETA = re.compile(r"^(\s*)[-*\u2022]\s+(.*)$")


def _inline(t: str) -> str:
    t = re.sub(r"`([^`]+)`", r"<code>\1</code>", E(t))
    return re.sub(r"\*\*([^*]+)\*\*", r"<b>\1</b>", t)


def _parrafos(texto: str) -> str:
    """Docstring a HTML: párrafos, listas y bloques de código.

    Las listas hay que reconocerlas expresamente. Un docstring escribe las
    continuaciones de una viñeta con más sangría, y una regla ingenua de
    "cuatro espacios = código" convierte la mitad de cada lista en un
    bloque monoespaciado sin sentido. Se vio en el docstring de `strain`,
    donde tres continuaciones salieron como código.
    """
    if not texto:
        return ""
    fuera = []
    for bloque in re.split(r"\n\s*\n", texto.strip("\n")):
        lineas = [l for l in bloque.split("\n") if l.strip()]
        if not lineas:
            continue
        # ¿lleva viñetas? entonces es una lista, con lo que venga antes
        # como texto de entrada
        i_v = next((i for i, l in enumerate(lineas) if _RE_VINETA.match(l)),
                   None)
        if i_v is not None:
            if i_v:
                fuera.append(f"<p>{_inline(' '.join(x.strip() for x in lineas[:i_v]))}</p>")
            items, actual = [], None
            for l in lineas[i_v:]:
                m = _RE_VINETA.match(l)
                if m:
                    if actual is not None:
                        items.append(actual)
                    actual = m.group(2).strip()
                elif actual is not None:
                    actual += " " + l.strip()
            if actual is not None:
                items.append(actual)
            fuera.append("<ul>" + "".join(f"<li>{_inline(x)}</li>"
                                          for x in items) + "</ul>")
            continue
        # ¿todo el bloque va sangrado? entonces es código
        if all(l.startswith("    ") for l in lineas):
            fuera.append(f"<pre>{E(bloque.strip(chr(10)))}</pre>")
            continue
        fuera.append(f"<p>{_inline(' '.join(x.strip() for x in lineas))}</p>")
    return "\n".join(fuera)


def _tabla(banderas: list, labels=None) -> str:
    if not banderas:
        return ""
    labels = labels or _labels("es")
    filas = []
    for b in banderas:
        ayuda = b.ayuda
        if labels.get("option_summaries"):
            for option in b.nombres.split(", "):
                if option in labels["option_summaries"]:
                    ayuda = labels["option_summaries"][option]
                    break
        extra = []
        if b.opciones:
            extra.append(f"{labels['one_of']}: " + ", ".join(b.opciones))
        if b.defecto is not None:
            extra.append(f"{labels['default']} {b.defecto}")
        cola = (f' <span class="extra">({"; ".join(extra)})</span>'
                if extra else "")
        filas.append(
            f'<tr><td><code>{E(b.nombres)}'
            + (f' {E(b.metavar)}' if b.metavar else "")
            + f'</code></td><td>{E(ayuda)}{cola}</td></tr>')
    return ('<div class="tw"><table class="flags"><tbody>'
            + "".join(filas) + "</tbody></table></div>")


def _html_recetas(labels=None, language="es") -> str:
    """Las mismas recetas del CLI, en la página. Una sola fuente.

    Si aquí hubiera una copia escrita a mano, se desincronizaría con el CLI
    a la tercera versión y nadie se daría cuenta. Se generan de RECETAS,
    traducidas por la misma tabla que usa ``olla-dft recetas --language en``.
    """
    from qekit.core import i18n
    from qekit.modules import recipes as rec
    labels = labels or _labels("es")
    # para el buscador de la página valen las palabras de los dos idiomas
    kw_en = i18n.load_table("recipes_en").get("keywords", {})
    B = []
    B.append(f'<section class="grupo" id="recetas">'
             f'<h2>{E(labels["recipes_title"])}</h2>'
             f'<p class="lema">{E(labels["recipes_lede"])}</p>')
    for r in rec.recetas(language):
        marca = (f'<span class="ok">{E(labels["full_run"])}</span>'
                 if r.verificada else
                 f'<span class="wip">{E(labels["not_run"])}</span>')
        palabras = list(r.palabras) + list(kw_en.get(r.clave, []))
        B.append(f'<article class="cmd receta" id="r-{E(r.clave)}" '
                 f'data-buscar="{E(r.clave)} {E(r.pregunta)} '
                 f'{E(" ".join(palabras))}">'
                 f'<h3><code>{__command_name__} recetas {E(r.clave)}</code></h3>'
                 f'<p class="corta">«{E(r.pregunta)}»</p>'
                 f'<p class="para">{E(r.para_que)}</p>'
                 f'<p class="meta2">{len(r.pasos)} {E(labels["steps"])} · {E(r.coste)} · '
                 f'{marca}</p><ol class="pasos">')
        for paso in r.pasos:
            flujo = ""
            for texto, de in paso.lee:
                flujo += (f'<li class="lee">{E(labels["reads"])} {E(texto)} '
                          f'<em>[{E(labels["from"])} {E(de)}]</em></li>')
            for f_ in paso.escribe:
                flujo += f'<li class="deja">{E(labels["leaves"])} <code>{E(f_)}</code></li>'
            B.append(f'<li><code class="cmdline">{E(paso.comando)}</code>'
                     + (' <span class="qe">pw.x</span>'
                        if paso.corre_qe else "")
                     + f'<p>{E(paso.hace)}</p>'
                     + (f'<ul class="flujo">{flujo}</ul>' if flujo else "")
                     + (f'<pre class="sal">{E(paso.salida)}</pre>'
                        if paso.salida else "")
                     + (f'<p class="ojo"><b>{E(labels["watch"])}</b> {E(paso.ojo)}</p>'
                        if paso.ojo else "")
                     + '</li>')
        B.append("</ol>")
        if r.despues:
            B.append(f'<p class="para"><b>{E(labels["after"])}</b></p><ul class="flujo">'
                     + "".join(f'<li>{E(t)} — <code>{E(c)}</code></li>'
                                for t, c in r.despues) + "</ul>")
        B.append("</article>")
    B.append("</section>")
    return "".join(B)


def generar(destino: str = "olla-dft-docs.html", language="es") -> str:
    labels = _labels(language)
    comandos = {c.nombre: c for c in extraer()}
    cortas = ayudas_cortas()

    sueltos = [n for n in comandos
               if not any(n in g[2] for g in GRUPOS)]
    grupos = list(GRUPOS)
    if sueltos:
        grupos.append(("Otros", "Comandos que no encajan en los grupos de "
                                "arriba.", sorted(sueltos)))

    from qekit.modules import recipes as _rec
    nav = [f'<div class="ng"><h4>{E(labels["start_here"])}</h4><ul>'
           + "".join(f'<li><a href="#r-{E(r.clave)}">{E(r.clave)}</a></li>'
                     for r in _rec.RECETAS) + "</ul></div>"]
    cuerpo = []
    for group_index, (titulo, lema, nombres) in enumerate(grupos):
        presentes = [n for n in nombres if n in comandos]
        if not presentes:
            continue
        if language == "en":
            if group_index < len(_GROUP_KEYS):
                titulo, lema = labels["groups"][_GROUP_KEYS[group_index]]
            else:
                titulo, lema = labels["other_group"]
        gid = re.sub(r"[^a-z]+", "-", titulo.lower())
        nav.append(f'<div class="ng"><h4>{E(titulo)}</h4><ul>'
                   + "".join(f'<li><a href="#c-{E(n)}">{E(n)}</a></li>'
                             for n in presentes) + "</ul></div>")
        cuerpo.append(f'<section class="grupo" id="{gid}">'
                      f'<h2>{E(titulo)}</h2><p class="lema">{E(lema)}</p>')
        for n in presentes:
            c = comandos[n]
            corta = (labels.get("command_summaries", {}).get(n)
                     if language == "en" else None) or cortas.get(n) or c.ayuda
            fisica = _parrafos(c.fisica)
            cuerpo.append(
                f'<article class="cmd" id="c-{E(n)}" data-buscar="{E(n)} '
                f'{E(corta)} {E(c.ayuda)}">'
                f'<h3><code>{__command_name__} {E(n)}</code></h3>'
                f'<p class="corta">{E(corta)}</p>'
                + (f'<div class="pos"><h5>{E(labels["arguments"])}</h5>'
                   f'{_tabla(c.posicionales, labels)}</div>' if c.posicionales else "")
                + (f'<details class="banderas"><summary>{len(c.banderas)} '
                   f'{E(labels["options"])}</summary>{_tabla(c.banderas, labels)}</details>'
                   if c.banderas else "")
                + (f'<details class="fisica"><summary>{E(labels["physics"])}</summary>'
                   f'<div class="fis">{fisica}</div></details>'
                   if fisica else "")
                + "</article>")
        cuerpo.append("</section>")

    total = len(comandos)
    doc = _PLANTILLA.replace("{{VERSION}}", E(__version__)) \
                    .replace("{{PRODUCT}}", E(__product_name__)) \
                    .replace("{{TOTAL}}", str(total)) \
                    .replace("{{NAV}}", "".join(nav)) \
                    .replace("{{CUERPO}}", "".join(cuerpo)) \
                    .replace("{{LANG}}", E(labels["lang"])) \
                    .replace("{{PAGE_TITLE}}", E(labels["page_title"])) \
                    .replace("{{SUBTITLE}}", E(labels["subtitle"].format(total=total))) \
                    .replace("{{SEARCH}}", E(labels["search"])) \
                    .replace("{{SEARCH_LABEL}}", E(labels["search_label"])) \
                    .replace("{{NO_MATCHES}}", E(labels["no_matches"])) \
                    .replace("{{NAVIGATION}}", E(labels["navigation"])) \
                    .replace("{{FOOTER}}", E(labels["footer"])) \
                    .replace("{{RECETAS}}", _html_recetas(labels, language))
    p = Path(destino)
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text(doc, encoding="utf-8")
    return str(p)


_PLANTILLA = """<!doctype html>
<html lang="{{LANG}}"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1">
<title>{{PRODUCT}} {{VERSION}} — {{PAGE_TITLE}}</title>
<style>
:root{
  --bg:#F3F4F1; --card:#fff; --line:#D8DCD7; --soft:#EEF0EC;
  --tx:#1A1E1B; --tx2:#525C55; --tx3:#7F887F; --acc:#12665F;
  --mono:ui-monospace,"SF Mono",Menlo,Consolas,monospace;
  --sans:system-ui,-apple-system,"Segoe UI",Roboto,sans-serif;
}
@media (prefers-color-scheme:dark){:root{
  --bg:#0E1211; --card:#171C1A; --line:#2B3431; --soft:#1D2422;
  --tx:#DFE6E2; --tx2:#9AA6A0; --tx3:#6F7A75; --acc:#63C6BC;
}}
*{box-sizing:border-box}
body{margin:0;background:var(--bg);color:var(--tx);font-family:var(--sans);
     font-size:15px;line-height:1.6}
a{color:var(--acc);text-decoration:none}
a:hover{text-decoration:underline}
code{font-family:var(--mono);font-size:.88em}
.wrap{display:flex;gap:0;max-width:1280px;margin:0 auto;align-items:flex-start}
nav{width:230px;flex:none;position:sticky;top:0;max-height:100vh;
    overflow-y:auto;padding:22px 14px 40px;border-right:1px solid var(--line)}
nav h4{font-size:11px;text-transform:uppercase;letter-spacing:.1em;
  color:var(--tx3);margin:16px 0 6px;font-weight:600}
nav ul{list-style:none;margin:0;padding:0}
nav li a{display:block;padding:2px 6px;border-radius:4px;font-family:var(--mono);
  font-size:12.5px;color:var(--tx2)}
nav li a:hover{background:var(--soft);color:var(--acc);text-decoration:none}
main{flex:1;min-width:0;padding:22px 24px 80px}
header{margin-bottom:8px}
h1{font-size:26px;margin:0 0 4px}
.sub{color:var(--tx2);margin:0 0 16px;font-size:14px}
#q{width:100%;padding:9px 12px;border:1px solid var(--line);border-radius:6px;
   background:var(--card);color:var(--tx);font-size:14px;font-family:var(--sans)}
#q:focus{outline:2px solid var(--acc);outline-offset:1px}
.grupo{margin-top:34px}
.grupo h2{font-size:20px;margin:0 0 2px;padding-bottom:6px;
  border-bottom:2px solid var(--tx)}
.lema{color:var(--tx2);font-size:13.5px;margin:6px 0 14px}
.cmd{background:var(--card);border:1px solid var(--line);border-radius:6px;
  padding:14px 16px;margin-bottom:8px}
.cmd h3{margin:0 0 3px;font-size:15px}
.cmd h3 code{color:var(--acc);font-size:15px}
.corta{margin:0 0 10px;color:var(--tx2);font-size:13.5px}
h5{margin:10px 0 4px;font-size:11px;text-transform:uppercase;
  letter-spacing:.08em;color:var(--tx3)}
.tw{overflow-x:auto}
table.flags{border-collapse:collapse;width:100%;font-size:13px}
table.flags td{padding:4px 10px 4px 0;vertical-align:top;
  border-bottom:1px solid var(--soft)}
table.flags td:first-child{white-space:nowrap;width:1%}
table.flags td:first-child code{color:var(--tx)}
.extra{color:var(--tx3);font-size:12px}
details{margin-top:8px}
summary{cursor:pointer;font-size:12.5px;color:var(--tx2);
  padding:3px 0;user-select:none}
summary:hover{color:var(--acc)}
.fis{border-left:3px solid var(--acc);padding:2px 0 2px 14px;margin-top:6px;
  font-size:13.5px;color:var(--tx2)}
.fis p{margin:0 0 .7em}
.fis ul{margin:.3em 0 .8em;padding-left:1.15em}
.fis li{margin-bottom:.35em}
.fis pre{background:var(--soft);padding:8px 11px;border-radius:4px;
  overflow-x:auto;font-family:var(--mono);font-size:12px;margin:.5em 0}
.oculto{display:none!important}
.sr-only{position:absolute;width:1px;height:1px;padding:0;margin:-1px;
  overflow:hidden;clip:rect(0,0,0,0);white-space:nowrap;border:0}
.search-status{min-height:1.5em;color:var(--tx2);font-size:13px}
footer{margin-top:40px;padding-top:14px;border-top:1px solid var(--line);
  color:var(--tx3);font-size:12px;font-family:var(--mono)}
@media(max-width:860px){nav{display:none}main{padding:18px 16px 60px}}

/* --- recetas --- */
.receta .para{color:var(--tx2);margin:.4rem 0 .2rem;max-width:70ch}
.receta .meta2{font-family:var(--mono);font-size:.78rem;color:var(--tx3);
  margin:.2rem 0 .8rem}
.receta .ok{color:var(--acc);font-weight:600}
.receta .wip{color:var(--tx3)}
ol.pasos{counter-reset:paso;list-style:none;padding-left:0;margin:0}
ol.pasos>li{counter-increment:paso;position:relative;padding:.7rem 0 .7rem 2.4rem;
  border-top:1px solid var(--line)}
ol.pasos>li::before{content:counter(paso);position:absolute;left:0;top:.75rem;
  width:1.6rem;height:1.6rem;border-radius:50%;background:var(--soft);
  color:var(--acc);font-family:var(--mono);font-size:.8rem;display:flex;
  align-items:center;justify-content:center;font-weight:600}
ol.pasos code.cmdline{background:var(--soft);padding:.25rem .5rem;
  border-radius:4px;display:inline-block;font-size:.86rem}
ol.pasos>li>p{margin:.45rem 0;color:var(--tx2);max-width:70ch}
.qe{font-family:var(--mono);font-size:.7rem;color:var(--tx3);
  border:1px solid var(--line);border-radius:10px;padding:.05rem .4rem;
  margin-left:.4rem}
ul.flujo{list-style:none;padding-left:0;margin:.35rem 0;font-size:.83rem;
  font-family:var(--mono);color:var(--tx3)}
ul.flujo li{padding:.1rem 0}
ul.flujo li.lee::before{content:"← lee  ";color:var(--acc)}
ul.flujo li.deja::before{content:"→ deja ";color:var(--tx3)}
pre.sal{background:var(--soft);border-left:2px solid var(--acc);
  padding:.5rem .7rem;margin:.5rem 0;font-size:.8rem;overflow-x:auto;
  white-space:pre-wrap}
p.ojo{font-size:.85rem;color:var(--tx2);border-left:2px solid var(--line);
  padding-left:.7rem;margin:.5rem 0;max-width:70ch}
</style></head><body>
<div class="wrap">
<nav aria-label="{{NAVIGATION}}">{{NAV}}</nav>
<main>
<header>
  <h1>{{PRODUCT}} {{VERSION}}</h1>
  <p class="sub">{{SUBTITLE}}</p>
  <div role="search"><label class="sr-only" for="q">{{SEARCH_LABEL}}</label>
  <input id="q" type="search" placeholder="{{SEARCH}}"
         autocomplete="off" aria-describedby="search-status"></div>
  <p id="search-status" class="search-status" aria-live="polite"></p>
</header>
{{RECETAS}}
{{CUERPO}}
<footer>{{FOOTER}}</footer>
</main></div>
<script>
(function(){
  var q=document.getElementById('q');
  var status=document.getElementById('search-status');
  var cmds=[].slice.call(document.querySelectorAll('.cmd'));
  var grupos=[].slice.call(document.querySelectorAll('.grupo'));
  function filter(){
    var t=q.value.trim().toLowerCase();
    var visible=0;
    cmds.forEach(function(c){
      var hay=!t||c.dataset.buscar.toLowerCase().indexOf(t)>=0
             ||c.textContent.toLowerCase().indexOf(t)>=0;
      c.classList.toggle('oculto',!hay);
      if(hay)visible++;
    });
    grupos.forEach(function(g){
      var vis=g.querySelectorAll('.cmd:not(.oculto)').length;
      g.classList.toggle('oculto',vis===0);
    });
    if(status)status.textContent=t&&!visible?'{{NO_MATCHES}}':'';
  }
  q.addEventListener('input',filter);
  document.addEventListener('keydown',function(e){
    if(e.key==='/'&&document.activeElement!==q){e.preventDefault();q.focus();}
  });
})();
</script></body></html>
"""
