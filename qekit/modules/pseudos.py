# Olla-DFT — command-line toolkit for Quantum ESPRESSO
# Copyright (C) 2026 Jorge Enrique González Sevilla
# SPDX-License-Identifier: AGPL-3.0-or-later
# This program is free software: you can redistribute it and/or modify it
# under the terms of the GNU Affero General Public License as published by the Free
# Software Foundation, either version 3 of the License, or (at your option)
# any later version. See the LICENSE file for details.

"""Elegir pseudopotencial con criterio, no por orden alfabético.

EL PROBLEMA
-----------
En una carpeta de pseudopotenciales suele haber varios para el mismo
elemento: `Fe.pbe-nd-rrkjus.UPF`, `Fe.pbe-spn-kjpaw_psl.1.0.0.UPF`,
`Fe.rel-pbe-spn-kjpaw_psl.1.0.0.UPF`... Hasta ahora Olla-DFT tomaba el
primero por orden alfabético. Eso funciona por casualidad, y cuando falla
lo hace en silencio: un pseudo escalar-relativista con `lspinorb` da un
desdoblamiento espín-órbita de CERO sin ningún error, y un ultrasuave con
`epsilon.x` da un espectro entero que está mal.

La elección depende de lo que vayas a hacer, y esa información Olla-DFT ya la
tiene. Este módulo la usa.

LO QUE MIRA
-----------
- **El funcional.** Mezclar funcionales entre elementos de la misma
  estructura invalida la energía total. Es lo primero que se comprueba.
- **El tipo.** Norma conservada (NC), ultrasuave (US) o PAW. Los NC son
  caros pero los únicos que sirven para `epsilon.x`; los US son baratos y
  necesitan un `ecutrho` mucho mayor; PAW da mejores densidades.
- **Lo relativista.** Solo un pseudo 'full' sirve para espín-órbita.
- **Los electrones de valencia.** Más valencia (semicore) es más caro y
  más transferible; para metales de transición y DFT+U suele hacer falta.
- **El cutoff sugerido.** El que declara el propio archivo. Un pseudo que
  pide 90 Ry cuesta el doble que uno que pide 45.

CÓMO PUNTÚA
-----------
No hay un "mejor pseudopotencial": hay uno adecuado para lo que vas a
hacer. El módulo aplica REQUISITOS DUROS (que descartan) y PREFERENCIAS
(que ordenan), las dos declaradas en una tabla que se lee, y explica por
qué quedó fuera cada uno. La última palabra es siempre del usuario.
"""

import re
from dataclasses import dataclass, field
from pathlib import Path

from qekit.core import pseudo as ps
from qekit.core.errors import ErrorDeUso

#: Requisitos duros por tarea: descartan un pseudo, no lo penalizan.
#: Cada entrada es (nombre legible, comprobación, por qué).
TAREAS = {
    "optics": {
        "nombre": "óptica con epsilon.x",
        "tipo": ("NC",),
        "razon_tipo": "epsilon.x SOLO funciona con norma conservada. Con "
                      "ultrasuaves o PAW devuelve un espectro sin quejarse, "
                      "y está mal.",
    },
    "soc": {
        "nombre": "espín-órbita",
        "relativista": ("full",),
        "razon_rel": "el acoplamiento espín-órbita necesita un pseudo "
                     "TOTALMENTE relativista. Uno escalar ya promedió el "
                     "SOC y daría un desdoblamiento cero disfrazado de "
                     "resultado.",
    },
    "xanes": {
        "nombre": "XANES con xspectra.x",
        "gipaw": True,
        "razon_gipaw": "xspectra.x reconstruye la función de onda de todos "
                       "los electrones con GIPAW: el pseudo tiene que traer "
                       "esa información.",
    },
    "hubbard": {
        "nombre": "DFT+U",
        "prefiere_semicore": True,
    },
    "fonones": {
        "nombre": "fonones (DFPT)",
        "prefiere_tipo": ("NC", "US"),
        "razon_pref": "la DFPT con PAW es más frágil y más cara; con norma "
                      "conservada o ultrasuave va mejor.",
    },
    "general": {"nombre": "cálculo general"},
}



#: Por debajo de este Z el espin-orbita es despreciable en la practica.
Z_SOC = 19


def _es_ligero(elemento: str) -> bool:
    from qekit.core.atomconf import Z_DE
    return Z_DE.get((elemento or "").capitalize(), 999) < Z_SOC

@dataclass
class Candidato:
    ruta: str = ""
    nombre: str = ""
    elemento: str = ""
    tipo: str = None
    funcional: str = None
    relativista: str = None
    z_valence: float = None
    ecutwfc: float = None
    ecutrho: float = None
    gipaw: bool = False
    tamano_kb: int = 0
    descartado: str = ""          # motivo, si lo hay
    puntos: float = 0.0
    notas: list = field(default_factory=list)

    @property
    def ok(self) -> bool:
        return not self.descartado


_RE_FUNC = re.compile(r'functional\s*=\s*"([^"]*)"', re.I)
_RE_FUNC_V1 = re.compile(r"^\s*(.+?)\s+Exchange-Correlation functional",
                         re.M)


#: Nombre corto <- las cuatro piezas con que QE escribe el mismo funcional.
NOMBRE_CORTO = {
    "SLA PW PBX PBC": "PBE",
    "SLA PZ NOGX NOGC": "PZ",
    "SLA PW PSX PSC": "PBESOL",
    "SLA B88 LYP BLYP": "BLYP",
    "SLA PW RPB PBC": "REVPBE",
    "SLA PW PBE": "PBE",
}


def _funcional(head: str):
    """El funcional, normalizado al nombre corto cuando se reconoce.

    QE escribe unas veces 'PBE' y otras 'SLA PW PBX PBC': son lo mismo.
    Mostrar las dos formas en la misma tabla hace pensar que dos pseudos
    son incompatibles cuando no lo son.
    """
    m = _RE_FUNC.search(head) or _RE_FUNC_V1.search(head)
    if not m:
        return None
    piezas = m.group(1).upper().split()
    if not piezas:
        return None
    crudo = " ".join(piezas)
    if crudo in NOMBRE_CORTO:
        return NOMBRE_CORTO[crudo]
    # UPF v1 escribe "SLA PW PBX PBC   PBE  Exchange-Correlation functional":
    # las cuatro piezas Y el nombre corto al final. Si la ultima pieza ya es
    # un nombre reconocible, esa es la respuesta.
    if piezas[-1] in NOMBRE_CORTO.values():
        return piezas[-1]
    if len(piezas) > 4 and " ".join(piezas[:4]) in NOMBRE_CORTO:
        return NOMBRE_CORTO[" ".join(piezas[:4])]
    return crudo


def leer(ruta) -> Candidato:
    """Lo que se puede saber de un UPF sin cargarlo entero."""
    p = Path(ruta)
    head = p.read_text(errors="ignore")[:30000]
    c = Candidato(ruta=str(p), nombre=p.name,
                  tipo=ps.pseudo_type(p), funcional=_funcional(head),
                  relativista=ps.relativistic(p), z_valence=ps.z_valence(p),
                  gipaw="PP_GIPAW" in p.read_text(errors="ignore"),
                  tamano_kb=int(p.stat().st_size / 1024))
    c.ecutwfc, c.ecutrho = ps.suggested_cutoffs(p)
    m = re.search(r'element\s*=\s*"?\s*(\w+)', head, re.I) or \
        re.search(r"^\s*(\w+)\s+Element", head, re.M)
    c.elemento = (m.group(1).strip() if m else p.name.split(".")[0])
    return c


def candidatos(elemento: str, pseudo_dir: str) -> list:
    """Todos los UPF de un elemento en la carpeta, leídos."""
    pdir = Path(pseudo_dir).expanduser()
    if not pdir.is_dir():
        raise ErrorDeUso(
            f"la carpeta de pseudopotenciales '{pseudo_dir}' no existe.\n"
            "Se fija con:  olla-dft config set pseudo_dir /ruta/a/tus/pseudos")
    return [leer(f) for f in ps.find_for_element(elemento, pdir)]


# ----------------------------------------------------------------------
# Selección
# ----------------------------------------------------------------------
def evaluar(cands: list, tarea: str = "general",
            funcional: str = None, prefiere_ligero: bool = False) -> list:
    """Descarta los que no sirven y ordena el resto. No elige por ti."""
    if tarea not in TAREAS:
        raise ErrorDeUso(
            f"tarea '{tarea}' desconocida. Opciones: "
            + ", ".join(sorted(TAREAS)))
    reglas = TAREAS[tarea]

    for c in cands:
        c.descartado, c.puntos, c.notas = "", 0.0, []

        tipos = reglas.get("tipo")
        if tipos and c.tipo and c.tipo not in tipos:
            c.descartado = (f"es {c.tipo} y hace falta "
                            f"{'/'.join(tipos)}: " + reglas["razon_tipo"])
            continue
        rel = reglas.get("relativista")
        if rel and c.relativista and c.relativista not in rel:
            if _es_ligero(c.elemento):
                # En un elemento ligero el espin-orbita es despreciable y
                # no existen pseudos relativistas para casi ninguno.
                # Descartarlo dejaria sin opciones a un calculo que es
                # perfectamente valido.
                c.notas.append(
                    f"es {c.relativista}, pero en {c.elemento} (Z pequeno) "
                    "el espin-orbita es\n    despreciable: se puede usar "
                    "junto a pseudos relativistas de los\n    elementos "
                    "pesados.")
                c.puntos -= 0.5
            else:
                c.descartado = (f"es {c.relativista}: " + reglas["razon_rel"])
                continue
        if reglas.get("gipaw") and not c.gipaw:
            c.descartado = "no trae GIPAW: " + reglas["razon_gipaw"]
            continue
        if funcional and c.funcional and \
                not _mismo_funcional(c.funcional, funcional):
            c.descartado = (f"su funcional es {c.funcional} y el pedido es "
                            f"{funcional}: mezclar funcionales invalida la "
                            "energía total.")
            continue

        # preferencias
        pref = reglas.get("prefiere_tipo")
        if pref and c.tipo in pref:
            c.puntos += 2.0
            c.notas.append(f"tipo {c.tipo}: " + reglas.get("razon_pref", ""))
        if reglas.get("prefiere_semicore") and c.z_valence:
            c.puntos += 0.15 * c.z_valence
            c.notas.append(f"{c.z_valence:g} electrones de valencia: más "
                           "semicore es más transferible para DFT+U")
        if c.ecutwfc:
            # un cutoff bajo es dinero; se premia, pero poco
            c.puntos += max(0.0, (90.0 - c.ecutwfc) / 30.0)
        else:
            c.notas.append("no declara cutoff sugerido: habrá que converger "
                           "a ciegas")
            c.puntos -= 0.5
        if prefiere_ligero and c.tipo in ("US", "PAW"):
            c.puntos += 1.0
            c.notas.append("ultrasuave/PAW: menos ondas planas, más barato")
        if c.gipaw:
            c.puntos += 0.3
            c.notas.append("trae GIPAW: sirve también para XANES y RMN")
        if c.relativista == "full":
            c.puntos += 0.2

    return sorted(cands, key=lambda c: (not c.ok, -c.puntos, c.nombre))


def elegir(elemento: str, pseudo_dir: str, tarea: str = "general",
           funcional: str = None, prefiere_ligero: bool = False) -> tuple:
    """(el mejor candidato, la lista entera evaluada)."""
    cands = candidatos(elemento, pseudo_dir)
    if not cands:
        raise ErrorDeUso(
            f"no hay ningún pseudopotencial de {elemento} en "
            f"'{pseudo_dir}'.\nSe descargan de pseudo-dojo.org o de "
            "quantum-espresso.org/pseudopotentials.")
    ev = evaluar(cands, tarea, funcional, prefiere_ligero)
    buenos = [c for c in ev if c.ok]
    if not buenos:
        raise ErrorDeUso(
            f"hay {len(ev)} pseudopotencial(es) de {elemento} pero ninguno "
            f"sirve para {TAREAS[tarea]['nombre']}:\n" +
            "\n".join(f"  {c.nombre}: {c.descartado}" for c in ev))
    return buenos[0], ev


def _mismo_funcional(a: str, b: str) -> bool:
    """Compara funcionales tolerando las distintas formas de escribirlos.

    'PBE' y 'SLA PW PBX PBC' son el mismo funcional escrito de dos formas:
    la corta y la lista de los cuatro trozos. Compararlos como cadenas
    descartaría pseudos perfectamente válidos.
    """
    ALIAS = {
        "PBE": {"PBE", "SLA PW PBX PBC"},
        "PZ": {"PZ", "SLA PZ NOGX NOGC", "LDA"},
        "PBESOL": {"PBESOL", "SLA PW PSX PSC"},
        "BLYP": {"BLYP", "SLA B88 LYP BLYP"},
        "REVPBE": {"REVPBE", "SLA PW RPB PBC"},
    }
    na = " ".join(a.upper().split())
    nb = " ".join(b.upper().split())
    if na == nb:
        return True
    for grupo in ALIAS.values():
        if na in grupo and nb in grupo:
            return True
    return False


def coherencia(elegidos: dict) -> list:
    """Avisos si los pseudos elegidos para distintos elementos no casan."""
    avisos = []
    funcs = {c.funcional for c in elegidos.values() if c.funcional}
    if len(funcs) > 1:
        # puede que sean el mismo escrito distinto
        base = list(funcs)[0]
        if not all(_mismo_funcional(base, f) for f in funcs):
            avisos.append(
                "FUNCIONALES DISTINTOS entre elementos: "
                + "; ".join(f"{k}={v.funcional}" for k, v in elegidos.items()
                            if v.funcional)
                + ".\nLa energía total de una estructura con pseudos de "
                  "funcionales distintos no\nsignifica nada. Esto no da "
                  "error en QE: sale un número perfectamente\nplausible y "
                  "equivocado.")
    tipos = {c.tipo for c in elegidos.values() if c.tipo}
    if "NC" in tipos and tipos & {"US", "PAW"}:
        avisos.append(
            "Se mezclan norma conservada con ultrasuave/PAW. QE lo permite, "
            "pero el\necutrho lo manda el más exigente: usa el dual del "
            "ultrasuave (8-12x) para\ntodos, no el del NC (4x).")
    cut = [c.ecutwfc for c in elegidos.values() if c.ecutwfc]
    if cut and max(cut) / min(cut) > 2.5:
        peor = max(elegidos.values(), key=lambda c: c.ecutwfc or 0)
        avisos.append(
            f"Los cutoffs sugeridos van de {min(cut):.0f} a {max(cut):.0f} "
            f"Ry. Manda el mayor ({peor.nombre}),\nasí que ese elemento "
            "decide el coste de todo el cálculo. Si hay otro pseudo\nde ese "
            "elemento más blando, puede ahorrar mucho tiempo.")
    return avisos


# ----------------------------------------------------------------------
# Reporte
# ----------------------------------------------------------------------
def report(elemento: str, evaluados: list, tarea: str = "general") -> str:
    lines = [f"--- Pseudopotenciales de {elemento} ---",
             f"Para: {TAREAS[tarea]['nombre']}", ""]
    buenos = [c for c in evaluados if c.ok]
    malos = [c for c in evaluados if not c.ok]

    if buenos:
        lines.append(f"{'':2s} {'archivo':<42s} {'tipo':>5s} {'func':>12s} "
                     f"{'rel':>7s} {'zval':>5s} {'ecut':>6s} {'rho':>6s}")
        for i, c in enumerate(buenos):
            marca = "->" if i == 0 else "  "
            lines.append(
                f"{marca} {c.nombre[:42]:<42s} {(c.tipo or '?'):>5s} "
                f"{(c.funcional or '?')[:12]:>12s} "
                f"{(c.relativista or '?')[:7]:>7s} "
                f"{(f'{c.z_valence:g}' if c.z_valence else '?'):>5s} "
                f"{(f'{c.ecutwfc:.0f}' if c.ecutwfc else '-'):>6s} "
                f"{(f'{c.ecutrho:.0f}' if c.ecutrho else '-'):>6s}")
        lines += ["", f"Recomendado: {buenos[0].nombre}"]
        for n in buenos[0].notas:
            if n:
                lines.append(f"  - {n}")

    if malos:
        lines += ["", "Descartados:"]
        for c in malos:
            lines.append(f"  {c.nombre}")
            lines.append(f"      {c.descartado}")

    lines += ["",
              "Esto es una recomendación, no una verdad: el pseudopotencial "
              "adecuado depende\ndel sistema. Se fuerza uno concreto con "
              "--pseudo " + elemento + "=archivo.UPF, y sea cual sea\nel que "
              "elijas, hay que converger el cutoff con 'olla-dft converge'."]
    return "\n".join(lines)


def report_coherencia(elegidos: dict) -> str:
    avisos = coherencia(elegidos)
    if not avisos:
        return ("Los pseudopotenciales elegidos son coherentes entre sí "
                "(mismo funcional,\ntipos compatibles, cutoffs del mismo "
                "orden).")
    return "\n\n".join(avisos)
