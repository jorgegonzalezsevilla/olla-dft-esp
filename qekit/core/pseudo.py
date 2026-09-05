# Olla-DFT — command-line toolkit for Quantum ESPRESSO
# Copyright (C) 2026 Jorge Enrique González Sevilla
# SPDX-License-Identifier: AGPL-3.0-or-later
# This program is free software: you can redistribute it and/or modify it
# under the terms of the GNU Affero General Public License as published by the Free
# Software Foundation, either version 3 of the License, or (at your option)
# any later version. See the LICENSE file for details.

"""Manejo de pseudopotenciales UPF.

- Busca el archivo UPF de cada elemento dentro de pseudo_dir.
- Lee los cutoffs sugeridos del encabezado del UPF (si el archivo los trae),
  para proponer ecutwfc/ecutrho automáticamente, como haría VASPKIT con
  los POTCAR.
"""

import re
from pathlib import Path

# patrones UPF v2:  wfc_cutoff="30.0"  rho_cutoff="240.0"
_RE_WFC_V2 = re.compile(r'wfc_cutoff\s*=\s*"([\d.eE+\-]+)"')
_RE_RHO_V2 = re.compile(r'rho_cutoff\s*=\s*"([\d.eE+\-]+)"')
_RE_TYPE_V2 = re.compile(r'pseudo_type\s*=\s*"([A-Za-z]+)"')
_RE_REL_V2 = re.compile(r'relativistic\s*=\s*"([A-Za-z]+)"')
_RE_ZVAL_V2 = re.compile(r'z_valence\s*=\s*"([\d.eE+\-]+)"')
_RE_ZVAL_V1 = re.compile(r"([\d.]+)\s+Z valence")
# patrones UPF v1 (texto libre en PP_INFO)
_RE_WFC_V1 = re.compile(
    r"[Ss]uggested\s+(?:minimum\s+)?cutoff\s+for\s+wavefunctions?\s*:?\s*([\d.]+)"
)
_RE_RHO_V1 = re.compile(
    r"[Ss]uggested\s+(?:minimum\s+)?cutoff\s+for\s+charge\s+density\s*:?\s*([\d.]+)"
)


def find_for_element(symbol: str, pseudo_dir: Path) -> list:
    """Archivos UPF candidatos para un elemento, ordenados alfabéticamente.

    Acepta nombres tipo 'Zn.UPF', 'Zn.pbe-dn-rrkjus_psl.1.0.0.UPF',
    'zn_pbe_v1.uspp.F.UPF', etc.: el nombre debe empezar con el símbolo
    (sin importar mayúsculas) seguido de un separador no alfabético.
    """
    if not pseudo_dir.is_dir():
        return []
    matches = []
    pattern = re.compile(rf"^{re.escape(symbol)}(?![a-zA-Z])", re.IGNORECASE)
    for path in sorted(pseudo_dir.iterdir()):
        if path.suffix.lower() != ".upf":
            continue
        if pattern.match(path.name):
            matches.append(path)
    return matches


def suggested_cutoffs(upf_path: Path) -> tuple:
    """(ecutwfc, ecutrho) sugeridos en el encabezado del UPF, en Ry.

    Devuelve (None, None) si el archivo no los declara o son cero.
    Solo se lee el inicio del archivo (el encabezado), no todo el UPF.
    """
    try:
        head = upf_path.read_text(errors="ignore")[:20000]
    except OSError:
        return (None, None)

    def _grab(pat_v2, pat_v1):
        m = pat_v2.search(head) or pat_v1.search(head)
        if m:
            try:
                val = float(m.group(1))
                return val if val > 1.0 else None
            except ValueError:
                return None
        return None

    return (_grab(_RE_WFC_V2, _RE_WFC_V1), _grab(_RE_RHO_V2, _RE_RHO_V1))


def pseudo_type(upf_path: Path):
    """Tipo del pseudopotencial: 'NC', 'US', 'PAW', 'SL' o None.

    Se lee del encabezado UPF v2; para UPF v1 se buscan las palabras clave
    en el PP_INFO. epsilon.x solo funciona con norma conservada (NC).
    """
    try:
        head = upf_path.read_text(errors="ignore")[:20000]
    except OSError:
        return None
    m = _RE_TYPE_V2.search(head)
    if m:
        t = m.group(1).upper()
        return {"NC": "NC", "US": "US", "USPP": "US", "PAW": "PAW",
                "SL": "SL", "1/R": "SL"}.get(t, t)
    low = head.lower()
    if "ultrasoft" in low:
        return "US"
    if "paw" in low:
        return "PAW"
    if "norm-conserving" in low or "norm conserving" in low:
        return "NC"
    return None


def relativistic(upf_path: Path):
    """'full', 'scalar', 'nonrelativistic' o None.

    El acoplamiento espín-órbita SOLO tiene sentido con pseudos
    totalmente relativistas ('full'): un pseudo escalar-relativista ya
    promedió el SOC y con lspinorb daría un desdoblamiento cero
    disfrazado de resultado.
    """
    try:
        head = upf_path.read_text(errors="ignore")[:20000]
    except OSError:
        return None
    m = _RE_REL_V2.search(head)
    if m:
        v = m.group(1).lower()
        # UPF v2 escribe "no" para no relativista; devolverlo tal cual
        # produce mensajes como "es no:" que no dicen nada.
        return "nonrelativistic" if v in ("no", "none", "n") else v
    low = head.lower()
    # UPF v1 lo escribe en prosa y de varias formas: "Full-Relativistic",
    # "Fully Relativistic", "full relativistic". Reconocer solo una deja
    # pseudos relativistas marcados como desconocidos, y entonces el filtro
    # de espin-orbita los descarta o los deja pasar por error.
    if re.search(r"full[- ]?relativistic|fully[- ]?relativistic", low):
        return "full"
    if re.search(r"scalar[- ]?relativistic", low):
        return "scalar"
    if re.search(r"non[- ]?relativistic", low):
        return "nonrelativistic"
    return None


def z_valence(upf_path: Path):
    """Electrones de valencia declarados en el UPF (o None)."""
    try:
        head = upf_path.read_text(errors="ignore")[:20000]
    except OSError:
        return None
    m = _RE_ZVAL_V2.search(head) or _RE_ZVAL_V1.search(head)
    if m:
        try:
            return float(m.group(1))
        except ValueError:
            return None
    return None



def _elegir(symbol, candidates, forzado=None, tarea=None, funcional=None):
    """Cual de los UPF de un elemento se usa.

    Por orden: el que el usuario forzo, el que el selector recomienda para
    la tarea, y en ultimo lugar el primero por orden alfabetico, que es lo
    que se hacia antes y funcionaba por casualidad.
    """
    if forzado:
        for c in candidates:
            if c.name == forzado:
                return c
        from qekit.core.errors import ErrorDeUso
        raise ErrorDeUso(
            f"pediste '{forzado}' para {symbol} y no esta en la carpeta de "
            "pseudopotenciales.\nLos que hay: "
            + ", ".join(c.name for c in candidates))
    if tarea and len(candidates) > 1:
        try:
            from qekit.modules import pseudos as _pz
            cands = [_pz.leer(c) for c in candidates]
            ev = _pz.evaluar(cands, tarea, funcional)
            buenos = [c for c in ev if c.ok]
            if buenos:
                return Path(buenos[0].ruta)
        except Exception:                              # noqa: BLE001
            pass
    return candidates[0]

def resolve(symbols: list, pseudo_dir: str, exclude=None, tarea: str = None,
            forzados: dict = None, funcional: str = None) -> dict:
    """Asigna un UPF a cada elemento de la lista.

    Devuelve un dict:
      simbolo -> {"filename": str, "found": bool, "alternatives": [...],
                  "ecutwfc": float|None, "ecutrho": float|None}
    Si no encuentra el archivo, usa '<Símbolo>.UPF' como marcador y
    found=False para que la CLI avise al usuario.
    """
    pdir = Path(pseudo_dir).expanduser()
    # `exclude` saca de la baraja archivos que NO deben elegirse solos.
    # El caso que lo motivó: un pseudo con hueco de core en la misma
    # carpeta ganaba por orden alfabético y se usaba para los átomos
    # normales; los dos tipos quedaban idénticos y XPS devolvía ceros.
    fuera = {Path(x).name for x in (exclude or [])}
    forzados = {k: Path(v).name for k, v in (forzados or {}).items()}
    result = {}
    for symbol in dict.fromkeys(symbols):  # únicos, en orden
        candidates = [c for c in find_for_element(symbol, pdir)
                      if c.name not in fuera]
        if candidates:
            chosen = _elegir(symbol, candidates, forzados.get(symbol),
                             tarea, funcional)
            wfc, rho = suggested_cutoffs(chosen)
            result[symbol] = {
                "filename": chosen.name,
                "found": True,
                "alternatives": [p.name for p in candidates[1:]],
                "ecutwfc": wfc,
                "ecutrho": rho,
                "z_valence": z_valence(chosen),
                "type": pseudo_type(chosen),
                "relativistic": relativistic(chosen),
            }
        else:
            result[symbol] = {
                "filename": f"{symbol}.UPF",
                "found": False,
                "alternatives": [],
                "ecutwfc": None,
                "ecutrho": None,
                "z_valence": None,
                "type": None,
                "relativistic": None,
            }
    return _coherencia_de_funcional(result, pdir, forzados)


#: Orden de preferencia entre funcionales cuando cubren lo mismo.
PREFERENCIA_FUNCIONAL = {"PBE": 0, "PBESOL": 1, "REVPBE": 2, "PZ": 3,
                         "BLYP": 4}


def _coherencia_de_funcional(result: dict, pdir, forzados: dict) -> dict:
    """Si los pseudos elegidos mezclan funcionales, busca una combinacion
    que no lo haga.

    Es la correccion silenciosa mas valiosa del selector. Por orden
    alfabetico, un NiO acababa con Ni de PBE y O de BLYP: QE no se queja y
    la energia total no significa nada. Aqui se busca el funcional que
    cubra mas elementos y se reeligen los que se puedan; los forzados a
    mano NO se tocan, porque el usuario mando.
    """
    from qekit.modules import pseudos as _pz

    encontrados = {k: v for k, v in result.items() if v["found"]}
    if len(encontrados) < 2:
        return result
    try:
        opciones = {k: [_pz.leer(f) for f in find_for_element(k, pdir)]
                    for k in encontrados}
    except Exception:                                  # noqa: BLE001
        return result

    actuales = {}
    for k, cands in opciones.items():
        for c in cands:
            if c.nombre == result[k]["filename"]:
                actuales[k] = c.funcional
                break
    presentes = {f for f in actuales.values() if f}
    if len(presentes) <= 1:
        return result

    # que funcional cubre mas elementos
    cobertura = {}
    for k, cands in opciones.items():
        for f in {c.funcional for c in cands if c.funcional}:
            cobertura[f] = cobertura.get(f, 0) + 1
    if not cobertura:
        return result
    # A igual cobertura manda la preferencia, no el alfabeto: PBE antes
    # que PZ. Sin esto el desempate alfabetico elegia LDA por sistema.
    mejor = max(cobertura,
                key=lambda f: (cobertura[f],
                               -PREFERENCIA_FUNCIONAL.get(f, 99), f))
    if cobertura[mejor] < len(encontrados):
        return result           # ninguno los cubre todos: se deja como esta

    for k, cands in opciones.items():
        if k in (forzados or {}):
            continue
        if actuales.get(k) == mejor:
            continue
        for c in cands:
            if c.funcional == mejor:
                p = Path(c.ruta)
                wfc, rho = suggested_cutoffs(p)
                result[k].update({
                    "filename": p.name, "ecutwfc": wfc, "ecutrho": rho,
                    "z_valence": z_valence(p), "type": pseudo_type(p),
                    "relativistic": relativistic(p),
                    "alternatives": [x.nombre for x in cands
                                     if x.nombre != p.name],
                })
                break
    return result


def recommend_cutoffs(pseudos: dict, default_wfc: float, dual: float) -> tuple:
    """(ecutwfc, ecutrho) para el sistema: el máximo sobre los elementos.

    Si ningún UPF declara cutoffs, se usan default_wfc y default_wfc*dual.
    """
    wfcs = [p["ecutwfc"] for p in pseudos.values() if p["ecutwfc"]]
    rhos = [p["ecutrho"] for p in pseudos.values() if p["ecutrho"]]
    ecutwfc = max(wfcs) if wfcs else default_wfc
    ecutrho = max(rhos) if rhos else ecutwfc * dual
    # ecutrho nunca por debajo de 4*ecutwfc (mínimo físico)
    ecutrho = max(ecutrho, 4.0 * ecutwfc)
    return (round(ecutwfc, 1), round(ecutrho, 1))
