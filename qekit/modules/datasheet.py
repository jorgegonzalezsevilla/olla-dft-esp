# Olla-DFT — command-line toolkit for Quantum ESPRESSO
# Copyright (C) 2026 Jorge Enrique González Sevilla
# SPDX-License-Identifier: AGPL-3.0-or-later
# This program is free software: you can redistribute it and/or modify it
# under the terms of the GNU Affero General Public License as published by the Free
# Software Foundation, either version 3 of the License, or (at your option)
# any later version. See the LICENSE file for details.

"""Ficha del material y párrafo de métodos.

Cada módulo de Olla-DFT escribe su propio reporte. Lo que falta al final es
la SÍNTESIS: un documento por material con todo junto y, sobre todo, con
la procedencia de cada número — de qué cálculo salió, con qué parámetros.

Y el subproducto que más tiempo ahorra: el párrafo de metodología. Todos
los parámetros ya están guardados en los XML y en los encabezados de
procedencia, así que redactarlo a mano solo sirve para equivocarse. Aquí
se genera del propio cálculo, junto con la lista de códigos a citar.

Se escribe en Markdown y en HTML. El HTML es autocontenido (sin archivos
externos) para que se pueda abrir en cualquier máquina.
"""

from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

from qekit import __version__
from qekit.core import provenance, qeout

CITAS = {
    "qe": ("P. Giannozzi et al., J. Phys.: Condens. Matter 21, 395502 "
           "(2009); J. Phys.: Condens. Matter 29, 465901 (2017)"),
    "dfpt": ("S. Baroni, S. de Gironcoli, A. Dal Corso, P. Giannozzi, "
             "Rev. Mod. Phys. 73, 515 (2001)"),
    "spglib": "A. Togo, I. Tanaka, arXiv:1808.01590 (2018)",
    "seekpath": ("Y. Hinuma, G. Pizzi, Y. Kumagai, F. Oba, I. Tanaka, "
                 "Comput. Mater. Sci. 128, 140 (2017)"),
    "ase": ("A. Hjorth Larsen et al., J. Phys.: Condens. Matter 29, "
            "273002 (2017)"),
    "mace": ("I. Batatia et al., arXiv:2401.00096 (2023) — MACE-MP-0"),
    "epsilon": ("A. Dal Corso, S. Baroni, R. Resta, Phys. Rev. B 49, 5323 "
                "(1994) — respuesta dieléctrica"),
}


@dataclass
class Ficha:
    formula: str = ""
    resultados: dict = field(default_factory=dict)   # sección -> filas
    procedencia: dict = field(default_factory=dict)  # sección -> ruta
    parametros: dict = field(default_factory=dict)
    codigos: list = field(default_factory=list)
    avisos: list = field(default_factory=list)


def _fila(nombre, valor, unidad="", nota=""):
    return {"magnitud": nombre, "valor": valor, "unidad": unidad,
            "nota": nota}


def recoger(project=".") -> Ficha:
    """Recorre la carpeta y junta lo que cada módulo haya dejado."""
    p = Path(project)
    f = Ficha()

    # --- estructura y parámetros, del primer XML que aparezca ---------
    try:
        xml = qeout.find_xml(str(p))
        r = qeout.read_xml(xml)
        from ase import Atoms
        f.formula = Atoms(symbols=r.symbols).get_chemical_formula()
        f.parametros = {
            "funcional": r.functional,
            "ecutwfc_Ry": r.ecutwfc, "ecutrho_Ry": r.ecutrho,
            "malla_k": "x".join(map(str, r.kgrid)) if r.kgrid else None,
            "pseudos": dict(r.pseudo_files),
            "ocupaciones": r.occupations_kind,
            "smearing": r.smearing or None,
            "degauss_Ry": r.degauss,
            "nspin": r.nspin,
        }
        f.resultados["Estructura"] = [
            _fila("volumen de celda", round(r.volume, 4), "Å³"),
            _fila("átomos por celda", len(r.symbols), ""),
            _fila("operaciones de simetría", r.n_sym, ""),
        ]
        f.procedencia["Estructura"] = str(xml)
        f.codigos.append("qe")
    except Exception as exc:                           # noqa: BLE001
        f.avisos.append(f"no se pudo leer ningún XML de pw.x: {exc}")

    # --- resultados de cada módulo, por su archivo característico -----
    def _buscar(patron):
        hits = list(p.rglob(patron))
        return hits[0] if hits else None

    eos = _buscar("EOS.txt")
    if eos:
        filas = []
        for linea in Path(eos).read_text(errors="ignore").splitlines():
            for clave, nombre, unidad in (("B0", "módulo volumétrico B₀", "GPa"),
                                          ("a0", "parámetro de red a₀", "Å"),
                                          ("V0", "volumen de equilibrio", "Å³")):
                if linea.strip().startswith(clave):
                    tok = [t for t in linea.replace("=", " ").split()
                           if _es_num(t)]
                    if tok:
                        filas.append(_fila(nombre, float(tok[0]), unidad))
        if filas:
            f.resultados["Ecuación de estado"] = filas
            f.procedencia["Ecuación de estado"] = str(eos)

    el = _buscar("ELASTIC_C.dat")
    if el:
        C = np.loadtxt(el, comments="#")
        if C.shape == (6, 6):
            from qekit.modules import elastic
            m = elastic.moduli(C)
            f.resultados["Elásticas"] = [
                _fila("C₁₁", round(float(C[0, 0]), 2), "GPa"),
                _fila("C₁₂", round(float(C[0, 1]), 2), "GPa"),
                _fila("C₄₄", round(float(C[3, 3]), 2), "GPa"),
                _fila("módulo volumétrico (VRH)", round(m.B_hill, 2), "GPa"),
                _fila("módulo de corte (VRH)", round(m.G_hill, 2), "GPa"),
                _fila("razón de Poisson", round(m.nu, 4), ""),
                _fila("estable (Born)", "sí" if m.stable else "NO", ""),
            ]
            f.procedencia["Elásticas"] = str(el)

    fon = _buscar("FONONES_TERMO.dat")
    if fon:
        d = np.loadtxt(fon, comments="#")
        cab = Path(fon).read_text(errors="ignore")
        zpe = None
        for linea in cab.splitlines():
            if "ZPE" in linea:
                tok = [t for t in linea.replace("=", " ").split()
                       if _es_num(t)]
                if tok:
                    zpe = float(tok[0]) * 1000.0
        i300 = int(np.argmin(np.abs(d[:, 0] - 300.0)))
        filas = [_fila("C_v (300 K)", round(d[i300, 4] * 1000, 4),
                       "meV/K por celda")]
        if zpe is not None:
            filas.insert(0, _fila("energía de punto cero", round(zpe, 2),
                                  "meV por celda"))
        f.resultados["Fonones"] = filas
        f.procedencia["Fonones"] = str(fon)
        f.codigos.append("dfpt")

    op = _buscar("OPTICS.dat")
    if op:
        d = np.loadtxt(op, comments="#")
        f.resultados["Ópticas"] = [
            _fila("ε₁(0)", round(float(d[0, 1]), 3), ""),
            _fila("n(0)", round(float(d[0, 3]), 4), ""),
        ]
        f.procedencia["Ópticas"] = str(op)
        f.codigos.append("epsilon")

    mef = _buscar("MASA_EFECTIVA.dat")
    if mef:
        filas = []
        for linea in Path(mef).read_text(errors="ignore").splitlines():
            if linea.startswith("#") or not linea.strip():
                continue
            t = linea.split()
            if len(t) >= 3 and _es_num(t[2]):
                filas.append(_fila(f"m* {t[0]} ({' '.join(t[6:])})",
                                   float(t[2]), "mₑ"))
        if filas:
            f.resultados["Masa efectiva"] = filas[:6]
            f.procedencia["Masa efectiva"] = str(mef)

    der = _buscar("DERIVED.dat")
    if der:
        filas = []
        for linea in Path(der).read_text(errors="ignore").splitlines():
            if linea.startswith("#") or not linea.strip():
                continue
            t = linea.split()
            if len(t) >= 2 and _es_num(t[1]):
                nom = t[0].replace("_", " ")
                uni = t[2] if len(t) > 2 else ""
                filas.append(_fila(nom, round(float(t[1]), 4), uni))
        if filas:
            f.resultados["Derivadas termoelásticas"] = filas
            f.procedencia["Derivadas termoelásticas"] = str(der)

    gam = _buscar("GAMMA.dat")
    if gam:
        try:
            d = np.loadtxt(gam, comments="#")
            d = d[None, :] if d.ndim == 1 else d
            cab = Path(gam).read_text(errors="ignore")
            g = None
            for linea in cab.splitlines():
                if linea.startswith("#") and "gamma" in linea.lower() \
                        and "J/m" in linea:
                    tok = [t for t in linea.replace("=", " ").split()
                           if _es_num(t)]
                    if tok:
                        g = float(tok[0])
            if g is not None:
                f.resultados["Superficie"] = [
                    _fila("energía de superficie γ", round(g, 4), "J/m²",
                          "del ajuste de Fiorentini–Methfessel, no de una "
                          "losa suelta")]
                f.procedencia["Superficie"] = str(gam)
        except Exception:                                   # noqa: BLE001
            pass

    wfd = _buscar("WF.dat")
    esm = _buscar("ESM.dat")
    filas_sup = []
    if wfd:
        phi = _cab_num(wfd, "Phi_eV")
        if phi is not None:
            filas_sup.append(_fila("función trabajo Φ", round(phi, 3), "eV",
                                   "meseta del potencial planar"))
    if esm:
        try:
            d = np.loadtxt(esm, comments="#")
            d = d[None, :] if d.ndim == 1 else d
            i = int(np.argmin(np.abs(d[:, 0])))
            filas_sup.append(_fila("función trabajo Φ (ESM)",
                                   round(float(d[i, 4]), 3), "eV",
                                   "nivel de vacío fijado a cero por la "
                                   "condición de contorno"))
        except Exception:                                   # noqa: BLE001
            pass
    if filas_sup:
        f.resultados["Función trabajo"] = filas_sup
        f.procedencia["Función trabajo"] = str(esm or wfd)

    kap = _buscar("KAPPA.dat")
    if kap:
        try:
            d = np.loadtxt(kap, comments="#")
            d = d[None, :] if d.ndim == 1 else d
            i = int(np.argmin(np.abs(d[:, 0] - 300.0)))
            f.resultados["Conductividad térmica de red"] = [
                _fila(f"κ_L ({d[i, 0]:.0f} K)", round(float(d[i, -1]), 2),
                      "W/m·K",
                      "RTA, solo tres fonones: subestima entre un 10 y un "
                      "15 %")]
            f.procedencia["Conductividad térmica de red"] = str(kap)
            f.codigos.append("phono3py")
        except Exception:                                   # noqa: BLE001
            pass

    ber = _buscar("BERRY.dat")
    if ber:
        try:
            d = np.loadtxt(ber, comments="#")
            d = d[None, :] if d.ndim == 1 else d
            filas = [_fila("P (proyección sobre R)",
                           round(float(d[-1, -1]), 5), "C/m²",
                           "definida módulo el cuanto: solo la DIFERENCIA "
                           "a lo largo de un camino es física")]
            if len(d) > 1:
                filas.append(_fila("ΔP a lo largo del camino",
                                   round(float(d[-1, -1] - d[0, -1]), 5),
                                   "C/m²"))
            f.resultados["Polarización"] = filas
            f.procedencia["Polarización"] = str(ber)
        except Exception:                                   # noqa: BLE001
            pass

    wan = _buscar("WANNIER_centros.dat")
    if wan:
        try:
            cab = Path(wan).read_text(errors="ignore")
            om = None
            for linea in cab.splitlines():
                if linea.startswith("#") and "Omega =" in linea:
                    om = float(linea.split("Omega =")[1].split()[0])
            d = np.loadtxt(wan, comments="#")
            d = d[None, :] if d.ndim == 1 else d
            filas = [_fila("funciones de Wannier", len(d), "")]
            if om is not None:
                filas.append(_fila("dispersión total Ω", round(om, 4), "Å²"))
                filas.append(_fila("Ω por función", round(om / len(d), 4),
                                   "Å²"))
            f.resultados["Funciones de Wannier"] = filas
            f.procedencia["Funciones de Wannier"] = str(wan)
        except Exception:                                   # noqa: BLE001
            pass

    ech = _buscar("ECHEM.dat")
    if ech:
        eta = _cab_num(ech, "sobrepotencial")
        if eta is not None:
            f.resultados["Electroquímica"] = [
                _fila("sobrepotencial η", round(eta, 3), "V")]
            f.procedencia["Electroquímica"] = str(ech)

    form = _buscar("FORMACION.dat")
    if form:
        try:
            d = np.loadtxt(form, comments="#")
            d = d[None, :] if d.ndim == 1 else d
            f.resultados["Defectos"] = [
                _fila("E_f más baja en el rango de E_F",
                      round(float(np.nanmin(d[:, 1:])), 4), "eV",
                      "depende de μ y de E_F: mira la envolvente completa")]
            f.procedencia["Defectos"] = str(form)
        except Exception:                                   # noqa: BLE001
            pass

    if _buscar("MLIP_PROCEDENCIA.json"):
        f.codigos.append("mace")
        f.avisos.append(
            "hay estructuras producidas por un potencial aprendido en esta "
            "carpeta. Revisa\nque no se hayan mezclado sus energías con las "
            "de DFT (olla-dft audit lo comprueba).")

    f.codigos = list(dict.fromkeys(f.codigos + ["spglib", "seekpath", "ase"]))
    return f


def _cab_num(path, clave):
    """Un número escrito como '# clave = valor' en la cabecera de un .dat."""
    for linea in Path(path).read_text(errors="ignore").splitlines():
        if not linea.startswith("#") or clave not in linea:
            continue
        tok = [t for t in linea.split("=")[-1].split() if _es_num(t)]
        if tok:
            return float(tok[0])
    return None


def _es_num(t):
    try:
        float(t)
        return True
    except ValueError:
        return False


def metodos(f: Ficha) -> str:
    """Párrafo de metodología redactado a partir de los parámetros reales."""
    p = f.parametros
    if not p:
        return ("No hay parámetros que redactar: no se encontró ningún XML "
                "de pw.x en la carpeta.")
    pseudos = ", ".join(f"{k} ({v})" for k, v in (p.get("pseudos") or {}).items())
    ocup = p.get("ocupaciones") or "?"
    if ocup == "smearing" and p.get("degauss_Ry"):
        ocup = (f"ensanchamiento de tipo {p.get('smearing')} con degauss = "
                f"{p['degauss_Ry']:g} Ry")
    elif ocup == "fixed":
        ocup = "ocupaciones fijas"
    texto = (
        f"Los cálculos de primeros principios se realizaron con Quantum "
        f"ESPRESSO [1], en el marco de la teoría del funcional de la "
        f"densidad. Se empleó el funcional de intercambio y correlación "
        f"{p.get('funcional') or '?'} y los pseudopotenciales {pseudos or '?'}. "
        f"Las funciones de onda y la densidad de carga se expandieron en "
        f"ondas planas con energías de corte de {p.get('ecutwfc_Ry') or '?'} "
        f"y {p.get('ecutrho_Ry') or '?'} Ry respectivamente. "
        f"La zona de Brillouin se muestreó con una malla uniforme centrada "
        f"en Γ de {p.get('malla_k') or '?'} (K_POINTS automatic, sin "
        f"desplazamiento), con {ocup}. "
    )
    if p.get("nspin") == 2:
        texto += "Los cálculos se realizaron con polarización de espín. "
    if "dfpt" in f.codigos:
        texto += ("Las propiedades vibracionales se obtuvieron por teoría "
                  "del funcional de la densidad perturbativa [2]. ")
    if "mace" in f.codigos:
        texto += ("Las geometrías de partida se pre-relajaron con el "
                  "potencial interatómico aprendido MACE-MP-0; los "
                  "resultados reportados provienen en todos los casos de "
                  "los cálculos DFT posteriores. ")
    texto += (f"El análisis de simetría y la generación de los caminos de "
              f"alta simetría se hicieron con spglib y seekpath. "
              f"El pre y post-proceso se realizó con Olla-DFT {__version__}.")
    return texto


def citas(f: Ficha) -> list:
    return [CITAS[c] for c in f.codigos if c in CITAS]


def markdown(f: Ficha) -> str:
    lines = [f"# Ficha del material: {f.formula or '?'}", "",
             f"*Generado por Olla-DFT {__version__} — "
             f"{provenance.fields()['generado']}*", ""]
    for seccion, filas in f.resultados.items():
        lines += [f"## {seccion}", "",
                  "| magnitud | valor | unidad |", "|---|---|---|"]
        for r in filas:
            lines.append(f"| {r['magnitud']} | {r['valor']} | "
                         f"{r['unidad']} |")
        if seccion in f.procedencia:
            lines += ["", f"*Fuente: `{f.procedencia[seccion]}`*"]
        lines.append("")
    lines += ["## Parámetros del cálculo", "",
              "| parámetro | valor |", "|---|---|"]
    for k, v in f.parametros.items():
        if isinstance(v, dict):
            v = ", ".join(f"{a}: {b}" for a, b in v.items())
        lines.append(f"| {k} | {v if v is not None else '—'} |")
    lines += ["", "## Métodos (borrador)", "", metodos(f), "",
              "## Referencias", ""]
    for i, c in enumerate(citas(f), start=1):
        lines.append(f"{i}. {c}")
    if f.avisos:
        lines += ["", "## Avisos", ""]
        for a in f.avisos:
            lines.append(f"- {a}")
    lines += ["", "---", "",
              "*El párrafo de métodos es un BORRADOR generado de los "
              "parámetros reales del\ncálculo. Revísalo antes de usarlo: "
              "sabe qué se hizo, no por qué se hizo.*"]
    return "\n".join(lines)


def html(f: Ficha) -> str:
    import html as _h
    filas_html = []
    for seccion, filas in f.resultados.items():
        filas_html.append(f"<h2>{_h.escape(seccion)}</h2><table>"
                          "<tr><th>magnitud</th><th>valor</th>"
                          "<th>unidad</th></tr>")
        for r in filas:
            filas_html.append(
                f"<tr><td>{_h.escape(str(r['magnitud']))}</td>"
                f"<td class='v'>{_h.escape(str(r['valor']))}</td>"
                f"<td>{_h.escape(str(r['unidad']))}</td></tr>")
        filas_html.append("</table>")
        if seccion in f.procedencia:
            filas_html.append(
                f"<p class='src'>Fuente: <code>"
                f"{_h.escape(f.procedencia[seccion])}</code></p>")
    params = "".join(
        f"<tr><td>{_h.escape(k)}</td><td>"
        f"{_h.escape(', '.join(f'{a}: {b}' for a, b in v.items()) if isinstance(v, dict) else str(v))}"
        "</td></tr>" for k, v in f.parametros.items())
    refs = "".join(f"<li>{_h.escape(c)}</li>" for c in citas(f))
    avisos = ""
    if f.avisos:
        avisos = "<h2>Avisos</h2><ul>" + "".join(
            f"<li>{_h.escape(a)}</li>" for a in f.avisos) + "</ul>"
    return f"""<!doctype html><html lang="es"><head><meta charset="utf-8">
<title>Ficha — {_h.escape(f.formula or '?')}</title><style>
body{{font-family:system-ui,sans-serif;max-width:52rem;margin:2rem auto;
padding:0 1rem;line-height:1.5;color:#1a1a1a}}
h1{{border-bottom:2px solid #333;padding-bottom:.3rem}}
h2{{margin-top:1.8rem;color:#333}}
table{{border-collapse:collapse;width:100%;margin:.5rem 0}}
th,td{{border:1px solid #ccc;padding:.35rem .6rem;text-align:left}}
th{{background:#f0f0f0}} td.v{{font-variant-numeric:tabular-nums;
text-align:right}}
.src{{font-size:.85em;color:#666;margin-top:-.2rem}}
.met{{background:#f7f7f7;padding:1rem;border-left:3px solid #888}}
</style></head><body>
<h1>Ficha del material: {_h.escape(f.formula or '?')}</h1>
<p class="src">Generado por Olla-DFT {__version__} —
{provenance.fields()['generado']}</p>
{''.join(filas_html)}
<h2>Parámetros del cálculo</h2><table>
<tr><th>parámetro</th><th>valor</th></tr>{params}</table>
<h2>Métodos (borrador)</h2>
<div class="met"><p>{_h.escape(metodos(f))}</p></div>
<h2>Referencias</h2><ol>{refs}</ol>
{avisos}
<hr><p class="src">El párrafo de métodos es un borrador generado de los
parámetros reales del cálculo. Revísalo antes de usarlo: sabe qué se hizo,
no por qué se hizo.</p></body></html>"""


def escribir(f: Ficha, outdir: str = ".", nombre: str = None) -> list:
    out = Path(outdir); out.mkdir(parents=True, exist_ok=True)
    base = nombre or f"ficha_{f.formula or 'material'}"
    md = out / f"{base}.md"
    ht = out / f"{base}.html"
    md.write_text(markdown(f), encoding="utf-8")
    ht.write_text(html(f), encoding="utf-8")
    return [str(md), str(ht)]
