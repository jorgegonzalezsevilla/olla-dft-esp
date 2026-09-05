# Olla-DFT — command-line toolkit for Quantum ESPRESSO
# Copyright (C) 2026 Jorge Enrique González Sevilla
# SPDX-License-Identifier: AGPL-3.0-or-later
# This program is free software: you can redistribute it and/or modify it
# under the terms of the GNU Affero General Public License as published by the Free
# Software Foundation, either version 3 of the License, or (at your option)
# any later version. See the LICENSE file for details.

"""Alineamiento de bandas: dónde queda el VBM de uno respecto al del otro.

Dos materiales calculados por separado NO comparten cero de energías. Cada
cálculo periódico fija su potencial hasta una constante arbitraria (el
término G=0 de Hartree), así que restar directamente los VBM de dos corridas
da un número perfectamente formado que no significa nada. Es probablemente
el error más caro de la literatura de heterouniones, porque el resultado
tiene el aspecto correcto.

Hay dos maneras honestas de ponerlos en la misma escala:

  --mode vacio      Cada lado es una LOSA con vacío, y se referencian los dos
                    al nivel de vacío de su propio cálculo. Es lo que se hace
                    en materiales 2D y superficies, y no necesita construir
                    la interfaz. El precio es que el resultado es el de las
                    superficies aisladas: ignora la reorganización de carga
                    al ponerlas en contacto.

  --mode interfaz   El método riguroso (Van de Walle-Martin): se calculan los
                    dos BULTOS y además la INTERFAZ, y se usa el potencial
                    electrostático macroscópicamente promediado de la
                    interfaz como puente entre las dos escalas:

                        ΔE_v = (E_v − V̄)_A − (E_v − V̄)_B + (V̄_A − V̄_B)_interfaz

                    Los dos primeros términos son propiedades de bulto; el
                    tercero es lo único que aporta la interfaz, y es donde
                    entra el dipolo del contacto.

El tipo de alineamiento (I, II o III) sale de comparar los dos gaps una vez
alineados, y es lo que decide si una heterounión sirve para separar cargas.
"""

from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

from qekit.core import provenance, qeout
from qekit.core.errors import ErrorDeUso, FaltanDatos
from qekit.core import style as qstyle

# Por debajo de esto los dos materiales están, a efectos prácticos, alineados.
# DFT con funcionales semilocales no distingue offsets más finos: el error
# típico frente al experimento es de varias décimas de eV.
TOL_ALINEADOS = 0.05

TIPOS = {
    "I": "tipo I (anidado): los dos portadores caen en el mismo material",
    "II": "tipo II (escalonado): electrón y hueco se separan en materiales "
          "distintos",
    "III": "tipo III (roto): los gaps no se solapan; es un semimetal de "
           "interfaz",
    "=": "sin offset apreciable: los dos quedan, a efectos prácticos, "
         "en el mismo sitio",
}


@dataclass
class Lado:
    nombre: str = ""
    ruta: str = ""
    vbm: float = None            # eV, en la escala del propio cálculo
    cbm: float = None
    gap: float = None
    fermi: float = None
    referencia: float = None     # nivel de vacío o V̄ de bulto, misma escala
    ref_tipo: str = ""
    planitud: float = None       # eV, variación en la meseta de vacío
    es_metal: bool = False

    @property
    def vbm_rel(self):
        """VBM medido desde la referencia (negativo: por debajo del vacío)."""
        if self.vbm is None or self.referencia is None:
            return None
        return self.vbm - self.referencia

    @property
    def cbm_rel(self):
        if self.cbm is None or self.referencia is None:
            return None
        return self.cbm - self.referencia


@dataclass
class Alineamiento:
    a: Lado = None
    b: Lado = None
    modo: str = "vacio"
    delta_v: float = None        # eV, offset de la banda de valencia
    delta_c: float = None        # eV, offset de la banda de conducción
    puente: float = None         # V̄_A − V̄_B de la interfaz (modo interfaz)
    tipo: str = ""
    avisos: list = field(default_factory=list)


# ----------------------------------------------------------------------
# Leer un lado
# ----------------------------------------------------------------------
def _potencial(ruta, eje: int = 2, rerun: bool = False,
               pw_cmd: str = None, nproc: int = None):
    from qekit.modules import fields

    cube = Path(ruta) / "potencial.cube"
    if not cube.exists() or rerun:
        cube = fields.run_pp(str(ruta), "potential", "potencial",
                             pw_cmd=pw_cmd, nproc=nproc)
    return fields.read_cube(str(cube))


def leer_lado(ruta, nombre: str = None, modo: str = "vacio", eje: int = 2,
              rerun: bool = False, pw_cmd: str = None,
              nproc: int = None) -> Lado:
    """VBM, CBM y la referencia de energías de un cálculo."""
    from qekit.modules import fields

    ruta = Path(ruta)
    try:
        qe = qeout.read_xml(str(ruta))
    except Exception as exc:                                # noqa: BLE001
        raise ErrorDeUso(
            f"no pude leer un resultado de QE en '{ruta}': {exc}") from None

    lado = Lado(nombre=nombre or ruta.name, ruta=str(ruta),
                vbm=qe.homo, cbm=qe.lumo, fermi=qe.fermi)
    if qe.homo is not None and qe.lumo is not None:
        lado.gap = qe.lumo - qe.homo
    if qe.homo is None:
        raise ErrorDeUso(
            f"'{ruta}' no da un VBM. En un metal no hay banda de valencia que "
            f"alinear; y si es un aislante, al cálculo le faltan bandas "
            f"vacías (nbnd) o no usó occupations='fixed'.")
    if lado.cbm is None:
        lado.es_metal = True

    cube = _potencial(ruta, eje, rerun, pw_cmd, nproc)
    if modo == "vacio":
        wf = fields.work_function(cube, qe.fermi if qe.fermi is not None
                                  else qe.homo, axis=eje,
                                  positions=qe.positions)
        lado.referencia = wf.v_vacuum
        lado.ref_tipo = "nivel de vacío"
        lado.planitud = wf.flatness
    else:
        z, prof = fields.planar_average(cube, axis=eje)
        lado.referencia = float(np.mean(prof)) * fields.RY_EV
        lado.ref_tipo = "potencial medio de la celda"
        lado.planitud = None
    return lado


def puente_interfaz(ruta, eje: int = 2, ancho: float = None,
                    rerun: bool = False, pw_cmd: str = None,
                    nproc: int = None) -> dict:
    """V̄ de cada mitad de la interfaz, con promedio macroscópico.

    El promedio planar oscila con la periodicidad de cada red; lo que hay que
    comparar es el promedio MACROSCÓPICO, que es el planar convolucionado con
    una ventana del tamaño del periodo de cada lado. Sin ese segundo promedio
    el "puente" depende de en qué plano atómico se decida mirar, y cambia
    varias décimas de eV con esa elección.
    """
    from qekit.modules import fields

    cube = _potencial(ruta, eje, rerun, pw_cmd, nproc)
    z, prof = fields.planar_average(cube, axis=eje)
    prof = np.asarray(prof) * fields.RY_EV
    L = float(z[-1] - z[0]) if len(z) > 1 else 1.0
    ventana = ancho if ancho else L / 8.0
    macro = fields.macroscopic_average(np.asarray(z), prof, ventana)

    n = len(macro)
    cuarto = max(2, n // 8)
    izq = float(np.mean(macro[cuarto:2 * cuarto]))
    der = float(np.mean(macro[n // 2 + cuarto:n // 2 + 2 * cuarto]))
    return {"V_a": izq, "V_b": der, "delta": izq - der,
            "z": np.asarray(z), "planar": prof, "macro": macro,
            "ventana": ventana}


# ----------------------------------------------------------------------
# Alinear
# ----------------------------------------------------------------------
def alinear(a: Lado, b: Lado, modo: str = "vacio",
            puente: float = None) -> Alineamiento:
    al = Alineamiento(a=a, b=b, modo=modo, puente=puente)
    if a.vbm_rel is None or b.vbm_rel is None:
        raise FaltanDatos("falta la referencia de energías de alguno de los "
                          "dos lados.")
    al.delta_v = a.vbm_rel - b.vbm_rel + (puente or 0.0)
    if a.cbm_rel is not None and b.cbm_rel is not None:
        al.delta_c = a.cbm_rel - b.cbm_rel + (puente or 0.0)

    if a.gap and b.gap and al.delta_v is not None and al.delta_c is not None:
        # POSICIONES en la escala de B, no offsets. El VBM de B es el cero,
        # así que el VBM de A está en ΔE_v; pero su CBM está en gap_B + ΔE_c,
        # no en ΔE_c. Confundirlos clasificaba como tipo III dos materiales
        # perfectamente alineados, porque ΔE_c ≈ 0 quedaba "por debajo" del
        # VBM de B.
        v_a, c_a = al.delta_v, b.gap + al.delta_c
        v_b, c_b = 0.0, b.gap
        if (abs(al.delta_v) < TOL_ALINEADOS
                and abs(al.delta_c) < TOL_ALINEADOS):
            al.tipo = "="
        elif (v_a <= v_b and c_a >= c_b) or (v_b <= v_a and c_b >= c_a):
            al.tipo = "I"
        elif c_a <= v_b or c_b <= v_a:
            al.tipo = "III"
        else:
            al.tipo = "II"

    for lado in (a, b):
        if lado.planitud is not None and lado.planitud > 0.05:
            al.avisos.append(
                f"la meseta de vacío de {lado.nombre} varía "
                f"{lado.planitud:.3f} eV. O falta vacío, o la losa tiene "
                f"dipolo neto: usa --dipole al generarla. El nivel de vacío "
                f"es la referencia de todo esto, así que ese error entra "
                f"entero en el offset.")
        if lado.es_metal:
            al.avisos.append(
                f"{lado.nombre} no tiene CBM: es un metal, o le faltan bandas "
                f"vacías. Solo se puede dar el offset de valencia.")
    if modo == "vacio":
        al.avisos.append(
            "Modo vacío: son las dos superficies AISLADAS. Al ponerlas en "
            "contacto se transfiere carga y aparece un dipolo de interfaz que "
            "desplaza el offset, típicamente entre 0.1 y 0.5 eV. Para incluirlo "
            "hace falta calcular la interfaz (--mode interfaz).")
    return al


def report(al: Alineamiento) -> str:
    a, b = al.a, al.b
    L = ["--- Alineamiento de bandas ---",
         f"Modo: {'nivel de vacío' if al.modo == 'vacio' else 'potencial de la interfaz'}",
         "",
         f"  {'':14s} {'VBM':>10s} {'CBM':>10s} {'gap':>9s}   referencia",
         "  " + "-" * 62]
    for lado in (a, b):
        cbm = (f"{lado.cbm_rel:10.4f}" if lado.cbm_rel is not None
               else f"{'—':>10s}")
        gap = f"{lado.gap:9.4f}" if lado.gap else f"{'—':>9s}"
        L.append(f"  {lado.nombre:14s} {lado.vbm_rel:10.4f} {cbm} {gap}   "
                 f"{lado.ref_tipo}")
    L.append("  (eV respecto a la referencia de cada lado)")

    if al.puente is not None:
        L += ["", f"Puente de la interfaz V̄_A − V̄_B = {al.puente:+.4f} eV"]

    L += ["", f"Offset de valencia   ΔE_v = {al.delta_v:+.4f} eV"]
    if al.delta_c is not None:
        L.append(f"Offset de conducción ΔE_c = {al.delta_c:+.4f} eV")
    L.append(f"  Positivo quiere decir que la banda de {a.nombre} queda por "
             f"ENCIMA de la de {b.nombre}.")

    if al.tipo:
        L += ["", f"Alineamiento {TIPOS[al.tipo]}"]
        if al.tipo == "II":
            # el electrón cae al CBM más bajo y el hueco sube al VBM más alto
            quien_e = a.nombre if (al.delta_c or 0) < 0 else b.nombre
            quien_h = a.nombre if (al.delta_v or 0) > 0 else b.nombre
            L.append(f"  El electrón se va a {quien_e} y el hueco a "
                     f"{quien_h}: es el alineamiento que sirve para separar "
                     f"cargas\n  (fotovoltaica, fotocatálisis).")
        elif al.tipo == "I":
            L.append("  Los dos portadores acaban en el mismo material: sirve "
                     "para confinar\n  y emitir luz, no para separar cargas.")
        elif al.tipo == "=":
            L.append(f"  Los dos offsets son menores que "
                     f"{TOL_ALINEADOS:g} eV. Decir de qué tipo es esta unión\n"
                     f"  sería leer ruido: con funcionales semilocales el "
                     f"error frente al\n  experimento es de varias décimas.")

    if al.avisos:
        L.append("")
        for aviso in al.avisos:
            L.append(f"AVISO: {aviso}")
    return "\n".join(L)


def export(al: Alineamiento, outdir: str = ".") -> list:
    out = Path(outdir); out.mkdir(parents=True, exist_ok=True)
    f = out / "ALINEAMIENTO.txt"
    f.write_text(report(al) + "\n")
    d = out / "ALINEAMIENTO.dat"
    lines = [provenance.header(
        "alineamiento de bandas",
        {"modo": al.modo, "dEv_eV": al.delta_v, "dEc_eV": al.delta_c,
         "tipo": al.tipo, "puente_eV": al.puente}),
        f"# {'lado':>14s} {'VBM_rel':>12s} {'CBM_rel':>12s} {'gap':>10s}"]
    for lado in (al.a, al.b):
        lines.append(
            f"{lado.nombre:>16s} {lado.vbm_rel:12.5f} "
            + (f"{lado.cbm_rel:12.5f}" if lado.cbm_rel is not None
               else f"{'nan':>12s}")
            + (f" {lado.gap:10.5f}" if lado.gap else f" {'nan':>10s}"))
    d.write_text("\n".join(lines) + "\n")
    return [str(d), str(f)]


def posiciones_en_escala_de_b(al: Alineamiento) -> dict:
    """Bordes de banda de los dos lados en la escala de B (VBM de B = 0).

    Es la misma convención que usa `alinear` para clasificar el tipo: el VBM
    de A está en ΔE_v, pero su CBM está en gap_B + ΔE_c, NO en ΔE_c (ΔE_c es
    un offset entre los dos CBM, no una posición). La figura dibujaba antes
    ΔE_c a secas y por eso el CBM de A salía por debajo de su propio VBM en
    cuanto los dos CBM estaban cerca. Se factoriza aquí para que el reporte,
    el export y la figura no puedan volver a discrepar.
    """
    v_b = 0.0
    c_b = al.b.gap if al.b.gap else 1.0
    v_a = al.delta_v
    if al.delta_c is not None:
        c_a = c_b + al.delta_c
    else:
        # un metal o sin bandas vacías: se dibuja el gap de A como caja
        c_a = v_a + (al.a.gap or 1.0)
    return {"v_a": v_a, "c_a": c_a, "v_b": v_b, "c_b": c_b}


def plot(al: Alineamiento, outfile: str = "alineamiento", formats="pdf,png",
         theme: str = None, size: str = None, family: str = None,
         background: str = None, palette=None, usetex: bool = None,
         width="single", journal: str = "generic", aspect: float = 0.80,
         mono: bool = False, dpi: int = None) -> list:
    """El diagrama de cajas que se publica: dos gaps, uno al lado del otro."""
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError as exc:                              # pragma: no cover
        raise RuntimeError("matplotlib no está instalado.") from exc
    if al.delta_v is None:
        raise FaltanDatos("no hay offsets que graficar.")

    st = qstyle.apply(theme, size=size, family=family, background=background,
                      palette=palette, usetex=usetex, mono=mono)
    fig, ax = qstyle.new_figure(width, journal, aspect)
    cols = qstyle.palette(2, mono=mono)

    # todo en la escala de B: su VBM en cero (misma convención que `alinear`)
    pos = posiciones_en_escala_de_b(al)
    v_a, c_a, v_b, c_b = pos["v_a"], pos["c_a"], pos["v_b"], pos["c_b"]

    for i, (nombre, v, c, col) in enumerate(
            [(al.a.nombre, v_a, c_a, cols[0]),
             (al.b.nombre, v_b, c_b, cols[1])]):
        x0 = i * 1.2
        ax.add_patch(plt.Rectangle((x0, -6.0), 1.0, 6.0 + v, color=col,
                                   alpha=0.30, lw=0))
        ax.add_patch(plt.Rectangle((x0, c), 1.0, 6.0 - c, color=col,
                                   alpha=0.30, lw=0))
        ax.plot([x0, x0 + 1.0], [v, v], color=col, lw=st["line"] * 1.6)
        ax.plot([x0, x0 + 1.0], [c, c], color=col, lw=st["line"] * 1.6)
        ax.annotate(nombre, xy=(x0 + 0.5, c + 0.25), ha="center",
                    fontsize=st["legend"], color=qstyle.INK)
        ax.annotate(f"{v:+.2f}", xy=(x0 + 0.5, v - 0.28), ha="center",
                    fontsize=st["legend"] * 0.9, color=qstyle.INK_SOFT)
        ax.annotate(f"{c:+.2f}", xy=(x0 + 0.5, c + 0.06), ha="center",
                    fontsize=st["legend"] * 0.9, color=qstyle.INK_SOFT)

    lo = min(v_a, v_b) - 1.2
    hi = max(c_a, c_b) + 1.2
    ax.set_ylim(lo, hi)
    ax.set_xlim(-0.25, 2.45)
    ax.set_xticks([])
    ax.set_ylabel("energía respecto al VBM de "
                  + qstyle.tex_safe(al.b.nombre) + " (eV)")
    ax.axhline(0.0, color=qstyle.INK_FAINT, lw=st["axis_line"],
               dashes=[3.5, 2.0])
    if al.tipo:
        ax.set_title(f"alineamiento tipo {al.tipo}", fontsize=st["legend"])
    written = qstyle.save(fig, outfile, formats, dpi=dpi, modulo="alineamiento")
    plt.close(fig)
    return written
