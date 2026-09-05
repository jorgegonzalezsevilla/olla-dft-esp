# Olla-DFT — command-line toolkit for Quantum ESPRESSO
# Copyright (C) 2026 Jorge Enrique González Sevilla
# SPDX-License-Identifier: AGPL-3.0-or-later
# This program is free software: you can redistribute it and/or modify it
# under the terms of the GNU Affero General Public License as published by the Free
# Software Foundation, either version 3 of the License, or (at your option)
# any later version. See the LICENSE file for details.

"""Formato de intercambio con el resto de la suite.

Olla-DFT calcula desde primeros principios lo mismo que varias apps de la
suite miden en el laboratorio, y esos pares se comparan directamente:

    olla-dft xrd      ->  patrón calculado como REFERENCIA de fases para la
                       app de DRX (que hoy solo tiene bases de datos
                       externas; esto le da una referencia del material
                       propio, ya relajado)
    olla-dft optics   ->  Tauc y α(E) de DFT contra el Tauc experimental de
                       la app de UV-Vis
    olla-dft phonons --gamma ->  frecuencias y actividad IR contra el
                       espectro medido en las apps de FTIR y Raman

Para que ese puente exista hace falta UN formato, no uno por app. Aquí
está: un JSON con esquema versionado, en UTF-8, con las unidades escritas
en el nombre de cada campo (nada de "energia" a secas).

Regla de compatibilidad: dentro de una misma versión de esquema solo se
AÑADEN campos, nunca se renombran ni se quitan. Un lector que ignore lo
que no conoce seguirá funcionando. Si algo tiene que cambiar de forma,
sube `qekit_suite_schema`.

NOTA: escribir este archivo es la mitad del puente; la otra mitad es el
lector en cada app, que hay que añadir allá.
"""

import json
from pathlib import Path

import numpy as np

from qekit.core import provenance

SCHEMA = 1


def _clean(x):
    """Convierte arreglos y escalares de numpy a tipos nativos de JSON."""
    if isinstance(x, np.ndarray):
        return [_clean(v) for v in x.tolist()]
    if isinstance(x, (np.floating, np.integer)):
        x = x.item()
    if isinstance(x, float):
        if not np.isfinite(x):
            return None
        return round(x, 8)
    if isinstance(x, dict):
        return {k: _clean(v) for k, v in x.items()}
    if isinstance(x, (list, tuple)):
        return [_clean(v) for v in x]
    return x


def material_block(atoms=None) -> dict:
    """Identidad del material, para que la app sepa qué está comparando."""
    if atoms is None:
        return {}
    d = {"formula": atoms.get_chemical_formula(),
         "n_atomos": len(atoms),
         "volumen_A3": float(atoms.get_volume()),
         "celda_A": [[float(c) for c in fila] for fila in atoms.cell.array]}
    try:
        from qekit.core import structure as st
        ds = st.symmetry_dataset(atoms)
        d["grupo_espacial"] = ds.international
        d["grupo_espacial_numero"] = int(ds.number)
    except Exception:
        pass
    return _clean(d)


def envelope(tipo: str, datos: dict, atoms=None, notas: str = "") -> dict:
    """Arma el JSON completo con su cabecera de procedencia."""
    doc = {
        "qekit_suite_schema": SCHEMA,
        "tipo": tipo,
        "generado_por": _clean(provenance.fields(tipo)),
        "material": material_block(atoms),
        "datos": _clean(datos),
    }
    if notas:
        doc["notas"] = notas
    return doc


def write(doc: dict, path) -> str:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(doc, ensure_ascii=False, indent=2) + "\n")
    return str(path)


# ----------------------------------------------------------------------
# Constructores por módulo
# ----------------------------------------------------------------------
def from_xrd(pattern, atoms=None) -> dict:
    """Patrón de polvos calculado, como referencia para la app de DRX."""
    picos = [{"dos_theta_grados": float(p.two_theta),
              "d_A": float(p.d),
              "intensidad_rel": float(p.intensity),
              "hkl": p.label,
              "hkl_equivalentes": [list(h) for h in p.hkls]}
             for p in sorted(pattern.peaks, key=lambda q: q.two_theta)]
    # los hkl están referidos a la celda en que se indexó el patrón: el
    # bloque "material" tiene que describir ESA celda, no la de entrada, o
    # la app leería índices de una celda y parámetros de otra
    if atoms is not None and pattern.basis == "conventional":
        try:
            from qekit.core import structure as st
            atoms = st.conventional(atoms)
        except Exception:
            pass
    datos = {
        "lambda_A": float(pattern.wavelength),
        "base_indices_hkl": pattern.basis,
        "picos": picos,
    }
    if pattern.two_theta is not None:
        datos["perfil"] = {
            "dos_theta_grados": pattern.two_theta,
            "intensidad_rel": pattern.intensity,
            "fwhm_grados": pattern.fwhm,
            "tamano_cristalito_nm": pattern.size_nm,
        }
    return envelope(
        "drx_patron_calculado", datos, atoms,
        notas="Patrón simulado desde una estructura relajada por DFT. Las "
              "intensidades son |F|² con Lorentz-polarización; sin "
              "orientación preferencial ni corrección de rugosidad, así que "
              "difieren de un difractograma medido de polvo real.")


def from_optics(run, atoms=None, tauc_kind: str = "direct") -> dict:
    """ε(ω), α(E) y Tauc de DFT, para comparar con la app de UV-Vis."""
    from qekit.modules import optics as opt
    d = opt.derived(run)
    gap, pend, ventana, curva = opt.tauc_gap(run, tauc_kind)
    datos = {
        "energia_eV": run.energies,
        "eps1": run.eps1,
        "eps2": run.eps2,
        "n": d["n"],
        "k": d["k"],
        "alpha_cm-1": d["alpha"],
        "reflectividad": d["R"],
        "scissor_eV": float(run.scissor or 0.0),
        "ensanchamiento_eV": float(run.intersmear),
        "tauc": {
            "tipo": tauc_kind,
            "exponente": 2 if tauc_kind == "direct" else 0.5,
            "gap_eV": gap,
            "ventana_ajuste_eV": list(ventana) if ventana else None,
            "curva": curva,
        },
    }
    return envelope(
        "optica_dft", datos, atoms,
        notas="RPA de partícula independiente: sin campos locales ni "
              "excitones. Para comparar contra un espectro UV-Vis medido, "
              "revisa el scissor aplicado. En un semiconductor de gap "
              "indirecto epsilon.x no incluye transiciones asistidas por "
              "fonones, así que el borde calculado es el gap DIRECTO.")


def from_phonons_gamma(run, atoms=None) -> dict:
    """Frecuencias y actividad IR en Γ, para las apps de FTIR y Raman."""
    modos = []
    for i, (w, act) in enumerate(run.gamma_freqs, start=1):
        m = {"modo": i,
             "omega_cm-1": float(w),
             "omega_THz": float(w) * 0.0299792458}
        if act is not None:
            m["actividad_IR"] = float(act)
        modos.append(m)
    return envelope(
        "fonones_gamma", {"modos": modos,
                          "con_tensor_dielectrico": bool(run.epsil)}, atoms,
        notas="Frecuencias armónicas en Γ. La actividad IR es la del cálculo "
              "DFPT (unidades de QE), proporcional a la intensidad medida "
              "pero no igual a ella: compara posiciones de banda y "
              "actividad relativa, no valores absolutos. Sin anarmonicidad, "
              "así que las frecuencias suelen salir altas frente al FTIR.")


def from_raman(run, atoms=None, laser_nm: float = 532.0,
               T: float = 300.0, fwhm: float = 5.0) -> dict:
    """Actividades Raman y espectro simulado, para la app de Raman."""
    from qekit.modules import phonons
    w, inten, picos = phonons.raman_spectrum(
        run, laser_nm=laser_nm, T=T, fwhm=fwhm)
    modos = []
    for d in run.modes:
        m = {"modo": d["modo"], "omega_cm-1": d["omega_cm1"],
             "omega_THz": d["omega_thz"], "actividad_IR": d.get("ir")}
        if "raman" in d:
            m["actividad_Raman_A4_amu"] = d["raman"]
            m["factor_despolarizacion"] = d.get("depol")
        modos.append(m)
    datos = {
        "modos": modos,
        "espectro_simulado": {
            "laser_nm": float(laser_nm),
            "temperatura_K": float(T),
            "fwhm_cm-1": float(fwhm),
            "desplazamiento_cm-1": w,
            "intensidad_rel": inten,
        },
        "picos": [{"omega_cm-1": a, "intensidad": b} for a, b in picos],
    }
    return envelope(
        "raman_dft", datos, atoms,
        notas="Frecuencias armonicas en Gamma. La intensidad simulada "
              "aplica (wL-w)^4/w y el factor de Bose a la actividad "
              "calculada, que es lo que hace comparable el calculo con un "
              "espectro medido: las actividades crudas NO son intensidades. "
              "Sin anarmonicidad, asi que las frecuencias suelen salir "
              "algo altas frente a la medida.")


def from_xps(res, atoms=None) -> dict:
    """Corrimientos de nivel de core, para la app de XPS."""
    atomos = []
    for i, sh in enumerate(res.shifts):
        d = {"indice": i + 1, "shift_eV": float(sh)}
        if i < len(res.symbols):
            d["especie"] = res.symbols[i]
        for k, v in res.contributions.items():
            if k != "TOTAL":
                d[f"contrib_{k.lower().replace(' ', '_')}_eV"] = float(v[i])
        atomos.append(d)
    return envelope(
        "xps_core_dft", {"atomos": atomos,
                         "todos_equivalentes": bool(res.equivalentes)},
        atoms,
        notas="Aproximacion de ESTADO INICIAL: comparable en corrimientos "
              "RELATIVOS entre sitios, no en energias de enlace absolutas. "
              "La relajacion frente al hueco de core (estado final) no esta "
              "incluida y puede valer varias decimas de eV.")
