# Olla-DFT — command-line toolkit for Quantum ESPRESSO
# Copyright (C) 2026 Jorge Enrique González Sevilla
# SPDX-License-Identifier: GPL-3.0-or-later
# This program is free software: you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by the Free
# Software Foundation, either version 3 of the License, or (at your option)
# any later version. See the LICENSE file for details.

"""Configuraciones electrónicas atómicas para generar pseudopotenciales.

`ld1.x` no adivina nada: hay que darle la configuración electrónica de
todos los estados, cuáles quedan en el core y cuáles son de valencia, y
un radio de corte por canal. Escribir eso a mano para cada elemento es
justo la barrera que hace que la gente no genere pseudopotenciales.

Este módulo construye esa configuración por llenado de Aufbau, con las
excepciones conocidas, y aplica una regla explícita para separar core de
valencia. Ninguna de esas dos cosas es "la verdad": son un punto de
partida razonable que el usuario puede sobrescribir entero. Por eso todo
lo que se decide aquí se REPORTA, y no se esconde.
"""

from qekit.core.errors import ErrorDeUso

SIMBOLOS = (
    "H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar K Ca Sc Ti V Cr Mn Fe "
    "Co Ni Cu Zn Ga Ge As Se Br Kr Rb Sr Y Zr Nb Mo Tc Ru Rh Pd Ag Cd In "
    "Sn Sb Te I Xe Cs Ba La Ce Pr Nd Pm Sm Eu Gd Tb Dy Ho Er Tm Yb Lu Hf "
    "Ta W Re Os Ir Pt Au Hg Tl Pb Bi Po At Rn"
).split()

Z_DE = {s: i + 1 for i, s in enumerate(SIMBOLOS)}

#: Orden de llenado (Madelung), como (n, l).
ORDEN = [(1, 0), (2, 0), (2, 1), (3, 0), (3, 1), (4, 0), (3, 2), (4, 1),
         (5, 0), (4, 2), (5, 1), (6, 0), (4, 3), (5, 2), (6, 1), (7, 0)]

LETRA = "SPDF"

#: Excepciones al llenado de Aufbau, como {símbolo: {(n, l): ocupación}}.
#: Solo las que cambian el número de electrones de valencia; se aplican
#: sobre la configuración de Aufbau ya construida.
EXCEPCIONES = {
    "Cr": {(3, 2): 5, (4, 0): 1},
    "Cu": {(3, 2): 10, (4, 0): 1},
    "Nb": {(4, 2): 4, (5, 0): 1},
    "Mo": {(4, 2): 5, (5, 0): 1},
    "Ru": {(4, 2): 7, (5, 0): 1},
    "Rh": {(4, 2): 8, (5, 0): 1},
    "Pd": {(4, 2): 10, (5, 0): 0},
    "Ag": {(4, 2): 10, (5, 0): 1},
    "La": {(5, 2): 1, (4, 3): 0},
    "Ce": {(5, 2): 1, (4, 3): 1},
    "Gd": {(5, 2): 1, (4, 3): 7},
    "Pt": {(5, 2): 9, (6, 0): 1},
    "Au": {(5, 2): 10, (6, 0): 1},
}

#: Bordes de absorción / niveles de core, como {nombre: (n, l)}.
BORDES = {
    "K":   (1, 0),
    "L1":  (2, 0),
    "L23": (2, 1),
    "M1":  (3, 0),
    "M23": (3, 1),
    "M45": (3, 2),
}


def _capacidad(l: int) -> int:
    return 2 * (2 * l + 1)


def etiqueta(n: int, l: int) -> str:
    return f"{n}{LETRA[l]}"


def aufbau(z: int) -> list:
    """Configuración electrónica de Aufbau: lista de (n, l, ocupación)."""
    restantes = z
    conf = []
    for n, l in ORDEN:
        if restantes <= 0:
            break
        occ = min(_capacidad(l), restantes)
        conf.append((n, l, float(occ)))
        restantes -= occ
    if restantes > 0:
        raise ErrorDeUso(f"Z = {z} queda fuera de la tabla de llenado.")
    return conf


def configuracion(simbolo: str) -> list:
    """Configuración del elemento, con las excepciones conocidas aplicadas."""
    if simbolo not in Z_DE:
        raise ErrorDeUso(
            f"elemento '{simbolo}' fuera de la tabla (H..Rn). "
            "Da la configuración a mano con --config.")
    conf = aufbau(Z_DE[simbolo])
    for (n, l), occ in EXCEPCIONES.get(simbolo, {}).items():
        for i, (nn, ll, _) in enumerate(conf):
            if (nn, ll) == (n, l):
                conf[i] = (nn, ll, float(occ))
                break
        else:
            conf.append((n, l, float(occ)))
    # ordenar por n y luego por l: es como ld1.x espera leerlos
    return sorted(conf, key=lambda t: (t[0], t[1]))


def particion(simbolo: str, semicore: bool = False) -> tuple:
    """Separa la configuración en (core, valencia).

    REGLA, dicha en voz alta para que se pueda discutir:

    - valencia = la capa s y p de n máximo, más cualquier d o f
      parcialmente llena (una d llena a la mitad es química, no core);
    - con `semicore=True` entra además la capa (n-1)s(n-1)p, que en
      metales alcalinos y alcalinotérreos hace falta de verdad;
    - todo lo demás va al core.

    No hay una respuesta universal: la partición correcta depende del
    sistema. Por eso Olla-DFT la reporta y se puede sobrescribir.
    """
    conf = configuracion(simbolo)
    nmax = max(n for n, _, _ in conf)
    valencia, core = [], []
    for n, l, occ in conf:
        es_val = False
        if n == nmax and l <= 1:
            es_val = True
        elif l >= 2 and 0 < occ < _capacidad(l):
            es_val = True                      # d o f parcialmente llena
        elif l == 2 and n == nmax - 1 and occ == _capacidad(l):
            es_val = True                      # d llena de la fila anterior
        elif semicore and n == nmax - 1 and l <= 1:
            es_val = True
        (valencia if es_val else core).append((n, l, occ))
    return core, valencia


def canales_pseudo(simbolo: str, semicore: bool = False,
                   canal_vacio: bool = True, proyectores: int = 1) -> list:
    """Canales para la tarjeta de pseudización, como (etiqueta, n, l, occ).

    Se añade un canal DESOCUPADO de l más alto (típicamente el d) porque
    sin él el potencial local queda mal descrito y el pseudo produce
    estados fantasma. Se marca con ocupación negativa, que es como ld1.x
    entiende "calcúlalo pero no lo ocupes".
    """
    _, valencia = particion(simbolo, semicore=semicore)
    if not valencia:
        raise ErrorDeUso(f"no se pudo separar la valencia de {simbolo}.")
    nmax = max(n for n, _, _ in valencia)
    canales = [(etiqueta(n, l), n, l, occ) for n, l, occ in valencia]
    if canal_vacio:
        # Un canal por cada l hasta d. Sin el p, el sodio y los metales de
        # transición quedan con un potencial local malo; sin el d aparecen
        # estados fantasma. Los que faltan entran DESOCUPADOS: ocupación
        # negativa es como ld1.x entiende "calcúlalo pero no lo llenes".
        presentes = {l for _, l, _ in valencia}
        for l in (0, 1, 2):
            if l in presentes:
                continue
            # el número cuántico principal más bajo que admite ese l:
            # no existe un 2d, así que para el oxígeno el canal vacío es 3d
            n = max(nmax, l + 1)
            canales.append((etiqueta(n, l), n, l, -2.0))
    canales = sorted(canales, key=lambda t: (t[2], t[1]))
    if proyectores >= 2:
        # Un segundo proyector por canal, desocupado y a energia mas alta.
        # XSpectra lo pide explicitamente: con uno solo el espectro deja de
        # ser fiable a partir de unos 10 eV del borde. En ld1.x se hace
        # repitiendo el canal con ocupacion 0 y otra energia de referencia.
        # El segundo lleva n+1 en la etiqueta. No es cosmetico: ld1.x
        # guarda los orbitales de reconstruccion GIPAW indexados por
        # etiqueta, asi que dos proyectores con el mismo nombre colapsan en
        # uno y el pseudo acaba con un solo canal de reconstruccion — que es
        # justo lo que XSpectra desaconseja.
        dobles = []
        for etq, n, l, occ in canales:
            dobles.append((etq, n, l, occ))
            if occ >= 0:
                dobles.append((etiqueta(n + 1, l), n + 1, l, -1.0))
        canales = dobles
    return canales


def config_hueco(simbolo: str, borde: str = "K") -> tuple:
    """Configuración con un hueco en el nivel de core del borde pedido.

    Devuelve (configuración, etiqueta del nivel). El hueco es de UN
    electrón: es lo que corresponde a arrancar uno solo, y es lo que hace
    que z_valence del pseudo resultante sea exactamente una unidad mayor
    que la del pseudo normal.
    """
    if borde not in BORDES:
        raise ErrorDeUso(
            f"borde '{borde}' desconocido. Opciones: {', '.join(BORDES)}")
    n_h, l_h = BORDES[borde]
    conf = configuracion(simbolo)
    fuera = []
    for n, l, occ in conf:
        if (n, l) == (n_h, l_h):
            if occ < 1:
                raise ErrorDeUso(
                    f"{simbolo} no tiene electrones en {etiqueta(n_h, l_h)}: "
                    f"el borde {borde} no existe para este elemento.")
            fuera.append((n, l, occ - 1.0))
        else:
            fuera.append((n, l, occ))
    core, _ = particion(simbolo)
    if (n_h, l_h) not in [(n, l) for n, l, _ in core]:
        raise ErrorDeUso(
            f"{etiqueta(n_h, l_h)} no queda en el core de {simbolo} con la "
            "partición actual: un hueco en un nivel de VALENCIA no es un "
            "hueco de core y el pseudo resultante no sirve para XPS ni XANES.")
    return fuera, etiqueta(n_h, l_h)
