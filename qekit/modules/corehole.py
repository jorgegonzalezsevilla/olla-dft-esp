# Olla-DFT — command-line toolkit for Quantum ESPRESSO
# Copyright (C) 2026 Jorge Enrique González Sevilla
# SPDX-License-Identifier: AGPL-3.0-or-later
# This program is free software: you can redistribute it and/or modify it
# under the terms of the GNU Affero General Public License as published by the Free
# Software Foundation, either version 3 of the License, or (at your option)
# any later version. See the LICENSE file for details.

"""Pseudopotenciales con hueco de core: la pieza que faltaba para XPS y XANES.

POR QUÉ EXISTE ESTE MÓDULO
--------------------------
Dos módulos de Olla-DFT estaban bloqueados por la misma razón. `initial_state.x`
(XPS) necesita que el input declare DOS especies del mismo elemento —la
normal y una con un electrón menos en el core— y calcula el corrimiento a
partir de `zv(excitada) - zv(normal)`. `xspectra.x` (XANES) necesita lo
mismo, más la función de onda del nivel de core.

Esos pseudopotenciales no vienen en ninguna tabla estándar. Hay que
generarlos con `ld1.x`, y escribir a mano su configuración electrónica y
los radios de corte es exactamente la barrera que hace que la gente
abandone el cálculo. Este módulo lo hace: le das el elemento y el borde,
y salen los dos pseudos consistentes entre sí.

LO QUE HAY QUE ENTENDER ANTES DE USARLO
---------------------------------------
1. **Los dos pseudos tienen que generarse juntos.** Comparar una energía
   hecha con un pseudo de una tabla contra otra hecha con un pseudo de
   `ld1.x` no significa nada: la referencia de energía es distinta. Por
   eso `generar()` siempre escribe LOS DOS con los mismos parámetros.

2. **z_valence sube exactamente una unidad.** El silicio normal tiene
   z_valence = 4; con el hueco 1s, 5. Esa unidad ES el hueco. Si los dos
   pseudos no difieren en exactamente 1, algo se hizo mal y `verificar()`
   lo dice.

3. **Un pseudo recién generado no está validado.** Que `ld1.x` termine no
   quiere decir que el pseudo sirva: puede tener estados fantasma, puede
   necesitar un cutoff enorme, puede reproducir mal las derivadas
   logarítmicas. Olla-DFT revisa lo que se puede revisar automáticamente y
   dice claramente lo que NO revisó.

4. **El cutoff hay que volver a converger.** Un pseudo nuevo no hereda el
   ecutwfc del anterior. `olla-dft converge` sigue siendo obligatorio.

SOBRE LA FUNCIÓN DE ONDA DE CORE
--------------------------------
`xspectra.x` lee la función de onda del nivel de core de un archivo de
texto (`filecore`) con un formato propio: una línea de encabezado y luego
un BLOQUE por orbital de core, separados por líneas en blanco, en el
orden 1s, 2s, 2p, 3s, 3p, 3d.

Ese archivo NO es el `.wfc` que escribe `ld1.x` —ese sale en columnas, y
si se le pasa tal cual a `xspectra.x` el programa lee la columna
equivocada sin quejarse—. Sale del propio UPF, de las secciones
`PP_GIPAW_CORE_ORBITAL`, que solo están si el pseudo se generó con
`lgipaw_reconstruction=.true.`. Olla-DFT lo extrae aquí, en Python, y
verifica que el orbital pedido esté.
"""

import re
import shutil
import subprocess
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

from qekit.core import atomconf, provenance
from qekit.core.errors import ErrorDeUso

#: Radios de corte por defecto (bohr), por fila de la tabla periódica.
#: Son un punto de partida conservador, no un óptimo: un radio grande da
#: un pseudo blando pero transferible de forma dudosa, y uno chico da lo
#: contrario. Se sobrescriben con --rcut.
RCUT_FILA = {1: 1.0, 2: 1.3, 3: 1.7, 4: 2.0, 5: 2.2, 6: 2.4}

#: Energía de referencia (Ry) para pseudizar un canal desocupado.
E_CANAL_VACIO = 0.15

#: Energía (Ry) del SEGUNDO proyector de cada canal, cuando se piden dos.
E_SEGUNDO_PROYECTOR = 0.05

#: rcutus / rc para los pseudos ultrasuaves.
FACTOR_RCUTUS = 1.25

_RE_ATRIB = re.compile(r'(\w+)\s*=\s*"([^"]*)"')


def _fila(simbolo: str) -> int:
    z = atomconf.Z_DE[simbolo]
    for fila, limite in enumerate((2, 10, 18, 36, 54, 86), start=1):
        if z <= limite:
            return fila
    return 6


def rcut_sugerido(simbolo: str) -> float:
    return RCUT_FILA[_fila(simbolo)]


@dataclass
class Pseudo:
    """Un pseudopotencial generado, con lo que se pudo verificar de él."""
    ruta: str = ""
    elemento: str = ""
    z_valence: float = None
    mesh: int = None
    tipo: str = ""
    funcional: str = ""
    orbitales_core: list = field(default_factory=list)   # etiquetas GIPAW
    ecutwfc_sugerido: float = None
    ecutrho_sugerido: float = None


@dataclass
class Generacion:
    """El par de pseudos y todo lo que hubo que decidir para hacerlos."""
    elemento: str = ""
    borde: str = ""
    nivel: str = ""              # 1S, 2P, ...
    base: Pseudo = None
    hueco: Pseudo = None
    core: list = field(default_factory=list)
    valencia: list = field(default_factory=list)
    canales: list = field(default_factory=list)
    rcut: float = None
    proyectores: int = 1
    funcional: str = ""
    entradas: list = field(default_factory=list)   # los .in escritos
    avisos: list = field(default_factory=list)
    ok: bool = False


# ----------------------------------------------------------------------
# Lectura del UPF
# ----------------------------------------------------------------------

def _cabecera(texto: str) -> str:
    """El bloque PP_HEADER, venga cerrado o autocerrado.

    ld1.x de QE 6.6 lo escribe autocerrado (`<PP_HEADER .../>`) y otros
    conversores lo escriben con etiqueta de cierre. Buscar solo el cierre
    devuelve -1 y deja la cabecera en once caracteres: el pseudo parece
    ilegible cuando está perfecto.
    """
    m = re.search(r"<PP_HEADER\b(.*?)(?:/>|</PP_HEADER>)", texto, re.S)
    return m.group(1) if m else texto[:20000]

def leer_upf(ruta) -> Pseudo:
    """Lee del UPF lo que hace falta para saber si sirve."""
    ruta = Path(ruta)
    texto = ruta.read_text(errors="ignore")
    p = Pseudo(ruta=str(ruta))
    at = dict(_RE_ATRIB.findall(_cabecera(texto)))
    p.elemento = (at.get("element") or "").strip()
    for clave, destino in (("z_valence", "z_valence"),
                           ("mesh_size", "mesh"),
                           ("wfc_cutoff", "ecutwfc_sugerido"),
                           ("rho_cutoff", "ecutrho_sugerido")):
        val = at.get(clave)
        if val is None:
            continue
        try:
            num = float(val)
        except ValueError:
            continue
        setattr(p, destino, int(num) if destino == "mesh" else num)
    p.tipo = (at.get("pseudo_type") or "").strip()
    p.funcional = (at.get("functional") or "").strip()
    p.orbitales_core = re.findall(
        r'<PP_GIPAW_CORE_ORBITAL\.\d+[^>]*label="([^"]+)"', texto)
    if not p.elemento:
        m = re.search(r"Element:\s*(\w+)", texto)
        if m:
            p.elemento = m.group(1)
    for campo in ("ecutwfc_sugerido", "ecutrho_sugerido"):
        if getattr(p, campo) == 0.0:
            setattr(p, campo, None)
    return p


def core_wfc(upf, destino, orbital: str = None) -> str:
    """Extrae la función de onda de core al formato que lee `xspectra.x`.

    Equivale al `upf2plotcore.sh` de QE, pero en Python y comprobando lo
    que ese script no comprueba: que el UPF traiga orbitales GIPAW y que
    el número de puntos coincida con la malla declarada. Si el pseudo se
    generó sin `lgipaw_reconstruction`, aquí se detecta — en vez de que
    `xspectra.x` lea basura y devuelva un espectro plausible pero falso.
    """
    upf = Path(upf)
    texto = upf.read_text(errors="ignore")
    at = dict(_RE_ATRIB.findall(_cabecera(texto)))
    try:
        mesh = int(float(at["mesh_size"]))
    except (KeyError, ValueError):
        raise ErrorDeUso(
            f"{upf.name} no declara mesh_size; no parece un UPF v2.") from None

    m = re.search(r"<PP_R[\s>](.*?)</PP_R>", texto, re.S)
    if m is None:
        raise ErrorDeUso(f"{upf.name} no trae la malla radial <PP_R>.")
    r = np.fromstring(re.sub(r"<[^>]*>", " ", m.group(1)), sep=" ")

    bloques = re.findall(
        r'<PP_GIPAW_CORE_ORBITAL\.(\d+)([^>]*)>(.*?)</PP_GIPAW_CORE_ORBITAL\.\1>',
        texto, re.S)
    if not bloques:
        raise ErrorDeUso(
            f"{upf.name} no trae orbitales de core (PP_GIPAW_CORE_ORBITAL).\n"
            "Ese pseudopotencial NO sirve para XANES: hay que regenerarlo con "
            "lgipaw_reconstruction=.true.\n"
            "  olla-dft corehole <elemento> --edge K")

    etiquetas, ondas = [], []
    for _, atrib, cuerpo in bloques:
        a = dict(_RE_ATRIB.findall(atrib))
        etiquetas.append((a.get("label") or "?").strip().upper())
        ondas.append(np.fromstring(re.sub(r"<[^>]*>", " ", cuerpo), sep=" "))

    if orbital is not None and orbital.upper() not in etiquetas:
        raise ErrorDeUso(
            f"{upf.name} no trae el orbital {orbital}. "
            f"Tiene: {', '.join(etiquetas)}")

    n = min([mesh, len(r)] + [len(w) for w in ondas])
    destino = Path(destino)
    destino.parent.mkdir(parents=True, exist_ok=True)
    with open(destino, "w") as f:
        f.write(f"# funciones de onda de core de {upf.name}: "
                f"{len(ondas)} estados ({', '.join(etiquetas)})\n")
        for w in ondas:
            for i in range(n):
                f.write(f"{r[i]:20.12f} {w[i]:20.12f}\n")
            f.write("\n")
    return str(destino)


# ----------------------------------------------------------------------
# Escritura de los inputs de ld1.x
# ----------------------------------------------------------------------
def input_ld1(simbolo: str, config: list, canales: list, salida_upf: str,
              rcut: float, dft: str = "PBE", prefix: str = None,
              titulo: str = "", rel: int = 0, pseudotype: int = 2,
              beta: float = 0.3) -> str:
    """Arma un input de ld1.x para generar un pseudopotencial (iswitch=3)."""
    z = float(atomconf.Z_DE[simbolo])
    lineas = [" &input",
              f"    title='{titulo or simbolo}',",
              f"    prefix='{prefix or simbolo}',",
              f"    zed={z},",
              f"    rel={rel},",
              f"    beta={beta},",
              "    iswitch=3,",
              f"    dft='{dft}',",
              " /",
              f"{len(config)}"]
    for n, l, occ in config:
        lineas.append(f"{atomconf.etiqueta(n, l):<4s}{n:2d} {l:2d} "
                      f"{occ:6.2f}  1")
    # El canal local: el de l más alto, que en la práctica es el d vacío.
    lloc = max(l for _, _, l, _ in canales)
    vistos = set()
    lineas += [" &inputp",
               f"    pseudotype={pseudotype},",
               f"    lloc={lloc},",
               "    tm=.true.,",
               f"    file_pseudopw='{salida_upf}',",
               "    lgipaw_reconstruction=.true.,",
               "    author='Olla-DFT',",
               " /",
               f"{len(canales)}"]
    # El segundo número de cada línea es el n del PSEUDO-orbital, y los dos
    # proyectores del mismo canal comparten n. Numerarlos por línea (1,2,3,4)
    # en vez de por canal (1,1,2,2) hace que ld1.x falle con "Zero norm:
    # self consistency problem", que no dice nada de esto.
    indice, orden = {}, 0
    for etq, n, l, occ in canales:
        # Un canal DESOCUPADO no tiene autovalor ligado que buscar: hay que
        # decirle a ld1.x a qué energía pseudizarlo. Con la columna en 0.00
        # intenta tratarlo como estado ligado y se planta con "mismatched
        # all-electron/pseudo occupations", que no dice nada de esto.
        if occ < 0:
            # Estado no ligado: hay que decir a qué energía pseudizarlo.
            # El segundo proyector de un canal ya ocupado va cerca del
            # nivel; un canal enteramente vacío, algo más arriba.
            energia = E_SEGUNDO_PROYECTOR if l in vistos else E_CANAL_VACIO
        else:
            energia = 0.00
            vistos.add(l)
        # Cuarta columna = rcutus, el radio de la aumentación ultrasuave.
        # Solo se usa con pseudotype=3; con norma conservada los dos radios
        # son el mismo.
        rcus = rcut * FACTOR_RCUTUS if pseudotype == 3 else rcut
        if (etq, l) not in indice:
            orden += 1
            indice[(etq, l)] = orden
        lineas.append(f"{etq:<4s}{indice[(etq, l)]:2d} {l:2d} {occ:6.2f} "
                      f"{energia:5.2f}  {rcut:5.2f}  {rcus:5.2f}")
    return "\n".join(lineas) + "\n"


def _correr_ld1(entrada: Path, salida: Path, ld1_cmd: str = None) -> None:
    exe = ld1_cmd or shutil.which("ld1.x")
    if not exe:
        raise ErrorDeUso(
            "no se encontró ld1.x. Es parte de Quantum ESPRESSO pero no se "
            "compila por defecto:\n  cd <fuente de QE> && make ld1\n"
            "Después ponlo en el PATH, o pásalo con --ld1-cmd.")
    with open(entrada) as fin, open(salida, "w") as fout:
        proc = subprocess.run([exe], stdin=fin, stdout=fout,
                              stderr=subprocess.STDOUT,
                              cwd=str(entrada.parent))
    texto = salida.read_text(errors="ignore")
    if proc.returncode != 0 or "Error in routine" in texto:
        cola = [ln.strip() for ln in texto.splitlines() if ln.strip()][-8:]
        raise RuntimeError(
            f"ld1.x falló; revisa {salida}\n" + "\n".join("  " + c for c in cola))


# ----------------------------------------------------------------------
# Generación del par
# ----------------------------------------------------------------------

def _con_canales_vacios(config: list, canales: list) -> list:
    """Añade a la configuración de todos los electrones los canales vacíos.

    `ld1.x` exige que cada canal de la tarjeta de pseudización tenga su
    estado correspondiente en la configuración de todos los electrones; si
    falta, se planta con "no all electron for this ps". Los canales vacíos
    entran con ocupación negativa, que es como ld1.x marca un estado que
    hay que calcular sin ocupar.
    """
    presentes = {(n, l) for n, l, _ in config}
    fuera = list(config)
    for _, n, l, occ in canales:
        if occ <= 0 and (n, l) not in presentes:
            fuera.append((n, l, -1.0))
    return sorted(fuera, key=lambda t: (t[0], t[1]))

def generar(simbolo: str, borde: str = "K", outdir: str = "pseudos",
            dft: str = "PBE", rcut: float = None, rel: int = 0,
            semicore: bool = False, pseudotype: int = 2,
            proyectores: int = 1, solo_base: bool = False,
            correr: bool = True, ld1_cmd: str = None) -> Generacion:
    """Genera el par de pseudopotenciales: normal y con hueco de core.

    Los dos con los MISMOS parámetros, que es la única forma de que la
    diferencia entre ellos signifique algo.
    """
    simbolo = simbolo.strip().capitalize()
    if simbolo not in atomconf.Z_DE:
        raise ErrorDeUso(
            f"elemento '{simbolo}' fuera de la tabla que Olla-DFT conoce "
            "(H..Rn).")
    if proyectores >= 2 and pseudotype != 3:
        # Con norma conservada, dos proyectores en el mismo canal hacen que
        # ld1.x no resuelva la ecuación de Kohn-Sham pseudizada. El esquema
        # que sí admite dos por canal es el ultrasuave.
        pseudotype = 3
    g = Generacion(elemento=simbolo, borde=borde.upper(), funcional=dft)
    g.rcut = rcut if rcut is not None else rcut_sugerido(simbolo)
    g.core, g.valencia = atomconf.particion(simbolo, semicore=semicore)
    g.canales = atomconf.canales_pseudo(simbolo, semicore=semicore,
                                        proyectores=proyectores)
    g.proyectores = proyectores

    base_conf = _con_canales_vacios(atomconf.configuracion(simbolo), g.canales)
    hueco_conf = None
    if not solo_base:
        hueco_bruto, nivel = atomconf.config_hueco(simbolo, g.borde)
        hueco_conf = _con_canales_vacios(hueco_bruto, g.canales)
        g.nivel = nivel

    out = Path(outdir)
    out.mkdir(parents=True, exist_ok=True)
    upf_base = f"{simbolo}.qekit.UPF"
    upf_hueco = f"{simbolo}.hueco{(g.nivel or 'x').lower()}.UPF"

    pares = [("base", base_conf, upf_base, f"{simbolo} base")]
    if not solo_base:
        pares.append(("hueco", hueco_conf, upf_hueco,
                      f"{simbolo} con hueco de core en {g.nivel}"))
    for etiqueta, conf, nombre_upf, titulo in pares:
        texto = input_ld1(simbolo, conf, g.canales, nombre_upf, g.rcut,
                          dft=dft, prefix=f"{simbolo}_{etiqueta}",
                          titulo=titulo, rel=rel, pseudotype=pseudotype)
        entrada = out / f"ld1_{etiqueta}.in"
        entrada.write_text(texto)
        g.entradas.append(str(entrada))
        if correr:
            _correr_ld1(entrada, out / f"ld1_{etiqueta}.out", ld1_cmd)
            p = leer_upf(out / nombre_upf)
            setattr(g, etiqueta, p)

    if correr:
        g.avisos = verificar(g) if not solo_base else _verificar_base(g)
        g.ok = not any(a.startswith("FALLA") for a in g.avisos)
    return g


def _verificar_base(g: Generacion) -> list:
    """Comprobaciones cuando solo se pidió el pseudopotencial normal."""
    if g.base is None:
        return ["FALLA: no se generó el pseudopotencial."]
    return [f"z_valence = {g.base.z_valence:g}, malla de {g.base.mesh} "
            f"puntos, tipo {g.base.tipo}.",
            "Solo se genero el pseudopotencial NORMAL (--plain). Para XPS "
            "o XANES hace falta tambien el de hueco de core: quita "
            "--plain.",
            "NO verificado automaticamente: estados fantasma, derivadas "
            "logaritmicas y transferibilidad. Corre 'olla-dft converge' "
            "antes de usarlo en serio."]


def verificar(g: Generacion) -> list:
    """Lo que se puede comprobar solo, y lo que NO se comprobó."""
    avisos = []
    if g.base is None or g.hueco is None:
        return ["FALLA: falta uno de los dos pseudopotenciales."]

    dz = (g.hueco.z_valence or 0) - (g.base.z_valence or 0)
    if abs(dz - 1.0) > 1e-6:
        avisos.append(
            f"FALLA: z_valence pasa de {g.base.z_valence} a "
            f"{g.hueco.z_valence} (diferencia {dz:+g}); tiene que ser "
            "exactamente +1. Un hueco de core es UN electrón.")
    else:
        avisos.append(
            f"z_valence {g.base.z_valence:g} -> {g.hueco.z_valence:g}: "
            "la diferencia de +1 es el hueco. Correcto.")

    if g.base.mesh != g.hueco.mesh:
        avisos.append(
            f"FALLA: las mallas radiales no coinciden ({g.base.mesh} vs "
            f"{g.hueco.mesh}); xspectra.x interpola la función de core sobre "
            "la malla del pseudo y saldría corrida.")

    if g.nivel.upper() not in [o.upper() for o in g.hueco.orbitales_core]:
        avisos.append(
            f"FALLA: el UPF con hueco no trae el orbital {g.nivel} entre sus "
            "funciones de core; XANES no va a poder leerlo.")
    else:
        avisos.append(
            f"funciones de onda de core presentes "
            f"({', '.join(g.hueco.orbitales_core)}): sirve para XANES.")

    if (g.base.funcional or "").split() != (g.hueco.funcional or "").split():
        avisos.append("FALLA: los dos pseudos usan funcionales distintos.")

    if g.proyectores < 2:
        avisos.append(
            "UN proyector por canal. XSpectra recomienda dos: con uno solo "
            "el espectro XANES es fiable cerca del borde pero se degrada a "
            "partir de unos 10 eV por encima. Para XPS no importa. Con dos "
            "(--projectors 2) el pseudo sale ultrasuave y casi siempre hay "
            "que ajustar --rcut a mano hasta que ld1.x converja.")

    avisos.append(
        "NO verificado automáticamente: estados fantasma, derivadas "
        "logarítmicas y transferibilidad. Antes de publicar nada con este "
        "pseudopotencial, revisa ld1_base.out y corre 'olla-dft converge' — el "
        "cutoff del pseudo anterior NO sirve para este.")
    return avisos


def report(g: Generacion) -> str:
    lines = ["--- Pseudopotenciales con hueco de core ---",
             (f"Elemento: {g.elemento}   borde {g.borde} "
              f"(nivel {g.nivel})") if g.nivel else
             f"Elemento: {g.elemento}   (solo el pseudo normal)",
             f"Funcional: {g.funcional}   radio de corte: {g.rcut:.2f} bohr",
             ""]
    lines.append("Partición core / valencia (regla de Olla-DFT, "
                 "sobrescribible):")
    lines.append("  core     = " + " ".join(
        f"{atomconf.etiqueta(n, l)}{occ:g}" for n, l, occ in g.core))
    lines.append("  valencia = " + " ".join(
        f"{atomconf.etiqueta(n, l)}{occ:g}" for n, l, occ in g.valencia))
    lines.append("  canales pseudizados = " + " ".join(
        f"{e}({o:g})" if o >= 0 else f"{e}(vacío)" for e, _, _, o in g.canales))
    lines.append("")

    if g.base:
        lines.append(f"{'':16s}{'z_valence':>10s} {'malla':>7s} {'tipo':>6s}")
        for nombre, p in [("normal", g.base)] + \
                ([("con hueco", g.hueco)] if g.hueco else []):
            lines.append(f"  {nombre:14s}{p.z_valence:10.2f} {p.mesh:7d} "
                         f"{p.tipo:>6s}   {Path(p.ruta).name}")
        lines.append("")

    if g.avisos:
        lines.append("Verificación:")
        for a in g.avisos:
            marca = "  [FALLA] " if a.startswith("FALLA") else "  [ok]    "
            texto = a[7:] if a.startswith("FALLA: ") else a
            lines.append(marca + texto)
        lines.append("")

    lines += [
        "Cómo se usan:",
        "  XPS   -> olla-dft xps estructura.cif --core-hole "
        f"{g.elemento}={Path(g.hueco.ruta).name if g.hueco else '(falta)'}",
        f"  XANES -> olla-dft xanes estructura.cif --pseudo-dir {Path(g.entradas[0]).parent if g.entradas else '.'}",
        "",
        "Los DOS pseudos van juntos: el input declara las dos especies, la",
        "normal para todos los átomos y la del hueco solo como contraparte.",
    ]
    return "\n".join(lines)


def export(g: Generacion, outdir: str = ".") -> list:
    out = Path(outdir); out.mkdir(parents=True, exist_ok=True)
    f = out / "PSEUDOS_HUECO.txt"
    f.write_text(provenance.header_plain(
        "pseudopotenciales con hueco de core",
        {"elemento": g.elemento, "borde": g.borde, "funcional": g.funcional,
         "rcut_bohr": g.rcut},
        titulo="Generacion de pseudopotenciales") + "\n" + report(g) + "\n")
    return [str(f)]
