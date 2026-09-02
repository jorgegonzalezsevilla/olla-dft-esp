# Olla-DFT — command-line toolkit for Quantum ESPRESSO
# Copyright (C) 2026 Jorge Enrique González Sevilla
# SPDX-License-Identifier: GPL-3.0-or-later
# This program is free software: you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by the Free
# Software Foundation, either version 3 of the License, or (at your option)
# any later version. See the LICENSE file for details.

"""Recetas: sesiones completas, de la estructura al resultado.

POR QUÉ EXISTE
--------------
`olla-dft wizard` responde a "¿qué quiero saber?" y devuelve un plan. La
referencia (`olla-dft docs`) explica cada comando por separado. Falta lo de en
medio, que es lo que de verdad cuesta aprender: **cómo encajan unos con
otros**.

Un comando aislado no enseña nada. Lo que enseña es ver que el cutoff que
sale de `converge` es el que entra en `gen`, que el `.dat` que escribe
`elastic` es el que lee `derived`, y que `crosscheck` no sirve de nada hasta
que hay DOS módulos que hayan dejado resultados. Eso no está en la ayuda de
ningún comando porque no pertenece a ninguno.

Cada receta es una sesión entera con:

  - los comandos exactos, en orden;
  - **qué archivo deja cada paso y qué paso posterior lo lee** — es la parte
    que convierte una lista de comandos en un flujo de trabajo;
  - un extracto de la salida REAL, cuando la receta se ha corrido entera;
  - y el error típico de cada paso, que casi siempre no es del comando sino
    de la costura entre dos.

QUÉ IMPIDE QUE ESTO SE PUDRA
----------------------------
Una sección de ejemplos escrita a mano se queda obsoleta a la tercera
versión y nadie se entera. Aquí no: hay una prueba que recorre TODOS los
comandos de TODAS las recetas y comprueba, contra el propio árbol de
argparse, que el subcomando existe y que cada bandera que se usa existe en
él. Si alguien renombra una opción, la receta falla en pytest antes de que
la vea un usuario.

EN INGLÉS
---------
Las recetas se escriben UNA vez, en español. Con ``--language en`` cada
función recibe el idioma y traduce la estructura al vuelo con la tabla de
``qekit/data/i18n/recipes_en.json`` (``i18n.translate_data``): los comandos,
los nombres de archivo y los extractos de salida real no están en la tabla
y por eso quedan como están. Una prueba comprueba que la tabla cubre todas
las cadenas de todas las recetas, así que añadir una receta sin su
traducción falla en pytest.
"""

from dataclasses import dataclass, field
from functools import lru_cache
import re

from qekit.core import i18n
from qekit.core.errors import ErrorDeUso


@dataclass
class Paso:
    """Un comando dentro de una receta, con lo que produce y lo que consume."""
    comando: str
    hace: str
    escribe: list = field(default_factory=list)
    lee: list = field(default_factory=list)   # (texto, de qué paso viene)
    salida: str = ""                          # extracto REAL, o vacío
    ojo: str = ""                             # el error típico de este paso
    corre_qe: bool = False                    # este paso lanza pw.x y tarda


@dataclass
class Receta:
    clave: str
    titulo: str
    pregunta: str                # cómo lo diría alguien que empieza
    para_que: str                # qué tienes al final y qué NO tienes
    coste: str                   # minutos | horas | días
    pasos: list = field(default_factory=list)
    despues: list = field(default_factory=list)   # (texto, comando)
    verificada: bool = False     # los extractos vienen de una corrida real
    ver_tambien: list = field(default_factory=list)
    # palabras con las que alguien buscaría esto SIN saber la jerga. No es
    # decoración: la búsqueda por lenguaje llano solo funciona si las
    # palabras que usa la gente están escritas en algún sitio, y el título
    # técnico casi nunca las lleva.
    palabras: list = field(default_factory=list)


RECETAS = [
    Receta(
        clave="primero",
        palabras=["empezar", "instalar", "nuevo", "principio", "primera", "básico", "relajar", "convergencia"],
        titulo="La primera sesión",
        pregunta="acabo de instalarlo y no sé por dónde empezar",
        para_que="Al final tienes una estructura relajada con parámetros "
                 "convergidos, que es el punto de partida de TODO lo demás. "
                 "Todavía no tienes ninguna propiedad: esto es la base.",
        coste="minutos de tu tiempo, ~20 min de máquina",
        verificada=True,
        pasos=[
            Paso("olla-dft info Si.cif",
                 "Mira la estructura antes de calcular nada: grupo espacial, "
                 "sitios inequivalentes, distancias.",
                 escribe=["nada, solo imprime"],
                 salida="Fórmula: Si2   (2 átomos)\n"
                        "Grupo espacial: Fd-3m (N.º 227)\n"
                        "Elementos: Si",
                 ojo="Si el grupo espacial no es el que esperabas, para "
                     "aquí. Una celda mal centrada no da error: da "
                     "resultados equivocados durante horas."),
            Paso("olla-dft pseudos Si.cif --task scf",
                 "Elige el pseudopotencial por TAREA y explica cada "
                 "descarte.",
                 lee=[("la lista de elementos", "①")],
                 escribe=["nada, solo imprime"],
                 ojo="Elegir por orden alfabético es cómo acaba uno con el "
                     "Ni de PBE y el O de BLYP en el mismo cálculo."),
            Paso("olla-dft converge Si.cif --run -o conv",
                 "Barre cutoff y malla k hasta que la energía deje de "
                 "moverse.",
                 lee=[("los pseudos elegidos", "②")],
                 escribe=["conv/CONVERGENCIA.dat", "conv/convergencia.pdf"],
                 corre_qe=True,
                 ojo="No hay un cutoff universal: depende del "
                     "pseudopotencial. El número que salga de aquí es el que "
                     "usarás en TODOS los cálculos siguientes, y tienen que "
                     "ser el mismo o las energías no se pueden restar."),
            Paso("olla-dft gen Si.cif -p vc-relax -o relax --ecutwfc 30",
                 "Relaja posiciones y celda con TU funcional.",
                 lee=[("el cutoff convergido", "③")],
                 escribe=["relax/pw.in", "relax/run.sh"],
                 corre_qe=True,
                 ojo="Una estructura de una base de datos NO está en el "
                     "mínimo de tu funcional. Esa tensión residual contamina "
                     "las elásticas, los fonones y todo lo demás."),
        ],
        despues=[("saber si es metal o semiconductor", "olla-dft recetas bandas"),
                 ("ver si es estable y cuánto cuesta deformarlo",
                  "olla-dft recetas mecanicas"),
                 ("si no sabes qué quieres, pregúntaselo al asistente",
                  "olla-dft wizard Si.cif --ask 'quiero saber si absorbe luz'")],
    ),

    Receta(
        clave="bandas",
        palabras=["gap", "metal", "semiconductor", "aislante", "bandas", "electrónica", "conduce", "electricidad", "DOS"],
        titulo="Estructura de bandas, DOS y gap",
        pregunta="¿es metal o semiconductor? ¿cuánto vale su gap?",
        para_que="Al final tienes la figura de bandas + DOS lista para un "
                 "artículo y el gap con su carácter (directo o indirecto). "
                 "OJO: el gap de LDA/PBE es un 30-50 % menor que el "
                 "experimental. Eso NO es un fallo del cálculo.",
        coste="~15 min de máquina en un sistema pequeño",
        verificada=True,
        pasos=[
            Paso("olla-dft gen Si.cif -p all -o bandas --ecutwfc 30",
                 "Escribe de una vez los tres inputs encadenados: scf, nscf "
                 "y bands, con el camino de alta simetría ya puesto.",
                 lee=[("la estructura relajada y el cutoff",
                       "receta «primero»")],
                 escribe=["bandas/scf.in", "bandas/nscf.in",
                          "bandas/bands.in", "bandas/run.sh"],
                 ojo="El camino de alta simetría se calcula sobre la celda "
                     "primitiva estándar. Si tu celda es otra, las etiquetas "
                     "Γ, X, L no señalan lo que dicen; `olla-dft prim` la "
                     "estandariza."),
            Paso("bash bandas/run.sh",
                 "Corre los tres pasos en orden.",
                 lee=[("los inputs", "①")],
                 escribe=["bandas/out/*.xml", "las salidas de pw.x"],
                 corre_qe=True),
            Paso("olla-dft bands bandas -o figuras",
                 "Lee las bandas, detecta el gap y dibuja.",
                 lee=[("el XML del cálculo de bandas", "②")],
                 escribe=["figuras/BANDAS.dat", "figuras/bandas.pdf"],
                 salida="Gap indirecto: 0.524 eV   (Γ → punto sobre Γ-X)\n"
                        "Gap directo:   2.556 eV   en Γ",
                 ojo="Un gap de 0 eV con ocupaciones fijas no significa "
                     "«metal»: significa que el cálculo no era válido. "
                     "Los metales necesitan smearing."),
            Paso("olla-dft dos bandas -o figuras",
                 "Densidad de estados y proyección por orbital.",
                 lee=[("el nscf de malla densa", "②")],
                 escribe=["figuras/DOS.dat", "figuras/dos.pdf"],
                 ojo="Integrar la DOS hasta E_F tiene que devolver el número "
                     "de electrones de valencia. Si no, la malla del nscf es "
                     "demasiado pobre."),
            Paso("olla-dft plot bandas -o figuras",
                 "Bandas y DOS en una sola figura, con el eje de energía "
                 "compartido. Es para lo que existe este comando: no hay "
                 "bandera que activar.",
                 lee=[("los dos .dat anteriores", "③ y ④")],
                 escribe=["figuras/bandas_dos.pdf"]),
        ],
        despues=[("masas efectivas para transporte",
                  "olla-dft effmass bandas -o masas"),
                 ("óptica y color", "olla-dft optics ..."),
                 ("un modelo pequeño que interpole las bandas gratis",
                  "olla-dft recetas modelo")],
        ver_tambien=["primero", "termoelectrico"],
    ),

    Receta(
        clave="mecanicas",
        palabras=["duro", "dureza", "elástico", "estable", "rígido", "módulo", "Debye", "sonido", "presión", "mecánicas"],
        titulo="Estabilidad, dureza y lo que se deduce de ellas",
        pregunta="¿es estable? ¿qué tan duro es? ¿a qué temperatura de "
                 "Debye corresponde?",
        para_que="Al final tienes las constantes elásticas con el criterio "
                 "de Born, los módulos, las velocidades del sonido, la "
                 "temperatura de Debye y una estimación de la conductividad "
                 "térmica — todo desde el mismo cálculo.",
        coste="~1 h de máquina",
        verificada=True,
        pasos=[
            Paso("olla-dft eos Si.cif --run -o eos",
                 "Curva energía-volumen y su ajuste: parámetro de red de "
                 "equilibrio y módulo volumétrico.",
                 escribe=["eos/EOS.txt", "eos/eos.pdf"],
                 corre_qe=True,
                 salida="a0 = 5.402 Å      B0 = 94.2 GPa      B' = 4.2",
                 ojo="Hacen falta al menos 5 puntos y a AMBOS lados del "
                     "mínimo. Un ajuste con todos los puntos a un lado da un "
                     "B0 que parece razonable y está mal."),
            Paso("olla-dft elastic Si.cif --run -o elastic",
                 "Las 6×6 constantes elásticas por deformaciones finitas.",
                 escribe=["elastic/ELASTIC_C.dat"],
                 corre_qe=True,
                 salida="C11 = 159.9   C12 = 61.7   C44 = 76.6 GPa\n"
                        "Estable según el criterio de Born: sí",
                 ojo="La celda tiene que estar relajada ANTES. Un esfuerzo "
                     "residual entra directamente en las Cij."),
            Paso("olla-dft derived Si.cif --cij elastic/ELASTIC_C.dat -o derivadas",
                 "Post-proceso puro: velocidades del sonido, Debye, "
                 "Grüneisen y Slack. No corre nada nuevo.",
                 lee=[("las constantes elásticas", "②")],
                 escribe=["derivadas/DERIVED.dat"],
                 salida="Velocidades del sonido: longitudinal 8788 m/s  |"
                        "  transversal 5246 m/s\n"
                        "Temperatura de Debye (elástica): 637 K   "
                        "(experimental 645 K)",
                 ojo="La Debye elástica y la que sale de la DOS de fonones "
                     "NO son la misma definición. Que difieran un 20 % es "
                     "normal; que coincidan al 1 % sería sospechoso."),
            Paso("olla-dft crosscheck . -f Si.cif",
                 "Compara el B0 de la EOS con el que sale de la traza de las "
                 "Cij: dos rutas independientes para la misma cantidad.",
                 lee=[("EOS.txt y ELASTIC_C.dat", "① y ②")],
                 salida="[OK  ] módulo volumétrico B₀  (0.2 % de desvío)\n"
                        "         ajuste de la ecuación de estado: 94.2 GPa\n"
                        "         traza de las constantes elásticas: 94.4 GPa",
                 ojo="Este comando no sirve de nada hasta que hay DOS "
                     "módulos que hayan dejado resultados. Es la idea "
                     "entera: un número solo no se puede validar."),
        ],
        despues=[("fonones, que dan la otra Debye",
                  "olla-dft recetas vibra"),
                 ("conductividad térmica de verdad, no la de Slack",
                  "olla-dft recetas termoelectrico")],
        ver_tambien=["vibra"],
    ),

    Receta(
        clave="vibra",
        palabras=["fonón", "fonones", "vibración", "Raman", "infrarrojo", "calor específico", "termodinámica", "entropía"],
        titulo="Fonones, termodinámica y espectro Raman",
        pregunta="¿es dinámicamente estable? ¿cuánto vale su calor "
                 "específico? ¿dónde salen sus picos Raman?",
        para_que="Al final tienes la dispersión de fonones, la DOS, C_v(T), "
                 "la energía de punto cero y, si lo pides, las "
                 "intensidades Raman con su despolarización.",
        coste="horas: es el cálculo más caro de esta lista",
        verificada=True,
        pasos=[
            Paso("olla-dft phonons Si.cif -o fon --qgrid 4x4x4",
                 "Escribe la cadena entera: scf muy convergido → ph.x → q2r "
                 "→ matdyn (dispersión y DOS).",
                 escribe=["fon/1_scf.in", "fon/2_ph.in", "fon/3_q2r.in",
                          "fon/4_matdyn.in", "fon/5_dos.in"],
                 ojo="El scf de aquí necesita conv_thr mucho más apretado "
                     "que uno normal. Con 1e-6 los fonones salen con "
                     "frecuencias imaginarias que NO son inestabilidad: son "
                     "ruido."),
            Paso("bash fon/correr.sh",
                 "Los cinco pasos en orden. Es reanudable: si se corta, "
                 "vuelve a lanzarlo y sigue donde estaba.",
                 lee=[("los inputs", "①")],
                 corre_qe=True),
            Paso("olla-dft phonons Si.cif --collect -o fon",
                 "Dispersión, DOS y termodinámica.",
                 lee=[("las salidas de matdyn", "②")],
                 escribe=["fon/FONONES_BANDAS.dat", "fon/FONONES_DOS.dat",
                          "fon/FONONES_TERMO.dat", "fon/fonones.pdf"],
                 salida="Γ: 0.0  0.0  0.0  507.3  507.3  507.3 cm⁻¹\n"
                        "(experimental 520.7; sin frecuencias imaginarias)",
                 ojo="Frecuencias imaginarias en Γ que no sean las tres "
                     "acústicas = la estructura no está relajada, o la malla "
                     "de q es demasiado pequeña. Casi nunca es física."),
            Paso("olla-dft crosscheck fon -f Si.cif",
                 "La termodinámica ya la dejó el paso anterior en "
                 "FONONES_TERMO.dat. Lo que falta es comprobarla: aquí se "
                 "cruza C_v a alta temperatura contra Dulong-Petit y la "
                 "integral de la DOS contra 3N.",
                 lee=[("FONONES_DOS.dat y FONONES_TERMO.dat", "③")],
                 salida="[OK  ] C_v en el límite clásico  (0.4 % de desvío)\n"
                        "[OK  ] número de modos: 6.000 contra 3N = 6",
                 ojo="A alta temperatura C_v tiene que tender a 3N·k_B. Si "
                     "no llega, la DOS está mal normalizada o le falta "
                     "espectro; si se pasa, hay modos de más."),
        ],
        despues=[("expansión térmica: repite los fonones a varios volúmenes",
                  "olla-dft qha ..."),
                 ("conductividad térmica de red, que necesita el tercer "
                  "orden", "olla-dft recetas termoelectrico")],
        ver_tambien=["mecanicas", "termoelectrico"],
    ),

    Receta(
        clave="termoelectrico",
        palabras=["calor", "conduce el calor", "conductividad térmica", "termoeléctrico", "Seebeck", "ZT", "figura de mérito", "kappa", "térmica"],
        titulo="Un termoeléctrico completo: las tres piezas de ZT",
        pregunta="¿sirve como termoeléctrico? ¿cuánto vale su ZT?",
        para_que="ZT = S²σT/(κ_e + κ_L). Esta receta calcula las CUATRO "
                 "cantidades por caminos distintos y las junta. Es la que "
                 "mejor enseña cómo se relacionan los módulos: ninguno de "
                 "ellos da ZT solo.",
        coste="días si lo haces todo con DFT; horas si usas el potencial "
              "aprendido para explorar primero",
        verificada=False,
        pasos=[
            Paso("olla-dft gen Si.cif -p nscf -o denso --kspacing 0.15",
                 "El transporte necesita una malla de k MUCHO más densa que "
                 "unas bandas: se integran derivadas de la ocupación, no "
                 "energías.",
                 escribe=["denso/nscf.in"],
                 corre_qe=True,
                 ojo="Con una malla de bandas normal el Seebeck sale con "
                     "ruido de decenas de µV/K. Si S cambia al densificar, "
                     "no estaba convergido."),
            Paso("olla-dft transport denso --collect -o trans",
                 "σ/τ, Seebeck y κ_e/τ en la aproximación de tiempo de "
                 "relajación constante.",
                 lee=[("el nscf denso", "①")],
                 escribe=["trans/TRANSPORTE.dat"],
                 ojo="S y el número de Lorenz NO dependen de τ; σ y κ_e sí. "
                     "Publicar σ de una CRTA sin decir qué τ usaste es "
                     "publicar un número sin unidades."),
            Paso("olla-dft elph Si.cif -o elph --qgrid 2x2x2",
                 "Acoplamiento electrón-fonón: de aquí sale el τ de verdad, "
                 "el que le falta al paso anterior.",
                 escribe=["elph/1_scf.in", "elph/2_nscf.in", "elph/3_ph.in"],
                 corre_qe=True,
                 ojo="Son TRES pasos y el segundo es el que se olvida: sin "
                     "el nscf con la2F=.true. no hay archivo a2Fsave y ph.x "
                     "se muere sin decir por qué."),
            Paso("olla-dft kappa Si.cif --dim 2x2x2 --model mace",
                 "Primero explora barato: κ_L con un potencial aprendido "
                 "para elegir el tamaño de la supercelda.",
                 escribe=["kappa/KAPPA.dat", "kappa/KAPPA_recorrido.dat"],
                 salida="κ_L(300 K) = 50.8 W/m·K   (MACE, NO es DFT)\n"
                        "la mitad de κ la llevan fonones con Λ < 745 nm",
                 ojo="El valor absoluto de un potencial aprendido puede "
                     "fallar por un factor 2-3. La forma de κ(T) y la "
                     "convergencia con la supercelda sí salen bien, y para "
                     "eso es para lo que sirve."),
            Paso("olla-dft kappa Si.cif --dim 2x2x2 -o kL",
                 "Y ahora el de verdad: las mismas configuraciones con pw.x.",
                 lee=[("el tamaño de supercelda que eligió ④", "④")],
                 escribe=["kL/fc3/d0000/pw.in ... (57 carpetas)"],
                 corre_qe=True),
            Paso("olla-dft kappa Si.cif --collect -o kL --mesh 19 --isotopes",
                 "Resuelve la ecuación de Boltzmann de fonones.",
                 lee=[("las fuerzas de las 57 configuraciones", "⑤")],
                 escribe=["kL/KAPPA.dat"],
                 salida="κ_L(300 K) = 100.7 W/m·K   (96 con isótopos)\n"
                        "κ ∝ T^−1.16   ← el T⁻¹ de Umklapp",
                 ojo="κ tiene que converger en el tamaño de la supercelda Y "
                     "en la malla de q A LA VEZ."),
            Paso("olla-dft crosscheck . -f Si.cif",
                 "Cruza la κ del tercer orden contra la estimación de Slack "
                 "que salió de las elásticas. Dos rutas independientes.",
                 lee=[("KAPPA.dat y ELASTIC_C.dat", "⑥ y la receta "
                       "«mecanicas»")],
                 salida="[OK  ] conductividad térmica de red  (10.6 %)\n"
                        "         Boltzmann con fc3: 100.7 W/m/K\n"
                        "         modelo de Slack:    90.1 W/m/K",
                 ojo="Con la κ del potencial aprendido este cruce FALLA "
                     "(87 % de desvío). Eso no es un fallo del cruce: es el "
                     "cruce haciendo su trabajo."),
        ],
        despues=[("la ficha con todo junto", "olla-dft datasheet . -o ficha"),
                 ("si κ_L sale alta, mira el recorrido libre medio: dice si "
                  "nanoestructurar sirve", "kL/KAPPA_recorrido.dat")],
        ver_tambien=["mecanicas", "vibra", "fiarme"],
    ),

    Receta(
        clave="superficie",
        palabras=["superficie", "losa", "cara", "adsorción", "catalizador", "función trabajo", "electrodo", "corte", "interfaz"],
        titulo="De cristal a superficie, y de superficie a catalizador",
        pregunta="quiero cortar una cara, saber cuánto cuesta hacerla, qué "
                 "función trabajo tiene y si algo se adsorbe encima",
        para_que="Al final tienes γ por el ajuste de Fiorentini–Methfessel, "
                 "la función "
                 "trabajo sin ajustar mesetas, los sitios de adsorción no "
                 "equivalentes con su energía, y el diagrama de energía "
                 "libre de la reacción.",
        coste="~2 h de máquina",
        verificada=True,
        pasos=[
            Paso("olla-dft surface Al.cif --miller 1 1 1 --layers 5 "
                 "--vacuum 14 -o losa",
                 "Corta la cara (111) sobre la celda convencional, con vacío "
                 "y avisando si la superficie es polar.",
                 escribe=["losa/slab.cif"],
                 ojo="Una superficie polar no converge a nada útil sin "
                     "compensar el dipolo. El comando avisa; hacerle caso es "
                     "cosa tuya."),
            Paso("olla-dft gamma Al.cif --miller 1 1 1 --layers 4,5,6,7 "
                 "--run -o gamma",
                 "Energía de superficie por el ajuste lineal de "
                 "Fiorentini–Methfessel sobre varios grosores, no por "
                 "restar una energía de bulto importada.",
                 escribe=["gamma/GAMMA.dat"],
                 corre_qe=True,
                 salida="γ = 1.10 J/m²   (R² del ajuste = 1.000000)",
                 ojo="γ calculada con una E_bulto de otro cálculo NO "
                     "converge al engrosar la losa: deriva linealmente. Por "
                     "eso hacen falta varios grosores."),
            Paso("olla-dft esm losa/slab.cif --bc bc1 --run -o wf",
                 "Función trabajo con la condición de contorno de ESM: el "
                 "nivel de vacío vale cero por construcción, así que Φ = −E_F "
                 "sin ajustar ninguna meseta.",
                 lee=[("la losa cortada", "①")],
                 escribe=["wf/ESM.dat"],
                 corre_qe=True,
                 salida="Φ = 4.24 eV   (experimental 4.24-4.26)",
                 ojo="ESM mide z desde el CENTRO de la celda. Olla-DFT centra "
                     "la losa solo; si la pones tú a mano en c/2, cae sobre "
                     "la frontera y salen cientos de Ry de error sin "
                     "mensaje."),
            Paso("olla-dft adsorb losa/slab.cif --mol H --run -o sitios",
                 "Encuentra los sitios no equivalentes (top, puente, hueco) "
                 "y calcula la energía de adsorción en cada uno.",
                 lee=[("la losa", "①")],
                 escribe=["sitios/ADSORCION.dat"],
                 corre_qe=True,
                 ojo="Los sitios se buscan por huella de vecinos, no por "
                     "radio: con radio salían tres huecos donde solo hay "
                     "dos."),
            Paso("olla-dft echem --her -0.33 -o her",
                 "Del ΔE de adsorción al diagrama de energía libre: "
                 "potencial limitante y sobrepotencial.",
                 lee=[("la energía de adsorción del mejor sitio", "④")],
                 escribe=["her/ECHEM.dat", "her/echem.pdf"],
                 ojo="Sin las correcciones térmicas (ZPE − TΔS) el ΔG está "
                     "mal en décimas de eV, que es justo el tamaño de las "
                     "diferencias que se discuten."),
        ],
        despues=[("superficie cargada, a un potencial dado",
                  "olla-dft esm losa/slab.cif --bc bc3 --charge -0.04,0,0.04"),
                 ("barrera de la reacción", "olla-dft neb ...")],
        ver_tambien=["primero"],
    ),

    Receta(
        clave="defecto",
        palabras=["defecto", "vacante", "impureza", "dopante", "sustitución", "intersticial", "carga", "nivel", "trampa"],
        titulo="Un defecto cargado bien hecho",
        pregunta="¿cuánto cuesta hacer una vacante? ¿en qué estado de carga "
                 "está según dónde caiga el nivel de Fermi?",
        para_que="Al final tienes E_f(E_F) para cada estado de carga, la "
                 "envolvente y los niveles de transición ε(q/q'), que es lo "
                 "que se compara con una medida de DLTS.",
        coste="~3 h de máquina",
        verificada=True,
        pasos=[
            Paso("olla-dft eform Si.cif -k vacancy --supercell 2x2x2 "
                 "-q -1,0,1 -o def --run",
                 "Construye la supercelda con el defecto en cada estado de "
                 "carga, más la supercelda perfecta de referencia.",
                 escribe=["def/FORMACION.dat", "def/q_0/pw.in", "..."],
                 corre_qe=True,
                 ojo="Un número impar de electrones con ocupaciones fijas "
                     "hace que pw.x aborte con «the system is metallic». "
                     "Olla-DFT pone nspin=2 y tot_magnetization en TODOS los "
                     "estados de carga, no solo en el impar, para que sean "
                     "comparables entre sí."),
            Paso("olla-dft eform Si.cif --collect -o def --mu Si:-107.4",
                 "Energías de formación con la corrección de imagen de "
                 "Madelung, y los niveles de transición.",
                 lee=[("las energías de cada estado de carga", "①")],
                 escribe=["def/FORMACION.dat", "def/formacion.pdf"],
                 salida="E_f(q=0) = 3.28 eV   (LDA de la literatura 3.2-3.6)\n"
                        "ε(−1/0) = 0.556 eV sobre el VBM",
                 ojo="El potencial químico μ NO es opcional: E_f depende de "
                     "él linealmente. Un E_f sin decir con qué μ se calculó "
                     "no significa nada."),
            Paso("olla-dft crosscheck . -f Si.cif",
                 "Comprueba que la supercelda del defecto y la perfecta son "
                 "comparables antes de restarlas.",
                 lee=[("los dos cálculos", "①")]),
        ],
        despues=[("si el defecto tiene un nivel dentro del gap, mira sus "
                  "estados", "olla-dft dos def/q_0 -o pdos"),
                 ("y su firma en fotoluminiscencia",
                  "olla-dft optics ...")],
        ver_tambien=["bandas", "fiarme"],
    ),

    Receta(
        clave="modelo",
        palabras=["Wannier", "interpolar", "modelo", "enlace fuerte", "tight binding", "polarización", "ferroeléctrico", "Berry", "rápido"],
        titulo="De DFT a un modelo pequeño que interpola gratis",
        pregunta="quiero la banda en cualquier punto k sin volver a correr "
                 "DFT, y de paso la polarización",
        para_que="Al final tienes H(R): una matriz pequeña con la que la "
                 "energía en CUALQUIER k sale en microsegundos. Con ella se "
                 "hacen superficies de Fermi finas, transporte con mallas "
                 "imposibles y polarización eléctrica.",
        coste="~10 min de máquina para un sistema pequeño",
        verificada=True,
        pasos=[
            Paso("olla-dft wannier Si.cif -g 4x4x4 -p 'Si:sp3' "
                 "--bands 8 --exclude 5-8 -o wann",
                 "Escribe los cuatro pasos: scf, nscf de malla COMPLETA, el "
                 ".nnkp (que normalmente escribe wannier90) y "
                 "pw2wannier90.",
                 escribe=["wann/1_scf.in", "wann/2_nscf.in", "wann/Si.nnkp",
                          "wann/3_pw2wan.in", "wann/4_bands.in"],
                 ojo="El nscf va con nosym y noinv a la fuerza: hacen falta "
                     "TODOS los puntos de la malla, no la cuña irreducible. "
                     "No hace falta tener wannier90 instalado."),
            Paso("olla-dft wannier Si.cif -o wann --run",
                 "Los cuatro pasos en orden, incluido el cálculo de bandas "
                 "de DFT con el que se va a comparar.",
                 lee=[("los inputs", "①")],
                 corre_qe=True),
            Paso("olla-dft wannier Si.cif --collect -o wann",
                 "Proyección, minimización de la dispersión, H(R) e "
                 "interpolación — y la comparación contra las bandas de DFT "
                 "en puntos que NO estaban en la malla.",
                 lee=[("los solapes .amn/.mmn/.eig", "②")],
                 escribe=["wann/WANNIER_hr.dat", "wann/WANNIER_centros.dat",
                          "wann/wannier.pdf"],
                 salida="centros en el enlace Si–Si (a 0.000 Å)\n"
                        "contra DFT fuera de la malla: 275 meV máximo\n"
                        "sin gauge sería 962 meV ← el gauge vale 3.5×",
                 ojo="Que reproduzca la malla de partida es TRIVIAL: es "
                     "interpolación. Lo único que dice si el modelo sirve es "
                     "compararlo en puntos que no estaban, y para eso hace "
                     "falta el paso 4_bands."),
            Paso("olla-dft berry Si.cif --displace 2:0,0,0.16 --nlambda 5 "
                 "-o born --run",
                 "La misma fase de Berry que hay detrás de los centros de "
                 "Wannier, ahora por el camino de lberry: de aquí sale la "
                 "carga efectiva de Born.",
                 escribe=["born/BERRY.dat"],
                 corre_qe=True,
                 salida="Z* = 0.0000 e   (en un cristal homopolar vale cero "
                        "exactamente)",
                 ojo="Un valor de P suelto no significa nada: está definido "
                     "módulo el cuanto. Solo la DIFERENCIA a lo largo de un "
                     "camino es física."),
            Paso("olla-dft crosscheck . -f Si.cif",
                 "Cruza la fase electrónica de lberry contra la que dan los "
                 "centros de Wannier. Dos rutinas que no comparten una línea "
                 "de código.",
                 lee=[("WANNIER_centros.dat y BERRY.dat", "③ y ④")],
                 salida="[OK  ] fase electrónica de Berry\n"
                        "         lberry: 0        centros: −6.4·10⁻⁷"),
        ],
        despues=[("DOS en una malla de 24³ sin tocar pw.x",
                  "olla-dft wannier Si.cif --collect -o wann --dos 24"),
                 ("bandas entrelazadas (conducción, metales): hacen falta "
                  "ventanas",
                  "olla-dft wannier Si.cif --collect -o wann "
                  "--window -10:20 --frozen -10:6.4")],
        ver_tambien=["bandas", "termoelectrico"],
    ),

    Receta(
        clave="fiarme",
        palabras=["validar", "comprobar", "confiar", "verificar", "error", "correcto", "bien", "basura", "revisar"],
        titulo="Cómo saber que lo que calculaste está bien",
        pregunta="ya tengo resultados: ¿cómo sé que no son basura?",
        para_que="Tres comprobaciones distintas y complementarias: que el "
                 "código hace bien las cuentas, que tus cálculos convergieron "
                 "y son comparables entre sí, y que dos rutas independientes "
                 "dan lo mismo.",
        coste="minutos",
        verificada=True,
        pasos=[
            Paso("olla-dft selftest",
                 "Contrasta el propio código contra valores conocidos de la "
                 "literatura, cada uno con su fuente. No necesita ni "
                 "Quantum ESPRESSO ni tus datos.",
                 salida="12 pruebas: 12 bien, 0 fuera de tolerancia\n"
                        "[ ok ] Constante de Madelung  α = 2.8373\n"
                        "(con --full corre también las que necesitan pw.x)",
                 ojo="Si esto falla, el problema no está en tu cálculo: está "
                     "en la instalación."),
            Paso("olla-dft doctor mi_calculo/",
                 "Mira la traza de convergencia y distingue oscilación de "
                 "carga de convergencia lenta, que piden remedios "
                 "OPUESTOS.",
                 lee=[("la salida de pw.x", "tu cálculo")],
                 ojo="Terminar y converger no son lo mismo. Un scf que agota "
                     "electron_maxstep escribe «JOB DONE» igual que uno "
                     "bueno."),
            Paso("olla-dft audit calculo1/ calculo2/ calculo3/",
                 "Comprueba que dos cálculos son comparables antes de restar "
                 "sus energías: funcional, pseudos, cutoffs, ocupaciones.",
                 lee=[("los XML de cada cálculo", "tus cálculos")],
                 ojo="Restar energías incomparables es el error más caro de "
                     "la DFT, y es el único que no da ningún aviso."),
            Paso("olla-dft crosscheck . -f Si.cif",
                 "La misma cantidad por dos caminos independientes. Necesita "
                 "que ya haya DOS módulos con resultados.",
                 lee=[("todo lo que hayan dejado los demás módulos",
                       "cualquier receta anterior")],
                 salida="3 cruces  |  3 coinciden  |  0 NO",
                 ojo="Un cruce que falla NO dice cuál de los dos caminos "
                     "está mal: dice que uno lo está."),
            Paso("olla-dft datasheet . -o ficha",
                 "Junta todo lo calculado en una ficha con su metodología y "
                 "su procedencia, lista para pegar en el artículo.",
                 lee=[("los .dat de todos los módulos", "todos")],
                 escribe=["ficha/ficha_material.md",
                          "ficha/ficha_material.html"]),
        ],
        despues=[("si algo falló, quedó registrado con su traza y sus "
                   "versiones", "olla-dft report --stats")],
        ver_tambien=["mecanicas", "termoelectrico"],
    ),
]

RECETAS_POR_CLAVE = {r.clave: r for r in RECETAS}
CIRCULOS = "①②③④⑤⑥⑦⑧⑨⑩"


# ----------------------------------------------------------------------
# Idioma
# ----------------------------------------------------------------------
def _idioma(language):
    return language or i18n.get_language()


@lru_cache(maxsize=None)
def _tabla(language):
    """El mapa {es: en} de recipes_<idioma>.json; vacío para el español."""
    if language == "es":
        return {}
    return i18n.load_table(f"recipes_{language}").get("strings", {})


@lru_cache(maxsize=None)
def _keywords(language):
    """Palabras de búsqueda en otro idioma: clave de receta -> lista."""
    if language == "es":
        return {r.clave: list(r.palabras) for r in RECETAS}
    return i18n.load_table(f"recipes_{language}").get("keywords", {})


def t(texto, language=None):
    """Un texto fijo de la interfaz, en el idioma pedido (o el activo)."""
    return _tabla(_idioma(language)).get(texto, texto)


@lru_cache(maxsize=None)
def _recetas(language):
    if language == "es":
        return RECETAS
    tabla = _tabla(language)
    return [i18n.translate_data(r, tabla) for r in RECETAS]


def recetas(language=None) -> list:
    """Las recetas en el idioma pedido. En español es la lista original."""
    return _recetas(_idioma(language))


def _traducir(r, language):
    """La receta `r` en el idioma pedido (idempotente: ver translate_data)."""
    if language == "es":
        return r
    return i18n.translate_data(r, _tabla(language))


# ----------------------------------------------------------------------
# Búsqueda, informe y guion
# ----------------------------------------------------------------------
def buscar(texto, language=None):
    """Recetas que encajan con lo que alguien escribe en lenguaje llano.

    Las palabras de `palabras` y de la pregunta pesan el triple que las del
    cuerpo: quien busca "conduce el calor" quiere la receta de conductividad
    térmica, no todas las que mencionan la palabra calor de pasada.

    Se busca a la vez en español y en inglés (con las `keywords` de la tabla
    de traducción), así que "band gap" encuentra la receta de bandas aunque
    la interfaz esté en español. Cada palabra puntúa como mucho una vez: lo
    que vale es el idioma en el que mejor encaja, no la suma de los dos.
    Devuelve las recetas en el idioma pedido.
    """
    t_ = str(texto or "").lower()
    idiomas = [lang for lang in i18n.LANGUAGES
               if lang == "es" or _tabla(lang)]
    cortas = {w.lower() for lang in idiomas
              for ws in _keywords(lang).values() for w in ws}
    # las palabras de 3 letras solo cuentan si son jerga conocida (gap, DOS)
    pedidas = [w for w in re.findall(r"\w+", t_)
               if len(w) > 3 or w in cortas]
    if not pedidas:
        return []
    puntos = []
    for i, r in enumerate(RECETAS):
        fuertes, flojos, frases = [], [], []
        for lang in idiomas:
            rl = _recetas(lang)[i]
            kw = _keywords(lang).get(r.clave, [])
            fuertes.append((f"{r.clave} {rl.pregunta} " + " ".join(kw)).lower())
            flojos.append(f"{rl.titulo} {rl.para_que}".lower())
            frases += kw
        n = 0
        for w in pedidas:
            n += max((3 if w in f else 0) + (1 if w in fl else 0)
                     for f, fl in zip(fuertes, flojos))
        # y la frase entera, si aparece tal cual, manda
        if any(len(frase) > 8 and frase.lower() in t_ for frase in frases):
            n += 6
        if n:
            puntos.append((n, r))
    puntos.sort(key=lambda x: (-x[0], x[1].clave))
    language = _idioma(language)
    return [_traducir(r, language) for _n, r in puntos]


def obtener(clave, language=None):
    language = _idioma(language)
    r = RECETAS_POR_CLAVE.get(str(clave).strip().lower())
    if r is None:
        raise ErrorDeUso(
            t("no hay ninguna receta llamada '{clave}'. Las que hay:",
              language).format(clave=clave) + "\n  "
            + "\n  ".join(f"{x.clave:16s} {x.pregunta}"
                          for x in recetas(language))
            + "\n\n" + t("O búscala con tus palabras:  olla-dft recetas "
                         "--buscar 'quiero saber si conduce'", language))
    return _traducir(r, language)


def _envolver(texto, ancho=72, sangria="     "):
    import textwrap
    fuera = []
    for parrafo in str(texto).split("\n"):
        fuera += textwrap.wrap(parrafo, ancho,
                               initial_indent=sangria,
                               subsequent_indent=sangria) or [sangria.rstrip()]
    return "\n".join(fuera)


def listar(language=None) -> str:
    language = _idioma(language)

    def T(s):
        return t(s, language)

    L = [T("--- Recetas: sesiones completas, de la estructura al resultado ---"),
         "",
         T("Cada una enseña no solo los comandos, sino QUÉ ARCHIVO deja cada "
           "paso y qué"),
         T("paso posterior lo lee. Eso es lo que convierte una lista de "
           "comandos en un"),
         T("flujo de trabajo, y es lo que no está en la ayuda de ningún "
           "comando suelto."),
         ""]
    ancho = max(len(r.clave) for r in RECETAS)
    for r in recetas(language):
        marca = "✓" if r.verificada else " "
        L.append(f"  {marca} {r.clave:<{ancho}}  {r.pregunta}")
        L.append(f"    {'':<{ancho}}  "
                 + T("{n} pasos · {coste}").format(n=len(r.pasos),
                                                   coste=r.coste))
    L += ["",
          T("  ✓ = la receta se ha corrido entera y los extractos de salida "
            "son los reales."),
          "",
          T("Ver una:      olla-dft recetas <clave>"),
          T("Buscar:       olla-dft recetas --buscar 'quiero saber si conduce el "
            "calor'"),
          T("Guionizar:    olla-dft recetas bandas --script mi_sesion.sh")]
    return "\n".join(L)


def report(r: Receta, language=None) -> str:
    language = _idioma(language)
    r = _traducir(r, language)

    def T(s):
        return t(s, language)

    L = [f"--- {r.titulo} ---",
         T("La pregunta: «{pregunta}»").format(pregunta=r.pregunta),
         ""]
    L.append(_envolver(r.para_que, sangria=""))
    L += ["", T("Coste: {coste}").format(coste=r.coste)]
    if not r.verificada:
        L.append(T("NOTA: esta receta no se ha corrido entera de principio a "
                   "fin, así que no lleva\nextractos de salida en todos los "
                   "pasos. Los comandos sí están comprobados "
                   "contra\nel propio programa."))
    L += ["", "=" * 74, ""]

    productores = {}          # archivo -> paso que lo escribió
    for i, p in enumerate(r.pasos):
        n = CIRCULOS[i] if i < len(CIRCULOS) else f"({i + 1})"
        marca = T("  [corre pw.x]") if p.corre_qe else ""
        L.append(f" {n}  $ {p.comando}{marca}")
        L.append(_envolver(p.hace, sangria="     "))
        for texto, de in p.lee:
            L.append(T("     ← lee   {texto}   [de {de}]").format(texto=texto,
                                                                  de=de))
        for f in p.escribe:
            L.append(T("     → deja  {f}").format(f=f))
            productores.setdefault(f, n)
        if p.salida:
            L.append("")
            if language != "es":
                # la salida es la de verdad, y el programa habla español
                L.append("       " + T("(salida real de la corrida, en "
                                       "español)"))
            for linea in p.salida.split("\n"):
                L.append(f"       │ {linea}")
        if p.ojo:
            L.append("")
            L.append(_envolver(T("OJO:") + " " + p.ojo, sangria="     "))
        L.append("")

    if r.despues:
        L += ["=" * 74, "", T("Y después:")]
        for texto, cmd in r.despues:
            L.append(f"  · {texto}")
            L.append(f"      $ {cmd}")
        L.append("")
    if r.ver_tambien:
        L.append(T("Recetas relacionadas:") + " "
                 + ", ".join(f"olla-dft recetas {c}" for c in r.ver_tambien))
    return "\n".join(L)


def script(r: Receta, destino=None, language=None) -> str:
    """Convierte la receta en un guion de shell que se puede correr.

    Sale COMENTADO: cada comando lleva encima lo que hace y el aviso de su
    paso. La idea no es que lo lances a ciegas, sino que tengas el esqueleto
    con las rutas ya encadenadas y lo edites.
    """
    language = _idioma(language)
    r = _traducir(r, language)

    def T(s):
        return t(s, language)

    L = ["#!/bin/bash",
         f"# {r.titulo}",
         f"# {r.pregunta}",
         "#",
         T("# Generado por 'olla-dft recetas {clave} --script'.").format(
             clave=r.clave),
         T("# Cambia el archivo de estructura y las rutas por las tuyas."),
         "set -e", ""]
    for i, p in enumerate(r.pasos):
        L.append(T("# --- paso {n}: {hace}").format(n=i + 1, hace=p.hace))
        for f in p.escribe:
            L.append(T("#     deja: {f}").format(f=f))
        if p.ojo:
            for linea in _envolver(p.ojo, 68, "#     ").split("\n"):
                L.append(linea)
        if p.corre_qe:
            L.append(T("#     ESTE PASO CORRE pw.x Y TARDA"))
        L.append(p.comando)
        L.append("")
    texto = "\n".join(L)
    if destino:
        from pathlib import Path
        d = Path(destino)
        d.parent.mkdir(parents=True, exist_ok=True)
        d.write_text(texto, encoding="utf-8")
        try:
            d.chmod(0o755)
        except OSError:
            pass
    return texto
