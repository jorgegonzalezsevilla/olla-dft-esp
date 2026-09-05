# Olla-DFT — command-line toolkit for Quantum ESPRESSO
# Copyright (C) 2026 Jorge Enrique González Sevilla
# SPDX-License-Identifier: AGPL-3.0-or-later
# This program is free software: you can redistribute it and/or modify it
# under the terms of the GNU Affero General Public License as published by the Free
# Software Foundation, either version 3 of the License, or (at your option)
# any later version. See the LICENSE file for details.

"""Asistente guiado: de lo que quieres SABER a los archivos que hay que correr.

POR QUÉ EXISTE
--------------
Olla-DFT tiene cincuenta subcomandos. Eso no sirve de nada si para llegar al
que te interesa hay que saber ya cómo se llama lo que buscas. La barrera
real de la simulación no es la física: es el vocabulario. Alguien que
quiere saber "de qué color va a ser mi material" tiene que averiguar
primero que eso se llama función dieléctrica, que sale de epsilon.x, que
necesita pseudopotenciales de norma conservada y un nscf previo.

Este módulo va al revés. Se pregunta QUÉ QUIERES SABER, en el idioma en
que uno lo piensa, y de ahí sale:

  - qué cálculo corresponde y por qué,
  - qué hace falta antes (y si ya lo tienes),
  - los comandos exactos, en orden,
  - lo que te va a costar,
  - y el error concreto en el que cae todo el mundo con ese cálculo.

CÓMO ESTÁ HECHO
---------------
Un catálogo de METAS. Cada meta dice qué pregunta responde, con qué
subcomandos se hace, qué necesita antes, y qué se suele hacer mal. No hay
ningún modelo aprendido: es una tabla escrita a mano y auditable. Se le
pueden añadir metas sin tocar la lógica.

El glosario es la otra mitad. Cada término que aparece se explica AL VUELO
la primera vez, en una frase, sin mandar a nadie a leer un manual.

EN INGLÉS
---------
El catálogo y el glosario se escriben una vez, en español. Con
``--language en`` las funciones de informe reciben el idioma y traducen la
estructura al vuelo con ``qekit/data/i18n/wizard_en.json``
(``i18n.translate_data``); los comandos no están en la tabla y quedan como
están. Una prueba comprueba que la tabla cubre todas las cadenas de todas
las metas y del glosario.
"""

import re
from dataclasses import dataclass, field
from functools import lru_cache

import numpy as np

from qekit.core import i18n
from qekit.core.errors import ErrorDeUso

# ----------------------------------------------------------------------
# Glosario: cada término, en una frase
# ----------------------------------------------------------------------
GLOSARIO = {
    "pseudopotencial": (
        "El archivo (.UPF) que sustituye al núcleo y a los electrones de "
        "core de un elemento. Sin uno por cada elemento de tu estructura, "
        "no se puede calcular nada."),
    "cutoff": (
        "Hasta qué detalle se describe la función de onda, en Ry. Más alto "
        "es más exacto y más caro. No hay un valor universal: depende del "
        "pseudopotencial, y hay que converger."),
    "malla k": (
        "Cuántos puntos del espacio recíproco se muestrean. Más densa es "
        "más exacta y más cara. Un metal necesita mucha más que un "
        "aislante."),
    "scf": (
        "El cálculo básico: encuentra la densidad electrónica del estado "
        "fundamental para unas posiciones atómicas fijas. Casi todo lo "
        "demás parte de un scf."),
    "nscf": (
        "Un cálculo que reutiliza la densidad de un scf para obtener "
        "energías en muchos más puntos k, sin volver a resolver la "
        "autoconsistencia. Es lo que hace baratas las bandas y la DOS."),
    "relajación": (
        "Dejar que los átomos se muevan hasta que las fuerzas sean casi "
        "cero. Si tu estructura viene de un artículo o de una base de "
        "datos con OTRO funcional, hay que relajarla antes de nada."),
    "gap": (
        "La distancia en energía entre el último estado ocupado y el "
        "primero vacío. Si es cero, el material conduce; si no, es "
        "semiconductor o aislante."),
    "DOS": (
        "Densidad de estados: cuántos estados electrónicos hay a cada "
        "energía. La PDOS los separa por átomo y por orbital, que es lo "
        "que dice quién aporta qué al enlace."),
    "fonones": (
        "Las vibraciones del cristal. De ellas salen el espectro "
        "infrarrojo y Raman, la capacidad calorífica, la expansión "
        "térmica, y la comprobación de si tu estructura es estable."),
    "supercelda": (
        "Repetir la celda varias veces para meter un defecto, un dopante "
        "o una molécula sin que interactúe consigo misma a través de la "
        "periodicidad."),
    "vacío": (
        "Espacio vacío que se deja en una dirección para simular una "
        "superficie o una monocapa. Si es poco, las dos caras de la losa "
        "se ven entre sí."),
    "funcional": (
        "La receta de intercambio y correlación (PBE, PZ, PBEsol...). "
        "Números calculados con funcionales distintos NO se pueden "
        "comparar entre sí."),
    "convergencia": (
        "Que el resultado ya no cambie al hacer el cálculo más fino. No "
        "es opcional: un número sin converger no es un resultado."),
    "función de Wannier": (
        "Una base localizada en el espacio real, obtenida transformando las "
        "bandas de Bloch. Con ella el hamiltoniano es una matriz pequeña que "
        "decae con la distancia, y la banda en cualquier punto k sale de "
        "diagonalizarla: microsegundos en vez de otro cálculo de DFT."),
    "fase de Berry": (
        "La fase que acumula un estado electrónico al recorrer la zona de "
        "Brillouin. Suena abstracta y es literalmente una posición: de ella "
        "salen los centros de Wannier y la polarización eléctrica."),
    "carga de Born": (
        "Cuánta carga hay que mover para producir la polarización que "
        "aparece al desplazar un átomo. No es la carga iónica nominal: en "
        "un cristal homopolar como el silicio vale cero exactamente."),
    "fc3": (
        "La tercera derivada de la energía respecto de los desplazamientos "
        "atómicos. Es lo que hace que los fonones se dispersen entre sí y, "
        "por tanto, lo que hace finita la conductividad térmica."),
    "recorrido libre medio": (
        "La distancia que recorre un fonón antes de dispersarse. Si la "
        "mayoría del calor la llevan fonones de recorrido largo, hacer el "
        "grano más pequeño que esa distancia baja la conductividad; si no, "
        "no sirve de nada."),
    "ESM": (
        "Medio de apantallamiento efectivo. Sustituye las imágenes "
        "periódicas en la dirección del vacío por una condición de contorno "
        "explícita, que es lo que permite calcular una superficie cargada "
        "sin el fondo compensador que la llena todo."),
    "función trabajo": (
        "La energía que cuesta arrancar un electrón de la superficie y "
        "dejarlo en el vacío: el nivel de vacío menos el nivel de Fermi."),
    "sobrepotencial": (
        "Cuánto voltaje de más hay que aplicar, por encima del "
        "termodinámico, para que una reacción electroquímica corra. Es la "
        "cifra de mérito de un catalizador."),
    "potencial de deformación": (
        "Cuánto se mueve una banda por unidad de deformación. Es lo que "
        "relaciona la tensión mecánica con el cambio del gap."),
    "alineamiento de bandas": (
        "Dónde caen las bandas de un material respecto de las de otro "
        "cuando se ponen en contacto. Decide si una heterounión confina "
        "electrones y huecos juntos (tipo I) o separados (tipo II)."),
}


@dataclass
class Meta:
    """Una cosa que alguien puede querer saber."""
    clave: str
    pregunta: str                 # cómo lo diría alguien que no sabe la jerga
    nombre: str                   # cómo se llama en la literatura
    explica: str                  # qué es y qué NO es
    pasos: list = field(default_factory=list)     # (descripcion, comando)
    necesita: list = field(default_factory=list)  # claves de otras metas
    coste: str = "medio"          # bajo | medio | alto | muy alto
    error_tipico: str = ""
    terminos: list = field(default_factory=list)
    requiere: list = field(default_factory=list)  # condiciones de estructura


#: El catálogo. Añadir una meta es añadir una entrada aquí.
METAS = [
    Meta("estructura",
         "¿qué es esto que tengo? ¿está bien la estructura?",
         "análisis de simetría",
         "Lee tu archivo, dice el grupo espacial, los sitios inequivalentes "
         "y las distancias. Es gratis y es lo primero que hay que mirar: si "
         "la simetría no es la que esperabas, todo lo demás sale mal.",
         pasos=[("mirar la estructura", "olla-dft info {file}"),
                ("ver el camino de alta simetría", "olla-dft kpath {file}")],
         coste="bajo",
         error_tipico="Dar por buena una estructura de una base de datos sin "
                      "mirarla. Una celda mal centrada o con un átomo de más "
                      "no da error: da resultados equivocados.",
         terminos=[]),

    Meta("relajar",
         "¿cuáles son las posiciones y el parámetro de red correctos?",
         "relajación de la estructura",
         "Deja que los átomos y la celda se muevan hasta el mínimo de "
         "energía CON TU funcional. Una estructura experimental o de otra "
         "base de datos no está en el mínimo del tuyo, y esa tensión "
         "residual contamina todo lo demás.",
         pasos=[("convergencia de cutoff y malla k",
                 "olla-dft converge {file} --run"),
                ("relajar posiciones y celda",
                 "olla-dft gen {file} -p vc-relax -o relax")],
         coste="medio",
         error_tipico="Saltarse la convergencia. Una relajación con un "
                      "cutoff bajo da un parámetro de red que parece "
                      "razonable y está mal en un 2 %.",
         terminos=["cutoff", "malla k", "convergencia", "relajación"]),

    Meta("conduce",
         "¿es metal, semiconductor o aislante? ¿cuánto vale su gap?",
         "estructura de bandas y gap",
         "Calcula la energía de los electrones a lo largo del camino de "
         "alta simetría. De ahí sale el gap, si es directo o indirecto, y "
         "dónde están el máximo de la banda de valencia y el mínimo de la "
         "de conducción.\n"
         "OJO: la DFT con PBE o LDA SUBESTIMA el gap, típicamente entre un "
         "30 y un 50 %. El valor que salga no es el experimental, y eso no "
         "es un fallo del cálculo: es una limitación conocida del método.",
         pasos=[("scf, nscf y bandas", "olla-dft gen {file} -p all -o bandas"),
                ("(correr los inputs con pw.x)", None),
                ("analizar y graficar", "olla-dft bands bandas"),
                ("solo el número del gap", "olla-dft gap bandas")],
         necesita=["relajar"],
         coste="medio",
         error_tipico="Comparar el gap calculado con el experimental y "
                      "concluir que el cálculo está mal. Está subestimado "
                      "por construcción.",
         terminos=["scf", "nscf", "gap"]),

    Meta("quien_aporta",
         "¿qué átomo y qué orbital aportan a cada parte del espectro?",
         "DOS y PDOS",
         "La densidad de estados dice cuántos estados hay a cada energía; "
         "la PROYECTADA los reparte por átomo y por orbital. Es lo que "
         "contesta 'este pico de aquí, ¿es del oxígeno o del metal?'.",
         pasos=[("scf y nscf denso", "olla-dft gen {file} -p dos -o dos"),
                ("(correr pw.x, dos.x y projwfc.x)", None),
                ("analizar", "olla-dft dos dos --mode orbital")],
         necesita=["relajar"],
         coste="medio",
         error_tipico="Usar la misma malla k que en el scf. La DOS necesita "
                      "una mucho más densa o sale con picos falsos.",
         terminos=["DOS", "nscf", "malla k"]),

    Meta("color",
         "¿de qué color va a ser? ¿absorbe luz visible?",
         "propiedades ópticas",
         "Calcula la función dieléctrica y de ahí la absorción, el índice "
         "de refracción y la reflectividad. El gráfico de Tauc da el gap "
         "óptico, que es el que se mide en un UV-Vis.",
         pasos=[("preparar", "olla-dft optics {file} -o optica"),
                ("(correr pw.x y epsilon.x)", None),
                ("analizar", "olla-dft optics {file} --collect -o optica "
                             "--tauc indirect")],
         necesita=["relajar"],
         coste="alto",
         error_tipico="epsilon.x SOLO funciona con pseudopotenciales de "
                      "norma conservada. Con ultrasuaves o PAW da números "
                      "sin quejarse, y están mal.",
         terminos=["funcional", "gap"]),

    Meta("difractograma",
         "¿cómo se vería en un difractómetro de rayos X?",
         "difractograma de polvos simulado",
         "Calcula las posiciones e intensidades de los picos de difracción "
         "a partir de la estructura. Sirve para comparar con tu medida y "
         "confirmar que la fase es la que crees.",
         pasos=[("simular", "olla-dft xrd {file} --wavelength 1.5406"),
                ("comparar con tu medida", "olla-dft xrd {file} --exp "
                                           "medida.xy")],
         coste="bajo",
         error_tipico="Comparar intensidades sin tener en cuenta la "
                      "orientación preferente: una muestra prensada no da "
                      "las intensidades de un polvo ideal.",
         terminos=[]),

    Meta("estable",
         "¿esta estructura es estable? ¿no se va a deshacer?",
         "fonones y criterio de estabilidad",
         "Si alguna frecuencia sale imaginaria, la estructura NO es un "
         "mínimo: hay una deformación que baja la energía. Es la prueba "
         "más dura que se le puede hacer a una estructura propuesta.",
         pasos=[("preparar la DFPT", "olla-dft phonons {file} -o fonones "
                                     "--qgrid 2x2x2"),
                ("(correr ph.x, q2r.x y matdyn.x)", None),
                ("analizar", "olla-dft phonons {file} --collect -o fonones")],
         necesita=["relajar"],
         coste="muy alto",
         error_tipico="Calcular fonones sobre una estructura mal relajada. "
                      "Las frecuencias imaginarias que salen son de la "
                      "relajación que faltaba, no del material.",
         terminos=["fonones", "relajación"]),

    Meta("mecanicas",
         "¿es duro? ¿es frágil? ¿cuánto se deforma?",
         "constantes elásticas",
         "De las Cij salen el módulo de bulk, el de Young, la razón de "
         "Poisson y el criterio de Pugh (frágil o dúctil). También la "
         "velocidad del sonido y una estimación de la conductividad "
         "térmica.",
         pasos=[("preparar el barrido", "olla-dft elastic {file} -o elastico"),
                ("(correr los inputs)", None),
                ("analizar", "olla-dft elastic {file} --collect -o elastico"),
                ("cantidades derivadas", "olla-dft derived {file} --cij "
                                         "elastico/ELASTIC_C.dat")],
         necesita=["relajar"],
         coste="alto",
         error_tipico="Deformaciones demasiado grandes: se sale del régimen "
                      "lineal y las Cij salen mal. Y demasiado pequeñas: se "
                      "ahoga en el ruido numérico.",
         terminos=["relajación", "convergencia"]),

    Meta("vibra",
         "¿qué picos voy a ver en el Raman o en el infrarrojo?",
         "espectro Raman e IR",
         "Los modos vibracionales en el punto Gamma, con sus intensidades. "
         "En un cristal con centro de inversión, un modo activo en Raman "
         "NO lo es en IR y al revés: esa regla de exclusión mutua es la "
         "primera comprobación de que el cálculo está bien.",
         pasos=[("preparar", "olla-dft phonons {file} -o raman --gamma --raman"),
                ("(correr ph.x y dynmat.x)", None),
                ("analizar", "olla-dft phonons {file} --collect -o raman "
                             "--gamma --raman --laser 532")],
         necesita=["relajar"],
         coste="alto",
         error_tipico="Comparar frecuencias calculadas con medidas sin "
                      "recordar que la armónica sobreestima típicamente un "
                      "2-5 %.",
         terminos=["fonones"]),

    Meta("superficie",
         "¿cuánto cuesta cortar este cristal? ¿cuál es su función trabajo?",
         "superficies y función trabajo",
         "Corta una superficie con los índices que pidas, le pone vacío, y "
         "calcula la función trabajo desde el potencial electrostático.",
         pasos=[("cortar", "olla-dft surface {file} -m '1 1 1' -l 6 "
                           "--vacuum 20 -o losa.cif"),
                ("generar el scf", "olla-dft gen losa.cif -p scf -o super"),
                ("(correr pw.x)", None),
                ("función trabajo", "olla-dft wf super")],
         necesita=["relajar"],
         coste="alto",
         error_tipico="Poco vacío. Con menos de 15 Å las dos caras de la "
                      "losa se ven, y la función trabajo sale mal.",
         terminos=["vacío", "supercelda"]),

    Meta("defecto",
         "¿qué le pasa a mi material si le meto un dopante o le quito un "
         "átomo?",
         "defectos puntuales y desdoblamiento de bandas",
         "Construye la supercelda con el defecto, y desdobla las bandas "
         "para ver qué le pasó a la estructura electrónica del material "
         "original — que es lo que no se ve en las bandas plegadas de una "
         "supercelda.",
         pasos=[("crear el defecto", "olla-dft defect {file} -k substitution "
                                     "--site 0 --new-element X "
                                     "--supercell 2x2x2 -o defecto"),
                ("bandas de la supercelda", "olla-dft gen defecto/*.cif "
                                            "-p all -o bandas_def"),
                ("(correr pw.x)", None),
                ("desdoblar", "olla-dft unfold bandas_def {file}")],
         necesita=["relajar"],
         coste="muy alto",
         error_tipico="Supercelda pequeña: el defecto ve sus propias "
                      "imágenes y la banda que sale es la de una red de "
                      "defectos, no la de un defecto aislado.",
         terminos=["supercelda", "gap"]),

    Meta("laminar",
         "¿es un material laminar? ¿cuánto cuesta exfoliarlo?",
         "capas y energía de exfoliación",
         "Detecta las capas por conectividad, mide el espaciado basal y "
         "calcula lo que cuesta separar una monocapa del bulk.",
         pasos=[("detectar capas", "olla-dft layers {file} --slab mono.cif"),
                ("preparar la exfoliación", "olla-dft exfoliate {file} "
                                            "-o exfo --vdw grimme-d3"),
                ("(correr pw.x)", None),
                ("analizar", "olla-dft exfoliate {file} --collect -o exfo")],
         coste="alto",
         error_tipico="Calcular exfoliación sin corrección de dispersión. "
                       "PBE a secas no describe la interacción entre capas "
                       "y la energía sale casi cero.",
         terminos=["vacío", "funcional"]),

    Meta("carga",
         "¿cómo se reparte la carga? ¿quién le da electrones a quién?",
         "cargas atómicas y densidad de carga",
         "Cargas de Bader y de Löwdin, y la diferencia de densidad entre "
         "el compuesto y sus partes. Es lo que cuantifica la "
         "transferencia de carga en un enlace.",
         pasos=[("scf con proyecciones", "olla-dft gen {file} -p dos -o carga"),
                ("(correr pw.x, projwfc.x y pp.x)", None),
                ("analizar", "olla-dft charges {file} --lowdin "
                             "carga/projwfc.out --bader carga/densidad.cube")],
         necesita=["relajar"],
         coste="medio",
         error_tipico="Comparar cargas de Bader con cargas formales. No son "
                      "lo mismo ni tienen por qué parecerse.",
         terminos=["DOS"]),

    Meta("termoelectrico",
         "¿sirve como termoeléctrico? ¿cuánto vale su Seebeck?",
         "transporte termoeléctrico",
         "Coeficiente Seebeck, conductividad y factor de potencia en "
         "función del dopaje y la temperatura. Por omisión usa tiempo de "
         "relajación constante (CRTA), que da sigma/tau; con el módulo de "
         "electrón-fonón sale el tau de verdad.",
         pasos=[("nscf muy denso", "olla-dft transport {file} -o transporte "
                                   "--grid 24x24x24"),
                ("(correr pw.x)", None),
                ("analizar", "olla-dft transport {file} --collect "
                             "-o transporte"),
                ("(opcional) tau real", "olla-dft elph {file} -o elph")],
         necesita=["relajar"],
         coste="alto",
         error_tipico="Malla k insuficiente. El Seebeck es una derivada de "
                      "la estructura de bandas y necesita muchísimos más "
                      "puntos k que la energía total.",
         terminos=["nscf", "malla k"]),

    Meta("oxido",
         "mi material es un óxido de metal de transición y sale metálico "
         "cuando no debería",
         "DFT+U con U calculado",
         "El problema clásico: la DFT normal deslocaliza de más los "
         "electrones d y f, y un aislante sale metálico. DFT+U lo corrige, "
         "pero necesita un parámetro U — que casi todo el mundo copia de "
         "un artículo. Aquí se calcula para TU sistema.",
         pasos=[("calcular U autoconsistente",
                 "olla-dft hubbard {file} -o hubbard --cycle --nspin 2 "
                 "--mag 'M=0.5'"),
                ("usarlo", "olla-dft gen {file} --hubbard M=4.1 -o conU")],
         necesita=["relajar"],
         coste="muy alto",
         error_tipico="Usar un U de la literatura calculado con OTRO "
                      "esquema de proyección. No es el mismo número.",
         terminos=["funcional", "convergencia"]),

    Meta("superficie_quimica",
         "¿se adsorbe esta molécula? ¿cuál es la barrera de la reacción?",
         "adsorción, NEB y termoquímica",
         "La energía de adsorción electrónica, la barrera del camino de "
         "reacción, y las correcciones térmicas que hacen falta para "
         "compararlas con un experimento.",
         pasos=[("cortar la superficie", "olla-dft surface {file} -m '1 1 1' "
                                          "-l 5 --vacuum 20 -o losa.cif"),
                ("camino de reacción", "olla-dft neb inicial.cif final.cif "
                                       "-o camino --images 7"),
                ("(correr neb.x)", None),
                ("analizar", "olla-dft neb inicial.cif --collect -o camino"),
                ("correcciones térmicas", "olla-dft thermochem frecuencias.dat "
                                          "--phase gas --structure mol.xyz")],
         necesita=["relajar"],
         coste="muy alto",
         error_tipico="Dar la barrera electrónica como energía de "
                      "activación. Sin punto cero ni entropía puede estar a "
                      "más de 0.5 eV de la que se mide.",
         terminos=["vacío", "supercelda", "relajación"]),

    Meta("caracterizar",
         "¿qué le pasa a este átomo en concreto? ¿en qué estado de "
         "oxidación está?",
         "XPS y XANES",
         "Dos espectroscopías que miran UN átomo: los corrimientos de "
         "nivel de core (XPS) y la absorción de rayos X cerca del borde "
         "(XANES). Las dos necesitan un pseudopotencial con hueco de core, "
         "que Olla-DFT genera.",
         pasos=[("generar los pseudos", "olla-dft corehole EL --edge K "
                                        "-o pseudos"),
                ("XPS", "olla-dft xps {file} --core-hole EL=pseudos/"
                        "EL.huecols.UPF -o xps"),
                ("XANES", "olla-dft xanes {file} --element EL --core-hole "
                          "pseudos/EL.huecols.UPF --average -o xanes")],
         necesita=["relajar"],
         coste="alto",
         error_tipico="Un pseudopotencial con hueco de core NO viene en "
                      "las tablas estándar. Sin él, XPS devuelve una tabla "
                      "de ceros sin dar error.",
         terminos=["pseudopotencial", "supercelda"]),

    Meta("liquido",
         "¿cómo se mueven los átomos? ¿difunde algo?",
         "dinámica molecular",
         "De una trayectoria salen la función de distribución radial "
         "(estructura), el desplazamiento cuadrático medio (difusión) y la "
         "densidad de estados vibracional, esta última incluyendo "
         "anarmonicidad.",
         pasos=[("(correr pw.x con calculation='md' y nosym=.true.)", None),
                ("analizar", "olla-dft md md.out --skip 200")],
         necesita=["relajar"],
         coste="muy alto",
         error_tipico="Analizar el equilibrado. Los primeros picosegundos "
                      "no son la trayectoria de equilibrio y sesgan todo.",
         terminos=["supercelda"]),

    Meta("interfase",
         "quiero poner un material sobre otro",
         "heteroestructuras",
         "Busca la supercelda común que hace que las dos redes encajen con "
         "la menor deformación, y arma el apilamiento.",
         pasos=[("ver las candidatas", "olla-dft interface abajo.cif arriba.cif "
                                       "--list"),
                ("construir", "olla-dft interface abajo.cif arriba.cif "
                              "--index 0 -o interfase")],
         coste="alto",
         error_tipico="Aceptar la celda más pequeña sin mirar la "
                      "deformación. Por encima del 3 % ya no estás "
                      "modelando tu material.",
         terminos=["supercelda", "vacío"]),

    Meta("comparar",
         "quiero comparar dos materiales o dos cálculos entre sí",
         "auditoría de consistencia",
         "Antes de restar dos energías hay que estar seguro de que son "
         "comparables: mismo funcional, mismos pseudopotenciales, mismos "
         "cutoffs, mismas ocupaciones. Restar energías incomparables es el "
         "error más caro de la DFT porque no da ningún aviso.",
         pasos=[("auditar", "olla-dft audit calculo1/ calculo2/"),
                ("indexar", "olla-dft db calculo*/ --db mis_calculos.db"),
                ("estabilidad de fases", "olla-dft hull calculo*/ -o casco")],
         coste="bajo",
         error_tipico="Es el error mismo: comparar sin auditar.",
         terminos=["funcional", "cutoff"]),

    Meta("estirar",
         "¿qué le pasa al material si lo estiro o lo comprimo?",
         "barrido de deformación e ingeniería de bandas",
         "Aplica deformaciones controladas y sigue el gap, la energía y el "
         "momento. De ahí salen el potencial de deformación, el módulo "
         "biaxial y, si el gap cruza cero, la deformación a la que el "
         "material se vuelve metálico.\n"
         "No confundir con relajar: aquí la deformación se IMPONE y se "
         "mantiene; es lo que pasa en una lámina crecida sobre otro "
         "sustrato.",
         pasos=[("barrido biaxial de −4 % a +4 %",
                 "olla-dft strain {file} -m biaxial -r -4:4:9 -o deformacion "
                 "--run"),
                ("analizar", "olla-dft strain {file} --collect -o deformacion")],
         necesita=["relajar"],
         coste="medio",
         error_tipico="Comparar el gap de dos deformaciones con distinto "
                      "número de bandas. Si nbnd cambia, el LUMO puede no "
                      "estar ahí y el gap sale inventado.",
         terminos=["potencial de deformación", "gap", "relajación"]),

    Meta("juntar",
         "voy a poner dos materiales en contacto: ¿dónde quedan sus bandas?",
         "alineamiento de bandas",
         "Da los desplazamientos ΔE_v y ΔE_c y el tipo de heterounión (I, II "
         "o III), que es lo que decide si el par sirve para un LED, para una "
         "celda solar o para nada.\n"
         "OJO: los máximos de valencia de dos cálculos distintos NO son "
         "comparables tal cual. Cada celda tiene su propio cero de energía, "
         "arbitrario. Hay que referirlos al vacío o al potencial "
         "macroscópico de la interfase, y eso es justamente lo que hace este "
         "comando.",
         pasos=[("cada material con su vacío",
                 "olla-dft gen A.cif -p scf -o ladoA   # y lo mismo con B"),
                ("(correr los dos)", None),
                ("alinear", "olla-dft align ladoA ladoB -o alineamiento")],
         necesita=["relajar"],
         coste="medio",
         error_tipico="Restar los VBM crudos de dos cálculos. Cambiando "
                      "solo el vacío de 16 a 22 Å el VBM del hBN se mueve "
                      "0.60 eV sin que la física cambie nada.",
         terminos=["alineamiento de bandas", "vacío", "gap"]),

    Meta("interpolar",
         "quiero la banda en cualquier punto k sin volver a correr DFT",
         "funciones de Wannier",
         "Baja la estructura de bandas a un modelo pequeño en el espacio "
         "real. Una vez construido, la energía en CUALQUIER punto k sale de "
         "diagonalizar una matriz diminuta: mallas de 50³ puntos que en DFT "
         "serían imposibles se vuelven gratis. Es lo que necesitan el "
         "transporte, las superficies de Fermi finas y las masas efectivas "
         "por diferencias.\n"
         "No hace falta instalar wannier90: Olla-DFT escribe el .nnkp y hace la "
         "proyección y la minimización él mismo.",
         pasos=[("preparar los cuatro pasos",
                 "olla-dft wannier {file} -g 4x4x4 -p 'Si:sp3' -o wann"),
                ("correrlos", "olla-dft wannier {file} -o wann --run"),
                ("analizar y comparar con las bandas de DFT",
                 "olla-dft wannier {file} --collect -o wann")],
         necesita=["relajar", "conduce"],
         coste="medio",
         error_tipico="Transformar directamente las energías propias sin "
                      "gauge. La fase que devuelve un diagonalizador es "
                      "arbitraria y cambia de un punto k al siguiente, así "
                      "que lo que transformas no es suave y la "
                      "interpolación sale mal por un factor 4.",
         terminos=["función de Wannier", "malla k", "gap"]),

    Meta("ferroelectrico",
         "¿es ferroeléctrico? ¿cuánta polarización tiene?",
         "polarización por fase de Berry",
         "La polarización de un sólido periódico NO es la integral del "
         "dipolo en la celda: esa integral depende de dónde pongas los "
         "bordes. Lo observable es la fase de Berry, y está definida módulo "
         "un cuanto.\n"
         "Consecuencia práctica: un valor de P suelto no significa nada. "
         "Hay que dar la DIFERENCIA entre la estructura polar y una de "
         "referencia centrosimétrica, recorriendo un camino con pasos "
         "pequeños. De aquí salen también las cargas efectivas de Born.",
         pasos=[("camino adiabático de la referencia a la polar",
                 "olla-dft berry {file} -r referencia.cif --nlambda 7 "
                 "-o polarizacion --run"),
                ("cargas de Born desplazando un átomo",
                 "olla-dft berry {file} --displace 1:0,0,0.1 -o born --run")],
         necesita=["relajar"],
         coste="alto",
         error_tipico="Dar un valor de P de un solo cálculo. Está definido "
                      "módulo el cuanto, así que por sí solo no es "
                      "publicable; y si el camino tiene pasos grandes, el "
                      "seguimiento de la rama se salta un cuanto entero y "
                      "el resultado sigue pareciendo razonable.",
         terminos=["fase de Berry", "carga de Born"]),

    Meta("conduce_calor",
         "¿conduce bien el calor? ¿me sirve para un termoeléctrico?",
         "conductividad térmica de red",
         "Los fonones armónicos no conducen calor de forma finita: κ sale "
         "del término cúbico de la energía, el que permite que un fonón se "
         "parta en dos. Hace falta la fc3, que es cara porque necesita una "
         "derivada por cada TRÍO de átomos.\n"
         "Además de κ(T) sale la curva acumulada frente al recorrido libre "
         "medio, que es la que dice si nanoestructurar sirve para bajar κ o "
         "no.",
         pasos=[("explorar con un potencial aprendido, para elegir la "
                 "supercelda", "olla-dft kappa {file} --dim 2x2x2 --model mace"),
                ("preparar el cálculo de verdad",
                 "olla-dft kappa {file} --dim 2x2x2 -o kL"),
                ("(correr todas las configuraciones)", None),
                ("resolver la ecuación de Boltzmann",
                 "olla-dft kappa {file} --collect -o kL --mesh 19 --isotopes")],
         necesita=["relajar"],
         coste="muy alto",
         error_tipico="Converger solo una cosa. κ tiene que converger en el "
                      "tamaño de la supercelda Y en la malla de q a la vez, "
                      "y además la RTA subestima el resultado entre un 10 y "
                      "un 15 % en silicio.",
         terminos=["fc3", "recorrido libre medio", "fonones", "supercelda"]),

    Meta("electrodo",
         "quiero una superficie cargada, o a un potencial dado",
         "superficies con medio de apantallamiento efectivo (ESM)",
         "En una celda periódica normal, una superficie CARGADA se calcula "
         "con un fondo uniforme compensador que llena todo el vacío. Ese "
         "fondo no existe en ningún experimento y la energía que sale no "
         "converge a nada.\n"
         "ESM sustituye las imágenes en z por una condición de contorno. "
         "Con bc1 (losa neutra) el nivel de vacío vale cero por "
         "construcción, así que la función trabajo es directamente −E_F, y "
         "deja de hacer falta un vacío enorme.",
         pasos=[("función trabajo de la losa neutra",
                 "olla-dft esm {file} --bc bc1 -o superficie --run"),
                ("barrido de carga con electrodo al otro lado",
                 "olla-dft esm {file} --bc bc3 --charge -0.04,0,0.04 "
                 "-o cargada --run")],
         necesita=["superficie"],
         coste="medio",
         error_tipico="Dejar la losa centrada donde la deja el constructor "
                      "de estructuras. ESM mide z desde el CENTRO de la "
                      "celda: una losa en c/2 cae justo sobre su frontera y "
                      "el resultado son cientos de Ry de error, sin ningún "
                      "mensaje.",
         terminos=["ESM", "función trabajo", "vacío"]),

    Meta("catalizador",
         "¿sirve como catalizador para producir hidrógeno u oxígeno?",
         "electrodo de hidrógeno computacional (HER / OER)",
         "Convierte energías de adsorción en un diagrama de energía libre "
         "frente al potencial aplicado, y de ahí salen el potencial "
         "limitante y el sobrepotencial, que es la cifra que se compara "
         "entre catalizadores.\n"
         "El truco del método es evitar calcular el protón solvatado: se "
         "usa H⁺ + e⁻ ⇌ ½H₂ a potencial cero, que es exacto por definición "
         "del electrodo de hidrógeno.",
         pasos=[("sitios de adsorción y energías",
                 "olla-dft adsorb losa.cif --molecule H -o sitios --run"),
                ("diagrama de energía libre",
                 "olla-dft echem --her -0.33 -U 0 --ph 0 -o her")],
         necesita=["superficie_quimica"],
         coste="alto",
         error_tipico="Comparar energías de adsorción sin las correcciones "
                      "térmicas (ZPE − TΔS). Son décimas de eV y del mismo "
                      "tamaño que las diferencias que se discuten.",
         terminos=["sobrepotencial"]),

    Meta("vidrio",
         "quiero un vidrio o un material amorfo, no un cristal",
         "generación de estructuras amorfas por fundido y temple",
         "Empaqueta los átomos a la densidad que le pidas, los funde a alta "
         "temperatura y los enfría. Con un potencial aprendido cuesta "
         "minutos en vez de días.\n"
         "OJO con la velocidad de temple: un vidrio real se enfría a 1-100 "
         "K/s y una simulación va a 10¹²-10¹⁵ K/s. Son diez órdenes de "
         "magnitud, y la estructura sale más desordenada y menos densa que "
         "la real.",
         pasos=[("fundir y templar",
                 "olla-dft amorphous SiO2 -n 8 -d 2.2 --quench-steps 5000 "
                 "-o amorfo"),
                ("relajar con DFT antes de calcularle nada",
                 "olla-dft gen amorfo/amorfo.cif -p relax -o relajado")],
         coste="medio",
         error_tipico="Quedarse con UNA realización. Un amorfo no es una "
                      "estructura, es un conjunto: hay que generar varias "
                      "con semillas distintas y promediar.",
         terminos=["supercelda", "relajación"]),

    Meta("fiarme",
         "¿cómo sé que todo esto que he calculado está bien?",
         "validación cruzada y contraste con la literatura",
         "Tres cosas distintas, y las tres hacen falta:\n"
         "  · `selftest` contrasta el propio código contra valores "
         "conocidos de la literatura, cada uno con su fuente.\n"
         "  · `crosscheck` calcula la misma cantidad por dos rutas "
         "independientes y compara. Que coincidan es la evidencia más "
         "fuerte que se puede tener sin experimento.\n"
         "  · `doctor` y `audit` miran si cada cálculo convergió y si dos "
         "cálculos son comparables entre sí.",
         pasos=[("validar el código contra la literatura", "olla-dft selftest"),
                ("cruzar tus propios resultados por rutas distintas",
                 "olla-dft crosscheck . -f {file}"),
                ("comprobar que son comparables", "olla-dft audit calculo*/")],
         coste="bajo",
         error_tipico="Fiarse de que un cálculo terminó. Terminar y "
                      "converger no son lo mismo, y un scf que agota "
                      "electron_maxstep escribe JOB DONE igual.",
         terminos=["convergencia", "funcional"]),
]

METAS_POR_CLAVE = {m.clave: m for m in METAS}


# ----------------------------------------------------------------------
# Idioma
# ----------------------------------------------------------------------
def _idioma(language):
    return language or i18n.get_language()


@lru_cache(maxsize=None)
def _tabla(language):
    """El mapa {es: en} de wizard_<idioma>.json; vacío para el español."""
    if language == "es":
        return {}
    return i18n.load_table(f"wizard_{language}").get("strings", {})


@lru_cache(maxsize=None)
def _keywords(language):
    """Palabras de búsqueda en otro idioma: clave de meta -> lista."""
    if language == "es":
        return {}
    return i18n.load_table(f"wizard_{language}").get("keywords", {})


def t(texto, language=None):
    """Un texto fijo de la interfaz, en el idioma pedido (o el activo)."""
    return _tabla(_idioma(language)).get(texto, texto)


@lru_cache(maxsize=None)
def _metas(language):
    if language == "es":
        return METAS
    tabla = _tabla(language)
    return [i18n.translate_data(m, tabla) for m in METAS]


def metas(language=None) -> list:
    """El catálogo en el idioma pedido. En español es la lista original."""
    return _metas(_idioma(language))


def metas_por_clave(language=None) -> dict:
    return {m.clave: m for m in metas(language)}


@lru_cache(maxsize=None)
def _glosario(language):
    if language == "es":
        return GLOSARIO
    return i18n.translate_data(GLOSARIO, _tabla(language))


def glosario(language=None) -> dict:
    """El glosario en el idioma pedido: también los términos se traducen."""
    return _glosario(_idioma(language))


def _traducir(m, language):
    """La meta `m` en el idioma pedido (idempotente: ver translate_data)."""
    if language == "es":
        return m
    return i18n.translate_data(m, _tabla(language))


# ----------------------------------------------------------------------
# Diagnóstico de la estructura
# ----------------------------------------------------------------------
@dataclass
class Diagnostico:
    formula: str = ""
    natoms: int = 0
    grupo: str = ""
    elementos: list = field(default_factory=list)
    es_laminar: bool = False
    tiene_vacio: bool = False
    metales_transicion: list = field(default_factory=list)
    pseudos_faltan: list = field(default_factory=list)
    notas: list = field(default_factory=list)
    # las mismas notas sin formatear: (plantilla, campos), para traducirlas
    notas_src: list = field(default_factory=list)


_NOTA_VACIO = ("Hay {hueco:.1f} Å de vacío a lo largo del eje {eje}: esto "
               "parece una losa o una monocapa, no un sólido en volumen.")
_NOTA_TM = ("Hay metales de transición o tierras raras ({lista}). Si el "
            "material es un óxido aislante y te\nsale metálico, el problema "
            "es la autointeracción: mira la meta 'oxido'.")
_GRUPO_DESCONOCIDO = "no determinado"


def diagnosticar(atoms, pseudo_dir: str = None) -> Diagnostico:
    """Lo que se puede decir de una estructura sin calcular nada."""
    from qekit.core import pseudo as ps
    from qekit.core import structure as st
    from qekit.modules.hubbard import ORBITAL_HUBBARD

    d = Diagnostico(formula=atoms.get_chemical_formula(),
                    natoms=len(atoms))
    d.elementos = list(dict.fromkeys(atoms.get_chemical_symbols()))
    try:
        ds = st.symmetry_dataset(atoms)
        d.grupo = f"{ds.international} (N.º {ds.number})"
    except Exception:                                  # noqa: BLE001
        d.grupo = _GRUPO_DESCONOCIDO

    def nota(plantilla, **campos):
        d.notas.append(plantilla.format(**campos))
        d.notas_src.append((plantilla, campos))

    celda = np.array(atoms.get_cell())
    pos = atoms.get_positions()
    for eje in range(3):
        largo = float(np.linalg.norm(celda[eje]))
        proy = pos @ (celda[eje] / max(largo, 1e-9))
        hueco = largo - (proy.max() - proy.min())
        if hueco > 8.0:
            d.tiene_vacio = True
            nota(_NOTA_VACIO, hueco=hueco, eje="abc"[eje])
            break

    d.metales_transicion = [e for e in d.elementos if e in ORBITAL_HUBBARD]
    if d.metales_transicion:
        nota(_NOTA_TM, lista=", ".join(d.metales_transicion))

    if pseudo_dir:
        res = ps.resolve(d.elementos, pseudo_dir)
        d.pseudos_faltan = [e for e, v in res.items() if not v["found"]]
    return d


# ----------------------------------------------------------------------
# Plan
# ----------------------------------------------------------------------
def plan(meta_clave: str, archivo: str = "estructura.cif",
         vistos: set = None, profundidad: int = 0, language=None) -> list:
    """Los pasos de una meta, con los de sus prerrequisitos delante."""
    language = _idioma(language)
    if meta_clave not in METAS_POR_CLAVE:
        raise ErrorDeUso(
            t("meta '{meta}' desconocida. Disponibles:", language).format(
                meta=meta_clave)
            + " " + ", ".join(m.clave for m in METAS))
    vistos = vistos if vistos is not None else set()
    if meta_clave in vistos:
        return []
    vistos.add(meta_clave)
    m = metas_por_clave(language)[meta_clave]
    fuera = []
    for prev in m.necesita:
        fuera += plan(prev, archivo, vistos, profundidad + 1, language)
    for desc, cmd in m.pasos:
        fuera.append((m.clave, desc,
                      cmd.replace("{file}", archivo) if cmd else None))
    return fuera


def buscar(texto: str, n: int = 4, language=None) -> list:
    """Metas que encajan con lo que alguien escribió en sus palabras.

    Se busca a la vez en español y en inglés (más las `keywords` de la
    tabla de traducción, que pesan como la pregunta); cada palabra puntúa
    una sola vez, por el idioma en el que mejor encaja. Devuelve las metas
    en el idioma pedido.
    """
    palabras = [w for w in re.split(r"[^\wáéíóúñü]+", texto.lower()) if
                len(w) > 2]
    if not palabras:
        return []
    idiomas = [lang for lang in i18n.LANGUAGES
               if lang == "es" or _tabla(lang)]
    puntuados = []
    for i, m in enumerate(METAS):
        variantes = []
        for lang in idiomas:
            ml = _metas(lang)[i]
            kw = " ".join(_keywords(lang).get(m.clave, [])).lower()
            fuerte = f"{ml.pregunta.lower()} {kw}"
            variantes.append((fuerte, ml.nombre.lower(),
                              " ".join([fuerte, ml.nombre, ml.explica,
                                        m.clave]).lower()))
        p = sum(max(3 if w in fuerte else
                    2 if w in nombre else
                    1 if w in campo else 0
                    for fuerte, nombre, campo in variantes)
                for w in palabras)
        if p:
            puntuados.append((p, m))
    puntuados.sort(key=lambda t_: -t_[0])
    language = _idioma(language)
    return [_traducir(m, language) for _, m in puntuados[:n]]


def report_meta(m: Meta, archivo: str = "estructura.cif",
                glosario: bool = True, language=None) -> str:
    language = _idioma(language)
    m = _traducir(m, language)
    por_clave = metas_por_clave(language)
    glos = _glosario(language)

    def T(s):
        return t(s, language)

    lines = [f"--- {m.nombre} ---",
             T('La pregunta: "{pregunta}"').format(pregunta=m.pregunta),
             "", m.explica, ""]
    if m.necesita:
        lines.append(T("Antes hace falta:") + " " + ", ".join(
            por_clave[k].nombre for k in m.necesita))
        lines.append("")
    lines.append(T("Coste: {coste}").format(coste=m.coste))
    lines.append("")
    lines.append(T("Pasos:"))
    for i, (clave, desc, cmd) in enumerate(
            plan(m.clave, archivo, language=language), start=1):
        marca = "" if clave == m.clave else T("  [de {clave}]").format(
            clave=clave)
        if cmd:
            lines.append(f"  {i}. {desc}{marca}")
            lines.append(f"     $ {cmd}")
        else:
            lines.append(f"  {i}. {desc}{marca}")
    if m.error_tipico:
        lines += ["", T("El error en el que cae todo el mundo:"),
                  "  " + m.error_tipico]
    if glosario and m.terminos:
        lines += ["", T("Términos que salen aquí:")]
        for term in m.terminos:
            if term in glos:
                lines.append(f"  {term}: {glos[term]}")
    return "\n".join(lines)


def report_catalogo(language=None) -> str:
    language = _idioma(language)
    lines = [t("--- ¿Qué quieres saber? ---", language), ""]
    for m in metas(language):
        lines.append(f"  {m.clave:20s} {m.pregunta}")
    lines += ["",
              t("Se elige con:  olla-dft wizard estructura.cif --goal <clave>",
                language),
              t("O se busca con tus palabras:  olla-dft wizard estructura.cif "
                '--ask "quiero saber si absorbe luz"', language)]
    return "\n".join(lines)


def report_diagnostico(d: Diagnostico, language=None) -> str:
    language = _idioma(language)

    def T(s):
        return t(s, language)

    lines = [T("--- Lo que veo en tu estructura ---"),
             T("  Fórmula: {formula}   ({natoms} átomos)").format(
                 formula=d.formula, natoms=d.natoms),
             T("  Grupo espacial: {grupo}").format(grupo=T(d.grupo)),
             T("  Elementos: {elementos}").format(
                 elementos=", ".join(d.elementos))]
    if d.pseudos_faltan:
        lines += ["",
                  T("FALTAN PSEUDOPOTENCIALES para:") + " "
                  + ", ".join(d.pseudos_faltan),
                  "  " + T(GLOSARIO["pseudopotencial"]),
                  T("  Se descargan de pseudo-dojo.org o de "
                    "quantum-espresso.org/pseudopotentials,"),
                  T("  y se le dicen a Olla-DFT con:  olla-dft config set "
                    "pseudo_dir /ruta/a/tus/pseudos")]
    # las notas se rehacen desde la plantilla para poder traducirlas; si
    # alguien construyó el Diagnostico a mano solo con `notas`, van tal cual
    notas = ([T(tpl).format(**campos) for tpl, campos in d.notas_src]
             if d.notas_src else list(d.notas))
    for n in notas:
        lines += ["", "  " + n.replace("\n", "\n  ")]
    return "\n".join(lines)
