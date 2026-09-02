# Arquitectura

Cómo está organizado el código de Olla-DFT, las reglas que sigue cada módulo y
qué hace falta para añadir un comando.

El comando instalado es `olla-dft`; el paquete Python conserva su nombre
original, `qekit`, para que los scripts, proyectos y archivos de configuración
existentes sigan funcionando.

### Árbol de módulos

```
qekit/
├── __init__.py            # versión, nombre del producto, nombre del comando, autor
├── __main__.py            # python -m qekit
├── cli.py                 # menú interactivo + subcomandos planos (argparse), despacho, códigos de salida
├── config.py              # configuración persistente (config.ini en la carpeta de cada sistema), migración
├── core/                  # infraestructura compartida por todos los módulos
│   ├── atomconf.py        # configuraciones electrónicas atómicas para generar pseudopotenciales
│   ├── compat.py          # compatibilidad entre versiones de las dependencias (numpy, ASE)
│   ├── consola.py         # salida que no revienta: UTF-8 primero, transliteración a ASCII si no
│   ├── errors.py          # ErrorDeUso (error de uso, código 2) frente a falla del programa (código 1)
│   ├── i18n.py            # idioma de la interfaz (es/en): --language, OLLA_DFT_LANG, config; tablas JSON
│   ├── kpoints.py         # mallas k centradas en Γ por espaciado y caminos de alta simetría (seekpath)
│   ├── layers.py          # detección de capas en materiales laminares por conectividad periódica
│   ├── plataforma.py      # carpetas por sistema, nombres de binarios (.x/.exe), lanzador MPI, guiones portables
│   ├── provenance.py      # versión, línea de comandos y parámetros escritos en cada salida
│   ├── pseudo.py          # manejo de UPF: búsqueda, cutoffs y valencia del encabezado
│   ├── qeout.py           # lectura de salidas de Quantum ESPRESSO (XML de pw.x, etiquetas de alta simetría)
│   ├── runner.py          # ejecución por lotes de pw.x (--run): reanudable, con fallos de QE traducidos
│   ├── structure.py       # lectura de estructuras y simetría (ASE + spglib), celdas primitiva/convencional
│   ├── style.py           # estilo de figuras para publicación: tamaños físicos, tipografía, paleta
│   ├── themes.py          # plantillas visuales y paletas verificadas
│   └── wfc.py             # lector de los archivos binarios de función de onda de QE
├── data/
│   ├── atomic_scattering_params.json   # factores de forma de rayos X (de pymatgen, MIT) con su atribución
│   ├── i18n/              # cli_es/en, menu_es/en, onboarding_es/en, dashboard_es/en, recipes_en, wizard_en, docs_en
│   └── theory/            # electronica/mecanica/espectros .es.md/.en.md: la física detrás de cada comando
└── modules/               # un archivo por tarea; cada uno expone su física en el docstring
    ├── adsorb.py          # sitios de adsorción sobre una losa: enumerarlos, montarlos y compararlos
    ├── align.py           # alineamiento de bandas: dónde queda el VBM de uno respecto al del otro
    ├── amorphous.py       # sólidos amorfos por fundido y temple (MACE)
    ├── audit.py           # auditoría de consistencia entre cálculos, y base de datos local
    ├── ballistic.py       # transporte balístico: conductancia de Landauer con pwcond.x
    ├── bands.py           # post-proceso de bandas, análisis de gap, exportación y gráfica
    ├── berry.py           # polarización eléctrica por fase de Berry, y lo que se deduce de ella (Z*)
    ├── builder.py         # constructores de estructuras: superficies por índices de Miller y defectos puntuales
    ├── campaign.py        # campañas reproducibles de cálculos parametrizados
    ├── charges.py         # cargas atómicas (Bader, Löwdin) y diferencia de densidad
    ├── combined.py        # figura combinada bandas + DOS
    ├── compare.py         # comparación segura de corridas de Quantum ESPRESSO
    ├── converge.py        # pruebas de convergencia: ecutwfc, ecutrho y malla de puntos k
    ├── corehole.py        # pseudopotenciales con hueco de core con ld1.x, para XPS y XANES
    ├── cost.py            # cuánto va a tardar esto, antes de lanzarlo
    ├── crosscheck.py      # validación cruzada: la misma cantidad por caminos independientes
    ├── dashboard.py       # dashboard HTML autocontenido de un proyecto
    ├── datasheet.py       # ficha del material y párrafo de métodos
    ├── defects.py         # energía de formación de defectos cargados (corrección de Madelung)
    ├── derived.py         # cantidades derivadas de las constantes elásticas y de los fonones: Debye, velocidades del sonido, Slack
    ├── diagnose.py        # diagnóstico de un cálculo de pw.x: ¿sirve, y si no, por qué?
    ├── docs.py            # referencia HTML navegable, generada del propio código
    ├── dos.py             # post-proceso de DOS y DOS proyectada
    ├── dynamics.py        # trayectorias de dinámica molecular: g(r), MSD/difusión, VDOS
    ├── echem.py           # electrodo de hidrógeno computacional: HER, OER y el diagrama de Pourbaix
    ├── effmass.py         # masa efectiva por ajuste parabólico de las bandas
    ├── elastic.py         # constantes elásticas por el método de esfuerzo-deformación
    ├── elph.py            # acoplamiento electrón-fonón: lambda, alpha²F, Tc y un tau de verdad
    ├── environment.py     # bloqueo ligero del entorno para reproducibilidad local
    ├── eos.py             # ecuación de estado E-V: volumen de equilibrio y módulo de bulk
    ├── esm.py             # superficies cargadas: el medio de apantallamiento efectivo (ESM)
    ├── exfoliate.py       # energía de exfoliación de un material laminar
    ├── feedback.py        # registro local de incidencias: fallas, confusiones y errores colados; sin telemetría
    ├── fields.py          # densidad de carga, potencial electrostático y función trabajo (pp.x)
    ├── health.py          # diagnóstico de la instalación y del entorno de ejecución
    ├── hubbard.py         # U de Hubbard por respuesta lineal, con hp.x
    ├── inputgen.py        # generador de inputs para pw.x y post-proceso; run.sh y run.py
    ├── interface.py       # heteroestructuras: apilar dos materiales con la menor deformación posible
    ├── interop.py         # formato de intercambio con el resto de la suite
    ├── kappa.py           # conductividad térmica de red: el fonón que se dispersa contra otro fonón (phono3py)
    ├── mlip.py            # potenciales interatómicos aprendidos: pre-relajación y cribado
    ├── neb.py             # caminos de reacción y barreras de activación con neb.x
    ├── onboarding.py      # inicio guiado para quien nunca ha usado una CLI científica
    ├── optics.py          # propiedades ópticas con epsilon.x: ε(ω), n, k, absorción, Tauc
    ├── phonons.py         # fonones por DFPT: dispersión, DOS, termodinámica e IR
    ├── project.py         # Project Hub: proyectos, workflows reanudables y procedencia
    ├── pseudos.py         # elegir pseudopotencial con criterio, no por orden alfabético
    ├── qha.py             # aproximación cuasi-armónica: expansión térmica y a(T)
    ├── quality.py         # puerta de calidad científica de un proyecto
    ├── recipes.py         # recetas: sesiones completas, de la estructura al resultado
    ├── recommend.py       # recomendaciones a partir de TU propio historial de cálculos
    ├── report.py          # informe PDF compacto de un proyecto
    ├── results.py         # motor local de resultados normalizados y trazables
    ├── selftest.py        # comprobación contra la física conocida, no contra uno mismo
    ├── strain.py          # barrido de deformación: propiedades en función de la deformación aplicada
    ├── surfen.py          # energía de superficie: cortar un cristal y ver cuánto cuesta
    ├── sweep.py           # infraestructura común de los barridos (prepare / run / collect)
    ├── tddft.py           # absorción óptica con TDDFPT: la parte que epsilon.x no ve
    ├── theory.py          # fundamento físico de cada comando, consultable desde la terminal (olla-dft teoria)
    ├── thermo.py          # energías de formación, casco convexo y estabilidad de fases
    ├── thermochem.py      # termoquímica: de una energía DFT a una energía libre comparable
    ├── topology.py        # invariantes topológicos de un Hamiltoniano de Wannier: Chern, lazos de Wilson
    ├── tphonons.py        # fonones a temperatura electrónica: ¿se estabiliza el modo imaginario?
    ├── transport.py       # transporte electrónico en aproximación de tiempo de relajación constante
    ├── tuning.py          # recomendación adaptativa a partir de una serie de convergencia
    ├── uncertainty.py     # utilidades pequeñas para propagar incertidumbres experimentales o numéricas
    ├── unfold.py          # desdoblamiento de bandas: recuperar la dispersión de una supercelda
    ├── validation.py      # validaciones estructurales y de integridad antes de gastar CPU
    ├── wannier.py         # funciones de Wannier: bajar la estructura de bandas a un modelo pequeño
    ├── wizard.py          # asistente guiado: de lo que quieres SABER a los archivos que hay que correr
    ├── xanes.py           # XANES / NEXAFS con xspectra.x
    ├── xps.py             # corrimientos de niveles de core (XPS) en aproximación de estado inicial
    └── xrd.py             # difracción de polvos simulada a partir de la estructura cristalina
```

`tools/build_docs.py` genera `docs/COMANDOS.md` y `docs/TEORIA.md` a partir
del árbol de argparse y de `qekit/data/theory/`; esos dos archivos nunca se
editan a mano (`--all` escribe también el par en inglés).

### Reglas de diseño

**Subcomandos planos, ayuda agrupada.** Cada comando es un `olla-dft
<comando>` plano (fácil de automatizar), pero `olla-dft --help` y el menú
interactivo los presentan por tareas: la tabla `COMMAND_GROUPS` de `cli.py`
lista cada comando exactamente una vez y una prueba la protege contra olvidos.
Tres comandos con nombre en español tienen alias en inglés
(`sistema`/`system`, `recetas`/`recipes`, `teoria`/`theory`).

**prepare / --run / --collect.** Todo barrido (convergencia, EOS, elásticas,
deformación, energía de superficie, adsorción, defectos, fonones, ESM, Berry,
kappa…) sigue el mismo ciclo, compartido a través de `modules/sweep.py`:
`prepare` escribe una carpeta por cálculo más `run.sh` y `run.py`; `--run` los
ejecuta aquí con `core/runner.py` (reanudable, `-j N` puntos en paralelo,
presupuesto `--max-time`, `--timeout`, `--redo`); sin `--run` el comando
explica cómo lanzarlos; `--collect` lee los cálculos terminados y escribe el
informe sin reescribir los inputs, de modo que un cálculo que el usuario editó
a mano se describe tal como corrió. `--estimate` predice el coste con el
histórico local y sale. Todos esos comandos comparten los mismos grupos de
opciones (`ejecución`, `parámetros DFT`, `figura`), que argparse lista después
de las opciones propias del comando.

**Códigos de salida.** `core/errors.py` separa los errores por a quién le toca
arreglarlos. `ErrorDeUso` (subclase de `ValueError`) significa que el comando o
sus datos no encajan y el mensaje ya dice qué hacer: el programa hizo lo
correcto, sale con **2** (como argparse) y se anota como tipo "uso", sin traza.
Cualquier otra excepción es una falla del programa: sale con **1** y se archiva
completa —comando, traza, versiones— con `modules/feedback.py`. Los comandos
que corrieron bien y *encontraron* un problema (`doctor`, `crosscheck`,
`audit`, `selftest`, `project quality`, una `campaign` con una tarea fallida)
devuelven también **1**, para que un script
pueda detenerse en ellos; `tests/barrido_cli.sh` declara el código esperado en
cada línea. `Ctrl-C` devuelve 130; una tubería rota (`| head`) se cierra en
silencio con 0.

**Procedencia en todo.** `core/provenance.py` escribe la versión de Olla-DFT,
la fecha UTC, la línea de comandos exacta y los parámetros como comentarios
`#` al principio de cada `.dat`, y como metadatos en cada PDF/PNG (legibles con
`pdfinfo`, `exiftool` u `olla-dft info --figura`). Meses después, una cifra de
una tesis se puede rastrear hasta el cálculo que la produjo.

**Capa de idioma de la interfaz.** `core/i18n.py` decide el idioma en este
orden: la bandera global `--language en` (aceptada en cualquier posición), la
variable `OLLA_DFT_LANG`, la clave `language` de la configuración, español. Lo
que se traduce es la *interfaz*: la ayuda de cada comando y bandera
(`data/i18n/cli_en.json`), el menú interactivo (`menu_*.json`), el inicio guiado
(`onboarding_*.json`), las recetas y el asistente (`recipes_en`, `wizard_en`,
aplicados al catálogo en español con `translate_data`), el dashboard
(`dashboard_*.json`), la referencia HTML (`docs_en.json`) y la teoría
(`data/theory/*.en.md`). Lo que sigue en español es el informe científico que
cada comando imprime y escribe: los docstrings de los módulos, los
encabezados de los `.dat` y el texto de análisis no se traducen en tiempo de
ejecución.

**Salida de consola que no muere.** `core/consola.py` pone la consola en UTF-8
cuando puede y, si no, translitera cada símbolo fuera de ASCII a un
equivalente que se lee igual (`Å → A`, `α → alpha`, `→ → ->`); `--ascii` lo
fuerza. Una prueba recorre todos los informes reales y falla si a la tabla le
falta un símbolo.

**Sin telemetría.** `modules/feedback.py` guarda el registro de incidencias en
la carpeta local de configuración; nada sale de la máquina salvo que el
usuario lo exporte con `olla-dft report --export`.

### Añadir un comando

1. Escribe la física en un módulo nuevo `qekit/modules/<nombre>.py`. Su
   docstring no es de cortesía: `olla-dft docs` y la referencia HTML lo citan
   como "la física detrás" del comando, así que debe explicar qué responde y
   cuándo el resultado no vale. Un barrido expone `prepare`, `collect`,
   `report`, `export` y `plot`, y construye sus trabajos con
   `sweep.prepare_common` / `sweep.write_scf_job`.
2. En `qekit/cli.py`: añade el parser (dentro de `build_parser`, con
   `_calc_opts` si es un barrido, para que herede
   `--run/--collect/-j/--pseudo-dir…`), escribe `_cmd_<nombre>(args) -> int`
   (lanza `ErrorDeUso` ante datos que no encajan; usa `_run_or_explain` para
   el paso de ejecución) y regístralo en `_DISPATCH` y en `COMMAND_GROUPS`.
   Toda opción necesita su `help=` en español.
3. En `qekit/modules/docs.py`: añade el comando a un grupo de `GRUPOS` y, si su
   física vive en un módulo con otro nombre, a `MODULO_DE`
   (`tests/test_docs.py` falla con huérfanos).
4. En `qekit/data/i18n/docs_en.json`: añade el resumen de una línea en
   `command_summaries`. En `qekit/data/i18n/cli_en.json`: añade la traducción
   al inglés de cada texto de ayuda nuevo en `help` (la clave es el texto en
   español). Los catálogos en español (`recetas`, `wizard`) llevan sus
   entradas `{es: en}` en el `_en.json` correspondiente.
5. En `qekit/data/theory/<area>.es.md` y `<area>.en.md`: añade una sección
   encabezada ``### `olla-dft <nombre>` — título`` con los apartados
   obligatorios (*Qué responde*, *Fundamento para no expertos*, *Cómo lo
   calcula Olla-DFT*, *Límites y trampas*, *Referencias*, y sus equivalentes en
   inglés); `tests/test_teoria.py` exige una para cada comando científico y
   comprueba la paridad es/en.
6. Añade una prueba en `tests/` (con una salida real de QE en `tests/datos/`
   si el comando lee una, y una referencia congelada en `tests/referencias.py`
   si produce un número validado contra experimento).
7. Corre `python tools/build_docs.py` para regenerar `docs/COMANDOS.md`
   y `docs/TEORIA.md`, y después
   `python -m pytest -q` y `python -m pyflakes qekit tests`.

### Pruebas

- `tests/` — 977 pruebas de pytest que corren sin Quantum ESPRESSO en
  menos de un minuto (`python -m pytest -q`). Cubren las funciones, el
  árbol de argparse (cada comando en un grupo, cada comando con sección de
  teoría, cada comando de README y de receta válido), las tres plataformas
  (`sys.platform` simulado, salida forzada a cp1252), las tablas de i18n (la
  tabla inglesa cubre todo el catálogo) y la licencia.
- `tests/datos/` — salidas reales de Quantum ESPRESSO que leen las pruebas:
  bandas, DOS y masa efectiva del Si, fonones del Si, espectros de
  `epsilon.x`, XANES, TDDFPT del etileno, `pwcond.x` de un hilo de Al,
  electrón-fonón del Al, una trayectoria de MD, un archivo de desdoblamiento y
  un UPF con hueco de core.
- `tests/referencias.py` — valores de referencia congelados, cada uno validado
  una vez contra experimento, una ficha PDF u otra implementación; desde
  entonces son detectores de regresiones, no documentación.
- `tests/barrido_cli.sh` — barrido de regresión a nivel de comando sobre
  salidas de QE ya calculadas (`OLLA_DFT_REG=/ruta bash tests/barrido_cli.sh`),
  con los `RuntimeWarning` de numpy convertidos en errores para que un NaN
  silencioso no pase, y el código de salida esperado (0, 1 o 2) declarado en
  cada línea.
- Las pruebas que necesitan binarios de QE llevan la marca `qe`.
