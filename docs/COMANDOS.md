# Referencia de comandos de Olla-DFT

Los 80 subcomandos de `olla-dft`, agrupados por área, con sus opciones. Generado del propio código con `python tools/build_docs.py`; la misma información sale en la terminal con `olla-dft COMANDO --help` y, navegable, con `olla-dft docs`.

## Índice

- **Primeros pasos**: [`start`](#start), [`wizard`](#wizard), [`recetas`](#recetas), [`teoria`](#teoria), [`docs`](#docs), [`sistema`](#sistema), [`selftest`](#selftest), [`update`](#update)
- **Estructuras e inputs**: [`gen`](#gen), [`info`](#info), [`kpath`](#kpath), [`prim`](#prim), [`conv`](#conv), [`supercell`](#supercell), [`convert`](#convert)
- **Estructura electrónica**: [`bands`](#bands), [`dos`](#dos), [`plot`](#plot), [`gap`](#gap), [`fermi`](#fermi), [`effmass`](#effmass), [`wannier`](#wannier), [`unfold`](#unfold), [`topology`](#topology), [`hubbard`](#hubbard)
- **Espectros y respuesta**: [`optics`](#optics), [`tddft`](#tddft), [`xanes`](#xanes), [`xps`](#xps), [`corehole`](#corehole), [`charge`](#charge), [`charges`](#charges), [`wf`](#wf), [`berry`](#berry)
- **Fonones, transporte y temperatura**: [`phonons`](#phonons), [`elph`](#elph), [`transport`](#transport), [`ballistic`](#ballistic), [`kappa`](#kappa), [`qha`](#qha), [`thermochem`](#thermochem), [`md`](#md), [`derived`](#derived)
- **Mecánica y estabilidad**: [`converge`](#converge), [`eos`](#eos), [`elastic`](#elastic), [`strain`](#strain), [`layers`](#layers), [`xrd`](#xrd), [`exfoliate`](#exfoliate), [`gamma`](#gamma)
- **Superficies, defectos y química**: [`surface`](#surface), [`defect`](#defect), [`interface`](#interface), [`adsorb`](#adsorb), [`eform`](#eform), [`align`](#align), [`esm`](#esm), [`echem`](#echem), [`neb`](#neb), [`amorphous`](#amorphous)
- **Automatización y calidad**: [`doctor`](#doctor), [`audit`](#audit), [`crosscheck`](#crosscheck), [`cost`](#cost), [`db`](#db), [`hull`](#hull), [`mlip`](#mlip), [`suggest`](#suggest), [`datasheet`](#datasheet), [`report`](#report), [`compare`](#compare), [`tune`](#tune), [`results`](#results), [`campaign`](#campaign), [`pseudos`](#pseudos)
- **Proyecto**: [`project`](#project), [`resilient`](#resilient)
- **Apariencia y configuración**: [`templates`](#templates), [`config`](#config)

## Primeros pasos

### `start`

inicio guiado para crear un proyecto sin conocer la CLI

**Uso:** `olla-dft start [-h] [--project PROJECT] [--structure STRUCTURE] [--goal GOAL] [--name NAME] [--non-interactive] [--no-validate] [--language {es,en}]`

**Opciones:**

| Opción | Descripción |
|---|---|
| `--project` | carpeta del proyecto (default: `.`) |
| `--structure` | CIF, POSCAR o input de pw.x |
| `--goal` | relax, gap, dos, phonons, optics o scf |
| `--name` | nombre visible del proyecto |
| `--non-interactive` | no preguntar; requiere --structure en un proyecto nuevo |
| `--no-validate` | no ejecutar la validación inicial |
| `--language {es,en}` | idioma del inicio guiado (default: es) |

### `wizard`

asistente: dime QUE quieres saber y te digo que hay que correr, en orden y con los comandos

**Uso:** `olla-dft wizard [-h] [--goal GOAL] [--ask TEXTO] [--list] [--term TERM] [--no-glossary] [--pseudo-dir PSEUDO_DIR] [--pseudo EL=UPF] [file]`

**Argumentos:**

- `file` — tu estructura (opcional)

**Opciones:**

| Opción | Descripción |
|---|---|
| `--goal` | clave de la meta; salen con --list |
| `--ask TEXTO` | describelo con tus palabras, por ejemplo 'quiero saber si absorbe luz' |
| `--list` | listar todo lo que el asistente sabe hacer |
| `--term` | que significa un termino |
| `--no-glossary` | no explicar los términos técnicos al final de la respuesta |
| `--pseudo-dir` | carpeta con los pseudopotenciales UPF (si no se da, la de 'olla-dft config') |
| `--pseudo EL=UPF` | forzar un pseudopotencial concreto, por ejemplo Fe=Fe.rel-pbe.UPF. Se puede repetir. Sin esto, Olla-DFT elige con 'olla-dft pseudos' |

### `recetas`

sesiones completas de principio a fin: qué comando va después de cuál y qué archivo se pasan entre ellos

**Uso:** `olla-dft recetas [-h] [--buscar TEXTO] [--script [ARCHIVO]] [receta]`

**Argumentos:**

- `receta` — clave de la receta; sin nada, las lista todas

**Opciones:**

| Opción | Descripción |
|---|---|
| `--buscar TEXTO` | buscarla con tus palabras, sin saber la clave |
| `--script ARCHIVO` | escribir la receta como un guion de shell comentado, listo para editar |

### `teoria`

el fundamento físico de un comando: qué responde, las fórmulas que implementa, de qué módulo salen y de dónde sale cada dato

**Uso:** `olla-dft teoria [-h] [--all] [-o ARCHIVO.md] [comando]`

**Argumentos:**

- `comando` — comando a explicar; sin nada, el índice

**Opciones:**

| Opción | Descripción |
|---|---|
| `--all` | el documento completo (todas las áreas) |
| `-o, --output ARCHIVO.md` | guardarlo en Markdown en vez de imprimirlo |

### `docs`

referencia navegable de todos los subcomandos, generada del propio código

**Uso:** `olla-dft docs [-h] [-o OUTPUT] [--open] [--language {es,en}] [--both]`

**Opciones:**

| Opción | Descripción |
|---|---|
| `-o, --output` | archivo HTML de salida (default: `olla-dft-docs.html`) |
| `--open` | abrirla en el navegador al terminar |
| `--language {es,en}` | idioma de la interfaz de referencia (default: es) |
| `--both` | generar referencias en español e inglés por separado |

### `sistema`

qué ve Olla-DFT de esta máquina: codificación, dónde guarda la configuración, qué binarios de QE encuentra y cómo lanzar los cálculos aquí

**Uso:** `olla-dft sistema [-h]`

### `selftest`

comprobar Olla-DFT contra valores publicados, no contra sí mismo

**Uso:** `olla-dft selftest [-h] [--full] [--mlip] [--only ONLY] [--list] [--pseudo-dir PSEUDO_DIR] [--pw-cmd PW_CMD] [--nproc NPROC] [-j JOBS] [--keep CARPETA]`

**Opciones:**

| Opción | Descripción |
|---|---|
| `--full` | incluir las pruebas que corren pw.x de verdad (unos diez minutos) |
| `--mlip` | incluir por separado la prueba con potencial aprendido (requiere MACE) |
| `--only` | solo estas pruebas, separadas por coma |
| `--list` | listar las pruebas y sus referencias, sin correr nada |
| `--pseudo-dir` | pseudopotenciales para las de --full |
| `--pw-cmd` | ejecutable de pw.x para --run; de su ruta salen los demás binarios de QE |
| `--nproc NPROC` | número de procesos MPI para los cálculos que se lanzan con --run |
| `-j, --jobs JOBS` | pruebas simultáneas (default: 1) |
| `--keep CARPETA` | dejar los cálculos aquí en vez de borrarlos |

**Fundamento físico:** [`olla-dft teoria selftest`](TEORIA.md)

### `update`

comprobar si hay una versión nueva de Olla-DFT y, si la hay, instalarla con una confirmación; nunca se ejecuta solo

**Uso:** `olla-dft update [-h] [--check] [--yes] [--version TAG]`

**Opciones:**

| Opción | Descripción |
|---|---|
| `--check` | solo comprobar e informar, sin instalar nada |
| `--yes` | no preguntar; instalar directamente si hay versión nueva |
| `--version TAG` | instalar una versión concreta (p. ej. v1.0.1) en vez de la última |

## Estructuras e inputs

### `gen`

generar inputs de pw.x y post-proceso

**Uso:** `olla-dft gen [-h] [-p {scf,relax,vc-relax,nscf,bands,dos,all,md}] [-o OUTDIR] [-k {coarse,fine,gamma,medium,very-fine}] [--kspacing KSPACING] [--kgrid N N N] [--band-points BAND_POINTS] [--ecutwfc ECUTWFC] [--ecutrho ECUTRHO] [--insulator] [--primitive] [--pseudo-dir PSEUDO_DIR] [--pseudo EL=UPF] [--prefix PREFIX] [--nspin {1,2}] [--mag MAG] [--vdw {grimme-d2,grimme-d3,DFT-D,ts-vdw,xdm,mbd}] [--soc] [--hubbard EL=U] [--hubbard-style {legacy,card}] [--charge Q] [--dipole [EJE]] [--nosym] [--functional {b3lyp,gaupbe,hse,pbe0}] [--exx-grid NxNxN] [--exx-fraction EXX_FRACTION] [--dt FS] [--nstep NSTEP] [--thermostat {none,rescaling,berendsen,andersen,initial,reduce-history}] [-T TEMPERATURE] file`

**Argumentos:**

- `file` — estructura (CIF, POSCAR, entrada de pw.x, ...)

**Opciones:**

| Opción | Descripción |
|---|---|
| `-p, --preset {scf,relax,vc-relax,nscf,bands,dos,all,md}` | tipo de cálculo (default: scf) |
| `-o, --outdir` | carpeta de salida (default: `.`) |
| `-k, --klevel {coarse,fine,gamma,medium,very-fine}` | densidad de la malla k (gamma/coarse/medium/fine/very-fine) |
| `--kspacing KSPACING` | espaciado k en Å^-1 (anula --klevel) |
| `--kgrid N` | malla k explícita para scf/relax (tres enteros; anula --kspacing y --klevel) |
| `--band-points BAND_POINTS` | puntos por segmento del k-path |
| `--ecutwfc ECUTWFC` | cutoff de funciones de onda (Ry) |
| `--ecutrho ECUTRHO` | cutoff de densidad (Ry) |
| `--insulator` | occupations='fixed' (aislantes; default: smearing) |
| `--primitive` | reducir a la celda primitiva estandarizada antes de generar |
| `--pseudo-dir` | carpeta de pseudopotenciales (anula config) |
| `--pseudo EL=UPF` | forzar un pseudopotencial concreto, por ejemplo Fe=Fe.rel-pbe.UPF. Se puede repetir. Sin esto, Olla-DFT elige con 'olla-dft pseudos' |
| `--prefix` | prefix del cálculo (default: fórmula) |
| `--nspin {1,2}` | 2 activa la polarización de espín (default: `1`) |
| `--mag` | magnetización inicial: un número (0.5) o por elemento (Fe=0.7,O=0). Implica --nspin 2 |
| `--vdw {grimme-d2,grimme-d3,DFT-D,ts-vdw,xdm,mbd}` | corrección de dispersión (van der Waals) |
| `--soc` | acoplamiento espín-órbita: cálculo no colineal con lspinorb (exige pseudos totalmente relativistas) |
| `--hubbard EL=U` | U de Hubbard en eV por elemento, por ejemplo Ni=4.1. Se puede repetir. Para calcularlo en vez de suponerlo:  olla-dft hubbard --cycle |
| `--hubbard-style {legacy,card}` | legacy = lda_plus_u (QE <= 7.0), card = tarjeta HUBBARD (QE >= 7.1) (default: `legacy`) |
| `--charge Q` | carga total de la celda (tot_charge): +1 le quita un electrón, -1 se lo añade |
| `--dipole EJE` | corrección dipolar para losas polares; sin valor usa el eje c. Coloca el diente de sierra dentro del vacío |
| `--nosym` | desactivar la simetría (nosym y noinv) |
| `--functional {b3lyp,gaupbe,hse,pbe0}` | funcional híbrido: hse, pbe0, b3lyp o gaupbe. Cuesta entre uno y dos órdenes de magnitud más que PBE, y el reporte lo dice con números |
| `--exx-grid NxNxN` | malla q del intercambio exacto (default 1x1x1). Tiene que dividir la malla de k |
| `--exx-fraction EXX_FRACTION` | fracción de intercambio exacto, si quieres cambiar la del funcional |
| `--dt FS` | paso de tiempo de la MD en fs (default: 1.0) |
| `--nstep NSTEP` | pasos de la MD (default: 1000) |
| `--thermostat {none,rescaling,berendsen,andersen,initial,reduce-history}` | termostato de la MD; none = NVE (default) |
| `-T, --temperature TEMPERATURE` | temperatura objetivo de la MD en K (default: 300) |

**Fundamento físico:** [`olla-dft teoria gen`](TEORIA.md)

### `info`

información de estructura y simetría

**Uso:** `olla-dft info [-h] file`

**Argumentos:**

- `file` — estructura de entrada (CIF, POSCAR, input de pw.x...)

**Fundamento físico:** [`olla-dft teoria info`](TEORIA.md)

### `kpath`

camino de alta simetría (seekpath)

**Uso:** `olla-dft kpath [-h] file`

**Argumentos:**

- `file` — estructura de entrada (CIF, POSCAR, input de pw.x...)

**Fundamento físico:** [`olla-dft teoria kpath`](TEORIA.md)

### `prim`

celda primitiva estandarizada

**Uso:** `olla-dft prim [-h] [-o OUTPUT] file`

**Argumentos:**

- `file` — estructura de entrada (CIF, POSCAR, input de pw.x...)

**Opciones:**

| Opción | Descripción |
|---|---|
| `-o, --output` | archivo de estructura de salida (por omisión, un .cif con el nombre del comando) (default: `primitive.cif`) |

**Fundamento físico:** [`olla-dft teoria prim`](TEORIA.md)

### `conv`

celda convencional estandarizada

**Uso:** `olla-dft conv [-h] [-o OUTPUT] file`

**Argumentos:**

- `file` — estructura de entrada (CIF, POSCAR, input de pw.x...)

**Opciones:**

| Opción | Descripción |
|---|---|
| `-o, --output` | archivo de estructura de salida (por omisión, un .cif con el nombre del comando) (default: `conventional.cif`) |

**Fundamento físico:** [`olla-dft teoria conv`](TEORIA.md)

### `supercell`

construir supercelda

**Uso:** `olla-dft supercell [-h] [-o OUTPUT] file nx ny nz`

**Argumentos:**

- `file` — estructura de entrada (CIF, POSCAR, input de pw.x...)
- `nx` — repeticiones de la celda a lo largo de a
- `ny` — repeticiones de la celda a lo largo de b
- `nz` — repeticiones de la celda a lo largo de c

**Opciones:**

| Opción | Descripción |
|---|---|
| `-o, --output` | archivo de estructura de salida (por omisión, un .cif con el nombre del comando) (default: `supercell.cif`) |

**Fundamento físico:** [`olla-dft teoria supercell`](TEORIA.md)

### `convert`

convertir formato (CIF/POSCAR/XYZ)

**Uso:** `olla-dft convert [-h] [-o OUTPUT_FLAG] file [output]`

**Argumentos:**

- `file` — estructura de entrada (CIF, POSCAR, input de pw.x...)
- `output` — archivo de destino; el formato se deduce de la extensión (.cif, .vasp, .xyz...)

**Opciones:**

| Opción | Descripción |
|---|---|
| `-o, --output-flag` | archivo de salida (alternativa a darlo posicional) |

## Estructura electrónica

### `bands`

analizar y graficar la estructura de bandas

**Uso:** `olla-dft bands [-h] [-o OUTDIR] [--prefix PREFIX] [--ref {auto,fermi,vbm,none}] [--emin EMIN] [--emax EMAX] [--no-plot] [--dpi DPI] [--format FORMAT] [-t TEMPLATE] [--size {paper,poster,presentation}] [--font {sans,serif,latex}] [--usetex] [--palette PALETTE] [--background BACKGROUND] [--journal {acs,aps,elsevier,generic,iop,nature,rsc,wiley}] [--width WIDTH] [--aspect ASPECT] [--mono] [--dashes {auto,always,never}] [--title TITLE] [--gap-label] [--panel PANEL] [--fat SELECTOR] [--fat-scale FAT_SCALE] [--projwfc ARCHIVO] [path]`

**Argumentos:**

- `path` — carpeta del cálculo (o ruta al .xml)

**Opciones:**

| Opción | Descripción |
|---|---|
| `-o, --outdir` | carpeta de salida (default: `.`) |
| `--prefix` | prefix del cálculo (se detecta solo) |
| `--ref {auto,fermi,vbm,none}` | origen de energías (default: auto) |
| `--emin EMIN` | límite inferior del eje de energía (eV) (default: `-6.0`) |
| `--emax EMAX` | límite superior del eje de energía (eV) (default: `6.0`) |
| `--no-plot` | solo exportar datos, sin generar la gráfica |
| `--dpi DPI` | resolución de los formatos de mapa de bits (default: `600`) |
| `--format` | formatos separados por coma: pdf,png,svg,eps,tif (default: `pdf,png`) |
| `-t, --template` | plantilla visual: dark, journal, latex, latex-true, minimal, mono, mono-latex, poster, slides (o la ruta a un JSON propio) |
| `--size {paper,poster,presentation}` | escala tipográfica: paper / presentation / poster |
| `--font {sans,serif,latex}` | familia tipográfica (latex = Computer Modern) |
| `--usetex` | renderizar el texto con LaTeX de verdad |
| `--palette` | paleta: grayscale, okabe-ito, okabe-ito-dark, o colores hexadecimales separados por coma |
| `--background` | color de fondo, por ejemplo '#FFFFFF' o 'none' |
| `--journal {acs,aps,elsevier,generic,iop,nature,rsc,wiley}` | anchos de columna de la editorial (default: `generic`) |
| `--width` | ancho: single / onehalf / double, o un número en mm |
| `--aspect ASPECT` | relación alto/ancho de la figura |
| `--mono` | monocromo: tinta negra y patrones de línea (para revistas que cobran el color) |
| `--dashes {auto,always,never}` | patrones de línea como codificación secundaria (default: `auto`) |
| `--title` | título dentro de la figura (por defecto ninguno: en un artículo el texto va en el pie de figura) |
| `--gap-label` | anotar el valor del gap dentro de la gráfica |
| `--panel` | etiqueta de panel, por ejemplo '(a)' |
| `--fat SELECTOR` | fatbands: peso de un orbital sobre cada banda. Por ejemplo Ni-d, Si-p, O, d o atomo:3. Necesita la salida de projwfc.x del MISMO cálculo de bandas |
| `--fat-scale FAT_SCALE` | tamaño de los puntos de las fatbands (default: `55.0`) |
| `--projwfc ARCHIVO` | salida de projwfc.x (por omisión projwfc.out en la misma carpeta) |

**Fundamento físico:** [`olla-dft teoria bands`](TEORIA.md)

### `dos`

analizar y graficar DOS y PDOS

**Uso:** `olla-dft dos [-h] [--mode {orbital,element,total}] [-o OUTDIR] [--prefix PREFIX] [--ref {auto,fermi,vbm,none}] [--emin EMIN] [--emax EMAX] [--no-plot] [--dpi DPI] [--format FORMAT] [-t TEMPLATE] [--size {paper,poster,presentation}] [--font {sans,serif,latex}] [--usetex] [--palette PALETTE] [--background BACKGROUND] [--journal {acs,aps,elsevier,generic,iop,nature,rsc,wiley}] [--width WIDTH] [--aspect ASPECT] [--mono] [--dashes {auto,always,never}] [--title TITLE] [--gap-label] [--panel PANEL] [--dband EL[-ORB]] [--dband-emax eV] [path]`

**Argumentos:**

- `path` — carpeta del cálculo (o ruta al .xml)

**Opciones:**

| Opción | Descripción |
|---|---|
| `--mode {orbital,element,total}` | cómo descomponer la PDOS (default: `orbital`) |
| `-o, --outdir` | carpeta de salida (default: `.`) |
| `--prefix` | prefix del cálculo (se detecta solo) |
| `--ref {auto,fermi,vbm,none}` | origen de energías (default: auto) |
| `--emin EMIN` | límite inferior del eje de energía (eV) (default: `-6.0`) |
| `--emax EMAX` | límite superior del eje de energía (eV) (default: `6.0`) |
| `--no-plot` | solo exportar datos, sin generar la gráfica |
| `--dpi DPI` | resolución de los formatos de mapa de bits (default: `600`) |
| `--format` | formatos separados por coma: pdf,png,svg,eps,tif (default: `pdf,png`) |
| `-t, --template` | plantilla visual: dark, journal, latex, latex-true, minimal, mono, mono-latex, poster, slides (o la ruta a un JSON propio) |
| `--size {paper,poster,presentation}` | escala tipográfica: paper / presentation / poster |
| `--font {sans,serif,latex}` | familia tipográfica (latex = Computer Modern) |
| `--usetex` | renderizar el texto con LaTeX de verdad |
| `--palette` | paleta: grayscale, okabe-ito, okabe-ito-dark, o colores hexadecimales separados por coma |
| `--background` | color de fondo, por ejemplo '#FFFFFF' o 'none' |
| `--journal {acs,aps,elsevier,generic,iop,nature,rsc,wiley}` | anchos de columna de la editorial (default: `generic`) |
| `--width` | ancho: single / onehalf / double, o un número en mm |
| `--aspect ASPECT` | relación alto/ancho de la figura |
| `--mono` | monocromo: tinta negra y patrones de línea (para revistas que cobran el color) |
| `--dashes {auto,always,never}` | patrones de línea como codificación secundaria (default: `auto`) |
| `--title` | título dentro de la figura (por defecto ninguno: en un artículo el texto va en el pie de figura) |
| `--gap-label` | anotar el valor del gap dentro de la gráfica |
| `--panel` | etiqueta de panel, por ejemplo '(a)' |
| `--dband EL[-ORB]` | centro, anchura y llenado de una banda proyectada, por ejemplo Pt (usa d) o Ni-p. Es el descriptor que se correlaciona con la energía de adsorción |
| `--dband-emax eV` | corte superior de la integral, respecto al Fermi |

**Fundamento físico:** [`olla-dft teoria dos`](TEORIA.md)

### `plot`

gráfica combinada de bandas + DOS

**Uso:** `olla-dft plot [-h] [--mode {orbital,element,total}] [-o OUTDIR] [--prefix PREFIX] [--ref {auto,fermi,vbm,none}] [--emin EMIN] [--emax EMAX] [--no-plot] [--dpi DPI] [--format FORMAT] [-t TEMPLATE] [--size {paper,poster,presentation}] [--font {sans,serif,latex}] [--usetex] [--palette PALETTE] [--background BACKGROUND] [--journal {acs,aps,elsevier,generic,iop,nature,rsc,wiley}] [--width WIDTH] [--aspect ASPECT] [--mono] [--dashes {auto,always,never}] [--title TITLE] [--gap-label] [--panel PANEL] [path]`

**Argumentos:**

- `path` — carpeta del cálculo (o ruta al .xml)

**Opciones:**

| Opción | Descripción |
|---|---|
| `--mode {orbital,element,total}` | cómo descomponer la PDOS (default: `orbital`) |
| `-o, --outdir` | carpeta de salida (default: `.`) |
| `--prefix` | prefix del cálculo (se detecta solo) |
| `--ref {auto,fermi,vbm,none}` | origen de energías (default: auto) |
| `--emin EMIN` | límite inferior del eje de energía (eV) (default: `-6.0`) |
| `--emax EMAX` | límite superior del eje de energía (eV) (default: `6.0`) |
| `--no-plot` | solo exportar datos, sin generar la gráfica |
| `--dpi DPI` | resolución de los formatos de mapa de bits (default: `600`) |
| `--format` | formatos separados por coma: pdf,png,svg,eps,tif (default: `pdf,png`) |
| `-t, --template` | plantilla visual: dark, journal, latex, latex-true, minimal, mono, mono-latex, poster, slides (o la ruta a un JSON propio) |
| `--size {paper,poster,presentation}` | escala tipográfica: paper / presentation / poster |
| `--font {sans,serif,latex}` | familia tipográfica (latex = Computer Modern) |
| `--usetex` | renderizar el texto con LaTeX de verdad |
| `--palette` | paleta: grayscale, okabe-ito, okabe-ito-dark, o colores hexadecimales separados por coma |
| `--background` | color de fondo, por ejemplo '#FFFFFF' o 'none' |
| `--journal {acs,aps,elsevier,generic,iop,nature,rsc,wiley}` | anchos de columna de la editorial (default: `generic`) |
| `--width` | ancho: single / onehalf / double, o un número en mm |
| `--aspect ASPECT` | relación alto/ancho de la figura |
| `--mono` | monocromo: tinta negra y patrones de línea (para revistas que cobran el color) |
| `--dashes {auto,always,never}` | patrones de línea como codificación secundaria (default: `auto`) |
| `--title` | título dentro de la figura (por defecto ninguno: en un artículo el texto va en el pie de figura) |
| `--gap-label` | anotar el valor del gap dentro de la gráfica |
| `--panel` | etiqueta de panel, por ejemplo '(a)' |

**Fundamento físico:** [`olla-dft teoria plot`](TEORIA.md)

### `gap`

solo el reporte de band gap (rápido)

**Uso:** `olla-dft gap [-h] [--prefix PREFIX] [path]`

**Argumentos:**

- `path` — carpeta del cálculo (o ruta al .xml)

**Opciones:**

| Opción | Descripción |
|---|---|
| `--prefix` | prefix del cálculo (se detecta solo) |

**Fundamento físico:** [`olla-dft teoria gap`](TEORIA.md)

### `fermi`

exportar la superficie de Fermi en BXSF

**Uso:** `olla-dft fermi [-h] [-o OUTDIR]`

**Opciones:**

| Opción | Descripción |
|---|---|
| `-o, --outdir` | carpeta de salida (default: `transporte`) |

**Fundamento físico:** [`olla-dft teoria fermi`](TEORIA.md)

### `effmass`

masa efectiva por ajuste parabólico de bandas

**Uso:** `olla-dft effmass [-h] [-o OUTDIR] [--bands-dir BANDS_DIR] [--collect] [--run] [--half-width HALF_WIDTH] [--points POINTS] [--window WINDOW] [--min-points MIN_POINTS] [--pseudo-dir PSEUDO_DIR] [--pseudo EL=UPF] [--ecutwfc ECUTWFC] [--ecutrho ECUTRHO] [--pw-cmd PW_CMD] [--nproc NPROC] [--timeout TIMEOUT] file`

**Argumentos:**

- `file` — estructura de entrada (CIF, POSCAR, input de pw.x...)

**Opciones:**

| Opción | Descripción |
|---|---|
| `-o, --outdir` | carpeta de salida (default: `masa_efectiva`) |
| `--bands-dir` | carpeta con un cálculo de bandas ya hecho (de ahí salen VBM y CBM) |
| `--collect` | leer el cálculo fino ya corrido |
| `--run` | correr el cálculo fino en cuanto se prepare |
| `--half-width HALF_WIDTH` | semiancho de cada línea en Å⁻¹ (default: `0.06`) |
| `--points POINTS` | puntos k por línea (impar) (default: `21`) |
| `--window WINDOW` | semiancho del ajuste rápido sobre el camino, en Å⁻¹ a cada lado del extremo (por omisión, la mitad del límite parabólico: ±0.06) |
| `--min-points MIN_POINTS` | puntos mínimos del ajuste rápido (default: `7`) |
| `--pseudo-dir` | carpeta con los pseudopotenciales UPF (si no se da, la de 'olla-dft config') |
| `--pseudo EL=UPF` | forzar un pseudopotencial concreto, por ejemplo Fe=Fe.rel-pbe.UPF. Se puede repetir. Sin esto, Olla-DFT elige con 'olla-dft pseudos' |
| `--ecutwfc ECUTWFC` | cutoff de funciones de onda en Ry (si no se da, el recomendado por los UPF) |
| `--ecutrho ECUTRHO` | cutoff de densidad en Ry (si no se da, el recomendado por los UPF) |
| `--pw-cmd` | ejecutable de pw.x para --run; de su ruta salen los demás binarios de QE |
| `--nproc NPROC` | número de procesos MPI para los cálculos que se lanzan con --run |
| `--timeout TIMEOUT` | límite de tiempo en segundos para cada ejecución de pw.x |

**Fundamento físico:** [`olla-dft teoria effmass`](TEORIA.md)

### `wannier`

funciones de Wannier: interpolar bandas, centros y dispersión, sin necesitar wannier90

**Uso:** `olla-dft wannier [-h] [-o OUTDIR] [-g NxNxN] [-p SITIO:ORBITAL] [--bands BANDS] [--exclude 5-8] [--window MIN:MAX] [--frozen MIN:MAX] [--no-minimize] [--iterations ITERATIONS] [--points POINTS] [--dft-bands DIR] [--no-dft-bands] [--dos N] [--sigma SIGMA] [--run] [--collect] [--pw-cmd PW_CMD] [--pw2wan-cmd PW2WAN_CMD] [--nproc NPROC] [--timeout TIMEOUT] [--pseudo-dir PSEUDO_DIR] [--ecutwfc ECUTWFC] [--ecutrho ECUTRHO] [--kgrid NxNxN] [--insulator] [--dpi DPI] [--format FORMAT] [--no-plot] [-t TEMPLATE] [--font {sans,serif,latex}] [--usetex] [--palette PALETTE] [--background BACKGROUND] [--journal {acs,aps,elsevier,generic,iop,nature,rsc,wiley}] [--width WIDTH] [--mono] [file]`

**Argumentos:**

- `file` — estructura

**Opciones:**

| Opción | Descripción |
|---|---|
| `-o, --outdir` | carpeta de salida (default: `wannier`) |
| `-g, --grid NxNxN` | malla COMPLETA de puntos k (default 4x4x4). Es la que fija la calidad de la interpolación |
| `-p, --projections SITIO:ORBITAL` | orbitales de prueba: 'Si:sp3', 'O:p;Ti:d', 'f=0.125,0.125,0.125:s'. Varias separadas por ';'. Con 'auto' se ponen s y p en cada átomo (default: `auto`) |
| `--bands BANDS` | bandas del nscf (default: las que hagan falta) |
| `--exclude 5-8` | bandas que NO entran en la wannierización |
| `--window MIN:MAX` | ventana exterior de desenredado en eV: de qué bandas se puede elegir el subespacio. Hace falta cuando las bandas están enredadas con otras (conducción, metales) |
| `--frozen MIN:MAX` | ventana congelada en eV: las bandas de dentro se reproducen EXACTAS. Suele ser la valencia más el trozo de conducción que te importe |
| `--no-minimize` | quedarse en la gauge de proyección, sin minimizar la dispersión |
| `--iterations ITERATIONS` | pasos de minimización (default 500) |
| `--points POINTS` | puntos por tramo del camino interpolado (default: `30`) |
| `--dft-bands DIR` | carpeta con el cálculo de bandas de DFT con el que comparar; sin esto no hay validación de verdad |
| `--no-dft-bands` | con --run, saltarse el paso 4 de bandas |
| `--dos N` | además, DOS interpolada en una malla NxNxN |
| `--sigma SIGMA` | ensanchamiento de la DOS interpolada (eV) (default: `0.05`) |
| `--run` | lanzar los cuatro pasos en orden |
| `--collect` | analizar lo que ya está corrido |
| `--pw-cmd` | ejecutable de pw.x para --run; de su ruta salen los demás binarios de QE |
| `--pw2wan-cmd` | ejecutable de pw2wannier90.x (default: al lado de pw.x) |
| `--nproc NPROC` | número de procesos MPI para los cálculos que se lanzan con --run |
| `--timeout TIMEOUT` | límite de tiempo en segundos para cada ejecución de pw.x |
| `--pseudo-dir` | carpeta con los pseudopotenciales UPF (si no se da, la de 'olla-dft config') |
| `--ecutwfc ECUTWFC` | cutoff de funciones de onda en Ry (si no se da, el recomendado por los UPF) |
| `--ecutrho ECUTRHO` | cutoff de densidad en Ry (si no se da, el recomendado por los UPF) |
| `--kgrid NxNxN` | malla del scf inicial |
| `--insulator` | occupations='fixed' (aislantes; default: smearing) |
| `--dpi DPI` | resolución de las figuras de mapa de bits en puntos por pulgada (default: 600) |
| `--format` | formatos de figura separados por coma: pdf,png,svg,eps,tif (default: pdf,png) |
| `--no-plot` | solo exportar los datos, sin generar la figura |
| `-t, --template` | plantilla visual de la figura (la lista sale con 'olla-dft templates list') |
| `--font {sans,serif,latex}` | familia tipográfica: sans, serif o latex (Computer Modern) |
| `--usetex` | componer los textos de la figura con LaTeX real (necesita una instalación de LaTeX) |
| `--palette` | paleta de colores: grayscale, okabe-ito, okabe-ito-dark o colores hex separados por coma |
| `--background` | color de fondo de la figura, por ejemplo '#FFFFFF' o 'none' |
| `--journal {acs,aps,elsevier,generic,iop,nature,rsc,wiley}` | ancho de figura según la revista (default: generic) |
| `--width` | ancho de la figura: single, onehalf o double, o un número en milímetros |
| `--mono` | versión en escala de grises: tinta negra y patrones de línea |

**Fundamento físico:** [`olla-dft teoria wannier`](TEORIA.md)

### `unfold`

desdoblar las bandas de una supercelda sobre la zona de Brillouin primitiva

**Uso:** `olla-dft unfold [-h] [-o OUTDIR] [--prefix PREFIX] [--bands BANDS] [--spin {up,dw}] [--emin EMIN] [--emax EMAX] [--dpi DPI] [--format FORMAT] [--no-plot] [-t TEMPLATE] [--font {sans,serif,latex}] [--usetex] [--palette PALETTE] [--background BACKGROUND] [--journal {acs,aps,elsevier,generic,iop,nature,rsc,wiley}] [--width WIDTH] [--mono] path primitive`

**Argumentos:**

- `path` — carpeta del calculo de bandas de la supercelda
- `primitive` — estructura de la celda PRIMITIVA

**Opciones:**

| Opción | Descripción |
|---|---|
| `-o, --outdir` | carpeta de salida (default: `.`) |
| `--prefix` | prefix del cálculo (se detecta solo) |
| `--bands BANDS` | cuantas bandas desdoblar (desde la mas baja) |
| `--spin {up,dw}` | canal de espin a desdoblar si el calculo es lsda (se desdobla UN canal por corrida; por omision, up) (default: `up`) |
| `--emin EMIN` | límite inferior del eje de energía (eV) (default: `-6.0`) |
| `--emax EMAX` | límite superior del eje de energía (eV) (default: `6.0`) |
| `--dpi DPI` | resolución de las figuras de mapa de bits en puntos por pulgada (default: 600) |
| `--format` | formatos de figura separados por coma: pdf,png,svg,eps,tif (default: pdf,png) |
| `--no-plot` | solo exportar los datos, sin generar la figura |
| `-t, --template` | plantilla visual de la figura (la lista sale con 'olla-dft templates list') |
| `--font {sans,serif,latex}` | familia tipográfica: sans, serif o latex (Computer Modern) |
| `--usetex` | componer los textos de la figura con LaTeX real (necesita una instalación de LaTeX) |
| `--palette` | paleta de colores: grayscale, okabe-ito, okabe-ito-dark o colores hex separados por coma |
| `--background` | color de fondo de la figura, por ejemplo '#FFFFFF' o 'none' |
| `--journal {acs,aps,elsevier,generic,iop,nature,rsc,wiley}` | ancho de figura según la revista (default: generic) |
| `--width` | ancho de la figura: single, onehalf o double, o un número en milímetros |
| `--mono` | versión en escala de grises: tinta negra y patrones de línea |

**Fundamento físico:** [`olla-dft teoria unfold`](TEORIA.md)

### `topology`

Chern y lazos de Wilson de un modelo Wannier

**Uso:** `olla-dft topology [-h] (--occupied N | --fermi EV) [-g NxN] [--plane {xy,xz,yz}] [--fixed K] [--gap-tol EV] [-o OUTDIR] [--dpi DPI] [--format FORMAT] [--no-plot] [-t TEMPLATE] [--font {sans,serif,latex}] [--usetex] [--palette PALETTE] [--background BACKGROUND] [--journal {acs,aps,elsevier,generic,iop,nature,rsc,wiley}] [--width WIDTH] [--mono] MODELO`

**Argumentos:**

- `MODELO` — archivo *_hr.dat o carpeta que contenga WANNIER_hr.dat

**Opciones:**

| Opción | Descripción |
|---|---|
| `--occupied N` | número de bandas ocupadas del subespacio aislado |
| `--fermi EV` | nivel de Fermi; se rechaza si corta una banda |
| `-g, --grid NxN` | malla periódica de la sección 2D (default: 40x40) |
| `--plane {xy,xz,yz}` | plano orientado de la sección del BZ (default: xy) |
| `--fixed K` | coordenada fraccionaria perpendicular (default: 0) |
| `--gap-tol EV` | gap directo mínimo para aceptar el invariante (default: 1e-8) |
| `-o, --outdir` | carpeta de salida (default: `topology`) |
| `--dpi DPI` | resolución de las figuras de mapa de bits en puntos por pulgada (default: 600) |
| `--format` | formatos de figura separados por coma: pdf,png,svg,eps,tif (default: pdf,png) |
| `--no-plot` | solo exportar los datos, sin generar la figura |
| `-t, --template` | plantilla visual de la figura (la lista sale con 'olla-dft templates list') |
| `--font {sans,serif,latex}` | familia tipográfica: sans, serif o latex (Computer Modern) |
| `--usetex` | componer los textos de la figura con LaTeX real (necesita una instalación de LaTeX) |
| `--palette` | paleta de colores: grayscale, okabe-ito, okabe-ito-dark o colores hex separados por coma |
| `--background` | color de fondo de la figura, por ejemplo '#FFFFFF' o 'none' |
| `--journal {acs,aps,elsevier,generic,iop,nature,rsc,wiley}` | ancho de figura según la revista (default: generic) |
| `--width` | ancho de la figura: single, onehalf o double, o un número en milímetros |
| `--mono` | versión en escala de grises: tinta negra y patrones de línea |

**Fundamento físico:** [`olla-dft teoria topology`](TEORIA.md)

### `hubbard`

U de Hubbard por respuesta lineal (hp.x), en vez de copiarlo de un articulo

**Uso:** `olla-dft hubbard [-h] [-o OUTDIR] [--species SPECIES] [--qgrid QGRID] [--projection {atomic,ortho-atomic,norm-atomic,wannier,pseudo}] [--hubbard-style {legacy,card}] [--cycle] [--max-iter MAX_ITER] [--tol TOL] [--mixing MIXING] [--collect] [--pw-cmd PW_CMD] [--nproc NPROC] [--pseudo-dir PSEUDO_DIR] [--pseudo EL=UPF] [--ecutwfc ECUTWFC] [--ecutrho ECUTRHO] [--kspacing KSPACING] [--metal] [--nspin {1,2}] [--mag MAG] [--intersite] [--v-threshold eV] file`

**Argumentos:**

- `file` — estructura

**Opciones:**

| Opción | Descripción |
|---|---|
| `-o, --outdir` | carpeta de salida (default: `hubbard`) |
| `--species` | especies a perturbar, separadas por coma. Por omision, los metales de transicion y tierras raras de la estructura |
| `--qgrid` | malla de q de la respuesta lineal; equivale a una supercelda de nq1*nq2*nq3 celdas (default: `2x2x2`) |
| `--projection {atomic,ortho-atomic,norm-atomic,wannier,pseudo}` | esquema de proyeccion. El U SOLO vale con el mismo esquema con el que se calculo (default: `ortho-atomic`) |
| `--hubbard-style {legacy,card}` | sintaxis de DFT+U del scf: legacy = lda_plus_u (QE <= 7.0), card = tarjeta HUBBARD (QE >= 7.1, donde la sintaxis vieja es un error) (default: `legacy`) |
| `--cycle` | ciclo de autoconsistencia completo: scf -> hp.x -> scf con el U nuevo, hasta que deje de moverse |
| `--max-iter MAX_ITER` | iteraciones máximas del ciclo scf -> hp.x -> scf con --cycle (default: 6) |
| `--tol TOL` | cambio en eV por debajo del cual se da por convergido (default: `0.05`) |
| `--mixing MIXING` | amortiguacion del paso; bajalo a 0.5 si oscila (default: `1.0`) |
| `--collect` | leer los resultados de un cálculo ya corrido en lugar de preparar los inputs |
| `--pw-cmd` | ejecutable de pw.x para --run; de su ruta salen los demás binarios de QE |
| `--nproc NPROC` | número de procesos MPI para los cálculos que se lanzan con --run |
| `--pseudo-dir` | carpeta con los pseudopotenciales UPF (si no se da, la de 'olla-dft config') |
| `--pseudo EL=UPF` | forzar un pseudopotencial concreto, por ejemplo Fe=Fe.rel-pbe.UPF. Se puede repetir. Sin esto, Olla-DFT elige con 'olla-dft pseudos' |
| `--ecutwfc ECUTWFC` | cutoff de funciones de onda en Ry (si no se da, el recomendado por los UPF) |
| `--ecutrho ECUTRHO` | cutoff de densidad en Ry (si no se da, el recomendado por los UPF) |
| `--kspacing KSPACING` | espaciado de la malla k en Å^-1 |
| `--metal` | sistema metálico: ocupaciones con smearing en vez de fijas |
| `--nspin {1,2}` | 2 activa la polarización de espín (default: 1) |
| `--mag` | magnetización inicial: un número (0.5) o por elemento (Fe=0.7,O=0). Implica --nspin 2 |
| `--intersite` | además de las U, leer los V intersitio que hp.x ya escribe y generar la tarjeta HUBBARD de QE >= 7.1 |
| `--v-threshold eV` | V por debajo de esto no se lista ni se escribe (default: `0.01`) |

**Fundamento físico:** [`olla-dft teoria hubbard`](TEORIA.md)

## Espectros y respuesta

### `optics`

ε(ω), absorción y Tauc con epsilon.x (pseudos NC)

**Uso:** `olla-dft optics [-h] [-o OUTDIR] [--run] [--collect] [--pw-cmd PW_CMD] [--nproc NPROC] [-j N] [--redo] [--max-time T] [--estimate] [--timeout TIMEOUT] [--pseudo-dir PSEUDO_DIR] [--pseudo EL=UPF] [--ecutwfc ECUTWFC] [--ecutrho ECUTRHO] [--kspacing KSPACING] [--insulator] [--dpi DPI] [--format FORMAT] [--no-plot] [-t TEMPLATE] [--size {paper,poster,presentation}] [--font {sans,serif,latex}] [--usetex] [--palette PALETTE] [--background BACKGROUND] [--journal {acs,aps,elsevier,generic,iop,nature,rsc,wiley}] [--width WIDTH] [--aspect ASPECT] [--mono] [--wmax WMAX] [--smear SMEAR] [--metal] [--suite] [--tauc {direct,indirect}] [--scissor SCISSOR] file`

**Argumentos:**

- `file` — estructura (CIF, POSCAR, input de pw.x...)

**Opciones:**

| Opción | Descripción |
|---|---|
| `-o, --outdir` | carpeta del barrido (default: `opticas`) |
| `--run` | ejecutar los cálculos ahora, uno tras otro |
| `--collect` | solo analizar cálculos ya corridos |
| `--pw-cmd` | ejecutable de pw.x (anula la configuración) |
| `--nproc NPROC` | procesos MPI por cálculo |
| `-j, --jobs N` | cálculos simultáneos (default: 1). Sin --nproc, los hilos de la máquina se reparten entre ellos |
| `--redo` | rehacer también los cálculos que ya estaban terminados |
| `--max-time T` | presupuesto TOTAL de tiempo: 90m, 2h, 3600. Al agotarse no se lanzan más y el barrido queda reanudable |
| `--estimate` | estimar cuánto va a tardar el barrido y salir, usando el histórico de 'olla-dft db' |
| `--timeout TIMEOUT` | límite en segundos por cálculo |
| `--pseudo-dir` | carpeta de pseudopotenciales |
| `--pseudo EL=UPF` | forzar un pseudopotencial concreto, por ejemplo Fe=Fe.rel-pbe.UPF. Se puede repetir. Sin esto, Olla-DFT elige con 'olla-dft pseudos' |
| `--ecutwfc ECUTWFC` | cutoff de ondas (Ry) |
| `--ecutrho ECUTRHO` | cutoff de densidad (Ry) |
| `--kspacing KSPACING` | espaciado k en Å^-1 |
| `--insulator` | occupations='fixed' |
| `--dpi DPI` | resolución de las figuras de mapa de bits en puntos por pulgada (default: 600) |
| `--format` | formatos de figura separados por coma: pdf,png,svg,eps,tif (default: pdf,png) |
| `--no-plot` | solo exportar los datos, sin generar la figura |
| `-t, --template` | plantilla visual de la figura |
| `--size {paper,poster,presentation}` | tamaño de figura: paper, presentation o poster |
| `--font {sans,serif,latex}` | familia tipográfica: sans, serif o latex (Computer Modern) |
| `--usetex` | componer los textos de la figura con LaTeX real (necesita una instalación de LaTeX) |
| `--palette` | paleta de colores: grayscale, okabe-ito, okabe-ito-dark o colores hex separados por coma |
| `--background` | color de fondo de la figura, por ejemplo '#FFFFFF' o 'none' |
| `--journal {acs,aps,elsevier,generic,iop,nature,rsc,wiley}` | ancho de figura según la revista (default: generic) |
| `--width` | ancho de la figura: single, onehalf o double, o un número en milímetros |
| `--aspect ASPECT` | relación alto/ancho de la figura |
| `--mono` | versión en escala de grises: tinta negra y patrones de línea |
| `--wmax WMAX` | energía máxima del espectro (eV) (default: `20.0`) |
| `--smear SMEAR` | ensanchamiento interbanda (eV) (default: `0.1`) |
| `--metal` | sistema metálico (ocupaciones con smearing) |
| `--suite` | además, exportar JSON de intercambio para las otras apps de la suite |
| `--tauc {direct,indirect}` | tipo de transición para la gráfica de Tauc (default: `direct`) |
| `--scissor SCISSOR` | corrimiento rígido del gap en eV (gap experimental o GW menos el gap del cálculo); desplaza ε2 y rehace ε1 por Kramers-Kronig |

**Fundamento físico:** [`olla-dft teoria optics`](TEORIA.md)

### `tddft`

absorcion optica con TDDFPT: deja que el electron excitado y su hueco se vean

**Uso:** `olla-dft tddft [-h] [-o OUTDIR] [--method {lanczos,davidson}] [--iter ITER] [--pol {1,2,3,4}] [--states STATES] [--emin EMIN] [--emax EMAX] [--broadening BROADENING] [--scissor SCISSOR] [--extrapolation {no,constant,osc}] [--tamm-dancoff] [--rpa] [--gamma] [--gap GAP] [--compare OPTICS.dat] [--nbnd NBND] [--collect] [--pseudo-dir PSEUDO_DIR] [--pseudo EL=UPF] [--ecutwfc ECUTWFC] [--ecutrho ECUTRHO] [--kspacing KSPACING] [--metal] [--dpi DPI] [--format FORMAT] [--no-plot] [-t TEMPLATE] [--font {sans,serif,latex}] [--usetex] [--palette PALETTE] [--background BACKGROUND] [--journal {acs,aps,elsevier,generic,iop,nature,rsc,wiley}] [--width WIDTH] [--mono] [file]`

**Argumentos:**

- `file` — estructura

**Opciones:**

| Opción | Descripción |
|---|---|
| `-o, --outdir` | carpeta de salida (default: `tddft`) |
| `--method {lanczos,davidson}` | lanczos da el espectro entero; davidson da las primeras excitaciones una a una (default: `lanczos`) |
| `--iter ITER` | iteraciones de Lanczos: manda la resolucion (default: `500`) |
| `--pol {1,2,3,4}` | 1/2/3 = xx/yy/zz, 4 = tensor completo (default: `4`) |
| `--states STATES` | excitaciones a buscar (davidson) (default: `10`) |
| `--emin EMIN` | límite inferior del eje de energía (eV) |
| `--emax EMAX` | límite superior del eje de energía (eV) (default: `15.0`) |
| `--broadening BROADENING` | ensanchamiento en eV (default 0.05). Con --collect fija el umbral de deteccion del exciton; si se omite se lee de spectrum.in |
| `--scissor SCISSOR` | corrimiento rigido de las bandas vacias en eV (solo lanczos): compensa el gap subestimado |
| `--extrapolation {no,constant,osc}` | extrapolación de la cadena de Lanczos en el espectro: no, constant u osc (default: osc) |
| `--tamm-dancoff` | aproximacion de Tamm-Dancoff: mas barata, no exacta |
| `--rpa` | apagar el kernel xc, para ver cuanto aporta |
| `--gamma` | forzar K_POINTS gamma. Se detecta solo cuando la estructura es una molecula |
| `--gap GAP` | gap de particulas independientes en eV, para detectar si hay exciton ligado |
| `--compare OPTICS.dat` | superponer el espectro de 'olla-dft optics' |
| `--nbnd NBND` | número de bandas del scf (por omisión, las que decida pw.x) |
| `--collect` | leer los resultados de un cálculo ya corrido en lugar de preparar los inputs |
| `--pseudo-dir` | carpeta con los pseudopotenciales UPF (si no se da, la de 'olla-dft config') |
| `--pseudo EL=UPF` | forzar un pseudopotencial concreto, por ejemplo Fe=Fe.rel-pbe.UPF; se puede repetir |
| `--ecutwfc ECUTWFC` | cutoff de funciones de onda en Ry (si no se da, el recomendado por los UPF) |
| `--ecutrho ECUTRHO` | cutoff de densidad en Ry (si no se da, el recomendado por los UPF) |
| `--kspacing KSPACING` | espaciado de la malla k en Å^-1 |
| `--metal` | sistema metálico: ocupaciones con smearing en vez de fijas |
| `--dpi DPI` | resolución de las figuras de mapa de bits en puntos por pulgada (default: 600) |
| `--format` | formatos de figura separados por coma: pdf,png,svg,eps,tif (default: pdf,png) |
| `--no-plot` | solo exportar los datos, sin generar la figura |
| `-t, --template` | plantilla visual de la figura (la lista sale con 'olla-dft templates list') |
| `--font {sans,serif,latex}` | familia tipográfica: sans, serif o latex (Computer Modern) |
| `--usetex` | componer los textos de la figura con LaTeX real (necesita una instalación de LaTeX) |
| `--palette` | paleta de colores: grayscale, okabe-ito, okabe-ito-dark o colores hex separados por coma |
| `--background` | color de fondo de la figura, por ejemplo '#FFFFFF' o 'none' |
| `--journal {acs,aps,elsevier,generic,iop,nature,rsc,wiley}` | ancho de figura según la revista (default: generic) |
| `--width` | ancho de la figura: single, onehalf o double, o un número en milímetros |
| `--mono` | versión en escala de grises: tinta negra y patrones de línea |

**Fundamento físico:** [`olla-dft teoria tddft`](TEORIA.md)

### `xanes`

XANES/NEXAFS: absorcion de rayos X cerca del borde (xspectra.x)

**Uso:** `olla-dft xanes [-h] [-o OUTDIR] [--element ELEMENT] [--site SITE] [--edge EDGE] [--core-hole UPF] [--polarization POLARIZATION] [--average] [--emin EMIN] [--emax EMAX] [--broadening BROADENING] [--r-paw R_PAW] [--collect] [--pseudo-dir PSEUDO_DIR] [--pseudo EL=UPF] [--ecutwfc ECUTWFC] [--ecutrho ECUTRHO] [--kspacing KSPACING] [--metal] [--dpi DPI] [--format FORMAT] [--no-plot] [-t TEMPLATE] [--font {sans,serif,latex}] [--usetex] [--palette PALETTE] [--background BACKGROUND] [--journal {acs,aps,elsevier,generic,iop,nature,rsc,wiley}] [--width WIDTH] [--mono] file`

**Argumentos:**

- `file` — estructura

**Opciones:**

| Opción | Descripción |
|---|---|
| `-o, --outdir` | carpeta de salida (default: `xanes`) |
| `--element` | elemento que absorbe |
| `--site SITE` | cual atomo de ese elemento (desde 0) |
| `--edge` | borde: K, L1, L2, L3 o L23 (los que calcula xspectra.x; los bordes M no) (default: `K`) |
| `--core-hole UPF` | pseudopotencial con hueco de core (olla-dft corehole) |
| `--polarization` | direccion del campo electrico, por ejemplo '0 0 1' (default: `1 0 0`) |
| `--average` | tres direcciones ortogonales y promedio: es lo que corresponde a una muestra en polvo |
| `--emin EMIN` | límite inferior del eje de energía (eV) (default: `-10.0`) |
| `--emax EMAX` | límite superior del eje de energía (eV) (default: `30.0`) |
| `--broadening BROADENING` | ensanchamiento en eV (xgamma) (default: `0.8`) |
| `--r-paw R_PAW` | radio de la esfera PAW del absorbedor para xspectra.x, en bohr (default: 3.0) |
| `--collect` | leer los resultados de un cálculo ya corrido en lugar de preparar los inputs |
| `--pseudo-dir` | carpeta con los pseudopotenciales UPF (si no se da, la de 'olla-dft config') |
| `--pseudo EL=UPF` | forzar un pseudopotencial concreto, por ejemplo Fe=Fe.rel-pbe.UPF. Se puede repetir. Sin esto, Olla-DFT elige con 'olla-dft pseudos' |
| `--ecutwfc ECUTWFC` | cutoff de funciones de onda en Ry (si no se da, el recomendado por los UPF) |
| `--ecutrho ECUTRHO` | cutoff de densidad en Ry (si no se da, el recomendado por los UPF) |
| `--kspacing KSPACING` | espaciado de la malla k en Å^-1 |
| `--metal` | sistema metálico: ocupaciones con smearing en vez de fijas |
| `--dpi DPI` | resolución de las figuras de mapa de bits en puntos por pulgada (default: 600) |
| `--format` | formatos de figura separados por coma: pdf,png,svg,eps,tif (default: pdf,png) |
| `--no-plot` | solo exportar los datos, sin generar la figura |
| `-t, --template` | plantilla visual de la figura (la lista sale con 'olla-dft templates list') |
| `--font {sans,serif,latex}` | familia tipográfica: sans, serif o latex (Computer Modern) |
| `--usetex` | componer los textos de la figura con LaTeX real (necesita una instalación de LaTeX) |
| `--palette` | paleta de colores: grayscale, okabe-ito, okabe-ito-dark o colores hex separados por coma |
| `--background` | color de fondo de la figura, por ejemplo '#FFFFFF' o 'none' |
| `--journal {acs,aps,elsevier,generic,iop,nature,rsc,wiley}` | ancho de figura según la revista (default: generic) |
| `--width` | ancho de la figura: single, onehalf o double, o un número en milímetros |
| `--mono` | versión en escala de grises: tinta negra y patrones de línea |

**Fundamento físico:** [`olla-dft teoria xanes`](TEORIA.md)

### `xps`

corrimientos de nivel de core (estado inicial)

**Uso:** `olla-dft xps [-h] [-o OUTDIR] [--core-hole EL=UPF] [--collect] [--suite] [--pseudo-dir PSEUDO_DIR] [--pseudo EL=UPF] [--ecutwfc ECUTWFC] [--ecutrho ECUTRHO] [--kspacing KSPACING] [--metal] file`

**Argumentos:**

- `file` — estructura de entrada (CIF, POSCAR, input de pw.x...)

**Opciones:**

| Opción | Descripción |
|---|---|
| `-o, --outdir` | carpeta de salida (default: `xps`) |
| `--core-hole EL=UPF` | pseudopotencial con hueco de core, por ejemplo Si=Si.star1s.UPF. Se puede repetir. Sin esto, initial_state.x devuelve una tabla de ceros |
| `--collect` | leer los resultados de un cálculo ya corrido en lugar de preparar los inputs |
| `--suite` | además, exportar JSON de intercambio para las otras apps de la suite |
| `--pseudo-dir` | carpeta con los pseudopotenciales UPF (si no se da, la de 'olla-dft config') |
| `--pseudo EL=UPF` | forzar un pseudopotencial concreto, por ejemplo Fe=Fe.rel-pbe.UPF. Se puede repetir. Sin esto, Olla-DFT elige con 'olla-dft pseudos' |
| `--ecutwfc ECUTWFC` | cutoff de funciones de onda en Ry (si no se da, el recomendado por los UPF) |
| `--ecutrho ECUTRHO` | cutoff de densidad en Ry (si no se da, el recomendado por los UPF) |
| `--kspacing KSPACING` | espaciado de la malla k en Å^-1 |
| `--metal` | sistema metálico: ocupaciones con smearing en vez de fijas |

**Fundamento físico:** [`olla-dft teoria xps`](TEORIA.md)

### `corehole`

generar el par de pseudopotenciales normal + hueco de core (ld1.x), para XPS y XANES

**Uso:** `olla-dft corehole [-h] [-o OUTDIR] [--edge EDGE] [--functional FUNCTIONAL] [--rcut RCUT] [--rel {0,1,2}] [--semicore] [--pseudotype {1,2,3}] [--plain] [--only-inputs] [--projectors {1,2}] [--ld1-cmd LD1_CMD] [--core-wfc UPF] [--orbital ORBITAL] [--output OUTPUT] [element]`

**Argumentos:**

- `element` — simbolo del elemento, por ejemplo Si

**Opciones:**

| Opción | Descripción |
|---|---|
| `-o, --outdir` | carpeta de salida (default: `pseudos`) |
| `--edge` | borde/nivel del hueco: K (1s), L1 (2s), L23 (2p), M1, M23, M45 (default: `K`) |
| `--functional` | funcional del pseudopotencial; tiene que ser el mismo con el que vas a correr pw.x (default: `PBE`) |
| `--rcut RCUT` | radio de corte en bohr (por omision, uno por fila de la tabla periodica) |
| `--rel {0,1,2}` | 0 no relativista, 1 escalar, 2 completo |
| `--semicore` | mete la capa (n-1)s(n-1)p en la valencia |
| `--pseudotype {1,2,3}` | 1 y 2 conservan la norma, 3 es ultrasuave (default: `2`) |
| `--plain` | generar SOLO el pseudopotencial normal, sin el de hueco de core. Sirve para tener un pseudo consistente de un elemento que no lo admita |
| `--only-inputs` | escribir los inputs de ld1.x sin ejecutarlos |
| `--projectors {1,2}` | proyectores GIPAW por canal. XSpectra recomienda 2, pero con 2 el pseudo sale ultrasuave y suele haber que ajustar --rcut a mano (default: `1`) |
| `--ld1-cmd` | ruta a ld1.x |
| `--core-wfc UPF` | en vez de generar: extraer de un UPF la funcion de onda de core en el formato que lee xspectra.x |
| `--orbital` | orbital a verificar, por ejemplo 1S |
| `--output` | archivo de salida para --core-wfc |

**Fundamento físico:** [`olla-dft teoria corehole`](TEORIA.md)

### `charge`

densidad de carga / ELF / espín con pp.x

**Uso:** `olla-dft charge [-h] [-o OUTDIR] [--field {density,elf,spin,potential,vtotal}] [--axis AXIS] [--rerun] [--pw-cmd PW_CMD] [--nproc NPROC] [--dpi DPI] [--format FORMAT] [--no-plot] [-t TEMPLATE] [--font {sans,serif,latex}] [--usetex] [--palette PALETTE] [--background BACKGROUND] [--journal {acs,aps,elsevier,generic,iop,nature,rsc,wiley}] [--width WIDTH] [--mono] [path]`

**Argumentos:**

- `path` — carpeta del cálculo

**Opciones:**

| Opción | Descripción |
|---|---|
| `-o, --outdir` | carpeta de salida (default: `.`) |
| `--field {density,elf,spin,potential,vtotal}` | campo a calcular con pp.x: density, elf, spin, potential o vtotal (default: density) |
| `--axis` | eje del perfil planar (default: `c`) |
| `--rerun` | volver a correr pp.x aunque ya exista el archivo cube |
| `--pw-cmd` | ejecutable de pw.x para --run; de su ruta salen los demás binarios de QE |
| `--nproc NPROC` | número de procesos MPI para los cálculos que se lanzan con --run |
| `--dpi DPI` | resolución de las figuras de mapa de bits en puntos por pulgada (default: 600) |
| `--format` | formatos de figura separados por coma: pdf,png,svg,eps,tif (default: pdf,png) |
| `--no-plot` | solo exportar los datos, sin generar la figura |
| `-t, --template` | plantilla visual de la figura (la lista sale con 'olla-dft templates list') |
| `--font {sans,serif,latex}` | familia tipográfica: sans, serif o latex (Computer Modern) |
| `--usetex` | componer los textos de la figura con LaTeX real (necesita una instalación de LaTeX) |
| `--palette` | paleta de colores: grayscale, okabe-ito, okabe-ito-dark o colores hex separados por coma |
| `--background` | color de fondo de la figura, por ejemplo '#FFFFFF' o 'none' |
| `--journal {acs,aps,elsevier,generic,iop,nature,rsc,wiley}` | ancho de figura según la revista (default: generic) |
| `--width` | ancho de la figura: single, onehalf o double, o un número en milímetros |
| `--mono` | versión en escala de grises: tinta negra y patrones de línea |

**Fundamento físico:** [`olla-dft teoria charge`](TEORIA.md)

### `charges`

cargas de Löwdin/Bader y diferencia de densidad

**Uso:** `olla-dft charges [-h] [--lowdin LOWDIN] [--bader BADER] [--difference CUBE [CUBE ...]] [--pseudo-dir PSEUDO_DIR] [--axis {0,1,2}] [-o OUTDIR] [--dpi DPI] [--format FORMAT] [--no-plot] [-t TEMPLATE] [--font {sans,serif,latex}] [--usetex] [--palette PALETTE] [--background BACKGROUND] [--journal {acs,aps,elsevier,generic,iop,nature,rsc,wiley}] [--width WIDTH] [--mono] [file]`

**Argumentos:**

- `file` — estructura (para Bader)

**Opciones:**

| Opción | Descripción |
|---|---|
| `--lowdin` | salida de projwfc.x |
| `--bader` | cube de densidad (plot_num=0) |
| `--difference CUBE` | total.cube parte1.cube parte2.cube ... |
| `--pseudo-dir` | carpeta con los UPF del cálculo: de ahí sale Z_valencia para la columna 'neta' (anula config) |
| `--axis {0,1,2}` | eje del perfil planar de la diferencia de densidad: 0, 1 o 2 (default: 2) |
| `-o, --outdir` | carpeta de salida (default: `.`) |
| `--dpi DPI` | resolución de las figuras de mapa de bits en puntos por pulgada (default: 600) |
| `--format` | formatos de figura separados por coma: pdf,png,svg,eps,tif (default: pdf,png) |
| `--no-plot` | solo exportar los datos, sin generar la figura |
| `-t, --template` | plantilla visual de la figura (la lista sale con 'olla-dft templates list') |
| `--font {sans,serif,latex}` | familia tipográfica: sans, serif o latex (Computer Modern) |
| `--usetex` | componer los textos de la figura con LaTeX real (necesita una instalación de LaTeX) |
| `--palette` | paleta de colores: grayscale, okabe-ito, okabe-ito-dark o colores hex separados por coma |
| `--background` | color de fondo de la figura, por ejemplo '#FFFFFF' o 'none' |
| `--journal {acs,aps,elsevier,generic,iop,nature,rsc,wiley}` | ancho de figura según la revista (default: generic) |
| `--width` | ancho de la figura: single, onehalf o double, o un número en milímetros |
| `--mono` | versión en escala de grises: tinta negra y patrones de línea |

**Fundamento físico:** [`olla-dft teoria charges`](TEORIA.md)

### `wf`

función trabajo desde un cálculo con vacío

**Uso:** `olla-dft wf [-h] [-o OUTDIR] [--axis AXIS] [--rerun] [--pw-cmd PW_CMD] [--nproc NPROC] [--dpi DPI] [--format FORMAT] [--no-plot] [-t TEMPLATE] [--font {sans,serif,latex}] [--usetex] [--palette PALETTE] [--background BACKGROUND] [--journal {acs,aps,elsevier,generic,iop,nature,rsc,wiley}] [--width WIDTH] [--mono] [path]`

**Argumentos:**

- `path` — carpeta del cálculo

**Opciones:**

| Opción | Descripción |
|---|---|
| `-o, --outdir` | carpeta de salida (default: `.`) |
| `--axis` | eje del vacío: a/b/c (default c) |
| `--rerun` | volver a correr pp.x aunque exista el cube |
| `--pw-cmd` | ejecutable de pw.x para --run; de su ruta salen los demás binarios de QE |
| `--nproc NPROC` | número de procesos MPI para los cálculos que se lanzan con --run |
| `--dpi DPI` | resolución de las figuras de mapa de bits en puntos por pulgada (default: 600) |
| `--format` | formatos de figura separados por coma: pdf,png,svg,eps,tif (default: pdf,png) |
| `--no-plot` | solo exportar los datos, sin generar la figura |
| `-t, --template` | plantilla visual de la figura (la lista sale con 'olla-dft templates list') |
| `--font {sans,serif,latex}` | familia tipográfica: sans, serif o latex (Computer Modern) |
| `--usetex` | componer los textos de la figura con LaTeX real (necesita una instalación de LaTeX) |
| `--palette` | paleta de colores: grayscale, okabe-ito, okabe-ito-dark o colores hex separados por coma |
| `--background` | color de fondo de la figura, por ejemplo '#FFFFFF' o 'none' |
| `--journal {acs,aps,elsevier,generic,iop,nature,rsc,wiley}` | ancho de figura según la revista (default: generic) |
| `--width` | ancho de la figura: single, onehalf o double, o un número en milímetros |
| `--mono` | versión en escala de grises: tinta negra y patrones de línea |

**Fundamento físico:** [`olla-dft teoria wf`](TEORIA.md)

### `berry`

polarización por fase de Berry: ΔP espontánea, cargas de Born y ferroelectricidad

**Uso:** `olla-dft berry [-h] [-o OUTDIR] [--gdir {1,2,3}] [--nppstr NPPSTR] [--kperp NxN] [-r ARCHIVO] [--displace ATOMO:dx,dy,dz] [--nlambda NLAMBDA] [--run] [--collect] [--redo] [--pw-cmd PW_CMD] [--nproc NPROC] [--timeout TIMEOUT] [--pseudo-dir PSEUDO_DIR] [--ecutwfc ECUTWFC] [--ecutrho ECUTRHO] [--kgrid NxNxN] [--dpi DPI] [--format FORMAT] [--no-plot] [-t TEMPLATE] [--font {sans,serif,latex}] [--usetex] [--palette PALETTE] [--background BACKGROUND] [--journal {acs,aps,elsevier,generic,iop,nature,rsc,wiley}] [--width WIDTH] [--mono] file`

**Argumentos:**

- `file` — estructura (la polar, si hay camino)

**Opciones:**

| Opción | Descripción |
|---|---|
| `-o, --outdir` | carpeta de salida (default: `berry`) |
| `--gdir {1,2,3}` | dirección: vector de la red recíproca (default 3) |
| `--nppstr NPPSTR` | puntos por cuerda de k (default 9); sube hasta que la fase deje de moverse |
| `--kperp NxN` | malla perpendicular a la cuerda (default 6x6) |
| `-r, --reference ARCHIVO` | estructura de referencia, normalmente la centrosimétrica: se interpola un camino adiabático hasta la polar y ΔP es la polarización espontánea |
| `--displace ATOMO:dx,dy,dz` | camino de desplazamiento de un átomo, en Å; la pendiente de P da la carga efectiva de Born |
| `--nlambda NLAMBDA` | puntos del camino (default 5) |
| `--run` | correr los cálculos en cuanto se preparen los inputs |
| `--collect` | leer los resultados de un cálculo ya corrido en lugar de preparar los inputs |
| `--redo` | rehacer también los cálculos que ya estaban terminados |
| `--pw-cmd` | ejecutable de pw.x para --run; de su ruta salen los demás binarios de QE |
| `--nproc NPROC` | número de procesos MPI para los cálculos que se lanzan con --run |
| `--timeout TIMEOUT` | límite de tiempo en segundos para cada ejecución de pw.x |
| `--pseudo-dir` | carpeta con los pseudopotenciales UPF (si no se da, la de 'olla-dft config') |
| `--ecutwfc ECUTWFC` | cutoff de funciones de onda en Ry (si no se da, el recomendado por los UPF) |
| `--ecutrho ECUTRHO` | cutoff de densidad en Ry (si no se da, el recomendado por los UPF) |
| `--kgrid NxNxN` | malla de k del scf, por ejemplo 6x6x6 (por omisión, según el espaciado k) |
| `--dpi DPI` | resolución de las figuras de mapa de bits en puntos por pulgada (default: 600) |
| `--format` | formatos de figura separados por coma: pdf,png,svg,eps,tif (default: pdf,png) |
| `--no-plot` | solo exportar los datos, sin generar la figura |
| `-t, --template` | plantilla visual de la figura (la lista sale con 'olla-dft templates list') |
| `--font {sans,serif,latex}` | familia tipográfica: sans, serif o latex (Computer Modern) |
| `--usetex` | componer los textos de la figura con LaTeX real (necesita una instalación de LaTeX) |
| `--palette` | paleta de colores: grayscale, okabe-ito, okabe-ito-dark o colores hex separados por coma |
| `--background` | color de fondo de la figura, por ejemplo '#FFFFFF' o 'none' |
| `--journal {acs,aps,elsevier,generic,iop,nature,rsc,wiley}` | ancho de figura según la revista (default: generic) |
| `--width` | ancho de la figura: single, onehalf o double, o un número en milímetros |
| `--mono` | versión en escala de grises: tinta negra y patrones de línea |

**Fundamento físico:** [`olla-dft teoria berry`](TEORIA.md)

## Fonones, transporte y temperatura

### `phonons`

fonones DFPT: dispersión, DOS, termodinámica, IR

**Uso:** `olla-dft phonons [-h] [-o OUTDIR] [--run] [--collect] [--pw-cmd PW_CMD] [--nproc NPROC] [-j N] [--redo] [--max-time T] [--estimate] [--timeout TIMEOUT] [--pseudo-dir PSEUDO_DIR] [--pseudo EL=UPF] [--ecutwfc ECUTWFC] [--ecutrho ECUTRHO] [--kspacing KSPACING] [--insulator] [--dpi DPI] [--format FORMAT] [--no-plot] [-t TEMPLATE] [--size {paper,poster,presentation}] [--font {sans,serif,latex}] [--usetex] [--palette PALETTE] [--background BACKGROUND] [--journal {acs,aps,elsevier,generic,iop,nature,rsc,wiley}] [--width WIDTH] [--aspect ASPECT] [--mono] [--qgrid QGRID] [--gamma] [--raman] [--laser LASER] [--suite] [--tscan T1,T2,...] file`

**Argumentos:**

- `file` — estructura (CIF, POSCAR, input de pw.x...)

**Opciones:**

| Opción | Descripción |
|---|---|
| `-o, --outdir` | carpeta del barrido (default: `fonones`) |
| `--run` | ejecutar los cálculos ahora, uno tras otro |
| `--collect` | solo analizar cálculos ya corridos |
| `--pw-cmd` | ejecutable de pw.x (anula la configuración) |
| `--nproc NPROC` | procesos MPI por cálculo |
| `-j, --jobs N` | cálculos simultáneos (default: 1). Sin --nproc, los hilos de la máquina se reparten entre ellos |
| `--redo` | rehacer también los cálculos que ya estaban terminados |
| `--max-time T` | presupuesto TOTAL de tiempo: 90m, 2h, 3600. Al agotarse no se lanzan más y el barrido queda reanudable |
| `--estimate` | estimar cuánto va a tardar el barrido y salir, usando el histórico de 'olla-dft db' |
| `--timeout TIMEOUT` | límite en segundos por cálculo |
| `--pseudo-dir` | carpeta de pseudopotenciales |
| `--pseudo EL=UPF` | forzar un pseudopotencial concreto, por ejemplo Fe=Fe.rel-pbe.UPF. Se puede repetir. Sin esto, Olla-DFT elige con 'olla-dft pseudos' |
| `--ecutwfc ECUTWFC` | cutoff de ondas (Ry) |
| `--ecutrho ECUTRHO` | cutoff de densidad (Ry) |
| `--kspacing KSPACING` | espaciado k en Å^-1 |
| `--insulator` | occupations='fixed' |
| `--dpi DPI` | resolución de las figuras de mapa de bits en puntos por pulgada (default: 600) |
| `--format` | formatos de figura separados por coma: pdf,png,svg,eps,tif (default: pdf,png) |
| `--no-plot` | solo exportar los datos, sin generar la figura |
| `-t, --template` | plantilla visual de la figura |
| `--size {paper,poster,presentation}` | tamaño de figura: paper, presentation o poster |
| `--font {sans,serif,latex}` | familia tipográfica: sans, serif o latex (Computer Modern) |
| `--usetex` | componer los textos de la figura con LaTeX real (necesita una instalación de LaTeX) |
| `--palette` | paleta de colores: grayscale, okabe-ito, okabe-ito-dark o colores hex separados por coma |
| `--background` | color de fondo de la figura, por ejemplo '#FFFFFF' o 'none' |
| `--journal {acs,aps,elsevier,generic,iop,nature,rsc,wiley}` | ancho de figura según la revista (default: generic) |
| `--width` | ancho de la figura: single, onehalf o double, o un número en milímetros |
| `--aspect ASPECT` | relación alto/ancho de la figura |
| `--mono` | versión en escala de grises: tinta negra y patrones de línea |
| `--qgrid` | malla de q, por ejemplo 2x2x2 |
| `--gamma` | solo Γ con dynmat.x: frecuencias y actividades IR |
| `--raman` | además, tensores e intensidades Raman en Γ (lraman; solo pseudos de norma conservada, y es bastante más caro) |
| `--laser LASER` | longitud de onda del láser en nm para simular el espectro Raman (default: `532.0`) |
| `--suite` | además, exportar JSON de intercambio (solo con --gamma) para las apps de FTIR y Raman |
| `--tscan T1,T2,...` | barrido de temperatura ELECTRÓNICA en K: repite los fonones con smearing fermi-dirac a cada una y mira si un modo imaginario se estabiliza al calentar (ondas de densidad de carga, transiciones estructurales) |

**Fundamento físico:** [`olla-dft teoria phonons`](TEORIA.md)

### `elph`

acoplamiento electron-fonon: lambda, Tc y un tau de verdad para el transporte

**Uso:** `olla-dft elph [-h] [-o OUTDIR] [--qgrid QGRID] [--kgrid KGRID] [--kgrid-nscf KGRID_NSCF] [--nsigma NSIGMA] [--sigma SIGMA] [--degauss DEGAUSS] [--debye DEBYE] [--collect] [--pseudo-dir PSEUDO_DIR] [--pseudo EL=UPF] [--ecutwfc ECUTWFC] [--ecutrho ECUTRHO] [--dpi DPI] [--format FORMAT] [--no-plot] [-t TEMPLATE] [--font {sans,serif,latex}] [--usetex] [--palette PALETTE] [--background BACKGROUND] [--journal {acs,aps,elsevier,generic,iop,nature,rsc,wiley}] [--width WIDTH] [--mono] [file]`

**Argumentos:**

- `file` — estructura

**Opciones:**

| Opción | Descripción |
|---|---|
| `-o, --outdir` | carpeta de salida (default: `elph`) |
| `--qgrid` | malla de q de la DFPT, por ejemplo 2x2x2 (default: 2x2x2) |
| `--kgrid` | malla de k del scf |
| `--kgrid-nscf` | malla de k del nscf denso; por omision, el doble de la del scf redondeada a un multiplo de la de q |
| `--nsigma NSIGMA` | cuántos ensanchamientos barre ph.x para lambda (el_ph_nsigma; default: 10) |
| `--sigma SIGMA` | paso del barrido de ensanchamiento, en Ry (default: `0.005`) |
| `--degauss DEGAUSS` | smearing del scf en Ry (default: 0.02) |
| `--debye DEBYE` | temperatura de Debye en K, para marcar el regimen en el que vale la formula de tau |
| `--collect` | leer los resultados de un cálculo ya corrido en lugar de preparar los inputs |
| `--pseudo-dir` | carpeta con los pseudopotenciales UPF (si no se da, la de 'olla-dft config') |
| `--pseudo EL=UPF` | forzar un pseudopotencial concreto, por ejemplo Fe=Fe.rel-pbe.UPF. Se puede repetir. Sin esto, Olla-DFT elige con 'olla-dft pseudos' |
| `--ecutwfc ECUTWFC` | cutoff de funciones de onda en Ry (si no se da, el recomendado por los UPF) |
| `--ecutrho ECUTRHO` | cutoff de densidad en Ry (si no se da, el recomendado por los UPF) |
| `--dpi DPI` | resolución de las figuras de mapa de bits en puntos por pulgada (default: 600) |
| `--format` | formatos de figura separados por coma: pdf,png,svg,eps,tif (default: pdf,png) |
| `--no-plot` | solo exportar los datos, sin generar la figura |
| `-t, --template` | plantilla visual de la figura (la lista sale con 'olla-dft templates list') |
| `--font {sans,serif,latex}` | familia tipográfica: sans, serif o latex (Computer Modern) |
| `--usetex` | componer los textos de la figura con LaTeX real (necesita una instalación de LaTeX) |
| `--palette` | paleta de colores: grayscale, okabe-ito, okabe-ito-dark o colores hex separados por coma |
| `--background` | color de fondo de la figura, por ejemplo '#FFFFFF' o 'none' |
| `--journal {acs,aps,elsevier,generic,iop,nature,rsc,wiley}` | ancho de figura según la revista (default: generic) |
| `--width` | ancho de la figura: single, onehalf o double, o un número en milímetros |
| `--mono` | versión en escala de grises: tinta negra y patrones de línea |

**Fundamento físico:** [`olla-dft teoria elph`](TEORIA.md)

### `transport`

Seebeck, sigma/tau y factor de potencia (CRTA)

**Uso:** `olla-dft transport [-h] [-o OUTDIR] [--run] [--collect] [--pw-cmd PW_CMD] [--nproc NPROC] [-j N] [--redo] [--max-time T] [--estimate] [--timeout TIMEOUT] [--pseudo-dir PSEUDO_DIR] [--pseudo EL=UPF] [--ecutwfc ECUTWFC] [--ecutrho ECUTRHO] [--kspacing KSPACING] [--insulator] [--dpi DPI] [--format FORMAT] [--no-plot] [-t TEMPLATE] [--size {paper,poster,presentation}] [--font {sans,serif,latex}] [--usetex] [--palette PALETTE] [--background BACKGROUND] [--journal {acs,aps,elsevier,generic,iop,nature,rsc,wiley}] [--width WIDTH] [--aspect ASPECT] [--mono] [--grid GRID] [--temperatures TEMPERATURES] [--mu-span MU_SPAN] [--metal] [--nspin {1,2}] [--mag MAG] [--spin-resolved] file`

**Argumentos:**

- `file` — estructura (CIF, POSCAR, input de pw.x...)

**Opciones:**

| Opción | Descripción |
|---|---|
| `-o, --outdir` | carpeta del barrido (default: `transporte`) |
| `--run` | ejecutar los cálculos ahora, uno tras otro |
| `--collect` | solo analizar cálculos ya corridos |
| `--pw-cmd` | ejecutable de pw.x (anula la configuración) |
| `--nproc NPROC` | procesos MPI por cálculo |
| `-j, --jobs N` | cálculos simultáneos (default: 1). Sin --nproc, los hilos de la máquina se reparten entre ellos |
| `--redo` | rehacer también los cálculos que ya estaban terminados |
| `--max-time T` | presupuesto TOTAL de tiempo: 90m, 2h, 3600. Al agotarse no se lanzan más y el barrido queda reanudable |
| `--estimate` | estimar cuánto va a tardar el barrido y salir, usando el histórico de 'olla-dft db' |
| `--timeout TIMEOUT` | límite en segundos por cálculo |
| `--pseudo-dir` | carpeta de pseudopotenciales |
| `--pseudo EL=UPF` | forzar un pseudopotencial concreto, por ejemplo Fe=Fe.rel-pbe.UPF. Se puede repetir. Sin esto, Olla-DFT elige con 'olla-dft pseudos' |
| `--ecutwfc ECUTWFC` | cutoff de ondas (Ry) |
| `--ecutrho ECUTRHO` | cutoff de densidad (Ry) |
| `--kspacing KSPACING` | espaciado k en Å^-1 |
| `--insulator` | occupations='fixed' |
| `--dpi DPI` | resolución de las figuras de mapa de bits en puntos por pulgada (default: 600) |
| `--format` | formatos de figura separados por coma: pdf,png,svg,eps,tif (default: pdf,png) |
| `--no-plot` | solo exportar los datos, sin generar la figura |
| `-t, --template` | plantilla visual de la figura |
| `--size {paper,poster,presentation}` | tamaño de figura: paper, presentation o poster |
| `--font {sans,serif,latex}` | familia tipográfica: sans, serif o latex (Computer Modern) |
| `--usetex` | componer los textos de la figura con LaTeX real (necesita una instalación de LaTeX) |
| `--palette` | paleta de colores: grayscale, okabe-ito, okabe-ito-dark o colores hex separados por coma |
| `--background` | color de fondo de la figura, por ejemplo '#FFFFFF' o 'none' |
| `--journal {acs,aps,elsevier,generic,iop,nature,rsc,wiley}` | ancho de figura según la revista (default: generic) |
| `--width` | ancho de la figura: single, onehalf o double, o un número en milímetros |
| `--aspect ASPECT` | relación alto/ancho de la figura |
| `--mono` | versión en escala de grises: tinta negra y patrones de línea |
| `--grid` | malla del nscf, por ejemplo 16x16x16 |
| `--temperatures` | temperaturas en K separadas por comas (default: `300`) |
| `--mu-span MU_SPAN` | rango de potencial químico alrededor de E_F (eV) (default: `1.0`) |
| `--metal` | sistema metálico: ocupaciones con smearing en vez de fijas |
| `--nspin {1,2}` | 2 activa la polarización de espín en scf y nscf (necesario para --spin-resolved) (default: `1`) |
| `--mag` | magnetización inicial, por ejemplo Fe=0.7 (implica --nspin 2) |
| `--spin-resolved` | separar los dos canales de espín (modelo de dos corrientes) y dar la polarización de la conductividad y la termopotencia de espín |

**Fundamento físico:** [`olla-dft teoria transport`](TEORIA.md)

### `ballistic`

conductancia balistica de Landauer (pwcond.x), para nanocontactos y moleculas entre electrodos

**Uso:** `olla-dft ballistic [-h] [--scatterer SCATTERER] [-o OUTDIR] [--ikind {0,1}] [--emin EMIN] [--emax EMAX] [--points POINTS] [--nz1 NZ1] [--collect] [--pseudo-dir PSEUDO_DIR] [--pseudo EL=UPF] [--ecutwfc ECUTWFC] [--ecutrho ECUTRHO] [--kspacing KSPACING] [--dpi DPI] [--format FORMAT] [--no-plot] [-t TEMPLATE] [--font {sans,serif,latex}] [--usetex] [--palette PALETTE] [--background BACKGROUND] [--journal {acs,aps,elsevier,generic,iop,nature,rsc,wiley}] [--width WIDTH] [--mono] [file]`

**Argumentos:**

- `file` — electrodo: la celda periodica en z

**Opciones:**

| Opción | Descripción |
|---|---|
| `--scatterer` | region de dispersion (la molecula o el defecto). Sin esto solo salen las bandas complejas |
| `-o, --outdir` | carpeta de salida (default: `balistico`) |
| `--ikind {0,1}` | 0 = solo bandas complejas, 1 = conductancia con el mismo electrodo a los dos lados (default: 1 si hay --scatterer, 0 si no). Electrodos distintos (ikind=2 de pwcond.x) no están soportados |
| `--emin EMIN` | límite inferior del eje de energía (eV) (default: `-3.0`) |
| `--emax EMAX` | límite superior del eje de energía (eV) (default: `3.0`) |
| `--points POINTS` | número de energías del barrido de conductancia (default: 61) |
| `--nz1 NZ1` | subdivisiones en z de cada rebanada de pwcond.x (nz1; default: 3) |
| `--collect` | leer los resultados de un cálculo ya corrido en lugar de preparar los inputs |
| `--pseudo-dir` | carpeta con los pseudopotenciales UPF (si no se da, la de 'olla-dft config') |
| `--pseudo EL=UPF` | forzar un pseudopotencial concreto, por ejemplo Fe=Fe.rel-pbe.UPF; se puede repetir |
| `--ecutwfc ECUTWFC` | cutoff de funciones de onda en Ry (si no se da, el recomendado por los UPF) |
| `--ecutrho ECUTRHO` | cutoff de densidad en Ry (si no se da, el recomendado por los UPF) |
| `--kspacing KSPACING` | espaciado de la malla k en Å^-1 |
| `--dpi DPI` | resolución de las figuras de mapa de bits en puntos por pulgada (default: 600) |
| `--format` | formatos de figura separados por coma: pdf,png,svg,eps,tif (default: pdf,png) |
| `--no-plot` | solo exportar los datos, sin generar la figura |
| `-t, --template` | plantilla visual de la figura (la lista sale con 'olla-dft templates list') |
| `--font {sans,serif,latex}` | familia tipográfica: sans, serif o latex (Computer Modern) |
| `--usetex` | componer los textos de la figura con LaTeX real (necesita una instalación de LaTeX) |
| `--palette` | paleta de colores: grayscale, okabe-ito, okabe-ito-dark o colores hex separados por coma |
| `--background` | color de fondo de la figura, por ejemplo '#FFFFFF' o 'none' |
| `--journal {acs,aps,elsevier,generic,iop,nature,rsc,wiley}` | ancho de figura según la revista (default: generic) |
| `--width` | ancho de la figura: single, onehalf o double, o un número en milímetros |
| `--mono` | versión en escala de grises: tinta negra y patrones de línea |

**Fundamento físico:** [`olla-dft teoria ballistic`](TEORIA.md)

### `kappa`

conductividad térmica de red: fc3, ecuación de Boltzmann de fonones y recorrido libre medio

**Uso:** `olla-dft kappa [-h] [-o OUTDIR] [--dim NxNxN] [--dim-fc2 NxNxN] [--distance DISTANCE] [--mesh MESH] [--temps TEMPS] [--isotopes] [--grain UM] [--model MODEL] [--collect] [--force] [--metal] [--pseudo-dir PSEUDO_DIR] [--ecutwfc ECUTWFC] [--ecutrho ECUTRHO] [--kspacing KSPACING] [--dpi DPI] [--format FORMAT] [--no-plot] [-t TEMPLATE] [--font {sans,serif,latex}] [--usetex] [--palette PALETTE] [--background BACKGROUND] [--journal {acs,aps,elsevier,generic,iop,nature,rsc,wiley}] [--width WIDTH] [--mono] file`

**Argumentos:**

- `file` — estructura (celda primitiva)

**Opciones:**

| Opción | Descripción |
|---|---|
| `-o, --outdir` | carpeta de salida (default: `kappa`) |
| `--dim NxNxN` | supercelda de la fc3 (default 2x2x2). Es lo que más cuesta: el número de configuraciones crece deprisa |
| `--dim-fc2 NxNxN` | supercelda MAYOR solo para la parte armónica, que es barata y necesita más alcance |
| `--distance DISTANCE` | desplazamiento finito en Å (default 0.03) |
| `--mesh MESH` | malla de q para la ecuación de Boltzmann (default 13) |
| `--temps` | temperaturas: 100:800:8 o 300,500,700 (default: `100:800:8`) |
| `--isotopes` | añadir dispersión por isótopos con las abundancias naturales (en Si son ~10 %%) |
| `--grain UM` | tamaño de grano en µm: añade dispersión por fronteras |
| `--model` | calcular las fuerzas con un potencial aprendido (mace, chgnet, m3gnet) en vez de con pw.x: segundos en vez de horas, pero el valor absoluto puede estar lejos |
| `--collect` | leer las fuerzas ya calculadas y resolver |
| `--force` | escribir los inputs aunque sean muchísimos |
| `--metal` | sistema metálico (ocupaciones con smearing en los scf de la fc2/fc3). Sin esto se usa occupations='fixed', lo correcto para aislantes |
| `--pseudo-dir` | carpeta con los pseudopotenciales UPF (si no se da, la de 'olla-dft config') |
| `--ecutwfc ECUTWFC` | cutoff de funciones de onda en Ry (si no se da, el recomendado por los UPF) |
| `--ecutrho ECUTRHO` | cutoff de densidad en Ry (si no se da, el recomendado por los UPF) |
| `--kspacing KSPACING` | espaciado de la malla k en Å^-1 (default: `0.35`) |
| `--dpi DPI` | resolución de las figuras de mapa de bits en puntos por pulgada (default: 600) |
| `--format` | formatos de figura separados por coma: pdf,png,svg,eps,tif (default: pdf,png) |
| `--no-plot` | solo exportar los datos, sin generar la figura |
| `-t, --template` | plantilla visual de la figura (la lista sale con 'olla-dft templates list') |
| `--font {sans,serif,latex}` | familia tipográfica: sans, serif o latex (Computer Modern) |
| `--usetex` | componer los textos de la figura con LaTeX real (necesita una instalación de LaTeX) |
| `--palette` | paleta de colores: grayscale, okabe-ito, okabe-ito-dark o colores hex separados por coma |
| `--background` | color de fondo de la figura, por ejemplo '#FFFFFF' o 'none' |
| `--journal {acs,aps,elsevier,generic,iop,nature,rsc,wiley}` | ancho de figura según la revista (default: generic) |
| `--width` | ancho de la figura: single, onehalf o double, o un número en milímetros |
| `--mono` | versión en escala de grises: tinta negra y patrones de línea |

**Fundamento físico:** [`olla-dft teoria kappa`](TEORIA.md)

### `qha`

cuasi-armónica: expansión térmica y a(T)

**Uso:** `olla-dft qha [-h] [-o OUTDIR] [--natoms NATOMS] [--cells CELLS] [--cubic] [--structure CIF] [--tmax TMAX] [--dt DT] [--temp TEMP] [--dpi DPI] [--format FORMAT] [--no-plot] [-t TEMPLATE] [--font {sans,serif,latex}] [--usetex] [--palette PALETTE] [--background BACKGROUND] [--journal {acs,aps,elsevier,generic,iop,nature,rsc,wiley}] [--width WIDTH] [--mono] data`

**Argumentos:**

- `data` — tabla: V(A^3) E(eV) w1 w2 ... por volumen

**Opciones:**

| Opción | Descripción |
|---|---|
| `-o, --outdir` | carpeta de salida (default: `.`) |
| `--natoms NATOMS` | átomos por celda primitiva, para las magnitudes por átomo (default: 1) |
| `--cells CELLS` | celdas primitivas por supercelda de los modos (default: `1`) |
| `--cubic` | además, a(T). Sin --structure es solo V_prim^(1/3) |
| `--structure CIF` | estructura del material: con ella a(T) se convierte al parámetro de red CONVENCIONAL (factor 4 en fcc/diamante, 2 en bcc) y se detecta si es cúbica |
| `--tmax TMAX` | temperatura máxima de la malla de T en K (default: 1000) |
| `--dt DT` | paso: de integración en fs (amorphous) o de la malla de temperaturas en K (qha) (default: `5.0`) |
| `--temp TEMP` | temperatura de trabajo en K (default: `300.0`) |
| `--dpi DPI` | resolución de las figuras de mapa de bits en puntos por pulgada (default: 600) |
| `--format` | formatos de figura separados por coma: pdf,png,svg,eps,tif (default: pdf,png) |
| `--no-plot` | solo exportar los datos, sin generar la figura |
| `-t, --template` | plantilla visual de la figura (la lista sale con 'olla-dft templates list') |
| `--font {sans,serif,latex}` | familia tipográfica: sans, serif o latex (Computer Modern) |
| `--usetex` | componer los textos de la figura con LaTeX real (necesita una instalación de LaTeX) |
| `--palette` | paleta de colores: grayscale, okabe-ito, okabe-ito-dark o colores hex separados por coma |
| `--background` | color de fondo de la figura, por ejemplo '#FFFFFF' o 'none' |
| `--journal {acs,aps,elsevier,generic,iop,nature,rsc,wiley}` | ancho de figura según la revista (default: generic) |
| `--width` | ancho de la figura: single, onehalf o double, o un número en milímetros |
| `--mono` | versión en escala de grises: tinta negra y patrones de línea |

**Fundamento físico:** [`olla-dft teoria qha`](TEORIA.md)

### `thermochem`

ZPE, entropia y energia libre: de una energia DFT a una comparable con el experimento

**Uso:** `olla-dft thermochem [-h] [--phase {solido,adsorbato,gas,transicion}] [--structure STRUCTURE] [--temp TEMP] [--pressure PRESSURE] [--symmetry SYMMETRY] [--multiplicity MULTIPLICITY] [--floor FLOOR] [--energy ENERGY] [-o OUTDIR] freqs`

**Argumentos:**

- `freqs` — archivo de frecuencias en cm-1, o la lista separada por comas

**Opciones:**

| Opción | Descripción |
|---|---|
| `--phase {solido,adsorbato,gas,transicion}` | gas anade traslaciones y rotaciones; transicion exige exactamente una frecuencia imaginaria (default: `solido`) |
| `--structure` | estructura (necesaria para la fase gas) |
| `--temp TEMP` | temperatura de trabajo en K (default: `298.15`) |
| `--pressure PRESSURE` | en bar (default: `1.0`) |
| `--symmetry SYMMETRY` | numero de simetria del grupo puntual: 2 para H2O y O2, 3 para NH3, 12 para CH4 (default: `1`) |
| `--multiplicity MULTIPLICITY` | multiplicidad de espin del estado fundamental (default: `1`) |
| `--floor FLOOR` | sube los modos por debajo de este valor (cm-1); 100 es lo habitual |
| `--energy ENERGY` | E_DFT en eV, para dar G(T) |
| `-o, --outdir` | carpeta de salida |

**Fundamento físico:** [`olla-dft teoria thermochem`](TEORIA.md)

### `md`

analizar una trayectoria de dinamica molecular: g(r), difusion y espectro vibracional

**Uso:** `olla-dft md [-h] [-o OUTDIR] [--skip SKIP] [--rmax RMAX] [--bins BINS] [--dpi DPI] [--format FORMAT] [--no-plot] [-t TEMPLATE] [--font {sans,serif,latex}] [--usetex] [--palette PALETTE] [--background BACKGROUND] [--journal {acs,aps,elsevier,generic,iop,nature,rsc,wiley}] [--width WIDTH] [--mono] path`

**Argumentos:**

- `path` — salida de pw.x con calculation='md', o su carpeta

**Opciones:**

| Opción | Descripción |
|---|---|
| `-o, --outdir` | carpeta de salida (default: `.`) |
| `--skip SKIP` | pasos iniciales a descartar (equilibrado) |
| `--rmax RMAX` | corte de g(r) en A; por omision, media arista de la celda, que es hasta donde la normalizacion vale |
| `--bins BINS` | número de intervalos del histograma de g(r) (default: 200) |
| `--dpi DPI` | resolución de las figuras de mapa de bits en puntos por pulgada (default: 600) |
| `--format` | formatos de figura separados por coma: pdf,png,svg,eps,tif (default: pdf,png) |
| `--no-plot` | solo exportar los datos, sin generar la figura |
| `-t, --template` | plantilla visual de la figura (la lista sale con 'olla-dft templates list') |
| `--font {sans,serif,latex}` | familia tipográfica: sans, serif o latex (Computer Modern) |
| `--usetex` | componer los textos de la figura con LaTeX real (necesita una instalación de LaTeX) |
| `--palette` | paleta de colores: grayscale, okabe-ito, okabe-ito-dark o colores hex separados por coma |
| `--background` | color de fondo de la figura, por ejemplo '#FFFFFF' o 'none' |
| `--journal {acs,aps,elsevier,generic,iop,nature,rsc,wiley}` | ancho de figura según la revista (default: generic) |
| `--width` | ancho de la figura: single, onehalf o double, o un número en milímetros |
| `--mono` | versión en escala de grises: tinta negra y patrones de línea |

**Fundamento físico:** [`olla-dft teoria md`](TEORIA.md)

### `derived`

Debye, velocidades del sonido y Slack desde las Cij

**Uso:** `olla-dft derived [-h] [--cij CIJ] [--temp TEMP] [-o OUTDIR] file`

**Argumentos:**

- `file` — estructura

**Opciones:**

| Opción | Descripción |
|---|---|
| `--cij` | archivo con la matriz elástica (default: `ELASTIC_C.dat`) |
| `--temp TEMP` | temperatura de trabajo en K (default: `300.0`) |
| `-o, --outdir` | carpeta donde dejar DERIVED.dat (default: `.`) |

**Fundamento físico:** [`olla-dft teoria derived`](TEORIA.md)

## Mecánica y estabilidad

### `converge`

pruebas de convergencia de cutoffs y malla k

**Uso:** `olla-dft converge [-h] [-o OUTDIR] [--run] [--collect] [--pw-cmd PW_CMD] [--nproc NPROC] [-j N] [--redo] [--max-time T] [--estimate] [--timeout TIMEOUT] [--pseudo-dir PSEUDO_DIR] [--pseudo EL=UPF] [--ecutwfc ECUTWFC] [--ecutrho ECUTRHO] [--kspacing KSPACING] [--insulator] [--dpi DPI] [--format FORMAT] [--no-plot] [-t TEMPLATE] [--size {paper,poster,presentation}] [--font {sans,serif,latex}] [--usetex] [--palette PALETTE] [--background BACKGROUND] [--journal {acs,aps,elsevier,generic,iop,nature,rsc,wiley}] [--width WIDTH] [--aspect ASPECT] [--mono] [-k {ecutwfc,ecutrho,kmesh}] [--values VALUES] [--threshold THRESHOLD] file`

**Argumentos:**

- `file` — estructura (CIF, POSCAR, input de pw.x...)

**Opciones:**

| Opción | Descripción |
|---|---|
| `-o, --outdir` | carpeta del barrido (default: `convergencia`) |
| `--run` | ejecutar los cálculos ahora, uno tras otro |
| `--collect` | solo analizar cálculos ya corridos |
| `--pw-cmd` | ejecutable de pw.x (anula la configuración) |
| `--nproc NPROC` | procesos MPI por cálculo |
| `-j, --jobs N` | cálculos simultáneos (default: 1). Sin --nproc, los hilos de la máquina se reparten entre ellos |
| `--redo` | rehacer también los cálculos que ya estaban terminados |
| `--max-time T` | presupuesto TOTAL de tiempo: 90m, 2h, 3600. Al agotarse no se lanzan más y el barrido queda reanudable |
| `--estimate` | estimar cuánto va a tardar el barrido y salir, usando el histórico de 'olla-dft db' |
| `--timeout TIMEOUT` | límite en segundos por cálculo |
| `--pseudo-dir` | carpeta de pseudopotenciales |
| `--pseudo EL=UPF` | forzar un pseudopotencial concreto, por ejemplo Fe=Fe.rel-pbe.UPF. Se puede repetir. Sin esto, Olla-DFT elige con 'olla-dft pseudos' |
| `--ecutwfc ECUTWFC` | cutoff de ondas (Ry) |
| `--ecutrho ECUTRHO` | cutoff de densidad (Ry) |
| `--kspacing KSPACING` | espaciado k en Å^-1 |
| `--insulator` | occupations='fixed' |
| `--dpi DPI` | resolución de las figuras de mapa de bits en puntos por pulgada (default: 600) |
| `--format` | formatos de figura separados por coma: pdf,png,svg,eps,tif (default: pdf,png) |
| `--no-plot` | solo exportar los datos, sin generar la figura |
| `-t, --template` | plantilla visual de la figura |
| `--size {paper,poster,presentation}` | tamaño de figura: paper, presentation o poster |
| `--font {sans,serif,latex}` | familia tipográfica: sans, serif o latex (Computer Modern) |
| `--usetex` | componer los textos de la figura con LaTeX real (necesita una instalación de LaTeX) |
| `--palette` | paleta de colores: grayscale, okabe-ito, okabe-ito-dark o colores hex separados por coma |
| `--background` | color de fondo de la figura, por ejemplo '#FFFFFF' o 'none' |
| `--journal {acs,aps,elsevier,generic,iop,nature,rsc,wiley}` | ancho de figura según la revista (default: generic) |
| `--width` | ancho de la figura: single, onehalf o double, o un número en milímetros |
| `--aspect ASPECT` | relación alto/ancho de la figura |
| `--mono` | versión en escala de grises: tinta negra y patrones de línea |
| `-k, --kind {ecutwfc,ecutrho,kmesh}` | qué parámetro se barre (default: ecutwfc) |
| `--values` | valores separados por coma; para kmesh admite 8x8x8 |
| `--threshold THRESHOLD` | umbral de convergencia en meV/átomo (default: 1) |

**Fundamento físico:** [`olla-dft teoria converge`](TEORIA.md)

### `eos`

ecuación de estado E–V y módulo de bulk

**Uso:** `olla-dft eos [-h] [-o OUTDIR] [--run] [--collect] [--pw-cmd PW_CMD] [--nproc NPROC] [-j N] [--redo] [--max-time T] [--estimate] [--timeout TIMEOUT] [--pseudo-dir PSEUDO_DIR] [--pseudo EL=UPF] [--ecutwfc ECUTWFC] [--ecutrho ECUTRHO] [--kspacing KSPACING] [--insulator] [--dpi DPI] [--format FORMAT] [--no-plot] [-t TEMPLATE] [--size {paper,poster,presentation}] [--font {sans,serif,latex}] [--usetex] [--palette PALETTE] [--background BACKGROUND] [--journal {acs,aps,elsevier,generic,iop,nature,rsc,wiley}] [--width WIDTH] [--aspect ASPECT] [--mono] [--npoints NPOINTS] [--scale SCALE] [--span SPAN] [--equation {birch-murnaghan,murnaghan,vinet}] [--relax-ions] file`

**Argumentos:**

- `file` — estructura (CIF, POSCAR, input de pw.x...)

**Opciones:**

| Opción | Descripción |
|---|---|
| `-o, --outdir` | carpeta del barrido (default: `eos`) |
| `--run` | ejecutar los cálculos ahora, uno tras otro |
| `--collect` | solo analizar cálculos ya corridos |
| `--pw-cmd` | ejecutable de pw.x (anula la configuración) |
| `--nproc NPROC` | procesos MPI por cálculo |
| `-j, --jobs N` | cálculos simultáneos (default: 1). Sin --nproc, los hilos de la máquina se reparten entre ellos |
| `--redo` | rehacer también los cálculos que ya estaban terminados |
| `--max-time T` | presupuesto TOTAL de tiempo: 90m, 2h, 3600. Al agotarse no se lanzan más y el barrido queda reanudable |
| `--estimate` | estimar cuánto va a tardar el barrido y salir, usando el histórico de 'olla-dft db' |
| `--timeout TIMEOUT` | límite en segundos por cálculo |
| `--pseudo-dir` | carpeta de pseudopotenciales |
| `--pseudo EL=UPF` | forzar un pseudopotencial concreto, por ejemplo Fe=Fe.rel-pbe.UPF. Se puede repetir. Sin esto, Olla-DFT elige con 'olla-dft pseudos' |
| `--ecutwfc ECUTWFC` | cutoff de ondas (Ry) |
| `--ecutrho ECUTRHO` | cutoff de densidad (Ry) |
| `--kspacing KSPACING` | espaciado k en Å^-1 |
| `--insulator` | occupations='fixed' |
| `--dpi DPI` | resolución de las figuras de mapa de bits en puntos por pulgada (default: 600) |
| `--format` | formatos de figura separados por coma: pdf,png,svg,eps,tif (default: pdf,png) |
| `--no-plot` | solo exportar los datos, sin generar la figura |
| `-t, --template` | plantilla visual de la figura |
| `--size {paper,poster,presentation}` | tamaño de figura: paper, presentation o poster |
| `--font {sans,serif,latex}` | familia tipográfica: sans, serif o latex (Computer Modern) |
| `--usetex` | componer los textos de la figura con LaTeX real (necesita una instalación de LaTeX) |
| `--palette` | paleta de colores: grayscale, okabe-ito, okabe-ito-dark o colores hex separados por coma |
| `--background` | color de fondo de la figura, por ejemplo '#FFFFFF' o 'none' |
| `--journal {acs,aps,elsevier,generic,iop,nature,rsc,wiley}` | ancho de figura según la revista (default: generic) |
| `--width` | ancho de la figura: single, onehalf o double, o un número en milímetros |
| `--aspect ASPECT` | relación alto/ancho de la figura |
| `--mono` | versión en escala de grises: tinta negra y patrones de línea |
| `--npoints NPOINTS` | número de volúmenes (default: 9) |
| `--scale SCALE` | factor lineal en que centrar el barrido (lo devuelve 'olla-dft mlip scan') (default: `1.0`) |
| `--span SPAN` | variación relativa de volumen a cada lado (default: 0.10) |
| `--equation {birch-murnaghan,murnaghan,vinet}` | ecuación que se grafica (default: `birch-murnaghan`) |
| `--relax-ions` | relajar posiciones internas en cada volumen |

**Fundamento físico:** [`olla-dft teoria eos`](TEORIA.md)

### `elastic`

constantes elásticas y propiedades mecánicas

**Uso:** `olla-dft elastic [-h] [-o OUTDIR] [--run] [--collect] [--pw-cmd PW_CMD] [--nproc NPROC] [-j N] [--redo] [--max-time T] [--estimate] [--timeout TIMEOUT] [--pseudo-dir PSEUDO_DIR] [--pseudo EL=UPF] [--ecutwfc ECUTWFC] [--ecutrho ECUTRHO] [--kspacing KSPACING] [--insulator] [--dpi DPI] [--format FORMAT] [--no-plot] [-t TEMPLATE] [--size {paper,poster,presentation}] [--font {sans,serif,latex}] [--usetex] [--palette PALETTE] [--background BACKGROUND] [--journal {acs,aps,elsevier,generic,iop,nature,rsc,wiley}] [--width WIDTH] [--aspect ASPECT] [--mono] [--delta DELTA] [--npoints NPOINTS] [--2d] [--thickness A] [--ion-mode {auto,relax,fixed}] file`

**Argumentos:**

- `file` — estructura (CIF, POSCAR, input de pw.x...)

**Opciones:**

| Opción | Descripción |
|---|---|
| `-o, --outdir` | carpeta del barrido (default: `elastic`) |
| `--run` | ejecutar los cálculos ahora, uno tras otro |
| `--collect` | solo analizar cálculos ya corridos |
| `--pw-cmd` | ejecutable de pw.x (anula la configuración) |
| `--nproc NPROC` | procesos MPI por cálculo |
| `-j, --jobs N` | cálculos simultáneos (default: 1). Sin --nproc, los hilos de la máquina se reparten entre ellos |
| `--redo` | rehacer también los cálculos que ya estaban terminados |
| `--max-time T` | presupuesto TOTAL de tiempo: 90m, 2h, 3600. Al agotarse no se lanzan más y el barrido queda reanudable |
| `--estimate` | estimar cuánto va a tardar el barrido y salir, usando el histórico de 'olla-dft db' |
| `--timeout TIMEOUT` | límite en segundos por cálculo |
| `--pseudo-dir` | carpeta de pseudopotenciales |
| `--pseudo EL=UPF` | forzar un pseudopotencial concreto, por ejemplo Fe=Fe.rel-pbe.UPF. Se puede repetir. Sin esto, Olla-DFT elige con 'olla-dft pseudos' |
| `--ecutwfc ECUTWFC` | cutoff de ondas (Ry) |
| `--ecutrho ECUTRHO` | cutoff de densidad (Ry) |
| `--kspacing KSPACING` | espaciado k en Å^-1 |
| `--insulator` | occupations='fixed' |
| `--dpi DPI` | resolución de las figuras de mapa de bits en puntos por pulgada (default: 600) |
| `--format` | formatos de figura separados por coma: pdf,png,svg,eps,tif (default: pdf,png) |
| `--no-plot` | solo exportar los datos, sin generar la figura |
| `-t, --template` | plantilla visual de la figura |
| `--size {paper,poster,presentation}` | tamaño de figura: paper, presentation o poster |
| `--font {sans,serif,latex}` | familia tipográfica: sans, serif o latex (Computer Modern) |
| `--usetex` | componer los textos de la figura con LaTeX real (necesita una instalación de LaTeX) |
| `--palette` | paleta de colores: grayscale, okabe-ito, okabe-ito-dark o colores hex separados por coma |
| `--background` | color de fondo de la figura, por ejemplo '#FFFFFF' o 'none' |
| `--journal {acs,aps,elsevier,generic,iop,nature,rsc,wiley}` | ancho de figura según la revista (default: generic) |
| `--width` | ancho de la figura: single, onehalf o double, o un número en milímetros |
| `--aspect ASPECT` | relación alto/ancho de la figura |
| `--mono` | versión en escala de grises: tinta negra y patrones de línea |
| `--delta DELTA` | deformación máxima aplicada (default: 0.010 = 1 %%) |
| `--npoints NPOINTS` | deformaciones no nulas por componente, par (default: 4) |
| `--2d` | lámina: constantes en N/m (no en GPa), solo ε1, ε2 y ε6, y criterios de Born en 2D |
| `--thickness A` | espesor supuesto en Å para dar también el equivalente en GPa (convenio, no medida) |
| `--ion-mode {auto,relax,fixed}` | posiciones internas: auto = fijas en deformaciones normales y relajadas en cizallas (recomendado); relax = relajar todas; fixed = clamped-ion (default: `auto`) |

**Fundamento físico:** [`olla-dft teoria elastic`](TEORIA.md)

### `strain`

barrido de deformación: gap, energía y momento en función de la deformación aplicada

**Uso:** `olla-dft strain [-h] [-o OUTDIR] [--run] [--collect] [--pw-cmd PW_CMD] [--nproc NPROC] [-j N] [--redo] [--max-time T] [--estimate] [--timeout TIMEOUT] [--pseudo-dir PSEUDO_DIR] [--pseudo EL=UPF] [--ecutwfc ECUTWFC] [--ecutrho ECUTRHO] [--kspacing KSPACING] [--insulator] [--dpi DPI] [--format FORMAT] [--no-plot] [-t TEMPLATE] [--size {paper,poster,presentation}] [--font {sans,serif,latex}] [--usetex] [--palette PALETTE] [--background BACKGROUND] [--journal {acs,aps,elsevier,generic,iop,nature,rsc,wiley}] [--width WIDTH] [--aspect ASPECT] [--mono] [-m {biaxial,cizalla,hidrostatica,uniaxial-a,uniaxial-b,uniaxial-c}] [-r MIN:MAX:N] [--fixed-ions] [--relax-perp] [--nspin {1,2}] [--mag MAG] [--hubbard EL=U] [--vdw {grimme-d2,grimme-d3,DFT-D,ts-vdw,xdm,mbd}] file`

**Argumentos:**

- `file` — estructura (CIF, POSCAR, input de pw.x...)

**Opciones:**

| Opción | Descripción |
|---|---|
| `-o, --outdir` | carpeta del barrido (default: `strain`) |
| `--run` | ejecutar los cálculos ahora, uno tras otro |
| `--collect` | solo analizar cálculos ya corridos |
| `--pw-cmd` | ejecutable de pw.x (anula la configuración) |
| `--nproc NPROC` | procesos MPI por cálculo |
| `-j, --jobs N` | cálculos simultáneos (default: 1). Sin --nproc, los hilos de la máquina se reparten entre ellos |
| `--redo` | rehacer también los cálculos que ya estaban terminados |
| `--max-time T` | presupuesto TOTAL de tiempo: 90m, 2h, 3600. Al agotarse no se lanzan más y el barrido queda reanudable |
| `--estimate` | estimar cuánto va a tardar el barrido y salir, usando el histórico de 'olla-dft db' |
| `--timeout TIMEOUT` | límite en segundos por cálculo |
| `--pseudo-dir` | carpeta de pseudopotenciales |
| `--pseudo EL=UPF` | forzar un pseudopotencial concreto, por ejemplo Fe=Fe.rel-pbe.UPF. Se puede repetir. Sin esto, Olla-DFT elige con 'olla-dft pseudos' |
| `--ecutwfc ECUTWFC` | cutoff de ondas (Ry) |
| `--ecutrho ECUTRHO` | cutoff de densidad (Ry) |
| `--kspacing KSPACING` | espaciado k en Å^-1 |
| `--insulator` | occupations='fixed' |
| `--dpi DPI` | resolución de las figuras de mapa de bits en puntos por pulgada (default: 600) |
| `--format` | formatos de figura separados por coma: pdf,png,svg,eps,tif (default: pdf,png) |
| `--no-plot` | solo exportar los datos, sin generar la figura |
| `-t, --template` | plantilla visual de la figura |
| `--size {paper,poster,presentation}` | tamaño de figura: paper, presentation o poster |
| `--font {sans,serif,latex}` | familia tipográfica: sans, serif o latex (Computer Modern) |
| `--usetex` | componer los textos de la figura con LaTeX real (necesita una instalación de LaTeX) |
| `--palette` | paleta de colores: grayscale, okabe-ito, okabe-ito-dark o colores hex separados por coma |
| `--background` | color de fondo de la figura, por ejemplo '#FFFFFF' o 'none' |
| `--journal {acs,aps,elsevier,generic,iop,nature,rsc,wiley}` | ancho de figura según la revista (default: generic) |
| `--width` | ancho de la figura: single, onehalf o double, o un número en milímetros |
| `--aspect ASPECT` | relación alto/ancho de la figura |
| `--mono` | versión en escala de grises: tinta negra y patrones de línea |
| `-m, --mode {biaxial,cizalla,hidrostatica,uniaxial-a,uniaxial-b,uniaxial-c}` | qué se deforma (default: biaxial) |
| `-r, --range MIN:MAX:N` | rango en POR CIENTO, por ejemplo -5:5:11 (de -5 %% a +5 %% en 11 puntos) (default: `-5:5:11`) |
| `--fixed-ions` | no relajar las posiciones internas en cada deformación (más rápido y menos realista) |
| `--relax-perp` | dejar libre el eje perpendicular al plano deformado (relajación de Poisson); imprescindible en láminas |
| `--nspin {1,2}` | 2 activa la polarización de espín (default: `1`) |
| `--mag` | magnetización inicial (implica --nspin 2) |
| `--hubbard EL=U` | U de Hubbard en eV por elemento |
| `--vdw {grimme-d2,grimme-d3,DFT-D,ts-vdw,xdm,mbd}` | corrección de dispersión |

**Fundamento físico:** [`olla-dft teoria strain`](TEORIA.md)

### `layers`

detectar capas, espaciado basal y hueco interlaminar

**Uso:** `olla-dft layers [-h] [--tol TOL] [--wavelength WAVELENGTH] [--slab ARCHIVO] [--vacuum VACUUM] file`

**Argumentos:**

- `file` — estructura de entrada (CIF, POSCAR, input de pw.x...)

**Opciones:**

| Opción | Descripción |
|---|---|
| `--tol TOL` | tolerancia de enlace sobre radios covalentes (Å) (default: `0.45`) |
| `--wavelength` | radiación para las reflexiones basales (default CuKa) |
| `--slab ARCHIVO` | además, escribir la monocapa con vacío a este archivo |
| `--vacuum VACUUM` | vacío de la monocapa en Å (default 20) |

**Fundamento físico:** [`olla-dft teoria layers`](TEORIA.md)

### `xrd`

difractograma de polvos simulado

**Uso:** `olla-dft xrd [-h] [-o OUTDIR] [--suite] [--basis {conventional,input}] [--wavelength WAVELENGTH] [--tt-min TT_MIN] [--tt-max TT_MAX] [--fwhm FWHM] [--size SIZE] [--biso BISO] [--exp EXP] [--dpi DPI] [--format FORMAT] [--no-plot] [-t TEMPLATE] [--size-preset {paper,poster,presentation}] [--font {sans,serif,latex}] [--usetex] [--palette PALETTE] [--background BACKGROUND] [--journal {acs,aps,elsevier,generic,iop,nature,rsc,wiley}] [--width WIDTH] [--aspect ASPECT] [--mono] file`

**Argumentos:**

- `file` — estructura de entrada (CIF, POSCAR, input de pw.x...)

**Opciones:**

| Opción | Descripción |
|---|---|
| `-o, --outdir` | carpeta de salida (default: `.`) |
| `--suite` | además, exportar JSON de intercambio para las otras apps de la suite |
| `--basis {conventional,input}` | celda en que se indexan los hkl: 'conventional' (por defecto, los índices de las fichas PDF) o 'input' (la celda del archivo tal cual) (default: `conventional`) |
| `--wavelength` | radiación: AgKa, CoKa, CrKa, CuKa, CuKa1, FeKa, MoKa o λ en Å (default: `CuKa`) |
| `--tt-min TT_MIN` | 2θ mínimo (°) (default: `5.0`) |
| `--tt-max TT_MAX` | 2θ máximo (°) (default: `70.0`) |
| `--fwhm FWHM` | anchura instrumental (° 2θ, default 0.15) |
| `--size SIZE` | tamaño de cristalito en nm (anchura por Scherrer) |
| `--biso BISO` | factor de temperatura global B (Å²) |
| `--exp` | difractograma experimental (2θ, I) para superponer |
| `--dpi DPI` | resolución de las figuras de mapa de bits en puntos por pulgada (default: 600) |
| `--format` | formatos de figura separados por coma: pdf,png,svg,eps,tif (default: pdf,png) |
| `--no-plot` | solo exportar los datos, sin generar la figura |
| `-t, --template` | plantilla visual de la figura (la lista sale con 'olla-dft templates list') |
| `--size-preset {paper,poster,presentation}` | escala tipográfica de la figura |
| `--font {sans,serif,latex}` | familia tipográfica: sans, serif o latex (Computer Modern) |
| `--usetex` | componer los textos de la figura con LaTeX real (necesita una instalación de LaTeX) |
| `--palette` | paleta de colores: grayscale, okabe-ito, okabe-ito-dark o colores hex separados por coma |
| `--background` | color de fondo de la figura, por ejemplo '#FFFFFF' o 'none' |
| `--journal {acs,aps,elsevier,generic,iop,nature,rsc,wiley}` | ancho de figura según la revista (default: generic) |
| `--width` | ancho de la figura: single, onehalf o double, o un número en milímetros |
| `--aspect ASPECT` | relación alto/ancho de la figura |
| `--mono` | versión en escala de grises: tinta negra y patrones de línea |

**Fundamento físico:** [`olla-dft teoria xrd`](TEORIA.md)

### `exfoliate`

energía de exfoliación (bulk vs monocapa)

**Uso:** `olla-dft exfoliate [-h] [-o OUTDIR] [--run] [--collect] [--pw-cmd PW_CMD] [--nproc NPROC] [-j N] [--redo] [--max-time T] [--estimate] [--timeout TIMEOUT] [--pseudo-dir PSEUDO_DIR] [--pseudo EL=UPF] [--ecutwfc ECUTWFC] [--ecutrho ECUTRHO] [--kspacing KSPACING] [--insulator] [--dpi DPI] [--format FORMAT] [--no-plot] [-t TEMPLATE] [--size {paper,poster,presentation}] [--font {sans,serif,latex}] [--usetex] [--palette PALETTE] [--background BACKGROUND] [--journal {acs,aps,elsevier,generic,iop,nature,rsc,wiley}] [--width WIDTH] [--aspect ASPECT] [--mono] [--vacuum VACUUM] [--vdw {grimme-d2,grimme-d3,DFT-D,ts-vdw,xdm,mbd}] [--tol TOL] [--relax-slab] file`

**Argumentos:**

- `file` — estructura (CIF, POSCAR, input de pw.x...)

**Opciones:**

| Opción | Descripción |
|---|---|
| `-o, --outdir` | carpeta del barrido (default: `exfoliacion`) |
| `--run` | ejecutar los cálculos ahora, uno tras otro |
| `--collect` | solo analizar cálculos ya corridos |
| `--pw-cmd` | ejecutable de pw.x (anula la configuración) |
| `--nproc NPROC` | procesos MPI por cálculo |
| `-j, --jobs N` | cálculos simultáneos (default: 1). Sin --nproc, los hilos de la máquina se reparten entre ellos |
| `--redo` | rehacer también los cálculos que ya estaban terminados |
| `--max-time T` | presupuesto TOTAL de tiempo: 90m, 2h, 3600. Al agotarse no se lanzan más y el barrido queda reanudable |
| `--estimate` | estimar cuánto va a tardar el barrido y salir, usando el histórico de 'olla-dft db' |
| `--timeout TIMEOUT` | límite en segundos por cálculo |
| `--pseudo-dir` | carpeta de pseudopotenciales |
| `--pseudo EL=UPF` | forzar un pseudopotencial concreto, por ejemplo Fe=Fe.rel-pbe.UPF. Se puede repetir. Sin esto, Olla-DFT elige con 'olla-dft pseudos' |
| `--ecutwfc ECUTWFC` | cutoff de ondas (Ry) |
| `--ecutrho ECUTRHO` | cutoff de densidad (Ry) |
| `--kspacing KSPACING` | espaciado k en Å^-1 |
| `--insulator` | occupations='fixed' |
| `--dpi DPI` | resolución de las figuras de mapa de bits en puntos por pulgada (default: 600) |
| `--format` | formatos de figura separados por coma: pdf,png,svg,eps,tif (default: pdf,png) |
| `--no-plot` | solo exportar los datos, sin generar la figura |
| `-t, --template` | plantilla visual de la figura |
| `--size {paper,poster,presentation}` | tamaño de figura: paper, presentation o poster |
| `--font {sans,serif,latex}` | familia tipográfica: sans, serif o latex (Computer Modern) |
| `--usetex` | componer los textos de la figura con LaTeX real (necesita una instalación de LaTeX) |
| `--palette` | paleta de colores: grayscale, okabe-ito, okabe-ito-dark o colores hex separados por coma |
| `--background` | color de fondo de la figura, por ejemplo '#FFFFFF' o 'none' |
| `--journal {acs,aps,elsevier,generic,iop,nature,rsc,wiley}` | ancho de figura según la revista (default: generic) |
| `--width` | ancho de la figura: single, onehalf o double, o un número en milímetros |
| `--aspect ASPECT` | relación alto/ancho de la figura |
| `--mono` | versión en escala de grises: tinta negra y patrones de línea |
| `--vacuum VACUUM` | vacío de la monocapa en Å (default 20) |
| `--vdw {grimme-d2,grimme-d3,DFT-D,ts-vdw,xdm,mbd}` | corrección de dispersión para ambos cálculos |
| `--tol TOL` | tolerancia de enlace para detectar las capas (Å) (default: `0.45`) |
| `--relax-slab` | relajar las posiciones de la monocapa |

**Fundamento físico:** [`olla-dft teoria exfoliate`](TEORIA.md)

### `gamma`

energía de superficie y de escisión por el ajuste lineal de Fiorentini–Methfessel, con la convergencia contra el grosor de la losa

**Uso:** `olla-dft gamma [-h] [-o OUTDIR] [--run] [--collect] [--pw-cmd PW_CMD] [--nproc NPROC] [-j N] [--redo] [--max-time T] [--estimate] [--timeout TIMEOUT] [--pseudo-dir PSEUDO_DIR] [--pseudo EL=UPF] [--ecutwfc ECUTWFC] [--ecutrho ECUTRHO] [--kspacing KSPACING] [--insulator] [--dpi DPI] [--format FORMAT] [--no-plot] [-t TEMPLATE] [--size {paper,poster,presentation}] [--font {sans,serif,latex}] [--usetex] [--palette PALETTE] [--background BACKGROUND] [--journal {acs,aps,elsevier,generic,iop,nature,rsc,wiley}] [--width WIDTH] [--aspect ASPECT] [--mono] [-m MILLER] [-l LAYERS] [--vacuum VACUUM] [--fix N] [--relax] [--no-bulk] [--no-reduce] [--vdw {grimme-d2,grimme-d3,DFT-D,ts-vdw,xdm,mbd}] [--dipole] [--nspin {1,2}] [--mag MAG] file`

**Argumentos:**

- `file` — estructura (CIF, POSCAR, input de pw.x...)

**Opciones:**

| Opción | Descripción |
|---|---|
| `-o, --outdir` | carpeta del barrido (default: `gamma`) |
| `--run` | ejecutar los cálculos ahora, uno tras otro |
| `--collect` | solo analizar cálculos ya corridos |
| `--pw-cmd` | ejecutable de pw.x (anula la configuración) |
| `--nproc NPROC` | procesos MPI por cálculo |
| `-j, --jobs N` | cálculos simultáneos (default: 1). Sin --nproc, los hilos de la máquina se reparten entre ellos |
| `--redo` | rehacer también los cálculos que ya estaban terminados |
| `--max-time T` | presupuesto TOTAL de tiempo: 90m, 2h, 3600. Al agotarse no se lanzan más y el barrido queda reanudable |
| `--estimate` | estimar cuánto va a tardar el barrido y salir, usando el histórico de 'olla-dft db' |
| `--timeout TIMEOUT` | límite en segundos por cálculo |
| `--pseudo-dir` | carpeta de pseudopotenciales |
| `--pseudo EL=UPF` | forzar un pseudopotencial concreto, por ejemplo Fe=Fe.rel-pbe.UPF. Se puede repetir. Sin esto, Olla-DFT elige con 'olla-dft pseudos' |
| `--ecutwfc ECUTWFC` | cutoff de ondas (Ry) |
| `--ecutrho ECUTRHO` | cutoff de densidad (Ry) |
| `--kspacing KSPACING` | espaciado k en Å^-1 |
| `--insulator` | occupations='fixed' |
| `--dpi DPI` | resolución de las figuras de mapa de bits en puntos por pulgada (default: 600) |
| `--format` | formatos de figura separados por coma: pdf,png,svg,eps,tif (default: pdf,png) |
| `--no-plot` | solo exportar los datos, sin generar la figura |
| `-t, --template` | plantilla visual de la figura |
| `--size {paper,poster,presentation}` | tamaño de figura: paper, presentation o poster |
| `--font {sans,serif,latex}` | familia tipográfica: sans, serif o latex (Computer Modern) |
| `--usetex` | componer los textos de la figura con LaTeX real (necesita una instalación de LaTeX) |
| `--palette` | paleta de colores: grayscale, okabe-ito, okabe-ito-dark o colores hex separados por coma |
| `--background` | color de fondo de la figura, por ejemplo '#FFFFFF' o 'none' |
| `--journal {acs,aps,elsevier,generic,iop,nature,rsc,wiley}` | ancho de figura según la revista (default: generic) |
| `--width` | ancho de la figura: single, onehalf o double, o un número en milímetros |
| `--aspect ASPECT` | relación alto/ancho de la figura |
| `--mono` | versión en escala de grises: tinta negra y patrones de línea |
| `-m, --miller` | índices de Miller de la cara, por ejemplo '1 1 1' (default: `1 0 0`) |
| `-l, --layers` | grosores a calcular, separados por coma (default: 3,4,5,6). Hacen falta al menos dos |
| `--vacuum VACUUM` | vacío en Å (default: 20) |
| `--fix N` | congelar N capas del fondo al relajar |
| `--relax` | relajar las posiciones (γ baja entre un 5 y un 20 %%) |
| `--no-bulk` | no calcular el bulto aparte; solo el ajuste lineal E_losa(N) = 2γA + N·E_bulto |
| `--no-reduce` | no reducir la celda superficial a la mínima (por omisión sí se reduce: mismo γ, mucho más barato) |
| `--vdw {grimme-d2,grimme-d3,DFT-D,ts-vdw,xdm,mbd}` | corrección de dispersión (van der Waals): grimme-d2, grimme-d3, DFT-D, ts-vdw, xdm o mbd |
| `--dipole` | corrección dipolar, para losas polares |
| `--nspin {1,2}` | 2 activa la polarización de espín (default: 1) |
| `--mag` | magnetización inicial (implica --nspin 2) |

**Fundamento físico:** [`olla-dft teoria gamma`](TEORIA.md)

## Superficies, defectos y química

### `surface`

cortar una superficie (hkl) con vacío

**Uso:** `olla-dft surface [-h] [-m MILLER] [-l LAYERS] [--vacuum VACUUM] [--fix FIX] [-o OUTPUT] file`

**Argumentos:**

- `file` — estructura de entrada (CIF, POSCAR, input de pw.x...)

**Opciones:**

| Opción | Descripción |
|---|---|
| `-m, --miller` | índices de Miller, por ejemplo '1 1 1' o 1,1,1 (default: `1 0 0`) |
| `-l, --layers LAYERS` | número de capas atómicas de la losa (default: 6) |
| `--vacuum VACUUM` | vacío total en Å (default: `15.0`) |
| `--fix FIX` | planos atómicos del fondo a congelar |
| `-o, --output` | archivo de salida (CIF/POSCAR) |

**Fundamento físico:** [`olla-dft teoria surface`](TEORIA.md)

### `defect`

crear un defecto puntual

**Uso:** `olla-dft defect [-h] [-k {vacancy,substitution,interstitial}] [--site SITE] [--new-element NEW_ELEMENT] [--supercell SUPERCELL] [--position POSITION] [-o OUTDIR] file`

**Argumentos:**

- `file` — estructura de entrada (CIF, POSCAR, input de pw.x...)

**Opciones:**

| Opción | Descripción |
|---|---|
| `-k, --kind {vacancy,substitution,interstitial}` | tipo de defecto: vacancy, substitution o interstitial (default: vacancy) |
| `--site SITE` | índice del átomo afectado (base 0) |
| `--new-element` | especie que entra |
| `--supercell` | por ejemplo 3x3x3 |
| `--position` | x,y,z fraccionarias (intersticial) |
| `-o, --outdir` | carpeta de salida (default: `defecto`) |

**Fundamento físico:** [`olla-dft teoria defect`](TEORIA.md)

### `interface`

heteroestructura: apilar dos materiales con la menor deformacion de red posible

**Uso:** `olla-dft interface [-h] [-o OUTDIR] [--name NAME] [--max-index MAX_INDEX] [--tol TOL] [--max-atoms MAX_ATOMS] [--index INDEX] [--top TOP] [--list] [--separation SEPARATION] [--vacuum VACUUM] [--strain {first,second,both}] [--shift SHIFT] file1 file2`

**Argumentos:**

- `file1` — material de abajo (el sustrato)
- `file2` — material de arriba

**Opciones:**

| Opción | Descripción |
|---|---|
| `-o, --outdir` | carpeta de salida (default: `.`) |
| `--name` | nombre base de los archivos de la heteroestructura (default: heteroestructura) |
| `--max-index MAX_INDEX` | mayor coeficiente entero de la supercelda; subirlo encuentra celdas mas giradas pero tarda mucho mas (default: `4`) |
| `--tol TOL` | deformacion maxima aceptada (0.05 = 5 %%) (default: `0.05`) |
| `--max-atoms MAX_ATOMS` | máximo de átomos admitido en la supercelda de la interfaz (default: 200) |
| `--index INDEX` | cual de las candidatas construir |
| `--top TOP` | cuántas candidatas listar, de menor a mayor deformación (default: 10) |
| `--list` | solo listar las candidatas, sin construir nada |
| `--separation SEPARATION` | distancia inicial entre capas en A; por omision, de los radios de van der Waals |
| `--vacuum VACUUM` | vacío sobre la heteroestructura en Å (default: 20) |
| `--strain {first,second,both}` | quien se deforma: el de abajo, el de arriba, o los dos a medias (default: `second`) |
| `--shift` | desplazamiento lateral del material de arriba, en fracciones de la celda comun |

**Fundamento físico:** [`olla-dft teoria interface`](TEORIA.md)

### `adsorb`

sitios de adsorción sobre una losa y su energía

**Uso:** `olla-dft adsorb [-h] [-o OUTDIR] [--run] [--collect] [--pw-cmd PW_CMD] [--nproc NPROC] [-j N] [--redo] [--max-time T] [--estimate] [--timeout TIMEOUT] [--pseudo-dir PSEUDO_DIR] [--pseudo EL=UPF] [--ecutwfc ECUTWFC] [--ecutrho ECUTRHO] [--kspacing KSPACING] [--insulator] [--dpi DPI] [--format FORMAT] [--no-plot] [-t TEMPLATE] [--size {paper,poster,presentation}] [--font {sans,serif,latex}] [--usetex] [--palette PALETTE] [--background BACKGROUND] [--journal {acs,aps,elsevier,generic,iop,nature,rsc,wiley}] [--width WIDTH] [--aspect ASPECT] [--mono] --mol MOLECULA [--sites SITES] [--height HEIGHT] [--face {top,bottom}] [--rotations ROTATIONS] [--anchor ANCHOR] [--fixed-ions] [--vdw {grimme-d2,grimme-d3,DFT-D,ts-vdw,xdm,mbd}] [--dipole] [--nspin {1,2}] [--mag MAG] file`

**Argumentos:**

- `file` — estructura (CIF, POSCAR, input de pw.x...)

**Opciones:**

| Opción | Descripción |
|---|---|
| `-o, --outdir` | carpeta del barrido (default: `adsorb`) |
| `--run` | ejecutar los cálculos ahora, uno tras otro |
| `--collect` | solo analizar cálculos ya corridos |
| `--pw-cmd` | ejecutable de pw.x (anula la configuración) |
| `--nproc NPROC` | procesos MPI por cálculo |
| `-j, --jobs N` | cálculos simultáneos (default: 1). Sin --nproc, los hilos de la máquina se reparten entre ellos |
| `--redo` | rehacer también los cálculos que ya estaban terminados |
| `--max-time T` | presupuesto TOTAL de tiempo: 90m, 2h, 3600. Al agotarse no se lanzan más y el barrido queda reanudable |
| `--estimate` | estimar cuánto va a tardar el barrido y salir, usando el histórico de 'olla-dft db' |
| `--timeout TIMEOUT` | límite en segundos por cálculo |
| `--pseudo-dir` | carpeta de pseudopotenciales |
| `--pseudo EL=UPF` | forzar un pseudopotencial concreto, por ejemplo Fe=Fe.rel-pbe.UPF. Se puede repetir. Sin esto, Olla-DFT elige con 'olla-dft pseudos' |
| `--ecutwfc ECUTWFC` | cutoff de ondas (Ry) |
| `--ecutrho ECUTRHO` | cutoff de densidad (Ry) |
| `--kspacing KSPACING` | espaciado k en Å^-1 |
| `--insulator` | occupations='fixed' |
| `--dpi DPI` | resolución de las figuras de mapa de bits en puntos por pulgada (default: 600) |
| `--format` | formatos de figura separados por coma: pdf,png,svg,eps,tif (default: pdf,png) |
| `--no-plot` | solo exportar los datos, sin generar la figura |
| `-t, --template` | plantilla visual de la figura |
| `--size {paper,poster,presentation}` | tamaño de figura: paper, presentation o poster |
| `--font {sans,serif,latex}` | familia tipográfica: sans, serif o latex (Computer Modern) |
| `--usetex` | componer los textos de la figura con LaTeX real (necesita una instalación de LaTeX) |
| `--palette` | paleta de colores: grayscale, okabe-ito, okabe-ito-dark o colores hex separados por coma |
| `--background` | color de fondo de la figura, por ejemplo '#FFFFFF' o 'none' |
| `--journal {acs,aps,elsevier,generic,iop,nature,rsc,wiley}` | ancho de figura según la revista (default: generic) |
| `--width` | ancho de la figura: single, onehalf o double, o un número en milímetros |
| `--aspect ASPECT` | relación alto/ancho de la figura |
| `--mono` | versión en escala de grises: tinta negra y patrones de línea |
| `--mol MOLECULA` | adsorbato: nombre de la base de ASE (CO2, H2O, CO, NH3, O2...) o un archivo con la molécula |
| `--sites` | tipos de sitio a probar (default: los tres) |
| `--height HEIGHT` | altura inicial del adsorbato sobre el sitio, en Å (default: 2.0) |
| `--face {top,bottom}` | cara de la losa donde adsorber (default: `top`) |
| `--rotations ROTATIONS` | orientaciones a probar girando alrededor de la normal (default: 1) |
| `--anchor ANCHOR` | átomo de la molécula que se apoya en el sitio (índice desde 0; default: 0) |
| `--fixed-ions` | no relajar: solo scf en la geometría inicial |
| `--vdw {grimme-d2,grimme-d3,DFT-D,ts-vdw,xdm,mbd}` | corrección de dispersión (casi obligatoria en fisisorción) |
| `--dipole` | corrección dipolar: la losa con adsorbato en una sola cara es polar |
| `--nspin {1,2}` | 2 activa la polarización de espín (default: 1) |
| `--mag` | magnetización inicial (implica --nspin 2) |

**Fundamento físico:** [`olla-dft teoria adsorb`](TEORIA.md)

### `eform`

energía de formación de defectos cargados, niveles de transición y diagrama E_f vs ε_F

**Uso:** `olla-dft eform [-h] [-o OUTDIR] [--run] [--collect] [--pw-cmd PW_CMD] [--nproc NPROC] [-j N] [--redo] [--max-time T] [--estimate] [--timeout TIMEOUT] [--pseudo-dir PSEUDO_DIR] [--pseudo EL=UPF] [--ecutwfc ECUTWFC] [--ecutrho ECUTRHO] [--kspacing KSPACING] [--insulator] [--dpi DPI] [--format FORMAT] [--no-plot] [-t TEMPLATE] [--size {paper,poster,presentation}] [--font {sans,serif,latex}] [--usetex] [--palette PALETTE] [--background BACKGROUND] [--journal {acs,aps,elsevier,generic,iop,nature,rsc,wiley}] [--width WIDTH] [--aspect ASPECT] [--mono] [-k {vacancy,substitution,interstitial}] [--site SITE] [--new-element NEW_ELEMENT] [--position POSITION] [--supercell SUPERCELL] [-q CHARGES] [--epsilon EPSILON] [--correction {ninguna,makov-payne,lany-zunger}] [--mu EL=eV] [--align POT_DEF POT_PERF] [--dv DV] [--fixed-ions] [--vdw {grimme-d2,grimme-d3,DFT-D,ts-vdw,xdm,mbd}] [--nspin {1,2}] [--mag MAG] file`

**Argumentos:**

- `file` — estructura (CIF, POSCAR, input de pw.x...)

**Opciones:**

| Opción | Descripción |
|---|---|
| `-o, --outdir` | carpeta del barrido (default: `formacion`) |
| `--run` | ejecutar los cálculos ahora, uno tras otro |
| `--collect` | solo analizar cálculos ya corridos |
| `--pw-cmd` | ejecutable de pw.x (anula la configuración) |
| `--nproc NPROC` | procesos MPI por cálculo |
| `-j, --jobs N` | cálculos simultáneos (default: 1). Sin --nproc, los hilos de la máquina se reparten entre ellos |
| `--redo` | rehacer también los cálculos que ya estaban terminados |
| `--max-time T` | presupuesto TOTAL de tiempo: 90m, 2h, 3600. Al agotarse no se lanzan más y el barrido queda reanudable |
| `--estimate` | estimar cuánto va a tardar el barrido y salir, usando el histórico de 'olla-dft db' |
| `--timeout TIMEOUT` | límite en segundos por cálculo |
| `--pseudo-dir` | carpeta de pseudopotenciales |
| `--pseudo EL=UPF` | forzar un pseudopotencial concreto, por ejemplo Fe=Fe.rel-pbe.UPF. Se puede repetir. Sin esto, Olla-DFT elige con 'olla-dft pseudos' |
| `--ecutwfc ECUTWFC` | cutoff de ondas (Ry) |
| `--ecutrho ECUTRHO` | cutoff de densidad (Ry) |
| `--kspacing KSPACING` | espaciado k en Å^-1 |
| `--insulator` | occupations='fixed' |
| `--dpi DPI` | resolución de las figuras de mapa de bits en puntos por pulgada (default: 600) |
| `--format` | formatos de figura separados por coma: pdf,png,svg,eps,tif (default: pdf,png) |
| `--no-plot` | solo exportar los datos, sin generar la figura |
| `-t, --template` | plantilla visual de la figura |
| `--size {paper,poster,presentation}` | tamaño de figura: paper, presentation o poster |
| `--font {sans,serif,latex}` | familia tipográfica: sans, serif o latex (Computer Modern) |
| `--usetex` | componer los textos de la figura con LaTeX real (necesita una instalación de LaTeX) |
| `--palette` | paleta de colores: grayscale, okabe-ito, okabe-ito-dark o colores hex separados por coma |
| `--background` | color de fondo de la figura, por ejemplo '#FFFFFF' o 'none' |
| `--journal {acs,aps,elsevier,generic,iop,nature,rsc,wiley}` | ancho de figura según la revista (default: generic) |
| `--width` | ancho de la figura: single, onehalf o double, o un número en milímetros |
| `--aspect ASPECT` | relación alto/ancho de la figura |
| `--mono` | versión en escala de grises: tinta negra y patrones de línea |
| `-k, --kind {vacancy,substitution,interstitial}` | tipo de defecto (default: `vacancy`) |
| `--site SITE` | índice del átomo afectado en la supercelda (base 0) |
| `--new-element` | especie que entra |
| `--position` | x,y,z fraccionarias (intersticial) |
| `--supercell` | tamaño de la supercelda (default: 2x2x2) |
| `-q, --charges` | estados de carga separados por coma, por ejemplo -2,-1,0,1,2 (default: `0`) |
| `--epsilon EPSILON` | constante dieléctrica del material, para apantallar la corrección de imagen |
| `--correction {ninguna,makov-payne,lany-zunger}` | esquema de corrección de tamaño finito (default: `lany-zunger`) |
| `--mu EL=eV` | potencial químico por elemento, en eV por átomo. Se puede repetir |
| `--align ('POT_DEF', 'POT_PERF')` | dos archivos cube de potencial electrostático (defecto y perfecto) para el término ΔV |
| `--dv DV` | alineamiento ΔV en eV, si ya lo tienes calculado |
| `--fixed-ions` | no relajar el defecto en cada estado de carga |
| `--vdw {grimme-d2,grimme-d3,DFT-D,ts-vdw,xdm,mbd}` | corrección de dispersión (van der Waals): grimme-d2, grimme-d3, DFT-D, ts-vdw, xdm o mbd |
| `--nspin {1,2}` | 2 activa la polarización de espín (default: 1) |
| `--mag` | magnetización inicial (implica --nspin 2) |

**Fundamento físico:** [`olla-dft teoria eform`](TEORIA.md)

### `align`

alineamiento de bandas entre dos materiales: offsets ΔE_v, ΔE_c y tipo I/II/III

**Uso:** `olla-dft align [-h] [--interface CARPETA] [--names NAMES] [--axis AXIS] [--window A] [--rerun] [-o OUTDIR] [--pw-cmd PW_CMD] [--nproc NPROC] [--dpi DPI] [--format FORMAT] [--no-plot] [-t TEMPLATE] [--font {sans,serif,latex}] [--usetex] [--palette PALETTE] [--background BACKGROUND] [--journal {acs,aps,elsevier,generic,iop,nature,rsc,wiley}] [--width WIDTH] [--mono] a b`

**Argumentos:**

- `a` — carpeta del cálculo del primer material
- `b` — carpeta del cálculo del segundo material

**Opciones:**

| Opción | Descripción |
|---|---|
| `--interface CARPETA` | carpeta de la interfaz; activa el método riguroso de Van de Walle-Martin |
| `--names` | nombres para el reporte, separados por coma (por omisión, los de las carpetas) |
| `--axis` | eje del perfil planar (default: `c`) |
| `--window A` | ventana del promedio macroscópico en Å (por omisión, un octavo de la celda) |
| `--rerun` | volver a correr pp.x aunque ya exista el cube |
| `-o, --outdir` | carpeta de salida (default: `alineamiento`) |
| `--pw-cmd` | ejecutable de pw.x para --run; de su ruta salen los demás binarios de QE |
| `--nproc NPROC` | número de procesos MPI para los cálculos que se lanzan con --run |
| `--dpi DPI` | resolución de las figuras de mapa de bits en puntos por pulgada (default: 600) |
| `--format` | formatos de figura separados por coma: pdf,png,svg,eps,tif (default: pdf,png) |
| `--no-plot` | solo exportar los datos, sin generar la figura |
| `-t, --template` | plantilla visual de la figura (la lista sale con 'olla-dft templates list') |
| `--font {sans,serif,latex}` | familia tipográfica: sans, serif o latex (Computer Modern) |
| `--usetex` | componer los textos de la figura con LaTeX real (necesita una instalación de LaTeX) |
| `--palette` | paleta de colores: grayscale, okabe-ito, okabe-ito-dark o colores hex separados por coma |
| `--background` | color de fondo de la figura, por ejemplo '#FFFFFF' o 'none' |
| `--journal {acs,aps,elsevier,generic,iop,nature,rsc,wiley}` | ancho de figura según la revista (default: generic) |
| `--width` | ancho de la figura: single, onehalf o double, o un número en milímetros |
| `--mono` | versión en escala de grises: tinta negra y patrones de línea |

**Fundamento físico:** [`olla-dft teoria align`](TEORIA.md)

### `esm`

superficies cargadas con medio de apantallamiento efectivo: función trabajo, capacitancia y potencial de carga cero

**Uso:** `olla-dft esm [-h] [-o OUTDIR] [--run] [--collect] [--pw-cmd PW_CMD] [--nproc NPROC] [-j N] [--redo] [--max-time T] [--estimate] [--timeout TIMEOUT] [--pseudo-dir PSEUDO_DIR] [--pseudo EL=UPF] [--ecutwfc ECUTWFC] [--ecutrho ECUTRHO] [--kspacing KSPACING] [--insulator] [--dpi DPI] [--format FORMAT] [--no-plot] [-t TEMPLATE] [--size {paper,poster,presentation}] [--font {sans,serif,latex}] [--usetex] [--palette PALETTE] [--background BACKGROUND] [--journal {acs,aps,elsevier,generic,iop,nature,rsc,wiley}] [--width WIDTH] [--aspect ASPECT] [--mono] [--bc {bc1,bc2,bc3}] [--charge CHARGE] [--field FIELD] [--esm-w WIDTH_ESM] [--nfit NFIT] file`

**Argumentos:**

- `file` — estructura (CIF, POSCAR, input de pw.x...)

**Opciones:**

| Opción | Descripción |
|---|---|
| `-o, --outdir` | carpeta del barrido (default: `esm`) |
| `--run` | ejecutar los cálculos ahora, uno tras otro |
| `--collect` | solo analizar cálculos ya corridos |
| `--pw-cmd` | ejecutable de pw.x (anula la configuración) |
| `--nproc NPROC` | procesos MPI por cálculo |
| `-j, --jobs N` | cálculos simultáneos (default: 1). Sin --nproc, los hilos de la máquina se reparten entre ellos |
| `--redo` | rehacer también los cálculos que ya estaban terminados |
| `--max-time T` | presupuesto TOTAL de tiempo: 90m, 2h, 3600. Al agotarse no se lanzan más y el barrido queda reanudable |
| `--estimate` | estimar cuánto va a tardar el barrido y salir, usando el histórico de 'olla-dft db' |
| `--timeout TIMEOUT` | límite en segundos por cálculo |
| `--pseudo-dir` | carpeta de pseudopotenciales |
| `--pseudo EL=UPF` | forzar un pseudopotencial concreto, por ejemplo Fe=Fe.rel-pbe.UPF. Se puede repetir. Sin esto, Olla-DFT elige con 'olla-dft pseudos' |
| `--ecutwfc ECUTWFC` | cutoff de ondas (Ry) |
| `--ecutrho ECUTRHO` | cutoff de densidad (Ry) |
| `--kspacing KSPACING` | espaciado k en Å^-1 |
| `--insulator` | occupations='fixed' |
| `--dpi DPI` | resolución de las figuras de mapa de bits en puntos por pulgada (default: 600) |
| `--format` | formatos de figura separados por coma: pdf,png,svg,eps,tif (default: pdf,png) |
| `--no-plot` | solo exportar los datos, sin generar la figura |
| `-t, --template` | plantilla visual de la figura |
| `--size {paper,poster,presentation}` | tamaño de figura: paper, presentation o poster |
| `--font {sans,serif,latex}` | familia tipográfica: sans, serif o latex (Computer Modern) |
| `--usetex` | componer los textos de la figura con LaTeX real (necesita una instalación de LaTeX) |
| `--palette` | paleta de colores: grayscale, okabe-ito, okabe-ito-dark o colores hex separados por coma |
| `--background` | color de fondo de la figura, por ejemplo '#FFFFFF' o 'none' |
| `--journal {acs,aps,elsevier,generic,iop,nature,rsc,wiley}` | ancho de figura según la revista (default: generic) |
| `--width` | ancho de la figura: single, onehalf o double, o un número en milímetros |
| `--aspect ASPECT` | relación alto/ancho de la figura |
| `--mono` | versión en escala de grises: tinta negra y patrones de línea |
| `--bc {bc1,bc2,bc3}` | bc1 vacío/vacío (losas neutras), bc2 metal/metal (condensador), bc3 vacío/metal (electrodo, el único junto con bc2 que admite carga neta) (default: `bc1`) |
| `--charge` | cargas netas en e, separadas por coma: -0.2,0,0.2 (default: `0`) |
| `--field FIELD` | campo aplicado en Ry/u.a. (solo con bc2) |
| `--esm-w WIDTH_ESM` | desplazamiento de la frontera de ESM en u.a. |
| `--nfit NFIT` | puntos de ajuste del potencial en la frontera (default: `4`) |

**Fundamento físico:** [`olla-dft teoria esm`](TEORIA.md)

### `echem`

electrodo de hidrógeno computacional: HER, OER, potencial limitante y sobrepotencial

**Uso:** `olla-dft echem [-h] [--her E_ads] [--oer OH=..,O=..,OOH=..] [--corrections X=eV] [-U POTENTIAL] [--ph PH] [-T TEMPERATURE] [-o OUTDIR] [--dpi DPI] [--format FORMAT] [--no-plot] [-t TEMPLATE] [--font {sans,serif,latex}] [--usetex] [--palette PALETTE] [--background BACKGROUND] [--journal {acs,aps,elsevier,generic,iop,nature,rsc,wiley}] [--width WIDTH] [--mono]`

**Opciones:**

| Opción | Descripción |
|---|---|
| `--her E_ads` | energía de adsorción de H en eV (reacción HER) |
| `--oer OH=..,O=..,OOH=..` | energías de adsorción de los tres intermedios de la OER, en eV y referidas al agua |
| `--corrections X=eV` | correcciones térmicas ZPE−TΔS por intermedio; sin esto se usan las estándar de la literatura |
| `-U, --potential POTENTIAL` | potencial aplicado en V frente al SHE (a pH 0 es el mismo que frente al RHE; el pH lo convierte) |
| `--ph PH` | pH |
| `-T, --temperature TEMPERATURE` | temperatura en K (default: 298.15) |
| `-o, --outdir` | carpeta de salida (default: `echem`) |
| `--dpi DPI` | resolución de las figuras de mapa de bits en puntos por pulgada (default: 600) |
| `--format` | formatos de figura separados por coma: pdf,png,svg,eps,tif (default: pdf,png) |
| `--no-plot` | solo exportar los datos, sin generar la figura |
| `-t, --template` | plantilla visual de la figura (la lista sale con 'olla-dft templates list') |
| `--font {sans,serif,latex}` | familia tipográfica: sans, serif o latex (Computer Modern) |
| `--usetex` | componer los textos de la figura con LaTeX real (necesita una instalación de LaTeX) |
| `--palette` | paleta de colores: grayscale, okabe-ito, okabe-ito-dark o colores hex separados por coma |
| `--background` | color de fondo de la figura, por ejemplo '#FFFFFF' o 'none' |
| `--journal {acs,aps,elsevier,generic,iop,nature,rsc,wiley}` | ancho de figura según la revista (default: generic) |
| `--width` | ancho de la figura: single, onehalf o double, o un número en milímetros |
| `--mono` | versión en escala de grises: tinta negra y patrones de línea |

**Fundamento físico:** [`olla-dft teoria echem`](TEORIA.md)

### `neb`

camino de reaccion y barrera de activacion (neb.x)

**Uso:** `olla-dft neb [-h] [-o OUTDIR] [--images IMAGES] [--no-ci] [--path-thr PATH_THR] [--nstep NSTEP] [--fix FIX] [--prefix PREFIX] [--collect] [--pseudo-dir PSEUDO_DIR] [--pseudo EL=UPF] [--ecutwfc ECUTWFC] [--ecutrho ECUTRHO] [--kspacing KSPACING] [--metal] [--nspin {1,2}] [--mag MAG] [--dpi DPI] [--format FORMAT] [--no-plot] [-t TEMPLATE] [--font {sans,serif,latex}] [--usetex] [--palette PALETTE] [--background BACKGROUND] [--journal {acs,aps,elsevier,generic,iop,nature,rsc,wiley}] [--width WIDTH] [--mono] file [final]`

**Argumentos:**

- `file` — estructura inicial (reactivo)
- `final` — estructura final (producto)

**Opciones:**

| Opción | Descripción |
|---|---|
| `-o, --outdir` | carpeta de salida (default: `neb`) |
| `--images IMAGES` | numero de imagenes de la cadena (default: `7`) |
| `--no-ci` | sin imagen trepadora; la barrera saldra SUBESTIMADA |
| `--path-thr PATH_THR` | umbral de fuerza del camino en eV/A (default: `0.05`) |
| `--nstep NSTEP` | pasos máximos de optimización del camino en neb.x (default: 50) |
| `--fix` | indices de atomos a congelar (base 0) |
| `--prefix` | prefix del cálculo (se detecta solo) |
| `--collect` | leer los resultados de un cálculo ya corrido en lugar de preparar los inputs |
| `--pseudo-dir` | carpeta con los pseudopotenciales UPF (si no se da, la de 'olla-dft config') |
| `--pseudo EL=UPF` | forzar un pseudopotencial concreto, por ejemplo Fe=Fe.rel-pbe.UPF. Se puede repetir. Sin esto, Olla-DFT elige con 'olla-dft pseudos' |
| `--ecutwfc ECUTWFC` | cutoff de funciones de onda en Ry (si no se da, el recomendado por los UPF) |
| `--ecutrho ECUTRHO` | cutoff de densidad en Ry (si no se da, el recomendado por los UPF) |
| `--kspacing KSPACING` | espaciado de la malla k en Å^-1 |
| `--metal` | sistema metálico: ocupaciones con smearing en vez de fijas |
| `--nspin {1,2}` | 2 activa la polarización de espín (default: 1) |
| `--mag` | magnetización inicial: un número (0.5) o por elemento (Fe=0.7,O=0). Implica --nspin 2 |
| `--dpi DPI` | resolución de las figuras de mapa de bits en puntos por pulgada (default: 600) |
| `--format` | formatos de figura separados por coma: pdf,png,svg,eps,tif (default: pdf,png) |
| `--no-plot` | solo exportar los datos, sin generar la figura |
| `-t, --template` | plantilla visual de la figura (la lista sale con 'olla-dft templates list') |
| `--font {sans,serif,latex}` | familia tipográfica: sans, serif o latex (Computer Modern) |
| `--usetex` | componer los textos de la figura con LaTeX real (necesita una instalación de LaTeX) |
| `--palette` | paleta de colores: grayscale, okabe-ito, okabe-ito-dark o colores hex separados por coma |
| `--background` | color de fondo de la figura, por ejemplo '#FFFFFF' o 'none' |
| `--journal {acs,aps,elsevier,generic,iop,nature,rsc,wiley}` | ancho de figura según la revista (default: generic) |
| `--width` | ancho de la figura: single, onehalf o double, o un número en milímetros |
| `--mono` | versión en escala de grises: tinta negra y patrones de línea |

**Fundamento físico:** [`olla-dft teoria neb`](TEORIA.md)

### `amorphous`

sólido amorfo por fundido y temple con un potencial aprendido

**Uso:** `olla-dft amorphous [-h] [-n UNITS] -d G_CM3 [--melt K] [--final K] [--melt-steps MELT_STEPS] [--quench-steps QUENCH_STEPS] [--anneal-steps ANNEAL_STEPS] [--dt FS] [--model MODEL] [--min-dist F] [--seed SEED] [--pack-only] [-o OUTDIR] formula`

**Argumentos:**

- `formula` — fórmula de la unidad, por ejemplo SiO2

**Opciones:**

| Opción | Descripción |
|---|---|
| `-n, --units UNITS` | unidades de fórmula en la celda (default: 8) |
| `-d, --density G_CM3` | densidad objetivo en g/cm³ |
| `--melt K` | temperatura de fundido (default: 3000 K) |
| `--final K` | temperatura final (default: 300 K) |
| `--melt-steps MELT_STEPS` | pasos de dinámica en la fase de fundido (default: 500) |
| `--quench-steps QUENCH_STEPS` | pasos del temple: son los que fijan la velocidad. El default (1000) es un temple de exploración a ~3e15 K/s y el reporte lo avisa; 27000 baja a 1e14 K/s |
| `--anneal-steps ANNEAL_STEPS` | pasos de dinámica del recocido a la temperatura final (default: 200) |
| `--dt FS` | paso: de integración en fs (amorphous) o de la malla de temperaturas en K (qha) (default: `1.0`) |
| `--model` | potencial interatómico (default: `mace`) |
| `--min-dist F` | factor sobre la suma de radios covalentes al empaquetar (default: 0.75) |
| `--seed SEED` | semilla; cambia para generar otra realización |
| `--pack-only` | solo empaquetar, sin dinámica |
| `-o, --outdir` | carpeta de salida (default: `amorfo`) |

**Fundamento físico:** [`olla-dft teoria amorphous`](TEORIA.md)

## Automatización y calidad

### `doctor`

diagnosticar un cálculo: convergencia, fuerzas y por qué no converge

**Uso:** `olla-dft doctor [-h] [--system] [--project PROJECT] [--json] [--prefix PREFIX] [-o OUTDIR] [--no-plot] [--dpi DPI] [--format FORMAT] [-t TEMPLATE] [--font {sans,serif,latex}] [--usetex] [--palette PALETTE] [--background BACKGROUND] [--journal {acs,aps,elsevier,generic,iop,nature,rsc,wiley}] [--width WIDTH] [--mono] [path]`

**Argumentos:**

- `path` — carpeta del cálculo o archivo de salida

**Opciones:**

| Opción | Descripción |
|---|---|
| `--system` | revisar instalación, recursos, QE y pseudopotenciales |
| `--project` | además, revisar la puerta de calidad de este proyecto |
| `--json` | imprimir el diagnóstico como JSON |
| `--prefix` | prefix del cálculo (se detecta solo) |
| `-o, --outdir` | carpeta de salida (default: `.`) |
| `--no-plot` | solo exportar los datos, sin generar la figura |
| `--dpi DPI` | resolución de las figuras de mapa de bits en puntos por pulgada (default: 600) |
| `--format` | formatos de figura separados por coma: pdf,png,svg,eps,tif (default: pdf,png) |
| `-t, --template` | plantilla visual de la figura (la lista sale con 'olla-dft templates list') |
| `--font {sans,serif,latex}` | familia tipográfica: sans, serif o latex (Computer Modern) |
| `--usetex` | componer los textos de la figura con LaTeX real (necesita una instalación de LaTeX) |
| `--palette` | paleta de colores: grayscale, okabe-ito, okabe-ito-dark o colores hex separados por coma |
| `--background` | color de fondo de la figura, por ejemplo '#FFFFFF' o 'none' |
| `--journal {acs,aps,elsevier,generic,iop,nature,rsc,wiley}` | ancho de figura según la revista (default: generic) |
| `--width` | ancho de la figura: single, onehalf o double, o un número en milímetros |
| `--mono` | versión en escala de grises: tinta negra y patrones de línea |

**Fundamento físico:** [`olla-dft teoria doctor`](TEORIA.md)

### `audit`

verificar que un conjunto de cálculos sea comparable antes de restar energías

**Uso:** `olla-dft audit [-h] [--index] [--db DB] paths [paths ...]`

**Argumentos:**

- `paths` — carpetas o archivos XML de los cálculos

**Opciones:**

| Opción | Descripción |
|---|---|
| `--index` | además, registrarlos en la base de datos |
| `--db` | archivo SQLite del índice de cálculos (default: olla-dft.db) |

**Fundamento físico:** [`olla-dft teoria audit`](TEORIA.md)

### `crosscheck`

cruzar la misma cantidad por rutas independientes

**Uso:** `olla-dft crosscheck [-h] [-f FILE] [--gap-bandas GAP_BANDAS] [--gap-tauc GAP_TAUC] [project]`

**Argumentos:**

- `project` — carpeta del proyecto

**Opciones:**

| Opción | Descripción |
|---|---|
| `-f, --file` | estructura (para masas y volumen) |
| `--gap-bandas GAP_BANDAS` | gap de la estructura de bandas en eV, para cruzarlo con el de Tauc |
| `--gap-tauc GAP_TAUC` | gap de la extrapolación de Tauc en eV, para cruzarlo con el de bandas |

**Fundamento físico:** [`olla-dft teoria crosscheck`](TEORIA.md)

### `cost`

qué sabe Olla-DFT de la velocidad de tu máquina

**Uso:** `olla-dft cost [-h] [--db DB]`

**Opciones:**

| Opción | Descripción |
|---|---|
| `--db` | base de cálculos (default: `olla-dft.db`) |

**Fundamento físico:** [`olla-dft teoria cost`](TEORIA.md)

### `db`

índice local de cálculos

**Uso:** `olla-dft db [-h] [--db DB] [-q QUERY] [--export EXPORT] [--formula FORMULA] [--calculation CALCULATION] [--gap-min GAP_MIN] [--gap-max GAP_MAX] [--limit LIMIT] [paths ...]`

**Argumentos:**

- `paths` — carpetas a registrar

**Opciones:**

| Opción | Descripción |
|---|---|
| `--db` | archivo SQLite del índice de cálculos (default: olla-dft.db) |
| `-q, --query` | consulta SQL (solo SELECT) |
| `--export` | exportar todo a un JSON |
| `--formula` | filtrar por fórmula, por ejemplo Si |
| `--calculation` | filtrar por tipo: scf, relax, nscf... |
| `--gap-min GAP_MIN` | gap mínimo en eV |
| `--gap-max GAP_MAX` | gap máximo en eV |
| `--limit LIMIT` | máximo de filas filtradas (default: `100`) |

**Fundamento físico:** [`olla-dft teoria db`](TEORIA.md)

### `hull`

energías de formación y casco convexo

**Uso:** `olla-dft hull [-h] [-o OUTDIR] [--elements ELEMENTS] [--threshold THRESHOLD] [--force] [--dpi DPI] [--format FORMAT] [--no-plot] [-t TEMPLATE] [--font {sans,serif,latex}] [--usetex] [--palette PALETTE] [--background BACKGROUND] [--journal {acs,aps,elsevier,generic,iop,nature,rsc,wiley}] [--width WIDTH] [--mono] paths [paths ...]`

**Argumentos:**

- `paths` — carpetas o archivos XML de los cálculos

**Opciones:**

| Opción | Descripción |
|---|---|
| `-o, --outdir` | carpeta de salida (default: `.`) |
| `--elements` | orden de los elementos, por ejemplo Zn,Al |
| `--threshold THRESHOLD` | umbral de metaestabilidad en eV/átomo (default: `0.025`) |
| `--force` | construir el casco aunque la auditoría falle |
| `--dpi DPI` | resolución de las figuras de mapa de bits en puntos por pulgada (default: 600) |
| `--format` | formatos de figura separados por coma: pdf,png,svg,eps,tif (default: pdf,png) |
| `--no-plot` | solo exportar los datos, sin generar la figura |
| `-t, --template` | plantilla visual de la figura (la lista sale con 'olla-dft templates list') |
| `--font {sans,serif,latex}` | familia tipográfica: sans, serif o latex (Computer Modern) |
| `--usetex` | componer los textos de la figura con LaTeX real (necesita una instalación de LaTeX) |
| `--palette` | paleta de colores: grayscale, okabe-ito, okabe-ito-dark o colores hex separados por coma |
| `--background` | color de fondo de la figura, por ejemplo '#FFFFFF' o 'none' |
| `--journal {acs,aps,elsevier,generic,iop,nature,rsc,wiley}` | ancho de figura según la revista (default: generic) |
| `--width` | ancho de la figura: single, onehalf o double, o un número en milímetros |
| `--mono` | versión en escala de grises: tinta negra y patrones de línea |

**Fundamento físico:** [`olla-dft teoria hull`](TEORIA.md)

### `mlip`

potencial aprendido: pre-relajar y cribar antes de gastar DFT

**Uso:** `olla-dft mlip [-h] [-o OUTPUT] [--model {mace,chgnet,m3gnet}] [--size SIZE] [--device DEVICE] [--fmax FMAX] [--steps STEPS] [--fixed-cell] [--span SPAN] [--npoints NPOINTS] [--supercell SUPERCELL] {relax,scan,phonons} file`

**Argumentos:**

- `action` {relax,scan,phonons} — acción a realizar (ver la lista de arriba)
- `file` — estructura de entrada (CIF, POSCAR, input de pw.x...)

**Opciones:**

| Opción | Descripción |
|---|---|
| `-o, --output` | estructura de salida (relax) |
| `--model {mace,chgnet,m3gnet}` | potencial aprendido: mace, chgnet o m3gnet (default: mace) |
| `--size` | tamaño del modelo MACE (small/medium/large) (default: `small`) |
| `--device` | dispositivo donde corre el potencial: cpu o cuda (default: cpu) |
| `--fmax FMAX` | fuerza objetivo en eV/Å (default: `0.01`) |
| `--steps STEPS` | pasos máximos de la relajación (default: 300) |
| `--fixed-cell` | no relajar la celda, solo las posiciones |
| `--span SPAN` | rango del barrido de volumen (scan) (default: `0.1`) |
| `--npoints NPOINTS` | puntos del barrido de volumen (default: 15) |
| `--supercell` | supercelda del cribado, ej. 2x2x2 |

**Fundamento físico:** [`olla-dft teoria mlip`](TEORIA.md)

### `suggest`

sugerir parámetros a partir de tus cálculos previos

**Uso:** `olla-dft suggest [-h] [--db DB] file`

**Argumentos:**

- `file` — estructura de entrada (CIF, POSCAR, input de pw.x...)

**Opciones:**

| Opción | Descripción |
|---|---|
| `--db` | archivo SQLite del índice de cálculos (default: olla-dft.db) |

**Fundamento físico:** [`olla-dft teoria suggest`](TEORIA.md)

### `datasheet`

ficha del material y párrafo de métodos

**Uso:** `olla-dft datasheet [-h] [-o OUTDIR] [--name NAME] [--methods] [project]`

**Argumentos:**

- `project` — carpeta del proyecto (default: .)

**Opciones:**

| Opción | Descripción |
|---|---|
| `-o, --outdir` | carpeta de salida (default: `.`) |
| `--name` | nombre base de los archivos |
| `--methods` | solo el párrafo de metodología y las citas |

### `report`

registro local de fallas y confusiones

**Uso:** `olla-dft report [-h] [--show SHOW] [--close CLOSE] [--note NOTE] [--stats] [--export EXPORT] [--only-open] [--attach ATTACH] [description ...]`

**Argumentos:**

- `description` — qué pasó (si se omite, lista las incidencias)

**Opciones:**

| Opción | Descripción |
|---|---|
| `--show` | ver una incidencia por su id |
| `--close` | marcar una incidencia como resuelta |
| `--note` | nota al cerrar |
| `--stats` | qué subcomandos fallan más |
| `--export` | empaquetar todo en un archivo JSON |
| `--only-open` | listar solo las incidencias abiertas |
| `--attach` | adjuntar un archivo (se copia al registro local) |

### `compare`

comparar corridas sin restar energías incompatibles

**Uso:** `olla-dft compare [-h] [--reference REFERENCE] [-o OUTPUT] paths [paths ...]`

**Argumentos:**

- `paths` — carpetas o XML de las corridas

**Opciones:**

| Opción | Descripción |
|---|---|
| `--reference REFERENCE` | índice de la corrida de referencia (default: 0) |
| `-o, --output` | guardar comparación en JSON |

### `tune`

recomendar el siguiente punto de una convergencia

**Uso:** `olla-dft tune [-h] [--threshold THRESHOLD] [-o OUTPUT] file`

**Argumentos:**

- `file` — CONVERGENCIA.dat

**Opciones:**

| Opción | Descripción |
|---|---|
| `--threshold THRESHOLD` | umbral en meV/átomo (default: 1) |
| `-o, --output` | guardar recomendación en JSON |

**Fundamento físico:** [`olla-dft teoria tune`](TEORIA.md)

### `results`

ingerir, consultar y exportar resultados normalizados del proyecto

**Uso:** `olla-dft results [-h] [--project PROJECT] [--db DB] [--tag TAG] [--formula FORMULA] [--calculation CALCULATION] [--status {invalid,not_converged,parsed_no_energy,parsed,converged}] [--review-status {unreviewed,accepted,rejected}] [--note NOTE] [--limit LIMIT] [--json] [-o OUTPUT] {ingest,list,show,review,export,explore} [target] [extra_paths ...]`

**Argumentos:**

- `action` {ingest,list,show,review,export,explore} — acción a realizar (ver la lista de arriba)
- `target` — ruta de entrada para ingest, o id para show
- `extra_paths` — más carpetas/XML para ingest

**Opciones:**

| Opción | Descripción |
|---|---|
| `--project` | carpeta del proyecto (default: .) |
| `--db` | SQLite alternativa; por defecto .qekit/results.sqlite3 |
| `--tag` | etiqueta de procedencia para ingest |
| `--formula` | filtrar por fórmula |
| `--calculation` | filtrar por tipo de cálculo |
| `--status {invalid,not_converged,parsed_no_energy,parsed,converged}` | filtrar por estado: invalid, not_converged, parsed_no_energy, parsed o converged |
| `--review-status {unreviewed,accepted,rejected}` | en review, estado de la revisión humana |
| `--note` | en review, nota que acompaña la decisión |
| `--limit LIMIT` | máximo de registros: list=100, explore=10000 |
| `--json` | en list, imprimir JSON |
| `-o, --output` | archivo de salida: export=JSON, explore=HTML interactivo |

### `campaign`

crear matrices reproducibles de tareas parametrizadas

**Uso:** `olla-dft campaign [-h] [--project PROJECT] [--command CAMPAIGN_COMMAND] [--axis AXIS] [--goal GOAL] [--convergence-file CONVERGENCE_FILE] [--adaptive] [--threshold THRESHOLD] [--execute] [--force] [--parallel PARALLEL] [--retries RETRIES] [--timeout TIMEOUT] [--cancel-file CANCEL_FILE] [-o OUTPUT] {create,list,status,export,run,extend} [target]`

**Argumentos:**

- `action` {create,list,status,export,run,extend} — acción a realizar (ver la lista de arriba)
- `target` — nombre o id de campaña

**Opciones:**

| Opción | Descripción |
|---|---|
| `--project` | carpeta del proyecto (default: .) |
| `--command` | plantilla del comando olla-dft; campos: {eje}, {index}, {id}, {structure} |
| `--axis` | eje nombre=v1,v2; se puede repetir (default: `[]`) |
| `--goal` | objetivo científico de la campaña |
| `--convergence-file` | CONVERGENCIA.dat para tomar una recomendación |
| `--adaptive` | añadir el siguiente valor recomendado al eje de convergencia |
| `--threshold THRESHOLD` | umbral de convergencia al extender (meV/átomo) |
| `--execute` | en run, ejecutar los puntos seleccionados |
| `--force` | en run, ignorar la caché de tareas |
| `--parallel PARALLEL` | en run, puntos independientes simultáneos (default: 1) |
| `--retries RETRIES` | reintentos por punto fallido (default: 0) |
| `--timeout TIMEOUT` | tiempo máximo por intento, en segundos |
| `--cancel-file` | marker de cancelación cooperativa personalizado |
| `-o, --output` | archivo JSON de export |

### `pseudos`

comparar los pseudopotenciales disponibles y elegir con criterio, no por orden alfabetico

**Uso:** `olla-dft pseudos [-h] [--element ELEMENT] [--task TASK] [--functional FUNCTIONAL] [--cheap] [--pseudo-dir PSEUDO_DIR] [--pseudo EL=UPF] [file]`

**Argumentos:**

- `file` — estructura

**Opciones:**

| Opción | Descripción |
|---|---|
| `--element` | elementos separados por coma |
| `--task` | para que es: general, optics, soc, xanes, hubbard, fonones. Cada tarea descarta los que no sirven (default: `general`) |
| `--functional` | exigir un funcional concreto (PBE, PZ, PBEsol...) |
| `--cheap` | preferir ultrasuave/PAW, que necesitan menos ondas planas |
| `--pseudo-dir` | carpeta con los pseudopotenciales UPF (si no se da, la de 'olla-dft config') |
| `--pseudo EL=UPF` | forzar un pseudopotencial concreto, por ejemplo Fe=Fe.rel-pbe.UPF. Se puede repetir. Sin esto, Olla-DFT elige con 'olla-dft pseudos' |

**Fundamento físico:** [`olla-dft teoria pseudos`](TEORIA.md)

## Proyecto

### `project`

gestionar un proyecto reproducible: fuentes, workflow, calidad y dashboard

**Uso:** `olla-dft project [-h] [--project PROJECT] [--name NAME] [--command TASK_COMMANDS] [--execute] [--force] [--parallel PARALLEL] [--retries RETRIES] [--timeout TIMEOUT] [--cancel-file CANCEL_FILE] [--reason REASON] [--selftest] [--advanced] [-o OUTPUT] [--pdf] [--theme {auto,light,dark}] [--language {es,en}] [--both] [--verify-environment] [--other OTHER] [--json] {init,add,plan,show,status,validate,run,dashboard,report,export,ingest,environment,diff,cancel,resume} [target]`

**Argumentos:**

- `action` {init,add,plan,show,status,validate,run,dashboard,report,export,ingest,environment,diff,cancel,resume} — acción sobre el proyecto
- `target` — directorio, archivo, objetivo, perfil o tarea según acción

**Opciones:**

| Opción | Descripción |
|---|---|
| `--project` | proyecto desde el que trabajar (default: .) |
| `--name` | nombre al inicializar |
| `--command` | tarea olla-dft personalizada; se puede repetir con plan |
| `--execute` | ejecutar run/submit; por omisión solo simula o escribe |
| `--force` | en run, ignorar la caché y volver a preparar todas las tareas |
| `--parallel PARALLEL` | en run, tareas independientes simultáneas (default: 1) |
| `--retries RETRIES` | reintentos por tarea fallida (default: 0) |
| `--timeout TIMEOUT` | tiempo máximo por intento, en segundos |
| `--cancel-file` | marker de cancelación cooperativa personalizado |
| `--reason` | en cancel, motivo opcional |
| `--selftest` | en validate, ejecutar la validación rápida contra referencias físicas |
| `--advanced` | en validate, revisar estructura, comandos, unidades y colisiones |
| `-o, --output` | salida para dashboard, report o export |
| `--pdf` | en report, generar un informe PDF autocontenido |
| `--theme {auto,light,dark}` | tema del dashboard (default: auto) |
| `--language {es,en}` | idioma del dashboard (default: es) |
| `--both` | generar dashboard español e inglés en archivos separados |
| `--verify-environment` | en environment, comprobar el bloqueo guardado |
| `--other` | en diff, snapshot o proyecto de comparación |
| `--json` | en diff, imprimir JSON |

### `resilient`

cálculos QE recuperables ante interrupciones del servidor

**Uso:** `olla-dft resilient [-h] [--state STATE] [--pw-cmd PW_CMD] [--runtime-id RUNTIME_ID] [--checkpoint-seconds CHECKPOINT_SECONDS] [--grace-seconds GRACE_SECONDS] [--max-failures MAX_FAILURES] [--threads THREADS] [--keep KEEP] [--max-segments MAX_SEGMENTS] [--resume] [--user USER] [-o OUTPUT] {init,run,status,pause,service} target`

**Argumentos:**

- `action` {init,run,status,pause,service} — acción a realizar (ver la lista de arriba)
- `target` — input para init; directorio persistente del trabajo para las demás acciones

**Opciones:**

| Opción | Descripción |
|---|---|
| `--state` | directorio nuevo del trabajo en un disco persistente conservado |
| `--pw-cmd` | comando de QE o MPI con paralelismo fijo (default: `pw.x`) |
| `--runtime-id` | identificador de la imagen inmutable del entorno |
| `--checkpoint-seconds CHECKPOINT_SECONDS` |  (default: `900`) |
| `--grace-seconds GRACE_SECONDS` |  (default: `300`) |
| `--max-failures MAX_FAILURES` |  (default: `3`) |
| `--threads THREADS` |  (default: `1`) |
| `--keep KEEP` | guardados íntegros que conservar (mínimo 2) (default: `2`) |
| `--max-segments MAX_SEGMENTS` | detenerse tras este número de segmentos guardados; 0 significa sin límite |
| `--resume` | retirar una pausa explícita antes de continuar |
| `--user` | usuario sin privilegios del servicio systemd generado |
| `-o, --output` | archivo de servicio generado; se instala por separado |

## Apariencia y configuración

### `templates`

listar, ver o exportar plantillas

**Uso:** `olla-dft templates [-h] [-o OUTPUT] [{list,show,export}] [name]`

**Argumentos:**

- `action` {list,show,export} — list (default), show o export
- `name` — nombre de la plantilla

**Opciones:**

| Opción | Descripción |
|---|---|
| `-o, --output` | archivo JSON de salida (export) |

### `config`

ver o cambiar la configuración

**Uso:** `olla-dft config [-h] [{show,set}] [key] [value]`

**Argumentos:**

- `action` {show,set} — acción a realizar (ver la lista de arriba)
- `key` — clave de configuración, por ejemplo pseudo_dir, nproc o language
- `value` — valor que se asigna a la clave

---

*Olla-DFT 1.2.0*
