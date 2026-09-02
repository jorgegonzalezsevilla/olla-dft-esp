# Ejemplos de Olla-DFT

Cada subcarpeta trae salidas REALES de Quantum ESPRESSO, no maquetas: los
datos, las figuras y un `README.md` con los comandos `olla-dft` exactos que
las produjeron y la comparación con el experimento. Los comandos de todos
los README se validan contra la CLI en `tests/test_examples.py`, así que no
pueden quedarse obsoletos sin que una prueba falle.

| Carpeta | Qué demuestra |
|---|---|
| [`demo_Si/`](demo_Si/) | bandas, DOS y PDOS del silicio: `gen`, `run.sh`, `plot` |
| [`demo_Fe/`](demo_Fe/) | hierro bcc con polarización de espín: `gen --mag`, `dos` |
| [`demo_calculo/`](demo_calculo/) | convergencia, ecuación de estado y constantes elásticas: `converge`, `eos`, `elastic` |
| [`demo_propiedades/`](demo_propiedades/) | funciones ópticas, fonones y masa efectiva del silicio: `optics`, `phonons`, `effmass` |
| [`demo_derivadas/`](demo_derivadas/) | Debye/Slack desde las Cij, cuasi-armónica y ficha del material: `derived`, `qha`, `datasheet` |
| [`demo_laminar/`](demo_laminar/) | capas, difractograma y exfoliación de grafito: `layers`, `xrd`, `exfoliate` |
| [`demo_espectros_avanzados/`](demo_espectros_avanzados/) | XANES, U de Hubbard, electrón-fonón, desdoblamiento de bandas y VDOS: `xanes`, `hubbard`, `elph`, `unfold`, `md` |
| [`demo_tddft_balistico/`](demo_tddft_balistico/) | TDDFPT del etileno y conductancia balística de un hilo de Al: `tddft`, `ballistic` |
| [`plantillas/`](plantillas/) | la misma figura en todas las plantillas visuales: `templates`, `plot -t` |

### Estructuras sueltas

Archivos de estructura para probar cualquier comando sin buscar uno propio:

| Archivo | Qué es |
|---|---|
| `grafito.cif` | grafito hexagonal (4 átomos, 2 capas); lo usa `demo_laminar/` |
| `hbn.cif` | nitruro de boro hexagonal, otro material laminar |
| `ZnO.cif` | óxido de cinc wurtzita (4 átomos), el ejemplo del README principal |
| `POSCAR_NaCl` | cloruro de sodio en formato POSCAR (celda convencional de 8 átomos): Olla-DFT lee VASP además de CIF |

Algunos comandos de las demos citan la estructura de partida con la que se
corrieron (`Si.cif`, `NiO.cif`, `Al.cif`, `c2h4.cif`...) y que no se
incluye: sirve cualquier CIF equivalente (Materials Project, COD...)
o uno de los archivos de esta carpeta.

Por ejemplo:

    olla-dft info POSCAR_NaCl
    olla-dft prim POSCAR_NaCl -o NaCl_prim.cif
    olla-dft layers hbn.cif
    olla-dft kpath ZnO.cif
    olla-dft gen ZnO.cif -p all -o ZnO_run --insulator

Para sesiones guiadas de principio a fin, `olla-dft recetas`.
