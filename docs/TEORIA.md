# Fundamento físico de Olla-DFT

## Estructura electrónica

Este capítulo documenta la física que Olla-DFT implementa de verdad —no la que promete el manual de Quantum ESPRESSO— en los comandos que preparan un cálculo de estructura electrónica (`gen`, `kpath`, `info`, `prim`, `conv`, `supercell`) y en los que leen sus resultados (`bands`, `gap`, `dos`, `plot`, `effmass`, `fermi`, `unfold`, `wannier`, `topology`, `berry`, `hubbard`, `align`). Cada sección dice qué pregunta contesta el comando, qué fórmulas usa el código (citando la función de Python que las contiene), de qué archivo de QE o de qué constante sale cada número, y en qué casos el resultado no vale. Las constantes y los valores por omisión se han leído del código fuente de la versión 0.35.0 (`qekit/config.py`, `qekit/core/qeout.py`, y cada módulo).

---

### `olla-dft gen` — generar los inputs de pw.x y del post-proceso

**Qué responde.** Traduce una estructura (CIF, POSCAR, XYZ con celda, input o output de pw.x) en un juego de archivos de entrada coherentes para pw.x, dos.x, projwfc.x y bands.x, eligiendo por ti los cutoffs, la malla de puntos k, el camino de alta simetría, el número de bandas y el tratamiento de la ocupación.

**Fundamento para no expertos.** Un cálculo DFT de ondas planas necesita cuatro decisiones numéricas antes de empezar: (1) cuántas ondas planas usar para describir los orbitales (el *cutoff* `ecutwfc`, una energía cinética máxima en Rydberg: cuanto más alta, más fina la descripción y más caro el cálculo), (2) cuántas para la densidad (`ecutrho`, que en pseudopotenciales ultrasuaves o PAW debe ser bastante mayor), (3) cuántos puntos k muestrear en la zona de Brillouin (la "resolución" con la que se integra sobre el cristal infinito), y (4) cómo repartir los electrones entre bandas cerca del nivel de Fermi (en un aislante la ocupación es fija; en un metal se "difumina" con un *smearing* para que la suma sobre k converja). Olla-DFT toma (1) y (2) de la cabecera del propio pseudopotencial, (3) de un espaciado en el espacio recíproco (la misma idea que el `KSPACING` de VASP) y (4) de si el usuario declara que el sistema es aislante.

Un pseudopotencial es la "versión suavizada" de un átomo: sustituye al núcleo y a los electrones internos por un potencial efectivo, y solo se resuelven explícitamente los electrones de valencia. Cada archivo UPF trae recomendaciones de cutoff y el número de electrones de valencia, y Olla-DFT las lee. Para la estructura de bandas hace falta además un *camino de alta simetría*: una ruta de segmentos rectos entre puntos especiales de la zona de Brillouin (Γ, X, L…). Ese camino lo decide la biblioteca seekpath con la convención de Hinuma et al., y está referido a una celda primitiva estandarizada, por lo que `gen` usa esa celda cuando el preset incluye bandas.

**Fórmulas.**

Malla de puntos k a partir de un espaciado (`qekit/core/kpoints.py: kgrid_from_spacing`):

$$
n_i = \max\!\left(1,\ \left\lceil \frac{|\mathbf{b}_i|}{\Delta k} \right\rceil\right), \qquad \mathbf{b}_i = 2\pi\,(\mathbf{A}^{-1})^{\mathsf T}_{i}
$$

- $n_i$: número de puntos de la malla a lo largo del vector recíproco $i$ (adimensional).
- $\mathbf{b}_i$: vector de la red recíproca **incluyendo el factor $2\pi$**, en Å⁻¹; $\mathbf{A}$ es la matriz de la celda con los vectores en filas (Å).
- $\Delta k$: espaciado pedido en Å⁻¹. Por omisión `kspacing = 0.20` (scf) y `kspacing_nscf = 0.12` (nscf/DOS), leídos de `qekit/config.py: DEFAULTS`. Los niveles `--klevel` son `coarse 0.30`, `medium 0.20`, `fine 0.15`, `very-fine 0.10` y `gamma` (solo Γ).

Si a lo largo de un eje el hueco más ancho entre átomos supera `VACIO_MINIMO = 8.0` Å (`kpoints.direcciones_con_vacio`), ese $n_i$ se fuerza a 1. La malla se escribe sin desplazamiento (`0 0 0`, centrada en Γ); si es $1\times1\times1$ se escribe `K_POINTS gamma`.

Espesor de la celda y hueco de vacío (`qekit/modules/inputgen.py: espesor_celda`, `hueco_vacio`):

$$
d_i = \frac{V}{|\mathbf{a}_j \times \mathbf{a}_k|}, \qquad h_{\text{Å}} = d_i \cdot \max_m \left(f^{(i)}_{m+1} - f^{(i)}_m\right)
$$

- $d_i$: altura de la celda a lo largo de la normal al plano de los otros dos vectores (Å); $V$ es el volumen (Å³).
- $h_{\text{Å}}$: el mayor hueco entre coordenadas fraccionarias ordenadas $f^{(i)}$ a lo largo del eje $i$ (con el hueco que cruza el borde periódico incluido), convertido a Å.

Cutoffs recomendados (`qekit/core/pseudo.py: recommend_cutoffs`):

$$
E_{\text{wfc}} = \max_s E^{\text{UPF}}_{\text{wfc},s}, \qquad
E_{\rho} = \max\!\left(\max_s E^{\text{UPF}}_{\rho,s},\ 4\,E_{\text{wfc}}\right)
$$

- $E^{\text{UPF}}_{\text{wfc},s}$, $E^{\text{UPF}}_{\rho,s}$: cutoffs sugeridos en la cabecera del UPF de la especie $s$ (Ry), leídos por `pseudo.suggested_cutoffs` (atributos `wfc_cutoff`/`rho_cutoff` en UPF v2 o el texto "Suggested minimum cutoff for wavefunctions/charge density" en UPF v1). Valores $\le 1$ se ignoran.
- Si ningún UPF declara cutoff se usan `ecutwfc = 60.0` Ry y `ecutrho = dual × ecutwfc` con `dual = 8` (config). El suelo $4E_{\text{wfc}}$ es el mínimo físico para ondas planas.

Número de bandas estimado para nscf/bands (`inputgen._estimate_nbnd`):

$$
n_{\text{bnd}} = \left\lceil 1.25\cdot\frac{N_{\text{el}}}{2} + 4 \right\rceil, \qquad N_{\text{el}} = \sum_{\text{átomos}} Z^{\text{UPF}}_{\text{val}}
$$

Con `--nspin 2` se amplía a $\lfloor 1.2\,n_{\text{bnd}}\rfloor + 2$. Si algún UPF no declara `z_valence`, no se escribe `nbnd` y pw.x usa su valor por defecto.

Paso de tiempo de la MD (`inputgen.build_pw_input`): $\mathrm{dt}_{\text{Ry}} = \mathrm{dt}_{\text{fs}} / 0.048378$, porque pw.x pide `dt` en unidades atómicas de Rydberg (`_FS_POR_UA = 4.8378e-2` fs).

Corrección dipolar (`inputgen._region_vacio`): el máximo de la sierra se coloca en el centro del hueco de vacío, `emaxpos = centro`, y su bajada ocupa `eopreg = clip(hueco/3, 0.02, 0.20)` (fracciones del eje). Se exige $h_{\text{Å}} \ge 5$ Å.

Coste estimado de un híbrido (`inputgen.generate`), medido en silicio de 2 átomos: $\text{factor} \approx 3 + 2.6\,n_q$, con $n_q = n_{q1}n_{q2}n_{q3}$ la malla del intercambio exacto. Es una regla empírica, no una fórmula.

**Cómo lo calcula Olla-DFT.**

1. `qekit/cli.py: _cmd_gen` lee la estructura con `qekit/core/structure.py: load` (ASE; para `POSCAR/CONTCAR` fuerza `format="vasp"`, y si el archivo trae varias imágenes se queda con la última). Combina `--klevel`/`--kspacing`, `--mag` (que activa `nspin=2`), `--hubbard`, `--soc`, `--functional`, `--exx-grid` y las opciones de MD en un `GenOptions`.
2. `qekit/modules/inputgen.py: generate` decide la celda de trabajo: si el preset es `bands` o `all`, llama a `kpoints.get_kpath` (seekpath) y usa la **celda primitiva estandarizada** que devuelve; con `--primitive` usa `structure.primitive` (spglib); si no, la celda tal cual.
3. `qekit/core/pseudo.py: resolve` busca un UPF por elemento en `pseudo_dir` (archivo cuyo nombre empieza por el símbolo seguido de un carácter no alfabético, extensión `.upf` sin distinguir mayúsculas). Respeta `--pseudo El=archivo`, y `_coherencia_de_funcional` reelige los que haga falta para que todos los pseudos compartan funcional (preferencia `PBE > PBESOL > REVPBE > PZ > BLYP`).
4. `pseudo.recommend_cutoffs` fija `ecutwfc`/`ecutrho`; el usuario puede anularlos con `--ecutwfc`/`--ecutrho`.
5. `kpoints.kgrid_from_spacing` produce las mallas scf y nscf; `_estimate_nbnd` el número de bandas.
6. `inputgen.build_pw_input` escribe `&CONTROL` (con `tprnfor`, `tstress`, `outdir='./out'`), `&SYSTEM` (`ibrav=0`, cutoffs, ocupaciones, espín, SOC, Hubbard, `tot_charge`, dipolo, híbrido, `nosym`/`noinv`), `&ELECTRONS` (`conv_thr=1e-8`, `mixing_beta=0.4`, `electron_maxstep=200`), `&IONS`/`&CELL` según el preset (BFGS para relax, Verlet para MD con el termostato pedido, `press_conv_thr=0.05`), y las tarjetas `ATOMIC_SPECIES` (masas de `ase.data.atomic_masses`), `ATOMIC_POSITIONS crystal`, `CELL_PARAMETERS angstrom`, la tarjeta `HUBBARD` si `--hubbard-style card`, y `K_POINTS`.
7. Para bandas, `kpoints.kpath_card` escribe `K_POINTS crystal_b` con `band_points` puntos por tramo (20 por omisión) y `KPATH.txt` con las etiquetas; `build_bandsx_input` escribe `bands_pp.in` (`lsym=.true.`). Para DOS, `build_dos_input` y `build_projwfc_input` escriben `dos.in` y `projwfc.in` con `DeltaE = 0.02` eV.
8. `build_run_script` y `build_run_python_script` escriben `run.sh` (con `set -e -o pipefail` y `mpirun -np $NP`) y `run.py` con el orden `pw.x → dos.x/projwfc.x → pw.x (bands) → bands.x`.

**De dónde sale cada dato.**

| Dato | Origen | Detalle |
|---|---|---|
| Celda y posiciones | archivo del usuario (CIF/POSCAR/…) | `structure.load` vía `ase.io.read` |
| Celda primitiva y k-path | biblioteca seekpath | `kpoints.get_kpath` con `symprec = 1e-4` Å (`structure.SYMPREC`) |
| Cutoffs sugeridos | cabecera del UPF (`wfc_cutoff`, `rho_cutoff`) | `pseudo.suggested_cutoffs`; solo lee los primeros 20 000 caracteres |
| Electrones de valencia | `z_valence` del UPF | `pseudo.z_valence`, usado en `_estimate_nbnd` |
| Tipo y relativismo del pseudo | `pseudo_type`, `relativistic` del UPF | `pseudo.pseudo_type`, `pseudo.relativistic` |
| `ecutwfc`, `dual`, `kspacing`, `degauss`, `smearing`, `nproc` | `~/.config/qekit/config.ini` o `config.DEFAULTS` | 60 Ry, 8, 0.20 Å⁻¹, 0.01 Ry, `cold`, 4 |
| Masas atómicas | `ase.data.atomic_masses` | tarjeta `ATOMIC_SPECIES` |
| Parámetros de híbridos | tabla `inputgen.HIBRIDOS` | HSE: `exx_fraction 0.25`, `screening_parameter 0.106` bohr⁻¹; PBE0 0.25; B3LYP 0.20; Gau-PBE 0.24 |
| Orbital de Hubbard (tarjeta) | número atómico (`inputgen._orbital_hubbard`) | 3d (Z 21–30), 4d (39–48), 5d (72–80), 4f (57–71), 5f (89–103), `2p` para el resto |
| Conversión fs → u.a. Rydberg | constante `inputgen._FS_POR_UA` | 4.8378e-2 fs |

**Límites y trampas.**

- Los cutoffs "automáticos" son los que **sugiere el UPF**, no una convergencia: el reporte dice "(automático)". Si el UPF no los declara, se cae a 60 Ry / 480 Ry sin avisar más allá del reporte.
- `--soc` escribe `noncolin=.true.` y `lspinorb=.true.` solo si todos los pseudos declaran `relativistic="full"`: `inputgen.generate` llama a `sweep.check_soc_pseudos` antes de escribir nada y, si alguno es escalar-relativista o no lo declara, se detiene con "el acoplamiento espín-órbita necesita pseudopotenciales TOTALMENTE RELATIVISTAS (relativistic='full'), y estos no lo son". El motivo, citado en el propio error: con pseudos escalares lspinorb "devuelve un desdoblamiento espín-órbita de cero que parece un resultado válido". `--soc` y `--nspin 2` se rechazan juntos.
- `--hubbard-style legacy` (por omisión) escribe `lda_plus_u` y `Hubbard_U(i)`, sintaxis retirada en QE ≥ 7.1; para esas versiones hace falta `--hubbard-style card`, que escribe la tarjeta `HUBBARD (ortho-atomic)`. El orbital de la tarjeta se deduce solo del número atómico.
- Con híbridos, `nqx` **tiene que dividir** la malla de k; si no, el comando se detiene con "la malla de intercambio exacto tiene que DIVIDIR la de k". El reporte avisa además de que con `1x1x1` "el resultado va a salir claramente sobrestimado" y de que pw.x no puede hacer un `calculation='bands'` con EXX.
- Sin `--mag`, `--nspin 2` arranca con magnetización cero y el reporte avisa: "sin magnetización inicial el cálculo suele converger a la solución no magnética".
- La corrección dipolar exige un hueco de vacío ≥ 5 Å; si no, error: "la corrección dipolar necesita vacío en la dirección …".
- La malla scf de un preset `bands` se calcula sobre la celda primitiva de seekpath, que puede no ser la que el usuario dio; el reporte lo dice ("AVISO: se usó la celda primitiva estandarizada").
- Para MD se fuerza `nosym`, y se avisa si hay menos de 20 átomos o `dt > 2` fs.
- `tot_charge` se compensa con un fondo uniforme; el reporte recuerda que la energía de una celda cargada no es comparable con la neutra.
- La malla es siempre uniforme y centrada en Γ (sin desplazamiento), como dice el docstring de `kpoints.py`: no es una malla Monkhorst-Pack desplazada, y con $n$ par contiene Γ donde la de MP no lo haría.

**Referencias.**

- Y. Hinuma, G. Pizzi, Y. Kumagai, F. Oba, I. Tanaka, *Comput. Mater. Sci.* **128**, 140 (2017) — convención de k-paths de seekpath. DOI 10.1016/j.commatsci.2016.10.015.
- A. Togo, I. Tanaka, "Spglib: a software library for crystal symmetry search", arXiv:1808.01590 (2018).
- N. Marzari, D. Vanderbilt, A. De Vita, M. C. Payne, *Phys. Rev. Lett.* **82**, 3296 (1999) — smearing "cold". DOI 10.1103/PhysRevLett.82.3296.
- J. Heyd, G. E. Scuseria, M. Ernzerhof, *J. Chem. Phys.* **118**, 8207 (2003) — HSE.
- L. Bengtsson, *Phys. Rev. B* **59**, 12301 (1999) — corrección dipolar en losas.
- G. Prandini et al., *npj Comput. Mater.* **4**, 72 (2018) — biblioteca SSSP de pseudopotenciales y cutoffs.

---

### `olla-dft kpath` — camino de alta simetría

**Qué responde.** Cuál es el camino estándar por la zona de Brillouin (y en qué celda está referido) que hay que usar para dibujar la estructura de bandas de esta estructura.

**Fundamento para no expertos.** La zona de Brillouin es la "celda unidad" del espacio de los vectores de onda; las bandas $E(\mathbf k)$ se dibujan a lo largo de una ruta que pasa por sus puntos de mayor simetría, porque ahí es donde las bandas se tocan, se cruzan o tienen sus extremos. Qué puntos y en qué orden depende del grupo espacial y de la forma de la celda, y hay varias convenciones incompatibles en la literatura. Olla-DFT delega la elección en seekpath (convención de Hinuma et al., la misma que usa Materials Cloud), que además **estandariza la celda**: las coordenadas de los puntos especiales solo valen en esa celda primitiva, no necesariamente en la que el usuario tiene en su CIF.

**Fórmulas.** No hay fórmulas propias: el comando llama a `seekpath.get_path` y transcribe su resultado. La única aritmética es el criterio de "misma celda" (`qekit/core/kpoints.py: get_kpath`):

$$
\text{cell\_changed} = \neg\left[N_{\text{prim}} = N_{\text{in}} \ \wedge\ \max_{ij} |A^{\text{prim}}_{ij} - A^{\text{in}}_{ij}| \le 10^{-5}\ \text{Å}\right]
$$

- $A^{\text{prim}}$, $A^{\text{in}}$: matrices de celda (Å) de la primitiva de seekpath y de la entrada; $N$: número de átomos.

**Cómo lo calcula Olla-DFT.**

1. `qekit/cli.py: _cmd_kpath` carga la estructura con `structure.load`.
2. `kpoints.get_kpath` convierte a la tupla de spglib (`structure.to_spglib_cell`), llama a `seekpath.get_path(..., symprec=1e-4)` y reconstruye la celda primitiva (`structure.from_spglib_cell`) a partir de `primitive_lattice`, `primitive_positions`, `primitive_types`.
3. `kpoints.kpath_text` imprime grupo espacial (`spacegroup_international`, `spacegroup_number`), el camino compactado (`Γ — X — U | K — Γ …`), las coordenadas fraccionarias de cada punto (`point_coords`) con las etiquetas "bonitas" de `pretty_label` (GAMMA→Γ, DELTA_0→Δ0), y un aviso si la celda cambió.

**De dónde sale cada dato.**

| Dato | Origen | Detalle |
|---|---|---|
| Camino y coordenadas | biblioteca seekpath (`get_path`) | claves `path`, `point_coords` |
| Grupo espacial | seekpath (spglib por debajo) | `spacegroup_international`, `spacegroup_number` |
| Tolerancia de simetría | constante `structure.SYMPREC` | 1e-4 Å |
| Celda primitiva | seekpath | `primitive_lattice/positions/types` |

**Límites y trampas.**

- Las coordenadas están **en la celda primitiva estandarizada**. Si `cell_changed` es verdadero, el texto avisa: "el k-path está referido a la celda primitiva estandarizada, que difiere de la celda de entrada. Usa esa celda primitiva en el cálculo de bandas". Usar esas coordenadas en la celda original da un camino equivocado sin ningún error.
- Con `symprec = 1e-4` Å, una estructura relajada con ruido numérico puede perder simetría y recibir un camino de un grupo espacial más bajo; no hay opción de línea de comandos para cambiar la tolerancia.
- seekpath no tiene en cuenta el magnetismo ni el espín-órbita al elegir el grupo espacial.

**Referencias.**

- Y. Hinuma, G. Pizzi, Y. Kumagai, F. Oba, I. Tanaka, *Comput. Mater. Sci.* **128**, 140 (2017). DOI 10.1016/j.commatsci.2016.10.015.
- W. Setyawan, S. Curtarolo, *Comput. Mater. Sci.* **49**, 299 (2010) — la convención alternativa que seekpath **no** usa.

---

### `olla-dft info` — estructura y simetría

**Qué responde.** Qué contiene el archivo de estructura: fórmula, parámetros de red, volumen, grupo espacial, grupo puntual, posiciones de Wyckoff y cuántos átomos tendría la celda primitiva.

**Fundamento para no expertos.** Antes de calcular nada conviene saber si la celda es la mínima posible (la primitiva) o una celda mayor (convencional o supercelda), porque el coste del cálculo crece con el número de átomos, y si la estructura tiene la simetría que uno cree. spglib compara cada átomo con las operaciones de simetría candidatas dentro de una tolerancia y devuelve el grupo espacial en notación internacional (por ejemplo `Fd-3m`, N.º 227 para el silicio), el símbolo de Hall, el grupo puntual y la letra de Wyckoff de cada sitio (una etiqueta de qué tipo de posición de simetría ocupa cada átomo).

**Fórmulas.** No hay fórmulas más allá de la geometría de la celda que ASE calcula (`atoms.cell.cellpar()` devuelve $a, b, c, \alpha, \beta, \gamma$ y `atoms.cell.volume` el volumen $V = |\det \mathbf A|$ en Å³).

**Cómo lo calcula Olla-DFT.**

1. `qekit/cli.py: _cmd_info` → `structure.load`.
2. `qekit/core/structure.py: info_text` llama a `symmetry_dataset` (`spglib.get_symmetry_dataset` con `symprec = 1e-4`) y a `primitive` (`spglib.standardize_cell(to_primitive=True)`) para contar los átomos de la primitiva.
3. Imprime fórmula, composición, número de átomos, volumen, parámetros de red, grupo espacial (`international`, `number`), símbolo de Hall, grupo puntual, átomos en la primitiva, posiciones de Wyckoff (conjunto ordenado de `ds.wyckoffs`) y los vectores de celda.

**De dónde sale cada dato.**

| Dato | Origen | Detalle |
|---|---|---|
| Grupo espacial, Hall, grupo puntual, Wyckoff | biblioteca spglib | `structure.symmetry_dataset` |
| Parámetros de red y volumen | ASE (`Cell.cellpar`, `Cell.volume`) | a, b, c en Å; ángulos en grados |
| Átomos en la primitiva | spglib `standardize_cell` | `structure.primitive` |
| Tolerancia | `structure.SYMPREC` | 1e-4 Å |

**Límites y trampas.**

- Si spglib no puede determinar la simetría, el comando falla con `RuntimeError("spglib no pudo determinar la simetría de la estructura")`.
- La tolerancia es fija (1e-4 Å); no hay `--symprec`.
- Solo se listan las letras de Wyckoff distintas, no cuántos átomos hay en cada una.

**Referencias.**

- A. Togo, I. Tanaka, "Spglib: a software library for crystal symmetry search", arXiv:1808.01590 (2018).
- International Tables for Crystallography, Vol. A (IUCr) — notación de grupos espaciales y posiciones de Wyckoff.

---

### `olla-dft prim` — celda primitiva estandarizada

**Qué responde.** Cuál es la celda más pequeña que reproduce el cristal por traslación, escrita en la orientación estándar de spglib.

**Fundamento para no expertos.** Muchos CIF vienen en la celda convencional (la que hace visible la simetría, por ejemplo el cubo de 8 átomos del silicio), pero el cálculo solo necesita la celda primitiva (2 átomos en el silicio). Reducirla ahorra un factor igual al cociente de átomos en el coste, sin cambiar la física. "Estandarizada" quiere decir que spglib la reorienta y la expresa con la elección de vectores que fija la convención de la International Tables, de modo que dos entradas equivalentes den la misma salida.

**Fórmulas.** No hay aritmética propia: es `spglib.standardize_cell(cell, to_primitive=True, symprec=1e-4)`.

**Cómo lo calcula Olla-DFT.**

1. `qekit/cli.py: _cmd_prim` → `structure.load`.
2. `structure.primitive` → spglib → `structure.from_spglib_cell` (un `Atoms` con `pbc=True`).
3. `structure.convert` escribe el resultado según la extensión de `-o` (por omisión `primitive.cif`): CIF, POSCAR/`.vasp` (con `direct=True, sort=True`) o cualquier formato que ASE deduzca.

**De dónde sale cada dato.**

| Dato | Origen | Detalle |
|---|---|---|
| Celda primitiva | spglib `standardize_cell(to_primitive=True)` | `structure.primitive` |
| Tolerancia | `structure.SYMPREC` | 1e-4 Å |
| Formato de salida | extensión del archivo | `structure.convert` |

**Límites y trampas.**

- La celda primitiva de spglib **no es** necesariamente la misma que la de seekpath que usa `gen -p bands` (seekpath aplica su propia estandarización adicional); para bandas hay que dejar que `gen` elija.
- Al escribir POSCAR se **reordenan** los átomos por especie (`sort=True`); si el usuario tenía un orden concreto (por ejemplo para un `--displace`), se pierde.
- Si spglib falla: `RuntimeError("spglib no pudo estandarizar la celda")`.

**Referencias.**

- A. Togo, I. Tanaka, arXiv:1808.01590 (2018).

---

### `olla-dft conv` — celda convencional estandarizada

**Qué responde.** Cuál es la celda convencional (la que muestra la simetría completa del sistema cristalino) de la estructura.

**Fundamento para no expertos.** Es la operación inversa a `prim`: partir de cualquier celda y obtener la celda "de libro" (cúbica para el silicio, hexagonal para el grafito), útil para construir superficies, superceldas o para comparar con datos de difracción, aunque tenga más átomos de los estrictamente necesarios para el cálculo.

**Fórmulas.** No hay aritmética propia: `spglib.standardize_cell(cell, to_primitive=False, symprec=1e-4)`.

**Cómo lo calcula Olla-DFT.**

1. `qekit/cli.py: _cmd_conv` → `structure.load`.
2. `structure.conventional` → spglib → `from_spglib_cell`.
3. `structure.convert` escribe `-o` (por omisión `conventional.cif`).

**De dónde sale cada dato.**

| Dato | Origen | Detalle |
|---|---|---|
| Celda convencional | spglib `standardize_cell(to_primitive=False)` | `structure.conventional` |
| Tolerancia | `structure.SYMPREC` | 1e-4 Å |

**Límites y trampas.**

- Igual que `prim`: tolerancia fija, reordenación de átomos en POSCAR, error si spglib no puede estandarizar.
- Para una celda de baja simetría (triclínica) la "convencional" coincide con la primitiva y el comando no cambia nada.

**Referencias.**

- A. Togo, I. Tanaka, arXiv:1808.01590 (2018).

---

### `olla-dft supercell` — construir una supercelda

**Qué responde.** La estructura repetida $n_x \times n_y \times n_z$ veces a lo largo de sus tres vectores de celda.

**Fundamento para no expertos.** Una supercelda es varias celdas pegadas y tratadas como una sola unidad periódica. Hace falta para poner un defecto, un dopante o una molécula adsorbida a una concentración baja, para dinámica molecular, o para calcular fonones por desplazamientos finitos. El precio es que la zona de Brillouin se encoge en el mismo factor y las bandas se "pliegan" (véase `unfold`).

**Fórmulas.** `atoms.repeat((nx, ny, nz))` de ASE: la nueva celda es $\mathbf a'_i = n_i \mathbf a_i$ y cada átomo se copia en las $n_x n_y n_z$ traslaciones $\sum_i m_i \mathbf a_i$ con $0 \le m_i < n_i$.

**Cómo lo calcula Olla-DFT.**

1. `qekit/cli.py: _cmd_supercell` → `structure.load`.
2. `structure.supercell` valida que los tres factores sean ≥ 1 (si no, `ErrorDeUso("los factores de la supercelda deben ser >= 1")`) y llama a `Atoms.repeat`.
3. `structure.convert` escribe `-o` (por omisión `supercell.cif`).

**De dónde sale cada dato.**

| Dato | Origen | Detalle |
|---|---|---|
| Factores $n_x, n_y, n_z$ | parámetros del usuario (posicionales) | enteros ≥ 1 |
| Repetición | ASE `Atoms.repeat` | solo múltiplos diagonales |

**Límites y trampas.**

- Solo superceldas **diagonales** (múltiplos de cada vector); no se pueden construir matrices generales $\mathbf A' = \mathbf M \mathbf a$ como las que `unfold` sí sabe reconocer.
- No se reduce por simetría ni se comprueba que la supercelda sea "razonable" (por ejemplo, cúbica).

**Referencias.**

- A. H. Larsen et al., "The atomic simulation environment—a Python library for working with atoms", *J. Phys.: Condens. Matter* **29**, 273002 (2017). DOI 10.1088/1361-648X/aa680e.

---
### `olla-dft bands` — estructura de bandas y band gap

**Qué responde.** Cómo varía la energía de cada estado electrónico a lo largo del camino de alta simetría, si el sistema tiene gap o es metálico, dónde están el máximo de la banda de valencia (VBM) y el mínimo de la de conducción (CBM), si el gap es directo o indirecto, y —con `--fat`— de qué orbital atómico "está hecha" cada banda.

**Fundamento para no expertos.** En un cristal los electrones no tienen energías discretas sino *bandas*: para cada vector de onda $\mathbf k$ (una "dirección y longitud de onda" del electrón) hay una lista de energías permitidas $\varepsilon_n(\mathbf k)$. Dibujarlas a lo largo de un camino en la zona de Brillouin da la figura clásica de "espaguetis". El *band gap* es la distancia entre la banda ocupada más alta y la vacía más baja. Si el máximo de una y el mínimo de la otra están en el mismo $\mathbf k$ el gap es *directo* (el material absorbe y emite luz eficientemente); si no, *indirecto*. Si alguna banda cruza el nivel de Fermi (la energía hasta la que se llenan los estados), es un metal y no hay gap.

Las *fatbands* ("bandas gordas") responden a otra pregunta: qué peso tiene cada orbital atómico (el $d$ del níquel, el $p$ del oxígeno) en cada estado. projwfc.x proyecta cada función de onda sobre orbitales atómicos y escribe los pesos; Olla-DFT los dibuja como puntos de tamaño proporcional al peso encima de las bandas.

**Fórmulas.**

Conversión de unidades del XML de pw.x (`qekit/core/qeout.py: read_xml`):

$$
E_{\text{eV}} = E_{\text{Ha}} \cdot 27.211386245988, \qquad
\mathbf{k}_{\text{Å}^{-1}} = \mathbf{k}_{2\pi/a} \cdot \frac{2\pi}{a_{\text{bohr}}\cdot 0.529177210903}, \qquad
\mathbf{k}_{\text{frac}} = \mathbf{k}_{\text{Å}^{-1}}\, \mathbf{B}^{-1}
$$

- $a_{\text{bohr}}$: `alat` del XML (bohr). $\mathbf B$: matriz con los vectores recíprocos $\mathbf b_i = 2\pi(\mathbf A^{-1})^{\mathsf T}_i$ en filas (Å⁻¹).

Distancia acumulada en el eje x (`qekit/modules/bands.py: _build_kdist`):

$$
x_0 = 0, \qquad x_i = x_{i-1} + \begin{cases} 0 & i \in \text{breaks} \\ |\mathbf{k}_i - \mathbf{k}_{i-1}| & \text{si no} \end{cases}
$$

- `breaks`: índices donde dos puntos especiales aparecen consecutivos (una discontinuidad `U|K` del camino), detectados por `_detect_breaks`.

Clasificación metal / aislante (`bands.analyze_gap`), con `CROSS_TOL = 1e-6` eV y referencia $E_{\text{ref}}$:

$$
\text{cruza}_n = \left[\min_{\mathbf k}\varepsilon_n < E_{\text{ref}} - \delta\right] \wedge \left[\max_{\mathbf k}\varepsilon_n > E_{\text{ref}} + \delta\right]
$$

Si alguna banda cruza, es metal. Si no, $n_v = \max\{n : \max_{\mathbf k}\varepsilon_n \le E_{\text{ref}}+\delta\}$ y $n_c = \min\{n : \min_{\mathbf k}\varepsilon_n > E_{\text{ref}}-\delta\}$ (forzando $n_c \ge n_v+1$), y

$$
E_{\text{VBM}} = \max_{\mathbf k}\varepsilon_{n_v}(\mathbf k), \quad
E_{\text{CBM}} = \min_{\mathbf k}\varepsilon_{n_c}(\mathbf k), \quad
E_g = E_{\text{CBM}} - E_{\text{VBM}}, \quad
E_g^{\text{dir}} = \min_{\mathbf k}\left[\varepsilon_{n_c}(\mathbf k) - \varepsilon_{n_v}(\mathbf k)\right]
$$

El gap es directo si $\arg\max\varepsilon_{n_v} = \arg\min\varepsilon_{n_c}$ (mismo índice de punto k, no misma coordenada).

Referencia $E_{\text{ref}}$ (en este orden): `<fermi_energy>` del XML; si no existe, `<highestOccupiedLevel>`; si tampoco y `nspin = 1`, el punto medio $\tfrac12[\max_{\mathbf k}\varepsilon_{n_{occ}-1} + \min_{\mathbf k}\varepsilon_{n_{occ}}]$ con $n_{occ} = \mathrm{round}(N_{\text{el}}/2)$; en último caso, la mediana de todas las energías.

Cero de energías en la gráfica y en los datos exportados (`bands.reference_energy`): `--ref auto` usa el VBM si hay gap y $E_F$ si es metal; `fermi`, `vbm`, `none` como indican.

Peso de un selector en las fatbands (`bands.peso_de`):

$$
w_{n}(\mathbf k) = \sum_{i \in \text{selector}} |\langle \phi_i | \psi_{n\mathbf k}\rangle|^2
$$

- $|\langle \phi_i | \psi_{n\mathbf k}\rangle|^2$: los coeficientes `psi = 0.498*[# 1] + …` de la salida de texto de projwfc.x para el estado atómico $i$, leídos por `bands.leer_proyecciones`. **No se normaliza**: la parte que falta hasta 1 es la de la función de onda que no cae en ninguna esfera atómica, y `report_fat` la cuantifica como $1 - \langle\sum_i w_i\rangle$.

**Cómo lo calcula Olla-DFT.**

1. `qekit/cli.py: _cmd_bands` → `bands.load(path, prefix)`.
2. `qeout.find_xml` localiza el XML (`./out/*.xml`, `./*.xml` o `*.save/data-file-schema.xml`, comprobando que contenga "espresso"); `qeout.read_xml` lee `<atomic_structure>` (celda en bohr), `<band_structure>` (`nbnd`, `nelec`, `lsda`, `noncolin`, `fermi_energy`, `highestOccupiedLevel`, `lowestUnoccupiedLevel`), y cada `<ks_energies>` (`k_point` con `weight`, `eigenvalues`, `occupations`). Con `lsda`, la lista de autovalores por k se parte en dos mitades (up/down).
3. Las etiquetas vienen de `KPATH.txt` (`qeout.read_kpath_labels`) o, si no existe, de la tarjeta `K_POINTS crystal_b` de `bands.in` con comentarios `! G` (`qeout.read_crystal_b_card`). `qeout.match_labels_to_kpoints` las asigna a índices con tolerancia `1e-3` en fraccionarias, avanzando siempre hacia adelante y tolerando traslaciones de red recíproca.
4. `bands.analyze_gap` por canal de espín; `bands.gap_report` imprime el resumen y el recordatorio de que "los funcionales GGA/LDA subestiman el gap sistemáticamente (típicamente 30–50 %)".
5. Con `--fat`: `bands.leer_proyecciones` lee `projwfc.out` (o `proj.out`, `projwfc_bands.out`) del mismo cálculo de bandas; `comprobar_compatibilidad` exige el mismo número de puntos k; `peso_de` suma los estados que encajan con el selector (`Ni`, `Ni-d`, `d`, `atomo:3`).
6. `bands.export` escribe `BAND.dat` (o `BAND_up.dat`/`BAND_dw.dat`), `KLABELS.dat` y `BAND_GAP.txt`; `bands.plot` dibuja con matplotlib (bandas en tinta, espín ↓ a trazos, VBM como círculo y CBM como cuadrado, fatbands como `scatter` con `s = w · fat_scale`, `fat_scale = 55`).

**De dónde sale cada dato.**

| Dato | Origen | Detalle |
|---|---|---|
| Autovalores, k-points, pesos | `prefix.xml` de pw.x (`<ks_energies>`) | `qeout.read_xml`; Ha → eV |
| Energía de Fermi | `<fermi_energy>` en `prefix.xml` | solo si el scf usó smearing |
| HOMO / LUMO | `<highestOccupiedLevel>` / `<lowestUnoccupiedLevel>` | ocupaciones fijas |
| Número de electrones, bandas, espín | `<nelec>`, `<nbnd>`, `<lsda>`, `<noncolin>` | `nbnd` se recalcula a partir de la longitud de la lista de autovalores |
| Celda y `alat` | `<atomic_structure alat=…>` y `<cell>` | bohr → Å con 0.529177210903 |
| Etiquetas de alta simetría | `KPATH.txt` de `olla-dft gen` o `bands.in` | tolerancia de emparejamiento 1e-3 |
| Pesos orbitales (fatbands) | `projwfc.out` (bloques `psi = …`) | `bands.leer_proyecciones`; `state #` da átomo, elemento y $l$ |
| Hartree en eV, bohr en Å | constantes `qeout.HARTREE_EV`, `qeout.BOHR_ANG` | 27.211386245988 eV; 0.529177210903 Å (CODATA 2018) |

**Límites y trampas.**

- En un cálculo `bands` con ocupaciones fijas el XML puede no traer `<fermi_energy>`; entonces la referencia es `<highestOccupiedLevel>`, que en QE se hereda del scf. Si el scf fue con smearing, $E_F$ puede caer en mitad del gap o dentro de una banda plana; una banda que toque $E_F \pm 10^{-6}$ eV por un solo punto se clasifica como metal.
- Si `nbnd` solo cubre las bandas ocupadas, el reporte dice "No hay bandas de conducción en el cálculo (aumenta nbnd para obtener el gap)".
- Con `nspin = 2` cada canal se analiza por separado con la **misma** referencia; el reporte no calcula el gap global entre canales distintos (por ejemplo VBM up y CBM down).
- La gráfica y `--ref auto` usan siempre el análisis del canal 0 (espín up) para decidir el cero y marcar los extremos.
- "Directo" se decide comparando **índices** de punto k; dos puntos k equivalentes por simetría en índices distintos cuentan como indirecto.
- Fatbands: si projwfc.x se corrió sobre el nscf de la DOS y no sobre las bandas, `comprobar_compatibilidad` se detiene: "las bandas tienen N puntos k y las proyecciones M. No son del mismo cálculo". Si más del 10 % del peso medio no cae en esferas atómicas, `report_fat` avisa: "De media, un X % de cada función de onda NO cae dentro de ninguna esfera atómica".
- Los pesos de estados con $l>3$ se etiquetan `l4`, `l5`… y no se pueden seleccionar por letra.
- Los cálculos con SOC (`noncolin`) se leen como un solo canal; los pesos de projwfc con $j$ (`p_j1.5`) se agrupan solo por la letra del orbital.

**Referencias.**

- P. Giannozzi et al., *J. Phys.: Condens. Matter* **29**, 465901 (2017) — Quantum ESPRESSO (formato XML, projwfc.x). DOI 10.1088/1361-648X/aa8f79.
- J. P. Perdew, M. Levy, *Phys. Rev. Lett.* **51**, 1884 (1983) y L. J. Sham, M. Schlüter, *Phys. Rev. Lett.* **51**, 1888 (1983) — por qué DFT subestima el gap.
- CODATA 2018, E. Tiesinga et al., *Rev. Mod. Phys.* **93**, 025010 (2021) — constantes.

---

### `olla-dft gap` — solo el reporte de band gap

**Qué responde.** Lo mismo que la parte de análisis de `bands` —metal o no, VBM, CBM, gap fundamental y gap directo mínimo por canal de espín— sin exportar datos ni dibujar.

**Fundamento para no expertos.** Es la pregunta más frecuente que se le hace a un cálculo de bandas ("¿cuánto vale el gap?") desacoplada de la figura. Sirve igual sobre un `scf`, un `nscf` de malla o un `bands` de camino: lee cualquier XML de pw.x con autovalores. Con una malla, el gap que sale es el de los puntos muestreados, que puede ser mayor que el verdadero si el extremo cae entre puntos.

**Fórmulas.** Exactamente las de `bands.analyze_gap` descritas en `olla-dft bands` (clasificación con `CROSS_TOL = 1e-6` eV, $E_g = E_{\text{CBM}} - E_{\text{VBM}}$, $E_g^{\text{dir}} = \min_{\mathbf k}[\varepsilon_{n_c} - \varepsilon_{n_v}]$).

**Cómo lo calcula Olla-DFT.**

1. `qekit/cli.py: _cmd_gap` → `bands.load(path, prefix)` (lee el XML y, si existen, `KPATH.txt` o `bands.in` para etiquetar los puntos).
2. `bands.gap_report` recorre `range(nspin)` llamando a `analyze_gap`, imprime ruta del XML, tipo de cálculo, `nbnd`, `nk`, `nelec`, $E_F$ si existe, y por canal el resultado.

**De dónde sale cada dato.**

| Dato | Origen | Detalle |
|---|---|---|
| Autovalores y k | `prefix.xml` de pw.x | `qeout.read_xml` |
| Referencia | `<fermi_energy>` → `<highestOccupiedLevel>` → conteo de electrones → mediana | `bands.analyze_gap` |
| Etiqueta del punto k del extremo | `KPATH.txt` / `bands.in`, o las coordenadas fraccionarias | `bands._label_for` |

**Límites y trampas.**

- Sobre una malla scf/nscf reducida por simetría, "directo" solo puede detectarse si VBM y CBM caen en el **mismo índice** de la lista de puntos irreducibles.
- Un XML sin `<output>` (cálculo no terminado) produce `FaltanDatos("… no contiene una sección <output>")`.
- `gap_report` llama dos veces a `analyze_gap` por canal (una para el reporte y otra para decidir si imprimir el recordatorio GGA); es solo coste, no afecta al resultado.

**Referencias.**

- Las mismas que `olla-dft bands`.

---

### `olla-dft dos` — densidad de estados total y proyectada

**Qué responde.** Cuántos estados electrónicos hay por unidad de energía (DOS), cómo se reparten entre elementos y orbitales (PDOS), cuánto vale la DOS en el nivel de Fermi, y —con `--dband`— el centro, la anchura y el llenado de una banda proyectada (el "centro de la banda d" de la catálisis).

**Fundamento para no expertos.** La DOS es el histograma de las energías de todos los estados de la zona de Brillouin: donde es alta hay muchos estados, donde es cero hay gap. dos.x la calcula a partir de los autovalores del nscf (con el método de los tetraedros que `gen -p dos` pide); projwfc.x descompone cada estado en orbitales atómicos y da una PDOS por átomo y orbital. Sumando los archivos por elemento y por letra orbital ($s, p, d, f$) se obtiene la descomposición química que se publica.

El *centro de la banda d* es la energía media de la PDOS $d$ de un metal de transición respecto al nivel de Fermi. Es un descriptor empírico: cuanto más cerca del Fermi (menos negativo), más fuerte adsorbe la superficie (modelo de Hammer y Nørskov).

**Fórmulas.**

Columnas de los archivos (`qekit/modules/dos.py: read_dos_file`, `read_pdos_file`): del `<prefix>.dos` se toman $E$, DOS (una o dos columnas según espín) y la DOS integrada; del `pdos_atm#N(El)_wfc#M(l)` se toma la columna `ldos` (ya sumada sobre $m$), o `ldosup`/`ldosdw` con espín, según el número de columnas $1 + n_s(1 + (2l+1))$.

PDOS agregada (`dos.load`): $\rho_{\text{El},l}(E) = \sum_{\text{átomos } a \in \text{El}} \sum_{\text{wfc con } l} \text{ldos}_{a,l}(E)$; si dos.x y projwfc.x usan mallas de energía distintas, las proyecciones se interpolan linealmente (`np.interp`, cero fuera del rango) sobre la malla de la DOS total. Sin `<prefix>.dos`, la DOS total se define como $\sum_{\text{El},l}\rho_{\text{El},l}$.

DOS en el nivel de Fermi (`dos.report`): $\rho(E_F) = \sum_s \rho_s(E_{i^*})$ con $i^* = \arg\min_i |E_i - E_F|$; se llama "compatible con gap" si $\rho(E_F) < 10^{-3}$ estados/eV.

Momentos de una banda proyectada (`dos.momentos`), con $e = E - E_F$ y $\rho$ la PDOS del selector (por canal de espín, integrando con la regla del trapecio `np.trapezoid`):

$$
N = \int \rho(e)\,de, \qquad
\varepsilon_c = \frac{1}{N}\int e\,\rho(e)\,de, \qquad
W = \sqrt{\frac{1}{N}\int (e-\varepsilon_c)^2\rho(e)\,de}, \qquad
f = \frac{1}{N}\int_{e \le 0}\rho(e)\,de
$$

- $N$: estados integrados (adimensional); $\varepsilon_c$: centro (eV respecto a $E_F$); $W$: anchura rms (eV); $f$: llenado (fracción).
- Con dos canales, el valor reportado es el promedio de cada magnitud ponderado por $N_s$; el "desdoblamiento de intercambio" es $\varepsilon_c^{\uparrow} - \varepsilon_c^{\downarrow}$.
- Cola relativa: $\max(\rho)$ en los últimos $\max(3, n/50)$ puntos dividido por $\max(\rho)$ global; si supera 0.05 la banda está cortada por arriba.

**Cómo lo calcula Olla-DFT.**

1. `qekit/cli.py: _cmd_dos` → `dos.load(path, prefix)`.
2. `dos.load` intenta leer el XML con `qeout.read_xml` para tomar $E_F$ y el prefix; busca `<prefix>.dos` (o `*.dos`) y todos los `*pdos_atm#*`; parsea el nombre con la expresión `pdos_atm#(\d+)\(([A-Za-z]+)\)_wfc#(\d+)\(([A-Za-z])…\)` para obtener elemento y letra orbital (también con `p_j1.5` de SOC).
3. Si el XML no da $E_F$, lo toma del comentario `EFermi = …` de la cabecera del `.dos`.
4. Ordena las proyecciones por elemento (orden de aparición) y orbital $s,p,d,f$; `by_element` suma orbitales.
5. `dos.report` imprime rango de energías, $E_F$, origen del cero (`reference_energy`: Fermi salvo `--ref none`), canales, proyecciones y $\rho(E_F)$.
6. `dos.export` escribe `DOS.dat` (E, DOS[_up/_dw], DOS_integrada) y `PDOS.dat` (por elemento_orbital y totales por elemento); `dos.plot` / `dos.draw` dibujan total con relleno, PDOS por `--mode orbital|element|total`, espín ↓ reflejado hacia abajo.
7. Con `--dband El[-orb]`: `dos.momentos` exige $E_F$ (si no, `ErrorDeUso("no se encontró la energía de Fermi…")`) y la clave (El, orb); `report_momentos` imprime centro, anchura, llenado y avisos.

**De dónde sale cada dato.**

| Dato | Origen | Detalle |
|---|---|---|
| DOS total, integrada | `<prefix>.dos` de dos.x | columnas E, dos[, dosup, dosdw], int |
| PDOS por átomo/orbital | `<prefix>.pdos.pdos_atm#N(El)_wfc#M(l)` de projwfc.x | columna `ldos` (o `ldosup`, `ldosdw`) |
| Energía de Fermi | `<fermi_energy>` de `prefix.xml`; si no, `EFermi` en la cabecera del `.dos` | `dos.load` |
| Malla de energía | la de dos.x (`DeltaE = 0.02` eV en `gen`) | las PDOS se interpolan a ella |
| Umbral de "gap" | constante en `dos.report` | 1e-3 estados/eV |
| Umbral de cola | constante en `dos.momentos`/`report_momentos` | 5 % del pico |

**Límites y trampas.**

- Un archivo `pdos_atm#` con distinto número de puntos que el primero se salta, y `dos.load` lo anota en `DOSData.avisos`; el reporte lo imprime: "se han SALTADO N archivo(s) de projwfc.x cuya malla de energía no coincide con la del primero (… puntos), así que la PDOS está incompleta… Casi siempre es que hay archivos de dos corridas de projwfc.x mezclados en la misma carpeta". Los datos exportados siguen sin ese orbital: hay que apartar los archivos viejos y volver a cargar.
- La DOS total definida como suma de PDOS (cuando falta `.dos`) omite la parte de las funciones de onda que no cae en esferas atómicas: será menor que la real.
- `momentos` integra **todo** el rango disponible salvo `--dband-emax`; si la PDOS no ha decaído, avisa: "al final del rango todavía queda un X % del pico de PDOS. La banda está CORTADA por arriba, así que el centro sale más bajo de lo que debería". El texto recomienda "Vuelve a correr projwfc.x con un Emax mayor".
- El centro de banda d "es una correlación empírica dentro de una misma familia de metales, no una ley" (texto del reporte).
- Con SOC, la letra del orbital se toma de `p_j1.5` → `p`; las componentes $j$ se suman.
- `--ref vbm` no existe para la DOS: `reference_energy` solo distingue `none` del resto (siempre Fermi).

**Referencias.**

- P. E. Blöchl, O. Jepsen, O. K. Andersen, *Phys. Rev. B* **49**, 16223 (1994) — método de tetraedros (dos.x, `tetrahedra_opt`).
- B. Hammer, J. K. Nørskov, *Surf. Sci.* **343**, 211 (1995); *Adv. Catal.* **45**, 71 (2000) — modelo del centro de banda d.
- P. Giannozzi et al., *J. Phys.: Condens. Matter* **29**, 465901 (2017) — projwfc.x.

---

### `olla-dft plot` — figura combinada bandas + DOS

**Qué responde.** Produce la figura estándar de un artículo de estructura electrónica: bandas a la izquierda y DOS girada a la derecha, compartiendo el eje de energía y el mismo cero.

**Fundamento para no expertos.** Las bandas dicen *dónde* en el espacio k están los estados; la DOS dice *cuántos* hay a cada energía. Ponerlas lado a lado con el mismo cero (el VBM si hay gap, el Fermi si es metal) permite leer de un vistazo qué orbitales forman cada banda.

**Fórmulas.** Ninguna propia: `qekit/modules/combined.py: plot` solo dibuja; el cero se toma del análisis de bandas (`bands.reference_energy`) y se aplica a los dos paneles.

**Cómo lo calcula Olla-DFT.**

1. `qekit/cli.py: _cmd_plot` carga `bands.load` y `dos.load` sobre la misma carpeta.
2. Imprime `bands.gap_report`.
3. `combined.plot` crea dos ejes con proporción `ratio = 2.6`, llama a la lógica de dibujo de bandas y a `dos.draw(vertical=True)` con el mismo desplazamiento de energía.

**De dónde sale cada dato.**

| Dato | Origen | Detalle |
|---|---|---|
| Bandas y gap | `prefix.xml` (véase `bands`) | `bands.load` |
| DOS/PDOS | `.dos` y `pdos_atm#` (véase `dos`) | `dos.load` |
| Cero de energía | `bands.reference_energy(bs, ref)` | el mismo para ambos paneles |

**Límites y trampas.**

- Las bandas y la DOS suelen venir de cálculos distintos (camino vs. malla) con **el mismo scf**; si vienen de scf distintos, sus $E_F$ no coinciden y el panel derecho queda desplazado sin ningún aviso.
- El cero de la DOS es el de las bandas (VBM en `auto` con gap), aunque `olla-dft dos` por separado usaría el Fermi.

**Referencias.**

- Las de `bands` y `dos`.

---
### `olla-dft effmass` — masa efectiva por ajuste parabólico

**Qué responde.** Cuánto "pesa" un electrón en el fondo de la banda de conducción y un hueco en el techo de la de valencia: la masa efectiva $m^*/m_e$ en cada dirección, que gobierna movilidades, densidades de estados efectivas y niveles de excitones.

**Fundamento para no expertos.** Cerca de un extremo, una banda se parece a una parábola, igual que la energía cinética de una partícula libre $E = \hbar^2 k^2/2m$. La curvatura de esa parábola define una masa: una banda muy curvada (muy "abierta") corresponde a un portador ligero y rápido; una banda plana, a uno pesado. La masa puede depender de la dirección (en el silicio el electrón tiene una masa longitudinal de ~0.92 y una transversal de ~0.19), así que hay que ajustar parábolas a lo largo de rectas concretas en el espacio k.

Olla-DFT trabaja en dos etapas. Primero ajusta sobre las bandas que ya tienes (rápido, pero con pocos puntos y solo en las direcciones del camino). Después escribe un cálculo `bands` dedicado con líneas muy finas que cruzan el VBM y el CBM en tres direcciones (para un valle fuera de Γ: la radial Γ→k₀, "longitudinal", y dos perpendiculares, "transversales"; en Γ: [100], [110], [111]) y, cuando termina, ajusta una parábola por línea.

**Fórmulas.**

Ajuste cuadrático y masa (`qekit/modules/effmass.py: from_bands`, `collect_fine`, `_mass_from_quadratic`):

$$
E(k) \approx a\,k^2 + b\,k + c, \qquad \frac{m^*}{m_e} = \frac{\hbar^2/m_e}{2a}, \qquad \frac{\hbar^2}{m_e} = 7.6199682\ \text{eV·Å}^2
$$

- $k$: distancia al extremo a lo largo de la línea (Å⁻¹, con signo); $a$ en eV·Å²; el ajuste es `np.polyfit(x, y, 2)`.
- El signo se conserva: $a<0$ (curvatura hacia abajo) da $m^*<0$, que es la convención del reporte para un hueco.
- El término lineal $b$ se ajusta pero **no** se usa: el extremo se supone en $k=0$.

Calidad del ajuste (`effmass._r2`): $R^2 = 1 - \sum(y - \hat y)^2 / \sum (y - \bar y)^2$.

Direcciones del valle (`effmass.valley_directions`): si $|\mathbf k_0| < 10^{-6}$ Å⁻¹ (extremo en Γ), $\{[100], [110]/\sqrt2, [111]/\sqrt3\}$; si no, $\hat e_1 = \mathbf k_0/|\mathbf k_0|$ y dos perpendiculares construidas con productos vectoriales.

Puntos de la línea fina (`effmass.prepare`): $\mathbf k_j = \mathbf k_0 + t_j\,\hat e$, $t_j \in [-h, h]$ con `half_width = 0.06` Å⁻¹ y `npts = 21` (se fuerza impar), convertidos a fraccionarias con $\mathbf k_{\text{frac}} = \mathbf k\,\mathbf B^{-1}$.

Identificación de la banda de valencia en el cálculo fino (`collect_fine`): $n_v = \mathrm{round}(N_{\text{el}}/2) - 1$ (índice base 0) si `nspin = 1`; si no, la última banda cuyo máximo queda por debajo de $E_F$ (o del HOMO, o de la mediana).

**Cómo lo calcula Olla-DFT.**

1. `qekit/cli.py: _cmd_effmass` exige `--bands-dir` (una carpeta con un cálculo de bandas) salvo con `--collect`. Carga la estructura y `bands.load(bands_dir)`.
2. `effmass.from_bands`: `bands.analyze_gap` da VBM/CBM; para huecos toma la banda del VBM y las que en ese k están a menos de `DEGEN_TOL = 0.05` eV por debajo; para electrones la del CBM y las degeneradas por encima. Para cada banda relocaliza el extremo (`argmax`/`argmin`), delimita el tramo sin cruzar puntos especiales ni discontinuidades (`_segment_bounds`), reúne los puntos con $|k - k_0| \le$ `--window` (semiancho; por omisión `WINDOW_DEFAULT = PARABOLIC_MAX/2 = 0.06` Å⁻¹, ampliando hasta `--min-points` = 7) a ambos lados si el extremo es interior o a un lado si está en un punto especial (`_collect_window`), y ajusta.
3. Avisa si hay menos de 5 puntos ("solo N puntos: el ajuste no es confiable; haz el cálculo dedicado (effmass sin --collect y luego --collect)") o si el **tramo total** ajustado ($k_{\max} - k_{\min}$, `MassFit.window`) supera `PARABOLIC_MAX = 0.12` Å⁻¹ ("tramo ajustado de X Å⁻¹ (límite parabólico 0.12): el camino no tiene puntos más finos; haz el cálculo dedicado"). Con el semiancho por omisión un ajuste centrado en el extremo mide justo 0.12, así que el aviso solo salta cuando hubo que ensanchar por falta de puntos.
4. `effmass.prepare` reduce a la primitiva (`structure.primitive`), resuelve pseudos y cutoffs con `sweep.prepare_common`, escribe `masa.in` (`calculation='bands'`, `K_POINTS crystal` con las 6 líneas, `nbnd` igual al del cálculo de bandas, `occupations='fixed'` porque `insulator=True`) y `scf.in` (malla `sweep.default_grid`), y guarda `masa_meta.json` con la descripción de cada línea.
5. Con `--run`, `runner.run_all` lanza pw.x sobre `scf.in` y `masa.in`; con `--collect`, `effmass.collect_fine` lee `out/*.xml`, corta la lista de k en trozos de `npts`, calcula $t_j = \pm|\mathbf k_j - \mathbf k_c|$ y ajusta todas las bandas degeneradas con el extremo.
6. `effmass.report` imprime la tabla (portador, banda, $m^*/m_e$, $R^2$, puntos, Δk, dirección) y `export` escribe `MASA_EFECTIVA.dat`.

**De dónde sale cada dato.**

| Dato | Origen | Detalle |
|---|---|---|
| Autovalores y k cartesianos | `prefix.xml` de pw.x | `qeout.read_xml` (bandas previas y cálculo fino) |
| VBM, CBM y sus k | `bands.analyze_gap` | ver `olla-dft bands` |
| $\hbar^2/m_e$ | constante `effmass.HBAR2_OVER_ME` | 7.6199682 eV·Å² |
| Número de electrones | `<nelec>` del XML | para identificar la valencia en `collect_fine` |
| Ventana, mínimo de puntos, semiancho, puntos por línea | parámetros del usuario | `--window` (por omisión `effmass.WINDOW_DEFAULT` = 0.06 Å⁻¹), `--min-points 7`, `--half-width 0.06`, `--points 21` |
| Límite parabólico | `effmass.PARABOLIC_MAX` | 0.12 Å⁻¹ de tramo total (holgura `_TOL_VENTANA = 1e-6`) |
| Tolerancia de degeneración | `effmass.DEGEN_TOL` | 0.05 eV |
| Descripción de las líneas | `masa_meta.json` escrito por `prepare` | `effmass.load_meta` |

**Límites y trampas.**

- No hay un modo "solo ajuste rápido": con `--bands-dir` el comando hace siempre el ajuste sobre el camino y a continuación prepara el cálculo fino en `--outdir` (o lo corre con `--run`); con `--collect` lee el fino. Es lo que describe el docstring del módulo.
- En un metal el comando se detiene: "El sistema es metálico: no hay un extremo de banda aislado que ajustar".
- El reporte advierte que el cálculo **no incluye espín-órbita**: "cerca de Γ hay un triplete degenerado, no el par hueco pesado / hueco ligero del modelo de Luttinger".
- Un $R^2$ de 1.0000 con 3 o 4 puntos "no dice nada — una parábola pasa exacta por tres puntos cualesquiera" (texto del reporte).
- `--window` es un **semiancho**: el tramo ajustado mide hasta el doble, y es ese tramo el que se compara con `PARABOLIC_MAX`. Subir `--window` por encima de 0.06 Å⁻¹ dispara el aviso de régimen no parabólico aunque el ajuste sea "bueno" en $R^2$.
- `--collect` usa el **primer** XML de `out/`; si hay varios prefixes, puede leer el equivocado.
- El cálculo fino se escribe siempre con `occupations='fixed'` y sin espín; para sistemas magnéticos hay que editar `masa.in` a mano.
- Las líneas transversales para un valle fuera de Γ se eligen con un producto vectorial arbitrario: no son necesariamente ejes cristalográficos, y en un valle no elipsoidal las dos transversales pueden diferir.

**Referencias.**

- N. W. Ashcroft, N. D. Mermin, *Solid State Physics* (1976), cap. 12 — definición de masa efectiva.
- J. M. Luttinger, W. Kohn, *Phys. Rev.* **97**, 869 (1955) — modelo k·p de bandas de valencia degeneradas.
- Valores de referencia del silicio: M. Cardona, F. H. Pollak, *Phys. Rev.* **142**, 530 (1966).

---

### `olla-dft fermi` — superficie de Fermi en formato BXSF

**Qué responde.** Qué bandas cruzan el nivel de Fermi y cómo es la superficie $\varepsilon_n(\mathbf k) = E_F$ de cada una, escrita en un archivo BXSF que XCrySDen o FermiSurfer dibujan en 3D.

**Fundamento para no expertos.** En un metal los estados se llenan hasta una energía $E_F$; el conjunto de puntos k cuya energía es exactamente $E_F$ forma una superficie en la zona de Brillouin, la superficie de Fermi. Su forma determina la conductividad, las oscilaciones cuánticas y muchas inestabilidades (ondas de densidad de carga, superconductividad). Para dibujarla hace falta $\varepsilon_n(\mathbf k)$ en una malla **completa y uniforme** de la zona de Brillouin, que es lo que produce `olla-dft transport` (un nscf con `nosym`, `noinv`).

**Fórmulas.**

Bandas que cruzan $E_F$ (`qekit/modules/transport.py: crossing_bands`), con `tol = 1e-6` eV:

$$
\min_{\mathbf k}\varepsilon_n(\mathbf k) < E_F - \delta \quad\wedge\quad \max_{\mathbf k}\varepsilon_n(\mathbf k) > E_F + \delta
$$

Reconstrucción de la malla (`transport.load`): las fraccionarias se llevan a $[0,1)$ con $f \leftarrow f - \lfloor f + 10^{-6}\rfloor$, se redondean a 6 decimales y $n_i$ es el número de valores distintos en cada eje; se exige $n_1 n_2 n_3 = N_k$.

Rejilla BXSF (`transport.export_bxsf`): se escriben $(n_i + 1)$ puntos por eje repitiendo el primer plano al final (`np.pad(..., mode="wrap")`), en orden C (último índice más rápido), con los vectores recíprocos $\mathbf b_i = 2\pi(\mathbf A^{-1})^{\mathsf T}_i$ en Å⁻¹ y las energías en eV.

**Cómo lo calcula Olla-DFT.**

1. `qekit/cli.py: _cmd_fermi` busca `out/*.xml` dentro de `--outdir` (por omisión `transporte`).
2. `transport.load` lee el XML con `qeout.read_xml`; rechaza un XML de tipo `scf` ("es de un cálculo SCF, no del nscf de malla densa"); reconstruye la malla y reordena las energías con `np.lexsort`. Calcula además velocidades de banda por diferencias finitas (no se usan aquí) y avisa si la malla es menor que 24×24×24 o tiene menos de 12 000 puntos.
3. `transport.crossing_bands` lista las bandas metálicas; si no hay, imprime "Ninguna banda cruza E_F: el sistema no es metálico y no tiene superficie de Fermi".
4. `transport.export_bxsf` escribe `superficie_fermi.bxsf` con `Fermi Energy`, la rejilla y un bloque `BAND:` por banda.

**De dónde sale cada dato.**

| Dato | Origen | Detalle |
|---|---|---|
| Autovalores en la malla | `prefix.xml` del nscf de `olla-dft transport` | `qeout.read_xml`; solo el canal de espín 0 |
| Energía de Fermi | `<fermi_energy>` del XML | `run.fermi`; sin ella, `ErrorDeUso("no hay nivel de Fermi…")` |
| Celda | `<atomic_structure>` del XML | vectores recíprocos con $2\pi$ |
| Tolerancia de cruce | argumento `tol` de `crossing_bands` | 1e-6 eV |

**Límites y trampas.**

- Solo funciona sobre la carpeta de `olla-dft transport` (mismo nscf de malla completa); un camino de bandas o una malla reducida por simetría falla con "los N puntos k no forman una malla uniforme".
- Solo se exporta el canal de espín 0 (`transport.load(spin=0)`); un metal ferromagnético necesitaría dos archivos y el comando no los produce.
- $E_F$ se toma tal cual del XML del nscf; con `occupations='fixed'` no existe y el comando falla.
- El nivel de Fermi de un nscf denso no se recalcula: es el que heredó del scf (malla más gruesa).
- Los vectores recíprocos se escriben en Å⁻¹ con el factor $2\pi$; el visualizador debe interpretarlos en esas unidades.

**Referencias.**

- A. Kokalj, *Comput. Mater. Sci.* **28**, 155 (2003) — XCrySDen y el formato BXSF. DOI 10.1016/S0927-0256(03)00104-6.
- M. Kawamura, *Comput. Phys. Commun.* **239**, 197 (2019) — FermiSurfer. DOI 10.1016/j.cpc.2019.01.017.

---

### `olla-dft unfold` — desdoblamiento de bandas de una supercelda

**Qué responde.** Qué parte de cada estado de la supercelda "pertenece" a cada punto k de la celda primitiva: el peso espectral que permite ver la banda del material original (y cuánto la difumina un defecto, un dopante o el desorden) a partir de un cálculo en supercelda.

**Fundamento para no expertos.** Una supercelda de $N$ celdas primitivas tiene una zona de Brillouin $N$ veces menor, así que sus bandas salen *plegadas*: donde la primitiva tenía una banda, hay $N$ ramas amontonadas. Cada estado de la supercelda es una suma de ondas planas $e^{i(\mathbf K + \mathbf G)\cdot\mathbf r}$, y cada onda plana tiene un vector de onda bien definido. Preguntar "¿cuánto de este estado vive en el punto $\mathbf k$ de la primitiva?" tiene respuesta exacta: la suma de $|C(\mathbf G)|^2$ sobre las ondas planas cuyo $\mathbf K + \mathbf G$ coincide con $\mathbf k$ módulo la red recíproca primitiva. Si la supercelda es perfecta, cada estado tiene peso 1 en un solo $\mathbf k$ y se recupera la banda primitiva; si hay un defecto, el peso se reparte y la banda se ve borrosa. Esa borrosidad es el resultado físico.

**Fórmulas.**

Matriz de supercelda (`qekit/modules/unfold.py: matriz_supercelda`): $\mathbf M = \mathbf A_{\text{sc}}\,\mathbf a_{\text{prim}}^{-1}$, redondeada a enteros; se acepta si $\max|\mathbf M - \mathrm{round}(\mathbf M)| \le 10^{-3}$. Si falla por orientación, `_m_por_metricas` busca $\mathbf M$ entera tal que $\mathbf G_{\text{sc}} = \mathbf M\,\mathbf G_p\,\mathbf M^{\mathsf T}$ con $\mathbf G = \mathbf X\mathbf X^{\mathsf T}$ el tensor métrico (invariante bajo rotaciones), fila a fila entre vectores enteros con la longitud correcta. Después la primitiva se **rederiva** como $\mathbf a = \mathbf M^{-1}\mathbf A_{\text{sc}}$ para que ambas compartan ejes. $N = |\det\mathbf M|$.

Coordenadas: como $\mathbf B_{\text{sc}} = \mathbf M^{-\mathsf T}\mathbf b_{\text{prim}}$, un vector con coordenadas $\mathbf c_{\text{sc}}$ en la base recíproca de la supercelda tiene coordenadas $\mathbf c_p = \mathbf c_{\text{sc}}\mathbf M^{-\mathsf T}$ en la primitiva, y $\mathbf k_{\text{prim}} = \mathbf k_{\text{sc}}\mathbf M^{-\mathsf T}$ (`desdoblar`).

Peso espectral (`unfold.pesos_de_k`), con $\mathbf m_0 = \mathbf k_{\text{prim}}\mathbf M^{\mathsf T} - \mathbf k_{\text{sc}}$ (debe ser entero a `TOL_ENTERO = 1e-4`; si no, el peso es 0 porque ese $\mathbf k$ no se pliega sobre este $\mathbf K$):

$$
P_{n}(\mathbf k) = \frac{\sum_{\mathbf G \in S}\ \sum_{\sigma}|C_{n\sigma}(\mathbf G)|^2}{\sum_{\mathbf G}\ \sum_{\sigma}|C_{n\sigma}(\mathbf G)|^2}, \qquad
S = \left\{\mathbf G : (\mathbf G - \mathbf m_0)\,\mathbf M^{-\mathsf T} \in \mathbb{Z}^3\right\}
$$

- $C_{n\sigma}(\mathbf G)$: coeficientes de ondas planas de la banda $n$ (componente de espinor $\sigma$ si `npol = 2`) leídos de `wfc<N>.dat`; $\mathbf G$ dado por sus índices de Miller en la base recíproca de la supercelda.
- El denominador normaliza por si los coeficientes no están normalizados; $P_n \in [0,1]$.

Distancia en el eje x (`unfold._distancias`): suma de $|\Delta\mathbf k|$ con $\mathbf k = \mathbf k_{\text{frac}}\,\mathbf b_{\text{prim}}$; un salto mayor que 5 veces la mediana de los pasos no nulos se cuenta como cero (cambio de rama).

**Cómo lo calcula Olla-DFT.**

1. `qekit/cli.py: _cmd_unfold` carga la estructura primitiva y llama a `unfold.desdoblar(path, celda_primitiva, bandas=range(--bands), spin=--spin)`; `--spin` es `up` (por omisión) o `dw`.
2. `desdoblar` lee el XML (`qeout.read_xml`), localiza la carpeta `.save` (`_carpeta_save`) y los archivos de función de onda del canal pedido (`qekit/core/wfc.py: buscar_wfc(save, spin)`, ordenados por número de k): si `wfc.es_lsda` detecta `wfcup*`/`wfcdw*`, devuelve solo los `wfc{up|dw}<N>.dat` de ese canal; si no, los `wfc<N>.dat` de un cálculo sin espín. Sin ellos: "El cálculo no guardó las funciones de onda: eso pasa con disk_io='nowf' o 'low'" (o, con lsda, "falta el canal '…'").
3. `matriz_supercelda` obtiene $\mathbf M$; los k del cálculo se pasan a coordenadas primitivas.
4. Para cada punto k, `wfc.leer_wfc` lee el archivo Fortran sin formato: registro 1 (`ik`, `xk`, `ispin`, `gamma_only`, `scalef`), registro 2 (`ngw`, `igwx`, `npol`, `nbnd`), registro 3 (`b1,b2,b3`), registro 4 (índices de Miller) y un registro por banda con `npol·igwx` complejos (solo se materializan las bandas pedidas).
5. `pesos_de_k` calcula $P_n(\mathbf k)$; las energías salen del XML, del mismo canal de espín que las funciones de onda (`res.eigenvalues[0]` para `up`, `[1]` para `dw`).
6. `unfold.report` imprime $N$, $\mathbf M$, la distribución del peso (media, fracción > 0.9, fracción < 0.1) y avisos; `export` escribe `UNFOLD.dat` (distancia, $E - E_F$, peso) y `UNFOLD.txt`; `plot` dibuja un `scatter` con tamaño $= 60\,P$ para pesos > 0.005.

**De dónde sale cada dato.**

| Dato | Origen | Detalle |
|---|---|---|
| Coeficientes $C(\mathbf G)$, índices de Miller, `npol` | `out/<prefix>.save/wfc<N>.dat` de pw.x | `wfc.leer_wfc` (formato Fortran secuencial, little-endian) |
| Autovalores, k fraccionarios, $E_F$, celda de la supercelda | `prefix.xml` | `qeout.read_xml` |
| Canal de espín | `--spin up|dw` (usuario) | `wfc.buscar_wfc`, `wfc.es_lsda`; en un cálculo sin espín no cambia nada |
| Celda primitiva | archivo del usuario | se rederiva como $\mathbf M^{-1}\mathbf A_{\text{sc}}$ |
| Tolerancia de enteros | `unfold.TOL_ENTERO` | 1e-4 (y 1e-3 para aceptar $\mathbf M$) |

**Límites y trampas.**

- Hace falta `disk_io='medium'` o `'high'` en el cálculo de bandas de la supercelda; `olla-dft gen` no lo pone por defecto.
- Se desdobla **un solo canal de espín por corrida**. En un cálculo `lsda` el reporte avisa: "el cálculo es de espín polarizado (lsda) y aquí solo se ha desdoblado el canal 'up' (wfcup<N>.dat y sus energías). El otro canal no se mezcla ni se suma: para verlo repite el desdoblamiento con --spin dw". Los dos canales nunca se combinan en una sola figura.
- Si la supercelda está relajada y la primitiva no (o al revés), $\mathbf M$ no sale entera: "la celda de la supercelda no es un múltiplo entero de la primitiva (error …)".
- Si casi todos los pesos valen 1, el reporte avisa: "la supercelda parece PERFECTA (sin defecto ni desorden). En ese caso el desdoblamiento reproduce exactamente las bandas primitivas — que es la comprobación de que funciona, pero no un resultado nuevo".
- Los k del cálculo se interpretan como k de la supercelda y se convierten a la primitiva; no se genera un camino primitivo ni se comprueba que los k de la supercelda sean los pliegues correctos del camino deseado.
- No se tienen en cuenta las funciones de onda de pseudopotenciales ultrasuaves/PAW más allá de la parte de ondas planas (no se incluye el término de aumento $S$); el peso es el de la parte suave.

**Referencias.**

- V. Popescu, A. Zunger, *Phys. Rev. B* **85**, 085201 (2012) — peso espectral de desdoblamiento. DOI 10.1103/PhysRevB.85.085201.
- P. V. C. Medeiros, S. Stafström, J. Björk, *Phys. Rev. B* **89**, 041407(R) (2014) — desdoblamiento con ondas planas (BandUP). DOI 10.1103/PhysRevB.89.041407.
- W. Ku, T. Berlijn, C.-C. Lee, *Phys. Rev. Lett.* **104**, 216401 (2010).

---
### `olla-dft wannier` — funciones de Wannier e interpolación de bandas

**Qué responde.** Construye, a partir de un cálculo DFT en una malla gruesa de puntos k, un modelo pequeño $H_{mn}(\mathbf R)$ en una base de funciones localizadas (funciones de Wannier) con el que las bandas se pueden evaluar en **cualquier** punto k sin volver a correr pw.x; da además dónde está centrada cada función, cuánto se extiende (su dispersión $\Omega$), y cuánto se parece la banda interpolada a la de DFT.

**Fundamento para no expertos.** Los estados de Bloch $\psi_{n\mathbf k}$ están deslocalizados por todo el cristal. Su transformada de Fourier en k da funciones $|\mathbf R n\rangle$ localizadas alrededor de una celda $\mathbf R$: las funciones de Wannier. En esa base el hamiltoniano es una matriz pequeña $H(\mathbf R)$ que decae con $|\mathbf R|$, y volver a transformar a k da la banda en cualquier punto (una "interpolación de Fourier" que es exacta en los puntos de partida). La dificultad es que cada $\psi_{n\mathbf k}$ está definida salvo una fase (y, si hay bandas degeneradas, salvo una rotación unitaria entre ellas): esa libertad se llama *gauge*. Con un gauge arbitrario las funciones de Wannier no están localizadas y la interpolación es basura. Marzari y Vanderbilt propusieron elegir el gauge que minimiza la dispersión total $\Omega$ (la suma de las "anchuras" cuadráticas); un buen punto de partida es proyectar sobre orbitales atómicos de prueba y ortonormalizar.

Cuando las bandas que interesan se cruzan con otras (metales, bandas de conducción) no existe un grupo aislado que transformar: hay que elegir en cada k un subespacio de $J$ estados que "se conecte suavemente" con el de sus vecinos. Eso es el *desenredado* de Souza, Marzari y Vanderbilt, con una ventana *exterior* (de dónde se puede elegir) y opcionalmente una *congelada* (estados que deben conservarse exactos). Olla-DFT implementa las dos cosas en Python, usando solo los solapes y proyecciones que calcula `pw2wannier90.x` (incluido con QE), sin necesitar wannier90; si el usuario tiene wannier90, también lee su `seedname_hr.dat`.

**Fórmulas.**

Malla completa y vectores $\mathbf b$ (`qekit/modules/wannier.py: malla_completa`, `capas_b`, `residuo_completitud`): $\mathbf k_{ijk} = (i/n_1, j/n_2, k/n_3)$ en el orden de QE (último índice más rápido). Las capas de vecinos $\mathbf b = (h_1/n_1, h_2/n_2, h_3/n_3)\,\mathbf B$ se añaden por distancia hasta cumplir, por mínimos cuadrados sobre las 6 componentes independientes,

$$
\sum_{\mathbf b} w_{\mathbf b}\, b_\alpha b_\beta = \delta_{\alpha\beta}, \qquad \text{residuo } = \left\|\textstyle\sum_{\mathbf b} w_{\mathbf b}\,\mathbf b\otimes\mathbf b - \mathbf 1\right\|_\infty < 10^{-5}
$$

- $w_{\mathbf b}$: peso de cada capa (Å²); capas que no añaden rango (SVD) o con $|w| < 10^{-8}$ se descartan.

Gauge de proyección (`wannier.gauge_proyeccion`), con $A_{mn}(\mathbf k) = \langle\psi_{m\mathbf k}|g_n\rangle$ del `.amn` y la SVD $A = u\,s\,v^\dagger$:

$$
U(\mathbf k) = A\,(A^\dagger A)^{-1/2} = u\,v^\dagger
$$

- $U$: matriz $N_b\times J$ con columnas ortonormales (Löwdin); se reporta el menor valor singular $s_{\min}$ (aviso si $< 0.2$).

Dispersión invariante y desenredado (`wannier.omega_I`, `gauge_desenredo`), con $M^{\mathbf k,\mathbf b}_{mn} = \langle u_{m\mathbf k}|u_{n,\mathbf k+\mathbf b}\rangle$ del `.mmn`:

$$
\Omega_I = \frac{1}{N_k}\sum_{\mathbf k}\sum_{\mathbf b} w_{\mathbf b}\left[J - \left\|U^\dagger(\mathbf k)\,M^{\mathbf k,\mathbf b}\,U(\mathbf k+\mathbf b)\right\|_F^2\right]
$$

$$
Z(\mathbf k) = \sum_{\mathbf b} w_{\mathbf b}\, M^{\mathbf k,\mathbf b}\,U(\mathbf k+\mathbf b)\,U^\dagger(\mathbf k+\mathbf b)\,M^{\mathbf k,\mathbf b\,\dagger}
$$

- En cada iteración $Z$ se mezcla con la anterior ($Z \leftarrow \mu Z_{\text{new}} + (1-\mu)Z_{\text{old}}$, $\mu = 0.5$ inicial, se reduce a la mitad si $\Omega_I$ sube), se restringe a las bandas de la ventana exterior, se proyectan fuera las congeladas ($Q Z Q$ con $Q = 1 - P_{\text{frozen}}$) y se toman los $J - N_{\text{frozen}}$ autovectores de mayor autovalor. Máximo 200 pasos, tolerancia $10^{-10}$ Å². Al final se reproyecta sobre los orbitales de prueba dentro del subespacio ($U \leftarrow U\,\mathrm{polar}(U^\dagger A)$) para tener un gauge de partida suave.

Hamiltoniano en espacio real e interpolación (`wannier.hamiltoniano_k`, `a_reales`, `interpolar`, `celda_wigner_seitz`):

$$
H(\mathbf k) = U^\dagger(\mathbf k)\,\mathrm{diag}\big(\varepsilon_n(\mathbf k)\big)\,U(\mathbf k), \qquad
H(\mathbf R) = \frac{1}{N_k}\sum_{\mathbf k} e^{-2\pi i\,\mathbf k\cdot\mathbf R}\,H(\mathbf k), \qquad
H^{\text{int}}(\mathbf k) = \sum_{\mathbf R}\frac{e^{2\pi i\,\mathbf k\cdot\mathbf R}}{\deg(\mathbf R)}\,H(\mathbf R)
$$

- $\mathbf k$ y $\mathbf R$ en coordenadas fraccionarias (por eso el $2\pi$ explícito). $\mathbf R$ recorre los vectores de la celda de Wigner-Seitz de la superred $n_1\times n_2\times n_3$; $\deg(\mathbf R)$ es el número de imágenes equidistantes (tolerancia $10^{-5}$ Å²). Las bandas son los autovalores de $\tfrac12(H^{\text{int}} + H^{\text{int}\dagger})$.

Centros y dispersión (`wannier.dispersion`), ecuaciones 31 y 34–36 de Marzari-Vanderbilt, con $M^W = U^\dagger(\mathbf k) M^{\mathbf k,\mathbf b} U(\mathbf k+\mathbf b)$ y $\phi_n = \operatorname{Im}\ln M^W_{nn}$:

$$
\bar{\mathbf r}_n = -\frac{1}{N_k}\sum_{\mathbf k,\mathbf b} w_{\mathbf b}\,\mathbf b\,\phi_n, \qquad
\Omega_n = \frac{1}{N_k}\sum_{\mathbf k,\mathbf b} w_{\mathbf b}\left[\left(1-|M^W_{nn}|^2\right) + \phi_n^2\right] - |\bar{\mathbf r}_n|^2
$$

$$
\Omega_I = \frac{1}{N_k}\sum_{\mathbf k,\mathbf b} w_{\mathbf b}\Big[J - \sum_{mn}|M^W_{mn}|^2\Big], \quad
\Omega_{OD} = \frac{1}{N_k}\sum_{\mathbf k,\mathbf b} w_{\mathbf b}\sum_{m\ne n}|M^W_{mn}|^2, \quad
\Omega_D = \frac{1}{N_k}\sum_{\mathbf k,\mathbf b} w_{\mathbf b}\sum_n\left(\phi_n + \mathbf b\cdot\bar{\mathbf r}_n\right)^2
$$

- $\bar{\mathbf r}_n$ en Å, $\Omega$ en Å²; $\Omega = \sum_n\Omega_n = \Omega_I + \Omega_D + \Omega_{OD}$ (el reporte imprime la suma para comprobarlo).

Minimización (`wannier._gradiente`, `_rotar`, `minimizar`), ec. 52–57 de Marzari-Vanderbilt, con $R_{mn} = M_{mn}M^*_{nn}$, $T_{mn} = (M_{mn}/M_{nn})\,q_n$, $q_n = \phi_n + \mathbf b\cdot\bar{\mathbf r}_n$, $\mathcal A(B) = (B - B^\dagger)/2$, $\mathcal S(B) = (B + B^\dagger)/2i$:

$$
G(\mathbf k) = -\frac{4}{N_k}\sum_{\mathbf b} w_{\mathbf b}\left[\mathcal A(R^{\mathbf k,\mathbf b}) - \mathcal S(T^{\mathbf k,\mathbf b})\right], \qquad
U(\mathbf k) \leftarrow U(\mathbf k)\,\exp\!\left(-\Delta t\,G(\mathbf k)\right), \qquad
\Delta t_0 = \frac{\alpha}{4\sum_{\mathbf b} w_{\mathbf b}},\ \alpha = 2
$$

- Si el paso sube $\Omega$ se parte por la mitad hasta 12 veces; máximo 500 pasos (`--iterations`); parada cuando $|\Delta\Omega| < 10^{-10}$. Se comprueba que $\Omega_I$ no cambie (`deriva_I`).

DOS interpolada (`wannier.dos_interpolada`): $\rho(E) = \frac{1}{N_k\,\sigma\sqrt{2\pi}}\sum_{\mathbf k,n} e^{-(E-\varepsilon_n(\mathbf k))^2/2\sigma^2}$ en una malla $N^3$ (`--dos N`), $\sigma$ = `--sigma` 0.05 eV; integra a $J$ estados por celda, **sin** factor 2 de espín. La cabecera de `WANNIER_dos.dat` lo declara con `wannier.DOS_UNIDADES`: "estados/eV/celda, sin factor de espín: integra a num_wann (x2 para comparar con dos.x sin espín)".

**Cómo lo calcula Olla-DFT.**

1. *Preparar* (`qekit/cli.py: _cmd_wannier` → `wannier.prepare`): traduce `--projections` (`Si:sp3`, `O:p;Ti:d`, `f=0.25,0.25,0.25:s`, o `auto` = $s$ y $p$ en cada átomo) a orbitales $(l, m_r)$ de la tabla `ORBITALES` (convención de wannier90); escribe `1_scf.in`, `2_nscf.in` (malla completa `--grid` 4×4×4 por omisión, `K_POINTS crystal`, `nosym`, `noinv`, `conv_thr 1e-10`, `nbnd = --bands` o $J$ + excluidas), el `<prefix>.nnkp` (`escribir_nnkp`: red real y recíproca, k-points, proyecciones, vecinos `nnkpts` con sus $\mathbf G$, `exclude_bands`), `3_pw2wan.in` (`write_amn`, `write_mmn`, `write_unk=.false.`), `<prefix>.win` (por si se prefiere wannier90) y `4_bands.in` (bandas DFT sobre el camino de seekpath, 30 puntos por tramo, `outdir='./out_bandas'`).
2. *Correr* (`--run` → `wannier.correr`): `pw.x` (scf), `pw.x` (nscf), `pw2wannier90.x`, `pw.x` (bandas), en ese orden, deteniéndose en el primer fallo.
3. *Recoger* (`--collect` → `wannier.collect`): lee `.eig` (`leer_eig`), `.amn` (`leer_amn`), `.mmn` (`leer_mmn`, con $m$ corriendo más rápido → `reshape(order="F")`); recalcula capas y vecinos desde el `.nnkp` (`_leer_nnkp`); si hay más bandas que funciones o se dieron `--window`/`--frozen`, `gauge_desenredo`; si no, `gauge_proyeccion`. `dispersion` antes y después de `minimizar` (salvo `--no-minimize`). `celda_wigner_seitz`, `hamiltoniano_k`, `a_reales`; comprueba que `interpolar` reproduce la malla (`error_malla` < `TOL_EXACTA = 1e-6` eV) y, si hay bandas DFT (`out_bandas`, `--dft-bands`), compara en puntos que no estaban en la malla. Como control negativo, repite la interpolación con $U = 1$ (`E_sin_gauge`).
4. Si existe un `*_hr.dat` de wannier90 (distinto del propio `WANNIER_hr.dat`), `leer_hr` lo usa directamente y se salta la localización.
5. `wannier.report` imprime malla, vecinos y residuo, ventanas, $\Omega_I$, $\Omega$ troceada, centros con asignación a átomo o enlace (`asignar`, ventana de enlace 0.5–3.2 Å), decaimiento de $H(\mathbf R)$, exactitud en la malla y error frente a DFT; `export` escribe `WANNIER_hr.dat` (formato de wannier90), `WANNIER_centros.dat`, `WANNIER_bandas.dat`, `WANNIER.txt` y opcionalmente `WANNIER_dos.dat`; `plot` dibuja Wannier sobre DFT y la traza de $\Omega$.

**De dónde sale cada dato.**

| Dato | Origen | Detalle |
|---|---|---|
| Energías $\varepsilon_n(\mathbf k)$ | `seedname.eig` de pw2wannier90.x | eV absolutos; `wannier.leer_eig` |
| Proyecciones $A_{mn}(\mathbf k)$ | `seedname.amn` | `wannier.leer_amn` |
| Solapes $M^{\mathbf k,\mathbf b}_{mn}$ | `seedname.mmn` | `wannier.leer_mmn` |
| Celda, malla, bandas excluidas | `seedname.nnkp` (escrito por Olla-DFT) | `wannier._leer_nnkp` |
| $H(\mathbf R)$ externo | `seedname_hr.dat` de wannier90 | `wannier.leer_hr` |
| Bandas DFT de validación | `out_bandas/*.xml` del paso 4 | `qeout.read_xml`, canal 0 |
| Orbitales de prueba $(l, m_r)$ | tabla `wannier.ORBITALES` | Tabla 3.1/3.2 del manual de wannier90 |
| Tolerancias | `TOL_COMPLETITUD 1e-5`, `TOL_PESO 1e-8`, `TOL_EXACTA 1e-6` eV | constantes del módulo |
| Camino de alta simetría | seekpath (`wannier.camino_denso`) | 30 puntos por tramo (`--points`) |

**Límites y trampas.**

- Las ventanas `--window` y `--frozen` se comparan con las energías **absolutas** del `.eig` (no relativas a $E_F$), como en wannier90.
- Con desenredado y sin ventana congelada, nada tiene que reproducirse exacto; el reporte avisa: "Sin ventana congelada no hay ninguna banda que la interpolación tenga que reproducir exactamente… Si quieres que la valencia salga exacta, pásala en --frozen".
- Las proyecciones `auto` ($s$ y $p$ por átomo) fallan en metales de transición (faltan las $d$) y en enlaces muy covalentes; el reporte lo dice siempre.
- Si $H(\mathbf R)$ en el borde de la superred supera el 5 % de $H(0)$: "H(R) apenas ha decaído al borde de la superred: la base no está localizada".
- Sin `--dft-bands` (o `out_bandas`) el reporte avisa: "No has comparado con bandas de DFT. Que la interpolación reproduzca la malla es trivial".
- Solo se lee el canal de espín 0 de las bandas DFT; el flujo no está pensado para `nspin = 2` ni SOC (pw2wannier90 sí los soporta, pero `prepare` no escribe `nspin`).
- La minimización es descenso por gradiente con búsqueda de línea, no el gradiente conjugado de wannier90: puede necesitar más pasos y puede parar en un mínimo local.
- Mallas muy anisótropas pueden no admitir capas que cumplan la completitud: "no encuentro un conjunto de capas de vecinos que cumpla la condición de completitud con esta malla".
- `--collect` sin la estructura como primer argumento falla: "para analizar hace falta la estructura".
- La DOS interpolada no lleva el factor 2 de espín (la cabecera del archivo lo dice y pide multiplicar por 2 para comparar con dos.x sin espín) y solo vale dentro del rango de energías cubierto por las funciones de Wannier.

**Referencias.**

- N. Marzari, D. Vanderbilt, *Phys. Rev. B* **56**, 12847 (1997) — funciones de Wannier maximalmente localizadas. DOI 10.1103/PhysRevB.56.12847.
- I. Souza, N. Marzari, D. Vanderbilt, *Phys. Rev. B* **65**, 035109 (2001) — desenredado. DOI 10.1103/PhysRevB.65.035109.
- N. Marzari, A. A. Mostofi, J. R. Yates, I. Souza, D. Vanderbilt, *Rev. Mod. Phys.* **84**, 1419 (2012) — revisión. DOI 10.1103/RevModPhys.84.1419.
- G. Pizzi et al., *J. Phys.: Condens. Matter* **32**, 165902 (2020) — Wannier90 v3 (formatos `.nnkp`, `.amn`, `.mmn`, `_hr.dat`). DOI 10.1088/1361-648X/ab51ff.
- P.-O. Löwdin, *J. Chem. Phys.* **18**, 365 (1950) — ortonormalización simétrica.

---

### `olla-dft topology` — número de Chern y lazos de Wilson

**Qué responde.** Si el subespacio ocupado de un modelo de Wannier, en una sección bidimensional de la zona de Brillouin, tiene número de Chern distinto de cero (un invariante topológico entero) y cómo evolucionan los centros híbridos de Wannier (lazos de Wilson) a lo largo de esa sección.

**Fundamento para no expertos.** Además de sus energías, las bandas tienen una "geometría": al recorrer un lazo cerrado en el espacio k, los estados ocupados acumulan una fase (fase de Berry) que no depende de cómo se elijan las fases de cada estado. Sumando esa fase sobre toda una sección 2D de la zona de Brillouin sale un número entero, el número de Chern, que no cambia con deformaciones suaves del sistema: es *topológico*. Un Chern no nulo implica corrientes de borde sin disipación (efecto Hall cuántico anómalo). El *lazo de Wilson* es la versión "por rebanadas": para cada $k_2$ se calcula el producto de solapes a lo largo de $k_1$; las fases de sus autovalores son las posiciones (módulo 1) de las funciones de Wannier híbridas, y su "bombeo" al variar $k_2$ es otra forma de ver el Chern.

**Fórmulas.**

Malla y estados (`qekit/modules/topology.py: kmesh`, `analyze`): $\mathbf k_{ij}$ con $k_a = i/n_1$, $k_b = j/n_2$ y la tercera coordenada fija en `--fixed` (mod 1), en el plano `--plane` (`xy`, `xz`, `yz`); los autovectores $|u_n(\mathbf k)\rangle$ salen de `wannier.interpolar(..., vectores=True)`.

Enlaces unitarios y curvatura de Berry discreta (`topology._unitary_overlap`, `invariants_from_vectors`), con $V(\mathbf k)$ la matriz $N_w\times N_{\text{occ}}$ de autovectores ocupados:

$$
O_\mu(\mathbf k) = V^\dagger(\mathbf k)\,V(\mathbf k+\hat\mu), \qquad
Q_\mu = u\,v^\dagger \ \text{(parte unitaria de } O_\mu = u\,s\,v^\dagger), \qquad
U_\mu(\mathbf k) = \frac{\det Q_\mu(\mathbf k)}{|\det Q_\mu(\mathbf k)|}
$$

$$
F_{12}(\mathbf k) = \arg\!\left[U_1(\mathbf k)\,U_2(\mathbf k+\hat 1)\,U_1^*(\mathbf k+\hat 2)\,U_2^*(\mathbf k)\right], \qquad
C = \frac{1}{2\pi}\sum_{\mathbf k} F_{12}(\mathbf k)
$$

- $\hat\mu$: paso de malla en la dirección $\mu$ (periódico). $F_{12} \in (-\pi, \pi]$ por plaqueta (rad); $C$ se redondea al entero más próximo y se reporta el residuo $|C - \mathrm{round}(C)|$.
- Se reporta también el menor valor singular de todos los $O_\mu$ (`min_overlap`); si $< 10^{-6}$, aviso de malla demasiado gruesa.

Lazos de Wilson (`invariants_from_vectors`):

$$
W(k_2) = \prod_{i=0}^{n_1-1} Q_1(k_1^{(i)}, k_2), \qquad
x_n(k_2) = \frac{\arg\lambda_n\!\left[W(k_2)\right]}{2\pi} \bmod 1
$$

- $x_n$: centros híbridos de Wannier ordenados, en fracciones del vector de red a lo largo de la dirección 1.

Gaps de la sección: $E_g^{\text{dir}} = \min_{\mathbf k}[\varepsilon_{N_{\text{occ}}+1} - \varepsilon_{N_{\text{occ}}}]$, $E_g^{\text{ind}} = \min_{\mathbf k}\varepsilon_{N_{\text{occ}}+1} - \max_{\mathbf k}\varepsilon_{N_{\text{occ}}}$. Se exige $E_g^{\text{dir}} > $ `--gap-tol` (1e-8 eV).

**Cómo lo calcula Olla-DFT.**

1. `qekit/cli.py: _cmd_topology` exige exactamente una de `--occupied N` o `--fermi EV`, y `--grid` de al menos 3×3 (40×40 por omisión).
2. `topology.resolve_model` acepta un `*_hr.dat` o una carpeta con `WANNIER_hr.dat` (o un único `*_hr.dat`; si hay varios, error "indica el archivo exacto").
3. `wannier.leer_hr` lee $H(\mathbf R)$, $\mathbf R$ y degeneraciones; `wannier.interpolar` diagonaliza en la malla de la sección.
4. Con `--fermi`, cuenta los estados con $\varepsilon < E_F$ en cada k; si el número varía, error: "el nivel de Fermi corta bandas… El sistema es metálico en esta sección y el Chern de 'las ocupadas' no está definido".
5. `invariants_from_vectors` calcula curvatura, Chern y lazos de Wilson.
6. `topology.report` imprime gaps, Chern discreto y entero, residuo y solape mínimo; `export` escribe `TOPOLOGY_curvature.dat` (flujo por plaqueta), `TOPOLOGY_wilson.dat` (centros vs. $k_2$) y `TOPOLOGY.txt`; `plot` dibuja el mapa de flujo y los centros.

**De dónde sale cada dato.**

| Dato | Origen | Detalle |
|---|---|---|
| $H(\mathbf R)$, $\mathbf R$, $\deg(\mathbf R)$ | `WANNIER_hr.dat` (Olla-DFT) o `seedname_hr.dat` (wannier90) | `wannier.leer_hr` |
| Autovectores en la malla | `wannier.interpolar` | coordenadas fraccionarias, fase $e^{2\pi i\mathbf k\cdot\mathbf R}$ |
| Ocupación | `--occupied` o `--fermi` (usuario) | nunca se adivina |
| Tolerancia de gap | `--gap-tol` | 1e-8 eV |

**Límites y trampas.**

- Solo secciones 2D: para un material 3D hay que barrer `--fixed` a mano; el número de Chern de una sección es el invariante de un aislante de Chern 2D o de una rebanada.
- "La señal cambia al invertir la orientación del plano" (texto del reporte): el signo de $C$ depende del orden `(a, b)` del plano elegido.
- No se calcula $\mathbb Z_2$: "no se asigna un Z2 automático sin comprobar simetría de reversión temporal". Con simetría de inversión temporal el Chern es siempre 0; los lazos de Wilson exportados permiten leer el $\mathbb Z_2$ a ojo, pero el código no lo hace.
- Si el gap directo se cierra en la malla: "el subespacio ocupado no está aislado… El número de Chern no está definido".
- Si el Chern discreto no cierra a entero con $10^{-6}$: "refina la malla y revisa la localización del modelo Wannier".
- El resultado hereda todos los defectos del modelo de Wannier (proyecciones malas, $H(\mathbf R)$ no decaído).

**Referencias.**

- T. Fukui, Y. Hatsugai, H. Suzuki, *J. Phys. Soc. Jpn.* **74**, 1674 (2005) — Chern discreto en malla. DOI 10.1143/JPSJ.74.1674.
- R. Yu, X. L. Qi, A. Bernevig, Z. Fang, X. Dai, *Phys. Rev. B* **84**, 075119 (2011) — lazos de Wilson y centros híbridos.
- A. A. Soluyanov, D. Vanderbilt, *Phys. Rev. B* **83**, 235401 (2011) — centros de Wannier híbridos e invariantes.
- D. Vanderbilt, *Berry Phases in Electronic Structure Theory* (Cambridge, 2018).
- X.-L. Qi, Y.-S. Wu, S.-C. Zhang, *Phys. Rev. B* **74**, 085308 (2006) — modelo de prueba usado en `tests/test_topology.py`.

---
### `olla-dft berry` — polarización por fase de Berry, cargas de Born

**Qué responde.** Cuánto cambia la polarización eléctrica de un cristal aislante al pasar de una estructura de referencia (normalmente la centrosimétrica) a la polar —la polarización espontánea de un ferroeléctrico—, y cuánta carga efectiva "se mueve" al desplazar un átomo (carga efectiva de Born $Z^*$).

**Fundamento para no expertos.** La polarización de un sólido periódico **no** se puede calcular como el momento dipolar de la celda: ese número depende de dónde se corten los bordes. King-Smith y Vanderbilt demostraron que lo que sí está bien definido es una fase geométrica (fase de Berry) acumulada por los estados ocupados al recorrer la zona de Brillouin a lo largo de una dirección. Esa fase está definida módulo $2\pi$, así que la polarización está definida módulo un "cuanto" $e\mathbf R/\Omega$: solo las **diferencias** entre dos estructuras conectadas por un camino son medibles, igual que en el experimento (se mide la carga que circula al cambiar la estructura, no $P$). pw.x calcula esa fase con `lberry = .true.` sobre "cuerdas" de puntos k paralelas a un vector recíproco; Olla-DFT prepara las cuerdas correctamente, sigue la rama de la fase a lo largo del camino, y comprueba la parte iónica contra su fórmula exacta.

**Fórmulas.**

Cuerdas de k (`qekit/modules/berry.py: cuerdas`): para cada punto $(i/n_\perp^{(1)}, j/n_\perp^{(2)})$ de la malla perpendicular (`--kperp` 6×6), `nppstr` puntos (9 por omisión) a lo largo de $\mathbf b_{\text{gdir}}$ con coordenada $l/(n_{\text{pp}}-1)$, $l = 0,\dots,n_{\text{pp}}-1$: el último punto es el primero más $\mathbf G$.

Fase iónica (`berry.fase_ionica`), en las unidades de QE (el cuanto vale `MOD_TOT` = 2 si todas las valencias son pares, 1 si alguna es impar; `berry.modulo_de`):

$$
\varphi_{\text{ion}} = \sum_a \left[Z_a f_a^{(g)}\right]_{\bmod\, m_a}\Big|_{\bmod\, m}, \qquad m_a = \begin{cases}1 & Z_a \text{ impar}\\ 2 & Z_a \text{ par}\end{cases}
$$

- $Z_a$: carga de valencia del pseudopotencial del átomo $a$ (electrones); $f_a^{(g)}$: coordenada fraccionaria a lo largo de `gdir`. El plegado por ion y el plegado final reproducen lo que hace pw.x; el plegado a $[-m/2, m/2)$ usa el `NINT` de Fortran (`berry._nint`, medio se aleja de cero), de modo que medio cuanto sale $-1$ como en QE.

Fase electrónica desde centros de Wannier (`berry.desde_wannier`), como comprobación independiente:

$$
\varphi_{\text{el}} = -f_s\sum_n \bar r_n^{(g)}, \qquad f_s = 2
$$

- $\bar r_n^{(g)}$: coordenada fraccionaria del centro de Wannier $n$ a lo largo de `gdir`; $f_s$ es el factor de espín. La fase total es $\varphi_{\text{el}} + \varphi_{\text{ion}}$ plegada.

Polarización y cuanto (`berry.polarizacion`, `berry.cuanto`):

$$
P_g = \varphi\,\frac{|\mathbf R_g|}{\Omega}, \qquad
\Delta P_{\text{cuanto}} = m\,\frac{|\mathbf R_g|}{\Omega}, \qquad
1\ e/\text{Å}^2 = 16.02176634\ \text{C/m}^2
$$

- $\varphi$: fase total en unidades de QE (adimensional, cuanto $m$); $\mathbf R_g$: vector de red `gdir` (Å); $\Omega$: volumen (Å³). $P_g$ es la **proyección** de $P\Omega/e$ sobre $\mathbf R_g$, no el módulo de $\mathbf P$.

Seguimiento de la rama (`berry.desenrollar`): $\tilde\varphi_0 = \varphi_0$, $\tilde\varphi_i = \varphi_i + m\cdot\mathrm{round}\big((\tilde\varphi_{i-1} - \varphi_i)/m\big)$; se avisa si algún salto $|\tilde\varphi_i - \tilde\varphi_{i-1}| > 0.25\,m$ (`FRACCION_SOSPECHOSA`).

Carga efectiva de Born (`berry.analizar`), con $\mathbf u$ el desplazamiento total (Å) y $\mathbf B_g$ el vector recíproco `gdir` (con $2\pi$):

$$
Z^*_{g} = \frac{2\pi\,\dfrac{d\tilde\varphi}{d\lambda}}{\mathbf u\cdot\mathbf B_g}
$$

- $d\tilde\varphi/d\lambda$: pendiente del ajuste lineal de la fase seguida frente a $\lambda \in [0,1]$ (`np.polyfit` grado 1 si hay más de 2 puntos; diferencia finita si no). Es la componente $Z^*_{g,\hat u}$ del tensor. Si $\mathbf u\perp\mathbf B_g$, no se calcula.

Camino adiabático (`berry._interpolar_estructuras`): posiciones interpoladas en fraccionarias por imagen mínima, $f(\lambda) = f_a + \lambda\,[(f_b - f_a) - \mathrm{round}(f_b - f_a)]$, y celda interpolada linealmente.

**Cómo lo calcula Olla-DFT.**

1. `qekit/cli.py: _cmd_berry` carga la estructura polar, opcionalmente `--reference` (centrosimétrica) o `--displace ATOMO:dx,dy,dz` (Å, átomo en base 1), y `--kperp`.
2. `berry.prepare` construye la lista de estructuras (`--nlambda` 5 puntos de $\lambda$; un solo punto si no hay camino), resuelve pseudos y cutoffs (`sweep.prepare_common(insulator=True)`) y, en cada `pNN/`, escribe `1_scf.in` y `2_berry.in` (`calculation='nscf'`, `occupations='fixed'`, `conv_thr 1e-10`, `nosym`, `noinv`, con `lberry`, `gdir` y `nppstr` insertados en `&CONTROL`), más `correr.sh`/`correr.py`.
3. `--run` → `berry.correr`: pw.x sobre scf y berry en cada punto, saltando los que ya tienen `JOB DONE` salvo `--redo`.
4. `--collect` → `berry.collect`: `leer_berry` extrae de `2_berry.out` `Ionic Phase`, `Electronic Phase`, `TOTAL PHASE`, `MOD_TOT`, `P = … (mod …) (e/Omega).bohr`, `direction of vector`, `Number of k-points per string`, `Number of different strings`; `valencias_de` lee la tabla "atomic species valence mass pseudopotential" de `1_scf.out`.
5. `berry.analizar`: desenrolla las fases, convierte a C/m², calcula $\Delta P$ y, si el camino es un desplazamiento, $Z^*$; `comprobar_ionica` compara la fase iónica de pw.x con $\sum Z_a f_a$ (aviso si difieren en más de $10^{-4}$).
6. `berry.report` imprime la tabla $\lambda$ / iónica / electrónica / total / seguida / $P$; `export` escribe `BERRY.dat` y `BERRY.txt`; `plot` dibuja $P(\lambda)$, los valores plegados de pw.x y una banda de anchura un cuanto.

**De dónde sale cada dato.**

| Dato | Origen | Detalle |
|---|---|---|
| Fases iónica, electrónica, total, `MOD_TOT` | `pNN/2_berry.out` de pw.x (`lberry`) | `berry.leer_berry`, expresión regular sobre el texto |
| Valencias $Z_a$ | tabla `atomic species / valence` de `1_scf.out` | `berry.valencias_de` |
| Celda, volumen, $\mathbf R_g$, $\mathbf B_g$ | estructura del usuario (última del camino) | `berry.cuanto`, `berry.analizar` |
| Conversión e/Å² → C/m² | constante `berry.E_A2_A_C_M2` | 16.02176634 |
| Umbral de salto sospechoso | `berry.FRACCION_SOSPECHOSA` | 0.25 del cuanto |
| Centros de Wannier (comprobación) | `olla-dft wannier` | `berry.desde_wannier`, solo por API/tests |

**Límites y trampas.**

- Solo para **aislantes**: el nscf se escribe con `occupations='fixed'`; en un metal la fase no está definida.
- Un solo punto no sirve: "Un solo punto. P está definida módulo el cuanto, así que este número por sí solo no significa nada".
- Si un paso mueve la fase más del 25 % del cuanto: "El seguimiento de la rama supone que el paso es pequeño; con saltos así, elegir la imagen más cercana es una apuesta. Sube --nlambda". Si $|\Delta P| > 0.9$ cuantos: "Comprueba con más puntos que no es un salto de rama disfrazado".
- Solo se calcula **una componente** (`--gdir`); para el vector $\mathbf P$ hacen falta tres corridas.
- Si pw.x se para con "Wrong k-strings", casi siempre faltó `nosym`/`noinv`; Olla-DFT los fuerza, pero un input editado a mano puede perderlos.
- En la figura, los marcadores "lo que escribe pw.x (plegado)" salen de `berry.polarizacion_plegada`: $P = \varphi_{\text{tot}}/m \cdot \Delta P_{\text{cuanto}}$ con el mismo `MOD_TOT` que usa `analizar`, de modo que coinciden con la rama seguida en $\lambda = 0$ y difieren de ella solo en múltiplos enteros del cuanto.
- `desde_wannier` (comprobación contra centros de Wannier) no está conectada a la CLI; solo se usa desde Python o en los tests.
- No se aplica ninguna corrección por espín polarizado ni SOC (el factor de espín es 2 fijo en `desde_wannier`; pw.x lo maneja internamente en `lberry`).

**Referencias.**

- R. D. King-Smith, D. Vanderbilt, *Phys. Rev. B* **47**, 1651 (1993) — teoría moderna de la polarización. DOI 10.1103/PhysRevB.47.1651.
- R. Resta, *Rev. Mod. Phys.* **66**, 899 (1994). DOI 10.1103/RevModPhys.66.899.
- N. A. Spaldin, "A beginner's guide to the modern theory of polarization", *J. Solid State Chem.* **195**, 2 (2012). DOI 10.1016/j.jssc.2012.05.010.
- D. Vanderbilt, *Berry Phases in Electronic Structure Theory* (Cambridge, 2018).

---

### `olla-dft hubbard` — U de Hubbard por respuesta lineal (hp.x)

**Qué responde.** Cuánto vale el parámetro $U$ de DFT+U para los orbitales localizados ($d$ o $f$) de tu sistema, calculado por respuesta lineal con `hp.x` en lugar de copiado de un artículo, y —con `--cycle`— su valor autoconsistente.

**Fundamento para no expertos.** Los funcionales semilocales (LDA, GGA) dejan que un electrón se "vea a sí mismo" (autointeracción), lo que deslocaliza en exceso los orbitales $d$ y $f$ y convierte a óxidos aislantes como el NiO en metales. DFT+U añade una penalización $U$ a la ocupación fraccionaria de esos orbitales. El valor de $U$ no es una propiedad del elemento sino del sistema y del *esquema de proyección* con el que se cuentan las ocupaciones. Cococcioni y de Gironcoli lo obtienen midiendo cómo responde la ocupación del orbital a una pequeña perturbación del potencial: la respuesta "desnuda" $\chi_0$ (sin dejar que el resto del sistema se reajuste) y la completa $\chi$. Su diferencia es la curvatura espuria que $U$ debe cancelar. `hp.x` hace ese cálculo con teoría de perturbaciones (DFPT) en una malla de vectores $\mathbf q$ que equivale a una supercelda. Como el $U$ obtenido depende del $U$ con el que se hizo el scf de partida, hay que iterar hasta que se estabilice.

**Fórmulas.**

Respuesta lineal (calculada por `hp.x`, no por Olla-DFT; `qekit/modules/hubbard.py`, docstring):

$$
U_I = \left(\chi_0^{-1} - \chi^{-1}\right)_{II}
$$

- $\chi_0$, $\chi$: matrices de respuesta de las ocupaciones $n_I$ del sitio de Hubbard $I$ a la perturbación $\alpha_J$ del potencial en el sitio $J$, sin y con reajuste autoconsistente (eV⁻¹). $U_I$ en eV.

Ciclo de autoconsistencia (`hubbard.ciclo`):

$$
U^{(k+1)}_s = (1 - \mu)\,U^{(k)}_s + \mu\,U^{\text{hp}}_s\!\left[U^{(k)}\right], \qquad
\text{convergido si } k \ge 1 \ \wedge\ \max_s\left|U^{\text{hp}}_s - U^{(k)}_s\right| < \text{tol}
$$

- $\mu$ = `--mixing` (1.0 por omisión), tol = `--tol` (0.05 eV), máximo `--max-iter` = 6 vueltas; $U^{(0)}_s$ = `U_SEMILLA` = $10^{-8}$ eV. El $U$ reportado por especie es la media sobre sus sitios (`HubbardRun.U`).

**Cómo lo calcula Olla-DFT.**

1. `qekit/cli.py: _cmd_hubbard` carga la estructura; `--species` o, por omisión, `hubbard.elementos_hubbard` (los de la tabla `ORBITAL_HUBBARD`: 3d Sc–Zn, 4d Y–Cd, 5d Hf–Hg, 4f La–Lu); `--qgrid` 2×2×2; `--hubbard-style legacy|card` (el mismo selector que `gen`).
2. `hubbard.prepare` escribe `scf.in` con `inputgen.build_pw_input(hubbard={s: U_semilla}, hubbard_style=…, conv_thr=1e-15)`. Con `legacy` (por omisión, QE ≤ 7.0): `lda_plus_u = .true.`, `Hubbard_U(i) = 1e-8` y `U_projection_type = 'ortho-atomic'` insertado en `&SYSTEM` (`_fijar_proyeccion`). Con `card` (QE ≥ 7.1): tarjeta `HUBBARD (<proyección>)` con `U El-orb 1e-8` y sin `U_projection_type`, que en esas versiones es un error. `--projection` admite `atomic`, `ortho-atomic`, `norm-atomic`, `wannier`, `pseudo`. Y `hp.in` (`build_hp_input`: `nq1..3`, `conv_thr_chi = 1e-8`, `iverbosity = 2`). Ocupaciones fijas salvo `--metal`; `--nspin 2` y `--mag` se pasan tal cual.
3. `--cycle` → `hubbard.ciclo`: por iteración crea `iterNN/`, corre pw.x (`runner.run_all`) y hp.x (`run_hp`, buscado junto a pw.x), lee `*.Hubbard_parameters.dat` (`collect` → `leer_parametros`, sección "Hubbard U parameters", columnas sitio/tipo/etiqueta/espín/nuevo tipo/nueva etiqueta/U) y mezcla.
4. `--collect` → `hubbard.collect` lee el primer `*.Hubbard_parameters.dat` de la carpeta; `--intersite` añade `leer_v` (sección "Hubbard V parameters", tabla átomo 1 / átomo 2 / distancia en bohr / V) y escribe `HUBBARD.card` con `tarjeta_hubbard` (`U El-orb valor` y `V El-orb El-orb i j valor`, con índices de la supercelda de hp.x y umbral `--v-threshold` 0.01 eV).
5. `hubbard.report` imprime la tabla de $U$ por sitio, la historia del ciclo y las advertencias; `export` escribe `HUBBARD_U.dat` y `HUBBARD_U.txt`, y sugiere la línea `olla-dft gen … --hubbard El=U`.

**De dónde sale cada dato.**

| Dato | Origen | Detalle |
|---|---|---|
| $U$ por sitio | `<prefix>.Hubbard_parameters.dat` de hp.x | `hubbard.leer_parametros` |
| $V$ intersitio, supercelda de vecinos | misma salida, sección "Hubbard V parameters" | `hubbard.leer_v` |
| Orbital corregido | tabla `hubbard.ORBITAL_HUBBARD` | por elemento; `3d` si no está en la tabla (`2p` para el segundo átomo de un V) |
| $U$ semilla | `hubbard.U_SEMILLA` | 1e-8 eV |
| `conv_thr` del scf, `conv_thr_chi` | constantes en `prepare` / `build_hp_input` | 1e-15 Ry, 1e-8 |
| Malla $\mathbf q$ | `--qgrid` | 2×2×2 (8 celdas) |

**Límites y trampas.**

- Olla-DFT **no calcula** $\chi$ ni $U$: los lee de hp.x. Sin hp.x compilado (`make hp`) el comando falla: "no se encontró hp.x junto a pw.x".
- Por omisión el scf usa la sintaxis `lda_plus_u`/`Hubbard_U(i)` (QE ≤ 7.0); con QE ≥ 7.1 hay que pedir `--hubbard-style card`. El docstring de `tarjeta_hubbard` advierte que la tarjeta "está probado contra la sintaxis documentada, no contra una corrida de QE 7.1, porque el QE de esta máquina es 6.6".
- Una sola vuelta da "U de PRIMERA ITERACIÓN. Depende del U que llevaba el scf de partida".
- Con `nq = 1×1×1`: "la perturbación ve sus propias imágenes periódicas y el U sale mal. Usa al menos 2x2x2".
- El $U$ "solo vale con la MISMA proyección"; el reporte lo repite en cada salida.
- Si el ciclo no converge en `--max-iter`: "Se hicieron N vueltas sin bajar de tol eV… si el número oscila arriba y abajo, baja --mixing a 0.5; si baja despacio pero siempre en el mismo sentido, sube --max-iter"; el comando devuelve código de salida 1.
- El orbital de la tarjeta HUBBARD para un elemento fuera de la tabla es `3d` (o `2p` como segundo átomo de un $V$), que puede ser incorrecto (por ejemplo `O-2p` está bien, `S` recibiría `2p`).
- Los índices de los pares $V$ están en la numeración de la **supercelda** de hp.x; la tarjeta los copia tal cual, como exige QE.

**Referencias.**

- M. Cococcioni, S. de Gironcoli, *Phys. Rev. B* **71**, 035105 (2005) — U por respuesta lineal. DOI 10.1103/PhysRevB.71.035105.
- I. Timrov, N. Marzari, M. Cococcioni, *Phys. Rev. B* **98**, 085127 (2018) — hp.x, DFPT para U. DOI 10.1103/PhysRevB.98.085127.
- I. Timrov, N. Marzari, M. Cococcioni, *Phys. Rev. B* **103**, 045141 (2021) — U y V autoconsistentes, ortho-atomic. DOI 10.1103/PhysRevB.103.045141.
- V. L. Campo Jr., M. Cococcioni, *J. Phys.: Condens. Matter* **22**, 055602 (2010) — DFT+U+V.
- S. L. Dudarev et al., *Phys. Rev. B* **57**, 1505 (1998) — formulación simplificada de DFT+U usada por QE.

---

### `olla-dft align` — alineamiento de bandas entre dos materiales

**Qué responde.** Dónde queda la banda de valencia (y la de conducción) de un material respecto a la del otro cuando se ponen en contacto: los *offsets* $\Delta E_v$ y $\Delta E_c$ y el tipo de heterounión (I anidada, II escalonada, III rota).

**Fundamento para no expertos.** Cada cálculo periódico fija el cero de su potencial de forma arbitraria (el término $G = 0$ del potencial de Hartree), así que restar directamente los VBM de dos cálculos distintos da un número sin significado. Hay dos maneras de ponerlos en una escala común. En el **modo vacío**, cada material se calcula como una losa con vacío y se mide su VBM respecto al nivel de vacío de su propio cálculo (su potencial de ionización); el offset es la diferencia. Ignora la carga que se transfiere al formar el contacto. En el **modo interfaz** (Van de Walle y Martin) se calculan los dos bultos y además la interfaz, y el potencial electrostático promediado macroscópicamente en cada lado de la interfaz sirve de puente entre las dos escalas: es el único término que sabe del dipolo de contacto.

**Fórmulas.**

Offsets (`qekit/modules/align.py: alinear`), con $E_v^{A}$ el VBM y $V^{A}_{\text{ref}}$ la referencia del cálculo de $A$:

$$
\Delta E_v = \left(E_v^{A} - V_{\text{ref}}^{A}\right) - \left(E_v^{B} - V_{\text{ref}}^{B}\right) + \Delta\bar V, \qquad
\Delta E_c = \left(E_c^{A} - V_{\text{ref}}^{A}\right) - \left(E_c^{B} - V_{\text{ref}}^{B}\right) + \Delta\bar V
$$

- Modo vacío: $V_{\text{ref}}$ = nivel de vacío (máximo del potencial planar, media de una ventana del 20 % alrededor; `fields.work_function`), $\Delta\bar V = 0$.
- Modo interfaz: $V_{\text{ref}}$ = potencial electrostático medio de la celda de bulto ($\langle V\rangle$ del promedio planar), y $\Delta\bar V = \bar V_A - \bar V_B$ medido en la interfaz (`align.puente_interfaz`).
- Todo en eV; el potencial de pp.x (`plot_num = 11`, $V_{\text{bare}} + V_H$) viene en Ry y se multiplica por `RY_EV = 13.605693122994`.

Puente de la interfaz (`align.puente_interfaz`): promedio planar $\bar V(z)$ del cube, promedio macroscópico móvil periódico de ventana $w$ (`fields.macroscopic_average`; $w$ = `--window` o $L/8$), y

$$
\bar V_A = \langle \bar{\bar V}\rangle_{z \in [L/8,\, L/4]}, \qquad \bar V_B = \langle \bar{\bar V}\rangle_{z \in [5L/8,\, 3L/4]}, \qquad \Delta\bar V = \bar V_A - \bar V_B
$$

- Se supone que el material $A$ ocupa la primera mitad de la celda de la interfaz y $B$ la segunda.

Tipo de alineamiento (`align.alinear`), en la escala de $B$ (VBM de $B$ en 0): $v_A = \Delta E_v$, $c_A = E_g^{B} + \Delta E_c$, $v_B = 0$, $c_B = E_g^{B}$:

- `=` si $|\Delta E_v| < 0.05$ y $|\Delta E_c| < 0.05$ eV (`TOL_ALINEADOS`);
- I si un gap contiene al otro ($v_A \le v_B \wedge c_A \ge c_B$, o al revés);
- III si $c_A \le v_B$ o $c_B \le v_A$;
- II en cualquier otro caso.

**Cómo lo calcula Olla-DFT.**

1. `qekit/cli.py: _cmd_align` recibe las carpetas `a` y `b`, `--interface CARPETA` (activa el modo interfaz), `--axis` (c por omisión), `--window`, `--names`.
2. `align.leer_lado` lee el XML (`qeout.read_xml`): VBM = `<highestOccupiedLevel>`, CBM = `<lowestUnoccupiedLevel>`, $E_F$; sin HOMO falla: "no da un VBM. En un metal no hay banda de valencia que alinear; y si es un aislante, al cálculo le faltan bandas vacías (nbnd) o no usó occupations='fixed'". Sin LUMO se marca `es_metal` y solo se da $\Delta E_v$.
3. `align._potencial` reutiliza `potencial.cube` o ejecuta `pp.x` (`fields.run_pp` con `plot_num = 11`, `output_format = 6`) y lo lee con `fields.read_cube`.
4. Modo vacío: `fields.work_function` da `v_vacuum` y la planitud de la meseta; modo interfaz: `fields.planar_average` y la media.
5. Con `--interface`, `align.puente_interfaz` calcula $\Delta\bar V$.
6. `align.alinear`, `report` (tabla VBM/CBM/gap respecto a la referencia, offsets, tipo y a qué material va cada portador en tipo II), `export` (`ALINEAMIENTO.dat`, `.txt`) y `plot` (diagrama de cajas). Las posiciones de las cajas salen de `align.posiciones_en_escala_de_b` —$v_A = \Delta E_v$, $c_A = E_g^{B} + \Delta E_c$, $v_B = 0$, $c_B = E_g^{B}$—, la misma convención con la que `alinear` clasifica el tipo, de modo que reporte, exportación y figura no pueden discrepar.

**De dónde sale cada dato.**

| Dato | Origen | Detalle |
|---|---|---|
| VBM, CBM, $E_F$ | `<highestOccupiedLevel>`, `<lowestUnoccupiedLevel>`, `<fermi_energy>` de `prefix.xml` | `align.leer_lado`; requiere ocupaciones fijas |
| Potencial electrostático | `potencial.cube` de pp.x (`plot_num = 11`) | `fields.read_cube`; Ry → eV con 13.605693122994 |
| Nivel de vacío y planitud | máximo del promedio planar, ventana del 20 % | `fields.work_function` |
| Ventana macroscópica | `--window` o $L/8$ | `align.puente_interfaz` |
| Umbral "alineados" | `align.TOL_ALINEADOS` | 0.05 eV |
| Umbral de planitud | constante en `alinear` | 0.05 eV |

**Límites y trampas.**

- Modo vacío: el reporte avisa siempre: "son las dos superficies AISLADAS. Al ponerlas en contacto se transfiere carga y aparece un dipolo de interfaz que desplaza el offset, típicamente entre 0.1 y 0.5 eV".
- Si la meseta de vacío varía más de 0.05 eV: "O falta vacío, o la losa tiene dipolo neto: usa --dipole al generarla. El nivel de vacío es la referencia de todo esto, así que ese error entra entero en el offset".
- El modo interfaz supone que $A$ está en la primera mitad de la celda y $B$ en la segunda, y toma dos ventanas fijas ($[L/8, L/4]$ y $[5L/8, 3L/4]$); una interfaz asimétrica o con capas de distinto espesor da un puente equivocado sin aviso.
- El VBM/CBM se leen de `highestOccupiedLevel`/`lowestUnoccupiedLevel`, que dependen de la malla k del scf; no se hace un análisis de bandas.
- Si $A$ no tiene CBM (metal o sin bandas vacías), la figura dibuja su gap como una caja de altura $E_g^{A}$ (o 1 eV si tampoco hay gap) sobre $v_A$: es un relleno visual, no un dato.
- Los offsets llevan el error sistemático del funcional; `TIPOS["="]` recuerda que "con funcionales semilocales el error frente al experimento es de varias décimas".

**Referencias.**

- C. G. Van de Walle, R. M. Martin, *Phys. Rev. B* **35**, 8154 (1987) — alineamiento por potencial macroscópico. DOI 10.1103/PhysRevB.35.8154.
- A. Baldereschi, S. Baroni, R. Resta, *Phys. Rev. Lett.* **61**, 734 (1988) — promedio macroscópico. DOI 10.1103/PhysRevLett.61.734.
- L. Kleinman, *Phys. Rev. B* **24**, 7412 (1981) — el cero arbitrario del potencial en cálculos periódicos.
- J. Tersoff, *Phys. Rev. B* **30**, 4874 (1984) — alineamiento y dipolos de interfaz.

## Mecánica, vibraciones, temperatura y transporte

Esta parte documenta la física que Olla-DFT implementa en los comandos que van de la energía total a las propiedades mecánicas, vibracionales, térmicas y de transporte: desde las pruebas de convergencia (`converge`, `tune`) y la ecuación de estado (`eos`), pasando por las constantes elásticas (`elastic`, `derived`), las deformaciones (`strain`), las superficies y los materiales laminares (`gamma`, `layers`, `xrd`, `exfoliate`), los fonones y todo lo que sale de ellos (`phonons`, `qha`, `thermochem`, `kappa`, `elph`), la dinámica molecular (`md`), el transporte electrónico difusivo y balístico (`transport`, `ballistic`) y el estimador de coste (`cost`). Cada sección se ha escrito leyendo el código de `qekit/modules/*.py` y `qekit/cli.py`, y solo recoge las fórmulas que el código ejecuta de verdad, con sus constantes y valores por omisión tal como están escritos. Cuando un docstring promete algo que el código no hace, se dice en "Límites y trampas". Nota sobre nombres de archivo: el módulo `qekit/modules/thermo.py` NO contiene la termodinámica armónica (esa vive en `phonons.thermodynamics`) sino el casco convexo de energías de formación del comando `hull`, que se documenta en otra parte.

---

### `olla-dft converge` — Convergencia de cutoffs y malla k

**Qué responde.** ¿A partir de qué `ecutwfc`, `ecutrho` o malla de puntos k la energía total deja de cambiar más de un umbral (por omisión 1 meV/átomo)? Es la primera pregunta con cualquier sistema nuevo.

**Fundamento para no expertos.** Un cálculo de ondas planas describe los electrones con una suma de ondas; `ecutwfc` fija cuántas ondas entran (la "resolución" de las funciones de onda), `ecutrho` la de la densidad de carga y la malla k cuántos puntos de la zona de Brillouin se muestrean. Con pocas ondas o pocos puntos el resultado es tosco; con demasiados, el cálculo cuesta más sin ganar nada. La prueba de convergencia repite el mismo cálculo subiendo el parámetro y mira cuándo la energía "se aplana", como quien afina el zoom de un microscopio hasta que la imagen deja de cambiar.

El criterio que usa Olla-DFT tiene una sutileza importante: compara cada punto contra el MÁS DENSO de la serie (el último), no contra el vecino anterior. Dos puntos contiguos pueden parecerse por casualidad en mitad de una curva que todavía no ha aplanado; compararlos entre sí es el error habitual.

**Fórmulas.** Diferencia por átomo respecto del punto más denso (`converge.ConvergenceRun.per_atom_diffs`):

$$
\Delta E_i = \frac{|E_i - E_{\mathrm{ref}}|}{N_{\mathrm{at}}} \times 1000
$$

- $E_i$: energía total del punto $i$, en eV por celda (leída del XML y convertida de Hartree con $27.211386245988$ eV).
- $E_{\mathrm{ref}}$: energía del último punto que terminó (el más denso).
- $N_{\mathrm{at}}$: átomos de la celda.
- $\Delta E_i$: en meV/átomo.

Índice de convergencia (`converge.ConvergenceRun.converged_index`): el primer $i$ tal que todos los $\Delta E_j$ con $j \ge i$ cumplen $\Delta E_j \le$ umbral (los puntos fallidos se ignoran). Malla k a partir de un espaciado (`kpoints.kgrid_from_spacing`):

$$
n_i = \left\lceil \frac{|\mathbf{b}_i|}{k_{\mathrm{spacing}}} \right\rceil, \qquad |\mathbf{b}_i| \text{ con el factor } 2\pi
$$

- $\mathbf{b}_i$: vectores recíprocos en Å⁻¹; $k_{\mathrm{spacing}}$ en Å⁻¹ (configuración `kspacing`, por omisión 0.20). Las direcciones con vacío ≥ 8 Å reciben un solo punto.

**Cómo lo calcula Olla-DFT.**
1. `qekit/cli.py: _cmd_converge` carga la estructura y llama a `qekit/modules/converge.py: prepare`.
2. `sweep.prepare_common` resuelve pseudopotenciales y cutoffs (`pseudo.recommend_cutoffs`: el máximo que declaran los UPF; si no declaran, `ecutwfc` de configuración (60 Ry) y `dual` (8); `ecutrho` nunca por debajo de $4\,\mathrm{ecutwfc}$).
3. Serie por omisión: `ecutwfc` = 30, 40, …, 100 Ry con `ecutrho = dual × ecutwfc`; `ecutrho` = 4, 6, 8, 10, 12 × ecutwfc; `kmesh` = mallas de los espaciados 0.40, 0.30, 0.25, 0.20, 0.15, 0.12 Å⁻¹ (sin repetir). Con `--values` se sustituye la serie (para `kmesh` admite `8x8x8` o espaciados).
4. Un `pw.in` (`calculation='scf'`, `conv_thr = 1e-8`, `tstress`/`tprnfor` activados) por punto vía `sweep.write_scf_job`, más `run.sh` y `run.py`.
5. Con `--run`, `runner.run_all` ejecuta `pw.x`; con `--collect`, `converge.collect` lee `out/*.xml` (`qeout.read_xml`, etiqueta `<total_energy><etot>`).
6. `converge.report` imprime la tabla, el punto de convergencia y la recomendación (`--ecutwfc N` o la malla); `converge.export` escribe `CONVERGENCIA.dat` y `.txt`; `converge.plot` dibuja $|\Delta E|$ en escala logarítmica con la banda del umbral.

**De dónde sale cada dato.**

| Dato | Origen | Detalle |
|---|---|---|
| Energía total | XML de pw.x (`output/total_energy/etot`, Hartree) | `qeout.read_xml` → eV |
| Convergencia del scf | XML (`convergence_info/scf_conv/convergence_achieved`) | un punto no convergido cuenta como fallido en `--run` |
| Umbral | parámetro `--threshold` | 1.0 meV/átomo por omisión |
| Cutoffs base | cabecera de los UPF o `olla-dft config` | `pseudo.recommend_cutoffs` |
| Malla k fija | `sweep.default_grid` | `kspacing` de configuración (0.20 Å⁻¹) |
| Ry ↔ eV | `qeout.RY_EV` | 13.605693122994 eV |

**Límites y trampas.** Solo mira la energía total; el reporte avisa: "la convergencia depende de la propiedad: la energía total converge antes que los esfuerzos o los fonones". Si solo el último punto cumple, dice: "Solo el último punto queda bajo … no hay margen para asegurar que ahí ya aplanó". Si ninguno cumple: "NO converge dentro de … Extiende la serie hacia valores más densos". El campo `energies` del dataclass está comentado como "eV por celda", pero la tabla se imprime en Ry (se divide por `RY_EV`): no es un error, solo una conversión de presentación. Con `--collect` no se reescriben los inputs (`sweep.set_write_inputs(False)`), así que el reporte describe lo que realmente corrió.

**Referencias.** Manual de Quantum ESPRESSO (`pw.x`, variables `ecutwfc`, `ecutrho`, `K_POINTS`). Monkhorst y Pack, *Phys. Rev. B* 13, 5188 (1976), DOI 10.1103/PhysRevB.13.5188.

---

### `olla-dft tune` — Recomendación adaptativa de convergencia

**Qué responde.** Dado un `CONVERGENCIA.dat` ya generado, ¿está convergida la serie, y si no, qué valor conviene probar a continuación?

**Fundamento para no expertos.** Es un post-proceso puro sobre la tabla de `converge`: aplica el mismo criterio ("desde este punto, toda la cola queda dentro del umbral") y, cuando no se cumple, propone el siguiente valor con un paso razonable, en vez de dejar que el usuario adivine.

**Fórmulas.** Criterio (`tuning.analyze`): índice $i$ mínimo con $|\Delta E_j| \le$ umbral para todo $j \ge i$. Estados: `ready` (existe $i$ y no es el último), `confirm` (solo el último cumple), `extend` (ninguno). Siguiente valor (`tuning._next_value`):

$$
v_{\mathrm{next}} = v_{\mathrm{last}} + \max\!\left(\mathrm{mediana}\{v_{k+1}-v_k > 0\},\; 0.10\,|v_{\mathrm{last}}|\right)
$$

- Con menos de dos valores, o sin pasos positivos: $v_{\mathrm{next}} = 1.25\,v_{\mathrm{last}}$ (o $v_{\mathrm{last}}+1$ si no es positivo).

**Cómo lo calcula Olla-DFT.**
1. `qekit/cli.py: _cmd_tune` → `qekit/modules/tuning.py: read` lee las filas numéricas (columna 1 valor, 2 energía en Ry, 3 $\Delta E$ en meV/átomo; ignora comentarios y NaN).
2. `tuning.analyze` aplica el criterio y elige el estado y el valor recomendado.
3. `tuning.report` lo imprime; con `-o`, `tuning.export` escribe un JSON (`CONVERGENCIA_RECOMENDACION.json` por omisión).

**De dónde sale cada dato.**

| Dato | Origen | Detalle |
|---|---|---|
| Valor, E, ΔE | `CONVERGENCIA.dat` (de `olla-dft converge`) | `tuning.read`; ΔE se toma en valor absoluto |
| Umbral | `--threshold` | 1.0 meV/átomo si se omite; debe ser > 0 |

**Límites y trampas.** No corre nada ni lee salidas de QE: solo la tabla. Usa la columna ΔE tal como está escrita, que `converge` calculó contra el punto más denso de ESA serie; si se añaden puntos después, hay que regenerar la tabla. El reporte recuerda: "La propiedad energía puede converger antes que fuerzas, fonones o tensores".

**Referencias.** Ninguna específica; es la misma lógica de `converge`.

---

### `olla-dft eos` — Ecuación de estado E–V y módulo de bulk

**Qué responde.** ¿Cuál es el volumen de equilibrio $V_0$, la energía mínima $E_0$, el módulo de bulk $B_0$ y su derivada $B_0'$ del cristal? Y, si es cúbico, el parámetro de red $a_0$.

**Fundamento para no expertos.** Se comprime y se estira la celda un poco alrededor del tamaño de partida, se calcula la energía a cada volumen y se ajusta una curva con forma de valle. El fondo del valle es el volumen de equilibrio; la "rigidez" del valle (su curvatura) es el módulo de bulk, que mide cuánto hay que apretar para reducir el volumen. Olla-DFT ajusta tres ecuaciones distintas; si las tres dan lo mismo, el ajuste es fiable, y si discrepan, suele faltar rango o sobran puntos ruidosos.

**Fórmulas.** Birch–Murnaghan de tercer orden (`eos.birch_murnaghan`), con $\eta = (V_0/V)^{2/3}$:

$$
E(V) = E_0 + \frac{9 V_0 B_0}{16}\left[(\eta-1)^3 B_0' + (\eta-1)^2 (6 - 4\eta)\right]
$$

Murnaghan (`eos.murnaghan`):

$$
E(V) = E_0 + \frac{B_0 V}{B_0'}\left[\frac{(V_0/V)^{B_0'}}{B_0'-1} + 1\right] - \frac{B_0 V_0}{B_0'-1}
$$

Vinet (`eos.vinet`), con $x = (V/V_0)^{1/3}$ y $\xi = \tfrac{3}{2}(B_0'-1)$:

$$
E(V) = E_0 + \frac{9 B_0 V_0}{\xi^2}\left[1 + \left(\xi(1-x) - 1\right) e^{\xi(1-x)}\right]
$$

- $V$, $V_0$: volúmenes en Å³; $E$, $E_0$: eV; $B_0$: eV/Å³ dentro del ajuste, convertido a GPa con `EV_A3_GPA = 160.21766208`; $B_0'$: adimensional.
- Semilla del ajuste (`eos.fit`): parábola $E = aV^2+bV+c$ → $V_0 = -b/2a$, $B_0 = 2aV_0$, $B_0' = 4$.
- RMSE: $\sqrt{\langle (E - E_{\mathrm{fit}})^2 \rangle}/N_{\mathrm{at}}$, en eV/átomo (se imprime en meV/át).
- Parámetro de red cúbico (`eos.fit`, campo `EOSFit.a0`, con $V_{\mathrm{conv}}/V_{\mathrm{prim}}$ medido en `prepare`): $a_0 = (V_0 \cdot V_{\mathrm{conv}}/V_{\mathrm{prim}})^{1/3}$.

**Cómo lo calcula Olla-DFT.**
1. `qekit/cli.py: _cmd_eos` → `qekit/modules/eos.py: prepare`. Exige `--npoints` ≥ 5 (por omisión 9); `--span` 0.10 (±10 % en VOLUMEN); `--scale` (factor lineal de centrado) 1.0.
2. Factores de volumen equiespaciados $f \in [c^3(1-s),\, c^3(1+s)]$; factor lineal $f^{1/3}$ aplicado a la celda con `set_cell(..., scale_atoms=True)`.
3. Se decide si es cúbico preguntando a spglib (`structure.symmetry_dataset`, grupo ≥ 195) y se guarda $V_{\mathrm{conv}}/V_{\mathrm{prim}}$ con `structure.conventional`.
4. Un `scf` (o `relax` con `--relax-ions`) por volumen, todos con la MISMA malla k (`sweep.default_grid`), escritos por `sweep.write_scf_job` en `V_<factor>/pw.in`.
5. `--run` ejecuta `pw.x`; `--collect`/`eos.collect` lee `etot` del XML.
6. `eos.fit_all` ajusta las tres ecuaciones con `scipy.optimize.curve_fit` (`maxfev=20000`); rechaza el ajuste si $V_0 \notin (0.6 V_{\min}, 1.4 V_{\max})$ o $B_0 \le 0$.
7. `eos.report` imprime tabla, los tres ajustes, el resultado de Birch–Murnaghan y la dispersión entre ecuaciones; `eos.export` escribe `EOS.dat` y `EOS.txt`; `eos.plot` dibuja $E - E_0$ con residuales.

**De dónde sale cada dato.**

| Dato | Origen | Detalle |
|---|---|---|
| Energía total por volumen | XML de pw.x (`etot`) | `qeout.read_xml` |
| Volumen de cada punto | celda escalada (ASE) | $|\det(\mathbf{a})|$ en Å³ |
| Simetría cúbica y celda convencional | spglib vía `structure` | `symmetry_dataset`, `conventional` |
| eV/Å³ → GPa | `eos.EV_A3_GPA` | 160.21766208 |
| Ajuste no lineal | `scipy.optimize.curve_fit` | biblioteca |

**Límites y trampas.** No relaja la forma de la celda: solo escala isotrópicamente (para no cúbicos, $c/a$ queda fijo salvo que se use `--relax-ions`, que relaja posiciones, no la celda). Avisa: "V₀ cae FUERA del rango calculado. Vuelve a correr el barrido centrado en ese volumen" y, si las tres ecuaciones difieren más del 5 %: "suele indicar que faltan puntos o que el rango de volúmenes es muy estrecho". `--relax-ions` usa `calculation='relax'`, por lo que los inputs llevan `forc_conv_thr = 1e-4`.

**Referencias.** F. Birch, *Phys. Rev.* 71, 809 (1947), DOI 10.1103/PhysRev.71.809. F. D. Murnaghan, *Proc. Natl. Acad. Sci. USA* 30, 244 (1944). P. Vinet, J. Ferrante, J. R. Smith y J. H. Rose, *J. Phys. C* 19, L467 (1986).

---

### `olla-dft elastic` — Constantes elásticas por esfuerzo–deformación

**Qué responde.** ¿Cuáles son las constantes elásticas $C_{ij}$ del cristal (o de la lámina, con `--2d`), los módulos de bulk, cizalla y Young, la razón de Poisson, y es la estructura mecánicamente estable?

**Fundamento para no expertos.** Un sólido deformado ligeramente responde con un esfuerzo (fuerza por unidad de área) proporcional a la deformación: es la ley de Hooke generalizada, y las constantes de proporcionalidad son las $C_{ij}$. Olla-DFT deforma la celda un ±1 % (y ±0.5 %) en cada una de las seis direcciones independientes (tres estiramientos y tres cizallas), pide a `pw.x` el tensor de esfuerzos en cada una y ajusta una recta. Como `pw.x` da los seis esfuerzos de una vez, cada deformación aporta seis ecuaciones, muchas menos corridas que con el método de energía.

En una lámina (grafeno, MoS₂) no tiene sentido estirar la dirección del vacío: se deforman solo las dos direcciones del plano y la cizalla en el plano, y las constantes se dan en N/m multiplicando por la altura de la celda, de modo que el vacío se cancela.

**Fórmulas.** Deformación aplicada (`elastic.strain_matrix`): $\mathbf{a}' = \mathbf{a}(\mathbf{I}+\boldsymbol{\varepsilon})$ con $\varepsilon_{ii}=\delta$ para las normales y $\varepsilon_{ij}=\varepsilon_{ji}=\delta/2$ para las cizallas (convenio de Voigt, $\varepsilon_4 = 2\varepsilon_{23}$). Ajuste (`elastic.fit`), con el signo invertido porque el tensor que escribe `pw.x` es el opuesto del de la elasticidad:

$$
C_{ij} = -\frac{\partial\,\sigma^{\mathrm{pw}}_i}{\partial \varepsilon_j}\Big|_{\text{mínimos cuadrados}}, \qquad \sigma^{\mathrm{pw}}_i \to \sigma^{\mathrm{pw}}_i - \sigma^{\mathrm{pw}}_i(\text{ref})
$$

Promedios de Voigt, Reuss y Hill (`elastic.moduli`), con $S = C^{-1}$:

$$
B_V = \frac{(C_{11}+C_{22}+C_{33}) + 2(C_{12}+C_{23}+C_{13})}{9}, \quad
G_V = \frac{(C_{11}+C_{22}+C_{33}) - (C_{12}+C_{23}+C_{13}) + 3(C_{44}+C_{55}+C_{66})}{15}
$$

$$
B_R = \frac{1}{(S_{11}+S_{22}+S_{33}) + 2(S_{12}+S_{23}+S_{13})}, \quad
G_R = \frac{15}{4(S_{11}+S_{22}+S_{33}) - 4(S_{12}+S_{23}+S_{13}) + 3(S_{44}+S_{55}+S_{66})}
$$

$$
B_H = \tfrac{1}{2}(B_V+B_R),\quad G_H = \tfrac{1}{2}(G_V+G_R),\quad
E = \frac{9 B_H G_H}{3B_H + G_H},\quad \nu = \frac{3B_H - 2G_H}{2(3B_H+G_H)},\quad
A^U = 5\frac{G_V}{G_R} + \frac{B_V}{B_R} - 6
$$

- $C_{ij}$, $B$, $G$, $E$ en GPa; $\nu$ y $A^U$ adimensionales; cociente de Pugh $B_H/G_H$ (umbral de ductilidad 1.75).
- Estabilidad (Born generalizado): todos los valores propios de $\tfrac{1}{2}(C+C^T)$ positivos.

Lámina (`elastic.constantes_2d`, `modulos_2d`, `born_2d`): $C^{2D}_{ij} = C_{ij}\,c\times 0.1$ (GPa·Å → N/m), con $c$ la altura de la celda;

$$
Y_x = \frac{C_{11}C_{22}-C_{12}^2}{C_{22}},\quad \nu_x = \frac{C_{12}}{C_{22}},\quad
K = \frac{C_{11}+C_{22}+2C_{12}}{4},\quad G = C_{66};\qquad
C_{11}>0,\; C_{66}>0,\; C_{11}C_{22}-C_{12}^2>0
$$

**Cómo lo calcula Olla-DFT.**
1. `qekit/cli.py: _cmd_elastic` → `qekit/modules/elastic.py: prepare`. En 3D la estructura se lleva SIEMPRE a la primitiva estandarizada de spglib (`structure.primitive`) para que los ejes cartesianos coincidan con los cristalofísicos; en `--2d` no (exige vacío en $c$ vía `kpoints.direcciones_con_vacio`).
2. Familia cristalina por número de grupo espacial (`elastic.crystal_family`: ≥195 cúbico, ≥168 hexagonal, ≥143 trigonal, ≥75 tetragonal, ≥16 ortorrómbico).
3. Deformaciones: `--delta` 0.010, `--npoints` 4 (par) → ±δ/2, ±δ. Componentes: 6 de Voigt, o (1, 2, 6) en 2D. Más una celda de referencia sin deformar.
4. `--ion-mode auto` (por omisión): `scf` (iones fijos) en ε1–ε3 y `relax` en ε4–ε6; `relax`: todas relajadas; `fixed`: todas fijas. Los deformados llevan `conv_thr = 1e-9`.
5. `pw.x` con `tstress = .true.`; `elastic.collect` lee `<stress>` del XML (Ha/bohr³ → GPa con `qeout.HA_BOHR3_GPA = 29421.026`).
6. `elastic.fit` ajusta columna a columna con `np.polyfit(..., 1)`; `elastic.symmetrize` promedia los equivalentes de la familia (cúbico, hexagonal con $C_{66}=(C_{11}-C_{12})/2$, tetragonal parcial); `elastic.moduli` calcula VRH y Born.
7. `elastic.report`/`_report_2d`, `elastic.export` (`ELASTIC_C.dat`, `ELASTIC.txt`), `elastic.plot` (rectas σ–ε).

**De dónde sale cada dato.**

| Dato | Origen | Detalle |
|---|---|---|
| Tensor de esfuerzos | XML de pw.x (`output/stress`, Ha/bohr³, orden Fortran) | `qeout.read_xml`; requiere `tstress=.true.` (siempre puesto) |
| Grupo espacial y familia | spglib (`structure.symmetry_dataset`) | `elastic.crystal_family` |
| Altura de la celda (2D) | `|a_3|` de la celda de entrada | `ElasticRun.altura` |
| GPa·Å → N/m | `elastic.GPA_A_NM` | 0.1 |
| Espesor supuesto (2D) | `--thickness` | solo para el equivalente en GPa, "convenio, no medida" |

**Límites y trampas.** La simetrización solo cubre cúbico, hexagonal y (parcialmente) tetragonal; trigonal, ortorrómbico y monoclínico/triclínico quedan con la matriz simetrizada $\tfrac{1}{2}(C+C^T)$ sin más. El criterio de Born es el general (valores propios), no las desigualdades específicas de cada familia. Avisa si el esfuerzo residual de la celda de referencia supera 0.5 GPa: "es alto. Relaja la celda con vc-relax antes de calcular las constantes elásticas". En 2D advierte que con `--ion-mode auto` la identidad $C_{66}=(C_{11}-C_{12})/2$ deja de cumplirse aunque la lámina sea isótropa. El equivalente en GPa de una lámina depende del espesor elegido: "Este espesor es un CONVENIO, no una medida". Sin al menos 3 esfuerzos leídos no ajusta.

**Referencias.** R. Hill, *Proc. Phys. Soc. A* 65, 349 (1952), DOI 10.1088/0370-1298/65/5/307. S. I. Ranganathan y M. Ostoja-Starzewski, *Phys. Rev. Lett.* 101, 055504 (2008), DOI 10.1103/PhysRevLett.101.055504 (índice $A^U$). F. Mouhat y F.-X. Coudert, *Phys. Rev. B* 90, 224104 (2014), DOI 10.1103/PhysRevB.90.224104 (criterios de Born). S. F. Pugh, *Philos. Mag.* 45, 823 (1954).

---

### `olla-dft strain` — Barrido de deformación: gap, energía y momento

**Qué responde.** ¿Cómo cambian la energía, el band gap, la presión y el momento magnético al deformar la celda (biaxial, uniaxial, hidrostática o cizalla)? ¿Cuál es el potencial de deformación $dE_{\mathrm{gap}}/d\varepsilon$ y a qué deformación se cierra el gap?

**Fundamento para no expertos.** Estirar o comprimir un cristal cambia las distancias entre átomos y con ellas la estructura electrónica: el gap puede abrirse, cerrarse o cambiar de tipo, y un material magnético puede perder su momento. El potencial de deformación es la pendiente de esa respuesta, y es lo que se compara con experimentos de presión o con láminas sobre sustratos que las estiran. Olla-DFT aplica cada deformación SIEMPRE sobre la celda original (no sobre la del punto anterior, que acumularía error) y relaja las posiciones internas en cada punto.

**Fórmulas.** Deformación (`strain.matriz`): $\mathbf{a}' = \mathbf{a}_0(\mathbf{I}+\boldsymbol{\varepsilon})$ con las componentes de Voigt de cada modo (`strain.MODOS`: biaxial (xx, yy), uniaxial-a/b/c, hidrostática (xx, yy, zz), cizalla xy con $\varepsilon_{xy}=\varepsilon_{yx}=\varepsilon/2$). Mínimo de energía por parábola local (`strain.minimo`, hasta 3 puntos a cada lado del mínimo muestreado): $\varepsilon^* = -b/2a$. Potencial de deformación (`strain.potencial_deformacion`):

$$
E_{\mathrm{gap}}(\varepsilon) \approx m\,\varepsilon + b, \qquad R^2 = 1 - \frac{\sum (y - \hat y)^2}{\sum (y-\bar y)^2}
$$

- $m$: en eV por unidad de deformación (fracción, no por ciento). Gap $= E_{\mathrm{LUMO}} - E_{\mathrm{HOMO}}$ del XML.
- Cierre del gap (`strain.cierre_de_gap`): interpolación lineal de la deformación en que el gap cruza 0.02 eV.

Módulo biaxial 2D (`strain.modulo_biaxial`, solo modo biaxial, puntos con $|\varepsilon|\le 0.03$):

$$
Y_{2D} = \frac{1}{A_0}\frac{d^2E}{d\varepsilon^2} \times 16.021766 \;\; [\mathrm{N/m}], \qquad \frac{d^2E}{d\varepsilon^2} = 2a
$$

- $A_0 = |\mathbf{a}_1\times\mathbf{a}_2|$ en Å²; $E$ en eV; 1 eV/Å² = 16.021766 N/m. Es la combinación $C_{11}+2C_{12}+C_{22}$, NO el módulo de Young.

**Cómo lo calcula Olla-DFT.**
1. `qekit/cli.py: _cmd_strain` → `qekit/modules/strain.py: prepare`. `--range MIN:MAX:N` en POR CIENTO (por omisión `-5:5:11`; rechaza $|\varepsilon| > 30$ %, exige $N \ge 3$ y añade $\varepsilon=0$ si falta).
2. Estima `nbnd` con bandas vacías (`inputgen._estimate_nbnd`: $\lceil 1.25\,N_{\mathrm{occ}} + 4\rceil$, ×1.2+2 si `nspin=2`) para que exista LUMO.
3. Tipo de cálculo: `relax` (por omisión), `scf` con `--fixed-ions`, `vc-relax` con `cell_dofree` (`z`, `shape` o `2Dxy` según el modo) con `--relax-perp`.
4. Un input por deformación (`sweep.write_scf_job`, admite `--nspin/--mag`, `--hubbard`, `--vdw`).
5. `strain.collect` lee de cada XML `etot`, `highestOccupiedLevel`, `lowestUnoccupiedLevel`, `stress` (presión = traza/3), `magnetization/total` y `convergence_achieved`.
6. `strain.report` imprime la tabla, el mínimo, el potencial de deformación, el cierre del gap, el momento y el módulo biaxial (si hay vacío en $c$); `strain.export` (`STRAIN.dat`, `.txt`); `strain.plot` (dos paneles).

**De dónde sale cada dato.**

| Dato | Origen | Detalle |
|---|---|---|
| Energía, HOMO, LUMO | XML de pw.x (`etot`, `highestOccupiedLevel`, `lowestUnoccupiedLevel`) | `qeout.read_xml` |
| Presión | XML (`stress`, traza/3, signo de QE) | `QEResult.pressure` en GPa |
| Momento magnético | XML (`magnetization/total`) | μ_B por celda |
| Convergencia | XML (`convergence_achieved`) | filas marcadas `<< SIN CONVERGER` |
| Área y volumen de referencia | celda de entrada | `StrainRun.area0`, `volume0` |
| eV/Å² → N/m | constante literal en `modulo_biaxial` | 16.021766 |

**Límites y trampas.** Si el HOMO existe pero no el LUMO (sin bandas vacías) la columna del gap queda vacía y avisa: "No hay gap en la tabla: los cálculos no tienen bandas vacías". Si el mínimo no cae en $\varepsilon=0$ (|ε*| > 0.3 %): "La estructura de partida no estaba relajada". Con $R^2 < 0.9$: "El gap no responde de forma lineal en este rango". Biaxial sin vacío: "si es material en bulto, quizá querías 'hidrostatica'". `--relax-perp` con hidrostática se rechaza. Los puntos no convergidos SÍ entran en la tabla (se leen del XML aunque el runner los marque fallidos) pero se avisa que "NO son comparables con el resto". El "gap" es el de la malla k del scf, no el gap fundamental de un camino de bandas.

**Referencias.** J. Bardeen y W. Shockley, *Phys. Rev.* 80, 72 (1950), DOI 10.1103/PhysRev.80.72 (potenciales de deformación). C. G. Van de Walle, *Phys. Rev. B* 39, 1871 (1989).

---

### `olla-dft gamma` — Energía de superficie y ajuste de Fiorentini–Methfessel

**Qué responde.** ¿Cuánta energía por unidad de área cuesta crear la superficie (hkl) de un cristal, $\gamma$ en J/m², y cómo converge con el grosor de la losa?

**Fundamento para no expertos.** Cortar un cristal deja átomos con menos vecinos: eso cuesta energía, y la energía de superficie es ese coste por unidad de área. Se calcula con una "losa" (unas capas atómicas con vacío encima y debajo) restando lo que valdrían los mismos átomos dentro del cristal. El problema es que la energía de bulto viene de OTRO cálculo, con otra malla k, y cualquier error residual por átomo se multiplica por el número de átomos: $\gamma$ no converge, deriva. La salida es ajustar una recta $E_{\mathrm{losa}}(N)$ sobre varios grosores: la ordenada al origen da $2\gamma A$ y la pendiente una energía de bulto consistente con las propias losas.

**Fórmulas.** Directa (`surfen.GammaRun.gamma_directo`):

$$
\gamma_{\mathrm{dir}}(N) = \frac{E_{\mathrm{losa}}(N) - N_{\mathrm{at}}\,E_{\mathrm{bulto}}}{2A}
$$

Ajuste (`surfen.ajustar`), mínimos cuadrados sobre los pares $(N_{\mathrm{at}}, E_{\mathrm{losa}})$:

$$
E_{\mathrm{losa}}(N_{\mathrm{at}}) = 2\gamma A + N_{\mathrm{at}}\,E_{\mathrm{bulto}}^{\mathrm{ajuste}}
$$

- $E_{\mathrm{losa}}$: eV por celda de losa; $N_{\mathrm{at}}$: átomos de la losa; $E_{\mathrm{bulto}}$: eV/átomo del cálculo aparte de la celda convencional; $A = |\mathbf{a}_1\times\mathbf{a}_2|$ en Å² (una cara); el 2 son las dos caras (`GammaRun.caras` siempre vale 2).
- $\gamma$ en eV/Å² → J/m² con `EV_A2_A_J_M2 = 16.021766`. Energía de escisión $= 2\gamma$.

**Cómo lo calcula Olla-DFT.**
1. `qekit/cli.py: _cmd_gamma` → `qekit/modules/surfen.py: prepare`. `--miller` (por omisión `1 0 0`), `--layers` 3,4,5,6 (al menos dos, mínimo 3 capas), `--vacuum` 20 Å.
2. `builder.surface` corta cada losa sobre la celda CONVENCIONAL con `ase.build.surface`, centra el vacío, detecta si es simétrica (perfil $z$ igual a su reflejo, tol 0.3 Å) y polar (composición de la capa superior ≠ inferior) y emite avisos.
3. Salvo `--no-reduce` o `--fix`, `surfen.reducir_losa` sustituye la losa por su primitiva de spglib si el eje $c$ no cambia (mismo $\gamma$, menos átomos).
4. Malla k fijada con la losa más pequeña y reutilizada en todas (`sweep.default_grid`); el bulto (celda convencional) lleva su propia malla. Cálculos `scf` o `relax` (`--relax`), con opciones `--vdw`, `--dipole` (`dipole_correction=3`), `--nspin/--mag`.
5. `surfen.collect` lee `etot` y `convergence_achieved` de cada XML y ajusta (`surfen.ajustar`).
6. `surfen.report` imprime la tabla de γ directa, la deriva, el ajuste con $R^2$ y la diferencia $E_{\mathrm{bulto}}^{\mathrm{ajuste}} - E_{\mathrm{bulto}}$; `surfen.export` (`GAMMA.dat`, `GAMMA.txt`); `surfen.plot`.

**De dónde sale cada dato.**

| Dato | Origen | Detalle |
|---|---|---|
| $E_{\mathrm{losa}}(N)$ | XML de cada `capasNN/` (`etot`) | `qeout.read_xml` |
| $E_{\mathrm{bulto}}$ | XML de `_bulto/` (`etot`) / átomos de la celda convencional | omitido con `--no-bulk` |
| Área $A$ | celda de la losa más fina | `np.cross(a1, a2)` |
| Simetría / polaridad | `builder.surface` (geometría) | tolerancia 0.3 Å |
| eV/Å² → J/m² | `surfen.EV_A2_A_J_M2` | 16.021766 |

**Límites y trampas.** Es el ajuste LINEAL de Fiorentini–Methfessel; el esquema incremental de Boettger (que toma $E_{\mathrm{bulto}}$ de la diferencia entre losas consecutivas) no está implementado, y el docstring lo aclara. Losa no simétrica: "γ es el PROMEDIO de sus dos caras, no el de una". Polar: "usa --dipole". Sin `--relax`: "Sin relajar: γ sale alta. La relajación superficial la baja entre un 5 y un 20 %". Si la deriva entre losas supera 0.05 J/m²: "No converge … Es el error residual de E_bulto multiplicado por el número de átomos, no física. El valor bueno es el del ajuste". $R^2 < 0.999$: "o falta convergencia en algún punto, o las losas finas todavía no tienen interior de bulto". Con `--fix` no se reduce la celda (las restricciones van sobre átomos concretos).

**Referencias.** V. Fiorentini y M. Methfessel, *J. Phys.: Condens. Matter* 8, 6525 (1996), DOI 10.1088/0953-8984/8/36/005. J. C. Boettger, *Phys. Rev. B* 49, 16798 (1994), DOI 10.1103/PhysRevB.49.16798.

---

### `olla-dft layers` — Detección de capas por conectividad

**Qué responde.** ¿Es laminar la estructura? ¿Cuántas capas hay por celda, en qué eje se apilan, cuál es el espaciado basal $d$ y el hueco interlaminar, y dónde caería el pico basal (00l) en un difractograma?

**Fundamento para no expertos.** No se mira la geometría "a ojo" sino los enlaces: dos átomos están enlazados si su distancia no supera la suma de radios covalentes más una tolerancia. Se construye la red de enlaces respetando la periodicidad y se separan las piezas conectadas. Una pieza que se repite en exactamente dos direcciones es una capa; en tres, un armazón 3D; en una, una cadena; en ninguna, una molécula. La dimensionalidad se lee de los "vectores de cierre": al recorrer los enlaces asignando a cada átomo un desplazamiento de celda respecto de un átomo raíz, cada enlace que "no cuadra" aporta un vector entero, y el rango del conjunto de esos vectores es el número de direcciones periódicas.

**Fórmulas.** Criterio de enlace (`layers.bonds`, con `ase.neighborlist.neighbor_list` y radios por átomo $r_i + \mathrm{tol}/2$):

$$
d_{ij} \le r^{\mathrm{cov}}_i + r^{\mathrm{cov}}_j + \mathrm{tol}, \qquad \mathrm{tol} = 0.45\ \text{Å (por omisión)}
$$

Dimensionalidad (`layers._components_and_dim`): $\dim = \operatorname{rango}\{\mathbf{d}\}$ con $\mathbf{d} = \mathbf{o}_a + \mathbf{S}_{ab} - \mathbf{o}_b \ne 0$. Espaciados (`layers.analyze`):

$$
d_{\mathrm{basal}} = \frac{P}{n_{\mathrm{capas}}}, \qquad P = |\mathbf{a}_{\mathrm{apil}}\cdot\hat{\mathbf{n}}|, \qquad
\mathrm{hueco} = \min_k\left(z^{\mathrm{inf}}_{k+1} - z^{\mathrm{sup}}_k\right)
$$

Reflexión basal (`layers.report`), Bragg: $2\theta = 2\arcsin\!\left(\lambda/(2 d_{\mathrm{basal}}/l)\right)$ para $l = 1, 2, 3$.

- $\hat{\mathbf{n}}$: normal unitaria al plano de los dos vectores de celda no apilados; $z^{\mathrm{sup/inf}}_k$: centro ± grosor/2 de la capa $k$ (sin radios de van der Waals); $\lambda$ en Å.

**Cómo lo calcula Olla-DFT.**
1. `qekit/cli.py: _cmd_layers` → `qekit/core/layers.py: analyze` (`--tol` 0.45 Å).
2. Enlaces con ASE; componentes conexas y rango con `np.linalg.matrix_rank` sobre los vectores de cierre.
3. Eje de apilamiento: la dirección fraccionaria fuera del plano generado por los vectores de cierre (SVD), tomando el vector de celda con mayor componente fuera del plano.
4. Cada capa se reconstruye contigua (BFS con desplazamientos cartesianos) para medir centro y grosor sin cortes de celda.
5. `layers.report` imprime capas, $d$, hueco, periodo y las reflexiones 00l con la λ de `--wavelength` (por omisión CuKα = 1.54184 Å, `xrd.wavelength_value`), rotuladas con el nombre real de la radiación (`xrd.wavelength_name`: "Cu Kα", "Mo Kα1", o "λ dada" si se pasó un número).
6. Con `--slab ARCHIVO`, `layers.make_slab` aísla la primera capa: la desenrolla a lo largo del eje de apilamiento (imagen mínima en fraccionarias respecto del primer átomo, para que una capa que cruce la frontera de la celda no quede partida), sustituye el vector de apilamiento por la normal con altura grosor + `--vacuum` (20 Å), la centra y `structure.convert` la escribe.

**De dónde sale cada dato.**

| Dato | Origen | Detalle |
|---|---|---|
| Radios covalentes | `ase.data.covalent_radii` | biblioteca |
| Lista de vecinos periódica | `ase.neighborlist.neighbor_list("ijS")` | biblioteca |
| Longitud de onda | `xrd.WAVELENGTHS` o valor en Å | CuKa 1.54184, MoKa 0.71073, CoKa 1.79026, … |
| Tolerancia | `--tol` | 0.45 Å (`layers.DEFAULT_TOL`) |

**Límites y trampas.** Ningún cálculo de QE: es geometría pura. Si no hay componentes 2D: "No se detectaron capas … puedes probar con --tol menor". El eje de apilamiento y la normal se calculan con la PRIMERA capa; si hay capas con orientaciones distintas no lo verá. `make_slab` solo desenrolla a lo largo del eje de apilamiento (no en el plano), que es lo único que afecta al centrado y al grosor.

**Referencias.** M. Ashton, J. Paul, S. B. Sinnott y R. G. Hennig, *Phys. Rev. Lett.* 118, 106101 (2017), DOI 10.1103/PhysRevLett.118.106101 (criterio topológico de dimensionalidad). W. H. Bragg y W. L. Bragg, *Proc. R. Soc. Lond. A* 88, 428 (1913).

---

### `olla-dft xrd` — Difractograma de polvos simulado

**Qué responde.** ¿Dónde salen los picos de difracción de rayos X de polvo de esta estructura, con qué intensidad relativa y con qué índices hkl? ¿Se parece al difractograma medido?

**Fundamento para no expertos.** Un cristal difracta rayos X en direcciones fijadas por la ley de Bragg: cada familia de planos con espaciado $d$ produce un pico a un ángulo $2\theta$. La intensidad depende de cómo interfieren las ondas dispersadas por cada átomo de la celda (factor de estructura), de cuánto dispersa cada elemento (factor de dispersión atómica, que decae con el ángulo), y de factores geométricos del experimento (Lorentz–polarización). Cristalitos pequeños ensanchan los picos (Scherrer). Olla-DFT calcula todo eso y superpone, si se le da, un difractograma experimental, para ver a simple vista si el modelo estructural es el correcto.

**Fórmulas.** (`xrd.compute`) Para cada $hkl$ con $|\mathbf{g}| = 1/d$ dentro de la esfera accesible:

$$
\sin\theta = \frac{\lambda |\mathbf{g}|}{2}, \qquad s^2 = \left(\frac{\sin\theta}{\lambda}\right)^2 = \frac{|\mathbf{g}|^2}{4}
$$

$$
f(s) = Z - 41.78214\, s^2 \sum_{i=1}^{4} a_i\, e^{-b_i s^2}, \qquad f \to f\, e^{-B_{\mathrm{iso}} s^2}
$$

$$
F(hkl) = \sum_j f_j(s)\, e^{2\pi i\,(hkl)\cdot\mathbf{r}_j}, \qquad
I \propto |F|^2 \cdot \frac{1+\cos^2 2\theta}{\sin^2\theta\,\cos\theta}
$$

- $\lambda$: Å; $\mathbf{g} = (hkl)\,\mathbf{B}$ con $\mathbf{B} = (\mathbf{A}^{-1})^T$ sin $2\pi$; $Z$: número atómico; $a_i, b_i$: coeficientes analíticos (archivo de datos tomado de pymatgen, valores de las *International Tables*); $B_{\mathrm{iso}}$: factor de temperatura global en Å² (`--biso`, 0 por omisión); $\mathbf{r}_j$: posiciones fraccionarias.
- Multiplicidad: sale sola al enumerar TODOS los hkl y fusionar los que coinciden en $2\theta$ (tolerancia 0.02°). Intensidades normalizadas a 100; se descartan picos < 0.1.

Perfil (`xrd.broaden`), pseudo-Voigt con $\eta = 0.5$ y anchura $w$ (FWHM en ° 2θ):

$$
y(x) = \sum_p I_p\left[(1-\eta)\, e^{-\frac{(x-x_p)^2}{2\sigma^2}} + \eta\,\frac{1}{1+\left(\frac{x-x_p}{w/2}\right)^2}\right], \quad \sigma = \frac{w}{2\sqrt{2\ln 2}}, \quad
w_{\mathrm{Scherrer}} = \frac{K\lambda}{L\cos\theta},\; K = 0.9
$$

- $L$: tamaño de cristalito (`--size` en nm, convertido a Å); sin `--size`, $w$ = `--fwhm` (0.15°).

**Cómo lo calcula Olla-DFT.**
1. `qekit/cli.py: _cmd_xrd` → `qekit/modules/xrd.py: compute`. Con `--basis conventional` (por omisión) la celda se estandariza a la convencional (`structure.conventional`) para que los hkl coincidan con las fichas PDF; `input` usa la celda tal cual.
2. Enumera hkl en la caja $|h_i| \le \lceil g_{\max}/|\mathbf{b}_i|\rceil + 1$, filtra $g_{\min} \le |\mathbf{g}| \le g_{\max}$ (rango `--tt-min` 5°, `--tt-max` 70°).
3. Factores $f_j(s)$ de `qekit/data/atomic_scattering_params.json` (`xrd.scattering_params`), factor de estructura vectorizado, LP, descarte de reflexiones extinguidas ($I < 10^{-8} I_{\max}$), fusión por 2θ y etiqueta hkl "más legible" (`Peak.label`, orientada por Friedel).
4. `xrd.broaden` genera el perfil continuo (paso 0.02°).
5. `xrd.read_experimental` lee `--exp` (dos columnas, ≥ 10 filas; resta el mínimo y normaliza a 100).
6. `xrd.report` (12 picos más intensos), `xrd.export` (`XRD.dat`, `XRD_HKL.dat`), `xrd.plot` (experimental desplazado +105 arriba del simulado), `--suite` (JSON de intercambio).

**De dónde sale cada dato.**

| Dato | Origen | Detalle |
|---|---|---|
| Coeficientes $a_i, b_i$ | `qekit/data/atomic_scattering_params.json` | de pymatgen (MIT), *International Tables* |
| $Z$ | `ase.data.atomic_numbers` | biblioteca |
| Longitud de onda | `xrd.WAVELENGTHS` | CuKa 1.54184 Å, CuKa1 1.54056, MoKa 0.71073, CoKa 1.79026, FeKa 1.93735, CrKa 2.29100, AgKa 0.56087 |
| Constante de Scherrer | `xrd.SCHERRER_K` | 0.9 |
| Celda convencional | spglib (`structure.conventional`) | si falla, se usa la de entrada |

**Límites y trampas.** No hay refinamiento Rietveld ni factor R: la comparación con `--exp` es puramente visual (superposición). No aplica factor de absorción, ni orientación preferente, ni corrección de dispersión anómala, ni doblete Kα1/Kα2 (una sola λ). El factor de temperatura es un único $B$ isotrópico para todos los átomos. El ensanchamiento de Scherrer ignora la deformación (strain broadening). La fórmula de $f(s)$ es la parametrización de pymatgen (misma que la de las *International Tables* en su forma $Z - 41.78214 s^2\sum a_i e^{-b_i s^2}$), válida para rayos X. Con `--basis input` sobre una primitiva, "los hkl NO son los de la ficha PDF".

**Referencias.** *International Tables for Crystallography*, Vol. C (factores de dispersión). P. Scherrer, *Nachr. Ges. Wiss. Göttingen* 2, 98 (1918). S. P. Ong et al., *Comput. Mater. Sci.* 68, 314 (2013), DOI 10.1016/j.commatsci.2012.10.028 (pymatgen, origen de los coeficientes). B. E. Warren, *X-ray Diffraction*, Dover (1990).

---

### `olla-dft exfoliate` — Energía de exfoliación

**Qué responde.** ¿Cuánto cuesta separar una capa del cristal laminar, en J/m² (y meV/Å², meV/átomo)? ¿Es exfoliable?

**Fundamento para no expertos.** Se compara la energía del cristal (por capa) con la de una monocapa aislada en vacío. La diferencia por unidad de área es la energía de exfoliación; los laminares típicos están entre 0.2 y 0.6 J/m² (grafito ≈ 0.35 J/m² experimental). La cohesión entre capas es sobre todo dispersión de van der Waals, que LDA y PBE describen mal: sin corrección de dispersión el número no es comparable con el experimento, y el módulo lo dice.

**Fórmulas.** (`exfoliate.report_result`)

$$
E_{\mathrm{exf}} = \frac{E_{\mathrm{mono}} - E_{\mathrm{bulk}}/N_{\mathrm{capas}}}{A}
$$

- $E_{\mathrm{mono}}$, $E_{\mathrm{bulk}}$: energías totales en eV; $N_{\mathrm{capas}}$: capas por celda de bulk detectadas por `layers.analyze`; $A = |\mathbf{a}_i\times\mathbf{a}_j|$ de los dos vectores no apilados, en Å². Conversión con `EV_A2_TO_J_M2 = 16.02176634`.

**Cómo lo calcula Olla-DFT.**
1. `qekit/cli.py: _cmd_exfoliate` → `qekit/modules/exfoliate.py: prepare`: `layers.analyze(atoms, tol)` (`--tol` 0.45 Å); sin capas, error de uso.
2. `layers.make_slab` construye la monocapa (primera capa) con `--vacuum` 20 Å.
3. Malla k del bulk por `sweep.default_grid`; la de la monocapa es la misma en el plano y 1 en el eje de apilamiento.
4. Dos `scf` (`bulk/pw.in`, `monocapa/pw.in`; `relax` para la monocapa con `--relax-slab`), ambos con el mismo `--vdw` (`grimme-d2`, `grimme-d3`, `DFT-D`, `ts-vdw`, `xdm`, `mbd`).
5. `exfoliate.collect` lee `etot` de ambos XML; `exfoliate.report_result` imprime el resultado.

**De dónde sale cada dato.**

| Dato | Origen | Detalle |
|---|---|---|
| $E_{\mathrm{bulk}}$, $E_{\mathrm{mono}}$ | XML de pw.x (`etot`) | `qeout.read_xml` |
| $N_{\mathrm{capas}}$, eje de apilamiento | `layers.analyze` | conectividad por enlaces |
| Área $A$ | celda de bulk | dos vectores no apilados |
| eV/Å² → J/m² | `exfoliate.EV_A2_TO_J_M2` | 16.02176634 |

**Límites y trampas.** Supone que todas las capas de la celda son equivalentes (divide $E_{\mathrm{bulk}}$ entre $N_{\mathrm{capas}}$ y usa solo la primera). Sin `--vdw`: "SIN corrección de van der Waals: PBE apenas liga las capas y LDA liga por cancelación de errores". Con pseudos que parecen LDA (nombre con `pz`, `lda`, `pw92`) y Grimme: "combinarlas cuenta la dispersión dos veces". Si sale negativa: "Casi siempre significa que falta la corrección vdW o que algún cálculo no está bien convergido". No hay relajación del bulk. No escribe archivos `.dat` (solo el reporte en pantalla).

**Referencias.** J. H. Jung, C.-H. Park y J. Ihm, *Nano Lett.* 18, 2759 (2018), DOI 10.1021/acs.nanolett.7b04201 (energía de exfoliación vs. de enlace interlaminar). S. Grimme, *J. Comput. Chem.* 27, 1787 (2006), DOI 10.1002/jcc.20495 (D2); S. Grimme et al., *J. Chem. Phys.* 132, 154104 (2010), DOI 10.1063/1.3382344 (D3).

---

### `olla-dft phonons` — Fonones DFPT: dispersión, DOS, termodinámica, IR, Raman y temperatura electrónica

**Qué responde.** ¿Cuáles son las frecuencias de vibración del cristal (en Γ o a lo largo de la zona de Brillouin), es dinámicamente estable (sin frecuencias imaginarias), qué energía de punto cero, energía libre, entropía y calor específico armónicos tiene, qué modos son activos en IR y Raman (`--raman`) y, con `--tscan`, se estabiliza un modo blando al subir la temperatura electrónica?

**Fundamento para no expertos.** Los átomos de un cristal vibran alrededor de sus posiciones de equilibrio como si estuvieran unidos por muelles. La teoría de perturbaciones de la densidad (DFPT, lo que hace `ph.x`) calcula la rigidez de esos muelles (las constantes de fuerza) a partir de la densidad electrónica, sin desplazar átomos a mano. De ahí salen las frecuencias de todas las ondas de vibración (fonones). Una frecuencia "imaginaria" (negativa en la salida) significa que la estructura no está en un mínimo: o no se relajó bien o es inestable. Con la densidad de estados de fonones se calcula la termodinámica armónica: incluso a 0 K los átomos vibran (energía de punto cero) y al calentar se poblan más modos (entropía, calor específico). Los modos en Γ son los que ve un espectrómetro de infrarrojo o Raman.

**Fórmulas.** Termodinámica armónica por celda desde la DOS $g(\omega)$ (`phonons.thermodynamics`), con $\epsilon = \hbar\omega$ en eV, $x = \epsilon/k_BT$ (acotado a 500), $n = 1/(e^x - 1)$, y $g$ renormalizada a $\int g\,d\omega = 3N_{\mathrm{at}}$ (solo $\omega > 1$ cm⁻¹):

$$
E_{\mathrm{ZPE}} = \int \tfrac{1}{2}\epsilon\, g\, d\omega, \qquad
F(T) = E_{\mathrm{ZPE}} + k_B T \int \ln\!\left(1 - e^{-x}\right) g\, d\omega
$$

$$
U(T) = \int \left(\tfrac{1}{2} + n\right)\epsilon\, g\, d\omega, \qquad
C_v(T) = k_B \int x^2 e^{x} n^2\, g\, d\omega, \qquad S(T) = \frac{U - F}{T}
$$

- $k_B$ = `KB_EV` = 8.617333262e-5 eV/K; cm⁻¹ → eV con `CM1_TO_EV` = 1.239841984e-4; cm⁻¹ → THz con `CM1_TO_THZ` = 0.0299792458. Integrales por regla del trapecio. $T$ = 0…1000 K en pasos de 10.

Espectro Raman Stokes (`phonons.raman_spectrum`) a partir de la actividad $A$ (Å⁴/amu) de `dynmat.x`:

$$
I(\omega) \propto \frac{(\omega_L - \omega)^4}{\omega}\,[n(\omega,T)+1]\,A(\omega), \qquad \omega_L = \frac{10^7}{\lambda_{\mathrm{láser}}[\mathrm{nm}]}\ \mathrm{cm^{-1}}
$$

convolucionado con lorentzianas de FWHM 5 cm⁻¹ a $T$ = 300 K; se excluyen $\omega \le 1$ cm⁻¹. Temperatura electrónica (`tphonons.degauss_de_T`):

$$
\mathrm{degauss} = k_B T, \qquad k_B = 6.333621\times10^{-6}\ \mathrm{Ry/K}, \qquad \text{smearing = fermi-dirac}
$$

Temperatura de estabilización (`tphonons.temperatura_de_estabilizacion`): interpolación lineal de $T$ donde el modo más blando (mínimo de las frecuencias con $|\omega| > 10$ cm⁻¹) cruza de negativo a no negativo.

**Cómo lo calcula Olla-DFT.**
1. `qekit/cli.py: _cmd_phonons` → `qekit/modules/phonons.py: prepare`. La estructura se lleva a la primitiva estandarizada (`structure.primitive`). Escribe `scf.in` (`conv_thr = 1e-12`, malla k por `kspacing`), `ph.in` (`tr2_ph = 1e-14`, `fildyn='dyn'`; `epsil=.true.` si `--insulator` o `--raman`; `lraman=.true.` y `trans=.true.` con `--raman`; `ldisp` con malla `nq` = `--qgrid` o `kgrid_from_spacing(atoms, 0.6)`).
2. Modo malla: `q2r.in` (`zasr='simple'`, `flfrc='fuerzas.fc'`), `matdyn_band.in` (camino de seekpath por `kpoints.get_kpath`, 30 puntos por tramo, `q_in_band_form`, `q_in_cryst_coord`), `matdyn_dos.in` (`dos=.true.`, malla 12×12×12, `fldos='fonones.dos'`). Modo Γ (`--gamma` o `--raman`): `dynmat.in` (`asr='simple'`).
3. `--raman` exige pseudopotenciales de norma conservada (`p["type"] == "NC"`), si no: error de uso.
4. `--run`: `runner.run_all` corre `pw.x`; `phonons.run_chain` ejecuta `ph.x` → (`dynmat.x` | `q2r.x` → `matdyn.x` ×2), saltando pasos cuyo `.out` ya diga `JOB DONE`.
5. `phonons.collect`: en Γ lee la tabla `# mode [cm-1] [THz] IR [Raman depol]` de `dynmat.out` (`read_dynmat_table`); en malla lee `bandas.freq` (`_read_flfrq`, formato `&plot`, q en cartesianas 2π/alat) y `fonones.dos`.
6. `phonons.report` / `report_gamma_activities` (regla de exclusión mutua, despolarización 0.75), `phonons.thermodynamics`, `phonons.export` (`FONONES_GAMMA.dat` o `FONONES_BANDAS.dat`, `FONONES_DOS.dat`, `FONONES_TERMO.dat`), y `phonons.plot` solo si `phonons.has_dispersion(run)` (hay `band_freqs` y `qdist`, es decir, no es una corrida en Γ).
7. `--tscan T1,T2,...` → `qekit/cli.py: _cmd_phonons_tscan` → `qekit/modules/tphonons.py: prepare`: una cadena completa por temperatura en `T00300/` etc., con `insulator=False`, `smearing='fermi-dirac'` y `degauss = k_B T`; `tphonons.collect`, `report` (tabla de modos imaginarios por T, monotonía, $T_{\mathrm{est}}$), `export` (`FONONES_T.dat`), `plot`.

**De dónde sale cada dato.**

| Dato | Origen | Detalle |
|---|---|---|
| Frecuencias en Γ, IR, Raman, depol | `dynmat.out` (tabla `# mode`) | `phonons.read_dynmat_table`; IR en (D/Å)²/amu, Raman en Å⁴/amu |
| Dispersión | `bandas.freq` de `matdyn.x` | `phonons._read_flfrq` |
| DOS de fonones | `fonones.dos` de `matdyn.x` | `np.loadtxt`, estados/cm⁻¹ |
| Camino de alta simetría | seekpath (`kpoints.get_kpath`) | etiquetas y discontinuidades |
| $k_B$, cm⁻¹→eV, cm⁻¹→THz | `phonons.KB_EV`, `CM1_TO_EV`, `CM1_TO_THZ` | CODATA |
| $k_B$ en Ry/K | `tphonons.KB_RY` | 6.333621e-6 |
| Umbral de imaginaria | literal −5 cm⁻¹ (`phonons.report`), `tphonons.UMBRAL_IMAGINARIO` = 10 | ruido numérico por debajo |

**Límites y trampas.** Es armónico: no hay expansión térmica ni anarmonicidad (para eso, `qha`). El docstring de `prepare` fija `insulator=True` por omisión, pero la CLI pasa `args.insulator`, que es `False` salvo `--insulator`: por omisión el scf lleva smearing y NO se activa `epsil` (sin separación LO–TO). Avisa: "hay frecuencias imaginarias (negativas). O la estructura no está relajada, o es inestable en Γ" y "la estructura debe estar relajada (vc-relax) con estos mismos cutoffs". La escala absoluta de la DOS de `matdyn` depende de la malla; se renormaliza a $3N$ en la termodinámica. Con `--raman` (que fuerza el modo Γ aunque no se pase `--gamma`) no se dibuja dispersión: la CLI decide con `has_dispersion(run)`, no con la bandera. La termodinámica ignora $\omega \le 1$ cm⁻¹. En `--tscan`, "esto es temperatura ELECTRÓNICA. Los iones siguen estando quietos"; si el número de imaginarias no baja monótonamente: "Suele querer decir que falta convergencia en la malla de k". Solo `smearing='fermi-dirac'` corresponde a una temperatura real.

**Referencias.** S. Baroni, S. de Gironcoli, A. Dal Corso y P. Giannozzi, *Rev. Mod. Phys.* 73, 515 (2001), DOI 10.1103/RevModPhys.73.515 (DFPT). M. Lazzeri y F. Mauri, *Phys. Rev. Lett.* 90, 036401 (2003), DOI 10.1103/PhysRevLett.90.036401 (Raman por 2n+1). D. Porezag y M. R. Pederson, *Phys. Rev. B* 54, 7830 (1996) (intensidades Raman). A. A. Maradudin et al., *Theory of Lattice Dynamics in the Harmonic Approximation*, Academic Press (1971).

---

### `olla-dft qha` — Aproximación cuasi-armónica

**Qué responde.** ¿Cómo se dilata el cristal con la temperatura ($V(T)$, $\alpha(T)$, $a(T)$), cuánto vale el parámetro de Grüneisen, $C_p$ frente a $C_v$ y $B(T)$?

**Fundamento para no expertos.** La aproximación armónica no da dilatación: si las frecuencias no dependen del volumen, el mínimo de la energía libre no se mueve. La QHA mantiene los modos armónicos pero deja que sus frecuencias cambien con el VOLUMEN. A cada temperatura se suma la energía estática $E(V)$ y la energía libre vibracional $F_{\mathrm{vib}}(V,T)$, y se busca el volumen que minimiza la suma. Al calentar, la energía libre vibracional favorece volúmenes mayores (los modos se ablandan) y el mínimo se desplaza: eso es la dilatación térmica.

**Fórmulas.** (`qha.f_vib`, `qha.cv_modos`, `qha.run`) Para cada volumen $V_i$ con sus modos $\omega_k$ (cm⁻¹, solo $\omega > 1$), $\epsilon_k = \hbar\omega_k$:

$$
F(V_i,T) = E(V_i) + \frac{1}{N_{\mathrm{celdas}}}\left[\sum_k \tfrac{1}{2}\epsilon_k + k_B T \sum_k \ln\!\left(1 - e^{-\epsilon_k/k_BT}\right)\right]
$$

Mínimo por parábola local (hasta 2 puntos a cada lado del mínimo muestreado): $V(T) = -b/2a$ (recortado al rango), $B(T) = 2a\,V(T)\times 160.21766208$ GPa.

$$
\alpha(T) = \frac{1}{V}\frac{dV}{dT}\ (\texttt{np.gradient}), \qquad
C_v = k_B\sum_k x_k^2 \frac{e^{x_k}}{(e^{x_k}-1)^2}\ \text{(interpolada en } V(T)), \qquad
C_p = C_v + \alpha^2 B V T
$$

$$
\gamma = -\frac{d\ln\langle\omega\rangle}{d\ln V}\ \text{(recta sobre } \ln V\text{, } \ln\bar\omega\text{)}, \qquad
a(T) = \begin{cases} (V \cdot V_{\mathrm{conv}}/V_{\mathrm{prim}})^{1/3} & \text{con } \texttt{--structure} \\ V_{\mathrm{prim}}^{1/3} & \text{sin ella} \end{cases}\ (\texttt{--cubic})
$$

- $E$ en eV, $V$ en Å³, $C_v$, $C_p$ en meV/K por celda, $\alpha$ en K⁻¹, $B$ en GPa; `KB_EV` = 8.617333262e-5, `CM1_EV` = 1.239841984e-4.

**Cómo lo calcula Olla-DFT.**
1. `qekit/cli.py: _cmd_qha` lee una TABLA (`data`): columnas $V$ (Å³), $E$ (eV), $\omega_1, \omega_2, \ldots$ (cm⁻¹) por volumen; los valores ≤ −1000 se descartan como relleno.
2. Con `--structure CIF`, `qha.factor_convencional` cuenta cuántas primitivas caben en la convencional ($N_{\mathrm{conv}}/N_{\mathrm{prim}}$ con spglib: 4 en fcc/diamante, 2 en bcc, 1 en cúbica simple) y `qha.es_cubico` (grupo ≥ 195) activa `--cubic` automáticamente.
3. `qekit/modules/qha.py: run` con $T$ = 0 … `--tmax` (1000) en pasos `--dt` (5), `--natoms` (1), `--cells` (celdas primitivas por supercelda de los modos, 1), `--cubic`, `factor_conv`.
4. Avisos si hay < 4 volúmenes o frecuencias < −5 cm⁻¹ en algún volumen.
5. `qha.report` (a `--temp` 300 K; $a(T)$ rotulado "parámetro de red (celda convencional)" o "V_prim^(1/3) (NO es el parámetro de red convencional)" según `QHAResult.a_convencional`), `qha.export` (`QHA.dat`), `qha.plot` (V, α, $C_v$/$C_p$).

**De dónde sale cada dato.**

| Dato | Origen | Detalle |
|---|---|---|
| $E(V)$ y $\omega_k(V)$ | tabla del usuario (`data`) | de `eos` + `phonons` (o `mlip phonons`) por volumen; Olla-DFT no la genera |
| $V_{\mathrm{conv}}/V_{\mathrm{prim}}$ y cubicidad | `--structure` vía spglib | `qha.factor_convencional`, `qha.es_cubico` |
| Constantes | `qha.KB_EV`, `CM1_EV`, `EV_A3_GPA` | CODATA; 160.21766208 |

**Límites y trampas.** No lanza ningún cálculo de QE: recibe la tabla. Usa frecuencias discretas (un conjunto de modos por volumen), no una DOS: la termodinámica se hace sobre la lista que se le pasa, así que la calidad depende de la supercelda/malla de esos modos. El Grüneisen es un promedio sobre la frecuencia media, no modo a modo. La QHA "vale hasta ~la mitad de la temperatura de fusión". Sin `--structure`, $a(T)$ es solo $V_{\mathrm{prim}}^{1/3}$ y el reporte avisa: "En fcc, bcc o diamante eso NO es el parámetro de red convencional (difieren en 4^(1/3) o 2^(1/3)). Pasa la estructura con --structure". Con una sola temperatura, $\alpha$ es NaN y se avisa.

**Referencias.** A. Togo, L. Chaput, I. Tanaka y G. Hug, *Phys. Rev. B* 81, 174301 (2010), DOI 10.1103/PhysRevB.81.174301. G. Grimvall, *Thermophysical Properties of Materials*, North-Holland (1999). P. Pavone et al., *Phys. Rev. B* 48, 3156 (1993) (Si, expansión negativa).

---

### `olla-dft derived` — Debye, velocidades del sonido y Slack desde las $C_{ij}$

**Qué responde.** A partir de una matriz elástica ya calculada: ¿cuánto valen la densidad, las velocidades del sonido, la temperatura de Debye elástica, la razón de Poisson, un Grüneisen aproximado y una estimación de la conductividad térmica de red?

**Fundamento para no expertos.** La rigidez de un sólido determina a qué velocidad viajan las ondas de sonido por él, y esa velocidad fija cuál es la vibración más rápida posible; la temperatura de Debye es esa frecuencia máxima expresada en kelvin. Con ella y un parámetro de anarmonicidad (Grüneisen) el modelo de Slack estima cuánto calor conduce la red. Todo es post-proceso: no cuesta ningún cálculo nuevo.

**Fórmulas.** (`derived.density`, `sound_velocities`, `debye_from_velocity`, `poisson_ratio`, `gruneisen_from_poisson`, `slack`, `cubic_directional`)

$$
\rho = \frac{\sum_i m_i}{V}, \qquad
v_l = \sqrt{\frac{B + \tfrac{4}{3}G}{\rho}}, \qquad v_t = \sqrt{\frac{G}{\rho}}, \qquad
v_m = \left[\frac{1}{3}\left(\frac{2}{v_t^3} + \frac{1}{v_l^3}\right)\right]^{-1/3}
$$

$$
\Theta_D = \frac{\hbar}{k_B}\left(6\pi^2 n\right)^{1/3} v_m, \qquad
\nu = \frac{3B - 2G}{2(3B+G)}, \qquad
\gamma = \frac{3(1+\nu)}{2(2-3\nu)}
$$

$$
\kappa_L = A\,\frac{\bar M\,\Theta_D^3\,\delta}{\gamma^2\, n_{\mathrm{at}}^{2/3}\, T}, \qquad
A = \frac{3.1\times10^{-6}}{1 - 0.514/\gamma + 0.228/\gamma^2}, \qquad
v_L^{[100]} = \sqrt{C_{11}/\rho},\ v_T^{[100]} = \sqrt{C_{44}/\rho}
$$

- $B$, $G$: promedios de Hill en GPa (×10⁹ a Pa); $\rho$ en kg/m³ (masas en amu × 1.66053906660e-27, $V$ en Å³ × 1e-30); $n$: átomos por m³; $\hbar$ = 1.054571817e-34 J·s, $k_B$ = 1.380649e-23 J/K; $\bar M$: masa media en amu; $\delta = (V/n_{\mathrm{at}})^{1/3}$ en Å; $T$ en K (`--temp`, 300); $\kappa_L$ en W/(m·K).

**Cómo lo calcula Olla-DFT.**
1. `qekit/cli.py: _cmd_derived` carga la estructura (masas y volumen) y `--cij` (`ELASTIC_C.dat` de `elastic`, matriz 6×6).
2. `elastic.moduli` → $B_H$, $G_H$; `qekit/modules/derived.py: analyze` calcula todo lo anterior.
3. `derived.cubic_directional` imprime $v_L$, $v_T$ a lo largo de [100] solo si la estructura es cúbica según spglib (`elastic.crystal_family`) o el tensor tiene forma cúbica (`derived.is_cubic_tensor`: $C_{11}=C_{22}=C_{33}$, $C_{12}=C_{13}=C_{23}$, $C_{44}=C_{55}=C_{66}$ y ceros fuera, con tolerancia 5 % o 2 GPa).
4. `derived.report`; `derived.export` escribe `DERIVED.dat`.

**De dónde sale cada dato.**

| Dato | Origen | Detalle |
|---|---|---|
| $C_{ij}$ | `ELASTIC_C.dat` (`--cij`) | `np.loadtxt`, GPa |
| Masas y volumen | estructura (ASE `get_masses`, `get_volume`) | amu, Å³ |
| $\hbar$, $k_B$, amu | `derived.HBAR`, `KB`, `AMU` | CODATA 2018 |
| Prefactor de Slack | literal 3.1e-6 y corrección $(1 - 0.514/\gamma + 0.228/\gamma^2)$ | Slack / Julian |

**Límites y trampas.** La $\Theta_D$ es la ELÁSTICA (límite acústico de baja temperatura): "La que sale de la DOS de fonones usa todo el espectro y da otro número; no son la misma cantidad" (existe `derived.debye_from_dos`, $\Theta_D = (\hbar/k_B)\sqrt{5\langle\omega^2\rangle/3}$, usada por `crosscheck`, no por este comando). El Grüneisen "viene de una correlación empírica con la razón de Poisson" (Belomestnykh) y Slack "es una estimación de orden de magnitud". Poisson negativo: aviso de material auxético. La κ de Slack se rotula con la temperatura realmente usada (`Termoelastico.T`, clave `kappa_Slack_<T>K` en `DERIVED.dat`). Si $G \le 0$ no hay velocidades.

**Referencias.** O. L. Anderson, *J. Phys. Chem. Solids* 24, 909 (1963), DOI 10.1016/0022-3697(63)90067-2 ($\Theta_D$ elástica). G. A. Slack, *Solid State Phys.* 34, 1 (1979); D. T. Morelli y G. A. Slack, en *High Thermal Conductivity Materials*, Springer (2006) (prefactor con corrección de Julian). V. N. Belomestnykh y E. P. Tesleva, *Tech. Phys.* 49, 1098 (2004) (Grüneisen–Poisson).

---

### `olla-dft thermochem` — ZPE, entropía y energía libre

**Qué responde.** ¿Cuánto hay que sumar a una energía DFT (a 0 K, sin vibraciones) para obtener una energía libre $G(T,p)$ comparable con el experimento, para un sólido, un adsorbato, un gas ideal o un estado de transición?

**Fundamento para no expertos.** Una energía DFT es electrónica y a 0 K. Lo que se mide es una energía libre a la temperatura y presión del laboratorio. Entre ambas hay tres términos: la energía de punto cero (los modos vibran incluso a 0 K), la corrección entálpica (los modos se poblan al calentar) y el término entrópico $-TS$, que en una molécula gaseosa incluye las entropías de traslación (Sackur–Tetrode) y rotación (rotor rígido) y puede valer del orden de 1 eV a 500 K. Olvidarlo puede cambiar el signo de una energía de adsorción.

**Fórmulas.** (`thermochem.zpe`, `H_vib`, `S_vib`, `Cv_vib`, `S_traslacional`, `S_rotacional`, `corregir`) Con $\epsilon_k = h c\,\tilde\nu_k$, $x_k = \epsilon_k/k_BT$ (acotado a 500):

$$
E_{\mathrm{ZPE}} = \tfrac{1}{2}\sum_k \epsilon_k, \quad
H_{\mathrm{vib}} = \sum_k \frac{\epsilon_k}{e^{x_k}-1}, \quad
S_{\mathrm{vib}} = k_B\sum_k\left[\frac{x_k}{e^{x_k}-1} - \ln\!\left(1-e^{-x_k}\right)\right], \quad
C_v = k_B\sum_k \frac{x_k^2 e^{x_k}}{(e^{x_k}-1)^2}
$$

$$
S_{\mathrm{tras}} = k_B\left[\ln\!\left(\frac{V}{\Lambda^3}\right) + \tfrac{5}{2}\right], \quad V = \frac{k_BT}{p}, \quad \Lambda = \frac{h}{\sqrt{2\pi m k_B T}}
$$

$$
S_{\mathrm{rot}}^{\mathrm{lineal}} = k_B\left[\ln\frac{T}{\sigma\Theta_r} + 1\right], \qquad
S_{\mathrm{rot}}^{\mathrm{no\,lineal}} = k_B\left[\tfrac{1}{2}\ln\frac{\pi T^3}{\sigma^2\Theta_A\Theta_B\Theta_C} + \tfrac{3}{2}\right], \qquad \Theta_i = \frac{\hbar^2}{2 I_i k_B}
$$

$$
G - E_{\mathrm{DFT}} = E_{\mathrm{ZPE}} + H_{\mathrm{corr}} - TS, \qquad
H_{\mathrm{corr}}^{\mathrm{gas}} = H_{\mathrm{vib}} + \left(\tfrac{3}{2} + n_{\mathrm{rot}} + 1\right)k_BT, \qquad S_{\mathrm{elec}} = k_B\ln(\text{multiplicidad})
$$

- $\tilde\nu$ en cm⁻¹ (`C_CM` = 2.99792458e10 cm/s, `H_EVS` = 4.135667696e-15 eV·s); $m$: masa molecular (amu → kg); $p$ en Pa (`--pressure` en bar); $\sigma$: número de simetría (`--symmetry`, 1); $I_i$: momentos principales de inercia (amu·Å² → kg·m²); lineal si $I_1 < 10^{-3} I_3$; $n_{\mathrm{rot}}$ = 1 (lineal), 1.5 (no lineal), 0 (átomo). Energía de adsorción (`thermochem.adsorcion`): $E_{\mathrm{ads}} = E_{\mathrm{slab+ads}} - E_{\mathrm{slab}} - nE_{\mathrm{gas}}$; $G_{\mathrm{ads}} = E_{\mathrm{ads}} + G^{\mathrm{corr}}_{\mathrm{ads}} - n\,G^{\mathrm{corr}}_{\mathrm{gas}}$.

**Cómo lo calcula Olla-DFT.**
1. `qekit/cli.py: _cmd_thermochem` lee las frecuencias (`_leer_frecuencias`: archivo de una o varias columnas —última columna—, o lista inline) y, para gas, la estructura con `ase.io.read`.
2. `qekit/modules/thermochem.py: limpiar_frecuencias`: separa imaginarias ($\tilde\nu < -1$), descarta $|\tilde\nu| \le 1$ (traslaciones/rotaciones residuales), sube al piso `--floor` (p. ej. 100 cm⁻¹) las blandas y emite avisos según `--phase` (`solido`, `adsorbato`, `gas`, `transicion`).
3. `thermochem.corregir` suma los términos a `--temp` (298.15 K) y `--pressure` (1 bar), con `--multiplicity`.
4. `thermochem.report` (con `G(T)` si se da `--energy`); con `-o`, escribe `TERMOQUIMICA.txt`.

**De dónde sale cada dato.**

| Dato | Origen | Detalle |
|---|---|---|
| Frecuencias | archivo (`FONONES_GAMMA.dat`, columna final) o lista | `_leer_frecuencias` |
| Masas y geometría (gas) | `--structure` vía ASE | momentos de inercia |
| $h$, $k_B$, $c$, amu, $\hbar$ | `thermochem.H_EVS`, `KB_EV`, `C_CM`, `AMU_KG`, `HBAR_JS`, `KB_J` | CODATA |
| Piso de modos blandos | `--floor` (`PISO_BLANDO` = 100 solo como referencia) | sin `--floor` no se sube ninguno |

**Límites y trampas.** Verificado en los tests contra NIST (H₂O, N₂, CH₄ al 0.5 %). Estado de transición sin imaginaria o con más de una: aviso explícito. Imaginarias en un mínimo: "la estructura NO es un mínimo … Se excluyen de las sumas". Modos subidos: "es una CORRECCIÓN, no un cálculo: dilo si publicas". Gas con número de modos ≠ $3N-6$ (o $3N-5$): "cuentan doble con los términos traslacional y rotacional". En fase gas no se incluye anarmonicidad ni rotores internos; el sólido no lleva término $pV$. `adsorcion` no aplica correcciones vibracionales a la losa limpia (asume que no cambian). El comando `thermochem` de la CLI no expone `adsorcion` (se usa desde `adsorb`/API).

**Referencias.** D. A. McQuarrie, *Statistical Mechanics*, University Science Books (2000). C. J. Cramer, *Essentials of Computational Chemistry*, Wiley (2004), cap. 10. O. Sackur, *Ann. Phys.* 36, 958 (1911); H. Tetrode, *Ann. Phys.* 38, 434 (1912).

---

### `olla-dft md` — Análisis de una trayectoria de dinámica molecular

**Qué responde.** De una salida de `pw.x` con `calculation='md'`: ¿qué estructura tiene el sistema (g(r), números de coordinación), difunden los átomos (MSD y coeficiente de difusión $D$) y cuál es su espectro vibracional (VDOS) incluyendo temperatura y anarmonicidad?

**Fundamento para no expertos.** Una dinámica molecular es una "película" de los átomos moviéndose. Tres funciones resumen la película: la función de distribución radial $g(r)$ dice cuántos vecinos hay a cada distancia (su primer pico es la distancia de enlace, su área el número de coordinación); el desplazamiento cuadrático medio (MSD) dice si los átomos se alejan de donde estaban (si crece linealmente, difunden; si se aplana, solo vibran); y la transformada de Fourier de la autocorrelación de velocidades da las frecuencias a las que vibran. Antes de creerse ninguna hay que descartar el tramo inicial de equilibrado y comprobar que la trayectoria es bastante larga.

**Fórmulas.** (`dynamics.rdf`, `coordinacion`, `msd`, `difusion`, `vdos`) Con distancias por imagen mínima en coordenadas fraccionarias:

$$
g(r) = \frac{h(r)}{N_{\mathrm{pasos}}\,\frac{N(N-1)}{2}\,\frac{4\pi r^2\,\Delta r}{V}}, \qquad
g_{AB}(r) = \frac{h_{AB}(r)}{N_{\mathrm{pasos}}\,N_{\mathrm{pares}}\,\frac{4\pi r^2\Delta r}{V}}, \qquad
n_{\mathrm{coord}} = \int_0^{r_{\min}} 4\pi r^2 \rho\, g(r)\, dr
$$

$$
\mathrm{MSD}(\tau) = \left\langle |\mathbf{r}_i(t+\tau) - \mathbf{r}_i(t)|^2 \right\rangle_{i,t}, \qquad
\mathrm{MSD} = 6 D \tau + b \;\Rightarrow\; D = \frac{m}{6}\times 10^{-1}\ [\mathrm{cm^2/s}]
$$

$$
C(t) = \frac{\sum_{i,\alpha}\langle v_{i\alpha}(0)v_{i\alpha}(t)\rangle}{C(0)}\ (\text{por FFT}), \qquad
\mathrm{VDOS}(\tilde\nu) = \left|\mathcal{F}\{C(t)\,w_{\mathrm{Hann}}(t)\}\right|, \quad \tilde\nu = \frac{f[\mathrm{fs^{-1}}]\times 10^{15}}{c[\mathrm{cm/s}]}
$$

- $r_{\max}$ = mitad de la arista menor de la celda (o `--rmax` si es menor); `--bins` 200; $N_{\mathrm{pares}} = N_AN_B$ o $N_A(N_A-1)/2$; $\rho = N/V$; $r_{\min}$: primer mínimo local con $g<1$ tras el primer máximo; el ajuste de $D$ usa solo el tramo 20–80 % de los retardos (retardos hasta $n/2$), pendiente $m$ en Å²/fs; velocidades por `np.gradient` de las posiciones DESDOBLADAS (sin saltos periódicos), sin ponderar por masa.

**Cómo lo calcula Olla-DFT.**
1. `qekit/cli.py: _cmd_md` → `qekit/modules/dynamics.py: leer_md`: lee de `pw.out` (texto) los bloques `ATOMIC_POSITIONS` (alat, bohr, angstrom o crystal → Å), la celda de `a(i) = (...)`·alat o `CELL_PARAMETERS`, `temperature = … K`, `!    total energy` y `Time step = … femto-seconds`; descarta `--skip` pasos.
2. `dynamics.analizar`: `rdf`, `coordinacion` por par, `desdoblar` + `msd` (total y por especie), `difusion`, `vdos` (≥ 8 pasos), deriva de temperatura entre mitades.
3. `dynamics.report`, `dynamics.export` (`MD_RDF.dat`, `MD_MSD.dat`, `MD_VDOS.dat`, `MD.txt`), `dynamics.plot` (tres paneles).

**De dónde sale cada dato.**

| Dato | Origen | Detalle |
|---|---|---|
| Posiciones por paso | `pw.out` (`ATOMIC_POSITIONS`) | `dynamics._leer_marcos`; unidades detectadas |
| Celda | `pw.out` (`a(1..3) = (...)` × alat, o `CELL_PARAMETERS`) | se asume constante |
| Paso de tiempo | `pw.out` (`Time step = … a.u., X femto-seconds`) | 1 fs si no aparece |
| Temperatura y energía | `pw.out` (`temperature =`, `!    total energy`) | K; Ry → eV con 13.605693122994 |
| bohr → Å | `dynamics.BOHR_A` | 0.529177210903 |

**Límites y trampas.** Celda constante: no sirve para `vc-md`. Si $g(r)$ está vacía hasta el corte: "la celda es demasiado pequeña para sacar estructura de ella: haz una supercelda". Menos de 2 ps: "Sirve para ver la estructura, no para un coeficiente de difusión". $R^2 < 0.95$: "el MSD NO es lineal: no hay difusión, o falta tiempo". Deriva de temperatura > 15 %: "todavía está equilibrando: descarta más pasos con --skip". El número de coordinación se integra hasta el PRIMER MÍNIMO, "una convención, no una medida". La VDOS no pondera por masa (no es la DOS de fonones) y toma el módulo del espectro, no la parte real; su resolución es $1/(N_{\mathrm{pasos}}\,\Delta t)$. `KB_RY` en `dynamics.py` no se usa.

**Referencias.** M. P. Allen y D. J. Tildesley, *Computer Simulation of Liquids*, Oxford (2017). A. Einstein, *Ann. Phys.* 17, 549 (1905). J.-P. Hansen e I. R. McDonald, *Theory of Simple Liquids*, Academic Press (2013).

---

### `olla-dft kappa` — Conductividad térmica de red (fc3 + BTE con phono3py)

**Qué responde.** ¿Cuánto calor conduce la red cristalina, $\kappa_L(T)$ en W/(m·K), qué exponente sigue con la temperatura y qué recorridos libres medios de fonón llevan el calor (para saber si nanoestructurar sirve)?

**Fundamento para no expertos.** En un cristal perfectamente armónico un fonón viajaría para siempre y la conductividad sería infinita. Lo que la hace finita es que un fonón puede partirse en dos o dos fundirse en uno: eso lo permite el término cúbico de la energía (las constantes de fuerza de tercer orden, fc3). Se obtienen desplazando dos átomos a la vez en una supercelda y calculando las fuerzas; con ellas, la ecuación de Boltzmann de fonones (en aproximación de tiempo de relajación, RTA) da $\kappa$. Es caro porque el número de configuraciones crece deprisa con la supercelda. Olla-DFT admite calcular las fuerzas con `pw.x` (el cálculo de verdad) o con un potencial aprendido (MACE, etc.) para explorar.

**Fórmulas.** Las resuelve phono3py (`kappa.resolver`); Olla-DFT post-procesa:

$$
\kappa_L^{\alpha\beta} = \frac{1}{NV}\sum_\lambda C_\lambda\, v_\lambda^\alpha v_\lambda^\beta\, \tau_\lambda, \qquad \tau_\lambda = \frac{1}{2\Gamma_\lambda}, \qquad \Lambda_\lambda = |\mathbf{v}_\lambda|\,\tau_\lambda
$$

$$
\bar\kappa = \frac{\kappa_{xx}+\kappa_{yy}+\kappa_{zz}}{3}, \qquad
\kappa \propto T^{-n}\ (n \text{ por recta en } \ln\kappa\text{–}\ln T,\ T \ge 200\ \mathrm{K}), \qquad
\kappa_{\mathrm{acum}}(\Lambda) = \frac{\sum_{\lambda:\Lambda_\lambda<\Lambda} w_\lambda C_\lambda \tfrac{|\mathbf{v}_\lambda|^2}{3}\tau_\lambda}{\sum_\lambda w_\lambda C_\lambda \tfrac{|\mathbf{v}_\lambda|^2}{3}\tau_\lambda}
$$

- $\Gamma_\lambda$: anchura de línea (THz) de phono3py; $\mathbf{v}_\lambda$: velocidad de grupo (THz·Å); $C_\lambda$: capacidad calorífica modal; $w_\lambda$: peso del punto q; $\Lambda$ en Å (se reporta en nm). Se descartan los modos con $\Gamma = 0$ (acústicos en Γ).

**Cómo lo calcula Olla-DFT.**
1. `qekit/cli.py: _cmd_kappa` → `qekit/modules/kappa.py: preparar`: `Phono3py(..., supercell_matrix=--dim (2x2x2), phonon_supercell_matrix=--dim-fc2, primitive_matrix="auto", symprec=1e-5)` y `generate_displacements(distance=--distance 0.03 Å)`.
2. `kappa.configuraciones` convierte las superceldas desplazadas a ASE (fc3 y, si hay, fc2).
3. Fuerzas: (a) `--model mace|chgnet|m3gnet` → `kappa.fuerzas_mlip`; (b) sin `--model` → `kappa.escribir_inputs` escribe un `scf` por configuración en `fc3/dNNNN/pw.in` (y `fc2/`), `conv_thr = 1e-10`, malla por `--kspacing` 0.35 Å⁻¹, `occupations='fixed'` salvo `--metal` (smearing), más `correr.sh`; se niega por encima de 150 configuraciones sin `--force`; (c) `--collect` → `kappa.leer_fuerzas` lee `<forces>` de cada XML (Ha/bohr → eV/Å) y exige TODAS.
4. `kappa.resolver`: `produce_fc3`, `produce_fc2`, simetrización, `mesh_numbers = --mesh (13)`, `init_phph_interaction`, `run_thermal_conductivity(temperatures=--temps 100:800:8, is_isotope=--isotopes, boundary_mfp=--grain µm ×1e4 Å o 1e6)`.
5. `kappa.recoger` guarda κ (Voigt 6), Γ, velocidades, $C_\lambda$, pesos; `kappa.report`, `export` (`KAPPA.dat`, `KAPPA_recorrido.dat`, `KAPPA.txt`), `plot` (κ(T) log-log con guía $T^{-1}$; acumulada vs Λ).

**De dónde sale cada dato.**

| Dato | Origen | Detalle |
|---|---|---|
| Fuerzas | XML de pw.x (`output/forces`) o potencial MLIP | `qeout.read_xml` / `mlip.calculator` |
| fc2, fc3, Γ, v, C, κ | biblioteca `phono3py` | `Phono3py.thermal_conductivity` (RTA) |
| Malla de recorridos | `kappa.RECORRIDOS` | `np.logspace(0, 7, 141)` Å |
| Isótopos | phono3py (abundancias naturales) | `--isotopes` |

**Límites y trampas.** "Es RTA, no la solución exacta de la ecuación de Boltzmann. La RTA subestima κ (≈10-15 % en silicio, mucho más en grafeno o diamante)". "Solo hay dispersión de tres fonones". Sin `--isotopes`: "El silicio natural conduce ~10 % menos que el isotópicamente puro". Con potencial aprendido: "el valor absoluto puede estar lejos: con MACE-MP pequeño el silicio da ~51 W/mK a 300 K donde el experimento son ~140". Supercelda ≤ 8 celdas: "κ tiene que converger en el tamaño de la supercelda Y en la malla de q a la vez". Por omisión los scf de la fc2/fc3 llevan `occupations='fixed'` ("lo correcto para aislantes"); en un metal hay que pasar `--metal`, o los scf no convergen. No hay opción para escribir/leer los `fc2.hdf5`/`fc3.hdf5` de phono3py: cada `--collect` reconstruye todo.

**Referencias.** A. Togo, L. Chaput e I. Tanaka, *Phys. Rev. B* 91, 094306 (2015), DOI 10.1103/PhysRevB.91.094306 (phono3py). J. M. Ziman, *Electrons and Phonons*, Oxford (1960). L. Lindsay, D. A. Broido y T. L. Reinecke, *Phys. Rev. B* 87, 165201 (2013) (RTA vs. solución exacta).

---

### `olla-dft elph` — Acoplamiento electrón-fonón: λ, ω_log, Tc y τ

**Qué responde.** ¿Cuánto se acoplan los electrones a los fonones ($\lambda$), cuál es la temperatura crítica superconductora de Allen–Dynes con sus correcciones de acoplamiento fuerte, y cuánto vale el tiempo de relajación por fonones $\tau(T)$ que le falta al transporte en CRTA?

**Fundamento para no expertos.** Los electrones que se mueven por un metal chocan con las vibraciones de la red; cuánto chocan lo mide $\lambda$, un número adimensional. En un superconductor convencional ese mismo acoplamiento es lo que empareja a los electrones, y la fórmula de Allen–Dynes convierte $\lambda$ y una frecuencia típica de fonón ($\omega_{\log}$) en una temperatura crítica. El mismo $\lambda$ da el tiempo entre choques a alta temperatura, $\tau$. `ph.x` calcula el acoplamiento para varios ensanchamientos numéricos; el valor bueno es el del "plató", donde deja de depender del ensanchamiento.

**Fórmulas.** (`elph.lambda_de_a2F`, `omega_log_de_a2F`, `omega_2`, `factores_correccion`, `allen_dynes`, `tau_elph`)

$$
\lambda = 2\int \frac{\alpha^2F(\omega)}{\omega}\,d\omega, \qquad
\omega_{\log} = \exp\!\left[\frac{2}{\lambda}\int \ln\omega\,\frac{\alpha^2F(\omega)}{\omega}\,d\omega\right], \qquad
\bar\omega_2 = \left[\frac{2}{\lambda}\int \omega\,\alpha^2F(\omega)\,d\omega\right]^{1/2}
$$

$$
T_c = f_1 f_2\,\frac{\omega_{\log}}{1.2}\exp\!\left[\frac{-1.04(1+\lambda)}{\lambda - \mu^*(1+0.62\lambda)}\right], \qquad
f_1 = \left[1 + \left(\frac{\lambda}{\Lambda_1}\right)^{3/2}\right]^{1/3}, \quad \Lambda_1 = 2.46(1+3.8\mu^*)
$$

$$
f_2 = 1 + \frac{(r-1)\lambda^2}{\lambda^2 + \Lambda_2^2}, \quad r = \frac{\bar\omega_2}{\omega_{\log}}, \quad \Lambda_2 = 1.82(1+6.3\mu^*)\,r, \qquad
\frac{1}{\tau} = \frac{2\pi\lambda k_B T}{\hbar}
$$

- $\omega$ en THz en `a2F.dos*`; $\omega_{\log}$, $\bar\omega_2$ en K (`THZ_K` = 47.9924 K/THz); $\mu^*$ = 0.10, 0.13, 0.16 (rango, no se calcula); $T_c$ = 0 si el denominador ≤ 0; $\hbar$ = `HBAR_EVS` = 6.582119569e-16 eV·s, $k_B$ = 8.617333262e-5 eV/K; $\tau$ en s. Plató (`elph.plato`): tramo más largo de ≥ 3 λ consecutivos que no difieren más del 5 % del primero; su punto medio.

**Cómo lo calcula Olla-DFT.**
1. Preparación (`qekit/cli.py: _cmd_elph` sin `--collect` → `qekit/modules/elph.py: prepare`): escribe `1_scf.in` (malla `--kgrid` o `kspacing`), `2_nscf.in` con `la2F = .true.` y malla `--kgrid-nscf` (por omisión $q_i\cdot\max(2, \lceil 2k_i/q_i\rceil)$, múltiplo de la de q), y `3_ph.in` con `electron_phonon='interpolated'`, `el_ph_sigma = --sigma (0.005 Ry)`, `el_ph_nsigma = --nsigma (10)`, `fildvscf='dvscf'`, `tr2_ph = 1e-12`, `ldisp` con `--qgrid` (2x2x2). Smearing `methfessel-paxton`, `--degauss` 0.02 Ry, `conv_thr = 1e-10`.
2. El usuario corre `pw.x` ×2, `ph.x` y, opcionalmente, `lambda.x` (Olla-DFT tiene `elph.build_lambda_input` pero la CLI no lo escribe).
3. `--collect`: `elph.leer_elph_ph` lee de `ph.out` los ensanchamientos (`Gaussian Broadening: X Ry`) y `DOS = … states/spin/Ry`; `elph.leer_lambda_out` lee `lambda.dat` (columnas σ, λ, ∫α²F, ⟨log ω⟩, N(E_F)) o el texto de `lambda.out`, toma el $\mu^*$ de `lambda.in` (última línea numérica) y rellena la columna $T_c$ por ensanchamiento: de la tabla final `lambda omega_log T_c` de `lambda.out` si existe y cuadra en tamaño, y si no, la calcula fila a fila con `allen_dynes(λ_i, ω_log,i, μ*, correcciones=False)`; `ElPhRun.Tc_fuente` dice cuál de las dos cosas hizo y el reporte lo imprime. `elph.leer_a2F` lee `a2F.dos*` (o `A2F.dat`) y de ahí λ, $\omega_{\log}$ (si no venían) y, siempre, $\bar\omega_2$ (`elph.omega_2`, calculada en la CLI), que es la que activa el factor $f_2$ en las $T_c$ del resumen.
4. `elph.plato`, `elph.report` (tabla por ensanchamiento, régimen, $T_c$ para tres $\mu^*$, τ a 100/300/500/800 K), `elph.export` (`ELPH.dat`, `A2F.dat`, `ELPH.txt`), `elph.plot`.

**De dónde sale cada dato.**

| Dato | Origen | Detalle |
|---|---|---|
| Ensanchamientos, N(E_F) | `ph.out` (`Gaussian Broadening`, `DOS =`) | `elph.leer_elph_ph` |
| λ, ⟨log ω⟩ por ensanchamiento | `lambda.dat` / `lambda.out` de `lambda.x` | `elph.leer_lambda_out` |
| $\alpha^2F(\omega)$ | `a2F.dos*` (o `A2F.dat`) | `elph.leer_a2F`, columna 1 THz, columna 2 α²F |
| THz → K | `elph.THZ_K` | 47.9924 |
| $\mu^*$ | literal (0.10, 0.13, 0.16) | empírico |
| $T_{\mathrm{Debye}}$ | `--debye` | solo para marcar el régimen de validez de τ |

**Límites y trampas.** El valor de λ de `ph.out` NO se lee (solo σ y N(E_F)); λ sale de `lambda.x` o del $\alpha^2F$. La columna "Tc(K)" de la tabla por ensanchamiento es la de `lambda.x` (Allen–Dynes SIN $f_1 f_2$, con el $\mu^*$ de `lambda.in`) o la recalculada igual por Olla-DFT; las $T_c$ con correcciones y para los tres $\mu^*$ son las del bloque "Temperatura crítica", y sin `a2F.dos*` no hay $\bar\omega_2$ y $f_2 = 1$. Sin plató: "la malla de k es insuficiente … Cualquier lambda que se reporte de aquí es arbitrario". $\mu^*$ "es empírico (0.10-0.16) y NO se calcula aquí". τ "vale por ENCIMA de la temperatura de Debye; por debajo sobreestima la dispersión". `lambda.x` con malla q gruesa deja $\omega_{\log}$ en NaN: se recalcula del $\alpha^2F$ si existe. El τ NO se inyecta automáticamente en `transport`: el docstring del módulo lo dice explícitamente y da la secuencia (`transport --collect` → `elph --collect` → $\sigma(T) = [\sigma/\tau](T)\cdot\tau(T)$ a mano sobre las columnas de `TRANSPORTE.dat`).

**Referencias.** P. B. Allen y R. C. Dynes, *Phys. Rev. B* 12, 905 (1975), DOI 10.1103/PhysRevB.12.905. W. L. McMillan, *Phys. Rev.* 167, 331 (1968), DOI 10.1103/PhysRev.167.331. P. B. Allen, *Phys. Rev. B* 3, 305 (1971) (τ a alta T). G. Grimvall, *The Electron–Phonon Interaction in Metals*, North-Holland (1981).

---

### `olla-dft transport` — Transporte electrónico en CRTA: Seebeck, σ/τ, κ_e/τ, Lorenz y espín

**Qué responde.** A partir de las bandas sobre una malla densa: ¿cuánto valen el coeficiente Seebeck $S$, la conductividad por unidad de tiempo de relajación $\sigma/\tau$, la conductividad térmica electrónica $\kappa_e/\tau$, el factor de potencia $S^2\sigma/\tau$ y la concentración de portadores en función del potencial químico y la temperatura? ¿Se cumple Wiedemann–Franz? ¿Cómo se reparte el transporte entre los dos canales de espín?

**Fundamento para no expertos.** Un electrón en una banda se mueve con velocidad $v = (1/\hbar)\,dE/dk$. A una temperatura dada solo los estados a unos $k_BT$ del potencial químico participan en el transporte (la "ventana" $-\partial f/\partial E$). Sumando velocidad por velocidad sobre esa ventana sale la conductividad; ponderando además por $(E-\mu)$ sale el Seebeck, que mide cuánto voltaje aparece por grado de diferencia de temperatura. La aproximación de tiempo de relajación constante (CRTA) supone que todos los electrones chocan con la misma frecuencia $1/\tau$: entonces $\tau$ se cancela en $S$ y en el número de Lorenz (predicciones reales) pero no en σ ni en κ_e, que se reportan divididas por τ.

**Fórmulas.** (`transport._fd_derivative`, `transport.compute`, `lorenz`, `cancelacion`, `TransporteEspin`) Con $x = (E-\mu)/k_BT$, $-\partial f/\partial E = \mathrm{sech}^2(x/2)/(4k_BT)$, pesos $w_k = 1/N_k$, $V$ el volumen de la celda:

$$
\mathbf{v}_{n\mathbf{k}} = \frac{1}{\hbar}\nabla_{\mathbf{k}}E_{n\mathbf{k}}\ (\text{diferencias finitas periódicas, } \texttt{np.gradient}), \qquad
\mathbf{s}_m = \sum_{n\mathbf{k}} w_k\, \mathbf{v}\otimes\mathbf{v}\,(E-\mu)^m\left(-\frac{\partial f}{\partial E}\right)
$$

$$
\frac{\boldsymbol\sigma}{\tau} = \frac{e}{V}\mathbf{s}_0, \qquad
\mathbf{S} = -\frac{1}{T}\,\mathbf{s}_1\mathbf{s}_0^{-1}, \qquad
\frac{\boldsymbol\kappa_e}{\tau} = \frac{e}{VT}\mathbf{s}_2 - \mathbf{S}\mathbf{S}\,\frac{\boldsymbol\sigma}{\tau}\,T, \qquad
\mathrm{PF} = \bar S^2\,\bar\sigma/\tau
$$

$$
n = \frac{N_{\mathrm{elec}} - 2\sum_{n\mathbf{k}} w_k f(E_{n\mathbf{k}})}{V}, \qquad
L = \frac{\bar\kappa_e}{\bar\sigma T}, \qquad L_0 = 2.44\times10^{-8}\ \mathrm{W\,\Omega/K^2}, \qquad
c = \frac{|\bar\kappa_e|}{|\bar\kappa_e + \bar S^2\bar\sigma T|}
$$

$$
\sigma_{\mathrm{tot}} = \sigma_\uparrow + \sigma_\downarrow, \qquad
S_{\mathrm{tot}} = \frac{S_\uparrow\sigma_\uparrow + S_\downarrow\sigma_\downarrow}{\sigma_\uparrow+\sigma_\downarrow}, \qquad
P = \frac{\sigma_\uparrow-\sigma_\downarrow}{\sigma_\uparrow+\sigma_\downarrow}, \qquad S_{\mathrm{espín}} = S_\uparrow - S_\downarrow
$$

- $e$ = 1.602176634e-19 C; $\hbar$ = 6.582119569e-16 eV·s; $k_B$ = 8.617333262e-5 eV/K; σ/τ en S/(m·s); $S$ en V/K (µV/K en el reporte); κ_e/τ en W/(m·K·s); $n$ en cm⁻³ (positivo = huecos); barra = traza/3. La "cancelación" $c$ mide qué fracción sobrevive a la resta $\kappa^0 - S^2\sigma T$.

**Cómo lo calcula Olla-DFT.**
1. `qekit/cli.py: _cmd_transport` → `qekit/modules/transport.py: prepare`: primitiva estandarizada; `scf.in` (malla por `--kspacing` o la de configuración) y `nscf.in` con `K_POINTS automatic` `--grid` (16x16x16), `nosym=.true.` (malla completa), `nbnd = 2 × nbnd estimado`; `--metal` desactiva `occupations='fixed'`; `--nspin 2` y `--mag EL=valor` (que implica `nspin=2`) escriben scf y nscf con polarización de espín, requisito de `--spin-resolved`.
2. `--run` ejecuta scf y nscf; `--collect` → `transport.load` lee el primer `out/*.xml`, reconstruye la rejilla a partir de las fraccionarias (rechaza lo que no sea una rejilla completa) y deriva $E(\mathbf{k})$ con `np.gradient` sobre la malla envuelta periódicamente, pasando a cartesianas con $\mathbf{B}^{-T}$.
3. `transport.compute` sobre $T$ = `--temperatures` (300) y 201 valores de µ en $E_F \pm$ `--mu-span` (1 eV).
4. `transport.report` (mejor $S$ tipo p y n, PF máximo), `report_lorenz`, con `--spin-resolved` carga `spin=1` y `report_espin`; `transport.export` (`TRANSPORTE.dat`), `transport.plot`.

**De dónde sale cada dato.**

| Dato | Origen | Detalle |
|---|---|---|
| Autovalores, k, pesos | XML del nscf (`ks_energies`, Hartree → eV) | `qeout.read_xml`; `weights` se sustituyen por $1/N_k$ |
| $E_F$, $N_{\mathrm{elec}}$, volumen, celda | XML (`fermi_energy`, `nelec`, `cell`) | `qeout.read_xml` |
| Constantes | `transport.E_CHARGE`, `HBAR_EVS`, `KB_EV`, `L0_SOMMERFELD` | CODATA; $L_0$ = 2.44e-8 |

**Límites y trampas.** Solo CRTA: "Para dar σ en S/m hace falta un τ que venga de un ajuste a una medida o de un cálculo de electrón-fonón — Olla-DFT no lo inventa". NO calcula ZT (haría falta κ_red y τ) ni acopla automáticamente el τ de `elph`. No interpola bandas (a diferencia de BoltzTraP): con malla < 24 por lado o < 12000 puntos avisa "INSUFICIENTE … sigma sale en picos aislados". El número de Lorenz dentro del gap sufre cancelación catastrófica: "NO TE FÍES DE ESTE NÚMERO … sobrevive el X %"; solo se resumen puntos con $c > 0.10$. `--spin-resolved` sobre un XML con `nspin = 1` se rechaza con la instrucción de volver a preparar con `--nspin 2 --mag EL=0.7 --run`. Modelo de dos corrientes: "Vale mientras la dispersión con inversión de espín sea lenta … deja de valer cerca de [la temperatura de Curie]".

**Referencias.** G. K. H. Madsen y D. J. Singh, *Comput. Phys. Commun.* 175, 67 (2006), DOI 10.1016/j.cpc.2006.03.007 (BoltzTraP, misma formulación CRTA). N. W. Ashcroft y N. D. Mermin, *Solid State Physics*, cap. 13 (Wiedemann–Franz). N. F. Mott, *Proc. R. Soc. A* 153, 699 (1936) (modelo de dos corrientes).

---

### `olla-dft ballistic` — Conductancia de Landauer con `pwcond.x`

**Qué responde.** ¿Cuántos canales de conducción tiene un electrodo a cada energía (bandas complejas) y qué probabilidad de transmisión $T(E)$ tiene un nanocontacto o una molécula entre electrodos? ¿Cuál es su conductancia en unidades de $G_0$?

**Fundamento para no expertos.** En un cristal macroscópico el electrón choca muchas veces (transporte difusivo, `transport`). En un contacto de pocos átomos lo cruza de un tirón: no hay conductividad, hay conductancia, y la da la fórmula de Landauer: $G = G_0 T(E_F)$, con $T$ la probabilidad de pasar sumada sobre todos los "carriles" abiertos. Como $T$ no puede superar el número de carriles, la conductancia sale cuantizada en escalones de $G_0 = 2e^2/h$, y ver esos escalones es la señal de que el cálculo está bien.

**Fórmulas.** (`ballistic.G0`, `CondRun`)

$$
G = G_0\,T(E_F), \qquad G_0 = \frac{2e^2}{h} = 7.748091729\times10^{-5}\ \mathrm{S}, \qquad R = \frac{1}{G} = \frac{12.906\ \mathrm{k\Omega}}{T(E_F)}, \qquad T(E) \le N_{\mathrm{canales}}(E)
$$

- $T(E_F)$: transmisión en la energía más próxima a $E - E_F = 0$ de la ventana. Límites de región (`ballistic.longitud_z`): $\mathrm{bdl} = |\mathbf{a}_3|_{\mathrm{electrodo}}/a$, $\mathrm{bds} = |\mathbf{a}_3|_{\mathrm{dispersor}}/a$ con $a = |\mathbf{a}_1|$ (unidades de alat): la frontera de cada región es el final de SU celda, no la altura del último átomo.

**Cómo lo calcula Olla-DFT.**
1. `qekit/cli.py: _cmd_ballistic` → `qekit/modules/ballistic.py: prepare`: `comprobar_geometria` exige $\mathbf{a}_3 \parallel z$, $\mathbf{a}_{1,2} \perp z$ y, con `--scatterer`, la misma celda en el plano.
2. Escribe `scf_electrodo.in` (y `scf_dispersor.in`) con `insulator=False` y prefijos `electr`/`disper`, y `cond.in` (`&inputcond`: `ikind` = `--ikind` (solo 0 o 1) o 1 si hay dispersor / 0 si no; `ikind=1` sin `--scatterer` y `ikind=2` se rechazan con explicación; `energy0 = --emax` (3), `denergy = -(emax-emin)/(n-1)` (`--emin` −3, `--points` 61), `ewind = 1`, `epsproj = 1e-3`, `nz1 = --nz1` (3), `bdl` = `longitud_z(electrodo)`, `bds` = `longitud_z(dispersor)`, un punto k (0, 0, 1)).
3. El usuario corre `pw.x` y `pwcond.x`; `--collect` → `ballistic.collect` lee `trans*.dat` (E, T) o las líneas `T_tot` de `cond*.out`, y `Nchannels of the left tip` por energía (máximo sobre k); `ikind` del `.out` o de `cond.in`.
4. `ballistic._avisar`, `report`, `export` (`BALISTICO.dat`, `.txt`), `plot` (T(E) y escalones de canales).

**De dónde sale cada dato.**

| Dato | Origen | Detalle |
|---|---|---|
| $T(E)$ | `trans.dat` de `pwcond.x` (o `T_tot` en `cond.out`) | `ballistic.collect` |
| Canales abiertos | `cond.out` (`Nchannels of the left tip`) | máximo por energía |
| $G_0$ | `ballistic.G0` | 7.748091729e-5 S (CODATA) |
| Límites bdl/bds | longitud de la celda en $z$ (`longitud_z`) | en unidades de alat |

**Límites y trampas.** `ikind=0` "NO es la conductancia. Es el número de canales abiertos, que acota la conductancia por arriba". Si $T > N$: "Eso es imposible: T <= N por construcción. Revisa que los límites bdl/bds…". Transmisiones negativas: "el cálculo no convergió o … la geometría de las regiones está mal cortada". Electrodos distintos (`ikind=2` de `pwcond.x`) no están soportados: `--ikind` solo admite 0 y 1, y pedir 2 da "no está implementado … hay que escribir a mano el tercer scf y 'prefixr' y 'bdr' en cond.in". Un solo punto k transversal por omisión (0, 0, 1). No hay `--run`: siempre se corre a mano.

**Referencias.** R. Landauer, *IBM J. Res. Dev.* 1, 223 (1957); M. Büttiker, *Phys. Rev. Lett.* 57, 1761 (1986), DOI 10.1103/PhysRevLett.57.1761. H. J. Choi y J. Ihm, *Phys. Rev. B* 59, 2267 (1999), DOI 10.1103/PhysRevB.59.2267; A. Smogunov, A. Dal Corso y E. Tosatti, *Phys. Rev. B* 70, 045417 (2004) (`pwcond.x`).

---

### `olla-dft cost` — Estimador de coste calibrado con tu historial

**Qué responde.** ¿Cuánto va a tardar este barrido en ESTA máquina, con qué incertidumbre? (`cost` muestra el modelo; `--estimate` en cualquier barrido lo aplica.)

**Fundamento para no expertos.** El tiempo de un cálculo de ondas planas escala de forma conocida con el número de puntos k, de ondas planas, de bandas y de iteraciones. Lo que no se sabe de antemano es la constante de proporcionalidad de cada máquina. Olla-DFT toma la forma de la física y ajusta la escala con los cálculos que el usuario ya indexó en `olla-dft db` (con sus tiempos de pared), y mide cuánto se equivoca dejando fuera un sistema y prediciéndolo con los demás.

**Fórmulas.** (`cost.n_ondas_planas`, `trabajo`, `iteraciones`, `_ajusta`, `estimar`)

$$
N_{\mathrm{PW}} = \frac{V\,E_{\mathrm{cut}}^{3/2}}{6\pi^2}\ (\text{bohr}^3,\ \mathrm{Ry}), \qquad
w_1 = n_k\, s\, N_{\mathrm{PW}}\, n_{\mathrm{b}}, \qquad
w_2 = n_k\, s\, N_{\mathrm{PW}}\, n_{\mathrm{b}}^2
$$

$$
t = t_0 + \left(C_1 w_1 + C_2 w_2\right)\, n_{\mathrm{scf}}\, n_{\mathrm{ion}}, \qquad
\text{ajuste NNLS con pesos } t^{-1/2},\qquad
[t/\mathrm{disp},\ t\cdot\mathrm{disp}],\ \mathrm{disp} = e^{\sigma(\ln(\mathrm{pred}/\mathrm{real}))}
$$

- $V$: volumen (Å³ → bohr³ con `A3_BOHR3`); $n_k$: puntos k irreducibles (spglib `get_ir_reciprocal_mesh`, o el `number of k points` real de un `pw.out`); $s$: 1, 2 o 4 (`nspin`, no colineal); $n_{\mathrm{b}}$: `nbnd` del input o $\max(4, 2N_{\mathrm{at}})$; $n_{\mathrm{scf}}$: mediana del historial (14 por omisión); $n_{\mathrm{ion}}$: mediana por tipo (`relax` 8, `vc-relax` 12 por omisión). El ajuste con tres coeficientes solo si hay ≥ 8 cálculos y $w_1^{\max}/w_1^{\min} \ge 5$; si no, $C_1$ = mediana geométrica de $t/w_1$.

**Cómo lo calcula Olla-DFT.**
1. `qekit/cli.py: _cmd_cost` → `qekit/modules/cost.py: calibrar(--db olla-dft.db)`: `cost.historial` consulta la tabla `calculos` (natoms, ecutwfc, kgrid, nspin, volumen, n_scf, nk, nbnd, n_bfgs, wall_s, calculation).
2. `cost._prepara` construye $w_1 n_{it}$, $w_2 n_{it}$ y $t$; `cost._ajusta` (`scipy.optimize.nnls`, o `lstsq` recortado a ≥ 0).
3. Validación fuera de muestra por sistema (`_clave_sistema`: natoms, ecutwfc, nk, calculation, nspin) si hay ≥ 4 sistemas y ≥ 8 cálculos restantes: sesgo y dispersión de $\ln(\mathrm{pred}/\mathrm{real})$; si no, residuo del ajuste.
4. `cost.report_modelo` imprime $t_0$, $C_1$, $C_2$, iteraciones, precisión y avisos. En un barrido, `cli._run_or_explain` llama a `cost.estimar_barrido` (lee cada `pw.in` con `descriptores_de_input`; reutiliza el $n_k$ real de un punto ya corrido o del historial con la misma fórmula y malla) y `cost.report`.

**De dónde sale cada dato.**

| Dato | Origen | Detalle |
|---|---|---|
| Historial de tiempos | `olla-dft.db` (SQLite, tabla `calculos`, `wall_s`) | de `olla-dft db` |
| $N_{\mathrm{PW}}$ | volumen y `ecutwfc` del `pw.in` | `cost.n_ondas_planas`; verificado contra QE en los tests |
| $n_k$ irreducibles | spglib o `pw.out` (`number of k points`) | `cost.k_irreducibles`, `nk_de_salida` |
| Peso del ajuste | `cost.EXP_PESO` | 0.5 (elegido sobre 63 cálculos reales) |

**Límites y trampas.** "Esta herramienta distingue diez minutos de seis horas … No es un cronómetro". Sin historial: "No hay con qué calibrar: la base de cálculos está vacía o no guarda tiempos". Historial poco variado (`extrapola_bien` exige ≥ 8 cálculos y rango ≥ 5): "predecir un sistema de otro tamaño puede irse por un factor dos". spglib y `pw.x` no siempre ven la misma simetría: "ahí se va un factor dos o tres". No modela paralelismo MPI (el tiempo con `-j N` es simplemente total/N), ni el coste de `ph.x`, ni el de la memoria.

**Referencias.** M. C. Payne et al., *Rev. Mod. Phys.* 64, 1045 (1992), DOI 10.1103/RevModPhys.64.1045 (escalado de los métodos de ondas planas). C. L. Lawson y R. J. Hanson, *Solving Least Squares Problems*, SIAM (1995) (NNLS).

## Espectros, superficies, química y control de calidad

Esta parte documenta la física que hay detrás de los comandos de Olla-DFT que van más allá de la energía total: espectros ópticos y de rayos X (`optics`, `tddft`, `corehole`, `xanes`, `xps`), análisis de la densidad electrónica y de superficies (`charges`, `charge`, `wf`, `esm`, `surface`, `adsorb`, `interface`), química de defectos y de reacciones (`defect`, `eform`, `echem`, `neb`, `hull`), generación de estructuras con potenciales aprendidos (`amorphous`, `mlip`) y las herramientas de control de calidad que vigilan que todo lo anterior sea comparable y creíble (`audit`, `db`, `doctor`, `crosscheck`, `selftest`, `suggest`, `pseudos`). Cada sección dice qué pregunta contesta el comando, qué fórmulas implementa realmente el código (con la función responsable), de qué archivo de Quantum ESPRESSO sale cada dato y dónde están los límites. Cuando la documentación interna del código promete algo que el código no hace, se dice en "Límites y trampas".

---

### `olla-dft optics` — Función dieléctrica, absorción y gap de Tauc

**Qué responde.** ¿Cómo responde el material a la luz? Da $\varepsilon(\omega)$, el índice de refracción $n$, el coeficiente de extinción $k$, el coeficiente de absorción $\alpha$, la reflectividad $R$ y un gap óptico extrapolado como se hace con un espectro UV-Vis.

**Fundamento para no expertos.** Cuando la luz atraviesa un sólido, el campo eléctrico de la onda empuja a los electrones. Si la energía del fotón coincide con la que necesita un electrón para saltar de una banda ocupada a una vacía, la luz se absorbe. La *función dieléctrica* $\varepsilon(\omega) = \varepsilon_1 + i\varepsilon_2$ resume esa respuesta: la parte imaginaria $\varepsilon_2$ cuenta cuántos saltos hay a cada energía (absorción) y la parte real $\varepsilon_1$ cuánto se polariza el material (refracción). Las dos no son independientes: la causalidad (la respuesta no puede adelantarse a la causa) las liga por las relaciones de Kramers–Kronig, así que conociendo una se puede reconstruir la otra.

`epsilon.x` calcula $\varepsilon_2$ sumando todas las transiciones *verticales* entre bandas de Kohn–Sham, como si cada electrón saltara solo, sin sentir al hueco que deja. Es la aproximación de partícula independiente (RPA sin campos locales). El gap que sale es el del funcional, normalmente demasiado pequeño; el "scissor" (tijera) es la corrección más simple: se suben rígidamente todas las transiciones en $\Delta$ y se rehace $\varepsilon_1$ por Kramers–Kronig para no romper la causalidad.

**Fórmulas.** Todas en `qekit/modules/optics.py`.

Promedio isótropo (`optics.collect`):
$$\varepsilon_{1,2}(\omega) = \tfrac{1}{3}\left[\varepsilon_{xx} + \varepsilon_{yy} + \varepsilon_{zz}\right]$$

Funciones ópticas derivadas (`optics.derived`):
$$|\varepsilon| = \sqrt{\varepsilon_1^2 + \varepsilon_2^2},\qquad n = \sqrt{\frac{|\varepsilon| + \varepsilon_1}{2}},\qquad k = \sqrt{\frac{|\varepsilon| - \varepsilon_1}{2}}$$
$$\alpha(E) = \frac{2\,k\,E}{\hbar c},\qquad R = \frac{(n-1)^2 + k^2}{(n+1)^2 + k^2}$$

- $E = \hbar\omega$: energía del fotón (eV). $\hbar c$ = `HBAR_C_EV_CM` = $1.9732698\times10^{-5}$ eV·cm, de modo que $\alpha$ sale en cm⁻¹. $n$, $k$, $R$ son adimensionales. Los radicandos negativos se truncan a cero (`np.maximum`).

Kramers–Kronig (`optics.kramers_kronig`):
$$\varepsilon_1(\omega) = 1 + \frac{2}{\pi}\,\mathcal{P}\!\int_0^{\omega_{\max}} \frac{\omega'\,\varepsilon_2(\omega')}{\omega'^2 - \omega^2}\,d\omega'$$
- $\mathcal{P}$: valor principal; se implementa quitando el punto $\omega'=\omega$ de la cuadratura trapezoidal sobre la malla uniforme de `epsilon.x`. La integral se trunca en `wmax`.

Scissor (`optics.scissor`):
$$\varepsilon_2'(E) = \varepsilon_2(E-\Delta)\left(\frac{E-\Delta}{E}\right)^2,\qquad \varepsilon_1' = \mathrm{KK}[\varepsilon_2']$$
- $\Delta$: corrimiento en eV (`--scissor`). El factor $((E-\Delta)/E)^2$ viene de $\varepsilon_2 \propto |p|^2/\omega^2$ con elementos de matriz $|p|^2$ intactos. Se aplica a cada componente cartesiana y luego se promedia.

Gráfica de Tauc (`optics.tauc_gap`):
$$y(E) = \left(\alpha E\right)^{1/r},\qquad r = \tfrac{1}{2}\ (\text{directa permitida}),\quad r = 2\ (\text{indirecta})$$
$$E_g^{\mathrm{opt}} = -\frac{b}{m}\quad\text{con}\quad y \approx m E + b\ \text{ajustada sobre el primer frente de absorción}$$

**Cómo lo calcula Olla-DFT.**
1. `optics.prepare` resuelve pseudos y cutoffs (`sweep.prepare_common`, tarea `optics`) y **se niega** si alguno no es de norma conservada (`epsilon.x` no tiene elementos de matriz para USPP/PAW).
2. Escribe `scf.in` (malla del `kspacing` de configuración, 0.20 Å⁻¹ por omisión), `nscf.in` con malla densa (`--kspacing`, por omisión 0.12 Å⁻¹), `nosym=.true.` y `nbnd = 3 ×` la estimación de `inputgen._estimate_nbnd` (`nbnd_factor=3.0`), y `epsilon.in` (`calculation='eps'`, `smeartype='gauss'`, `intersmear=--smear` (0.10 eV), `wmin=0`, `wmax=--wmax` (20 eV), `nw=800`).
3. Con `--run`: `pw.x` scf → `pw.x` nscf (`runner.run_all`) → `optics.run_epsilon` lanza `epsilon.x` (buscado junto a `pw.x`).
4. `optics.collect` lee `epsr_<prefix>.dat` y `epsi_<prefix>.dat` (columnas: energía, xx, yy, zz) y promedia.
5. Si `--scissor Δ ≠ 0`: `optics.scissor` desplaza $\varepsilon_2$ y rehace $\varepsilon_1$ por `optics.kramers_kronig`.
6. `optics.derived` obtiene $n, k, \alpha, R$; `optics.tauc_gap` ajusta el gap directo e indirecto; `optics.report` imprime $\varepsilon_1(0)$ (valor en $E \approx 0.05$ eV), $n(0)$, el máximo de $\varepsilon_2$ y los gaps.
7. `optics.export` escribe `OPTICS.dat` con las columnas de `optics.OPTICS_COLUMNS` (`E(eV)`, `eps1`, `eps2`, `n`, `k`, `alpha(1/cm)`, `R`), nombradas en la última línea de comentario para que `optics.read_optics_dat` las lea por nombre; `optics.plot` dibuja la figura de tres paneles.

Detalle del ajuste de Tauc (`optics.tauc_gap`): la curva se suaviza con un promedio móvil de ~0.05 eV; el piso de ruido es el máximo de $y$ en el 1 % inicial del espectro; el frente arranca donde $y$ supera $\max(2\cdot\text{piso}, 10^{-3}\cdot\mathrm{mediana}(y>0))$ y $E > 0.1$ eV; termina en el primer máximo local que triplique el arranque o a lo sumo `max_span` = 1.5 eV más arriba; se ajusta una recta en una ventana de `fit_window` = 0.6 eV centrada en el punto de pendiente máxima de ese tramo. Devuelve `None` si no hay absorción, si la pendiente no es positiva o si el corte cae fuera del rango.

**De dónde sale cada dato.**

| Dato | Origen | Detalle |
|---|---|---|
| $\varepsilon_1(\omega)$ por dirección | `epsr_<prefix>.dat` de `epsilon.x` | `optics.collect`, columnas 1–3 |
| $\varepsilon_2(\omega)$ por dirección | `epsi_<prefix>.dat` de `epsilon.x` | `optics.collect` |
| $\hbar c$ | constante `optics.HBAR_C_EV_CM` | $1.9732698\times10^{-5}$ eV·cm |
| Tipo de cada pseudo (NC/US/PAW) | cabecera UPF (`pseudo_type`) | vía `sweep.prepare_common` |
| $\Delta$ (scissor) | parámetro `--scissor` | eV; recomendado gap exp./GW − gap DFT |
| Ensanchamiento | parámetro `--smear` | `intersmear` gaussiano, 0.10 eV |
| Ventana y puntos | `--wmax` (20 eV), `nw=800` | fijos en `optics.prepare` |
| Malla nscf | `--kspacing` (0.12 Å⁻¹) | `sweep.default_grid` |

**Límites y trampas.**
- Es RPA de partícula independiente: sin campos locales ni excitones. El reporte lo recuerda: *"Recuerda: RPA de partícula independiente y gap del funcional…"*.
- Sin pseudos NC el comando aborta: *"epsilon.x solo funciona con pseudopotenciales de NORMA CONSERVADA…"*.
- `epsilon.x` no incluye transiciones asistidas por fonones: en un semiconductor indirecto $\varepsilon_2 = 0$ por debajo del gap directo y el ajuste "indirect" **no** da el gap indirecto real (docstring de `tauc_gap`).
- La integral de Kramers–Kronig se trunca en `wmax`: $\varepsilon_1(0)$ hereda un error si hay absorción fuerte por encima de 20 eV.
- El scissor sólo mueve el gap; no corrige intensidades ni añade excitones.
- Si el ajuste de Tauc falla, el reporte imprime *"no se pudo ajustar"* en vez de un número.

**Referencias.**
- J. Tauc, R. Grigorovici, A. Vancu, *Phys. Status Solidi* 15, 627 (1966) — gráfica de Tauc.
- Manual de `epsilon.x` (Quantum ESPRESSO, paquete PP): A. Benassi, *"epsilon.x: a post-processing tool for the calculation of the dielectric properties"*.
- M. Dressel, G. Grüner, *Electrodynamics of Solids* (Cambridge, 2002) — Kramers–Kronig y funciones ópticas.

---

### `olla-dft tddft` — Absorción óptica con TDDFPT (Lanczos/Davidson)

**Qué responde.** ¿Cambia el espectro de absorción cuando el electrón excitado y el hueco que deja se ven entre sí? ¿Dónde están las primeras excitaciones, cuáles son brillantes y cuáles oscuras, y hay un excitón ligado por debajo del gap?

**Fundamento para no expertos.** `optics` suma transiciones de un electrón a la vez. En la realidad el electrón excitado (carga −) y el hueco (carga +) se atraen; en moléculas y en aislantes de gap ancho esa atracción baja la energía del par y crea un pico de absorción **dentro** del gap: el excitón. La teoría del funcional de la densidad dependiente del tiempo en respuesta lineal (TDDFPT) incluye esa interacción a través del kernel de intercambio-correlación. Quantum ESPRESSO la resuelve de dos formas: con el algoritmo de **Lanczos** (`turbo_lanczos.x` + `turbo_spectrum.x`), que da el espectro entero sin calcular estados vacíos, o con **Davidson** (`turbo_davidson.x`), que da una a una las primeras N excitaciones con su energía y su fuerza de oscilador $f$. Una excitación con $f \approx 0$ existe pero no absorbe luz: es "oscura".

**Fórmulas.** En `qekit/modules/tddft.py`.

Conversión de unidades de los inputs (`tddft.build_lanczos_input`, `build_spectrum_input`, `build_davidson_input`):
$$E_{\mathrm{Ry}} = \frac{E_{\mathrm{eV}}}{\mathrm{RY\_EV}},\qquad \mathrm{RY\_EV} = 13.605693122994\ \mathrm{eV}$$

Longitud de onda (`tddft.report`):
$$\lambda\,(\mathrm{nm}) = \frac{1239.84}{E\,(\mathrm{eV})}$$

Borde de absorción (`TddftRun.onset`): primer máximo local de $dS/dE$ que supere el 20 % del máximo de la derivada (punto de inflexión de la primera subida).

Firma de excitón (`tddft._avisar`):
$$d = E_{\mathrm{onset}} - E_g^{\mathrm{IP}},\qquad d < -\max(0.10\ \mathrm{eV},\ 2\,\eta)\ \Rightarrow\ \text{excitón ligado}$$
- $E_g^{\mathrm{IP}}$: gap de partículas independientes dado por el usuario (`--gap`); $\eta$: ensanchamiento en eV, el de `--broadening` en `--collect` o, si se omite, el leído de `spectrum.in` (`epsil`) o `davidson.in` (`broadening`) por `tddft._broadening_de_inputs` (Ry → eV). `UMBRAL_EXCITON` = 0.10 eV; `BROADENING_DEFAULT` = 0.05 eV.

Anisotropía (`tddft._anisotropia`): $\max_E[\max_i S_i(E) - \min_i S_i(E)] / \max_{i,E} S_i(E)$ sobre las componentes $x,y,z$.

**Cómo lo calcula Olla-DFT.**
1. `tddft.prepare`: si el vacío mínimo (`_vacio_minimo`) supera 5 Å o se pasa `--gamma`, usa `K_POINTS gamma` (lo único que TDDFPT implementa); si no, una malla automática con aviso de que `turbo_*.x` se plantará. Escribe `scf.in` con `nosym` y `noinv`.
2. Lanczos: `lanczos.in` (`itermax=--iter` 500, `ipol=--pol` 4 → `n_ipol=3`, `ltammd` con `--tamm-dancoff`, `lrpa` con `--rpa`, `scissor=--scissor/RY_EV` si se pide un corrimiento rígido de las bandas vacías; `prepare` rechaza un scissor negativo o con `--method davidson`) y `spectrum.in` (`itermax0=itermax`, `itermax=4×itermax` para la extrapolación `--extrapolation` osc/constant/no, `epsil=--broadening/RY_EV`, `units=1` (eV), `start/end/increment`).
3. Davidson: `davidson.in` con `num_eign=--states` (10), `num_init=2N`, `num_basis_max=max(80, 8N)`, `residue_conv_thr=1e-4`, `p_nbnd_virt=15`, ventana y `broadening` en Ry, `reference` en el centro de la ventana.
4. El usuario corre `pw.x` → `turbo_lanczos.x` → `turbo_spectrum.x` (o `turbo_davidson.x`) a mano.
5. `tddft.collect --collect` (con `broadening` de `--broadening` o de los inputs): Lanczos lee el primer `*plot*.dat` (columnas: E en eV, S total, S_x, S_y, S_z) y de `lanczos.out` el `itermax` y el funcional. Davidson (`_collect_davidson`) lee `<prefix>.eigen` (energía en Ry → eV, fuerza total, fuerzas por dirección) y el `*plot*.dat` si existe.
6. `_picos` lista máximos locales por encima del 5 % del máximo; `_avisar` compara el borde con `--gap`; `report` marca como brillantes las excitaciones con $f > 0.01$ y cuenta las oscuras.
7. `export` escribe `TDDFT.dat`, `TDDFT_EXCITACIONES.dat`, `TDDFT.txt`; `plot` superpone opcionalmente el espectro de `optics` (`--compare OPTICS.dat`: la CLI lee la columna `alpha(1/cm)` por su nombre con `optics.read_optics_dat` y la normaliza al máximo del TDDFPT).

**De dónde sale cada dato.**

| Dato | Origen | Detalle |
|---|---|---|
| $S(E)$ y componentes | `<prefix>.plot_S.dat` (o `*plot*.dat`) de `turbo_spectrum.x` | `tddft.collect`; energía en eV por `units=1` |
| Excitaciones $(E, f, f_x, f_y, f_z)$ | `<prefix>.eigen` de `turbo_davidson.x` | `_collect_davidson`; E en Ry × `RY_EV` |
| `itermax`, funcional | `lanczos.out` | regex en `tddft.collect` |
| Gap IP | parámetro `--gap` | eV |
| Ensanchamiento $\eta$ | `--broadening` o `spectrum.in`/`davidson.in` | `tddft._broadening_de_inputs`, Ry × `RY_EV` |
| Scissor | `--scissor` (sólo Lanczos) | eV → Ry en `lanczos.in` |
| $\alpha(E)$ de `optics` | `OPTICS.dat`, columna `alpha(1/cm)` | `optics.read_optics_dat` |
| Ry → eV | `tddft.RY_EV` | 13.605693122994 |
| Vacío mínimo | geometría de la celda | `tddft._vacio_minimo` |

**Límites y trampas.**
- Con LDA/GGA el kernel adiabático **no** liga excitones en un sólido; el reporte lo dice: *"con LDA o GGA el kernel adiabático NO liga excitones en un SÓLIDO… En MOLÉCULAS sí mejora."*
- Sólo punto Γ: con malla k el reporte avisa *"OJO: TDDFPT solo tiene implementado el caso gamma y se plantará al leer el input"*.
- Molécula con < 6 Å de vacío y < 30 átomos: *"AVISO: solo hay X Å de vacío…"*.
- `--scissor` sólo existe en `turbo_lanczos.x`: con `--method davidson` o con valor negativo el comando aborta (*"--scissor solo existe en turbo_lanczos.x…"*). El scissor de TDDFPT no rehace nada por Kramers–Kronig: lo aplica el propio código de QE a las bandas vacías.
- `--compare` exige un `OPTICS.dat` con la columna `alpha(1/cm)`; si falta: *"'…' no tiene la columna 'alpha(1/cm)'; --compare espera el OPTICS.dat de 'olla-dft optics'."*
- Si ni `--broadening` ni los inputs dan el ensanchamiento, el umbral de excitón se queda en `UMBRAL_EXCITON` = 0.10 eV.
- El comando no lanza los ejecutables `turbo_*.x`: sólo escribe inputs y lee salidas.

**Referencias.**
- D. Rocca, R. Gebauer, Y. Saad, S. Baroni, *J. Chem. Phys.* 128, 154105 (2008) — TDDFPT Lanczos.
- O. B. Malcıoğlu, R. Gebauer, D. Rocca, S. Baroni, *Comput. Phys. Commun.* 182, 1744 (2011) — turboTDDFT.
- X. Ge, S. J. Binnie, D. Rocca, R. Gebauer, S. Baroni, *Comput. Phys. Commun.* 185, 2080 (2014) — turboTDDFT 2.0 (Davidson).

---

### `olla-dft corehole` — Pseudopotenciales con hueco de core (ld1.x)

**Qué responde.** ¿Cómo describir un átomo al que se le ha arrancado un electrón de una capa interna? Genera el par de pseudopotenciales (normal + con hueco de core) que necesitan `xps` y `xanes`, con la misma configuración y los mismos radios, y extrae la función de onda de core que lee `xspectra.x`.

**Fundamento para no expertos.** Un pseudopotencial reemplaza al núcleo y a los electrones internos ("core") por un potencial efectivo, para que el cálculo sólo trate los electrones de valencia. Para simular una espectroscopia de rayos X hay que quitar un electrón de ese core congelado: eso exige un pseudopotencial distinto, generado a propósito con el programa atómico `ld1.x`, en el que la ocupación del nivel de core (1s para el borde K, 2p para el L₂,₃, etc.) vale uno menos. Como el core tiene un electrón menos, la carga de valencia declarada `z_valence` sube exactamente en 1: esa unidad **es** el hueco. Los dos pseudos deben generarse juntos con los mismos parámetros, porque comparar energías hechas con pseudos de familias distintas no significa nada.

**Fórmulas.** Este módulo no evalúa fórmulas físicas; construye los inputs de `ld1.x` a partir de reglas explícitas en `qekit/core/atomconf.py` y `qekit/modules/corehole.py`:

- Configuración electrónica por llenado de Aufbau (`atomconf.aufbau`, orden de Madelung `ORDEN`) con las excepciones de `atomconf.EXCEPCIONES` (Cr, Cu, Nb, Mo, Ru, Rh, Pd, Ag, La, Ce, Gd, Pt, Au).
- Partición core/valencia (`atomconf.particion`): valencia = capa $s,p$ de $n_{\max}$ + cualquier $d$/$f$ parcialmente llena + $d$ llena de la fila anterior; con `--semicore`, además $(n-1)s,(n-1)p$.
- Hueco (`atomconf.config_hueco`): ocupación del nivel `BORDES[borde]` reducida en 1.0; se rechaza si el nivel no está en el core.
- Canales de pseudización (`atomconf.canales_pseudo`): la valencia más un canal **desocupado** (ocupación −2) por cada $l \le 2$ ausente, con $n = \max(n_{\max}, l+1)$; con `--projectors 2` un segundo proyector por canal con etiqueta $n+1$ y ocupación −1.
- Radio de corte por fila (`corehole.RCUT_FILA`): {1: 1.0, 2: 1.3, 3: 1.7, 4: 2.0, 5: 2.2, 6: 2.4} bohr; `rcutus = 1.25 · rcut` sólo si `pseudotype=3`.
- Energías de referencia de canales no ligados: `E_CANAL_VACIO` = 0.15 Ry; segundo proyector `E_SEGUNDO_PROYECTOR` = 0.05 Ry.

**Cómo lo calcula Olla-DFT.**
1. `corehole.generar`: valida el elemento (H..Rn), fuerza `pseudotype=3` si se piden 2 proyectores, obtiene partición, canales y `rcut` (o `--rcut`).
2. `corehole.input_ld1` escribe `ld1_base.in` y `ld1_hueco.in` (`iswitch=3`, `rel=--rel`, `beta=0.3`, `dft=--functional` (PBE), `tm=.true.`, `lloc` = mayor $l$ de los canales, `lgipaw_reconstruction=.true.`, `author='Olla-DFT'`). Los canales vacíos se añaden también a la configuración de todos los electrones (`_con_canales_vacios`) porque `ld1.x` exige que existan.
3. `corehole._correr_ld1` ejecuta `ld1.x < ld1_X.in > ld1_X.out` (salvo `--only-inputs`) y falla si aparece `Error in routine`.
4. `corehole.leer_upf` lee de cada UPF `element`, `z_valence`, `mesh_size`, `pseudo_type`, `functional`, `wfc_cutoff`, `rho_cutoff` y las etiquetas `PP_GIPAW_CORE_ORBITAL`.
5. `corehole.verificar` aplica las comprobaciones de la tabla siguiente; `report` y `export` (`PSEUDOS_HUECO.txt`) las listan. El código de salida es 1 si hay alguna `FALLA`.
6. Con `--core-wfc UPF`: `corehole.core_wfc` extrae las funciones de onda de core en el formato de `filecore` de `xspectra.x` (un bloque por orbital, separados por línea en blanco, orden del UPF) y verifica que el número de puntos coincida con `mesh_size`.

| Comprobación (`corehole.verificar`) | Criterio | Marca |
|---|---|---|
| Diferencia de `z_valence` | exactamente +1 (tolerancia 1e-6) | FALLA si no |
| Mallas radiales | `mesh_size` igual en los dos UPF | FALLA si no |
| Orbital del hueco | presente entre los `PP_GIPAW_CORE_ORBITAL` del UPF con hueco | FALLA si no |
| Funcional | mismo en los dos UPF | FALLA si no |
| Proyectores | aviso si sólo hay uno por canal (XSpectra recomienda dos) | aviso |
| Estados fantasma, derivadas logarítmicas, transferibilidad | **no se comprueban** | aviso explícito |

**De dónde sale cada dato.**

| Dato | Origen | Detalle |
|---|---|---|
| Configuración electrónica | `atomconf.aufbau` + `EXCEPCIONES` | regla, no dato experimental |
| Nivel del borde | `atomconf.BORDES` | K=1s, L1=2s, L23=2p, M1=3s, M23=3p, M45=3d |
| `z_valence`, `mesh_size`, tipo, funcional | cabecera `PP_HEADER` del UPF generado | `corehole.leer_upf` |
| Orbitales de core | secciones `PP_GIPAW_CORE_ORBITAL.n` del UPF | `leer_upf`, `core_wfc` |
| Malla radial | `PP_R` del UPF | `core_wfc` |
| Radio de corte | `RCUT_FILA` o `--rcut` | bohr |

**Límites y trampas.**
- Los bordes M (`M1`, `M23`, `M45`) existen en `atomconf.BORDES` y sirven para generar el pseudo con hueco (XPS), pero `xspectra.x` sólo implementa K, L1, L2, L3 y L23: `olla-dft xanes` los rechaza (`xanes.validar_borde`).
- El reporte avisa: *"NO verificado automáticamente: estados fantasma, derivadas logarítmicas y transferibilidad… el cutoff del pseudo anterior NO sirve para este."* Hay que volver a converger con `olla-dft converge`.
- Con `--projectors 2` el pseudo sale ultrasuave y *"casi siempre hay que ajustar --rcut a mano hasta que ld1.x converja"*.
- `ld1.x` no se compila por omisión en QE (`make ld1`).

**Referencias.**
- A. Dal Corso, *Comput. Mater. Sci.* 95, 337 (2014) — pslibrary y `ld1.x`.
- N. Troullier, J. L. Martins, *Phys. Rev. B* 43, 1993 (1991) — pseudización TM (`tm=.true.`).
- C. J. Pickard, F. Mauri, *Phys. Rev. B* 63, 245101 (2001) — reconstrucción GIPAW.

---

### `olla-dft xanes` — Absorción de rayos X cerca del borde (xspectra.x)

**Qué responde.** ¿Qué forma tiene el espectro XANES/NEXAFS de un átomo concreto en un borde concreto, con hueco de core y polarización dada, y cuánto depende de la dirección del campo?

**Fundamento para no expertos.** Un fotón de rayos X arranca un electrón de un nivel profundo (1s en el borde K) y lo manda a los estados vacíos. La regla de selección dipolar sólo permite estados finales con momento angular $l \pm 1$: desde 1s se ven los estados $p$ vacíos **de ese átomo**. El espectro es, en esencia, la densidad de estados vacíos proyectada sobre el absorbedor, y por eso es local, selectivo por elemento y sensible al estado de oxidación y a la coordinación. El hueco que deja el electrón atrae a los estados vacíos y corre el borde, así que el átomo que absorbe se describe con el pseudopotencial de hueco de core de `corehole`, y como el electrón arrancado se supone fuera del sistema, la celda lleva carga total +1 (aproximación de hueco de core completo, FCH). `xspectra.x` calcula la sección eficaz con el método de Lanczos y fracciones continuas sin construir los estados vacíos.

**Fórmulas.** En `qekit/modules/xanes.py`.

Promedio de polvo (`xanes.collect`):
$$\sigma(E) = \tfrac{1}{3}\left[\sigma_x(E) + \sigma_y(E) + \sigma_z(E)\right]$$

Distancia mínima entre imágenes del absorbedor (`xanes.distancia_imagen_minima`):
$$d_{\min} = \min_{(i,j,k)\neq 0,\ |i|,|j|,|k|\le 1}\left|i\,\mathbf{a} + j\,\mathbf{b} + k\,\mathbf{c}\right|$$

Borde operativo (`xanes.onset`): primera energía en la que $\sigma \ge 0.5\,\sigma_{\max}$. Anisotropía (`_anisotropia`): $\max_E[\mathrm{ptp}_i\,\sigma_i(E)]/\max\sigma$; se destaca si $> 0.1$.

**Cómo lo calcula Olla-DFT.**
1. `xanes.validar_borde` (también en `_cmd_xanes`) normaliza `--edge` y sólo admite `BORDES_XSPECTRA` = K, L1, L2, L3, L23; los bordes M se rechazan con un mensaje explícito. `BORDE_COREHOLE` dice con qué `--edge` de `corehole` se genera el hueco de cada borde (L2 y L3 comparten el hueco 2p = `L23`). `xanes.prepare` localiza el átomo `--element`/`--site`, lo mueve al **primer** lugar de la lista y lo declara como especie aparte con etiqueta de tres letras (`etiqueta_excitada`, p. ej. `Sih`; el límite de QE es `CHARACTER(LEN=3)`).
2. `sweep.prepare_common` (tarea `xanes`, excluyendo el UPF de hueco) y `inputgen.build_pw_input` escriben `scf.in` con `tot_charge = 1.0`; `_marcar_absorbedor` añade la especie excitada a `ATOMIC_SPECIES`, cambia la etiqueta del primer átomo y sube `ntyp` en 1.
3. `corehole.core_wfc` extrae `<El>.wfc` del UPF de hueco (secciones `PP_GIPAW_CORE_ORBITAL`).
4. `xanes.build_xspectra_input` escribe `xspectra_pol.in` (o `xspectra_x/y/z.in` con `--average`): `calculation='xanes_dipole'`, `edge=--edge`, `xiabs=1`, `xepsilon=--polarization`, `xniter=2000`, `xcheck_conv=10`, `xerror=0.001`, `&plot`: `xnepoint=1000`, `xgamma=--broadening` (0.8 eV), `xemin=-10`, `xemax=30`, `terminator`, `cut_occ_states=.true.`; `&pseudos`: `filecore`, `r_paw(1)=--r-paw` (3.0); `&cut_occ`: `cut_desmooth=0.1`, `cut_stepl=0.01`; malla k al final.
5. El reporte mide $d_{\min}$ y avisa si es menor que `DIST_MINIMA` = 8 Å.
6. El usuario corre `pw.x -in scf.in` y `xspectra.x -in xspectra_*.in`.
7. `xanes.collect --collect` lee todos los `xanes_*.dat` (columnas E − E_F, σ), promedia si hay varios, y lee el `xgamma` del comentario *"Broadening parameter (in eV)"* del primer archivo.
8. `report` da borde al 50 %, máximo principal, picos (> 5 % del máximo), anisotropía; `export` escribe `XANES.dat` y `XANES.txt`; `plot` la figura.

**De dónde sale cada dato.**

| Dato | Origen | Detalle |
|---|---|---|
| $\sigma(E)$ por polarización | `xanes_<dir>.dat` de `xspectra.x` | `xanes._leer_dat` |
| Ensanchamiento `xgamma` | cabecera del `xanes_*.dat` | regex *"Broadening parameter"* |
| Función de onda de core | UPF de hueco (`PP_GIPAW_CORE_ORBITAL`) | `corehole.core_wfc` |
| Carga total +1 | fija en `xanes.prepare` | `tot_charge=1.0` |
| Polarización | `--polarization` (1 0 0) o ejes `EJES` con `--average` | vector cartesiano (`xcoordcrys=.false.`) |
| Malla k | `--kspacing` → `sweep.default_grid` | también en `xspectra.in` |
| $d_{\min}$ | vectores de la celda | `distancia_imagen_minima` |

**Límites y trampas.**
- El eje de energía es relativo al nivel de Fermi, no energía de fotón: *"Para comparar con un experimento se alinea el borde y se compara la FORMA."*
- Aviso de supercelda: *"AVISO: X Å es poco. Con condiciones periódicas el hueco de core ve sus propias imágenes…"* (umbral 8 Å).
- Una sola polarización: *"UNA sola polarización. En un cristal anisótropo el espectro depende de la dirección…"*.
- El borde (`xanes.onset`) es el primer punto donde σ alcanza el 50 % del **máximo global**: un pre-borde débil antes de la línea blanca no cuenta como borde (así lo declara ahora el docstring).
- Bordes M: *"xspectra.x solo calcula bordes K y L (K, L1, L2, L3, L23); los bordes M no están implementados en QE, aunque 'olla-dft corehole' pueda generar el pseudo con ese hueco."*
- `distancia_imagen_minima` sólo mira las 26 celdas vecinas: para celdas muy oblicuas puede sobrestimar $d_{\min}$.
- Sin `--core-hole` el comando aborta: *"falta --core-hole con el UPF de hueco de core. Sin él se calcularía el espectro del estado fundamental…"* y sugiere `olla-dft corehole <El> --edge <BORDE_COREHOLE[borde]>`.

**Referencias.**
- M. Taillefumier, D. Cabaret, A.-M. Flank, F. Mauri, *Phys. Rev. B* 66, 195107 (2002) — XSpectra, Lanczos con fracciones continuas.
- C. Gougoussis, M. Calandra, A. P. Seitsonen, F. Mauri, *Phys. Rev. B* 80, 075102 (2009) — XSpectra con PAW/GIPAW.
- O. Bunău, M. Calandra, *Phys. Rev. B* 87, 205105 (2013) — bordes L₂,₃.

---

### `olla-dft xps` — Corrimientos de nivel de core en estado inicial (initial_state.x)

**Qué responde.** ¿Cuánto se desplaza la energía del nivel de core de cada átomo respecto de los demás de su misma especie? Es la contraparte teórica del corrimiento químico de un espectro XPS.

**Fundamento para no expertos.** En XPS se mide la energía necesaria para arrancar un electrón de core. Un átomo rodeado de vecinos electronegativos tiene su core más ligado (corrimiento positivo) que uno en un entorno metálico. La aproximación de **estado inicial** calcula sólo cómo cambia el potencial que siente el electrón de core *antes* de arrancarlo; ignora la relajación de los demás electrones alrededor del hueco (el *estado final*), que puede valer varias décimas de eV. Por eso lo que sale son corrimientos **relativos** entre sitios, no energías de enlace absolutas. `initial_state.x` necesita dos especies del mismo elemento en el input —la normal y una con hueco de core— porque define el corrimiento a partir de `delta_zv = zv(excitada) − zv(normal)`; si las dos son iguales devuelve ceros sin avisar.

**Fórmulas.** En `qekit/modules/xps.py`. El corrimiento lo calcula `initial_state.x`; Olla-DFT sólo lo lee y lo reordena:

$$\Delta_i = \text{shift}_i^{\mathrm{TOTAL}},\qquad \Delta_i^{\mathrm{rel}} = \Delta_i - \min_j \Delta_j,\qquad \text{dispersión} = \max_i\Delta_i - \min_i\Delta_i$$

Indicador de cancelación (`xps.report`):
$$\frac{\max_{c}\,\mathrm{ptp}(\text{contribución}_c)}{\text{dispersión}} > 20 \Rightarrow \text{aviso de cancelación numérica}$$

- $\Delta_i$: corrimiento del átomo $i$ en eV, leído de la línea `atom i type t shift = … Ry, = … eV` de la sección *TOTAL*. Las contribuciones $c$ (Fermi, local, no local, iónica, core-correction, Hubbard…) se leen de las secciones *"The X contribution to shift"*.

**Cómo lo calcula Olla-DFT.**
1. `xps.prepare` lee `--core-hole EL=archivo.UPF` (repetible). Para cada elemento: `_verificar_par` exige que el UPF normal y el de hueco sean archivos distintos y que `z_valence` difiera exactamente en +1 (`qekit.core.pseudo.z_valence`).
2. `inputgen.build_pw_input` escribe `scf.in` con las especies extra (`extra_species`) declaradas en `ATOMIC_SPECIES` **sin** ningún átomo que las use; `_copiar_pseudos` copia el UPF de hueco a `pseudo_dir`.
3. `xps.build_input` escribe `initial_state.in` con `excite(t_normal) = t_hueco` (índices base 1 en el orden de `ATOMIC_SPECIES`); se rechaza `excite(t)=t`.
4. `structure.symmetry_dataset` cuenta órbitas de átomos equivalentes; si hay una sola avisa que todo saldrá cero.
5. El usuario corre `pw.x -in scf.in` y `initial_state.x -in initial_state.in > initial_state.out`.
6. `xps.collect --collect` parsea `initial_state.out` con `_RE_SECCION` y `_RE_ATOMO`, toma la columna en eV, y marca `equivalentes=True` si todos los $|\Delta_i| < 10^{-6}$ eV.
7. `report` tabula corrimientos, relativo al mínimo, dispersión, descomposición por contribución y el aviso de cancelación; `export` escribe `XPS_CORE.dat`.

**De dónde sale cada dato.**

| Dato | Origen | Detalle |
|---|---|---|
| Corrimiento por átomo y contribución | `initial_state.out` | `xps.collect`, regex `atom N type T shift = X Ry, = Y eV` |
| `z_valence` normal y con hueco | cabeceras UPF | `pseudo.z_valence` en `_verificar_par` |
| Sitios inequivalentes | spglib vía `structure.symmetry_dataset` | `equivalent_atoms` |
| Etiqueta de la especie excitada | `xanes.etiqueta_excitada` | 3 caracteres |
| Símbolos por átomo | estructura de entrada | `atoms.get_chemical_symbols()` |

**Límites y trampas.**
- **Sólo estado inicial.** No hay ΔSCF ni estado final; el docstring del módulo lo declara ahora explícitamente: el UPF con hueco de core se usa *únicamente como la "especie excitada" que initial_state.x necesita para definir el corrimiento, no para relajar el sistema frente al hueco*. El reporte remite: *"las energías de enlace absolutas necesitan un ΔSCF con hueco de core."*
- Dispersión < 0.1 eV: *"Por debajo de ~0.1 eV el corrimiento no es concluyente: la relajación de estado final… es del mismo orden."*
- Sin `--core-hole` se escribe sólo `scf.in` y el reporte explica que `initial_state.x` devolvería ceros.
- Todos los átomos equivalentes: *"AVISO: todos los átomos son equivalentes por simetría, así que todos los corrimientos van a salir exactamente cero."*
- Cancelación grande: *"CUIDADO con la cancelacion… baja conv_thr (1e-10 o menos) y sube la malla k antes de creerte la tercera cifra."*
- Los mensajes de error remiten a `olla-dft corehole <El> --edge K` para generar el par consistente.

**Referencias.**
- E. Pehlke, M. Scheffler, *Phys. Rev. Lett.* 71, 2338 (1993) — estado inicial vs. estado final en corrimientos de core.
- L. Köhler, G. Kresse, *Phys. Rev. B* 70, 165405 (2004) — energías de enlace de core con hueco.
- Documentación de `initial_state.x` (Quantum ESPRESSO, paquete PP).

---

### `olla-dft charges` — Cargas de Löwdin, Bader on-grid y diferencia de densidad

**Qué responde.** ¿Cuánta carga electrónica "pertenece" a cada átomo, y dónde se acumula o se vacía la densidad al formarse un enlace o una adsorción?

**Fundamento para no expertos.** La densidad electrónica es continua; repartirla entre átomos exige una regla. **Löwdin** proyecta los estados sobre orbitales atómicos ortogonalizados (lo hace `projwfc.x`); es barata y depende de la base de orbitales del pseudopotencial. **Bader** no usa orbitales: divide el espacio en "cuencas" siguiendo la subida más empinada de la densidad desde cada punto hasta un máximo, como el agua de lluvia que baja por laderas hasta cada valle, pero al revés. La **diferencia de densidad** $\rho_{AB} - \rho_A - \rho_B$ muestra, punto a punto, qué cambió al juntar dos partes.

**Fórmulas.** En `qekit/modules/charges.py`.

Löwdin (`charges.read_lowdin`, `report_lowdin`):
$$q_i^{\mathrm{neta}} = Z_i^{\mathrm{val}} - Q_i^{\mathrm{Löwdin}}$$

Bader on-grid (`charges.bader`): para cada punto de la rejilla se elige el vecino $\nu$ (de 26) que maximiza la pendiente
$$s_\nu = \frac{\rho(\mathbf{r}+\mathbf{d}_\nu) - \rho(\mathbf{r})}{|\mathbf{d}_\nu|}$$
y se sigue la cadena hasta un máximo local (compresión de caminos, máx. 64 iteraciones). Cada máximo se asigna al átomo más cercano con imágenes periódicas. Entonces
$$Q_i = \sum_{\mathbf{r}\in\Omega_i}\rho(\mathbf{r})\,\Delta V,\qquad V_i = N_i\,\Delta V_{\mathrm{Å}^3},\qquad \Delta V_{\mathrm{Å}^3} = \frac{V_{\mathrm{celda}}}{n_1 n_2 n_3},\quad \Delta V = \frac{\Delta V_{\mathrm{Å}^3}}{a_0^3}$$
- $\rho$: densidad del `.cube` en e/bohr³ (lo que escribe `pp.x`, `density_units="e/bohr3"` por omisión; `"e/A3"` también se admite); `charges._voxel_volume` devuelve el volumen del vóxel en las unidades de la densidad (bohr³, con $a_0$ = `fields.BOHR` = 0.529177210903 Å) para que $\rho\,\Delta V$ sean electrones, y en Å³ para reportar los volúmenes de cuenca.

Diferencia de densidad (`charges.difference`, `report_difference`), con el mismo $\Delta V$ en bohr³:
$$\Delta\rho = \rho_{\mathrm{total}} - \sum_p \rho_p,\qquad Q_{\mathrm{neta}} = \sum \Delta\rho\,\Delta V,\qquad Q_{\mathrm{acum}} = \sum_{\Delta\rho>0}\Delta\rho\,\Delta V$$

**Cómo lo calcula Olla-DFT.**
1. Si se da la estructura, `charges.valence_from_pseudos` lee `z_valence` de los UPF de `--pseudo-dir` (o del `pseudo_dir` de la configuración) vía `pseudo.resolve`; si algún UPF no se puede leer, devuelve `None`, la CLI avisa (*"no pude leer z_valence de los UPF…"*) y la columna "neta" queda en `n/d`.
2. `--lowdin projwfc.out`: `charges.read_lowdin` busca `Atom #  i: total charge = q` y `Spilling Parameter:`; con la estructura pone símbolos y la carga neta $Z^{\mathrm{val}} - Q$.
3. `--bader densidad.cube` (necesita la estructura): `fields.read_cube` lee el cube (`plot_num=0` de `pp.x`), `charges.bader` reparte y compara la suma de cuencas con la integral total; `report_bader` compara además la integral con $\sum_i Z_i^{\mathrm{val}}$ y avisa si difieren en más del 5 %.
4. `--difference total.cube parte1.cube …`: `charges.difference` exige rejillas idénticas y resta; `report_difference` da carga neta, carga acumulada y extremos del perfil planar (`fields.planar_average`, eje `--axis`); `plot_difference` dibuja el perfil.

**De dónde sale cada dato.**

| Dato | Origen | Detalle |
|---|---|---|
| Cargas de Löwdin, spilling | salida de `projwfc.x` | regex `_RE_LOWDIN`, `_RE_SPILL` |
| $\rho(\mathbf{r})$ | `.cube` de `pp.x` (`plot_num=0`, `output_format=6`) | `fields.read_cube` |
| Posiciones atómicas (Bader) | estructura `file` | `atoms.positions` |
| $Z^{\mathrm{val}}$ por átomo | `z_valence` de los UPF (`--pseudo-dir` o configuración) | `charges.valence_from_pseudos` → `pseudo.resolve` |
| $a_0$ (bohr → Å) | `fields.BOHR` | 0.529177210903 |
| Eje del perfil | `--axis` (0/1/2) | `fields.planar_average` |

**Límites y trampas.**
- Bader on-grid: *"Hereda el sesgo de malla del método (centésimas de electrón); para números finos usa la variante near-grid del código `bader` de Henkelman."* Aviso si la suma de cuencas difiere de la integral en más de 1e-3 e: *"la malla del cube es demasiado gruesa."*
- Si la integral de la rejilla no coincide con $\sum Z^{\mathrm{val}}$ (más del 5 %): *"Revisa que el cube sea la densidad de valencia completa (plot_num=0) y que los UPF de --pseudo-dir sean los del cálculo."* Un cube que ya esté en e/Å³ hay que declararlo con `density_units="e/A3"` (sólo desde Python; la CLI asume e/bohr³).
- Löwdin: spilling > 0.05 → *"AVISO: por encima de ~0.05 la base atómica no describe bien los estados"*. Sirve para comparar átomos, no como carga absoluta.
- Sin `--pseudo-dir` legible la columna "neta" sale `n/d` con aviso; los UPF deben ser los del cálculo, porque $Z^{\mathrm{val}}$ depende del pseudo (semicore o no).
- El perfil de $\Delta\rho$ se dibuja en e/bohr³ (unidades del cube), no en e/Å³.
- `--difference` exige la misma celda, malla FFT y cutoffs: *"las rejillas no coinciden… la resta no significa nada."*

**Referencias.**
- R. F. W. Bader, *Atoms in Molecules: A Quantum Theory* (Oxford, 1990).
- G. Henkelman, A. Arnaldsson, H. Jónsson, *Comput. Mater. Sci.* 36, 354 (2006) — Bader on-grid.
- P.-O. Löwdin, *J. Chem. Phys.* 18, 365 (1950).

---

### `olla-dft charge` — Campos escalares de pp.x y perfil planar

**Qué responde.** ¿Cómo se distribuye a lo largo de un eje la densidad de carga, la densidad de espín, la ELF o el potencial electrostático de un cálculo terminado?

**Fundamento para no expertos.** `pp.x` extrae de la función de onda y la densidad ya calculadas un campo escalar en la rejilla 3D. Promediarlo sobre los planos perpendiculares a un eje da un "perfil" 1D fácil de leer: dónde están las capas de una losa, dónde se acumula espín, dónde el vacío.

**Fórmulas.** `fields.planar_average`:
$$\bar f(z_k) = \frac{1}{n_1 n_2}\sum_{i,j} f(i,j,k),\qquad z_k = k\,|\mathbf{h}_3|$$
- $\mathbf{h}_3$: paso de la rejilla a lo largo del eje elegido (Å). Los otros ejes se obtienen permutando.

**Cómo lo calcula Olla-DFT.**
1. `_cmd_charge`: si no existe `<nombre>.cube` (o con `--rerun`), `fields.run_pp` escribe `pp_<campo>.in` con `plot_num` de `fields.PLOTS` (density 0, vtotal 1, spin 6, elf 8, potential 11), `iflag=3`, `output_format=6`, y ejecuta `pp.x` (buscado junto a `pw.x`); exige `JOB DONE`.
2. `fields.read_cube` lee origen, ejes (bohr → Å si $n>0$) y valores.
3. `fields.planar_average` a lo largo de `--axis` (a/b/c); se escribe `PERFIL_PLANAR.dat` y la figura `perfil_<nombre>`.

**De dónde sale cada dato.**

| Dato | Origen | Detalle |
|---|---|---|
| Campo 3D | `<nombre>.cube` de `pp.x` | unidades de `pp.x`: e/bohr³ (densidad), Ry (potenciales) |
| `prefix` | XML del cálculo | `qeout.read_xml(...).prefix` |
| `plot_num` | tabla `fields.PLOTS` | 0, 1, 6, 8, 11 |
| Bohr → Å | `qeout.BOHR_ANG` | 0.529177210903 |

**Límites y trampas.**
- El perfil se exporta en las unidades crudas del cube (no se convierte Ry → eV aquí; sí en `wf`).
- Necesita `pp.x` compilado (`make pp`); si falta: *"no se encontró pp.x junto a pw.x…"*.
- El comando no interpreta el campo: sólo lo promedia y lo dibuja. El `.cube` se abre en VESTA para isosuperficies.

**Referencias.** Documentación de `pp.x` (INPUT_PP, Quantum ESPRESSO).

---

### `olla-dft wf` — Función trabajo desde el nivel de vacío

**Qué responde.** ¿Cuánta energía cuesta sacar un electrón de una superficie al vacío? $\Phi = V_{\mathrm{vac}} - E_F$.

**Fundamento para no expertos.** En una losa con vacío, el potencial electrostático se aplana lejos del material: esa meseta es el "nivel de vacío", la energía de un electrón en reposo fuera del sólido. La función trabajo es la distancia desde el nivel de Fermi (el último nivel ocupado) hasta esa meseta. Si la meseta no es plana, o el vacío es corto, o la losa tiene un dipolo neto que inclina el potencial.

**Fórmulas.** `fields.work_function`:
$$\bar V(z) = \mathrm{RY\_EV}\cdot\overline{V_{\mathrm{pp}}}(z),\qquad V_{\mathrm{vac}} = \frac{1}{2h+1}\sum_{k=-h}^{h}\bar V\big(z_{i^\ast + k}\big),\qquad \Phi = V_{\mathrm{vac}} - E_F$$
$$\text{planitud} = \max_{k}\bar V - \min_{k}\bar V\ \text{en la misma ventana}$$
- La ventana de índices $\{i^\ast + k\}$ la da `fields.vacuum_window` cuando se conocen las posiciones atómicas (la CLI las pasa desde el XML): es el 20 % central del hueco más ancho **sin átomos** a lo largo del eje (medido con periodicidad en coordenadas fraccionarias), con $h = \max(2, 0.1\,f_{\mathrm{hueco}} N_z)$. Sin posiciones, se cae al criterio ciego: $i^\ast = \arg\max_z \bar V$ y $h = \max(2, N_z/10)$ (±10 % de la celda alrededor del máximo). $E_F$ en eV del XML; `RY_EV` = 13.605693122994.

**Cómo lo calcula Olla-DFT.**
1. `_cmd_wf`: si no existe `potencial.cube`, `fields.run_pp(path, "potential", ...)` ejecuta `pp.x` con `plot_num=11` ($V_{\mathrm{bare}} + V_H$).
2. `fields.read_cube` y `qeout.read_xml` (para `fermi`, del tag `fermi_energy` en Ha → eV).
3. `fields.work_function(cube, E_F, axis, positions=qe.positions)` promedia en el plano, localiza la meseta de vacío (`vacuum_window`) y calcula $\Phi$ y la planitud; el reporte dice en qué tramo de $z$ se evaluó.
4. `report_wf`, `export_wf` (`WF.dat` con cabecera `Phi_eV`, `V_vacio_eV`, `E_Fermi_eV`, `planitud_eV` y el perfil) y `plot_profile`.

**De dónde sale cada dato.**

| Dato | Origen | Detalle |
|---|---|---|
| $V(\mathbf{r})$ | `potencial.cube` (`pp.x`, `plot_num=11`, Ry) | `fields.read_cube` |
| $E_F$ | XML de `pw.x`, tag `fermi_energy` | `qeout.read_xml`, Ha → eV |
| Ry → eV | `qeout.RY_EV` | 13.605693122994 |
| Posiciones atómicas (para `vacuum_window`) | XML de `pw.x` (`atomic_positions`) | `qeout.read_xml(...).positions` |
| Eje | `--axis` (c por omisión) | `_AXES` |

**Límites y trampas.**
- Aviso si la planitud > 0.05 eV: *"la meseta de vacío varía más de 0.05 eV. El vacío es insuficiente o hay un dipolo neto; aumenta el vacío (o usa una losa simétrica)…"*.
- Sin posiciones (uso desde Python) la meseta se busca a ciegas alrededor del máximo del potencial: *"con poco vacío la ventana puede pisar la cola del potencial atómico"* (docstring). La CLI siempre pasa las posiciones del XML.
- No aplica corrección dipolar por sí mismo: una losa polar da dos niveles de vacío distintos y este comando toma el más alto. Para losas polares hay que generar el cálculo con `--dipole` (`gen`, `eform`) o usar `esm`.
- Si el XML no trae `fermi_energy` (ocupaciones fijas): *"el XML no trae energía de Fermi (¿terminó el scf?)"*.

**Referencias.**
- N. D. Lang, W. Kohn, *Phys. Rev. B* 3, 1215 (1971) — función trabajo en el modelo de jellium.
- L. Bengtsson, *Phys. Rev. B* 59, 12301 (1999) — corrección dipolar en losas.

---

### `olla-dft esm` — Superficies cargadas con medio de apantallamiento efectivo

**Qué responde.** ¿Cuál es la función trabajo, la capacitancia y el potencial de carga cero de una losa (neutra o cargada) sin que las imágenes periódicas ni el fondo compensador contaminen el resultado?

**Fundamento para no expertos.** Una losa cargada en una celda periódica es un problema mal planteado: QE reparte un fondo uniforme de carga opuesta por todo el volumen, vacío incluido, y la energía depende del tamaño de la celda sin converger a nada. El **ESM** (Effective Screening Medium) sustituye la periodicidad en $z$ por una condición de contorno explícita: se resuelve la ecuación de Poisson dentro de la celda y se empalma con una solución analítica fuera. Tres variantes: `bc1` (vacío a ambos lados, losas neutras; el nivel de vacío vale cero por construcción), `bc2` (dos placas metálicas: un condensador, admite campo) y `bc3` (vacío/metal: un electrodo que recibe la contracarga). Con `bc2`/`bc3` la distancia al electrodo ya no es un parámetro de convergencia sino **física**: fija la capacitancia.

**Fórmulas.** En `qekit/modules/esm.py`.

Centrado (`esm.centrar`): $z_i \leftarrow z_i - \tfrac{1}{2}(z_{\min}+z_{\max})$ (ESM mide $z$ desde el centro de la celda).

Nivel de vacío (`esm.nivel_vacio`): promedio de $V_{\mathrm{tot}}(z)$ del `.esm1` en la región $|z| > t/2 + m$, con $t$ el espesor de la losa y un margen $m$ que empieza en `MARGEN_VACIO` = 2 Å y crece de 0.5 en 0.5 Å (hasta `margen_max` = 8 Å) hasta que la desviación típica del potencial baje de `tol` = 1e-3 eV; con `bc3` sólo el lado $z<0$.

$$\Phi = V_{\mathrm{vac}} - E_F$$

Capacitancia (`esm.capacitancia`), ajuste lineal $q = C' V + b$:
$$C = \frac{dq}{dV}\,\frac{1}{A}\cdot 1.602176634\times10^{3}\quad[\mu\mathrm{F/cm^2}],\qquad R^2 = 1 - \frac{\sum(q-\hat q)^2}{\sum(q-\bar q)^2}$$
- $q$ en e por celda, $V$ en V (eV/e), $A$ = área de la celda en Å² (`|(\mathbf a\times\mathbf b)_z|`); `E_A2_A_UF_CM2` = $1.602176634\times10^{3}$ convierte e/(Å²·V) a µF/cm².

Linealidad (`esm.linealidad`): $\max|P - \hat P| / (\max P - \min P) \le$ `tol` = 0.02.

Potencial de carga cero (`esm.potencial_de_carga_cero`): interpolación lineal de $\Phi(q)$ en $q = 0$.

Gran canónico (`esm.gran_canonico`, sólo biblioteca): $\Omega = E + q\,\Phi$.

**Cómo lo calcula Olla-DFT.**
1. `esm.comprobar`: rechaza `bc1` con carga (*"bc1 es vacío por los dos lados… la energía diverge"*) y celdas no ortogonales en $z$; avisa si vacío < `VACIO_MINIMO` = 6 Å, si la losa no estaba centrada y, con `bc2/bc3` y carga, de que el vacío es física.
2. `esm.prepare` centra la losa, calcula espesor, vacío y área, y escribe un `scf` por carga en `q00/`, `q01/`… (`inputgen.build_pw_input`, `conv_thr=1e-8`, smearing `mv` con `degauss=0.02`, malla $n_1\times n_2\times 1$, `tot_charge=q`) e inserta en `&SYSTEM`: `assume_isolated='esm'`, `esm_bc`, `esm_nfit=--nfit` (4), `esm_w=--esm-w` si ≠ 0, `esm_efield=--field` sólo con `bc2`. Escribe `run.sh`.
3. `--run` o a mano: `pw.x` en cada carpeta.
4. `esm.collect` lee de cada carpeta el XML (`total_energy`, `fermi`) y el `<prefix>.esm1` (`esm.leer_esm1`: z (Å), carga (e/Å), $V_H$, $V_{\mathrm{loc}}$, $V_{\mathrm{tot}}$ en eV); `nivel_vacio` y $\Phi$.
5. `esm.report`: tabla $q, E, E_F, V_{\mathrm{vac}}, \Phi$; con `bc1` comprueba $|V_{\mathrm{vac}}| < 10^{-3}$ eV; con varias cargas, capacitancia de $V_{\mathrm{vac}}(q)$ (voltaje de la celda) y, si $\Phi(q)$ es lineal, también de $\Phi(q)$ con PZC.
6. `export` (`ESM.dat`, `ESM_perfil_qNN.dat`, `ESM.txt`) y `plot` (perfiles y $q$ vs $\Phi$).

**De dónde sale cada dato.**

| Dato | Origen | Detalle |
|---|---|---|
| $V_{\mathrm{tot}}(z)$, carga$(z)$ | `<prefix>.esm1` escrito por `pw.x` con ESM | `esm.leer_esm1`, columnas 0–4 |
| $E$, $E_F$ | XML de `pw.x` | `qeout.read_xml` (Ha → eV) |
| Área $A$ | vectores $\mathbf a,\mathbf b$ de la celda | `esm.prepare` |
| Factor µF/cm² | `esm.E_A2_A_UF_CM2` | $e/(10^{-8}\,\mathrm{cm})^2$ |
| Cargas | `--charge` (lista) | e por celda |
| Campo | `--field` (Ry/u.a.) | sólo `bc2` |

**Límites y trampas.**
- *"Con bc2 o bc3 la capacitancia depende de la distancia al contraelectrodo: es una capacitancia DE ESTE MONTAJE, no una propiedad del material."*
- Las energías con carga neta no son comparables entre sí: *"la energía de ESM incluye la interacción con la carga imagen del electrodo, que crece como q²."*
- Si $\Phi(q)$ no es lineal: *"Φ(q) = V_vac − E_F NO es una recta… no doy un potencial de carga cero sobre ella."*
- `gran_canonico` (Ω = E + qΦ) existe en el módulo pero **ningún comando lo usa**; el "gran canónico" del título del módulo no está expuesto en la CLI.
- La losa se centra automáticamente; si el usuario ya la centró en $c/2$ (ASE) el aviso explica por qué se recentró.
- El cálculo siempre usa smearing (`insulator=False`): pensado para metales/electrodos.

**Referencias.**
- M. Otani, O. Sugino, *Phys. Rev. B* 73, 115407 (2006) — ESM.
- N. Bonnet, T. Morishita, O. Sugino, M. Otani, *Phys. Rev. Lett.* 109, 266101 (2012) — potencial constante con ESM.

---

### `olla-dft echem` — Electrodo de hidrógeno computacional: HER y OER

**Qué responde.** ¿Qué potencial hay que aplicar para que todos los pasos de la evolución de hidrógeno (HER) o de oxígeno (OER) sean cuesta abajo, y cuánto se aleja del potencial de equilibrio (sobrepotencial)?

**Fundamento para no expertos.** Calcular un protón solvatado es un problema durísimo. El truco del electrodo de hidrógeno computacional (CHE) es notar que, a 0 V frente al electrodo estándar de hidrógeno y pH 0, el par $\mathrm{H^+ + e^-}$ tiene la misma energía libre que $\tfrac12\mathrm{H_2(g)}$, que sí se calcula. Cada paso que libera un $(\mathrm{H^+ + e^-})$ se evalúa así, y el potencial $U$ y el pH entran después como términos que se suman. El paso con mayor $\Delta G$ es el "limitante": el potencial que lo hace exergónico es el potencial limitante, y su distancia al de equilibrio es el sobrepotencial. Es termodinámica de intermedios: no hay barreras cinéticas ni disolvente.

**Fórmulas.** En `qekit/modules/echem.py`.

Dependencia con $U$ y pH (`Echem.dG`):
$$\Delta G_i(U, \mathrm{pH}) = \Delta G_i(0,0) - eU - k_B T\ln 10\cdot\mathrm{pH} = \Delta G_i(0,0) - e\,U_{\mathrm{RHE}}$$
$$U_{\mathrm{RHE}} = U_{\mathrm{SHE}} + k_B T\ln 10\cdot\mathrm{pH}\quad(\text{`echem.u_rhe`; } 0.0592\,\mathrm{pH\ V\ a\ 298\ K})$$
- $k_B$ = `KB_EV` = $8.617333262\times10^{-5}$ eV/K; $T$ = `--temperature` (298.15 K); $U$ = `-U` en V **frente al SHE** (a pH 0 coincide con RHE); el término de pH es exactamente la conversión SHE → RHE, así que en la escala RHE los $\Delta G$ no dependen del pH. Un electrón por paso.

HER (`echem.her`):
$$\Delta G_{\mathrm{H^*}} = E_{\mathrm{ads}}(\mathrm{H}) + c_{\mathrm{H}},\qquad \text{pasos: } (+\Delta G_{\mathrm{H^*}},\ -\Delta G_{\mathrm{H^*}})$$
- $E_{\mathrm{ads}}(\mathrm{H})$: `--her`, referida a $\tfrac12\mathrm{H_2}$ (eV); $c_{\mathrm{H}}$ = ZPE − TΔS = 0.24 eV por omisión (`CORRECCIONES`).

OER (`echem.oer`), con $G_X = E_{\mathrm{ads}}(X) + c_X$:
$$\Delta G_1 = G_{\mathrm{OH}},\quad \Delta G_2 = G_{\mathrm{O}} - G_{\mathrm{OH}},\quad \Delta G_3 = G_{\mathrm{OOH}} - G_{\mathrm{O}},\quad \Delta G_4 = 4.92\ \mathrm{eV} - (\Delta G_1+\Delta G_2+\Delta G_3)$$
- $c_{\mathrm{OH}} = 0.35$, $c_{\mathrm{O}} = 0.05$, $c_{\mathrm{OOH}} = 0.40$ eV por omisión; `DG_AGUA_TOTAL` = 4.92 eV (experimental, $2\mathrm{H_2O} \to \mathrm{O_2} + 2\mathrm{H_2}$).

Potencial limitante y sobrepotencial (`Echem.U_limitante`, `Echem.sobrepotencial`):
$$U_L = \max_i \Delta G_i(0,0)/e,\qquad \eta = U_L - U_{\mathrm{eq}},\quad U_{\mathrm{eq}}^{\mathrm{OER}} = 1.229\ \mathrm{V},\ U_{\mathrm{eq}}^{\mathrm{HER}} = 0$$
- $\eta$ se devuelve **con signo**: positivo = en $U_{\mathrm{eq}}$ el paso limitante sigue cuesta arriba (con los perfiles de aquí nunca sale negativo; sólo podría con un `dG_total` distinto del experimental).

Relación de escala (`echem.escala_ooh_oh`, OER) y su límite (`echem.sobrepotencial_minimo_escala`):
$$\Delta_{\mathrm{esc}} = G_{\mathrm{OOH}} - G_{\mathrm{OH}}\ \text{(comparada con `ESCALA_OOH_OH` = 3.2 ± 0.2 eV)},\qquad \eta_{\min} = \frac{\Delta_{\mathrm{esc}}}{2} - \frac{\Delta G_{\mathrm{total}}}{4} = 0.37\ \mathrm{V}$$
- Si OOH* y OH* están separados por $\Delta_{\mathrm{esc}}$ fijo, los pasos 2 y 3 suman $\Delta_{\mathrm{esc}}$ y el peor no baja de $\Delta_{\mathrm{esc}}/2$ = 1.6 eV; frente a 4.92/4 = 1.23 V quedan ~0.37 V.

Rejilla tipo Pourbaix (`echem.pourbaix`, sólo biblioteca): $\Delta G_{\lim}(U,\mathrm{pH}) = \max_i\Delta G_i(0,0) - eU - k_BT\ln10\cdot\mathrm{pH}$ sobre $U\in[-0.5,2]$ V y pH $\in[0,14]$.

**Cómo lo calcula Olla-DFT.**
1. `_cmd_echem` exige exactamente una de `--her E` o `--oer OH=..,O=..,OOH=..`; `--corrections X=eV` sobrescribe las correcciones térmicas.
2. `echem.her` o `echem.oer` arman la lista de pasos $(\text{nombre}, \Delta G_i)$; `oer` avisa si $\Delta G_4 < 0$ y si se usaron correcciones de tabla.
3. Se fijan `U` (vs SHE) y `pH` del usuario; `echem.report` imprime también $U_{\mathrm{RHE}}$ (`Echem.U_rhe`) si pH ≠ 0, tabula $\Delta G(0)$ y $\Delta G(U,\mathrm{pH})$, el paso limitante, $U_L$ (vs RHE), $\eta$ con signo, el descriptor $\Delta G_{\mathrm{H^*}}$ (HER) o la relación de escala y el $\eta_{\min}$ que impone (OER).
4. `export` escribe `ECHEM.dat` y `ECHEM.txt`; `plot` dibuja el diagrama en escalera a $U = 0$, $U_{\mathrm{eq}}$ y $U_L$.

**De dónde sale cada dato.**

| Dato | Origen | Detalle |
|---|---|---|
| $E_{\mathrm{ads}}$ de H, OH, O, OOH | parámetros `--her`, `--oer` | eV, referidas a H₂O y ½H₂ (de `adsorb`) |
| Correcciones ZPE − TΔS | `echem.CORRECCIONES` o `--corrections` | H 0.24, OH 0.35, O 0.05, OOH 0.40 eV (Nørskov y col.) |
| $\Delta G$ total del agua | `echem.DG_AGUA_TOTAL` | 4.92 eV, experimental |
| $U_{\mathrm{eq}}$ | `echem.U_EQ_OER`, `U_EQ_HER` | 1.229 V, 0 V |
| $k_B$ | `echem.KB_EV` | $8.617333262\times10^{-5}$ eV/K (CODATA) |
| $\Delta_{\mathrm{esc}}$ universal | `echem.ESCALA_OOH_OH` | 3.2 eV (Man et al. 2011) |

**Límites y trampas.**
- *"El CHE es termodinámica de intermedios: NO hay barreras cinéticas, ni disolvente explícito, ni doble capa."*
- `-U` es frente al **SHE** (ayuda de la CLI: *"a pH 0 es el mismo que frente al RHE; el pH lo convierte"*); $U_L$ y $\eta$ están en la escala RHE. Para la HER $U_L = |\Delta G_{\mathrm{H^*}}|$, así que $\eta \ge 0$ siempre.
- Cuarto paso por diferencia: *"El cuarto paso sale NEGATIVO… o hay un error en las referencias, o tu superficie liga los intermedios muchísimo."*
- `pourbaix()` no está conectada a ningún comando: el "diagrama de Pourbaix" del título del módulo no se produce desde la CLI.

**Referencias.**
- J. K. Nørskov, J. Rossmeisl, A. Logadottir, L. Lindqvist, J. R. Kitchin, T. Bligaard, H. Jónsson, *J. Phys. Chem. B* 108, 17886 (2004) — CHE. DOI: 10.1021/jp047349j.
- J. K. Nørskov, T. Bligaard, A. Logadottir, J. R. Kitchin, J. G. Chen, S. Pandelov, U. Stimming, *J. Electrochem. Soc.* 152, J23 (2005) — volcán de HER.
- I. C. Man et al., *ChemCatChem* 3, 1159 (2011) — relación de escala de la OER.

---

### `olla-dft adsorb` — Sitios de adsorción y energía de adsorción

**Qué responde.** ¿En qué sitios no equivalentes de una superficie puede posarse una molécula, y cuánto gana (o pierde) el sistema al hacerlo en cada uno?

**Fundamento para no expertos.** Una molécula sobre un metal se posa encima de un átomo (*top*), sobre el punto medio entre dos (*bridge*) o sobre el centro de un triángulo de átomos (*hollow*; en fcc(111) hay dos: con o sin átomo debajo en la segunda capa). Muchos de esos sitios son copias por simetría, así que se agrupan por su "huella": la lista ordenada de distancias a sus 24 vecinos más cercanos contando todas las capas. La energía de adsorción es una resta de tres energías totales que sólo tiene sentido si los tres cálculos comparten celda, cutoffs, malla k y pseudos; por eso se generan juntos.

**Fórmulas.** `thermochem.adsorcion` (llamada desde `AdsorbRun.energias_ads`):
$$E_{\mathrm{ads}} = E(\text{losa}+\text{mol}) - E(\text{losa}) - n\,E(\text{mol})$$
- Todas en eV; $n$ = número de moléculas (`n_mol`, 1). Negativa = favorable.

Geometría tras relajar (`adsorb.collect`):
$$h = \min_{a\in\mathrm{ads}} z_a - \max_{s\in\mathrm{losa}} z_s,\qquad d_{\mathrm{contacto}} = \min_{a,s}|\mathbf r_a - \mathbf r_s|$$

Huella de un sitio (`adsorb._huella`): distancias ordenadas a los $k$ = `N_VECINOS_HUELLA` = 24 átomos más cercanos (con réplicas periódicas); dos sitios son el mismo si $\max|\Delta d| <$ `TOL_HUELLA` = 0.05 Å.

**Cómo lo calcula Olla-DFT.**
1. `adsorb.prepare` exige vacío en $c$ (`kpoints.direcciones_con_vacio`) y carga la molécula (`cargar_molecula`: archivo o base G2 de ASE).
2. `adsorb.sitios`: capa expuesta = átomos a menos de `TOL_CAPA` = 0.6 Å del $z$ extremo; *top* sobre cada uno; *bridge* entre pares a menos de `R_VECINO` = 3.6 Å; *hollow* en los baricentros de la triangulación de Delaunay (se descartan triángulos con lado > 1.6·3.6 Å); se llevan a la celda y se deduplican por huella; se etiquetan `top1`, `bridge1`, `hollow1`…
3. Con `--rotations N` y molécula poliatómica, cada sitio se repite con giros de $360k/N$ grados alrededor de $z$.
4. `sweep.prepare_common` se resuelve sobre la **unión** losa + molécula (mismos pseudos y cutoffs para todo). Se escriben `_losa/`, `_molecula/` (molécula centrada en la **misma** celda) y una carpeta por sitio (`adsorb.colocar`: átomo `--anchor` a `--height` = 2.0 Å sobre el sitio), todos `relax` salvo `--fixed-ions`; `run.sh`. Con `--dipole`, `dipole_correction=3` entra en los **tres** cálculos (`inputgen.build_pw_input`: `tefield`, `dipfield`, `edir=3`, `emaxpos`/`eopreg` en el centro del hueco de vacío por `inputgen._region_vacio`, `eamp=0`); si no hay ≥ 5 Å de vacío, aborta.
5. `--run`/`--collect`: `adsorb.collect` lee los XML (`qeout.read_xml`), energías, convergencia, altura y contacto.
6. `adsorb.report`: tabla ordenada por $E_{\mathrm{ads}}$, mejor sitio, diagnóstico por rangos (>0: no se pega; > −0.30 eV: fisisorción débil; < −2 eV: probable reacción/disociación o quimisorción atómica), diferencia con el segundo (< 50 meV: indistinguibles). `export`: `ADSORCION.dat/.txt`; `plot`: barras.

**De dónde sale cada dato.**

| Dato | Origen | Detalle |
|---|---|---|
| $E$(losa), $E$(mol), $E$(losa+mol) | XML de `pw.x` en `_losa/`, `_molecula/`, `<sitio>/` | `total_energy` (Ha → eV) |
| Posiciones relajadas | XML (`atomic_positions`) | altura y contacto |
| Molécula | archivo o `ase.build.molecule` | `--mol` |
| Radio de vecinos, tolerancias | `adsorb.R_VECINO`, `TOL_CAPA`, `N_VECINOS_HUELLA`, `TOL_HUELLA` | 3.6 Å, 0.6 Å, 24, 0.05 Å |
| Altura inicial | `--height` | 2.0 Å |
| Corrección vdW | `--vdw` | pasa a `inputgen.build_pw_input` |
| Corrección dipolar | `--dipole` | `dipole_correction=3` en losa, molécula y losa+molécula |

**Límites y trampas.**
- Sin `--vdw`: *"AVISO: sin corrección de van der Waals. En fisisorción… la energía sale cerca de cero y la geometría desligada."*
- Sin `--dipole` en cara `top`: *"Sugerencia: una molécula adsorbida en una sola cara deja la losa polar. Con --dipole se cancela el dipolo artificial a través del vacío."* La sierra se pone en los tres cálculos a propósito: *"si la referencia se calcula sin corregir, la resta arrastra el error."*
- $E_{\mathrm{ads}} > 0$ con iones fijos: *"lo más probable es que la altura inicial… no sea la de equilibrio y estés midiendo la repulsión."*
- La referencia es la molécula tal cual se pasó: con `--mol H` la referencia es el **átomo**, no ½H₂ (el reporte lo advierte para $|E_{\mathrm{ads}}| > 2$ eV).
- La molécula aislada se calcula con la misma malla k que la losa (consistencia deliberada, no una caja aparte).
- La enumeración de sitios es geométrica: no detecta sitios sobre segundas capas ni reconstrucciones.

**Referencias.**
- B. Hammer, J. K. Nørskov, *Adv. Catal.* 45, 71 (2000) — adsorción en superficies metálicas.
- S. Grimme, J. Antony, S. Ehrlich, H. Krieg, *J. Chem. Phys.* 132, 154104 (2010) — DFT-D3.

---

### `olla-dft surface` — Cortar una losa (hkl) con vacío

**Qué responde.** Dado un cristal, ¿cómo es la losa de superficie $(hkl)$ con $N$ capas y vacío, es simétrica, es polar, y cuánto vacío real queda?

**Fundamento para no expertos.** Una superficie se simula con una "losa": unas cuantas capas atómicas paralelas al plano $(hkl)$ y, encima, vacío suficiente para que la losa no vea su copia periódica. Si las dos caras no son iguales (losa *polar*), aparece un dipolo artificial a través del vacío que desplaza las funciones trabajo; QE lo corrige con `dipfield`. El vacío que importa es el que hay entre átomos, no entre bordes de celda.

**Fórmulas.** `builder.surface`:
$$t = z_{\max} - z_{\min},\qquad v_{\mathrm{real}} = c - t$$
- Simétrica: el perfil ordenado $z_i - \bar z$ coincide con su reflejo dentro de `tol` = 0.3 Å. Polar: la composición de la capa superior ≠ la de la inferior (átomos a menos de `tol` del extremo).

**Cómo lo calcula Olla-DFT.**
1. `structure.conventional` → `ase.build.surface(base, miller, layers, vacuum=vacuum/2, periodic=True)` y `slab.center(vacuum=vacuum/2, axis=2)`.
2. `builder.surface` calcula grosor, vacío real, número de planos atómicos (`_planos_z`, tolerancia 0.3 Å), simetría y polaridad; con `--fix N` marca los átomos de los $N$ planos inferiores de dos formas (`_fijar_capas`): el array `slab.arrays['qekit_fijo']` y una restricción `FixAtoms` de ASE. `inputgen.fixed_atoms` lee cualquiera de las dos y escribe `0 0 0` en la tercera columna de `ATOMIC_POSITIONS`.
3. Avisos: > 1.5 átomos por plano (celda múltiple de la mínima), vacío real < 10 Å, losa polar, < 4 capas, congelar todos los planos.
4. `report_slab` y, con `-o`, `structure.convert` escribe CIF/POSCAR/XYZ. Si hay átomos fijos y el formato no los conserva (`structure.conserva_fijos`: sólo POSCAR/CONTCAR/`.vasp` los guardan como *Selective dynamics*), la CLI avisa y recomienda `builder.FORMATO_CON_FIJOS` (POSCAR o `.vasp`) o `olla-dft gamma --fix`.

**De dónde sale cada dato.**

| Dato | Origen | Detalle |
|---|---|---|
| Losa | `ase.build.surface` sobre la celda convencional | `--miller`, `--layers` (6), `--vacuum` (15 Å) |
| Celda convencional | spglib vía `structure.conventional` | referencia de los índices hkl |
| Planos atómicos | alturas $z$ distintas (tol 0.3 Å) | `builder._planos_z` |

**Límites y trampas.**
- *"la losa es POLAR… Añade 'dipfield = .true.' y 'edir = 3' al input, o corta una losa simétrica."*
- `--fix` se pierde al exportar a CIF o XYZ: *"el CIF no tiene dónde ponerlo, así que al volver a cargarlo se relajaría todo. Escribe la losa en POSCAR (o .vasp)…"*. Sólo POSCAR conserva la restricción `FixAtoms`, que `inputgen.fixed_atoms` traduce a `0 0 0`.
- La detección de polaridad compara sólo composiciones de las capas extremas: una losa con terminaciones de igual composición pero geometría distinta no se marca.
- El corte sobre la celda convencional puede dar una celda superficial mayor que la mínima (se avisa).

**Referencias.**
- P. W. Tasker, *J. Phys. C* 12, 4977 (1979) — superficies polares.
- ASE: A. H. Larsen et al., *J. Phys.: Condens. Matter* 29, 273002 (2017).

---

### `olla-dft defect` — Construir un defecto puntual

**Qué responde.** ¿Cómo son la supercelda perfecta y la supercelda con una vacancia, una sustitución o un intersticial, y cuál es la fórmula de energía de formación que habrá que llenar?

**Fundamento para no expertos.** Un defecto puntual se modela repitiendo la celda primitiva $n_1\times n_2\times n_3$ veces y modificando un átomo. La supercelda debe ser grande para que el defecto no interactúe con sus imágenes periódicas. Este comando sólo construye las dos estructuras y escribe la fórmula con sus términos; `eform` hace el cálculo.

**Fórmulas.** `builder.formation_energy_text` escribe:
$$E_f = E(\text{defecto}) - E(\text{perfecto}) \pm \mu(\cdot)\ \ [+\,q(E_F + E_v) + E_{\mathrm{corr}}]$$
- vacancia: $+\mu(\text{especie que sale})$; sustitución: $+\mu(\text{sale}) - \mu(\text{entra})$; intersticial: $-\mu(\text{entra})$.

**Cómo lo calcula Olla-DFT.**
1. `structure.primitive` → `repeat(supercell)` (por omisión 2×2×2).
2. `builder.defect`: vacancia (`del d[site]`), sustitución (`d[site].symbol = new`), intersticial (posición fraccionaria `--position` de la supercelda; avisa si queda a < 1.0 Å de un vecino, distancia con imagen mínima).
3. Aviso si el lado más corto de la supercelda < 10 Å.
4. `report_defect` y escritura de `perfecto.cif` y `defecto.cif` en `--outdir`.

**De dónde sale cada dato.**

| Dato | Origen | Detalle |
|---|---|---|
| Celda primitiva | spglib vía `structure.primitive` | base de la supercelda |
| Sitio, especie, posición | `--site`, `--new-element`, `--position` | índices base 0 en la supercelda |

**Límites y trampas.** *"la supercelda mide X Å en su lado más corto: el defecto se ve con sus imágenes periódicas. Para energías de formación conviene ≥ 10-12 Å."* No relaja nada ni calcula energías; el índice `--site` se refiere a la supercelda repetida, no al cristal de entrada.

**Referencias.** C. Freysoldt, B. Grabowski, T. Hickel, J. Neugebauer, G. Kresse, A. Janotti, C. G. Van de Walle, *Rev. Mod. Phys.* 86, 253 (2014).

---

### `olla-dft eform` — Energía de formación de defectos cargados

**Qué responde.** ¿Cuánto cuesta formar el defecto en cada estado de carga, cómo varía con el nivel de Fermi, dónde están los niveles de transición de carga y cuál es la corrección por tamaño finito?

**Fundamento para no expertos.** Formar un defecto cuesta energía que depende de tres cosas: de dónde vienen o adónde van los átomos (potencial químico $\mu$, fijado por las condiciones de síntesis), de dónde vienen o adónde van los electrones (nivel de Fermi $\varepsilon_F$, medido desde el máximo de la banda de valencia) y de un artefacto: una celda cargada periódica interacciona con sus propias imágenes y con el fondo neutralizante que QE añade. Ese artefacto se corrige con la energía electrostática de una carga puntual en una red de cargas imagen (Makov–Payne) apantallada por la constante dieléctrica, o con la versión de Lany–Zunger que incluye un término de forma. El punto donde dos rectas $E_f(q)$ se cruzan es un nivel de transición: el nivel de Fermi al que el defecto cambia de carga.

**Fórmulas.** En `qekit/modules/defects.py`.

Energía de formación (`DefectRun.E_f`):
$$E_f[D^q](\varepsilon_F) = E[D^q] - E[\mathrm{perf}] - \sum_i n_i\mu_i + q\,(\varepsilon_{\mathrm{VBM}} + \varepsilon_F) + E_{\mathrm{corr}}(q) + q\,\Delta V$$
- $n_i$: átomos **añadidos** de la especie $i$ (−1 para la que sale); $\mu_i$: `--mu EL=eV` (para un cristal elemental, $\mu = E[\mathrm{perf}]/N$ automáticamente, `asignar_mu_elemental`); $\varepsilon_{\mathrm{VBM}}$ = `highestOccupiedLevel` de la supercelda perfecta (eV); $\varepsilon_F \in [0, E_g]$; $\Delta V$: alineamiento de potencial (`--dv` o `--align`).

Constante de Madelung por Ewald (`defects.madelung_xi`, `constante_madelung`):
$$\xi = \sum_{\mathbf R\neq 0}\frac{\mathrm{erfc}(\eta R)}{R} + \frac{4\pi}{V}\sum_{\mathbf G\neq 0}\frac{e^{-G^2/4\eta^2}}{G^2} - \frac{2\eta}{\sqrt\pi} - \frac{\pi}{\eta^2 V},\qquad \alpha_M = -\xi\,L,\quad L = V^{1/3}$$
- $\eta = \sqrt\pi / V^{1/3}$; cortes reales y recíprocos ajustados a `tol` = 1e-10. Da $\alpha_M = 2.8372974$ para la red cúbica simple.

Corrección de imagen (`defects.correccion_imagen`):
$$E_{\mathrm{MP}} = \frac{k_e\,q^2\,\alpha_M}{2\,\varepsilon\,L},\qquad E_{\mathrm{LZ}} = E_{\mathrm{MP}}\left[1 + c_{\mathrm{sh}}\left(1 - \frac{1}{\varepsilon}\right)\right]$$
- $k_e$ = `KE` = 14.399645 eV·Å; $\varepsilon$ = `--epsilon`; $c_{\mathrm{sh}}$ = `C_SHAPE` = −0.35 (valor único; LZ dan −0.369 sc, −0.343 fcc, −0.342 bcc). `--correction` ∈ {`ninguna`, `makov-payne`, `lany-zunger`}.

Alineamiento (`defects.alineamiento`): $\Delta V = f\,\langle \bar V_{\mathrm{def}}(z) - \bar V_{\mathrm{perf}}(z)\rangle$ promediado en el 25 % de la celda opuesto al punto de mayor $|\Delta V - \mathrm{mediana}|$, con su desviación típica; $f$ = `UNIDADES_POTENCIAL[unidades_cube]` convierte los cubes de `pp.x` (`plot_num=11`, Ry, `unidades_cube="Ry"` por omisión, $f$ = `RY_EV`) a eV; con `"eV"`, $f = 1$. El resultado (`dV`, `sigma`, `perfil`) está siempre en eV.

Niveles de transición (`defects.niveles_transicion`), una entrada por cada par de cargas consecutivas $a<b$ (ordenadas por $q$), con la bandera `dentro` = $0 \le \varepsilon \le E_g$; **no** se filtra por la envolvente inferior (para los niveles observables hay que cruzar con `envolvente`):
$$\varepsilon(a/b) = \frac{E_f(a, 0) - E_f(b, 0)}{b - a}$$

**Cómo lo calcula Olla-DFT.**
1. `defects.prepare`: exige `--epsilon` si hay cargas ≠ 0 y corrección ≠ `ninguna`; construye las celdas con `builder.defect`; resuelve pseudos sobre la unión de especies.
2. Paridad: si `--insulator` y algún estado de carga deja un número impar de electrones (`defects.electrones` con `z_valence` de los UPF), activa `nspin=2` en **todos** los estados con `tot_magnetization` 1 (impares) o 0 (pares).
3. Escribe `_perfecto/` (scf) y `qm1/`, `qp0/`, `qp1/`… (`relax` salvo `--fixed-ions`, `tot_charge=q`, mismo `nbnd` estimado) y `run.sh`.
4. `--run`/`--collect`: `defects.collect` lee energías, convergencia, `homo` (VBM) y `lumo` (gap) de la perfecta; `--mu`; `--align POT_DEF POT_PERF` o `--dv`.
5. `report`: tabla $q$, $E$, $E_{\mathrm{corr}}$, $E_f(\varepsilon_F=0)$, $E_f(\varepsilon_F=E_g)$; niveles de transición marcando los que caen fuera del gap; envolvente inferior (`envolvente`) y cargas estables al recorrer el gap. `export`: `FORMACION.dat` (tabla y $E_f(\varepsilon_F)$ en 51 puntos); `plot`: $E_f$ vs $\varepsilon_F$.

**De dónde sale cada dato.**

| Dato | Origen | Detalle |
|---|---|---|
| $E[D^q]$, $E[\mathrm{perf}]$ | XML de `pw.x` | `total_energy` |
| $\varepsilon_{\mathrm{VBM}}$, $E_g$ | XML de la supercelda perfecta | `highestOccupiedLevel`, `lowestUnoccupiedLevel` |
| $\mu_i$ | `--mu` o $E[\mathrm{perf}]/N$ (elemental) | eV/átomo |
| $\varepsilon$ | `--epsilon` | p. ej. $\varepsilon_1(0)$ de `optics` |
| $\alpha_M$ | suma de Ewald sobre la celda real | `defects.madelung_xi` |
| $k_e$, $c_{\mathrm{sh}}$ | `defects.KE`, `defects.C_SHAPE` | 14.399645 eV·Å, −0.35 |
| $\Delta V$ | dos `.cube` de potencial (`--align`, Ry → eV) o `--dv` (ya en eV) | `defects.alineamiento`, `UNIDADES_POTENCIAL` |
| Electrones por celda | `z_valence` de los UPF | `defects.electrones` |

**Límites y trampas.**
- Sin `--epsilon`: *"la constante dieléctrica es lo que apantalla la interacción del defecto con sus imágenes; sin ella la corrección sale ε veces de más."*
- Con `--correction ninguna`: *"SIN CORREGIR: las E_f de los estados cargados están sistemáticamente bajas, y el error crece con q²."*
- `--dv` se da directamente en eV (no se convierte); `--align` asume cubes de `pp.x` en Ry y el reporte lo dice: *"entra en E_f como q·ΔV = … eV por unidad de carga (el potencial de pp.x viene en Ry y se pasó a eV)"*. Si $\sigma_{\Delta V} > 0.3\,|\Delta V|$: *"el defecto todavía se nota en la zona 'lejana', o sea que la supercelda es pequeña."*
- Los niveles de transición listados incluyen cruces entre estados que nunca son los más estables; el reporte marca *"<< fuera del gap"* los que caen fuera de $[0, E_g]$, pero un nivel dentro del gap entre dos estados que no están en la envolvente tampoco es observable.
- Sin VBM (metal, sin bandas vacías): *"No pude leer el VBM… E_f de los estados cargados no está definida."*
- $\mu$ ausente en un compuesto: *"FALTA el potencial químico… las DIFERENCIAS entre cargas y los niveles de transición sí valen, el valor absoluto de E_f no."*
- La corrección sólo quita el término principal $\propto q^2/L$; lado < 10 Å con carga: aviso.

**Referencias.**
- G. Makov, M. C. Payne, *Phys. Rev. B* 51, 4014 (1995).
- S. Lany, A. Zunger, *Phys. Rev. B* 78, 235104 (2008); *Modelling Simul. Mater. Sci. Eng.* 17, 084002 (2009).
- C. Freysoldt, J. Neugebauer, C. G. Van de Walle, *Phys. Rev. Lett.* 102, 016402 (2009).
- C. Freysoldt et al., *Rev. Mod. Phys.* 86, 253 (2014). DOI: 10.1103/RevModPhys.86.253.

---

### `olla-dft interface` — Heteroestructuras y desajuste de red

**Qué responde.** ¿Qué supercelda común permite apilar dos materiales 2D (o dos losas) con la menor deformación posible, cuánto vale esa deformación y cómo queda la estructura inicial?

**Fundamento para no expertos.** Dos redes cristalinas casi nunca encajan. Para ponerlas en la misma celda periódica hay que buscar múltiplos enteros de los vectores de cada una que se parezcan y estirar una de las dos. Esa deformación es el número que decide si el cálculo describe el material o una versión estirada de él: 1 % es tolerable, 8 % ya es otro material.

**Fórmulas.** En `qekit/modules/interface.py`.

Superceldas candidatas (`_celdas_candidatas`): $\mathbf A' = M\mathbf a$, $\mathbf B' = N\mathbf b$ con $M, N \in \mathbb Z^{2\times2}$, $|M_{ij}|,|N_{ij}| \le$ `--max-index` (4), $\det > 0$, agrupadas por determinante (las áreas deben coincidir dentro de $2\cdot$`tol`).

Deformación (`_deformacion`):
$$\boldsymbol\epsilon = B'^{-1}A' - I,\qquad \epsilon_{\max} = \max_{ij}|\epsilon_{ij}| \le \texttt{--tol}\ (0.05)$$

Reducción de Lagrange–Gauss (`reducir_2d`) para no repetir la misma red con bases distintas; desempate por "simplicidad" de $M, N$ (`_simplicidad`: suma de |entradas|, máximo, negativos, no nulos).

Separación inicial (`separacion_vdw`): $d_0 = 0.85\,(r_1 + r_2)$ con radios de van der Waals de `R_VDW` (Bondi; 2.0 Å si falta).

Con `--strain both`: celda objetivo $= (w A' + v B')/(w+v)$ con $w = n_1\,|\det \mathbf a|$, $v = n_2\,|\det\mathbf b|$.

**Cómo lo calcula Olla-DFT.**
1. `interface.buscar`: enumera, filtra por átomos (`--max-atoms` 200) y deformación, deduplica por $(n_1, n_2, \text{forma reducida}, \epsilon_{\max})$, ordena por $(\epsilon_{\max}, N_{\mathrm{at}}, \text{simplicidad})$ y devuelve las `--top` (10) mejores. `--list` sólo las imprime.
2. `interface.emparejar` elige `--index` y `construir`: `ase.build.make_supercell` para cada material, se lleva la celda en el plano a la objetivo arrastrando posiciones fraccionarias (`_supercelda_deformada`), se apila el material 2 a `--separation` (o $d_0$) sobre el 1, se aplica `--shift` (fracciones de la celda común), se añade `--vacuum` (20 Å) y se centra.
3. Avisos: $\epsilon_{\max} > 3\,\%$, separación de vdW como punto de partida, registro no optimizado.
4. `export`: `<name>.cif` y `<name>.txt`.

**De dónde sale cada dato.**

| Dato | Origen | Detalle |
|---|---|---|
| Vectores en el plano | celdas de `file1`, `file2` (filas 0–1, columnas 0–1) | `interface._plano` |
| Radios de vdW | tabla `interface.R_VDW` | Å; `R_VDW_DEFECTO` = 2.0 |
| Límites de búsqueda | `--max-index`, `--tol`, `--max-atoms` | 4, 0.05, 200 |

**Límites y trampas.**
- Se reporta la **componente mayor** $\max|\epsilon_{ij}|$ de la matriz, no una norma ni un promedio: *"una deformación de 0 % en una dirección y 6 % en la otra no es '3 %'."*
- *"La deformación es del X %. Por encima de ~3 % no se está modelando el material sino una versión estirada de él."*
- La separación es un punto de partida: *"con un funcional sin corrección de dispersión la distancia de equilibrio saldrá demasiado grande."*
- *"El REGISTRO… no está optimizado. Dos apilamientos distintos pueden diferir en decenas de meV por átomo."*
- Se asume que $c$ es la normal y que la celda es una losa; la deformación real con `--strain both` no coincide con la $\boldsymbol\epsilon$ reportada (que es la de llevar B a A).

**Referencias.**
- A. Bondi, *J. Phys. Chem.* 68, 441 (1964) — radios de van der Waals.
- P. Lazić, *Comput. Phys. Commun.* 197, 324 (2015) — CellMatch, emparejamiento de redes.

---

### `olla-dft neb` — Barreras de reacción con neb.x

**Qué responde.** ¿Cuál es el camino de mínima energía entre reactivo y producto y cuánto vale la barrera de activación (directa e inversa)?

**Fundamento para no expertos.** Entre dos mínimos de energía hay un "puerto de montaña": el estado de transición. La banda elástica (NEB) tiende una cadena de imágenes entre reactivo y producto, unidas por muelles, y relaja cada imagen perpendicularmente al camino hasta que la cadena descansa en el valle. La imagen trepadora (CI) empuja la imagen más alta hasta el puerto exacto; sin ella la barrera sale subestimada.

**Fórmulas.** En `qekit/modules/neb.py`, `neb.collect`:
$$E_a^{\rightarrow} = E_{\max} - E_1,\qquad E_a^{\leftarrow} = E_{\max} - E_N,\qquad \Delta E = E_N - E_1$$
- Energías en eV relativas a la primera imagen (columna 2 de `<prefix>.dat`); si `neb.out` trae `activation energy (->)`/`(<-)`, se usan esos. Conversión a kJ/mol: × 96.485.

**Cómo lo calcula Olla-DFT.**
1. `neb.comprobar_extremos`: mismo número y **orden** de átomos, misma celda (tol 1e-4), estructuras no idénticas; si falla, aborta.
2. `neb.build_neb_input` escribe `neb.in`: `&PATH` con `string_method='neb'`, `nstep_path=--nstep` (50), `ds=1`, `opt_scheme='broyden'`, `num_of_images=--images` (7), `k_max=0.3`, `k_min=0.2`, `CI_scheme='auto'` (o `'no-CI'` con `--no-ci`), `path_thr=--path-thr` (0.05 eV/Å); motor `pw.x` recortado de `inputgen.build_pw_input` (sin posiciones ni celda); `FIRST_IMAGE`/`LAST_IMAGE` en Å con `0 0 0` en los átomos `--fix`; `CELL_PARAMETERS`.
3. El usuario corre `neb.x -inp neb.in > neb.out`.
4. `neb.collect --collect`: lee `<prefix>.dat` (s, E, F), `<prefix>.int` (interpolación), y de `*.out` barreras, convergencia (`convergence achieved`), iteraciones, `CI_scheme` y las imágenes con *"scf convergence NOT achieved on image"*.
5. `report`: barreras, tabla por imagen, aviso si el máximo interpolado cae a más de 0.4 pasos de cualquier imagen; `export` (`NEB.dat`, `NEB.txt`); `plot`.

**De dónde sale cada dato.**

| Dato | Origen | Detalle |
|---|---|---|
| $s$, $E$, $F$ por imagen | `<prefix>.dat` de `neb.x` | `neb.collect` |
| Curva interpolada | `<prefix>.int` de `neb.x` | opcional |
| Barreras, convergencia, iteraciones, CI | `neb.out` (regex) | prioridad sobre el cálculo propio |
| eV → kJ/mol | 96.485 en `neb.report` | — |

**Límites y trampas.**
- *"Esta barrera es ELECTRÓNICA, a 0 K y sin energía de punto cero."* Correcciones térmicas en `thermochem`.
- Sin CI: *"esta barrera es una COTA INFERIOR."* Pocas imágenes (< 5): aviso.
- Imágenes con scf no convergido: *"El scf NO convergió en la(s) imagen(es)…: por eso el perfil sale dentado."*
- Los extremos deben estar relajados con los mismos parámetros; el módulo no lo comprueba.

**Referencias.**
- G. Henkelman, B. P. Uberuaga, H. Jónsson, *J. Chem. Phys.* 113, 9901 (2000) — climbing-image NEB. DOI: 10.1063/1.1329672.
- G. Henkelman, H. Jónsson, *J. Chem. Phys.* 113, 9978 (2000) — tangente mejorada.

---

### `olla-dft amorphous` — Sólido amorfo por fundido y temple con MLIP

**Qué responde.** ¿Cómo generar una estructura amorfa de composición y densidad dadas, y qué coordinación y distancias de primer vecino tiene?

**Fundamento para no expertos.** Un vidrio no se dibuja: se fabrica calentando el material hasta fundirlo y enfriándolo tan rápido que no le da tiempo a cristalizar. En el ordenador el temple es millones de veces más rápido que en el laboratorio, así que el resultado es más desordenado y algo menos denso que el real. Aquí la dinámica se hace con un potencial interatómico aprendido (MACE por omisión), no con DFT, porque hacen falta miles de pasos; la estructura resultante es un punto de partida que luego debe relajarse con `pw.x`.

**Fórmulas.** En `qekit/modules/amorphous.py`.

Arista de la celda cúbica (`celda_para_densidad`) y densidad (`densidad_de`):
$$L = \left(\frac{\sum_i m_i\,u}{\rho}\right)^{1/3}\times 10^{8},\qquad \rho = \frac{\sum_i m_i\,u}{V}$$
- $m_i$ en uma; $u = 1.66053906660\times10^{-24}$ g; $\rho$ en g/cm³; $V$ en Å³ (× $10^{-24}$ cm³).

Velocidad de temple (`Protocolo.velocidad_temple`):
$$\dot T = \frac{T_{\mathrm{fundido}} - T_{\mathrm{final}}}{N_{\mathrm{temple}}\,\Delta t}$$
- Por omisión $(3000 - 300)\,\mathrm{K}/(1000\times 1\ \mathrm{fs}) = 2.7\times10^{15}$ K/s.

Coordinación (`coordinaciones`): $Z_{ab} = \frac{1}{N_a}\sum_{i\in a}\#\{j\in b: d_{ij} < 1.25\,(r_a^{\mathrm{cov}} + r_b^{\mathrm{cov}})\}$ con imagen mínima; distancia media de primer vecino con el mismo corte (`distancia_media`).

**Cómo lo calcula Olla-DFT.**
1. `formula_a_simbolos` expande `SiO2` × `--units` (8).
2. `empaquetar` coloca átomos al azar (semilla `--seed`) rechazando distancias < `--min-dist` × (suma de radios covalentes), `FACTOR_MINIMO` = 0.75; hasta 20000 intentos por átomo; error si no caben.
3. `fundir_y_templar` (salvo `--pack-only`): calculador `mlip.calculator(--model)`; velocidades de Maxwell–Boltzmann a `--melt` (3000 K); `Langevin` de ASE con `friction=0.02` y `--dt` (1 fs); `--melt-steps` (500) a $T_{\mathrm{fundido}}$; temple en 20 tramos de $N_{\mathrm{temple}}/20$ pasos bajando la temperatura del termostato linealmente hasta `--final` (300 K); `--anneal-steps` (200) a $T_{\mathrm{final}}$. Se registra $E$ y $T$ cada 10 pasos (`traza.dat`).
4. Avisos: temperatura final $> 2.5\,T_{\mathrm{final}} + 200$ K (el termostato no siguió la rampa) y $\dot T > 10^{13}$ K/s.
5. `report` (densidad, protocolo, coordinaciones, distancias, $T$ final) y `export` (`amorfo.cif`, `AMORFO.dat`, `AMORFO.txt`).

**De dónde sale cada dato.**

| Dato | Origen | Detalle |
|---|---|---|
| Masas y radios covalentes | `ase.data.atomic_masses`, `covalent_radii` | — |
| uma → g | constante local $1.66053906660\times10^{-24}$ | CODATA 2018 |
| Energías y fuerzas | potencial MLIP (`mlip.calculator`) | MACE-MP-0 small, CHGNet o M3GNet |
| Densidad objetivo | `--density` | g/cm³ |
| Protocolo | `--melt`, `--final`, `--melt-steps`, `--quench-steps`, `--anneal-steps`, `--dt` | K, pasos, fs |

**Límites y trampas.**
- *"Esta estructura viene de un potencial aprendido, NO de DFT… relájala con 'olla-dft gen -p relax'… y compara varias realizaciones (--seed distintas)."*
- El protocolo por omisión es de **exploración**: 2.7×10¹⁵ K/s, y el reporte lo avisa (*"Velocidad de temple X K/s. Un vidrio de verdad se enfría a 1-100 K/s"*). El docstring y la ayuda de `--quench-steps` lo dicen: 27 000 pasos bajan a 10¹⁴ K/s, diez veces más a 10¹³ K/s, que es donde desaparece el aviso.
- Dinámica NVT a volumen fijo: la densidad final es la impuesta, no se relaja.
- Con `friction=0.02` y rampas rápidas el sistema puede acabar líquido: *"El sistema acabó a X K, no a los Y K pedidos."*
- Requiere `torch` + el paquete del modelo (no son dependencias de Olla-DFT).

**Referencias.**
- I. Batatia et al., *MACE-MP-0* (arXiv:2401.00096, 2023).
- ASE Langevin: A. H. Larsen et al., *J. Phys.: Condens. Matter* 29, 273002 (2017).

---

### `olla-dft mlip` — Pre-relajación, barrido de volumen y cribado de fonones con un potencial aprendido

**Qué responde.** Antes de gastar DFT: ¿cuál es una geometría casi relajada, dónde está aproximadamente el mínimo $E(V)$ y tiene la estructura frecuencias imaginarias?

**Fundamento para no expertos.** Un potencial interatómico aprendido (MLIP) da energías y fuerzas miles de veces más barato que DFT. No sustituye a `pw.x` —está entrenado con datos PBE de Materials Project y describe *otra* superficie de energía— pero sirve para llegar al cálculo DFT con la geometría casi lista, para acotar el rango de una ecuación de estado y para detectar antes de la DFPT que una estructura no está en un mínimo.

**Fórmulas.** En `qekit/modules/mlip.py`.

Relajación (`mlip.relax`): BFGS de ASE hasta $f_{\max} <$ `--fmax` (0.01 eV/Å) o `--steps` (300), con `FrechetCellFilter` si se relaja la celda. Presión:
$$P = -\tfrac{1}{3}\,\mathrm{tr}\,\boldsymbol\sigma\times 160.21766208\ \ [\mathrm{GPa}]$$

Barrido de volumen (`mlip.volume_scan`): 15 escalas en $[1-s, 1+s]$, $s$ = `--span` (0.10); parábola $E = aV^2 + bV + c$:
$$V_0 = -\frac{b}{2a},\qquad B_0 \approx 2aV_0\times160.21766208\ \mathrm{GPa},\qquad \text{escala} = (V_0/V)^{1/3}$$

Fonones por diferencias finitas (`mlip.phonon_check`, `frequencies`): hessiano $H_{i\alpha,j\beta} = -\partial F_{j\beta}/\partial u_{i\alpha}$ centrado con $\delta$ = 0.01 Å en una supercelda `--supercell` (2×2×2), simetrizado; matriz dinámica $D = H/\sqrt{m_im_j}$;
$$\omega = \mathrm{sign}(\lambda)\sqrt{|\lambda|}\times 521.4708\ \mathrm{cm^{-1}}$$
- $\lambda$: autovalores de $D$ en eV/(Å²·uma); imaginarias si $\omega < -5$ cm⁻¹.

**Cómo lo calcula Olla-DFT.**
1. `mlip.calculator` carga MACE (`mace_mp(model=--size, default_dtype='float64')`), CHGNet o M3GNet; si falta el paquete explica qué instalar.
2. `relax`: fuerzas y presión inicial/final, desplazamiento máximo, cambio de volumen; avisos si no converge o si algún átomo se movió > 0.5 Å. Escribe la estructura (`relajado_mlip.cif`) y `MLIP_PROCEDENCIA.json` (`write_provenance`) para que `audit` sepa que no es DFT.
3. `scan`: `report_scan` sugiere `olla-dft eos --scale X --span 0.04`; avisa si el mínimo cae fuera del rango.
4. `phonons`: `report_phonon`; código de salida 1 si hay imaginarias.

**De dónde sale cada dato.**

| Dato | Origen | Detalle |
|---|---|---|
| $E$, $F$, $\sigma$ | calculador MLIP | `mace_mp`, `CHGNetCalculator`, `PESCalculator` |
| eV/Å³ → GPa | 160.21766208 | constante local |
| $\sqrt{\mathrm{eV/(Å^2\,uma)}}$ → cm⁻¹ | 521.4708 | `CONV` |
| Masas | `atoms.get_masses()` (ASE) | uma |

**Límites y trampas.**
- *"ESTO NO ES EL RESULTADO FINAL. El modelo está entrenado con datos PBE… no mezcles sus energías con las de QE."* Ejemplo del reporte: Si, MACE 5.464 Å vs LDA 5.402 Å.
- `phonon_check` diagonaliza la matriz dinámica **completa** de la supercelda: salen los modos de Γ de la primitiva y además los de los puntos q que la supercelda pliega sobre Γ (así lo declara el docstring). No es una dispersión.
- El $B_0$ del barrido es de una parábola: *"sirve para saber el orden de magnitud, no para reportarlo."*
- Sin `torch`/`mace-torch`: *"para usar 'mace' hace falta instalar 'mace-torch'… Ocupa algo más de 1 GB."*

**Referencias.**
- I. Batatia, D. P. Kovács, G. N. C. Simm, C. Ortner, G. Csányi, *NeurIPS* 35 (2022) — MACE.
- B. Deng et al., *Nat. Mach. Intell.* 5, 1031 (2023) — CHGNet.
- C. Chen, S. P. Ong, *Nat. Comput. Sci.* 2, 718 (2022) — M3GNet.

---

### `olla-dft audit` y `olla-dft db` — Comparabilidad entre cálculos e índice local

**Qué responde.** ¿Se pueden restar las energías totales de este conjunto de cálculos? Y `db`: ¿qué cálculos tengo, con qué parámetros y qué salió?

**Fundamento para no expertos.** Dos energías totales de QE sólo se pueden restar si vienen de la misma "receta": mismo funcional, mismos pseudopotenciales, mismos cutoffs y mismo tratamiento de las ocupaciones. Si no, la diferencia es un número perfectamente formado sin significado, y QE no avisa. La auditoría calcula una huella con esos parámetros y agrupa: más de un grupo = no comparables. La malla k se trata aparte como aviso, comparando la **densidad** de puntos k, que es lo comparable entre celdas distintas.

**Fórmulas.** `audit.kdensity`:
$$\rho_k = \frac{n_1 n_2 n_3}{(2\pi)^3 / V}\quad[\text{puntos}/\text{Å}^{-3}]$$

Huella (`qeout.QEResult.fingerprint` + `origen`): (origen, funcional, {elemento: UPF}, `ecutwfc`, `ecutrho`, `smearing`, `degauss`, `occupations`, `nspin`).

**Reglas implementadas.**

| Regla | Dónde | Efecto |
|---|---|---|
| Origen DFT vs MLIP entra en la huella | `audit.audit` (lee `MLIP_PROCEDENCIA.json` vía `mlip.read_provenance`) | grupos distintos: NO COMPARABLES |
| Funcional, pseudos, ecutwfc, ecutrho, smearing, degauss, ocupaciones, nspin | `_campos`/`ETIQUETAS` | se listan los que difieren |
| SCF no convergido | sólo `scf/relax/vc-relax/md/vc-md` con `converged=False` | "NO CONVERGIERON — sus energías no sirven" |
| `nscf`/`bands` | por tipo de cálculo | "Sin energía utilizable" |
| Densidad de k dispar | $\max\rho_k/\min\rho_k > 2$ | AVISO, no incompatibilidad |
| Carpeta sin XML propio pero con hijas | `audit.collect` | se auditan las hijas (barrido) |

**Cómo lo calcula Olla-DFT.**
1. `audit.collect(paths)`: para cada carpeta lee la marca MLIP y el XML (`qeout.read_xml`).
2. `audit.audit`: agrupa por huella, lista diferencias, no convergidos y sin energía.
3. `audit.report`; código de salida 1 si no son comparables. `--index` registra en `olla-dft.db`.
4. `db carpeta/…` indexa (`audit.index`, `INSERT OR REPLACE` por ruta absoluta); `db --query "SELECT …"` (sólo SELECT); `db --formula/--calculation/--gap-min/--gap-max` (`audit.search`); `db --export` (JSON); sin argumentos, `audit.summary`.

**De dónde sale cada dato.** Todo del XML de `pw.x` (`qeout.read_xml`): funcional, `pseudo_files`, cutoffs, smearing, ocupaciones, `nspin`, energía (Ha → eV), volumen, presión, fuerza máxima, `homo/lumo` → gap, magnetización, convergencia, pasos SCF, `nk` (puntos k usados), `nbnd`, pasos BFGS, tiempo de reloj; más `MLIP_PROCEDENCIA.json` si existe. Columnas de la tabla `calculos` en `audit.ESQUEMA`.

**Límites y trampas.**
- La huella no incluye la malla k ni la celda: *"un bulk y una losa necesitan mallas distintas por construcción."*
- Compara los **nombres** de los UPF, no su contenido: dos archivos distintos con el mismo nombre pasan.
- `hull` y `thermo.from_runs` se apoyan en esta auditoría y se niegan a mezclar orígenes.
- `db --query` sólo admite `SELECT`; bases antiguas se migran añadiendo `nk`, `nbnd`, `n_bfgs` (`_migrar`).

**Referencias.** Manual de Quantum ESPRESSO (esquema XML `qes`); K. Lejaeghere et al., *Science* 351, aad3000 (2016) — por qué los pseudos y cutoffs fijan la referencia de energía.

---

### `olla-dft hull` — Energías de formación y casco convexo

**Qué responde.** ¿Es estable cada fase frente a descomponerse en las demás, y cuánta energía por átomo está por encima del casco convexo?

**Fundamento para no expertos.** Se dibuja la energía de formación por átomo contra la composición. La curva más baja que envuelve todos los puntos desde abajo (casco convexo) une las fases estables; cualquier fase por encima gana energía descomponiéndose en las dos (o tres) del casco que la rodean, y esa distancia vertical es $E_{\mathrm{hull}}$. Es energía a 0 K sin entropía: una fase 25 meV/átomo por encima a veces se sintetiza igual.

**Fórmulas.** En `qekit/modules/thermo.py`.
$$E_f = \frac{E(\text{compuesto}) - \sum_i n_i\,\mu_i}{N},\qquad \mu_i = \min_{\text{fases puras de } i}\frac{E}{N}$$
$$E_{\mathrm{hull}} = E_f - E_{\mathrm{casco}}(\mathbf x)$$
- Binario (`_casco`): envolvente inferior por cadena monótona sobre $x$ y interpolación lineal. Ternario o más: `scipy.spatial.ConvexHull` en $(x_1,\dots,x_{n-1}, E_f)$, se conservan las facetas con normal hacia abajo en la energía (`eq[-2] < 0`), y $E_{\mathrm{casco}}$ se obtiene por coordenadas baricéntricas dentro de la faceta (`Delaunay.find_simplex`).

**Cómo lo calcula Olla-DFT.**
1. `audit.collect` + `audit.audit`; si no son comparables, imprime la auditoría y se niega salvo `--force`.
2. `thermo.from_runs`: descarta `nscf/bands`, sin energía o no convergidos; rechaza mezclar DFT y MLIP; fórmula con `ase.Atoms`; referencias elementales = mínima energía por átomo de las fases puras (aviso si falta alguna).
3. `_casco`; `report` con umbral de metaestabilidad `--threshold` (0.025 eV/átomo): ESTABLE / metaestable / inestable / fuera del dominio.
4. `export` (`CASCO_CONVEXO.dat`); `plot` sólo para binarios.

**De dónde sale cada dato.** Energías totales y símbolos del XML de cada carpeta (`qeout.read_xml`); orden de elementos de `--elements` o alfabético.

**Límites y trampas.**
- *"Esto es energía a 0 K, sin punto cero ni entropía."*
- Sin referencias elementales: *"hay que calcular cada elemento puro en su fase estable, con los mismos parámetros."*
- `--force` construye el casco con cálculos no comparables bajo la responsabilidad del usuario.
- Un elemento puro con varias fases: la más baja es la referencia; las otras salen con $E_f > 0$.
- La gráfica sólo para binarios.

**Referencias.** S. P. Ong, L. Wang, B. Kang, G. Ceder, *Chem. Mater.* 20, 1798 (2008); W. Sun et al., *Sci. Adv.* 2, e1600225 (2016) — escala de metaestabilidad.

---

### `olla-dft doctor` — Diagnóstico de convergencia de pw.x

**Qué responde.** ¿Sirve este cálculo y, si el SCF no convergió, es por oscilación de carga (mezclar menos) o por lentitud (mezclar más o más pasos)?

**Fundamento para no expertos.** El ciclo autoconsistente mezcla la densidad nueva con la vieja. Si mezcla demasiado, la carga "chapotea" de un lado a otro de la celda (oscilación, típica en losas y metales) y el error sube y baja; si mezcla poco, el error baja siempre pero despacio. Los dos remedios son opuestos, así que el módulo mira la **forma** de la curva de `estimated scf accuracy`.

**Reglas implementadas** (`diagnose._clasificar`, sólo si no convergió):

| Condición | Diagnóstico | Consejo |
|---|---|---|
| < 8 iteraciones | `pocos_datos` | subir `electron_maxstep` a ≥ 100 |
| (≥ 6 puntos y > 25 % de subidas tras las 2 primeras iteraciones) **o** el error se multiplica > 5× en una iteración | `oscilacion` | `mixing_beta = max(0.05, β/3)`, `mixing_mode='local-TF'`, `mixing_ndim=12` |
| bajó < 3 órdenes de magnitud en total | `estancada` | revisar `starting_magnetization`, smearing, distancias |
| resto, con β ≥ 0.6 | `lenta` | `electron_maxstep = 300` (no subir β) |
| resto, con β < 0.6 | `lenta` | `mixing_beta = min(0.7, max(1.75β, 0.3))`, `electron_maxstep = 300` |

Problemas del XML (`diagnose.diagnose`): SCF no convergido; fuerza residual > 0.05 eV/Å; $|P| > 1$ GPa en `scf/relax/vc-relax`; `Error in routine`. Relajación: aviso si la energía subió en más de $N/3$ pasos.

**Cómo lo calcula Olla-DFT.**
1. `qeout.find_xml` + `read_xml` (convergencia, pasos, error, fuerzas, presión, magnetización, tiempos).
2. `diagnose.find_stdout` busca el archivo con `Program PWSCF`; `read_scf_history` parte el stdout en ciclos SCF con `_ciclos_scf` (cada `iteration #  1` abre uno; en un `relax` hay uno por paso iónico), guarda `n_ciclos` y extrae **sólo del último ciclo** `estimated scf accuracy`, `total energy`, `beta`, `convergence has been achieved` / `convergence NOT achieved`; `read_trajectory` lee los `!    total energy`, `Total force`, `P=` de todo el archivo.
3. `report` y `plot` (precisión SCF en log y energía por paso iónico). Código 1 si hay problemas. `--system` delega en `health.check` (instalación).

**De dónde sale cada dato.** XML (`converged`, `n_scf_steps`, `scf_error`, `max_force` en eV/Å, `pressure` en GPa, `wall_time`) y stdout de `pw.x` (regex `_RE_ACC`, `_RE_ETOT`, `_RE_ITER`, `_RE_FORCE`, `_RE_PRESS`, `_RE_WARN`, `_RE_MAXSTEP`). $\beta$ por omisión 0.4 si no aparece.

**Límites y trampas.** En un `relax` sólo se diagnostica el último ciclo SCF (el reporte lo dice: *"en el último de N ciclos SCF (uno por paso iónico; se diagnostica solo el último)"*); un ciclo intermedio que oscilara no se ve. Los umbrales (0.05 eV/Å, 1 GPa) son fijos. No detecta problemas de simetría ni de pseudos.

**Referencias.** D. D. Johnson, *Phys. Rev. B* 38, 12807 (1988) — mezcla de Broyden; G. Kresse, J. Furthmüller, *Phys. Rev. B* 54, 11169 (1996) — oscilación de carga y `local-TF`.

---

### `olla-dft crosscheck` — La misma cantidad por dos caminos independientes

**Qué responde.** ¿Coinciden dos rutas físicamente independientes hacia la misma magnitud? Si no, algo está mal en una de ellas.

**Fundamento para no expertos.** Comparar contra la literatura detecta errores en un módulo, pero no un sesgo sistemático compartido. Calcular $B_0$ de la ecuación de estado y de las constantes elásticas, o el gap de bandas y el de Tauc, son rutas que no comparten código: si coinciden, es difícil que las dos estén mal igual.

**Cruces implementados** (`crosscheck.run`; desviación relativa $|b-a|/|a|$, o absoluta si $a = 0$):

| # | Magnitud | Ruta A | Ruta B | Tolerancia | Datos |
|---|---|---|---|---|---|
| 1 | $B_0$ | `EOS.txt` (línea con `B0` y `GPa`) | $B_{\mathrm{Hill}}$ de `ELASTIC_C.dat` (`elastic.moduli`) | 5 % | ambos archivos |
| 2 | $v_L[100]$, $v_T[100]$ | $\sqrt{C_{11}/\rho}$, $\sqrt{C_{44}/\rho}$ (`derived.cubic_directional`) | pendiente LA/TA en Γ de `FONONES_BANDAS.dat` | 10 % | Cij, bandas, masas, volumen |
| 3 | $\Theta_D$ | velocidades del sonido (`derived.debye_from_velocity`) | segundo momento de `FONONES_DOS.dat` | 30 % (definiciones distintas) | Cij, DOS, N |
| 4 | gap óptico | `--gap-bandas` | `--gap-tauc` | 6 % | parámetros |
| 5 | $C_v$ a 1500 K | $3Nk_B$ (Dulong–Petit) | $k_B\int x^2 e^x/(e^x-1)^2\,g(\omega)\,d\omega$ con $g$ normalizada a $3N$ (`_cv_alta_T`) | 3 % | DOS |
| 6 | número de modos | $3N$ | $\int g(\omega)\,d\omega$ | 5 % | DOS |
| 7 | $\kappa_L$ | `KAPPA.dat` a ~300 K | modelo de Slack desde Cij (`derived.slack`) | 60 % | KAPPA, Cij |
| 8 | fase de Berry | `BERRY.dat` (columna 3 en carga 0) | $-2\sum_n (\bar r_n\cdot b)/2\pi$ de `WANNIER_centros.dat`, misma rama mód. 2 | 0.05 | ambos, celda |
| 9 | función trabajo | `ESM.dat` (Φ en $q = 0$) | `WF.dat` (`Phi_eV`) | 5 % | ambos |
| 10 | $B_0$ (tercera vía) | `EOS.txt` | $-\tfrac{1}{3}\,dP/d\epsilon$ de `STRAIN.dat` (kbar → GPa × 0.1) | 10 % | STRAIN (hidrostático) |

Constantes: `KB_EV` = $8.617333262\times10^{-5}$ eV/K; cm⁻¹ → eV: $1.239841984\times10^{-4}$.

**Cómo lo calcula Olla-DFT.** `crosscheck._cargar` busca recursivamente los archivos de resultados en la carpeta del proyecto; con `-f estructura` toma masas, volumen, N y celda; `run` ejecuta cada cruce para el que haya datos; `report` marca OK/FALLA con el diagnóstico de qué mirar primero. Código 1 si alguno falla.

**Límites y trampas.** *"Un cruce que falla NO dice cuál de los dos caminos está mal."* El cruce 3 compara definiciones distintas de $\Theta_D$ (*"coincidir al 1 % sería sospechoso"*); el 10 sólo vale si el barrido fue hidrostático; el 2 acusa antes a la malla q que a las Cij. Los cruces 8–10 tragan cualquier excepción en silencio (`except Exception: pass`).

**Referencias.** R. Hill, *Proc. Phys. Soc. A* 65, 349 (1952); G. A. Slack, *Solid State Phys.* 34, 1 (1979); R. D. King-Smith, D. Vanderbilt, *Phys. Rev. B* 47, 1651 (1993).

---

### `olla-dft selftest` — Validación contra la física conocida

**Qué responde.** ¿Reproduce Olla-DFT valores medidos, publicados o exactos, y no sólo lo que él mismo dice?

**Fundamento para no expertos.** Las pruebas unitarias comparan el código consigo mismo. Aquí cada prueba calcula una magnitud con una respuesta conocida (constantes de Ewald, entropía de Sackur–Tetrode, $T_c$ del aluminio, invariantes topológicos…) y la contrasta con esa referencia y su fuente. `--quick` (por omisión) corre las que no necesitan `pw.x`; `--full` añade las que sí; `--mlip` la que necesita MACE.

**Pruebas y referencias** (`selftest.PRUEBAS`; desviación relativa, o absoluta si la referencia es 0):

| Clave | Magnitud | Referencia | Tol. | Fuente (según el código) | Función probada |
|---|---|---|---|---|---|
| `madelung` | $\alpha_M$ cúbica simple | 2.8372974 | 1e-5 | valor clásico de Ewald | `defects.constante_madelung` |
| `lorenz` | $L/L_0$ gas de electrones libres | 1.0 | 12 % | límite de Sommerfeld | `transport.compute`, `lorenz` |
| `npw` | ondas planas de Si a 30 Ry | 725 | 6 % | lo que reporta `pw.x` (V = 39.5 Å³) | `cost.n_ondas_planas` |
| `sackur` | $S_{\mathrm{trans}}$ de N₂ a 298 K | 150.4 J/(mol·K) | 1 % | Sackur–Tetrode, NIST-JANAF | `thermochem.S_traslacional` |
| `allen_dynes` | $T_c$ del Al (λ=0.44, ω_log=270 K, µ*=0.12) | 1.18 K | 12 % | Allen–Dynes 1975, exp. | `elph.allen_dynes` |
| `allen_dynes_mu` | $T_c(0.10)/T_c(0.12)$ | 1.56 | 5 % | exponencialidad en µ* | `elph.allen_dynes` |
| `born2d` | $Y_{2D}$ con C11=352, C12=60 N/m | 341.8 N/m | 1 % | $Y = C_{11} - C_{12}^2/C_{11}$ (grafeno DFT) | `elastic.modulos_2d` |
| `gap_invariante` | ΔE_v de un material consigo mismo | 0 eV | 1e-9 | identidad exacta | `align.alinear` |
| `ewald_escala` | $\lvert\alpha(3) - \alpha(30)\rvert$ | 0 | 1e-6 | invariancia de escala | `defects.constante_madelung` |
| `chern_qwz` | $C$ del modelo de Qi–Wu–Zhang (m=−1) | −1 | 1e-10 | PRB 74, 085308 (2006) | `topology.invariants_from_vectors` |
| `umklapp` (`--mlip`) | exponente $n$ en $\kappa\propto T^{-n}$ del Si | 1.0 | 25 % | ley de Umklapp por encima de $\Theta_D$ | `kappa.*` con MACE |
| `her_pt` | $\Delta G_{\mathrm{H^*}}$ con $E_{\mathrm{ads}} = -0.33$ eV | −0.09 eV | 5 % | Nørskov 2005, Pt(111) | `echem.her` |
| `oer_ruo2` | η con ΔG(OH,O,OOH)=(0.77, 2.16, 3.87) | 0.48 V | 10 % | Man et al. 2011 | `echem.oer` |
| `escala_oer` | ΔG(OOH) − ΔG(OH) del perfil de RuO₂ | 3.2 eV | 10 % | relación universal de escala | `echem.oer` + `echem.escala_ooh_oh` |
| `escala_eta_min` | $\eta_{\min}$ = Δ/2 − ΔG_total/4 | 0.37 V | 2 % | Man et al. 2011 | `echem.sobrepotencial_minimo_escala` |
| `fonon_si` (`--full`) | ω(Γ) óptico del Si | 520 cm⁻¹ | 10 % | Raman exp. 520.7 cm⁻¹ | `phonons.*` con `ph.x` |
| `wannier_si` (`--full`) | centro de Wannier Si–Si | 1.17563 Å | 2 % | $\sqrt3\,a/8$ con a = 5.43 Å | `wannier.*` |
| `condensador` (`--full`) | pendiente de $1/C$ vs $d$ / $(1/\varepsilon_0)$ | 1.0 | 6 % | electrostática del condensador plano | `esm.*` bc3 Al(111) |
| `born_si` (`--full`) | $Z^*$ del Si | 0 e | 0.05 | regla de suma acústica | `berry.*` |
| `gamma_al` (`--full`) | γ de Al(111) | 1.10 J/m² | 25 % | Vitos 1998 (1.20), exp. 1.14 | `surfen.*` |
| `bulk_si` (`--full`) | $B$ del Si por deformación | 95 GPa | 15 % | LDA 93–97 (Nielsen & Martin 1985), exp. 98 | `strain.*` |
| `sitio_h_al` (`--full`) | $E_{\mathrm{ads}}$(top) − $E_{\mathrm{ads}}$(hollow), H/Al(111) | 5.6 eV | 60 % | orden hueco < puente < top | `adsorb.*` |

**Cómo lo calcula Olla-DFT.** `selftest.ejecutar` filtra por `--only`, `--full`, `--mlip`; crea una carpeta temporal (`--keep` para conservarla); ejecuta cada `fn(ctx)` y mide el tiempo; `report` lista valor, referencia, desviación, tolerancia y fuente. Código 1 si alguna falla o da error. `--list` imprime la tabla sin correr nada.

**Límites y trampas.** *"Las que salen MAL no siempre son un fallo del código: una tolerancia ajustada, un pseudopotencial distinto o un cutoff bajo también las mueven."* Las de `--full` dependen de los pseudos de `--pseudo-dir` y de que `pw.x`/`ph.x` funcionen.

**Nota sobre `qekit/modules/uncertainty.py`.** No tiene comando propio. Ofrece `propagate(f, valores, sigmas)` — propagación en cuadratura con derivadas centradas, $\sigma_f^2 = \sum_i (\partial f/\partial x_i)^2\sigma_i^2$, paso relativo $10^{-6}$, entradas independientes — y `weighted_mean` — media ponderada con $w_i = 1/\sigma_i^2$ y $\sigma = (\sum w_i)^{-1/2}$. Ningún módulo de esta parte la invoca; sólo `validation`/`results` comprueban que las incertidumbres declaradas sean finitas y no negativas.

**Referencias.** P. B. Allen, R. C. Dynes, *Phys. Rev. B* 12, 905 (1975); X.-L. Qi, Y.-S. Wu, S.-C. Zhang, *Phys. Rev. B* 74, 085308 (2006); L. Vitos et al., *Surf. Sci.* 411, 186 (1998); O. H. Nielsen, R. M. Martin, *Phys. Rev. B* 32, 3792 (1985).

---

### `olla-dft suggest` — Parámetros a partir del historial propio

**Qué responde.** Según los cálculos que ya convergieron con estos elementos, ¿qué `ecutwfc`, dual, densidad de k y `electron_maxstep` conviene usar?

**Fundamento para no expertos.** Con unas decenas de cálculos no tiene sentido entrenar nada: se buscan los cálculos parecidos (comparten elementos, tamaño similar) y se mira qué les funcionó, diciendo cuántos casos respaldan cada número.

**Reglas implementadas** (`recommend.similares`, `recommend.sugerir`):

| Regla | Detalle |
|---|---|
| Similitud | sólo cálculos con `convergido`; puntaje = Jaccard de elementos $\lvert A\cap B\rvert/\lvert A\cup B\rvert$; × 0.5 si $N_{\mathrm{at}}$ difiere en más de un factor 2 |
| `ecutwfc` | **máximo** entre los similares (no la media), con rango |
| dual | máximo de `ecutrho/ecutwfc` |
| densidad de k | mediana de `kdensity` (puntos/Å⁻³) |
| `electron_maxstep = 300` | si la mediana de `n_scf` > 40 |
| `mixing_beta = 0.3` + `local-TF` | si la estructura es una losa (vacío en $c$ > 8 Å), regla general, 0 casos |
| Confianza | alta ≥ 8 casos, media ≥ 3, baja < 3 |
| Sin historial | remite a los cutoffs del propio UPF / SSSP |

**Cómo lo calcula Olla-DFT.** `_cmd_suggest` carga la estructura, lee `SELECT * FROM calculos` de `--db` (`olla-dft.db`), detecta si es losa y llama a `recommend.sugerir`; `report` imprime valor, número de casos y razón.

**De dónde sale cada dato.** Tabla `calculos` de `olla-dft.db` (`audit.index`): `formula`, `natoms`, `ecutwfc`, `ecutrho`, `kdensity`, `n_scf`, `convergido`.

**Límites y trampas.** *"No sustituyen a una prueba de convergencia: 'olla-dft converge' sigue siendo la forma de saberlo de verdad."* Con confianza "baja" el reporte marca *"UN SOLO CASO: tómalo como indicio"* si hay 1 caso y *"SOLO n CASOS"* si hay 2. No inventa cutoffs sin historial.

**Referencias.** G. Prandini, A. Marrazzo, I. E. Castelli, N. Mounet, N. Marzari, *npj Comput. Mater.* 4, 72 (2018) — SSSP.

---

### `olla-dft pseudos` — Elegir pseudopotenciales con criterio

**Qué responde.** De los UPF disponibles para cada elemento, ¿cuáles sirven para la tarea (óptica, espín-órbita, XANES, DFT+U, fonones) y cuál conviene?

**Fundamento para no expertos.** En una carpeta suele haber varios pseudopotenciales por elemento, de familias y funcionales distintos. Elegir el primero por orden alfabético falla en silencio: un pseudo escalar-relativista con `lspinorb` da un desdoblamiento de cero, un ultrasuave con `epsilon.x` da un espectro entero equivocado, y mezclar funcionales entre elementos invalida la energía total. El selector aplica requisitos duros (que descartan) y preferencias (que ordenan) y explica cada decisión.

**Reglas implementadas** (`pseudos.TAREAS`, `pseudos.evaluar`):

| Tarea | Requisito duro | Preferencia |
|---|---|---|
| `optics` | tipo ∈ {NC} | — |
| `soc` | relativista = `full` (salvo elementos con Z < 19: aviso y −0.5 puntos) | — |
| `xanes` | UPF con secciones `PP_GIPAW` | — |
| `hubbard` | — | +0.15 × `z_valence` (semicore) |
| `fonones` | — | +2.0 si tipo ∈ {NC, US} |
| `general` | — | — |
| todas | funcional igual al de `--functional` (alias PBE/`SLA PW PBX PBC`, PZ/LDA, PBEsol, BLYP, revPBE) | +max(0, (90 − ecutwfc)/30); −0.5 sin cutoff declarado; +1.0 US/PAW con `--cheap`; +0.3 con GIPAW; +0.2 si `full` |

Orden final: no descartados primero, luego puntos descendentes, luego nombre. Coherencia entre elementos (`pseudos.coherencia`): aviso si mezclan funcionales, si mezclan NC con US/PAW (manda el dual del ultrasuave) y si los cutoffs sugeridos difieren en más de 2.5×.

**Cómo lo calcula Olla-DFT.**
1. `pseudos.candidatos`: `pseudo.find_for_element` (archivos `.UPF` cuyo nombre empieza por el símbolo) y `pseudos.leer` (tipo, funcional normalizado por `_funcional`/`NOMBRE_CORTO`, relativista, `z_valence`, cutoffs sugeridos, GIPAW, tamaño).
2. `pseudos.evaluar` y `elegir`; `report` con tabla y descartados; `report_coherencia` si hay más de un elemento; imprime la línea `--pseudo EL=archivo` para reutilizar.
3. El mismo selector lo usa `sweep.prepare_common` en todos los comandos (`pseudo.resolve` → `_elegir` con la tarea) y `_coherencia_de_funcional` reelige para unificar funcional (preferencia PBE > PBEsol > revPBE > PZ > BLYP).

**De dónde sale cada dato.** Cabecera del UPF (primeros 20–30 kB): `pseudo_type`, `functional`, `relativistic`, `z_valence`, `wfc_cutoff`/`rho_cutoff` (o sus equivalentes v1), presencia de `PP_GIPAW`. `Z_SOC` = 19 en `pseudos.py`.

**Límites y trampas.** *"Esto es una recomendación, no una verdad… hay que converger el cutoff con 'olla-dft converge'."* El tipo/funcional se deduce por regex del encabezado: un UPF sin esos campos queda como `?` y no se descarta. Los cutoffs sugeridos que declara el UPF son un punto de partida, no una convergencia.

**Referencias.** M. J. van Setten et al., *Comput. Phys. Commun.* 226, 39 (2018) — PseudoDojo; A. Dal Corso, *Comput. Mater. Sci.* 95, 337 (2014) — pslibrary; G. Prandini et al., *npj Comput. Mater.* 4, 72 (2018) — SSSP.
