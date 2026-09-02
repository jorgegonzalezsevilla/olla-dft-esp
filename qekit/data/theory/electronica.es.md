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
