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
