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
