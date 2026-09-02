# Validación

Cada número de este documento sale de una corrida real de Quantum ESPRESSO
6.6, compilado desde el código fuente, preparada, ejecutada y post-procesada
con la propia Olla-DFT: generar los inputs, correr QE, leer los resultados.
Nada de esto es una maqueta. Los datos crudos están en el repositorio: las
carpetas de `examples/` llevan los inputs, los `.dat` y las figuras de cada
caso junto con los comandos `olla-dft` exactos que los produjeron, y
`tests/datos/` guarda las salidas reales de QE (bandas, DOS, fonones,
`epsilon.x`, XANES, TDDFPT, `pwcond.x`, electrón-fonón, MD) que lee la suite de
pruebas. Los valores de referencia que se validaron una vez contra experimento
o literatura están congelados en `tests/referencias.py`, donde hacen de
detectores de regresiones. Los inputs se verificaron además con el parser de
Quantum ESPRESSO incluido en ASE.

Se eligieron tres sistemas para cubrir los casos difíciles: un semiconductor
(Si), un metal (Al) y un metal con polarización de espín (Fe bcc); las
secciones restantes extienden la comprobación a espectros, fonones,
superficies, funciones de Wannier, polarización, transporte térmico y figuras.

## Silicio (semiconductor): bandas y DOS

`scf → nscf → dos.x → projwfc.x → bands → bands.x` corrieron sin editar un solo
archivo generado (`examples/demo_Si/`). Olla-DFT detecta un gap indirecto con
el VBM en Γ y el CBM sobre Γ→X —el mínimo Δ característico del Si— y un gap
directo en Γ; los dos concuerdan con los valores LDA conocidos. Integrando la
DOS hasta E_F se recuperan los ocho electrones de valencia.

| Magnitud | Olla-DFT | Referencia | Fuente de la referencia |
|---|---|---|---|
| Gap indirecto (VBM Γ, CBM sobre Γ→X) | 0.524 eV | valor LDA | resultado LDA conocido del Si |
| Gap directo en Γ | 2.56 eV | valor LDA | resultado LDA conocido del Si |
| Electrones de valencia por ∫DOS hasta E_F | 7.98 de 8 | 8 | conteo de electrones (0.2 % de error de malla) |

## Aluminio (metal)

El análisis lo clasifica correctamente como metálico por bandas que cruzan
E_F, con la forma √E de electrón casi libre. El origen de energías cambia solo
a la energía de Fermi, como corresponde a un metal.

| Magnitud | Olla-DFT | Referencia | Fuente de la referencia |
|---|---|---|---|
| Clasificación | metal (bandas cruzan E_F) | metal | — |
| DOS(E_F) | 0.44 estados/eV | forma √E de electrón libre | modelo de electrón casi libre |

## Hierro bcc (metal con polarización de espín)

Con `--mag Fe=0.7` el cálculo converge a 2.28 μB/celda (`examples/demo_Fe/`).
Integrando la DOS resuelta por espín, Olla-DFT recupera el momento de forma
independiente, lo que confirma que los dos canales se leen y separan
correctamente.

| Magnitud | Olla-DFT | Referencia | Fuente de la referencia |
|---|---|---|---|
| Momento magnético según pw.x | 2.28 μB/celda | 2.22 μB | experimento |
| Momento por integración de la DOS resuelta por espín | 2.27 μB | 2.28 μB | el valor de pw.x de arriba |

## Ópticas (Si)

El espectro de `epsilon.x` (`examples/demo_propiedades/`) pasa dos pruebas
independientes: la regla de suma f, ∫E·ε₂(E)dE = (π/2)ħω_p², se cumple con
factor 1.000 frente a la frecuencia de plasmón que el propio código reporta, y
la transformada de Kramers-Kronig de ε₂ reproduce el ε₁ que escribe
`epsilon.x`. Con el scissor de 0.65 eV (gap experimental menos gap LDA) el pico
de ε₂ cae exactamente en el punto crítico E₂ del silicio.

| Magnitud | Olla-DFT | Referencia | Fuente de la referencia |
|---|---|---|---|
| Factor de la regla de suma f | 1.000 | 1 | regla de suma exacta |
| Frecuencia de plasmón | 16.95 eV | 16.7 eV | experimento |
| ε₁ por Kramers-Kronig frente al ε₁ de `epsilon.x` | 0.3 % de diferencia | 0 | relación KK analítica |
| Pico de ε₂ con scissor de 0.65 eV | 4.30 eV | 4.30 eV | punto crítico E₂ del Si |

## Fonones (Si)

DFPT con malla de q 2×2×2 (`examples/demo_propiedades/`). Frecuencias en Γ, X,
L y W dentro del 1–6 % de los datos de neutrones, con las degeneraciones
correctas y sin frecuencias imaginarias; el C_v a 300 K coincide con el
experimental en 1 %.

| Magnitud | Olla-DFT | Referencia | Fuente de la referencia |
|---|---|---|---|
| Γ TO/LO | 508.9 cm⁻¹ | 517 cm⁻¹ | dispersión inelástica de neutrones |
| X TA | 140.8 cm⁻¹ | 150 cm⁻¹ | dispersión inelástica de neutrones |
| X LA/LO | 406.5 cm⁻¹ | 410 cm⁻¹ | dispersión inelástica de neutrones |
| X TO | 455.9 cm⁻¹ | 463 cm⁻¹ | dispersión inelástica de neutrones |
| L TA | 107.7 cm⁻¹ | 114 cm⁻¹ | dispersión inelástica de neutrones |
| L LA | 372.3 cm⁻¹ | 378 cm⁻¹ | dispersión inelástica de neutrones |
| L LO | 408.5 cm⁻¹ | 417 cm⁻¹ | dispersión inelástica de neutrones |
| L TO | 484.8 cm⁻¹ | 490 cm⁻¹ | dispersión inelástica de neutrones |
| C_v a 300 K (celda de 2 átomos) | 0.411 meV/K | 0.415 meV/K (20 J/(mol·K)) | experimento |

## Función trabajo (grafeno)

| Magnitud | Olla-DFT | Referencia | Fuente de la referencia |
|---|---|---|---|
| Φ (grafeno) | 4.54 eV | 4.6 eV | experimento |

## Difracción de rayos X

Contrastada contra las fichas PDF: Si (111)/(220)/(311)/(400) y NaCl
(200)/(220)/(222)/(400)/(111)/(311) con Δ2θ < 0.05°. Los índices se dan en la
celda convencional aunque la entrada sea la primitiva, que es lo que hace
comparables los hkl con la literatura (`examples/demo_laminar/`).

| Magnitud | Olla-DFT | Referencia | Fuente de la referencia |
|---|---|---|---|
| Posiciones de Si (111), (220), (311), (400) | Δ2θ < 0.05° | PDF 27-1402 | ficha PDF del ICDD |
| Posiciones de NaCl (200), (220), (222), (400), (111), (311) | Δ2θ < 0.05° | PDF 05-0628 | ficha PDF del ICDD |

## Sólidos amorfos (a-SiO₂)

Un fundido y temple de 24 átomos a 2.2 g/cm³ con MACE (1.8 ps) da exactamente
la red aleatoria continua de tetraedros SiO₄ que comparten vértice, que es la
estructura del vidrio de sílice, sin un solo enlace O–O. La coordinación se
mide con un corte por par tomado de los radios covalentes, no con un único
radio global; con un corte global de 3 Å los oxígenos aparecerían "enlazados"
entre sí, que es el error clásico de este análisis.

| Magnitud | Olla-DFT | Referencia | Fuente de la referencia |
|---|---|---|---|
| Coordinación Si–O | 4.00 | 4 | tetraedros SiO₄ |
| Coordinación O–Si | 2.00 | 2 | red que comparte vértices |
| Enlaces O–O | 0 | 0 | estructura del vidrio de sílice |
| Distancia Si–O | 1.690 Å | 1.61 Å | experimento |

## Funciones de Wannier (Si)

Las cuatro bandas de valencia del silicio se wannierizan sobre orbitales s
centrados en los enlaces. Los centros salen a 0.6788 Å de cada átomo en las
tres direcciones, o sea a 1.1756 Å = √3·a/8: el punto medio del enlace, con
cuatro decimales y sin que nada lo haya impuesto — es la fase de Berry la que
lo pone ahí. La dispersión converge con la malla y el valor de 6³ coincide con
Marzari y Vanderbilt. La interpolación reproduce las bandas de DFT en puntos
que NO estaban en la malla; transformar directamente las energías propias —sin
gauge— es entre 3.5 y 5.9 veces peor. Además Ω se parte en Ω_I + Ω_D + Ω_OD
exactamente, Ω_I no se mueve más de 10⁻¹² Å² al minimizar (es invariante de
gauge, y que se moviera sería la señal de un error), y el gradiente del
funcional de dispersión se contrastó contra su propia derivada numérica: ahí
apareció que le faltaba un factor 1/N_k, que con una malla 4×4×4 hace el paso
64 veces demasiado largo.

| Magnitud | Olla-DFT | Referencia | Fuente de la referencia |
|---|---|---|---|
| Distancia del centro de Wannier al átomo | 1.1756 Å | 1.1756 Å (√3·a/8) | punto medio del enlace, geometría |
| Dispersión por función, malla 4³ / 6³ / 8³ | 1.605 / 1.901 / 2.047 Å² | 1.93 Å² | Marzari y Vanderbilt |
| Error máximo de interpolación con gauge, 4³ / 6³ / 8³ | 275 / 108 / 48 meV | bandas DFT fuera de la malla | energías propias directas de pw.x |
| Error máximo sin gauge, 4³ / 6³ / 8³ | 962 / 473 / 281 meV | — | ídem |
| Cambio de Ω_I al minimizar | ≤ 10⁻¹² Å² | 0 | invariancia de gauge |

## Polarización por fase de Berry (Si y BN cúbico)

La parte iónica se contrasta contra su fórmula exacta, Σ Z_a·f_a, y coincide
con lo que escribe pw.x hasta el último decimal que imprime (5·10⁻⁶) en los
dos sistemas — incluido el plegado ion por ion módulo 1 que hace Quantum
ESPRESSO cuando alguna carga de valencia es impar, que es también lo que parte
el cuanto de polarización por la mitad. Al desplazar un átomo de silicio
0.16 Å la fase iónica se mueve 0.204 y la electrónica se mueve −0.204: la
carga efectiva de Born sale de una cancelación, no de que no pase nada. Y la
fase electrónica del silicio distorsionado, calculada con `lberry`, coincide
con la que dan los centros de Wannier del mismo sistema: dos rutas que no
comparten una línea de código y llegan al mismo número, porque son la misma
fase de Berry.

| Magnitud | Olla-DFT | Referencia | Fuente de la referencia |
|---|---|---|---|
| Fase iónica frente a Σ Z_a·f_a | coincide a 5·10⁻⁶ | salida de pw.x | último decimal que imprime pw.x |
| Z* del Si (átomo desplazado 0.16 Å) | 0 a 10⁻¹⁴ e | 0 | exacto en un cristal homopolar |
| Z* del BN cúbico, malla 6×6, 11 puntos por cuerda, 60 Ry | 1.94 e (2.01 con la malla gruesa) | 1.92 e | literatura |
| Fase electrónica, `lberry` frente a centros de Wannier | 0.0884 frente a 0.0892 | la misma fase | dos rutas independientes |

## Conductividad térmica de red (Si)

Las 57 configuraciones desplazadas de una supercelda 2×2×2 se calcularon con
pw.x y la ecuación de Boltzmann de fonones en RTA da κ a 300 K por debajo de lo
medido en la diferencia que se espera de la RTA (10–15 % por debajo de la
solución exacta) más el corte de la fc3. La dependencia con la temperatura es
el T⁻¹ de los procesos Umklapp, y la mitad de κ la llevan fonones con recorrido
libre medio mayor que 1.0 µm, que es justo lo que miden los experimentos de
espectroscopía de recorrido libre en silicio — es el número que dice por qué
nanoestructurar el silicio funciona tan bien para termoeléctricos. El mismo
cálculo con fuerzas de MACE en lugar de DFT tarda 8 segundos en vez de 40
minutos, reproduce el exponente y falla el valor absoluto por un factor 2: por
eso el informe lo dice cada vez que las fuerzas no vienen de DFT. La
convergencia con el tamaño de la supercelda se comprobó con ese camino barato
(2×2×2 → 3×3×3 mueve κ de 50.1 a 50.8 W/m·K, o sea nada), que es exactamente
para lo que sirve.

| Magnitud | Olla-DFT | Referencia | Fuente de la referencia |
|---|---|---|---|
| κ a 300 K, fuerzas DFT, RTA | 101 W/m·K (96 con isótopos naturales) | ~140 W/m·K | experimento |
| Exponente de temperatura, fuerzas DFT | κ ∝ T⁻¹·¹⁶ | T⁻¹ | procesos Umklapp |
| Recorrido libre medio que lleva la mitad de κ | 1.0 µm | ~1 µm | espectroscopía de recorrido libre |
| κ a 300 K, fuerzas MACE | 51 W/m·K | 101 W/m·K (DFT) | este trabajo |
| Exponente de temperatura, fuerzas MACE | κ ∝ T⁻¹·⁰⁶ | T⁻¹ | procesos Umklapp |
| Convergencia de supercelda (MACE), 2×2×2 → 3×3×3 | 50.1 → 50.8 W/m·K | convergido | este trabajo |

## Superficies con ESM (Al(111))

Con `esm_bc='bc1'` el nivel de vacío vale cero por construcción, así que la
función trabajo es directamente −E_F sin ajustar ninguna meseta, y deja de
depender del vacío — con 8, 12 y 16 Å la energía cambia 6·10⁻⁶ Ry y E_F
0.4 meV —, así que se puede trabajar con media celda. El módulo se niega a
correr `bc1` con carga neta, que es un problema mal planteado y que pw.x
calcula igualmente (dio −379 y −677 Ry para la misma losa con dos vacíos
distintos), y centra la losa en z = 0 antes de escribir nada, porque ESM mide
z desde el centro de la celda y una losa dejada donde ASE la pone cae sobre su
frontera. Con `bc3` y carga, comprueba que Φ(q) sea una recta antes de dar una
capacitancia: en Al(111) con ±0.04 e todavía no lo es (16 % de desviación), y
en ese caso lo dice en vez de dar el número.

| Magnitud | Olla-DFT | Referencia | Fuente de la referencia |
|---|---|---|---|
| Φ de Al(111), `bc1` | 4.24 eV | 4.24–4.26 eV | experimento |
| Cambio de energía, vacío 8/12/16 Å | 6·10⁻⁶ Ry | 0 | independencia del vacío en ESM |
| Cambio de E_F, vacío 8/12/16 Å | 0.4 meV | 0 | independencia del vacío en ESM |
| Linealidad de Φ(q) a ±0.04 e | 16 % de desviación (se avisa, no se da C) | recta | electrostática del condensador plano |

## Desenredado de Wannier (Si)

Con ocho funciones sp³ sacadas de doce bandas de DFT, la proyección sola deja
las bandas de valencia lejos de las de DFT; eligiendo el subespacio por el
método de Souza-Marzari-Vanderbilt el error baja un factor 4.3. Ω_I baja y
luego NO se mueve durante la minimización de la dispersión, que es lo que tiene
que pasar: es invariante de gauge. Las bandas dentro de la ventana congelada se
reproducen exactas, como promete el método.

| Magnitud | Olla-DFT | Referencia | Fuente de la referencia |
|---|---|---|---|
| Error de las bandas de valencia, solo proyección | 899 meV | — | bandas DFT |
| Error de las bandas de valencia, tras el desenredado | 208 meV | factor 4.3 mejor | bandas DFT |
| Ω_I antes → después del desenredado | 12.37 → 10.37 Å² | baja | Souza-Marzari-Vanderbilt |
| Cambio de Ω_I al minimizar | 5·10⁻¹⁴ Å² | 0 | invariancia de gauge |
| Error dentro de la ventana congelada | 3·10⁻¹³ eV | 0 | exacto por construcción |

## ESM cargado (Al(111)): capacitancia

La capacitancia se validó sin recurrir a ningún número de la literatura, con
electrostática: para un condensador plano 1/C = d/ε₀, y esa pendiente no
depende del material ni del funcional. Midiendo C con cuatro separaciones
distintas entre 4 y 11 Å sale una recta con pendiente 1/ε₀. Eso valida a la
vez la fórmula, el área y la conversión de unidades — y por tanto la rama
cargada del módulo.

| Magnitud | Olla-DFT | Referencia | Fuente de la referencia |
|---|---|---|---|
| Pendiente de 1/C frente a d (4 separaciones, 4–11 Å) | 1/ε₀ con 0.4 % de error | 1/ε₀ | electrostática del condensador plano |
| R² del ajuste lineal | 0.99998 | 1 | — |

## Figuras

Las figuras se comprobaron midiendo el PDF resultante (178.0 mm exactos para
doble columna, 86.0 mm para una columna), inspeccionando que las únicas fuentes
incrustadas sean las de la familia elegida y que lo estén como TrueType, y
revisando la legibilidad a tamaño de impresión real y no ampliadas. La paleta
se validó con un verificador de separación de color que simula protanopia y
deuteranopia.

| Magnitud | Olla-DFT | Referencia | Fuente de la referencia |
|---|---|---|---|
| Ancho de figura a doble columna | 178.0 mm | 178 mm | especificación de revista |
| Ancho de figura a una columna | 86.0 mm | 86 mm | especificación de revista |
| Fuentes incrustadas | solo la familia elegida, TrueType | — | inspección del PDF |
| Separación de color de las cuatro primeras series (ΔE en OKLab) | ≥ 11 | ≥ 8 (umbral seguro) | simulación de protanopia y deuteranopia |

## Comprobaciones incorporadas: `olla-dft selftest`

Las pruebas de pytest miran que el código haga lo que el código dice;
`selftest` compara Olla-DFT con el mundo. Cada prueba calcula una magnitud que
alguien ha medido o deducido y la contrasta con ese valor y su fuente.
`olla-dft selftest` corre las rápidas (sin Quantum ESPRESSO, segundos);
`--full` añade las que corren pw.x de verdad sobre sistemas pequeños (del
orden de diez minutos; hace falta un `pw.x` que funcione y `--pseudo-dir`);
`--mlip` añade, por separado, la que necesita MACE; `--list` imprime la tabla
de abajo sin correr nada. Las tolerancias son relativas salvo cuando la
referencia es cero, donde son absolutas.

| Clave | Prueba | Referencia | Tolerancia | Fuente | Necesita |
|---|---|---|---|---|---|
| `madelung` | Constante de Madelung, red cúbica simple, α_M | 2.8372974 | 1·10⁻⁵ | valor clásico de la suma de Ewald para una carga puntual en un fondo neutralizante | — |
| `lorenz` | Número de Lorenz de un gas de electrones libres, L/L₀ | 1.0 | 12 % | límite de Sommerfeld, L₀ = (π²/3)(k_B/e)² = 2.44·10⁻⁸ W·Ω/K² | — |
| `npw` | Ondas planas de Si a 30 Ry, N_PW | 725 | 6 % | lo que reporta pw.x para la celda primitiva de Si (V = 39.5 Å³) a 30 Ry | — |
| `sackur` | Entropía traslacional del N₂ a 298 K | 150.4 J/(mol·K) | 1 % | Sackur-Tetrode a 1 bar; tablas NIST-JANAF | — |
| `allen_dynes` | T_c de Allen-Dynes para el aluminio | 1.18 K | 12 % | T_c experimental del Al con λ = 0.44, ω_log = 270 K (Allen-Dynes 1975) y µ* = 0.12; µ* es un parámetro ajustado, con 0.10 la misma fórmula da 1.9 K | — |
| `allen_dynes_mu` | Sensibilidad de Allen-Dynes a µ*, T_c(0.10)/T_c(0.12) | 1.56 | 5 % | la fórmula es exponencial en µ*: subirlo de 0.10 a 0.12 baja T_c a dos tercios | — |
| `born2d` | Módulos de lámina de una hoja isótropa, Y_2D | 341.8 N/m | 1 % | C11 = 352, C12 = 60 N/m (grafeno, DFT), Y = C11 − C12²/C11 | — |
| `gap_invariante` | El alineamiento quita el cero arbitrario, ΔE_v de un material consigo mismo | 0 eV | 1·10⁻⁹ | identidad exacta | — |
| `ewald_escala` | La constante de Madelung no depende de la escala, \|α(L=3) − α(L=30)\| | 0 | 1·10⁻⁶ | invariancia exacta de la suma de Ewald bajo un cambio de unidades | — |
| `chern_qwz` | Chern del aislante de Qi-Wu-Zhang (banda inferior, m = −1) | −1 | 1·10⁻¹⁰ | Qi, Wu y Zhang, Phys. Rev. B 74, 085308 (2006) | — |
| `her_pt` | HER: el platino está en la cumbre del volcán, ΔG_H* | −0.09 eV | 5 % | Nørskov y col. 2005, Pt(111) | — |
| `oer_ruo2` | OER: sobrepotencial del RuO₂(110), η | 0.48 V | 10 % | Man et al. 2011 (ChemCatChem), ΔG(OH) = 0.77, ΔG(O) = 2.16, ΔG(OOH) = 3.87 eV | — |
| `escala_oer` | Relación de escala OOH−OH de la OER, ΔG(OOH) − ΔG(OH) | 3.2 eV | 10 % | relación universal de escala, 3.2 ± 0.2 eV en casi toda superficie de óxido | — |
| `escala_eta_min` | Límite de escala del sobrepotencial de la OER, η_min | 0.37 V | 2 % | Man et al. 2011 | — |
| `umklapp` | κ_L del silicio decae como 1/T, exponente n | 1.0 | 25 % | por encima de la temperatura de Debye los procesos Umklapp dan κ ∝ 1/T; se comprueba el exponente, no κ | MACE (`--mlip`, ~25 s) |
| `fonon_si` | Modo óptico del Si en Γ, ω(Γ) | 520 cm⁻¹ | 10 % | Raman experimental del silicio, 520.7 cm⁻¹ a 300 K | QE (~20 s) |
| `wannier_si` | Centro de Wannier del enlace Si–Si, \|r̄\| | 1.17563 Å | 2 % | centro del enlace de la estructura diamante a √3·a/8 con a = 5.43 Å | QE (~30 s) |
| `condensador` | ESM cargado: pendiente de 1/C frente a la distancia sobre 1/ε₀ | 1.0 | 6 % | electrostática del condensador plano, independiente del material, del pseudopotencial y del funcional | QE (~90 s) |
| `born_si` | Carga efectiva de Born del silicio, Z* | 0 e | 0.05 | cero exacto en un cristal homopolar por la regla de suma acústica | QE (~60 s) |
| `gamma_al` | Energía de superficie de Al(111), γ | 1.10 J/m² | 25 % | LDA de potencial completo (Vitos et al. 1998) da 1.20 J/m²; el experimento policristalino, 1.14 | QE (~60 s) |
| `bulk_si` | Módulo de bulto del Si por deformación, B | 95 GPa | 15 % | LDA da 93–97 GPa (Nielsen y Martin 1985); el experimento, 98 | QE (~50 s) |
| `sitio_h_al` | H sobre Al(111): el hueco gana al top, E_ads(top) − E_ads(hueco) | 5.6 eV | 60 % | el hidrógeno quimisorbe en el hueco de fcc(111); el orden hueco < puente < top es de manual | QE (~60 s) |

Las que salen mal no siempre son un fallo del código: una tolerancia ajustada,
un pseudopotencial distinto o un cutoff bajo también las mueven. Lo que sí
quieren decir es que ese número ha cambiado y hay que mirar por qué.

## Suite de pruebas

`tests/` contiene 977 pruebas de pytest que corren sin Quantum ESPRESSO
(`python -m pytest -q`, menos de un minuto). Leen las salidas reales de QE
de `tests/datos/`, comparan contra las referencias congeladas de
`tests/referencias.py`, validan cada comando `olla-dft` citado en los README de
`examples/` y en las recetas contra el árbol real de argparse, y lanzan el
programa entero con la salida forzada a cp1252 para asegurar que ningún
informe muere en una consola heredada de Windows. `tests/barrido_cli.sh` es el
barrido de regresión complementario a nivel de comando sobre salidas de QE ya
calculadas, con el código de salida esperado declarado en cada línea.
