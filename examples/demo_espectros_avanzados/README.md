# Espectros y módulos avanzados: XANES, U de Hubbard, electrón-fonón, desdoblamiento de bandas y VDOS

Resultados REALES de cinco módulos (`xanes`, `hubbard`, `elph`, `unfold`,
`md`), con lo que hay que mirar en cada uno.

### XANES (`xanes_Si.png`, `XANES.dat`)

Borde K del silicio. Las tres polarizaciones están dibujadas y no se ven
porque coinciden: el silicio es cúbico, y que coincidan al 0.04 % es la
comprobación de que el manejo de la polarización es correcto.

- Borde a +1.1 eV sobre E_F
- Máximo a +3.7 eV
- Estructuras a +10.4, +12.5 y +17.2 eV

Se reproduce con:

    olla-dft corehole Si --edge K -o pseudos --functional PZ --rcut 1.6
    olla-dft xanes Si8.cif --element Si --core-hole pseudos/Si.hueco1s.UPF --average -o xanes --ecutwfc 40
    # correr pw.x y luego xspectra.x sobre los tres xspectra_*.in
    olla-dft xanes Si8.cif --collect -o xanes --element Si --edge K

### U de Hubbard (`HUBBARD_U.dat`)

U de Hubbard del NiO por respuesta lineal (hp.x), proyección ortho-atomic.
El ciclo de autoconsistencia da:

    iter 0:  5.4429 eV      <- lo que reporta un cálculo de UNA vuelta
    iter 1:  3.9429
    iter 2:  4.1323
    iter 3:  4.1087 eV      <- el autoconsistente

1.33 eV de diferencia. Ese es el argumento entero del módulo.

    olla-dft hubbard NiO.cif -o hub --qgrid 2x2x2          # prepara scf + hp.x
    # correr pw.x y hp.x
    olla-dft hubbard NiO.cif --collect -o hub --qgrid 2x2x2 # lee el U de UNA vuelta
    olla-dft hubbard NiO.cif --cycle -o hub --qgrid 2x2x2   # ciclo scf -> hp.x -> scf hasta converger

### Electrón-fonón (`ELPH_Al_lambda.dat`)

Acoplamiento electrón-fonón del aluminio, malla q 2×2×2. La columna que
importa no es una: es la serie contra el ensanchamiento, y hay que leer el
PLATÓ (aquí en torno a 0.35). Que λ suba de 0.018 a 0.35 y luego se
estabilice es lo normal; si NO se estabilizara, la malla de k sería
insuficiente y cualquier número sería arbitrario.

De ahí sale τ(300 K) = 11.4 fs, que es lo que sustituye al τ constante de
la aproximación CRTA del módulo de transporte.

    olla-dft elph Al.cif --qgrid 2x2x2 -o elph        # prepara scf, nscf y ph.x
    # correr pw.x y ph.x
    olla-dft elph Al.cif --collect -o elph

### Desdoblamiento de bandas (`UNFOLD.dat`)

Desdoblamiento de una supercelda 2× de silicio SIN defecto. Los pesos salen
exactamente 1.0 o 0.0 (0.5/0.5 en las degeneraciones del borde de zona), y
la suma sobre las bandas da exactamente 4.000 en los 21 puntos k. Es el
teorema: con N = 2, la mitad del peso sobrevive.

Sobre una supercelda CON defecto los pesos salen repartidos, y esa
difuminación es el resultado físico.

    olla-dft unfold carpeta_bandas_supercelda/ prim.cif -o . --bands 8 --format png

### VDOS desde dinámica molecular (`MD_VDOS.dat`)

Densidad de estados vibracional del silicio a 900 K, desde una trayectoria
de dinámica molecular. A diferencia de los fonones armónicos, incluye la
anarmonicidad y la temperatura. La resolución en frecuencia la fija la
duración de la trayectoria, y el módulo la reporta: con 0.24 ps son
138 cm⁻¹, que es poco. Una trayectoria de 20 ps daría 1.7 cm⁻¹.

    olla-dft md md.out -o . --skip 50 --no-plot

### Archivos

| Archivo | Qué es |
|---|---|
| `XANES.dat` | σ(E) del borde K del Si: promedio y las tres polarizaciones |
| `xanes_Si.png` | figura del espectro XANES |
| `HUBBARD_U.dat` | U por sitio del NiO tal como lo escribe `hubbard --collect` (una vuelta, 5.44 eV) |
| `ELPH_Al_lambda.dat` | λ, ∫α²F, ⟨log ω⟩ y N(E_F) frente al ensanchamiento del Al |
| `UNFOLD.dat` | peso espectral desdoblado: distancia en el camino, E−E_F y peso |
| `MD_VDOS.dat` | VDOS del Si a 900 K desde la trayectoria (250 pasos, dt = 0.97 fs) |
