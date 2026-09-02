# Etileno con TDDFPT y conductancia balística de un hilo de Al

Excitaciones ópticas de una molécula con `tddft` (turbo_davidson.x) y
bandas complejas / canales de conducción de un hilo con `ballistic`
(pwcond.x).

### TDDFPT (`tddft_C2H4.png`, `TDDFT_etileno.eigen`)

Etileno con TDDFPT (turbo_davidson.x, funcional PZ). Seis excitaciones:

    n    E (eV)    f (osc.)   pol.
    1    6.4955     0.02893    x     <- la π->π*, la única realmente brillante
    2    7.1554     0.00000    z     <- oscura
    3    7.1726     0.00011    z
    4    7.2074     0.01237    x
    5    7.4331     0.00005    z
    6    7.4947     0.00001    z

Lo que hay que mirar: cuatro de las seis son OSCURAS. Existen como estados
excitados y no se ven en un espectro de absorción. Un cálculo que solo
reportara "el primer estado excitado está a 6.50 eV" estaría escondiendo
que los siguientes cuatro no aportan nada al espectro.

La π->π* experimental del etileno está a 7.66 eV. 6.50 eV es bajo, y es lo
que se espera: LDA subestima esta transición. El archivo `.eigen` trae las
energías en RYDBERG; leerlas como eV daría 0.48 eV, que es absurdo, y es un
error fácil de cometer.

Se reproduce con:

    olla-dft corehole C --plain -o ps --functional PZ --rcut 1.3
    olla-dft corehole H --plain -o ps --functional PZ --rcut 1.0
    olla-dft tddft c2h4.cif -o run --method davidson --states 6 --pseudo-dir ps --ecutwfc 40 --emax 12
    # correr pw.x y turbo_davidson.x
    olla-dft tddft --collect -o run --method davidson --gap 6.6

### Conductancia balística (`BALISTICO.dat`)

Estructura de bandas compleja de un hilo monoatómico de aluminio
(pwcond.x, ikind=0). La columna "canales" es el número de canales abiertos
a cada energía, que acota la conductancia por arriba: G ≤ canales × G0.

    por debajo de -0.3 eV:  0 canales  (no conduce)
    de -0.3 a  0.2 eV:      1 canal    (la banda s)
    cerca de  +0.3 eV:      3 canales  (entran las p degeneradas)

Ese salto de 1 a 3, y no a 2, es la firma de una degeneración doble: las
p_x y p_y del hilo entran juntas.

Con G0 = 2e²/h = 7.748e-5 S, un solo canal transmitiendo perfectamente son
12.906 kΩ. Para la conductancia de verdad hace falta una región de
dispersión en medio (`--scatterer`) e ikind=1.

    olla-dft ballistic Al_hilo.cif -o . --ikind 0      # prepara scf + pwcond.x
    # correr pw.x y pwcond.x
    olla-dft ballistic --collect -o . --no-plot

### Archivos

| Archivo | Qué es |
|---|---|
| `TDDFT_etileno.eigen` | salida cruda de turbo_davidson.x: energía (Ry) y fuerzas de oscilador total y por componente |
| `tddft_C2H4.png` | espectro de absorción con las seis excitaciones marcadas |
| `BALISTICO.dat` | canales abiertos frente a E−E_F del hilo de Al, tal como lo escribe `ballistic --collect` |
