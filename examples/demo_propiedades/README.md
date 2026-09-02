# Silicio: funciones ópticas, fonones y masa efectiva

Módulos `optics`, `phonons` y `effmass` corridos sobre silicio relajado
(LDA, norma conservada, QE 6.6). La estructura es `Si_relajado.cif`
(a0 = 5.402 Å, salida de la EOS de `demo_calculo`).

    olla-dft optics  Si_relajado.cif --ecutwfc 60 --run
    olla-dft optics  Si_relajado.cif --collect --scissor 0.65
    olla-dft phonons Si_relajado.cif --qgrid 2x2x2 --run
    olla-dft phonons Si_relajado.cif --collect

### Ópticas (`opticas_Si.png`, `OPTICS.dat`)

nscf de 14×14×14 sin simetría, 27 bandas, ensanchamiento 0.1 eV.

Sin scissor el cálculo da ε1(0) = 16.13, muy por encima del 11.7
experimental: es la consecuencia directa de que el gap LDA sea 0.52 eV en
vez de 1.17 eV. El espectro en sí es correcto y se comprobó dos veces:

- regla de suma f: ∫ E·ε2(E) dE = 451.5 eV² frente a (π/2)(ħωp)² = 451.5 eV²
  con la frecuencia de plasmón que reporta epsilon.x (16.95 eV;
  experimental 16.7). Factor 1.000.
- Kramers–Kronig de ε2 reproduce el ε1 de epsilon.x: 16.09 contra 16.13
  (0.3 %).

Con `--scissor 0.65` (1.17 − 0.52) el pico de ε2 se coloca en 4.30 eV, que
es justo el punto crítico E2 del silicio, y ε1(0) baja a 10.44.

### Fonones (`fonones_Si.png`, `FONONES_DOS.dat`, `FONONES_TERMO.dat`)

Malla de q 2×2×2 (8 puntos), DOS interpolada en 12×12×12.

Frecuencias (cm⁻¹) contra dispersión inelástica de neutrones:

    punto        Olla-DFT   experimental
    Γ     TO/LO  508.9      517
    X     TA     140.8      150
    X     LA/LO  406.5      410
    X     TO     455.9      463
    L     TA     107.7      114
    L     LA     372.3      378
    L     LO     408.5      417
    L     TO     484.8      490

Todas por debajo entre 1 y 6 %, que es lo esperado de LDA con una malla de
q tan pequeña, y con las degeneraciones correctas. Sin frecuencias
imaginarias: la estructura está bien relajada.

Termodinámica armónica por celda (2 átomos):
- energía de punto cero 122.25 meV
- a 300 K: C_v = 0.411 meV/K, S = 0.422 meV/K, F = 65.0 meV

El C_v experimental del Si a 300 K es 20 J/(mol·K), que por celda de dos
átomos son 0.415 meV/K: 1 % de diferencia.

### Masa efectiva (`MASA_EFECTIVA.dat`)

    olla-dft gen Si_relajado.cif --preset bands -o bandas    # y correrlo
    olla-dft effmass Si_relajado.cif --bands-dir bandas -o masa --run

Camino fino de 6 líneas, ±0.06 Å⁻¹, 21 puntos cada una.

    masa                 Olla-DFT   referencia
    electrón long.       0.949      0.916  (experimental)
    electrón transv.     0.193      0.190  (experimental)
    hueco pesado [100]   0.269      0.277  (Luttinger)
    hueco pesado [111]   0.670      0.718  (Luttinger)

Las dos transversales salen idénticas, como exige la simetría del valle.
El ajuste rápido sobre el camino de bandas normal NO sirve para publicar:
la ventana sale de 0.35 Å⁻¹, fuera del régimen parabólico, y el reporte lo
marca como no confiable.

### Archivos

| Archivo | Qué es |
|---|---|
| `Si_relajado.cif` | estructura de partida (silicio relajado, a0 = 5.402 Å) |
| `OPTICS.dat` | ε1, ε2, n, k, α y R frente a la energía (promedio isótropo) |
| `opticas_Si.png` | figura de las funciones ópticas |
| `FONONES_DOS.dat` | densidad de estados de fonones |
| `FONONES_TERMO.dat` | F, U, S y Cv frente a T en la aproximación armónica |
| `fonones_Si.png` | dispersión y DOS de fonones |
| `MASA_EFECTIVA.dat` | masas efectivas ajustadas con su ventana y su calidad |
