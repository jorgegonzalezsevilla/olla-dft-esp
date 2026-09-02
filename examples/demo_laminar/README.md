# Grafito: capas, difractograma y energía de exfoliación

Módulo de materiales laminares (`layers`, `xrd`, `exfoliate`), demostrado
con grafito. La estructura es `../grafito.cif`.

    olla-dft layers ../grafito.cif
      -> 2 capas, d basal = 3.356 Å, (002) esperado en 26.54° (Cu Kα)

    olla-dft xrd ../grafito.cif --size 18 --exp demo_experimental.xy
      -> difractograma simulado + comparación con un "experimental"
         (el "experimental" de aquí es sintético, generado solo para la
         demo: el mismo grafito con el espaciado basal un 2 % mayor,
         ruido y fondo)

    olla-dft exfoliate ../grafito.cif --run
      -> E_exf(LDA) = 25.8 meV/átomo = 0.157 J/m² (valor de literatura LDA)

El XRD está verificado numéricamente contra pymatgen (Si y grafito, todos
los picos con Δ2θ < 0.05° y ΔI < 1.5).

### Archivos

| Archivo | Qué es |
|---|---|
| `XRD_HKL.dat` | lista de reflexiones del grafito: 2θ, d, intensidad y (hkl), λ = 1.54184 Å |
| `demo_experimental.xy` | difractograma "experimental" sintético (2θ, I) para probar `--exp` |
| `xrd.pdf`, `xrd.png` | difractograma simulado frente al sintético |
