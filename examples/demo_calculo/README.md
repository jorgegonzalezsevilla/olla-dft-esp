# Silicio: convergencia, ecuación de estado y constantes elásticas

Los tres módulos de cálculo (`converge`, `eos`, `elastic`) corridos sobre
silicio con LDA y pseudopotencial de norma conservada.

    olla-dft converge Si.cif --kind ecutwfc --values 20,30,40,50,60,70 --run
    olla-dft eos Si.cif --ecutwfc 60 --run
    olla-dft elastic Si_eq.cif --ecutwfc 60 --run

Resultados: el cutoff converge en 50 Ry (1 meV/átomo); a0 = 5.402 Å y
B0 = 94.2 GPa; C11/C12/C44 = 159.9/61.7/76.6 GPa, con B = 94.45 GPa
calculado desde las Cij (0.25 % de diferencia con el de la EOS).
Experimento: 165.8/63.9/79.6 GPa y B ≈ 98 GPa.

### Archivos

| Archivo | Qué es |
|---|---|
| `CONVERGENCIA.dat`, `CONVERGENCIA.txt` | tabla E(ecutwfc) y el reporte de convergencia con el cutoff recomendado |
| `convergencia.png` | figura de la convergencia |
| `EOS.dat`, `EOS.txt` | puntos E(V) y el ajuste de Birch–Murnaghan (V0, B0, B0') |
| `eos.png` | figura de la ecuación de estado |
| `ELASTIC_C.dat`, `ELASTIC.txt` | matriz elástica Cij y el reporte (módulos de Voigt–Reuss–Hill, estabilidad) |
| `elastic.png` | figura esfuerzo–deformación de las deformaciones aplicadas |
