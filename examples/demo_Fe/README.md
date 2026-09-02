# Hierro bcc: DOS con polarización de espín

Hierro bcc ferromagnético: generación de los inputs con momento inicial y
densidad de estados resuelta por espín.

    olla-dft gen Fe.cif -p dos -o . --mag Fe=0.7
    olla-dft dos . -o . --journal aps --emin -10 --emax 6

pw.x converge a 2.28 μB/celda; integrando la DOS resuelta por espín que
exporta Olla-DFT se recuperan 2.27 μB de forma independiente.

### Archivos

| Archivo | Qué es |
|---|---|
| `scf.in`, `nscf.in` | inputs de pw.x escritos por `olla-dft gen` (scf con `starting_magnetization` y nscf denso para la DOS) |
| `Fe_dos.pdf`, `Fe_dos.png` | DOS por espín (arriba/abajo) con el estilo de la APS |
