#!/bin/bash
# Generado por Olla-DFT (olla-dft gen) — ejecuta los cálculos en orden.
set -e
NP=4

echo ">> pw.x < scf.in"
mpirun -np $NP pw.x -in scf.in | tee scf.out

echo ">> pw.x < nscf.in"
mpirun -np $NP pw.x -in nscf.in | tee nscf.out

echo ">> dos.x < dos.in"
mpirun -np $NP dos.x -in dos.in | tee dos.out

echo ">> projwfc.x < projwfc.in"
mpirun -np $NP projwfc.x -in projwfc.in | tee projwfc.out

echo ">> pw.x < bands.in"
mpirun -np $NP pw.x -in bands.in | tee bands.out

echo ">> bands.x < bands_pp.in"
mpirun -np $NP bands.x -in bands_pp.in | tee bands_pp.out

