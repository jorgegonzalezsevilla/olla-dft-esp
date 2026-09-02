# Silicio: bandas, DOS y PDOS de principio a fin

Ejemplo completo sobre silicio: generación de los inputs de Quantum
ESPRESSO, ejecución y figuras de bandas y densidad de estados listas para
una revista.

**Generación de los inputs:**

    olla-dft gen Si.cif -p all -o . --insulator

**Cálculo** (pw.x, dos.x, projwfc.x y bands.x en orden):

    ./run.sh

**Figuras:**

    olla-dft plot . -o . --gap-label            # -> Si_bandas_dos
    olla-dft plot . -o . --gap-label --mono     # -> versión monocroma

Resultado: gap indirecto de 0.524 eV, VBM en Γ y CBM sobre Γ–X. Los PDF
miden exactamente el ancho de columna pedido y llevan las fuentes
incrustadas como TrueType, listos para enviar a una revista.

### Archivos

| Archivo | Qué es |
|---|---|
| `scf.in`, `nscf.in`, `bands.in` | inputs de pw.x escritos por `olla-dft gen` (scf, nscf denso y camino de bandas) |
| `KPATH.txt` | el camino de alta simetría usado en `bands.in`, con etiquetas |
| `run.sh` | guion de ejecución generado; espera también `dos.in`, `projwfc.in` y `bands_pp.in`, que `gen -p all` escribe y no se incluyen aquí |
| `Si_bandas_dos.pdf`, `.png` | bandas + DOS en una sola figura, con el gap etiquetado |
| `Si_bandas_dos_mono.pdf` | la misma figura en monocromo (`--mono`) |
| `Si_bandas.pdf`, `Si_dos.pdf` | bandas y DOS por separado |
