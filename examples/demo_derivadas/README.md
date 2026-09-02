# Silicio: Debye y Slack desde las Cij, cuasi-armónica y ficha del material

Módulos derivados (`derived`, `qha`, `datasheet`): propiedades que salen de
resultados ya calculados, sin ningún cálculo nuevo de Quantum ESPRESSO.

Todo lo de esta carpeta es POST-PROCESO: no cuesta ningún cálculo nuevo de
Quantum ESPRESSO, sale de resultados que ya tenías.

### Qué correr

Temperatura de Debye, velocidades del sonido y conductividad térmica de
Slack a partir de las constantes elásticas:

    olla-dft derived ../demo_propiedades/Si_relajado.cif --cij ELASTIC_C.dat

Expansión térmica cuasi-armónica: V(T), α(T), a(T), Cv y Cp:

    olla-dft qha QHA_entrada_Si.dat --cells 8 --natoms 2 --cubic -o salida_qha

Ficha del material (reúne en un Markdown todo lo que hay en la carpeta de un
proyecto):

    olla-dft datasheet carpeta_del_proyecto/ -o .

### Qué mirar

En `derived`, que la temperatura de Debye elástica (~635 K para el silicio,
645 K experimental) NO es la misma que sale de la DOS de fonones: son dos
definiciones distintas, y el reporte lo dice.

En `qha`, que α(T) sale NEGATIVA por debajo de ~165 K. No es un error: el
silicio de verdad tiene expansión térmica negativa por debajo de ~120 K.
Que el modelo lo reproduzca —con el cruce por cero corrido, porque las
frecuencias vienen de un potencial aprendido y no de DFPT— es justo la
señal de que la cuasi-armónica está haciendo su trabajo.

### Archivos

| Archivo | Qué es |
|---|---|
| `ELASTIC_C.dat` | matriz elástica del silicio (LDA, celda primitiva), tal como la escribe `olla-dft elastic --collect` |
| `QHA_entrada_Si.dat` | tabla para `olla-dft qha`: una línea por volumen, con V (Å³), E (eV) y las 48 frecuencias (cm⁻¹) de una supercelda 2×2×2 (8 celdas primitivas). Se generó con MACE-MP-0 en segundos, no con DFPT |
| `qha_Si.png` | la figura que sale del comando `qha` de arriba |
| `ficha_Si.md` | ejemplo de la salida de `olla-dft datasheet` |
