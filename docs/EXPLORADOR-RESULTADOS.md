# Explorador de resultados

La versión 1.2.0 añade un explorador sin conexión a `project dashboard` y
`results explore`. Lee resultados existentes: no ejecuta QE ni modifica inputs.

```sh
olla-dft results ingest ./calculo --project ./mi-proyecto
olla-dft results explore --project ./mi-proyecto -o resultados.html
# Base de resultados normalizados independiente; no requiere manifiesto:
olla-dft --language es results explore --db results.sqlite3 -o resultados.html
```

1. Abre `resultados.html`. Por defecto se muestran cálculos convergidos.
2. Filtra por estado, cálculo, fórmula, etiqueta o ID; selecciona filas en la tabla.
3. Elige métricas X/Y o índice del registro; abre **Personalizar presentación**
   para ajustar título, color, puntos/barras, dimensiones, texto y rango Y.
4. Revisa el contador de exportación. **Exportar solo los puntos visibles** excluye
   valores ausentes, unidades distintas y valores fuera del rango. Al desmarcarlo,
   CSV/JSON/HTML incluyen todos los seleccionados que coinciden con los filtros.
   Las imágenes SVG/PNG siempre corresponden a la figura visible.
5. Descarga SVG editable, PNG a 1×/2×, CSV para hojas de cálculo, JSON con valores
   y procedencia, o HTML para reabrir selección y presentación sin conexión.

## Interpretación científica

Las exportaciones de datos no convierten unidades ni redondean valores. La tabla
acorta la presentación; enfoca un punto o abre sus detalles para ver la precisión
almacenada. Si un rango estrecho necesita un desplazamiento de eje, la figura lo
indica numéricamente y etiqueta el eje delta. El orden del registro no representa
convergencia ni una variable física. Las barras incluyen cero; usa puntos para
mostrar un rango estrecho de energía distinto de cero.

Cada opción de eje corresponde a una unidad exacta. Se cuenta y advierte la omisión
de registros con otra unidad. Los valores ausentes/no finitos son null en JSON o
celdas vacías con razón en CSV. La incertidumbre se conserva en CSV/JSON, pero esta
figura aún no dibuja barras de error. Convergencia y revisión humana son independientes.

La advertencia de composición/método señala una comparación exploratoria. Huellas
iguales no prueban equivalencia física. Las nuevas ingestas incluyen cutoffs (Ry),
malla/desplazamiento k, ocupaciones, espín y smearing. Registros anteriores pueden
carecer de estos parámetros: revisa input, pseudopotenciales y geometría antes de
interpretar diferencias. No se eliminan automáticamente IDs distintos como duplicados.

## Límites y archivos compartidos

El HTML es una instantánea: vuelve a generarlo cuando haya nuevos resultados.
`explore` carga hasta 10000 registros (`--limit` permite reducirlos) y muestra
cargados/total. La figura admite hasta 2000 puntos elegibles; filtra selecciones
mayores. No se muestrean en silencio. La tabla presenta 50 filas por página.
El botón **Ampliar gráfico** facilita leer los ejes en pantallas pequeñas.

Todos los recursos están dentro del archivo; no hay servicios externos, fuentes
remotas, seguimiento ni almacenamiento del navegador. El HTML descargado contiene
solo el alcance de exportación y su presentación. Al copiar un dashboard, conserva
su archivo hermano `.results.html`; el explorador se puede compartir solo.

Se omiten rutas de origen y notas de revisión, pero los títulos y etiquetas pueden
contener texto privado: revísalos antes de compartir. Las huellas SHA-256 identifican
contenido y método registrado; no son firmas digitales. CSV protege textos que
podrían ser fórmulas con un apóstrofo; los valores numéricos negativos se conservan.
Usa UTF-8 con BOM, comas, punto decimal y finales CRLF.

Esta exportación portable es independiente del contrato JSON existente de
`results export`, que conserva registros completos con procedencia local.


El eje de registro sigue el orden de la instantánea. CLI/dashboard conservan en
metadatos `ingested_desc_path_asc`: ingesta más reciente primero y ruta como desempate.
Las selecciones exportadas mantienen el número original; este orden no es físico.
SVG, JSON y HTML guardan versión de Olla-DFT, fecha de generación, esquema y orden.

En Excel configurado con coma decimal, importa con **Desde texto/CSV**, elige coma
como delimitador y una configuración regional que interprete punto decimal en las
columnas de métricas. Evita depender de la apertura con doble clic. En pandas usa
`pd.read_csv("olla-results.csv", encoding="utf-8-sig")`. Una celda numérica vacía
representa un ausente; `<métrica>.reason` explica el motivo y `<métrica>.unit` da
la unidad del valor y de su incertidumbre por fila, incluso al mezclar unidades.
