# Contribuir a Olla-DFT

## Preparar el entorno

```bash
git clone https://github.com/jorgegonzalezsevilla/olla-dft-esp
cd olla-dft-esp
python -m venv .venv && source .venv/bin/activate   # opcional
pip install -e ".[test]"
```

Eso instala el comando `olla-dft` en modo editable más pytest y pyflakes.
Quantum ESPRESSO no hace falta para desarrollar ni para correr la suite de
pruebas; `olla-dft sistema` dice qué tiene tu máquina. Los requisitos y las
notas por sistema están en [docs/PLATAFORMAS.md](docs/PLATAFORMAS.md).

## Antes de abrir un pull request

```bash
python -m pytest -q                 # la suite entera, sin QE, menos de un minuto
python -m pyflakes qekit tests      # no debe imprimir nada
python tools/build_docs.py          # regenera docs/COMANDOS.md y docs/TEORIA.md
```

Corre `build_docs.py` cada vez que añadas, renombres o cambies las opciones de
un comando, o toques `qekit/data/theory/`: `tests/test_docs.py` y
`tests/test_teoria.py` fallan si los archivos generados están viejos. No edites
esos cuatro archivos a mano. Opcional pero bienvenido: `olla-dft selftest`
(segundos, sin QE) y, si tienes `pw.x`, `olla-dft selftest --full --pseudo-dir
/ruta/a/upf` (unos diez minutos).

## Convenciones de estilo

- **Los identificadores y comentarios en español son el estilo de la casa.**
  Funciones, variables, docstrings, textos de ayuda e informes científicos se
  escriben en español, como el código existente (`ErrorDeUso`, `preparar`,
  `informe`). No traduzcas nombres existentes.
- **El inglés va en las tablas de i18n.** Cada texto de interfaz tiene su
  equivalente en inglés en `qekit/data/i18n/`: los textos de ayuda en
  `cli_en.json` (la clave es el texto en español), los resúmenes de una línea
  en `docs_en.json`, las etiquetas de menú, inicio guiado y dashboard en su
  `_en.json`, y la teoría en `qekit/data/theory/*.en.md`. Los informes
  científicos siguen en español. La capa se describe en
  [docs/ARQUITECTURA.md](docs/ARQUITECTURA.md).
- **Toda opción nueva necesita su texto de ayuda** (`help=` en español y su
  traducción en `cli_en.json`). Las opciones que se repiten entre comandos
  toman el suyo de la tabla `defaults` de `cli_es.json`.
- **Todo cambio de física necesita una prueba y una nota en la teoría.** Una
  fórmula, constante, valor por omisión o corrección nuevos llevan una prueba
  en `tests/` (con una salida real de QE en `tests/datos/` si lee una, y un
  valor congelado en `tests/referencias.py` si se validó contra experimento) y
  la actualización de la sección correspondiente en
  `qekit/data/theory/<area>.es.md` **y** `<area>.en.md` (la suite comprueba la
  paridad es/en y los apartados obligatorios). Nunca cambies un valor de
  `tests/referencias.py` solo para que pase una prueba: es un detector de
  regresiones y solo se actualiza cuando el valor nuevo se validó otra vez
  contra la fuente externa.
- Cada `.dat` y cada figura pasan por `core/provenance.py`; cada barrido pasa
  por `modules/sweep.py` y sigue prepare / `--run` / `--collect`.
- Los errores de uso lanzan `ErrorDeUso` (código 2, mensaje limpio, sin
  traza). Todo lo demás es una falla del programa (código 1) y se registra en
  local.
- La salida tiene que sobrevivir a una consola cp1252: si imprimes un símbolo
  nuevo fuera de ASCII, añádelo a `TRANSLITERACION` en `core/consola.py` o
  `test_portabilidad.py` te lo dirá.
- Mantén sincronizados `COMMAND_GROUPS` (cli.py), `docs.GRUPOS`/`MODULO_DE`
  (docs.py) y las secciones de teoría al añadir un comando; la lista paso a
  paso está en [docs/ARQUITECTURA.md](docs/ARQUITECTURA.md#añadir-un-comando).
- Sin emojis en código, salida ni documentación. Líneas de unas 79 columnas.

## Reportar errores

Cuando un comando falla de forma inesperada, Olla-DFT registra la incidencia
en tu máquina —comando exacto, traza, versiones de Python y de las
dependencias, si QE estaba disponible— e imprime su id. Empaquétalo todo con

```bash
olla-dft report --export incidencias.json
```

y adjunta ese archivo al issue en
https://github.com/jorgegonzalezsevilla/olla-dft-esp/issues, junto con la
estructura si no es confidencial (`--attach archivo.cif` la copia al
registro). `olla-dft report "qué pasó"` anota algo que no reventó pero
confundió; `olla-dft report --stats` dice qué comandos fallan más. Nunca se
manda nada solo: no hay telemetría, el registro vive en tu carpeta de
configuración y tú decides si lo compartes.

Si el problema es un número equivocado y no un fallo, di contra qué referencia
comparaste y de dónde sale; la salida de `olla-dft crosscheck` y de `olla-dft
selftest` ayuda.

## Licencia

Olla-DFT es software libre bajo la GNU General Public License, versión 3 (ver
[LICENSE](LICENSE)). Al contribuir aceptas que tu contribución se distribuya
bajo la misma licencia. Los componentes y datos de terceros están en
[THIRD_PARTY_NOTICES.md](THIRD_PARTY_NOTICES.md).

Copyright © 2026 Jorge Enrique González Sevilla.
