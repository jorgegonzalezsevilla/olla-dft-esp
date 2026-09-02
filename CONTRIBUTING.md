# Contribuir a Olla-DFT

## Cómo se desarrolla este proyecto

Olla-DFT lo escribe y mantiene una sola persona, Jorge Enrique González
Sevilla, y así seguirá. **No se aceptan pull requests**: cualquier PR abierto
en este repositorio se cerrará sin revisar. No es un juicio sobre tu trabajo;
es simplemente cómo se lleva el proyecto.

Lo que sí es muy bienvenido es tu retroalimentación. Si encuentras un fallo,
obtienes un número equivocado o te gustaría una función, **abre un issue** en
https://github.com/jorgegonzalezsevilla/olla-dft-esp/issues y descríbelo. El
autor lee todos los issues y decide qué entra en la siguiente versión.

## Reportar errores

Cuando un comando falla de forma inesperada, Olla-DFT registra la incidencia
en tu máquina (comando exacto, traza, versiones de Python y de las
dependencias, si QE estaba disponible) e imprime su id. Empaquétalo todo con

```bash
olla-dft report --export incidencias.json
```

y adjunta ese archivo al issue, junto con la estructura si no es confidencial
(`--attach archivo.cif` la copia al registro). `olla-dft report "qué pasó"`
anota algo que no reventó pero confundió; `olla-dft report --stats` dice qué
comandos fallan más. Nunca se manda nada solo: no hay telemetría, el registro
vive en tu carpeta de configuración y tú decides si lo compartes.

Si el problema es un número equivocado y no un fallo, di contra qué referencia
comparaste y de dónde sale; la salida de `olla-dft crosscheck` y de `olla-dft
selftest` ayuda.

## Pedir funciones

Abre un issue con:

- qué intentas calcular o automatizar, y por qué;
- qué hace Olla-DFT hoy y en qué se queda corto;
- si es posible, un ejemplo pequeño (archivo de estructura, línea de
  comandos, resultado esperado).

Los casos de uso claros son los que más a menudo llegan a una versión.

## Correr la suite de pruebas por tu cuenta

Eres libre de clonar, leer y correr el código bajo la GPL. Si quieres
comprobar que todo funciona en tu máquina:

```bash
git clone https://github.com/jorgegonzalezsevilla/olla-dft-esp
cd olla-dft-esp
python -m venv .venv && source .venv/bin/activate   # opcional
pip install -e ".[test]"
python -m pytest -q                 # sin QE, menos de un minuto
```

Quantum ESPRESSO no hace falta para correr la suite; `olla-dft sistema` dice
qué tiene tu máquina. Los requisitos y las notas por sistema están en
[docs/PLATAFORMAS.md](docs/PLATAFORMAS.md).

## Licencia

Olla-DFT es software libre bajo la GNU General Public License, versión 3 (ver
[LICENSE](LICENSE)). Los componentes y datos de terceros están en
[THIRD_PARTY_NOTICES.md](THIRD_PARTY_NOTICES.md).

Copyright © 2026 Jorge Enrique González Sevilla.
