# Olla-DFT: español integrado en la aplicación principal

**Este repositorio queda archivado.** Desde la versión 1.4.0, inglés y español
forman parte de una sola instalación de
[Olla-DFT](https://github.com/jorgegonzalezsevilla/olla-dft).

Abre `olla-dft` y elige **English** o **Español**. La aplicación recuerda tu
preferencia; puedes cambiarla con **l** en el menú o iniciar directamente con
`olla-dft --language es`.

## Migrar

Activa el entorno de Python donde usas Olla-DFT y ejecuta:

```bash
python -m pip install --upgrade "olla-dft @ git+https://github.com/jorgegonzalezsevilla/olla-dft.git@v1.4.0"
olla-dft config set language es
olla-dft --version
olla-dft
```

No necesitas dos paquetes. Se conserva la carpeta de configuración; tus
proyectos y el clon anterior no se borran. Si tienes una instalación editable
o cambios locales, consulta primero la
[guía de migración](https://github.com/jorgegonzalezsevilla/olla-dft/blob/main/docs/LANGUAGES.md#migrar-desde-la-edición-española).
No actualices el entorno de un cálculo activo como parte de esta migración.

- [Documentación en español](https://github.com/jorgegonzalezsevilla/olla-dft/blob/main/README.es.md)
- [Versiones y descargas](https://github.com/jorgegonzalezsevilla/olla-dft/releases)
- [Incidencias e ideas](https://github.com/jorgegonzalezsevilla/olla-dft/issues)

El historial y las versiones anteriores de este repositorio se conservan
como referencia. No se publicarán nuevas versiones aquí. Ambos idiomas
se citan como una sola aplicación en el
[registro principal de Zenodo](https://doi.org/10.5281/zenodo.22263121);
los DOI anteriores continúan disponibles.
