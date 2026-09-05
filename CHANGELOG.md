# Cambios

Todos los cambios relevantes de Olla-DFT. Las fechas son ISO 8601.

## 1.3.1 — 2026-09-04

- Cargar Matplotlib y buscar fuentes solo al dibujar figuras; los comandos de estructuras y la ayuda conservan sus opciones sin inicializar las gráficas.
- Separar el ajuste de ecuaciones de estado de las dependencias de preparación de estructuras y dibujo.
- Verificar en procesos nuevos los comandos en español e inglés, los tres ajustes EOS y la exportación de figuras tras diferir las importaciones.
- Conservar fórmulas, tolerancias, entradas generadas y competidores del benchmark.

## 1.3.0 — 2026-09-04

- El software propio pasa a AGPL-3.0-or-later; las versiones anteriores conservan sus condiciones GPL originales.
- Sincronizar avisos, metadatos del paquete, citas, GitHub y Zenodo; conservar las licencias de terceros y de los ejemplos científicos existentes.
- Incluir la licencia completa del programa y un enlace al código de la versión en el HTML interactivo exportado.
- Aclarar en el inventario el pseudopotencial con hueco de core generado por el proyecto para pruebas.
- Sin cambios en algoritmos científicos, parámetros de cálculo ni resultados del benchmark.

## 1.2.0 — 2026-09-04

- Añade un explorador de resultados sin conexión con ejes numéricos, filtros, selección y figuras personalizables.
- Exporta SVG, PNG, CSV, JSON y HTML interactivo según la selección; conserva unidades, precisión e incertidumbre en los datos.
- Sustituye la serie de energías unidas del dashboard por el explorador interactivo.
- Registra malla k y parámetros, advierte sobre métodos mezclados y declara límites de instantáneas y exportaciones.

- Respeta pausas durante restauración, recoge QE si falla el registro del PID y elimina PIDs de estados terminados.

- Añade `resilient init/run/status/pause/service` para cálculos pw.x recuperables.
- Conserva dos guardados completos verificados y restaura un espacio privado tras cada corte.
- Congela input, UPFs, MPI, hilos, bibliotecas y arquitectura; limita fallos consecutivos.
- Genera un servicio Linux para continuar al arrancar con el disco persistente montado.
- Mide cálculo, restauración y copia para comparar costes con Olla-Lungo sin cambiar la física.
- Recuperación SCF, relax y vc-relax validada localmente con QE 7.4; recuperación tras apagones físicos o pérdida del disco aún sin medir.

## 1.1.1 — 2026-09-04

- El análisis de gap distingue bandas insuficientes de metalicidad.
- Advierte si el cálculo no convergió y devuelve estado de error cuando el gap no está disponible o validado.
- Rechaza mallas y cutoffs no positivos y parámetros no finitos antes de escribir inputs.
- Añade regresiones de la auditoría de septiembre.

## 1.1.0 — 2026-09-03

Cambios motivados por las primeras corridas de la comparativa (olla-dft-bench):

- Arrancar cuesta unas 8 veces menos: `import qekit.cli` pasó de ~0.6 s a
  ~0.07 s. seekpath, ase.io, matplotlib, strain y defects se importan al
  primer uso y no en cada invocación.
- `gen --kgrid N N N`: malla k explícita para scf/relax (anula `--kspacing` y
  `--klevel`). Hasta ahora solo se podía dar un espaciado.
- `mixing_beta` es 0.7 (el valor por omisión de QE) con ocupaciones fijas y
  sigue en 0.4 con smearing. En la celda de Si de la comparativa el scf pasa
  de 14 a unas 7 iteraciones con la misma energía.

## 1.0.1 — 2026-09-03

- Nuevo `olla-dft update` (alias `actualizar`): consulta la última versión
  publicada en GitHub, muestra las novedades y los comandos exactos que
  ejecutaría, e instala solo tras confirmar (`--check` para solo mirar,
  `--yes` para no preguntar, `--version TAG` para una versión concreta).
  Funciona con instalaciones hechas con pip desde GitHub y con clones
  locales. Olla-DFT nunca busca actualizaciones por su cuenta.
- Enlace a la comparativa reproducible frente a ASE, pymatgen y seekpath
  (olla-dft-bench).

## 1.0.0 — 2026-09-02

Primera versión pública (GPL-3.0).

- Comando único `olla-dft`; el paquete Python conserva el nombre `qekit` para
  que sigan funcionando los guiones y las carpetas de proyecto (`.qekit/`).
  La configuración y los datos van ahora a una carpeta `olla-dft`
  (`~/.config/olla-dft` en Linux), migrada sola desde las carpetas
  `qekit`/`QEkit` anteriores.
- Interfaz bilingüe (español por defecto, inglés con `--language en`,
  `OLLA_DFT_LANG` o `config set language en`): ayuda de los 78 comandos y
  1 300 opciones, menú interactivo, inicio guiado, recetas, wizard, dashboard
  y referencia HTML. Alias en inglés `recipes`, `theory`, `system`.
- Nuevo `olla-dft teoria` / `theory`: el fundamento físico de cada comando
  científico (qué responde, fórmulas implementadas, procedimiento con la
  función y el binario de QE responsables, de dónde sale cada dato, límites y
  referencias), publicado también como `docs/TEORIA.md`.
- Ayuda agrupada en *opciones* / *ejecución* / *parámetros DFT* / *figura*;
  todas las opciones tienen ya texto de ayuda.
- Retirada la capa de plataforma experimental que no aportaba física:
  asistente LLM local, servidor web, plugins, monitor/envío HPC, preflight de
  release/SBOM y conectores a bases de datos externas.
- La auditoría de las fórmulas contra el código corrigió más de 30 fallos,
  entre ellos: cargas de Bader integradas con la unidad de volumen
  equivocada (≈6.75× de más); alineamiento del potencial de defectos sumado
  en Ry en vez de eV; `adsorb --dipole` no escribía la corrección dipolar;
  `surface --fix` no llegaba a `ATOMIC_POSITIONS`; `tddft --compare` leía la
  reflectividad en vez de la absorción; `--raman` sin `--gamma` reventaba
  tras toda la cadena DFPT; `kappa` forzaba ocupaciones fijas en metales;
  parámetro de red de QHA mal en celdas primitivas de 1 átomo; el
  desdoblamiento mezclaba canales de espín en corridas lsda; el factor de
  forma f₂ de Allen–Dynes nunca se aplicaba; la figura de alineamiento
  dibujaba mal el CBM; los marcadores plegados de la fase de Berry usaban el
  módulo equivocado; `gen --soc` no verificaba los pseudos relativistas;
  `doctor` mezclaba ciclos SCF de una relajación; la planitud de la función
  trabajo se medía fuera del vacío; `transport --spin-resolved` y
  `--kspacing` no estaban conectados; la columna Tc de `elph` siempre vacía;
  XANES aceptaba bordes que xspectra.x no calcula; y varios mensajes
  engañosos.
- Banderas nuevas: `hubbard --hubbard-style`, `unfold --spin`,
  `kappa --metal`, `qha --structure`, `transport --nspin/--mag`,
  `tddft --scissor`, `charges --pseudo-dir`; se retira `ballistic --ikind 2`
  (nunca estuvo implementado).
- Ejemplos renombrados por tema con README bilingüe cuyos comandos valida una
  prueba contra el parser; pruebas renombradas por tema; 977 pruebas.
- Archivos de repositorio para GitHub: licencia GPL-3.0, avisos de terceros,
  CONTRIBUTING, CITATION.cff, CI en Python 3.9–3.13, `.gitignore`.

## Antes de la 1.0

El proyecto creció en privado como *QEkit* (0.1–0.34) y *Olla-DFT* (0.35). Hitos, a título de referencia:

| Versión | Añadido |
|---|---|
| 0.1–0.4 | generación de inputs de pw.x, herramientas de estructura, post-proceso de bandas/DOS/gap, estilos y plantillas de publicación |
| 0.5–0.7 | convergencia, EOS, constantes elásticas; materiales laminares (capas, DRX, exfoliación); ópticas, densidad de carga, función trabajo, fonones DFPT |
| 0.8–0.12 | masa efectiva, procedencia, Raman, XPS, transporte, Bader/Löwdin, superficies y defectos, SOC y DFT+U, doctor/audit/db/hull, MLIP, registro de incidencias, crosscheck, derived, QHA, datasheet |
| 0.13–0.14 | pseudopotenciales con hueco de core, XANES, U de Hubbard autoconsistente, electrón-fonón, desdoblamiento, NEB, termoquímica, análisis de MD, interfaces, wizard, selección de pseudopotenciales, TDDFPT, transporte balístico |
| 0.15–0.20 | deformación, adsorción, elásticas 2D, corrección dipolar, defectos cargados, barridos en paralelo con presupuesto de tiempo, estimador de coste, energía de superficie, alineamiento de bandas, fatbands, Hubbard V, fonones a temperatura electrónica, híbridos, número de Lorenz, transporte por espín, selftest, CHE, referencia HTML |
| 0.21–0.28 | sólidos amorfos, funciones de Wannier y desenredado, fase de Berry, conductividad térmica de red, ESM, recetas, portabilidad (Linux/macOS/Windows, salida ASCII) |
| 0.29–0.35 | topología, workflow de proyecto, base de resultados, campañas, dashboard, inicio guiado bilingüe, cambio de nombre a Olla-DFT |
