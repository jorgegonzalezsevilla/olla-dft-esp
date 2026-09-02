# Cambios

Todos los cambios relevantes de Olla-DFT. Las fechas son ISO 8601.

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
