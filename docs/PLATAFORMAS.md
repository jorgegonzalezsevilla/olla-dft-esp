# Plataformas y requisitos

## Dónde corre

Olla-DFT es Python puro: corre igual en Linux, en macOS y en Windows. Lo que no
es portable es Quantum ESPRESSO, y eso cambia el consejo según el sistema.
Antes de nada:

```
olla-dft sistema
```

(también vale `olla-dft system`) dice qué ve Olla-DFT de tu máquina
—codificación de la consola, dónde guarda la configuración, qué binarios de
QE encuentra y dónde, si hay un lanzador MPI, cuántos núcleos ve— y qué hacer
con lo que falte. Es el primer comando que hay que correr en un sistema nuevo.

Olla-DFT sirve sin un solo binario de QE en la máquina: genera los inputs y
post-procesa salidas traídas de otro sitio. Solo `--run` necesita los binarios
aquí.

## Instalar Quantum ESPRESSO en cada sistema

**Linux.** Quantum ESPRESSO está empaquetado en todas las distribuciones
grandes: `sudo apt install quantum-espresso` en Debian y Ubuntu, `sudo dnf
install quantum-espresso` en Fedora y RHEL, `sudo zypper install
quantum-espresso` en openSUSE, `sudo pacman -S quantum-espresso` en Arch, o
`conda install -c conda-forge qe`. Compilarlo desde
https://www.quantum-espresso.org da mejor rendimiento. Si lo compilaste a mano
y no está en el PATH: `olla-dft config set pw_cmd /ruta/a/bin/pw.x`. Olla-DFT
mira además en `/usr/bin`, `/usr/local/bin`, `/opt/qe/bin`, `~/q-e/bin` y
`/usr/lib64/openmpi/bin` sin configurar nada.

**macOS.** `brew install quantum-espresso` o `sudo port install
quantum-espresso`. En Apple Silicon Homebrew instala en `/opt/homebrew/bin`,
que no siempre está en el PATH de un shell no interactivo; Olla-DFT mira ahí
de todos modos (y en `/usr/local/bin`, en `/opt/local/bin` de MacPorts y en
`~/q-e/bin`). Los guiones generados cuentan núcleos con `sysctl -n hw.ncpu`
además de con `nproc`, que en macOS no existe.

**Windows.** Hay tres caminos, de más a menos recomendable:

1. **WSL2.** `wsl --install`, y dentro de Ubuntu todo se comporta como en
   Linux (`sudo apt install quantum-espresso python3-pip` y luego Olla-DFT
   con pip). Es lo mejor probado y lo que menos sorpresas da.
2. **Binarios nativos de QE.** Se llaman `pw.exe`, no `pw.x`; Olla-DFT prueba
   las dos terminaciones y los busca también en `C:\Program Files\QE\bin`,
   `C:\Program Files\quantum-espresso\bin` y `C:\qe\bin`. Si están en otro
   sitio: `olla-dft config set pw_cmd C:\ruta\a\pw.exe`. Los barridos se
   lanzan con `python run.py` en vez de `run.sh`: **cada carpeta de trabajo
   lleva los dos guiones**, el `.sh` para POSIX y el `.py`, que hace lo mismo
   sin necesitar bash ni xargs y detecta solo si hay `mpiexec`.
3. **Olla-DFT aquí, QE en un clúster.** Genera los inputs en el portátil,
   córrelos donde sea y trae las salidas: todo el post-proceso funciona sin un
   solo binario de QE en tu máquina.

Los guiones `.sh` generados se escriben siempre con finales de línea POSIX,
también desde Windows, para que no fallen con `bad interpreter: /bin/bash^M`
al correrlos en WSL o en un clúster.

## MPI

El lanzador se detecta en este orden: `mpirun -np N` (OpenMPI), `mpiexec -n N`
(MPICH, MS-MPI), `srun -n N` (Slurm). Si no hay ninguno, los cálculos van en
serie, que en un portátil suele ser lo correcto; ahí, correr varios puntos a la
vez con `-j N` aprovecha mejor 2 o 4 núcleos que un solo `pw.x` con MPI.

## La codificación de la consola y `--ascii`

Los informes llevan Å, α, ε, →, ①, ✓, y la página de códigos heredada de
Windows (cp1252) no puede escribirlos: `print` lanza `UnicodeEncodeError` y el
comando **muere a media salida** con código 1, así que parece que falló el
cálculo. Se comprobó: `info`, `selftest` y `recetas` reventaban las tres.
Olla-DFT lo resuelve en dos escalones — primero intenta poner la consola en
UTF-8 (Windows 10 en adelante lo admite), y si no puede translitera a ASCII
conservando el significado:

```
  │ a0 = 5.402 Å   κ_L → 100.7 W/m·K   ✓          (UTF-8)
  | a0 = 5.402 A   kappa_L -> 100.7 W/m.K   ok    (--ascii)
```

Nunca sustituye por `?`: un informe que dice «el parámetro de red es 5.43 ?»
es peor que no imprimir nada, porque parece un dato corrupto. Con `--ascii` se
fuerza en cualquier sistema, que es lo que hay que usar al redirigir a un
archivo o al pegar la salida en un correo; como `--language`, se acepta en
cualquier posición de la línea de comandos. En Windows también se puede
activar UTF-8 con `chcp 65001` o poniendo `PYTHONUTF8=1`.

Hay pruebas para las tres plataformas: simulan `sys.platform`, lanzan el
programa entero con la salida forzada a cp1252, y una recorre **todos los
informes que Olla-DFT genera de verdad** comprobando que cada carácter fuera de
ASCII se puede transliterar. Añadir un símbolo nuevo a un informe sin enseñarle
a la tabla rompe pytest antes de que lo sufra nadie — y así se encontraron
cuatro (`Λ`, `⇌`, `∝` y el menos tipográfico `−`, que parece un guion y no lo
es) y un `mpirun` cableado en el módulo de XPS.

## Carpetas de configuración y de datos

La configuración (`config.ini`, plantillas) sigue la convención de cada
sistema; los datos grandes y reemplazables (modelos, histórico) van a una
carpeta aparte para no llenar copias de seguridad ni romper las convenciones
de cada escritorio.

| | Configuración | Datos |
|---|---|---|
| Linux | `~/.config/olla-dft` (respeta `XDG_CONFIG_HOME`) | `~/.local/share/olla-dft` (respeta `XDG_DATA_HOME`) |
| macOS | `~/Library/Application Support/olla-dft` | `~/Library/Application Support/olla-dft` |
| Windows | `%APPDATA%\olla-dft` | `%LOCALAPPDATA%\olla-dft` |

Las variables de entorno mandan sobre los valores por omisión:

| Variable | Efecto |
|---|---|
| `OLLA_DFT_CONFIG_DIR` | carpeta de configuración (correr desde un lápiz USB, o en un clúster con el HOME lleno; las pruebas la usan para no tocar nunca la configuración de verdad) |
| `OLLA_DFT_DATA_DIR` | carpeta de datos (poner modelos e histórico en un disco con más espacio) |
| `OLLA_DFT_LANG` | idioma de la interfaz, `es` o `en` (por debajo de `--language`, por encima de `olla-dft config set language`) |

Los nombres antiguos `QEKIT_CONFIG_DIR` y `QEKIT_DATA_DIR` se siguen
respetando. Una configuración dejada por una versión anterior —en
`~/.config/qekit` o en la carpeta `QEkit` de cada sistema— se copia sola la
primera vez que corre Olla-DFT, y `olla-dft sistema` lo dice cuando encuentra
una.

El registro local de incidencias (`olla-dft report`) vive dentro de la carpeta
de configuración. Nunca se manda nada a ningún sitio: no hay telemetría.

## Requisitos

- Python ≥ 3.9 con numpy ≥ 1.20, scipy ≥ 1.8, matplotlib ≥ 3.5, ASE ≥ 3.22,
  spglib ≥ 2.0 y seekpath ≥ 2.0 (se instalan solos con `pip install .` o
  `pip install -e .`).
- Quantum ESPRESSO para correr los cálculos: `pw.x`, `dos.x`, `projwfc.x`,
  `bands.x`, `pp.x`, `epsilon.x`, `ph.x`, `q2r.x`, `matdyn.x`, `dynmat.x` y
  `pw2wannier90.x` para `wannier`. Todos forman parte de una compilación
  normal de QE.
- Algunos módulos necesitan binarios que QE **no compila por omisión**. Desde
  la carpeta del código fuente de QE:

  ```bash
  make ld1        # olla-dft corehole  (pseudos con hueco de core)
  make xspectra   # olla-dft xanes
  make hp         # olla-dft hubbard
  make neb        # olla-dft neb
  make tddfpt     # olla-dft tddft     (turbo_lanczos.x, turbo_spectrum.x)
  make pwcond     # olla-dft ballistic
  ```

  Olla-DFT avisa con el comando exacto cuando falta alguno, en vez de fallar
  con un "command not found".
- Extras opcionales:

  | Instalar | Añade | Lo usan |
  |---|---|---|
  | `pip install "olla-dft[mlip]"` | torch, mace-torch (unos 1.2 GB) | `olla-dft mlip`, `olla-dft amorphous`, fuerzas MACE en `kappa`, `selftest --mlip` |
  | `pip install "olla-dft[kappa]"` | phono3py ≥ 3 | `olla-dft kappa` (ecuación de Boltzmann de fonones) |
  | `pip install "olla-dft[test]"` | pytest ≥ 7, pyflakes ≥ 3 | la suite de pruebas |
