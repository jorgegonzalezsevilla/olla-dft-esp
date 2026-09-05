# Recuperar QE después de un corte

El modo `resilient` ejecuta QE por tramos y guarda generaciones completas en
un disco persistente. Ante una interrupción, vuelve al último guardado válido.
Puede perder el trabajo realizado desde ese guardado. No recupera un disco
borrado. La validación publicada simula interrupciones de procesos locales.

Inicializa el trabajo una sola vez con su input original:

```sh
olla-dft resilient init scf.in --state /mnt/olla/calculo-001 \
  --pw-cmd '/opt/qe/bin/pw.x' --checkpoint-seconds 900 \
  --runtime-id IMAGEN_INMUTABLE
olla-dft resilient run /mnt/olla/calculo-001
```

Después de reiniciar la máquina, ejecuta el mismo `run`. Conserva exactamente
el mismo QE, bibliotecas, imagen, pseudopotenciales, hilos y paralelismo MPI.
Olla rechaza cambios detectados en el entorno o en los archivos originales.
No uses este directorio para otros cálculos ni edites sus guardados.

```sh
olla-dft resilient status /mnt/olla/calculo-001
olla-dft resilient pause /mnt/olla/calculo-001
olla-dft resilient run /mnt/olla/calculo-001 --resume
```

Una pausa manual requiere `--resume`. Un corte o SIGTERM no crea esa pausa
persistente. El código 0 significa convergencia verificada; 75 significa que
el trabajo se detuvo o alcanzó el límite solicitado, y 2 requiere revisar el
error. `JOB DONE` por sí solo no demuestra convergencia.

## Arrancar automáticamente al volver la máquina

Con Olla instalado en la imagen definitiva, genera un servicio para el usuario
que posee el disco; después instálalo en esa máquina:

```sh
olla-dft resilient service /mnt/olla/calculo-001 --user olla \
  --output olla-calculo.service
sudo install -m 644 olla-calculo.service /etc/systemd/system/
sudo systemctl daemon-reload
sudo systemctl enable --now olla-calculo.service
```

El servicio requiere que el disco esté montado. En un apagado intenta detener
QE limpiamente, pero si el aviso no alcanza sigue existiendo el guardado
anterior. El servicio se ejecuta en la máquina local; no puede encender un
equipo apagado.

Conserva el directorio de estado en un disco local con montaje estable y un
único escritor. Mantén el mismo entorno y las mismas rutas. El sistema de
archivos debe soportar bloqueo, sincronización y renombrado atómico; consulta
el contrato local. Reserva espacio para las dos generaciones, el intento activo
y la nueva copia durante su publicación, además de logs y margen.

## Alcance y validación

La implementación contempla pw.x SCF, relax y vc-relax; se probaron los tres
modos con casos pequeños de silicio y QE 7.4, un proceso y un hilo. Otros motores y
modos (NSCF, bandas, fonones, MD o NEB) requieren protocolos específicos y se
rechazan. Valida tus materiales, tamaño de cálculo y configuración MPI antes de producción.
No se cambian malla k, cortes, tolerancias ni parámetros físicos para reducir
coste. Sí se administran las rutas, el modo de reinicio y el límite por tramo.

El intervalo de 900 segundos es un punto de partida, no una optimización
universal. Un tramo demasiado corto puede gastar todo su tiempo iniciando QE.
Mide el tiempo de copiar, guardar, restaurar y repetir trabajo perdido para
elegir un intervalo adecuado para tu equipo y cálculo.

Fuentes: [reinicio de QE](https://www.quantum-espresso.org/Doc/pw_user_guide/node20.html),
[parámetros de pw.x](https://www.quantum-espresso.org/Doc/INPUT_PW.html).

Una pausa recibida durante la restauración se respeta antes de iniciar QE y
no suma un intento. Si falla la escritura del registro de arranque, Olla detiene
y recoge el proceso hijo. Los estados terminados no mantienen un PID activo;
después de un apagón, un estado `running` todavía puede ser histórico. No uses
ese número por sí solo para matar procesos: el bloqueo del trabajo controla
si corresponde recuperar.
