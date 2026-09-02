#!/bin/bash
# ----------------------------------------------------------------------
# Barrido de regresión de la interfaz de línea de comandos de olla-dft.
#
# Las pruebas de pytest miran las FUNCIONES; esto mira los COMANDOS, que
# es donde se cuelan los fallos que el usuario sí ve: una bandera que no
# llega al módulo, un archivo que se escribe a medias, un mensaje que
# aparece cuando no toca. La 0.12.0 salió de aquí con tres arreglos que
# ninguna prueba unitaria había visto.
#
# Cada línea declara el código de salida que ESPERA:
#     0  todo bien
#     1  el comando corrió y encontró un problema (doctor, crosscheck)
#     2  error de uso: la bandera o el dato no encajan
#
# Uso:
#     OLLA_DFT_REG=/ruta/con/datos bash tests/barrido_cli.sh
#
# Necesita salidas de QE ya calculadas; no las trae el repositorio porque
# pesan. Ajusta las rutas de arriba a las tuyas.
# ----------------------------------------------------------------------
REG="${OLLA_DFT_REG:-/tmp/reg012}"
cd "$REG" || { echo "no existe $REG; define OLLA_DFT_REG"; exit 1; }
# Equivale al comando `olla-dft`, pero convirtiendo los RuntimeWarning de
# numpy en errores para que un NaN silencioso no pase por OK.
Q="python3 -W error::RuntimeWarning -m qekit.cli"
LOG="$REG/barrido.log"
: > "$LOG"
ok=0; bad=0
run() {
  esp=0
  case "$1" in --rc*) esp="${1#--rc}"; shift;; esac
  desc="$1"; shift
  echo "### $desc :: $*" >> "$LOG"
  out=$($Q "$@" 2>&1); rc=$?
  echo "$out" >> "$LOG"
  echo "--- rc=$rc" >> "$LOG"; echo >> "$LOG"
  if [ "$rc" = "$esp" ]; then ok=$((ok+1)); printf 'OK    rc=%d  %s\n' "$rc" "$desc";
  else bad=$((bad+1)); printf 'FALLA rc=%d (esperaba %s)  %s\n' "$rc" "$esp" "$desc";
       echo "$out" | tail -6 | sed 's/^/      /'; fi
}

echo "=== ESTRUCTURA ==="
run "info Si"             info Si.cif
run "info hBN"            info hbn.cif
run "kpath"               kpath Si.cif
run "prim"                prim Si.cif -o o/prim.cif
run "conv"                conv Si.cif -o o/conv.cif
run "supercell 2x2x2"     supercell Si.cif 2 2 2 -o o/sc.cif
run "convert cif->POSCAR" convert Si.cif o/POSCAR_Si
run "convert cif->xyz"    convert Si.cif o/Si.xyz
run "surface (111)"       surface Si.cif -m "1 1 1" -l 4 --vacuum 15 -o o/slab.cif
run "defect vacancia"     defect Si.cif -k vacancy --site 0 --supercell 2x2x2 -o o/def
run "layers hBN"          layers hbn.cif
run "layers + slab"       layers hbn.cif --slab o/mono.cif

echo "=== GENERACION DE INPUTS ==="
run "gen scf"             gen Si.cif -p scf -o o/gen
run "gen all"             gen Si.cif -p all -o o/genall
run "gen relax insulator" gen Si.cif -p relax --insulator -k fine -o o/genrel
run "gen nspin=2 Fe"      gen Fe.cif -p scf --nspin 2 --mag 0.3 -o o/genfe
run "templates list"      templates list
run "templates show"      templates show journal
run "templates export"    templates export journal -o o/tpl.json
run "config show"         config show

echo "=== BANDAS / DOS ==="
run "bands"               bands bandas_si --no-plot -o o/bandas
run "bands figura"        bands bandas_si -o o/bandasfig --journal aps --format png
run "gap"                 gap bandas_si
run "dos"                 dos /tmp/test_si --no-plot -o o/dos
run "dos elemento"        dos /tmp/test_si --mode element --no-plot -o o/dosel
run "plot bandas+dos"     plot /tmp/test_si -o o/plot --format png
run "dos Fe spin"         dos /tmp/test_fe --no-plot -o o/dosfe

echo "=== PROPIEDADES ==="
run "xrd"                 xrd Si.cif -o o/xrd --no-plot
run "xrd suite"           xrd Si.cif -o o/xrds --suite --no-plot
run "eos collect"         eos Si.cif --collect -o /tmp/val/eos --no-plot
run "elastic collect"     elastic Si.cif --collect -o /tmp/val/elastic --no-plot
run "converge collect"    converge Si.cif --collect -o /tmp/val/conv_ecut --no-plot
run "effmass collect"     effmass Si.cif --collect -o /tmp/mef_cli
run "optics collect"      optics Si.cif --collect -o /tmp/opt_si --no-plot
run "optics tauc indir"   optics Si.cif --collect -o /tmp/opt_si --tauc indirect --no-plot
run "phonons collect"     phonons Si.cif --collect -o /tmp/fon_si --no-plot
run "phonons raman"       phonons Si.cif --collect -o /tmp/raman_si --gamma --raman --no-plot
run "transport collect"   transport Si.cif --collect -o /tmp/tr_si --no-plot
run "wf"                  wf /tmp/xps_si --no-plot -o o/wf
run "charges lowdin"      charges --lowdin /tmp/test_si/projwfc.out -o o/chg --no-plot
run "fermi"               fermi -o /tmp/tr_si
run "xps collect"         xps Si.cif --collect -o /tmp/xps_si

echo "=== DATOS Y CONTROL DE CALIDAD ==="
run --rc1 "doctor"        doctor /tmp/chk -o o/doc --no-plot
run "audit"               audit /tmp/chk /tmp/tr_si /tmp/raman_si
run "audit --index"       audit /tmp/chk /tmp/tr_si --index --db o/reg.db
run "db registrar"        db /tmp/chk /tmp/raman_si --db o/reg.db
run "db consulta"         db --db o/reg.db -q "SELECT formula, ecutwfc FROM calculos"
run "db export"           db --db o/reg.db --export o/db.json
run "hull"                hull /tmp/raman_si -o o/hull --no-plot
run "suggest"             suggest Si.cif --db o/reg.db
run "report stats"        report --stats
run "report listar"       report --only-open

echo "=== NUEVOS EN 0.12.0 ==="
run "derived Cij"         derived Si.cif --cij ELASTIC_C.dat
run "derived 500 K"       derived Si.cif --cij ELASTIC_C.dat --temp 500
run "qha"                 qha qha_tabla.dat -o o/qha --cells 8 --natoms 2 --cubic --no-plot
run "qha figura"          qha qha_tabla.dat -o o/qhafig --cells 8 --natoms 2 --cubic --format png
run --rc1 "crosscheck"    crosscheck /tmp/proy -f Si.cif
run --rc1 "crosscheck gaps" crosscheck /tmp/proy -f Si.cif --gap-bandas 0.61 --gap-tauc 2.56
run "datasheet"           datasheet /tmp/proy -o o/ficha --name Si_reg
run "datasheet metodos"   datasheet /tmp/proy --methods

echo "=== ERRORES DE USO (deben salir con codigo 2 y sin traza) ==="
run --rc2 "malla incompleta"   transport Si.cif --collect -o /tmp/tr_si --grid 1x2 --no-plot
run --rc2 "malla no entera"    transport Si.cif --collect -o /tmp/tr_si --grid 8x8xocho --no-plot
run --rc2 "plantilla=revista"  xrd Si.cif -o o/xrdt -t nature --format png
run --rc2 "columna inexistente" db --db o/reg.db -q "SELECT prefix FROM calculos"
run --rc2 "ecuacion desconocida" eos Si.cif --collect -o /tmp/val/eos --equation birch --no-plot

echo "=== NUEVOS EN 0.13.0 ==="
run "corehole (solo inputs)" corehole Si --edge K -o o/ch --functional PZ --rcut 1.6 --only-inputs
run "corehole core-wfc"     corehole --core-wfc ps_ch/Si.hueco1s.UPF --orbital 1S --output o/Si.wfc
run "xps prepara"           xps Si.cif -o o/xps --core-hole Si=ps_ch/Si.hueco1s.UPF --pseudo-dir ps_ch --ecutwfc 40 --metal
run "xanes prepara"         xanes gr.cif -o o/xan --element C --core-hole ps_ch/Si.hueco1s.UPF --pseudo-dir ps_ch --ecutwfc 40 --average
run "xanes collect"         xanes Si.cif --collect -o xanes_si --element Si --no-plot
run "xanes figura"          xanes Si.cif --collect -o xanes_si --element Si --format png
run "hubbard prepara"       hubbard NiO.cif -o o/hub --pseudo-dir ps_ch --ecutwfc 45 --qgrid 2x2x2 --nspin 2 --mag Ni=0.5 --metal
run "hubbard collect"       hubbard NiO.cif --collect -o hub_nio --qgrid 2x2x2
run "elph prepara"          elph Al.cif -o o/elph --qgrid 2x2x2 --pseudo-dir ps_ch --ecutwfc 20
run "elph collect"          elph --collect -o elph_al --debye 428 --no-plot
run "md analiza"            md md_si/md.out -o o/md --skip 50 --no-plot
run "md figura"             md md_si/md.out -o o/mdfig --skip 50 --format png
run "interface lista"       interface gr.cif hbn_mono.cif --list --max-index 2 --max-atoms 40
run "interface construye"   interface gr.cif hbn_mono.cif -o o/itf --max-index 2 --max-atoms 40 --name gr_hbn
run "unfold"                unfold unf unf/prim.cif -o o/unf --bands 8 --no-plot
run "unfold figura"         unfold unf unf/prim.cif -o o/unffig --bands 8 --format png
run --rc1 "neb collect"     neb unf/prim.cif --collect -o neb_h3 --no-plot
run "thermochem gas"        thermochem "1595,3657,3756" --phase gas --structure h2o.xyz --symmetry 2 --energy -14.22
run "thermochem solido"     thermochem "120,340,520,610,800,1100" --phase solido --floor 100
run "wizard catalogo"       wizard --list
run "wizard meta"           wizard Si.cif --goal conduce
run "wizard pregunta"       wizard NiO.cif --ask "mi oxido sale metalico"
run "wizard termino"        wizard --term "malla k"

echo "=== ERRORES DE USO NUEVOS (codigo 2) ==="
run --rc2 "corehole elemento raro" corehole Xx --edge K -o o/x
run --rc2 "borde inexistente"      corehole H --edge L23 -o o/x
run --rc2 "xanes sin core-hole"    xanes Si.cif -o o/x --element Si
run --rc2 "wizard meta inexistente" wizard Si.cif --goal noexiste
run --rc2 "thermochem fase mala"   thermochem "500,900" --phase liquido
run --rc2 "polarizacion incompleta" xanes Si.cif -o o/x --element Si --core-hole ps_ch/Si.hueco1s.UPF --polarization "1 0"

echo "=== NUEVOS EN 0.14.0 ==="
run "pseudos de un elemento"  pseudos --element Si --pseudo-dir /usr/share/espresso/pseudo
run "pseudos de una estructura" pseudos NiO.cif --pseudo-dir /usr/share/espresso/pseudo
run "pseudos para fonones"    pseudos --element Ni --task fonones --pseudo-dir /usr/share/espresso/pseudo
run "pseudos para SOC"        pseudos --element Ni,O --task soc --pseudo-dir /usr/share/espresso/pseudo
run --rc1 "pseudos sin opcion" pseudos --element Ni --task optics --pseudo-dir /usr/share/espresso/pseudo
run "gen con --pseudo"        gen NiO.cif -p scf -o o/gp --pseudo-dir /usr/share/espresso/pseudo --pseudo Ni=Ni.pbe-nd-rrkjus.UPF --pseudo O=O.pbe-rrkjus.UPF
run "corehole --plain"        corehole Si --plain -o o/plain --functional PZ --rcut 1.6 --only-inputs
run "tddft lanczos prepara"   tddft Si.cif -o o/td --method lanczos --iter 200 --pseudo-dir ps_ch --ecutwfc 30
run "tddft davidson prepara"  tddft Si.cif -o o/td2 --method davidson --states 5 --pseudo-dir ps_ch --ecutwfc 30
run "tddft collect davidson"  tddft --collect -o tddft_c2h4 --method davidson --gap 6.6 --no-plot
run "tddft figura"            tddft --collect -o tddft_c2h4 --method davidson --gap 6.6 --format png
run "ballistic prepara"       ballistic hilo.cif -o o/bal --pseudo-dir ps_ch --ecutwfc 15 --emin -3 --emax 3 --points 21
run "ballistic collect"       ballistic --collect -o balistico_al --no-plot
run "ballistic figura"        ballistic --collect -o balistico_al --format png

echo "=== ERRORES DE USO (codigo 2) (codigo 2) ==="
run --rc2 "pseudo forzado inexistente" gen NiO.cif -p scf -o o/x --pseudo-dir /usr/share/espresso/pseudo --pseudo Ni=noexiste.UPF
run --rc2 "pseudo mal escrito"    gen NiO.cif -p scf -o o/x --pseudo-dir /usr/share/espresso/pseudo --pseudo NiUPF
run --rc2 "tarea inexistente"     pseudos --element Si --task magia --pseudo-dir /usr/share/espresso/pseudo
run --rc2 "tddft metodo malo"     tddft Si.cif -o o/x --method montecarlo
run --rc2 "extrapolacion mala"    tddft Si.cif -o o/x --extrapolation magia
run --rc2 "ballistic sin z"       ballistic Si.cif -o o/x --pseudo-dir ps_ch
run --rc2 "dipolo sin vacio"      gen Si.cif -p scf -o o/x --dipole
run --rc2 "soc con nspin 2"       gen Si.cif -p scf -o o/x --soc --nspin 2
run --rc2 "hubbard sin igual"     gen NiO.cif -p scf -o o/x --hubbard Ni:4.6
run --rc2 "hubbard no numerico"   gen NiO.cif -p scf -o o/x --hubbard Ni=mucho
run --rc2 "dt sin preset md"      gen Si.cif -p scf -o o/x --dt 0.5
run --rc2 "strain modo malo"      strain Si.cif -m diagonal -o o/x
run --rc2 "strain rango corto"    strain Si.cif -r 0:5:2 -o o/x
run --rc2 "strain rango al reves" strain Si.cif -r 5:-5:11 -o o/x
run --rc2 "strain rango enorme"   strain Si.cif -r -50:50:11 -o o/x
run --rc2 "strain perp con hidro" strain Si.cif -m hidrostatica --relax-perp -o o/x
run --rc2 "adsorb sobre bulto"    adsorb Si.cif --mol CO -o o/x
run --rc2 "adsorb molecula rara"  adsorb al_slab.cif --mol unobtanio -o o/x
run --rc2 "adsorb sitio raro"     adsorb al_slab.cif --mol H --sites esquina -o o/x
run --rc2 "adsorb ancla fuera"    adsorb al_slab.cif --mol CO --anchor 9 -o o/x
run --rc2 "elastic 2d sobre bulto" elastic Si.cif --2d -o o/x
run --rc2 "eform carga sin epsilon" eform Si.cif -k vacancy -q 1 --supercell 2x2x2 -o o/x
run --rc2 "eform correccion mala"   eform Si.cif -k vacancy -q 0 --supercell 2x2x2 --correction magia -o o/x
run --rc2 "eform cargas no enteras" eform Si.cif -k vacancy -q 0.5 --supercell 2x2x2 -o o/x
run --rc2 "eform sustitucion sin el" eform Si.cif -k substitution -q 0 --supercell 2x2x2 -o o/x
run --rc2 "eform mu mal escrito"     eform Si.cif -k vacancy -q 0 --supercell 2x2x2 --mu Si:-107 -o o/x
run --rc2 "max-time sin sentido"    strain Si.cif -r -2:2:3 -o o/x --max-time manana
run --rc2 "max-time negativo"       strain Si.cif -r -2:2:3 -o o/x --max-time 0
run --rc2 "jobs cero"               strain Si.cif -r -2:2:3 -o o/x -j 0
run --rc2 "nproc cero"              strain Si.cif -r -2:2:3 -o o/x --nproc 0
run --rc2 "timeout negativo"        strain Si.cif -r -2:2:3 -o o/x --timeout -5
run --rc2 "gamma un solo grosor" gamma Si.cif -m "1 1 1" -l 4 -o o/x
run --rc2 "gamma una capa"       gamma Si.cif -m "1 1 1" -l 1,4 -o o/x
run --rc2 "gamma miller corto"   gamma Si.cif -m "1 1" -l 3,4 -o o/x
run --rc2 "gamma miller no num"  gamma Si.cif -m "1 1 x" -l 3,4 -o o/x
run --rc2 "gamma capas no num"   gamma Si.cif -m "1 1 1" -l tres,4 -o o/x
run --rc2 "dband sin ese elemento" dos /tmp/test_si --dband Au
run --rc2 "dband orbital raro"     dos /tmp/test_si --dband Si-g
run --rc2 "fat selector raro"      bands fat_si --fat Ni-d -o o/x --no-plot
run --rc2 "fat de otro calculo"    bands bandas_si --fat Si-p --projwfc /tmp/test_si/projwfc.out -o o/x --no-plot
run --rc2 "align carpeta mala"   align /tmp/no_existe hbn_b -o o/x --no-plot
run --rc2 "tscan una sola T"     phonons Si.cif --tscan 300 --gamma -o o/x
run --rc2 "tscan T absurda"      phonons Si.cif --tscan 300,99000 --gamma -o o/x
run --rc2 "tscan no numerico"    phonons Si.cif --tscan 300,frio --gamma -o o/x
run --rc2 "EXX no divide la k"   gen Si.cif -p scf --functional hse --exx-grid 7x7x7 -o o/x
run --rc2 "malla EXX incompleta" gen Si.cif -p scf --functional hse --exx-grid 2x2 -o o/x
run --rc2 "eform mu no numerico"     eform Si.cif -k vacancy -q 0 --supercell 2x2x2 --mu Si=mucho -o o/x
run --rc2 "eform intersticial sin pos" eform Si.cif -k interstitial --new-element H -q 0 --supercell 2x2x2 -o o/x
run --rc2 "defect sustitucion sin el"  defect Si.cif -k substitution -o o/x

echo "=== BANDERAS DE FISICA EN gen (0.15.0) ==="
run "gen con vdw"             gen Si.cif -p scf -o o/g_vdw --vdw grimme-d3
run "gen con hubbard"         gen NiO.cif -p scf -o o/g_u --hubbard Ni=4.6
run "gen hubbard tarjeta"     gen NiO.cif -p scf -o o/g_u2 --hubbard Ni=4.6 --hubbard-style card
run "gen con carga"           gen Si.cif -p scf -o o/g_q --charge 1
run "gen con soc"             gen Si.cif -p scf -o o/g_soc --soc
run "gen con nosym"           gen Si.cif -p scf -o o/g_ns --nosym
run "gen dipolo en losa"      gen al_slab.cif -p scf -o o/g_dip --dipole
run "gen dipolo eje 3"        gen al_slab.cif -p scf -o o/g_dip3 --dipole 3
run "gen md"                  gen al_slab.cif -p md -o o/g_md --nstep 200 --dt 1.0
run "gen md con termostato"   gen al_slab.cif -p md -o o/g_md2 --thermostat berendsen -T 500

echo "=== BARRIDO DE DEFORMACION (0.15.0) ==="
run "strain biaxial"          strain hbn_mono.cif -m biaxial -r -3:3:5 -o o/st1
run "strain hidrostatica"     strain Si.cif -m hidrostatica -r -2:2:5 -o o/st2
run "strain uniaxial"         strain Si.cif -m uniaxial-a -r -2:2:5 -o o/st3
run "strain cizalla"          strain Si.cif -m cizalla -r -1:1:3 -o o/st4
run "strain iones fijos"      strain Si.cif -r -2:2:3 --fixed-ions -o o/st5
run "strain relax perp"       strain hbn_mono.cif -r -2:2:3 --relax-perp -o o/st6
run "strain con espin"        strain Fe.cif -r -2:2:3 --nspin 2 --mag 2.2 -o o/st7

echo "=== ADSORCION (0.15.0) ==="
run "adsorb CO"               adsorb al_slab.cif --mol CO --vdw grimme-d3 -o o/ad1
run "adsorb solo huecos"      adsorb al_slab.cif --mol H --sites hollow -o o/ad2
run "adsorb con rotaciones"   adsorb al_slab.cif --mol CO2 --rotations 3 -o o/ad3
run "adsorb con dipolo"       adsorb al_slab.cif --mol H --dipole -o o/ad4
run "adsorb ancla explicita"  adsorb al_slab.cif --mol CO --anchor 1 -o o/ad5

echo "=== ELASTICO 2D (0.15.0) ==="
run "elastic 2d"              elastic hbn_mono.cif --2d --npoints 2 -o o/el2d
run "elastic 2d con espesor"  elastic hbn_mono.cif --2d --thickness 3.33 --npoints 2 -o o/el2dt

echo "=== DEFECTOS CARGADOS (0.16.0) ==="
run "eform neutro"            eform Si.cif -k vacancy --site 0 -q 0 --supercell 2x2x2 -o o/ef1
run "eform con cargas"        eform Si.cif -k vacancy --site 0 -q -1,0,1 --supercell 2x2x2 --epsilon 11.7 -o o/ef2
run "eform makov-payne"       eform Si.cif -k vacancy -q -1,1 --supercell 2x2x2 --epsilon 11.7 --correction makov-payne -o o/ef3
run "eform sin corregir"      eform Si.cif -k vacancy -q 1 --supercell 2x2x2 --correction ninguna -o o/ef4
run "eform sustitucion"       eform Si.cif -k substitution --new-element P -q 0,1 --supercell 2x2x2 --epsilon 11.7 -o o/ef5
run "eform intersticial"      eform Si.cif -k interstitial --new-element H --position 0.5,0.5,0.5 -q -1,0,1 --supercell 2x2x2 --epsilon 11.7 -o o/ef6
run "eform con mu explicito"  eform Si.cif -k vacancy -q 0 --supercell 2x2x2 --mu Si=-107.8 -o o/ef7
run "eform con dv"            eform Si.cif -k vacancy -q -1,1 --supercell 2x2x2 --epsilon 11.7 --dv -0.05 -o o/ef8

echo "=== EJECUCION: PARALELO, REANUDAR, PRESUPUESTO (0.17.0) ==="
run "sugerencia de -j"        strain Si.cif -r -2:2:5 -o o/rn1
run "run.sh paralelizable"    elastic Si.cif --npoints 2 -o o/rn2
run "flags de ejecucion"      eos Si.cif --npoints 5 -o o/rn3 -j 2 --max-time 2h --redo
run "estimar coste"           strain Si.cif -r -3:3:7 -o o/rn4 --estimate
run "estimar en paralelo"     elastic Si.cif --npoints 2 -o o/rn5 --estimate -j 4
run --rc1 "cost sin base"     cost --db /tmp/no_existe_esta_base.db

echo "=== ENERGIA DE SUPERFICIE (0.17.0) ==="
run "gamma prepara"           gamma Si.cif -m "1 1 1" -l 3,4,5 -o o/gm1
run "gamma sin bulto"         gamma Si.cif -m "1 0 0" -l 3,4 --no-bulk -o o/gm2
run "gamma sin reducir"       gamma Si.cif -m "1 1 1" -l 3,4 --no-reduce -o o/gm3
run "gamma relajando"         gamma Al.cif -m "1 1 1" -l 3,4 --relax --fix 1 -o o/gm4
run "gamma con dipolo"        gamma Al.cif -m "1 1 1" -l 3,4 --dipole -o o/gm5

echo "=== CENTRO DE BANDA Y FATBANDS (0.17.0) ==="
run "centro de banda d"       dos /tmp/test_fe --dband Fe
run "centro de banda p"       dos /tmp/test_si --dband Si-p
run "centro con corte"        dos /tmp/test_fe --dband Fe --dband-emax 5
run "fatbands s"              bands fat_si --fat Si-s -o o/fat1 --format png
run "fatbands p"              bands fat_si --fat Si-p -o o/fat2 --no-plot
run "fatbands por atomo"      bands fat_si --fat atomo:1 -o o/fat3 --no-plot

echo "=== ALINEAMIENTO DE BANDAS (0.17.0) ==="
run "align dos losas"         align hbn_a hbn_b -o o/al1 --no-plot
run "align con nombres"       align hbn_a hbn_b --names "A,B" -o o/al2 --format png

echo "=== FONONES A TEMPERATURA ELECTRONICA (0.17.0) ==="
run "tscan prepara"           phonons Si.cif --tscan 300,2000,6000 --gamma -o o/tp1
run "tscan collect"           phonons Si.cif --tscan 300,2000,6000 --gamma -o tp_si --collect --no-plot
run "tscan figura"            phonons Si.cif --tscan 300,2000,6000 --gamma -o tp_si --collect --format png

echo "=== HUBBARD V INTERSITIO (0.17.0) ==="
run "hubbard con intersitio"  hubbard NiO.cif --collect -o hub_nio --qgrid 2x2x2 --intersite
run "hubbard umbral de V"     hubbard NiO.cif --collect -o hub_nio --qgrid 2x2x2 --intersite --v-threshold 0.5

echo "=== FUNCIONALES HIBRIDOS (0.18.0) ==="
run "gen con HSE"             gen Si.cif -p scf --functional hse -o o/hse1
run "gen con PBE0"            gen Si.cif -p scf --functional pbe0 -o o/hse2
run "gen malla de EXX"        gen Si.cif -p scf --functional hse --exx-grid 1x1x1 -o o/hse3
run "gen fraccion de EXX"     gen Si.cif -p scf --functional hse --exx-fraction 0.4 -o o/hse4

echo "=== LORENZ Y ESPIN EN TRANSPORTE (0.18.0) ==="
run "transport con Lorenz"    transport Si.cif --collect -o /tmp/tr_si --no-plot
run --rc2 "espin sin nspin 2"  transport Si.cif --collect -o /tmp/tr_si --spin-resolved --no-plot

echo "=== VALIDACION FISICA (0.18.0) ==="
run "selftest rapido"         selftest
run "selftest lista"          selftest --list
run "selftest una sola"       selftest --only madelung,lorenz

echo "=== ELECTROQUIMICA (0.19.0) ==="
run "echem HER"               echem --her -0.33 -o o/ec1 --no-plot
run "echem HER figura"        echem --her -0.33 -o o/ec2 --format png
run "echem OER"               echem --oer OH=0.77,O=2.16,OOH=3.87 --corrections OH=0,O=0,OOH=0 -o o/ec3 --no-plot
run "echem con U y pH"        echem --her -0.33 -U 0.2 --ph 7 -o o/ec4 --no-plot

echo "=== DOCUMENTACION (0.19.0) ==="
run "docs genera la pagina"   docs -o o/olla-dft-docs.html
run --rc2 "echem sin reaccion" echem -o o/x --no-plot
run --rc2 "echem las dos"      echem --her -0.3 --oer OH=1 -o o/x --no-plot
run --rc2 "echem OER a medias" echem --oer OH=0.77 -o o/x --no-plot
run --rc2 "echem valor no num" echem --oer OH=mucho,O=1,OOH=2 -o o/x --no-plot
run --rc2 "selftest inventada" selftest --only no_existe
run --rc2 "selftest full sin pseudos" selftest --full

echo "=== PORTABILIDAD (0.35.0) ==="
run "sistema"                 sistema
run "salida ascii"            --ascii recetas
run "ascii detras"            recetas primero --ascii
run "ascii en un informe"     --ascii info Si.cif


echo "=== RECETAS (0.27.0) ==="
run "recetas lista"           recetas
run "recetas una"             recetas mecanicas
run "recetas la de empezar"   recetas primero
run "recetas sin verificar"   recetas termoelectrico
run "recetas busca"           recetas --buscar "quiero saber si conduce el calor"
run "recetas busca sin exito" recetas --buscar "zzzz qqqq"
run "recetas guion"           recetas bandas --script o/sesion.sh
run "recetas guion sin ruta"  recetas modelo --script o/modelo.sh
run --rc2 "recetas inventada"    recetas no_existe
run "derived escribe dat"     derived Si.cif --cij ELASTIC_C.dat -o o/dv


echo "=== EXFOLIACION (0.6.0, sin cubrir hasta la 0.26) ==="
run "exfoliate grafito"       exfoliate gr.cif -o o/ex1 --pseudo-dir /usr/share/espresso/pseudo --ecutwfc 20
run "exfoliate hBN"           exfoliate hbn.cif -o o/ex2 --pseudo-dir /usr/share/espresso/pseudo --ecutwfc 20
run "exfoliate con vdW"       exfoliate gr.cif -o o/ex3 --vdw grimme-d3 --pseudo-dir /usr/share/espresso/pseudo --ecutwfc 20
run --rc2 "exfoliate no laminar" exfoliate Si.cif -o o/x --pseudo-dir /usr/share/espresso/pseudo --ecutwfc 20
run --rc2 "exfoliate vdW inventado" exfoliate gr.cif -o o/x --vdw ninguno --pseudo-dir /usr/share/espresso/pseudo


echo "=== FUNCIONES DE WANNIER (0.21.0) ==="
run "wannier prepara"         wannier Si.cif -o o/w1 -g 2x2x2 -p Si:sp3 --pseudo-dir /usr/share/espresso/pseudo --ecutwfc 20
run "wannier con excluidas"   wannier Si.cif -o o/w2 -g 2x2x2 -p "f=0.125,0.125,0.125:s" --bands 8 --exclude 5-8 --pseudo-dir /usr/share/espresso/pseudo --ecutwfc 20
run "wannier proy automatica" wannier Si.cif -o o/w3 -g 2x2x2 --pseudo-dir /usr/share/espresso/pseudo --ecutwfc 20
run "wannier collect"         wannier Si.cif --collect -o wann --no-plot
run "wannier collect sin min" wannier Si.cif --collect -o wann --no-minimize --no-plot
run "wannier collect con dos" wannier Si.cif --collect -o wann --dos 8 --no-plot
run "wannier desenreda"       wannier Si.cif --collect -o wann12 --window -10:20 --frozen -10:6.4 --no-plot
run "wannier sin congelar"    wannier Si.cif --collect -o wann12 --window -10:20 --no-plot
run --rc2 "wannier ventana al reves"  wannier Si.cif --collect -o wann12 --window 20:-10 --no-plot
run --rc2 "wannier ventana mal"       wannier Si.cif --collect -o wann12 --window abc --no-plot
run --rc2 "wannier ventana estrecha"  wannier Si.cif --collect -o wann12 --window -10:-9 --no-plot
run --rc2 "wannier congela de mas"    wannier Si.cif --collect -o wann12 --window -10:20 --frozen -10:19 --no-plot
run --rc2 "wannier sin estructura"   wannier --collect -o wann
run --rc2 "wannier orbital inventado" wannier Si.cif -o o/x -g 2x2x2 -p Si:sp9 --pseudo-dir /usr/share/espresso/pseudo
run --rc2 "wannier mas wf que bandas" wannier Si.cif -o o/x -g 2x2x2 -p Si:sp3 --bands 4 --pseudo-dir /usr/share/espresso/pseudo
run --rc2 "wannier rango al reves"    wannier Si.cif -o o/x -g 2x2x2 -p Si:s --exclude 8-5 --pseudo-dir /usr/share/espresso/pseudo
run --rc2 "wannier collect sin datos" wannier Si.cif --collect -o o/w1 --no-plot


echo "=== FASE DE BERRY (0.23.0) ==="
run "berry prepara"           berry Si.cif -o o/b1 --nppstr 7 --kperp 2x2 --pseudo-dir /usr/share/espresso/pseudo --ecutwfc 20
run "berry camino"            berry Si.cif -o o/b2 -r Si.cif --nlambda 3 --nppstr 7 --kperp 2x2 --pseudo-dir /usr/share/espresso/pseudo --ecutwfc 20
run "berry desplazamiento"    berry Si.cif -o o/b3 --displace 2:0,0,0.1 --nlambda 3 --nppstr 7 --kperp 2x2 --pseudo-dir /usr/share/espresso/pseudo --ecutwfc 20
run "berry gdir 1"            berry Si.cif -o o/b4 --gdir 1 --nppstr 7 --kperp 2x2 --pseudo-dir /usr/share/espresso/pseudo --ecutwfc 20
run "berry collect"           berry Si.cif --collect -o berry --displace 2:0,0,0.16 --nlambda 5 --nppstr 9 --kperp 6x6 --pseudo-dir /usr/share/espresso/pseudo --ecutwfc 25 --no-plot
run "berry collect figura"    berry Si.cif --collect -o berry --displace 2:0,0,0.16 --nlambda 5 --nppstr 9 --kperp 6x6 --pseudo-dir /usr/share/espresso/pseudo --ecutwfc 25 --format png
run --rc2 "berry gdir invalido"   berry Si.cif -o o/x --gdir 5 --pseudo-dir /usr/share/espresso/pseudo
run --rc2 "berry cuerda corta"    berry Si.cif -o o/x --nppstr 2 --kperp 2x2 --pseudo-dir /usr/share/espresso/pseudo
run --rc2 "berry kperp mal"       berry Si.cif -o o/x --kperp 6 --pseudo-dir /usr/share/espresso/pseudo
run --rc2 "berry displace mal"    berry Si.cif -o o/x --displace 2 --pseudo-dir /usr/share/espresso/pseudo
run --rc2 "berry displace corto"  berry Si.cif -o o/x --displace 2:0,0 --pseudo-dir /usr/share/espresso/pseudo
run --rc2 "berry atomo inexistente" berry Si.cif -o o/x --displace 9:0,0,0.1 --nppstr 5 --kperp 2x2 --pseudo-dir /usr/share/espresso/pseudo
run --rc2 "berry sin datos"       berry Si.cif --collect -o o/b1 --no-plot


echo "=== CONDUCTIVIDAD TERMICA DE RED (0.24.0) ==="
run "kappa prepara"           kappa Si.cif -o o/k1 --dim 2x2x2 --pseudo-dir /usr/share/espresso/pseudo --ecutwfc 20
run "kappa fc2 aparte"        kappa Si.cif -o o/k2 --dim 2x2x2 --dim-fc2 3x3x3 --pseudo-dir /usr/share/espresso/pseudo --ecutwfc 20
run "kappa con mlip"          kappa Si.cif -o o/k3 --dim 2x2x2 --mesh 7 --temps 300,500 --model mace --no-plot
run "kappa isotopos"          kappa Si.cif -o o/k4 --dim 2x2x2 --mesh 7 --temps 300 --model mace --isotopes --no-plot
run "kappa grano"             kappa Si.cif -o o/k5 --dim 2x2x2 --mesh 7 --temps 300 --model mace --grain 0.1 --no-plot
run --rc2 "kappa dim mal"        kappa Si.cif -o o/x --dim 2x2 --pseudo-dir /usr/share/espresso/pseudo
run --rc2 "kappa temps mal"      kappa Si.cif -o o/x --temps 100:800 --model mace
run --rc2 "kappa temps sin numero" kappa Si.cif -o o/x --temps abc --model mace
run --rc2 "kappa demasiado grande" kappa Si.cif -o o/x --dim 4x4x4 --pseudo-dir /usr/share/espresso/pseudo
run --rc2 "kappa sin fuerzas"    kappa Si.cif --collect -o o/k1 --no-plot
run --rc2 "kappa modelo raro"    kappa Si.cif -o o/x --dim 2x2x2 --model inventado --no-plot


echo "=== SUPERFICIES CARGADAS CON ESM (0.25.0) ==="
run "esm bc1 prepara"         esm Al111.cif -o o/s1 --bc bc1 --pseudo-dir /usr/share/espresso/pseudo --ecutwfc 20
run "esm bc3 con cargas"      esm Al111.cif -o o/s2 --bc bc3 --charge -0.05,0,0.05 --pseudo-dir /usr/share/espresso/pseudo --ecutwfc 20
run "esm bc2 con campo"       esm Al111.cif -o o/s3 --bc bc2 --field 0.001 --pseudo-dir /usr/share/espresso/pseudo --ecutwfc 20
run "esm collect bc1"         esm Al111.cif --collect -o esm1 --bc bc1 --pseudo-dir /usr/share/espresso/pseudo --ecutwfc 25 --ecutrho 100 --kspacing 0.25 --no-plot
run "esm collect figura"      esm Al111.cif --collect -o esm1 --bc bc1 --pseudo-dir /usr/share/espresso/pseudo --ecutwfc 25 --ecutrho 100 --kspacing 0.25 --format png
run "esm collect bc3"         esm Al111.cif --collect -o esm3 --bc bc3 --charge -0.04,-0.02,0,0.02,0.04 --pseudo-dir /usr/share/espresso/pseudo --ecutwfc 25 --ecutrho 100 --kspacing 0.25 --no-plot
run --rc2 "esm bc1 con carga"    esm Al111.cif -o o/x --bc bc1 --charge 0.1 --pseudo-dir /usr/share/espresso/pseudo
run --rc2 "esm bc inventada"     esm Al111.cif -o o/x --bc bc9 --pseudo-dir /usr/share/espresso/pseudo
run --rc2 "esm carga no numerica" esm Al111.cif -o o/x --bc bc3 --charge mucho --pseudo-dir /usr/share/espresso/pseudo
run --rc2 "esm sin carga"        esm Al111.cif -o o/x --bc bc3 --charge "" --pseudo-dir /usr/share/espresso/pseudo
run --rc2 "esm sin datos"        esm Al111.cif --collect -o o/s1 --bc bc1 --pseudo-dir /usr/share/espresso/pseudo --no-plot
run --rc2 "esm celda torcida"    esm torcida.vasp -o o/x --bc bc1 --pseudo-dir /usr/share/espresso/pseudo


echo "=== SOLIDOS AMORFOS (0.21.0) ==="
run "amorfo solo empaqueta"   amorphous SiO2 -n 4 -d 2.2 --pack-only -o o/am1
run "amorfo otra semilla"     amorphous SiO2 -n 4 -d 2.2 --pack-only --seed 7 -o o/am2
run "amorfo binario"          amorphous GeTe -n 6 -d 6.0 --pack-only -o o/am3
run "amorfo min-dist"         amorphous SiO2 -n 4 -d 2.2 --pack-only --min-dist 0.6 -o o/am4
run "amorfo fundido y temple" amorphous SiO2 -n 2 -d 2.2 --melt 2500 --melt-steps 5 --quench-steps 10 --anneal-steps 5 -o o/am5
run --rc2 "amorfo sin densidad"  amorphous SiO2 -n 4 --pack-only -o o/x
run --rc2 "amorfo densidad cero" amorphous SiO2 -n 4 -d 0 --pack-only -o o/x
run --rc2 "amorfo elemento raro" amorphous XxO2 -n 4 -d 2.2 --pack-only -o o/x
run --rc2 "amorfo no cabe"       amorphous SiO2 -n 4 -d 40 --pack-only -o o/x

echo "=== ML ==="
run "mlip relax"          mlip relax Si.cif -o o/si_mlip.cif --steps 30
run "mlip scan"           mlip scan Si.cif --npoints 5 --span 0.04

echo
echo "==================================="
printf 'total OK=%d  FALLAS=%d\n' "$ok" "$bad"
echo "log completo: $LOG"
exit $bad
