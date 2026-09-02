# Ficha del material: Si2

*Generado por Olla-DFT 0.11.0 — 2026-08-28 19:36 UTC*

## Estructura

| magnitud | valor | unidad |
|---|---|---|
| volumen de celda | 39.4019 | Å³ |
| átomos por celda | 2 |  |
| operaciones de simetría | 48 |  |

*Fuente: `/tmp/proy/out/Si.xml`*

## Ecuación de estado

| magnitud | valor | unidad |
|---|---|---|
| parámetro de red a₀ | 5.402 | Å |
| módulo volumétrico B₀ | 94.2 | GPa |
| volumen de equilibrio | 40.05 | Å³ |

*Fuente: `/tmp/proy/EOS.txt`*

## Elásticas

| magnitud | valor | unidad |
|---|---|---|
| C₁₁ | 159.9 | GPa |
| C₁₂ | 61.7 | GPa |
| C₄₄ | 76.6 | GPa |
| módulo volumétrico (VRH) | 94.43 | GPa |
| módulo de corte (VRH) | 64.09 | GPa |
| razón de Poisson | 0.2233 |  |
| estable (Born) | sí |  |

*Fuente: `/tmp/proy/ELASTIC_C.dat`*

## Fonones

| magnitud | valor | unidad |
|---|---|---|
| energía de punto cero | 122.25 | meV por celda |
| C_v (300 K) | 0.4113 | meV/K por celda |

*Fuente: `/tmp/proy/FONONES_TERMO.dat`*

## Ópticas

| magnitud | valor | unidad |
|---|---|---|
| ε₁(0) | 16.132 |  |
| n(0) | 4.0165 |  |

*Fuente: `/tmp/proy/OPTICS.dat`*

## Masa efectiva

| magnitud | valor | unidad |
|---|---|---|
| m* hueco ([100]) | -0.2693 | mₑ |
| m* hueco ([100]) | -0.2693 | mₑ |
| m* hueco ([100]) | -0.178 | mₑ |
| m* hueco ([110]) | -2.7152 | mₑ |
| m* hueco ([110]) | -0.266 | mₑ |
| m* hueco ([110]) | -0.1136 | mₑ |

*Fuente: `/tmp/proy/MASA_EFECTIVA.dat`*

## Parámetros del cálculo

| parámetro | valor |
|---|---|
| funcional | PZ |
| ecutwfc_Ry | 60.0 |
| ecutrho_Ry | 480.0 |
| malla_k | 11x11x11 |
| pseudos | Si: Si.pz-vbc.UPF |
| ocupaciones | fixed |
| smearing | — |
| degauss_Ry | — |
| nspin | 1 |

## Métodos (borrador)

Los cálculos de primeros principios se realizaron con Quantum ESPRESSO [1], en el marco de la teoría del funcional de la densidad. Se empleó el funcional de intercambio y correlación PZ y los pseudopotenciales Si (Si.pz-vbc.UPF). Las funciones de onda y la densidad de carga se expandieron en ondas planas con energías de corte de 60.0 y 480.0 Ry respectivamente. La zona de Brillouin se muestreó con una malla de Monkhorst-Pack de 11x11x11, con ocupaciones fijas. Las propiedades vibracionales se obtuvieron por teoría del funcional de la densidad perturbativa [2]. El análisis de simetría y la generación de los caminos de alta simetría se hicieron con spglib y seekpath. El pre y post-proceso se realizó con Olla-DFT 0.11.0.

## Referencias

1. P. Giannozzi et al., J. Phys.: Condens. Matter 21, 395502 (2009); J. Phys.: Condens. Matter 29, 465901 (2017)
2. S. Baroni, S. de Gironcoli, A. Dal Corso, P. Giannozzi, Rev. Mod. Phys. 73, 515 (2001)
3. A. Dal Corso, S. Baroni, R. Resta, Phys. Rev. B 49, 5323 (1994) — respuesta dieléctrica
4. A. Togo, I. Tanaka, arXiv:1808.01590 (2018)
5. Y. Hinuma, G. Pizzi, Y. Kumagai, F. Oba, I. Tanaka, Comput. Mater. Sci. 128, 140 (2017)
6. A. Hjorth Larsen et al., J. Phys.: Condens. Matter 29, 273002 (2017)

---

*El párrafo de métodos es un BORRADOR generado de los parámetros reales del
cálculo. Revísalo antes de usarlo: sabe qué se hizo, no por qué se hizo.*