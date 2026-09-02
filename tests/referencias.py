"""Valores de referencia congelados.

Cada número aquí se validó UNA VEZ contra un cálculo real de Quantum
ESPRESSO o contra una fuente independiente (experimento, ficha PDF, otra
implementación). A partir de ese momento su función es distinta: no es
documentación, es un detector de regresiones. Si un cambio en el código
mueve uno de estos números, la prueba falla y hay que entender por qué
antes de actualizarlo.

Actualizar un valor de aquí para "que pase la prueba" anula el propósito
del archivo. Solo se cambia cuando se sabe por qué cambió y el nuevo valor
también se validó contra la fuente externa.

Cada entrada dice contra QUÉ se validó y con qué tolerancia.
"""

# --- Silicio, LDA (Si.pz-vbc.UPF), ecutwfc 60 Ry ----------------------
SI_GAP_INDIRECTO = 0.4987      # eV, camino de bandas; ±0.005
SI_GAP_DIRECTO = 2.5694        # eV, en Gamma;         ±0.005
SI_A0 = 5.402                  # A, de la EOS;         ±0.005   (exp 5.431)
SI_B0 = 94.2                   # GPa, EOS;             ±0.5     (exp ~98)
SI_C11, SI_C12, SI_C44 = 159.9, 61.7, 76.6   # GPa; ±1.0  (exp 165.8/63.9/79.6)

# Fonones, malla q 2x2x2. Referencia externa: dispersion inelastica de
# neutrones (Nilsson & Nelin 1972). LDA queda 1-6 % por debajo.
SI_FONON_GAMMA_TO = 508.9      # cm-1; ±1.0   (exp 517)
SI_FONON_X_TA = 140.8          # cm-1; ±1.0   (exp 150)
SI_FONON_X_LA = 406.5          # cm-1; ±1.0   (exp 410)
SI_FONON_L_TA = 107.7          # cm-1; ±1.0   (exp 114)
SI_FONON_L_TO = 484.8          # cm-1; ±1.0   (exp 490)
SI_ZPE_MEV = 122.25            # meV/celda; ±0.5
SI_CV_300K = 0.411             # meV/K por celda de 2 atomos; ±0.005
#   contraste externo: 20 J/(mol*K) -> 0.415 meV/K por celda

# Opticas (epsilon.x, nscf 14x14x14, 27 bandas, ensanchamiento 0.1 eV)
SI_EPS1_0_SIN_SCISSOR = 16.13  # ±0.05
SI_EPS1_0_SCISSOR_065 = 10.44  # ±0.05
SI_PICO_EPS2_SCISSOR = 4.30    # eV; ±0.05. Coincide con el punto critico E2
SI_PLASMON = 16.9538           # eV, reportado por epsilon.x; exp 16.7
# La extrapolacion de Tauc sobre el espectro real cae en el gap DIRECTO,
# porque epsilon.x no incluye transiciones asistidas por fonones.
SI_TAUC_DIRECTO = 2.562        # eV; ±0.05. Gap directo del calculo: 2.5694
SI_TAUC_SCISSOR_065 = 3.214    # eV; ±0.05

# Masa efectiva (camino fino, +-0.06 A^-1, 21 puntos)
SI_ME_LONGITUDINAL = 0.949     # ±0.01   (exp 0.916)
SI_ME_TRANSVERSAL = 0.193      # ±0.005  (exp 0.190)
SI_MH_100_PESADO = -0.269      # ±0.01   (Luttinger 0.277)
SI_MH_111_PESADO = -0.670      # ±0.01   (Luttinger 0.718)

# --- Grafito / grafeno, LDA ------------------------------------------
GRAFITO_ESPACIADO = 3.356      # A; ±0.005
HBN_ESPACIADO = 3.331          # A; ±0.005
GRAFITO_EXFOLIACION = 25.8     # meV/atomo; ±0.5 (valor de literatura LDA)
GRAFENO_FUNCION_TRABAJO = 4.539  # eV; ±0.01 (exp ~4.6)

# --- Difraccion: contra fichas PDF ------------------------------------
# Si, PDF 27-1402. NaCl, PDF 05-0628. Tolerancia 0.06 grados en 2theta.
SI_XRD = [(28.465, "(111)"), (47.343, "(220)"),
          (56.171, "(311)"), (69.193, "(400)")]
NACL_XRD = [(27.390, "(111)"), (31.731, "(200)"), (45.488, "(220)"),
            (53.917, "(311)"), (56.524, "(222)"), (66.289, "(400)")]
GRAFITO_XRD_002 = 26.564       # grados; PDF 41-1487 da 26.54
