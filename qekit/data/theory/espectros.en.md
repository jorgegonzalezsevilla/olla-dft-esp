## Spectra, surfaces, chemistry and quality control

This part documents the physics behind the Olla-DFT commands that go beyond the total energy: optical and X-ray spectra (`optics`, `tddft`, `corehole`, `xanes`, `xps`), analysis of the electron density and of surfaces (`charges`, `charge`, `wf`, `esm`, `surface`, `adsorb`, `interface`), defect and reaction chemistry (`defect`, `eform`, `echem`, `neb`, `hull`), structure generation with machine-learned potentials (`amorphous`, `mlip`) and the quality-control tools that check that all of the above is comparable and credible (`audit`, `db`, `doctor`, `crosscheck`, `selftest`, `suggest`, `pseudos`). Each section states which question the command answers, which formulas the code actually implements (with the responsible function), which Quantum ESPRESSO file each datum comes from, and where the limits are. Whenever the code's internal documentation promises something the code does not do, it is said under "Limits and pitfalls".

---

### `olla-dft optics` — Dielectric function, absorption and Tauc gap

**What it answers.** How does the material respond to light? It gives $\varepsilon(\omega)$, the refractive index $n$, the extinction coefficient $k$, the absorption coefficient $\alpha$, the reflectivity $R$ and an optical gap extrapolated the way it is done with a UV-Vis spectrum.

**Background for non-experts.** When light crosses a solid, the electric field of the wave pushes the electrons. If the photon energy matches what an electron needs to jump from an occupied band to an empty one, the light is absorbed. The *dielectric function* $\varepsilon(\omega) = \varepsilon_1 + i\varepsilon_2$ summarises that response: the imaginary part $\varepsilon_2$ counts how many jumps exist at each energy (absorption) and the real part $\varepsilon_1$ how much the material polarises (refraction). The two are not independent: causality (the response cannot precede the cause) ties them through the Kramers–Kronig relations, so knowing one lets you rebuild the other.

`epsilon.x` computes $\varepsilon_2$ by summing all *vertical* transitions between Kohn–Sham bands, as if each electron jumped alone, without feeling the hole it leaves behind. This is the independent-particle approximation (RPA without local fields). The gap that comes out is the functional's, usually too small; the "scissor" is the simplest correction: all transitions are shifted rigidly by $\Delta$ and $\varepsilon_1$ is rebuilt by Kramers–Kronig so as not to break causality.

**Formulas.** All in `qekit/modules/optics.py`.

Isotropic average (`optics.collect`):
$$\varepsilon_{1,2}(\omega) = \tfrac{1}{3}\left[\varepsilon_{xx} + \varepsilon_{yy} + \varepsilon_{zz}\right]$$

Derived optical functions (`optics.derived`):
$$|\varepsilon| = \sqrt{\varepsilon_1^2 + \varepsilon_2^2},\qquad n = \sqrt{\frac{|\varepsilon| + \varepsilon_1}{2}},\qquad k = \sqrt{\frac{|\varepsilon| - \varepsilon_1}{2}}$$
$$\alpha(E) = \frac{2\,k\,E}{\hbar c},\qquad R = \frac{(n-1)^2 + k^2}{(n+1)^2 + k^2}$$

- $E = \hbar\omega$: photon energy (eV). $\hbar c$ = `HBAR_C_EV_CM` = $1.9732698\times10^{-5}$ eV·cm, so that $\alpha$ comes out in cm⁻¹. $n$, $k$, $R$ are dimensionless. Negative radicands are clipped to zero (`np.maximum`).

Kramers–Kronig (`optics.kramers_kronig`):
$$\varepsilon_1(\omega) = 1 + \frac{2}{\pi}\,\mathcal{P}\!\int_0^{\omega_{\max}} \frac{\omega'\,\varepsilon_2(\omega')}{\omega'^2 - \omega^2}\,d\omega'$$
- $\mathcal{P}$: principal value; implemented by removing the point $\omega'=\omega$ from the trapezoidal quadrature on the uniform `epsilon.x` grid. The integral is truncated at `wmax`.

Scissor (`optics.scissor`):
$$\varepsilon_2'(E) = \varepsilon_2(E-\Delta)\left(\frac{E-\Delta}{E}\right)^2,\qquad \varepsilon_1' = \mathrm{KK}[\varepsilon_2']$$
- $\Delta$: shift in eV (`--scissor`). The factor $((E-\Delta)/E)^2$ comes from $\varepsilon_2 \propto |p|^2/\omega^2$ with untouched matrix elements $|p|^2$. It is applied to each Cartesian component and then averaged.

Tauc plot (`optics.tauc_gap`):
$$y(E) = \left(\alpha E\right)^{1/r},\qquad r = \tfrac{1}{2}\ (\text{allowed direct}),\quad r = 2\ (\text{indirect})$$
$$E_g^{\mathrm{opt}} = -\frac{b}{m}\quad\text{with}\quad y \approx m E + b\ \text{fitted on the first absorption edge}$$

**How Olla-DFT computes it.**
1. `optics.prepare` resolves pseudopotentials and cutoffs (`sweep.prepare_common`, task `optics`) and **refuses** if any is not norm-conserving (`epsilon.x` has no matrix elements for USPP/PAW).
2. Writes `scf.in` (grid from the configured `kspacing`, 0.20 Å⁻¹ by default), `nscf.in` with a dense grid (`--kspacing`, default 0.12 Å⁻¹), `nosym=.true.` and `nbnd = 3 ×` the estimate of `inputgen._estimate_nbnd` (`nbnd_factor=3.0`), and `epsilon.in` (`calculation='eps'`, `smeartype='gauss'`, `intersmear=--smear` (0.10 eV), `wmin=0`, `wmax=--wmax` (20 eV), `nw=800`).
3. With `--run`: `pw.x` scf → `pw.x` nscf (`runner.run_all`) → `optics.run_epsilon` launches `epsilon.x` (looked up next to `pw.x`).
4. `optics.collect` reads `epsr_<prefix>.dat` and `epsi_<prefix>.dat` (columns: energy, xx, yy, zz) and averages.
5. If `--scissor Δ ≠ 0`: `optics.scissor` shifts $\varepsilon_2$ and rebuilds $\varepsilon_1$ with `optics.kramers_kronig`.
6. `optics.derived` yields $n, k, \alpha, R$; `optics.tauc_gap` fits the direct and indirect gaps; `optics.report` prints $\varepsilon_1(0)$ (value at $E \approx 0.05$ eV), $n(0)$, the maximum of $\varepsilon_2$ and the gaps.
7. `optics.export` writes `OPTICS.dat` with the columns of `optics.OPTICS_COLUMNS` (`E(eV)`, `eps1`, `eps2`, `n`, `k`, `alpha(1/cm)`, `R`), named in the last comment line so that `optics.read_optics_dat` can read them by name; `optics.plot` draws the three-panel figure.

Detail of the Tauc fit (`optics.tauc_gap`): the curve is smoothed with a ~0.05 eV moving average; the noise floor is the maximum of $y$ in the first 1 % of the spectrum; the edge starts where $y$ exceeds $\max(2\cdot\text{floor}, 10^{-3}\cdot\mathrm{median}(y>0))$ and $E > 0.1$ eV; it ends at the first local maximum that triples the onset value or at most `max_span` = 1.5 eV higher; a straight line is fitted in a window of `fit_window` = 0.6 eV centred on the steepest point of that stretch. Returns `None` if there is no absorption, if the slope is not positive or if the intercept falls outside the range.

**Where each datum comes from.**

| Datum | Source | Detail |
|---|---|---|
| $\varepsilon_1(\omega)$ per direction | `epsr_<prefix>.dat` from `epsilon.x` | `optics.collect`, columns 1–3 |
| $\varepsilon_2(\omega)$ per direction | `epsi_<prefix>.dat` from `epsilon.x` | `optics.collect` |
| $\hbar c$ | constant `optics.HBAR_C_EV_CM` | $1.9732698\times10^{-5}$ eV·cm |
| Type of each pseudo (NC/US/PAW) | UPF header (`pseudo_type`) | via `sweep.prepare_common` |
| $\Delta$ (scissor) | parameter `--scissor` | eV; recommended exp./GW gap − DFT gap |
| Broadening | parameter `--smear` | Gaussian `intersmear`, 0.10 eV |
| Window and points | `--wmax` (20 eV), `nw=800` | fixed in `optics.prepare` |
| nscf grid | `--kspacing` (0.12 Å⁻¹) | `sweep.default_grid` |

**Limits and pitfalls.**
- It is independent-particle RPA: no local fields, no excitons. The report reminds you: *"Recuerda: RPA de partícula independiente y gap del funcional…"*.
- Without NC pseudos the command aborts: *"epsilon.x solo funciona con pseudopotenciales de NORMA CONSERVADA…"*.
- `epsilon.x` does not include phonon-assisted transitions: in an indirect semiconductor $\varepsilon_2 = 0$ below the direct gap and the "indirect" fit does **not** give the true indirect gap (`tauc_gap` docstring).
- The Kramers–Kronig integral is truncated at `wmax`: $\varepsilon_1(0)$ inherits an error if there is strong absorption above 20 eV.
- The scissor only moves the gap; it neither corrects intensities nor adds excitons.
- If the Tauc fit fails, the report prints *"no se pudo ajustar"* instead of a number.

**References.**
- J. Tauc, R. Grigorovici, A. Vancu, *Phys. Status Solidi* 15, 627 (1966) — Tauc plot.
- `epsilon.x` manual (Quantum ESPRESSO, PP package): A. Benassi, *"epsilon.x: a post-processing tool for the calculation of the dielectric properties"*.
- M. Dressel, G. Grüner, *Electrodynamics of Solids* (Cambridge, 2002) — Kramers–Kronig and optical functions.

---

### `olla-dft tddft` — Optical absorption with TDDFPT (Lanczos/Davidson)

**What it answers.** Does the absorption spectrum change when the excited electron and the hole it leaves see each other? Where are the first excitations, which are bright and which dark, and is there a bound exciton below the gap?

**Background for non-experts.** `optics` sums one-electron transitions one at a time. In reality the excited electron (charge −) and the hole (charge +) attract each other; in molecules and wide-gap insulators that attraction lowers the energy of the pair and creates an absorption peak **inside** the gap: the exciton. Time-dependent density functional theory in linear response (TDDFPT) includes that interaction through the exchange-correlation kernel. Quantum ESPRESSO solves it in two ways: with the **Lanczos** algorithm (`turbo_lanczos.x` + `turbo_spectrum.x`), which gives the whole spectrum without computing empty states, or with **Davidson** (`turbo_davidson.x`), which gives the first N excitations one by one with their energy and oscillator strength $f$. An excitation with $f \approx 0$ exists but does not absorb light: it is "dark".

**Formulas.** In `qekit/modules/tddft.py`.

Unit conversion of the inputs (`tddft.build_lanczos_input`, `build_spectrum_input`, `build_davidson_input`):
$$E_{\mathrm{Ry}} = \frac{E_{\mathrm{eV}}}{\mathrm{RY\_EV}},\qquad \mathrm{RY\_EV} = 13.605693122994\ \mathrm{eV}$$

Wavelength (`tddft.report`):
$$\lambda\,(\mathrm{nm}) = \frac{1239.84}{E\,(\mathrm{eV})}$$

Absorption onset (`TddftRun.onset`): first local maximum of $dS/dE$ exceeding 20 % of the maximum of the derivative (inflection point of the first rise).

Exciton signature (`tddft._avisar`):
$$d = E_{\mathrm{onset}} - E_g^{\mathrm{IP}},\qquad d < -\max(0.10\ \mathrm{eV},\ 2\,\eta)\ \Rightarrow\ \text{bound exciton}$$
- $E_g^{\mathrm{IP}}$: independent-particle gap supplied by the user (`--gap`); $\eta$: broadening in eV, from `--broadening` at `--collect` or, if omitted, read from `spectrum.in` (`epsil`) or `davidson.in` (`broadening`) by `tddft._broadening_de_inputs` (Ry → eV). `UMBRAL_EXCITON` = 0.10 eV; `BROADENING_DEFAULT` = 0.05 eV.

Anisotropy (`tddft._anisotropia`): $\max_E[\max_i S_i(E) - \min_i S_i(E)] / \max_{i,E} S_i(E)$ over the $x,y,z$ components.

**How Olla-DFT computes it.**
1. `tddft.prepare`: if the minimum vacuum (`_vacio_minimo`) exceeds 5 Å or `--gamma` is given, it uses `K_POINTS gamma` (the only case TDDFPT implements); otherwise an automatic grid with a warning that `turbo_*.x` will stop. Writes `scf.in` with `nosym` and `noinv`.
2. Lanczos: `lanczos.in` (`itermax=--iter` 500, `ipol=--pol` 4 → `n_ipol=3`, `ltammd` with `--tamm-dancoff`, `lrpa` with `--rpa`, `scissor=--scissor/RY_EV` if a rigid shift of the empty bands is requested; `prepare` rejects a negative scissor or one combined with `--method davidson`) and `spectrum.in` (`itermax0=itermax`, `itermax=4×itermax` for the `--extrapolation` osc/constant/no, `epsil=--broadening/RY_EV`, `units=1` (eV), `start/end/increment`).
3. Davidson: `davidson.in` with `num_eign=--states` (10), `num_init=2N`, `num_basis_max=max(80, 8N)`, `residue_conv_thr=1e-4`, `p_nbnd_virt=15`, window and `broadening` in Ry, `reference` at the window centre.
4. The user runs `pw.x` → `turbo_lanczos.x` → `turbo_spectrum.x` (or `turbo_davidson.x`) by hand.
5. `tddft.collect --collect` (with `broadening` from `--broadening` or from the inputs): Lanczos reads the first `*plot*.dat` (columns: E in eV, total S, S_x, S_y, S_z) and, from `lanczos.out`, the `itermax` and the functional. Davidson (`_collect_davidson`) reads `<prefix>.eigen` (energy in Ry → eV, total strength, strengths per direction) and the `*plot*.dat` if present.
6. `_picos` lists local maxima above 5 % of the maximum; `_avisar` compares the onset with `--gap`; `report` flags as bright the excitations with $f > 0.01$ and counts the dark ones.
7. `export` writes `TDDFT.dat`, `TDDFT_EXCITACIONES.dat`, `TDDFT.txt`; `plot` optionally overlays the `optics` spectrum (`--compare OPTICS.dat`: the CLI reads the `alpha(1/cm)` column by name with `optics.read_optics_dat` and normalises it to the TDDFPT maximum).

**Where each datum comes from.**

| Datum | Source | Detail |
|---|---|---|
| $S(E)$ and components | `<prefix>.plot_S.dat` (or `*plot*.dat`) from `turbo_spectrum.x` | `tddft.collect`; energy in eV because of `units=1` |
| Excitations $(E, f, f_x, f_y, f_z)$ | `<prefix>.eigen` from `turbo_davidson.x` | `_collect_davidson`; E in Ry × `RY_EV` |
| `itermax`, functional | `lanczos.out` | regex in `tddft.collect` |
| IP gap | parameter `--gap` | eV |
| Broadening $\eta$ | `--broadening` or `spectrum.in`/`davidson.in` | `tddft._broadening_de_inputs`, Ry × `RY_EV` |
| Scissor | `--scissor` (Lanczos only) | eV → Ry in `lanczos.in` |
| $\alpha(E)$ from `optics` | `OPTICS.dat`, column `alpha(1/cm)` | `optics.read_optics_dat` |
| Ry → eV | `tddft.RY_EV` | 13.605693122994 |
| Minimum vacuum | cell geometry | `tddft._vacio_minimo` |

**Limits and pitfalls.**
- With LDA/GGA the adiabatic kernel does **not** bind excitons in a solid; the report says so: *"con LDA o GGA el kernel adiabático NO liga excitones en un SÓLIDO… En MOLÉCULAS sí mejora."*
- Γ point only: with a k grid the report warns *"OJO: TDDFPT solo tiene implementado el caso gamma y se plantará al leer el input"*.
- Molecule with < 6 Å of vacuum and < 30 atoms: *"AVISO: solo hay X Å de vacío…"*.
- `--scissor` only exists in `turbo_lanczos.x`: with `--method davidson` or with a negative value the command aborts (*"--scissor solo existe en turbo_lanczos.x…"*). The TDDFPT scissor rebuilds nothing by Kramers–Kronig: QE's own code applies it to the empty bands.
- `--compare` requires an `OPTICS.dat` with the `alpha(1/cm)` column; if missing: *"'…' no tiene la columna 'alpha(1/cm)'; --compare espera el OPTICS.dat de 'olla-dft optics'."*
- If neither `--broadening` nor the inputs give the broadening, the exciton threshold stays at `UMBRAL_EXCITON` = 0.10 eV.
- The command does not launch the `turbo_*.x` executables: it only writes inputs and reads outputs.

**References.**
- D. Rocca, R. Gebauer, Y. Saad, S. Baroni, *J. Chem. Phys.* 128, 154105 (2008) — TDDFPT Lanczos.
- O. B. Malcıoğlu, R. Gebauer, D. Rocca, S. Baroni, *Comput. Phys. Commun.* 182, 1744 (2011) — turboTDDFT.
- X. Ge, S. J. Binnie, D. Rocca, R. Gebauer, S. Baroni, *Comput. Phys. Commun.* 185, 2080 (2014) — turboTDDFT 2.0 (Davidson).

---

### `olla-dft corehole` — Core-hole pseudopotentials (ld1.x)

**What it answers.** How to describe an atom from which an electron of an inner shell has been removed? It generates the pair of pseudopotentials (normal + core-hole) that `xps` and `xanes` need, with the same configuration and the same radii, and extracts the core wavefunction that `xspectra.x` reads.

**Background for non-experts.** A pseudopotential replaces the nucleus and the inner ("core") electrons by an effective potential, so that the calculation only treats the valence electrons. To simulate an X-ray spectroscopy one must remove an electron from that frozen core: this requires a different pseudopotential, generated on purpose with the atomic program `ld1.x`, in which the occupation of the core level (1s for the K edge, 2p for L₂,₃, etc.) is one less. Since the core has one electron less, the declared valence charge `z_valence` rises by exactly 1: that unit **is** the hole. The two pseudos must be generated together with the same parameters, because comparing energies made with pseudos from different families means nothing.

**Formulas.** This module evaluates no physical formulas; it builds `ld1.x` inputs from explicit rules in `qekit/core/atomconf.py` and `qekit/modules/corehole.py`:

- Electronic configuration by Aufbau filling (`atomconf.aufbau`, Madelung order `ORDEN`) with the exceptions in `atomconf.EXCEPCIONES` (Cr, Cu, Nb, Mo, Ru, Rh, Pd, Ag, La, Ce, Gd, Pt, Au).
- Core/valence partition (`atomconf.particion`): valence = $s,p$ shell of $n_{\max}$ + any partially filled $d$/$f$ + filled $d$ of the previous row; with `--semicore`, also $(n-1)s,(n-1)p$.
- Hole (`atomconf.config_hueco`): occupation of the `BORDES[edge]` level reduced by 1.0; rejected if the level is not in the core.
- Pseudisation channels (`atomconf.canales_pseudo`): the valence plus one **unoccupied** channel (occupation −2) for every missing $l \le 2$, with $n = \max(n_{\max}, l+1)$; with `--projectors 2` a second projector per channel labelled $n+1$ with occupation −1.
- Cutoff radius per row (`corehole.RCUT_FILA`): {1: 1.0, 2: 1.3, 3: 1.7, 4: 2.0, 5: 2.2, 6: 2.4} bohr; `rcutus = 1.25 · rcut` only if `pseudotype=3`.
- Reference energies of unbound channels: `E_CANAL_VACIO` = 0.15 Ry; second projector `E_SEGUNDO_PROYECTOR` = 0.05 Ry.

**How Olla-DFT computes it.**
1. `corehole.generar`: validates the element (H..Rn), forces `pseudotype=3` if 2 projectors are requested, obtains partition, channels and `rcut` (or `--rcut`).
2. `corehole.input_ld1` writes `ld1_base.in` and `ld1_hueco.in` (`iswitch=3`, `rel=--rel`, `beta=0.3`, `dft=--functional` (PBE), `tm=.true.`, `lloc` = highest $l$ among the channels, `lgipaw_reconstruction=.true.`, `author='Olla-DFT'`). The empty channels are also added to the all-electron configuration (`_con_canales_vacios`) because `ld1.x` requires them to exist.
3. `corehole._correr_ld1` runs `ld1.x < ld1_X.in > ld1_X.out` (unless `--only-inputs`) and fails if `Error in routine` appears.
4. `corehole.leer_upf` reads from each UPF `element`, `z_valence`, `mesh_size`, `pseudo_type`, `functional`, `wfc_cutoff`, `rho_cutoff` and the `PP_GIPAW_CORE_ORBITAL` labels.
5. `corehole.verificar` applies the checks in the following table; `report` and `export` (`PSEUDOS_HUECO.txt`) list them. The exit code is 1 if there is any `FALLA`.
6. With `--core-wfc UPF`: `corehole.core_wfc` extracts the core wavefunctions in the `filecore` format of `xspectra.x` (one block per orbital, separated by a blank line, in UPF order) and verifies that the number of points matches `mesh_size`.

| Check (`corehole.verificar`) | Criterion | Flag |
|---|---|---|
| Difference in `z_valence` | exactly +1 (tolerance 1e-6) | FALLA otherwise |
| Radial meshes | `mesh_size` equal in both UPFs | FALLA otherwise |
| Hole orbital | present among the `PP_GIPAW_CORE_ORBITAL` of the core-hole UPF | FALLA otherwise |
| Functional | same in both UPFs | FALLA otherwise |
| Projectors | warning if only one per channel (XSpectra recommends two) | warning |
| Ghost states, logarithmic derivatives, transferability | **not checked** | explicit warning |

**Where each datum comes from.**

| Datum | Source | Detail |
|---|---|---|
| Electronic configuration | `atomconf.aufbau` + `EXCEPCIONES` | rule, not experimental data |
| Edge level | `atomconf.BORDES` | K=1s, L1=2s, L23=2p, M1=3s, M23=3p, M45=3d |
| `z_valence`, `mesh_size`, type, functional | `PP_HEADER` of the generated UPF | `corehole.leer_upf` |
| Core orbitals | `PP_GIPAW_CORE_ORBITAL.n` sections of the UPF | `leer_upf`, `core_wfc` |
| Radial mesh | `PP_R` of the UPF | `core_wfc` |
| Cutoff radius | `RCUT_FILA` or `--rcut` | bohr |

**Limits and pitfalls.**
- The M edges (`M1`, `M23`, `M45`) exist in `atomconf.BORDES` and serve to generate the core-hole pseudo (XPS), but `xspectra.x` only implements K, L1, L2, L3 and L23: `olla-dft xanes` rejects them (`xanes.validar_borde`).
- The report warns: *"NO verificado automáticamente: estados fantasma, derivadas logarítmicas y transferibilidad… el cutoff del pseudo anterior NO sirve para este."* The cutoff must be reconverged with `olla-dft converge`.
- With `--projectors 2` the pseudo comes out ultrasoft and *"casi siempre hay que ajustar --rcut a mano hasta que ld1.x converja"*.
- `ld1.x` is not built by default in QE (`make ld1`).

**References.**
- A. Dal Corso, *Comput. Mater. Sci.* 95, 337 (2014) — pslibrary and `ld1.x`.
- N. Troullier, J. L. Martins, *Phys. Rev. B* 43, 1993 (1991) — TM pseudisation (`tm=.true.`).
- C. J. Pickard, F. Mauri, *Phys. Rev. B* 63, 245101 (2001) — GIPAW reconstruction.

---

### `olla-dft xanes` — X-ray absorption near the edge (xspectra.x)

**What it answers.** What is the shape of the XANES/NEXAFS spectrum of a given atom at a given edge, with a core hole and a given polarisation, and how much does it depend on the field direction?

**Background for non-experts.** An X-ray photon knocks an electron out of a deep level (1s at the K edge) and sends it to the empty states. The dipole selection rule only allows final states with angular momentum $l \pm 1$: from 1s one sees the empty $p$ states **of that atom**. The spectrum is, in essence, the density of empty states projected on the absorber, which is why it is local, element-selective and sensitive to oxidation state and coordination. The hole left by the electron attracts the empty states and shifts the edge, so the absorbing atom is described with the core-hole pseudopotential from `corehole`, and since the ejected electron is assumed to leave the system, the cell carries total charge +1 (full core hole approximation, FCH). `xspectra.x` computes the cross-section with the Lanczos method and continued fractions without building the empty states.

**Formulas.** In `qekit/modules/xanes.py`.

Powder average (`xanes.collect`):
$$\sigma(E) = \tfrac{1}{3}\left[\sigma_x(E) + \sigma_y(E) + \sigma_z(E)\right]$$

Minimum distance between absorber images (`xanes.distancia_imagen_minima`):
$$d_{\min} = \min_{(i,j,k)\neq 0,\ |i|,|j|,|k|\le 1}\left|i\,\mathbf{a} + j\,\mathbf{b} + k\,\mathbf{c}\right|$$

Operational onset (`xanes.onset`): first energy at which $\sigma \ge 0.5\,\sigma_{\max}$. Anisotropy (`_anisotropia`): $\max_E[\mathrm{ptp}_i\,\sigma_i(E)]/\max\sigma$; highlighted if $> 0.1$.

**How Olla-DFT computes it.**
1. `xanes.validar_borde` (also in `_cmd_xanes`) normalises `--edge` and only accepts `BORDES_XSPECTRA` = K, L1, L2, L3, L23; M edges are rejected with an explicit message. `BORDE_COREHOLE` says which `--edge` of `corehole` generates the hole for each edge (L2 and L3 share the 2p hole = `L23`). `xanes.prepare` locates the `--element`/`--site` atom, moves it to the **first** position of the list and declares it as a separate species with a three-letter label (`etiqueta_excitada`, e.g. `Sih`; the QE limit is `CHARACTER(LEN=3)`).
2. `sweep.prepare_common` (task `xanes`, excluding the core-hole UPF) and `inputgen.build_pw_input` write `scf.in` with `tot_charge = 1.0`; `_marcar_absorbedor` adds the excited species to `ATOMIC_SPECIES`, changes the label of the first atom and increments `ntyp` by 1.
3. `corehole.core_wfc` extracts `<El>.wfc` from the core-hole UPF (`PP_GIPAW_CORE_ORBITAL` sections).
4. `xanes.build_xspectra_input` writes `xspectra_pol.in` (or `xspectra_x/y/z.in` with `--average`): `calculation='xanes_dipole'`, `edge=--edge`, `xiabs=1`, `xepsilon=--polarization`, `xniter=2000`, `xcheck_conv=10`, `xerror=0.001`; `&plot`: `xnepoint=1000`, `xgamma=--broadening` (0.8 eV), `xemin=-10`, `xemax=30`, `terminator`, `cut_occ_states=.true.`; `&pseudos`: `filecore`, `r_paw(1)=--r-paw` (3.0); `&cut_occ`: `cut_desmooth=0.1`, `cut_stepl=0.01`; k grid at the end.
5. The report measures $d_{\min}$ and warns if it is below `DIST_MINIMA` = 8 Å.
6. The user runs `pw.x -in scf.in` and `xspectra.x -in xspectra_*.in`.
7. `xanes.collect --collect` reads all `xanes_*.dat` (columns E − E_F, σ), averages if there are several, and reads `xgamma` from the *"Broadening parameter (in eV)"* comment in the first file.
8. `report` gives the 50 % onset, main maximum, peaks (> 5 % of the maximum), anisotropy; `export` writes `XANES.dat` and `XANES.txt`; `plot` the figure.

**Where each datum comes from.**

| Datum | Source | Detail |
|---|---|---|
| $\sigma(E)$ per polarisation | `xanes_<dir>.dat` from `xspectra.x` | `xanes._leer_dat` |
| Broadening `xgamma` | header of `xanes_*.dat` | regex *"Broadening parameter"* |
| Core wavefunction | core-hole UPF (`PP_GIPAW_CORE_ORBITAL`) | `corehole.core_wfc` |
| Total charge +1 | fixed in `xanes.prepare` | `tot_charge=1.0` |
| Polarisation | `--polarization` (1 0 0) or the `EJES` axes with `--average` | Cartesian vector (`xcoordcrys=.false.`) |
| k grid | `--kspacing` → `sweep.default_grid` | also in `xspectra.in` |
| $d_{\min}$ | cell vectors | `distancia_imagen_minima` |

**Limits and pitfalls.**
- The energy axis is relative to the Fermi level, not photon energy: *"Para comparar con un experimento se alinea el borde y se compara la FORMA."*
- Supercell warning: *"AVISO: X Å es poco. Con condiciones periódicas el hueco de core ve sus propias imágenes…"* (threshold 8 Å).
- Single polarisation: *"UNA sola polarización. En un cristal anisótropo el espectro depende de la dirección…"*.
- The onset (`xanes.onset`) is the first point where σ reaches 50 % of the **global maximum**: a weak pre-edge before the white line does not count as the onset (the docstring now states this).
- M edges: *"xspectra.x solo calcula bordes K y L (K, L1, L2, L3, L23); los bordes M no están implementados en QE, aunque 'olla-dft corehole' pueda generar el pseudo con ese hueco."*
- `distancia_imagen_minima` only looks at the 26 neighbouring cells: for very oblique cells it can overestimate $d_{\min}$.
- Without `--core-hole` the command aborts: *"falta --core-hole con el UPF de hueco de core. Sin él se calcularía el espectro del estado fundamental…"* and suggests `olla-dft corehole <El> --edge <BORDE_COREHOLE[edge]>`.

**References.**
- M. Taillefumier, D. Cabaret, A.-M. Flank, F. Mauri, *Phys. Rev. B* 66, 195107 (2002) — XSpectra, Lanczos with continued fractions.
- C. Gougoussis, M. Calandra, A. P. Seitsonen, F. Mauri, *Phys. Rev. B* 80, 075102 (2009) — XSpectra with PAW/GIPAW.
- O. Bunău, M. Calandra, *Phys. Rev. B* 87, 205105 (2013) — L₂,₃ edges.

---

### `olla-dft xps` — Initial-state core-level shifts (initial_state.x)

**What it answers.** By how much does the core-level energy of each atom shift relative to the others of its species? It is the theoretical counterpart of the chemical shift in an XPS spectrum.

**Background for non-experts.** XPS measures the energy needed to eject a core electron. An atom surrounded by electronegative neighbours has its core more tightly bound (positive shift) than one in a metallic environment. The **initial-state** approximation computes only how the potential felt by the core electron changes *before* it is removed; it ignores the relaxation of the other electrons around the hole (the *final state*), which can amount to several tenths of an eV. That is why what comes out are **relative** shifts between sites, not absolute binding energies. `initial_state.x` needs two species of the same element in the input — the normal one and one with a core hole — because it defines the shift from `delta_zv = zv(excited) − zv(normal)`; if both are the same it returns zeros without warning.

**Formulas.** In `qekit/modules/xps.py`. The shift is computed by `initial_state.x`; Olla-DFT only reads and rearranges it:

$$\Delta_i = \text{shift}_i^{\mathrm{TOTAL}},\qquad \Delta_i^{\mathrm{rel}} = \Delta_i - \min_j \Delta_j,\qquad \text{spread} = \max_i\Delta_i - \min_i\Delta_i$$

Cancellation indicator (`xps.report`):
$$\frac{\max_{c}\,\mathrm{ptp}(\text{contribution}_c)}{\text{spread}} > 20 \Rightarrow \text{numerical-cancellation warning}$$

- $\Delta_i$: shift of atom $i$ in eV, read from the line `atom i type t shift = … Ry, = … eV` of the *TOTAL* section. The contributions $c$ (Fermi, local, non-local, ionic, core-correction, Hubbard…) are read from the *"The X contribution to shift"* sections.

**How Olla-DFT computes it.**
1. `xps.prepare` reads `--core-hole EL=file.UPF` (repeatable). For each element: `_verificar_par` requires the normal and the core-hole UPF to be different files and `z_valence` to differ by exactly +1 (`qekit.core.pseudo.z_valence`).
2. `inputgen.build_pw_input` writes `scf.in` with the extra species (`extra_species`) declared in `ATOMIC_SPECIES` **without** any atom using them; `_copiar_pseudos` copies the core-hole UPF into `pseudo_dir`.
3. `xps.build_input` writes `initial_state.in` with `excite(t_normal) = t_hole` (1-based indices in `ATOMIC_SPECIES` order); `excite(t)=t` is rejected.
4. `structure.symmetry_dataset` counts orbits of equivalent atoms; if there is only one it warns that everything will come out zero.
5. The user runs `pw.x -in scf.in` and `initial_state.x -in initial_state.in > initial_state.out`.
6. `xps.collect --collect` parses `initial_state.out` with `_RE_SECCION` and `_RE_ATOMO`, takes the eV column, and sets `equivalentes=True` if all $|\Delta_i| < 10^{-6}$ eV.
7. `report` tabulates shifts, shift relative to the minimum, spread, decomposition per contribution and the cancellation warning; `export` writes `XPS_CORE.dat`.

**Where each datum comes from.**

| Datum | Source | Detail |
|---|---|---|
| Shift per atom and contribution | `initial_state.out` | `xps.collect`, regex `atom N type T shift = X Ry, = Y eV` |
| `z_valence` normal and core-hole | UPF headers | `pseudo.z_valence` in `_verificar_par` |
| Inequivalent sites | spglib via `structure.symmetry_dataset` | `equivalent_atoms` |
| Excited species label | `xanes.etiqueta_excitada` | 3 characters |
| Symbols per atom | input structure | `atoms.get_chemical_symbols()` |

**Limits and pitfalls.**
- **Initial state only.** There is no ΔSCF and no final state; the module docstring now states it explicitly: the core-hole UPF is used *only as the "excited species" that initial_state.x needs to define the shift, not to relax the system around the hole*. The report refers you on: *"las energías de enlace absolutas necesitan un ΔSCF con hueco de core."*
- Spread < 0.1 eV: *"Por debajo de ~0.1 eV el corrimiento no es concluyente: la relajación de estado final… es del mismo orden."*
- Without `--core-hole` only `scf.in` is written and the report explains that `initial_state.x` would return zeros.
- All atoms equivalent: *"AVISO: todos los átomos son equivalentes por simetría, así que todos los corrimientos van a salir exactamente cero."*
- Large cancellation: *"CUIDADO con la cancelacion… baja conv_thr (1e-10 o menos) y sube la malla k antes de creerte la tercera cifra."*
- The error messages point to `olla-dft corehole <El> --edge K` to generate the consistent pair.

**References.**
- E. Pehlke, M. Scheffler, *Phys. Rev. Lett.* 71, 2338 (1993) — initial vs final state in core-level shifts.
- L. Köhler, G. Kresse, *Phys. Rev. B* 70, 165405 (2004) — core-level binding energies with a core hole.
- `initial_state.x` documentation (Quantum ESPRESSO, PP package).

---

### `olla-dft charges` — Löwdin charges, on-grid Bader and density difference

**What it answers.** How much electronic charge "belongs" to each atom, and where does the density accumulate or deplete when a bond or an adsorption forms?

**Background for non-experts.** The electron density is continuous; sharing it among atoms requires a rule. **Löwdin** projects the states onto orthogonalised atomic orbitals (done by `projwfc.x`); it is cheap and depends on the orbital basis of the pseudopotential. **Bader** uses no orbitals: it divides space into "basins" by following the steepest ascent of the density from each point to a maximum, like rainwater running down slopes into each valley, but reversed. The **density difference** $\rho_{AB} - \rho_A - \rho_B$ shows, point by point, what changed when the two parts were joined.

**Formulas.** In `qekit/modules/charges.py`.

Löwdin (`charges.read_lowdin`, `report_lowdin`):
$$q_i^{\mathrm{net}} = Z_i^{\mathrm{val}} - Q_i^{\mathrm{Löwdin}}$$

On-grid Bader (`charges.bader`): for each grid point the neighbour $\nu$ (out of 26) that maximises the slope
$$s_\nu = \frac{\rho(\mathbf{r}+\mathbf{d}_\nu) - \rho(\mathbf{r})}{|\mathbf{d}_\nu|}$$
is chosen and the chain is followed up to a local maximum (path compression, max. 64 iterations). Each maximum is assigned to the nearest atom with periodic images. Then
$$Q_i = \sum_{\mathbf{r}\in\Omega_i}\rho(\mathbf{r})\,\Delta V,\qquad V_i = N_i\,\Delta V_{\mathrm{Å}^3},\qquad \Delta V_{\mathrm{Å}^3} = \frac{V_{\mathrm{cell}}}{n_1 n_2 n_3},\quad \Delta V = \frac{\Delta V_{\mathrm{Å}^3}}{a_0^3}$$
- $\rho$: density from the `.cube` in e/bohr³ (what `pp.x` writes, `density_units="e/bohr3"` by default; `"e/A3"` is also accepted); `charges._voxel_volume` returns the voxel volume in the units of the density (bohr³, with $a_0$ = `fields.BOHR` = 0.529177210903 Å) so that $\rho\,\Delta V$ is a number of electrons, and in Å³ to report the basin volumes.

Density difference (`charges.difference`, `report_difference`), with the same $\Delta V$ in bohr³:
$$\Delta\rho = \rho_{\mathrm{total}} - \sum_p \rho_p,\qquad Q_{\mathrm{net}} = \sum \Delta\rho\,\Delta V,\qquad Q_{\mathrm{acc}} = \sum_{\Delta\rho>0}\Delta\rho\,\Delta V$$

**How Olla-DFT computes it.**
1. If the structure is given, `charges.valence_from_pseudos` reads `z_valence` from the UPFs in `--pseudo-dir` (or the configured `pseudo_dir`) via `pseudo.resolve`; if any UPF cannot be read it returns `None`, the CLI warns (*"no pude leer z_valence de los UPF…"*) and the "neta" column stays `n/d`.
2. `--lowdin projwfc.out`: `charges.read_lowdin` looks for `Atom #  i: total charge = q` and `Spilling Parameter:`; with the structure it adds symbols and the net charge $Z^{\mathrm{val}} - Q$.
3. `--bader density.cube` (needs the structure): `fields.read_cube` reads the cube (`plot_num=0` from `pp.x`), `charges.bader` partitions and compares the sum of basins with the total integral; `report_bader` also compares the integral with $\sum_i Z_i^{\mathrm{val}}$ and warns if they differ by more than 5 %.
4. `--difference total.cube part1.cube …`: `charges.difference` requires identical grids and subtracts; `report_difference` gives net charge, accumulated charge and the extrema of the planar profile (`fields.planar_average`, axis `--axis`); `plot_difference` draws the profile.

**Where each datum comes from.**

| Datum | Source | Detail |
|---|---|---|
| Löwdin charges, spilling | `projwfc.x` output | regex `_RE_LOWDIN`, `_RE_SPILL` |
| $\rho(\mathbf{r})$ | `.cube` from `pp.x` (`plot_num=0`, `output_format=6`) | `fields.read_cube` |
| Atomic positions (Bader) | structure `file` | `atoms.positions` |
| $Z^{\mathrm{val}}$ per atom | `z_valence` from the UPFs (`--pseudo-dir` or configuration) | `charges.valence_from_pseudos` → `pseudo.resolve` |
| $a_0$ (bohr → Å) | `fields.BOHR` | 0.529177210903 |
| Profile axis | `--axis` (0/1/2) | `fields.planar_average` |

**Limits and pitfalls.**
- On-grid Bader: *"Hereda el sesgo de malla del método (centésimas de electrón); para números finos usa la variante near-grid del código `bader` de Henkelman."* Warning if the sum of basins differs from the integral by more than 1e-3 e: *"la malla del cube es demasiado gruesa."*
- If the grid integral does not match $\sum Z^{\mathrm{val}}$ (by more than 5 %): *"Revisa que el cube sea la densidad de valencia completa (plot_num=0) y que los UPF de --pseudo-dir sean los del cálculo."* A cube already in e/Å³ must be declared with `density_units="e/A3"` (Python only; the CLI assumes e/bohr³).
- Löwdin: spilling > 0.05 → *"AVISO: por encima de ~0.05 la base atómica no describe bien los estados"*. Useful to compare atoms, not as an absolute charge.
- Without a readable `--pseudo-dir` the "neta" column is `n/d` with a warning; the UPFs must be those of the calculation, because $Z^{\mathrm{val}}$ depends on the pseudo (semicore or not).
- The $\Delta\rho$ profile is plotted in e/bohr³ (cube units), not e/Å³.
- `--difference` requires the same cell, FFT grid and cutoffs: *"las rejillas no coinciden… la resta no significa nada."*

**References.**
- R. F. W. Bader, *Atoms in Molecules: A Quantum Theory* (Oxford, 1990).
- G. Henkelman, A. Arnaldsson, H. Jónsson, *Comput. Mater. Sci.* 36, 354 (2006) — on-grid Bader.
- P.-O. Löwdin, *J. Chem. Phys.* 18, 365 (1950).

---

### `olla-dft charge` — pp.x scalar fields and planar profile

**What it answers.** How are the charge density, the spin density, the ELF or the electrostatic potential of a finished calculation distributed along an axis?

**Background for non-experts.** `pp.x` extracts from the already computed wavefunctions and density a scalar field on the 3D grid. Averaging it over the planes perpendicular to an axis gives a 1D "profile" that is easy to read: where the layers of a slab are, where spin accumulates, where the vacuum is.

**Formulas.** `fields.planar_average`:
$$\bar f(z_k) = \frac{1}{n_1 n_2}\sum_{i,j} f(i,j,k),\qquad z_k = k\,|\mathbf{h}_3|$$
- $\mathbf{h}_3$: grid step along the chosen axis (Å). The other axes are obtained by permutation.

**How Olla-DFT computes it.**
1. `_cmd_charge`: if `<name>.cube` does not exist (or with `--rerun`), `fields.run_pp` writes `pp_<field>.in` with the `plot_num` from `fields.PLOTS` (density 0, vtotal 1, spin 6, elf 8, potential 11), `iflag=3`, `output_format=6`, and runs `pp.x` (looked up next to `pw.x`); requires `JOB DONE`.
2. `fields.read_cube` reads origin, axes (bohr → Å if $n>0$) and values.
3. `fields.planar_average` along `--axis` (a/b/c); `PERFIL_PLANAR.dat` and the figure `perfil_<name>` are written.

**Where each datum comes from.**

| Datum | Source | Detail |
|---|---|---|
| 3D field | `<name>.cube` from `pp.x` | `pp.x` units: e/bohr³ (density), Ry (potentials) |
| `prefix` | XML of the calculation | `qeout.read_xml(...).prefix` |
| `plot_num` | table `fields.PLOTS` | 0, 1, 6, 8, 11 |
| Bohr → Å | `qeout.BOHR_ANG` | 0.529177210903 |

**Limits and pitfalls.**
- The profile is exported in the raw cube units (no Ry → eV conversion here; `wf` does it).
- Needs `pp.x` compiled (`make pp`); if missing: *"no se encontró pp.x junto a pw.x…"*.
- The command does not interpret the field: it only averages and plots it. The `.cube` opens in VESTA for isosurfaces.

**References.** `pp.x` documentation (INPUT_PP, Quantum ESPRESSO).

---

### `olla-dft wf` — Work function from the vacuum level

**What it answers.** How much energy does it take to remove an electron from a surface into vacuum? $\Phi = V_{\mathrm{vac}} - E_F$.

**Background for non-experts.** In a slab with vacuum, the electrostatic potential flattens far from the material: that plateau is the "vacuum level", the energy of an electron at rest outside the solid. The work function is the distance from the Fermi level (the last occupied level) to that plateau. If the plateau is not flat, either the vacuum is short or the slab has a net dipole that tilts the potential.

**Formulas.** `fields.work_function`:
$$\bar V(z) = \mathrm{RY\_EV}\cdot\overline{V_{\mathrm{pp}}}(z),\qquad V_{\mathrm{vac}} = \frac{1}{2h+1}\sum_{k=-h}^{h}\bar V\big(z_{i^\ast + k}\big),\qquad \Phi = V_{\mathrm{vac}} - E_F$$
$$\text{flatness} = \max_{k}\bar V - \min_{k}\bar V\ \text{in the same window}$$
- The index window $\{i^\ast + k\}$ is given by `fields.vacuum_window` when the atomic positions are known (the CLI passes them from the XML): it is the central 20 % of the widest **atom-free** gap along the axis (measured periodically in fractional coordinates), with $h = \max(2, 0.1\,f_{\mathrm{gap}} N_z)$. Without positions it falls back to the blind criterion: $i^\ast = \arg\max_z \bar V$ and $h = \max(2, N_z/10)$ (±10 % of the cell around the maximum). $E_F$ in eV from the XML; `RY_EV` = 13.605693122994.

**How Olla-DFT computes it.**
1. `_cmd_wf`: if `potencial.cube` does not exist, `fields.run_pp(path, "potential", ...)` runs `pp.x` with `plot_num=11` ($V_{\mathrm{bare}} + V_H$).
2. `fields.read_cube` and `qeout.read_xml` (for `fermi`, from the `fermi_energy` tag in Ha → eV).
3. `fields.work_function(cube, E_F, axis, positions=qe.positions)` averages in the plane, locates the vacuum plateau (`vacuum_window`) and computes $\Phi$ and the flatness; the report states over which $z$ range it was evaluated.
4. `report_wf`, `export_wf` (`WF.dat` with header `Phi_eV`, `V_vacio_eV`, `E_Fermi_eV`, `planitud_eV` and the profile) and `plot_profile`.

**Where each datum comes from.**

| Datum | Source | Detail |
|---|---|---|
| $V(\mathbf{r})$ | `potencial.cube` (`pp.x`, `plot_num=11`, Ry) | `fields.read_cube` |
| $E_F$ | `pw.x` XML, tag `fermi_energy` | `qeout.read_xml`, Ha → eV |
| Ry → eV | `qeout.RY_EV` | 13.605693122994 |
| Atomic positions (for `vacuum_window`) | `pw.x` XML (`atomic_positions`) | `qeout.read_xml(...).positions` |
| Axis | `--axis` (c by default) | `_AXES` |

**Limits and pitfalls.**
- Warning if flatness > 0.05 eV: *"la meseta de vacío varía más de 0.05 eV. El vacío es insuficiente o hay un dipolo neto; aumenta el vacío (o usa una losa simétrica)…"*.
- Without positions (use from Python) the plateau is searched blindly around the potential maximum: *"con poco vacío la ventana puede pisar la cola del potencial atómico"* (docstring). The CLI always passes the positions from the XML.
- It does not apply a dipole correction by itself: a polar slab gives two different vacuum levels and this command takes the higher one. For polar slabs the calculation must be generated with `--dipole` (`gen`, `eform`) or `esm` must be used.
- If the XML has no `fermi_energy` (fixed occupations): *"el XML no trae energía de Fermi (¿terminó el scf?)"*.

**References.**
- N. D. Lang, W. Kohn, *Phys. Rev. B* 3, 1215 (1971) — work function in the jellium model.
- L. Bengtsson, *Phys. Rev. B* 59, 12301 (1999) — dipole correction in slabs.

---

### `olla-dft esm` — Charged surfaces with the effective screening medium

**What it answers.** What are the work function, the capacitance and the potential of zero charge of a slab (neutral or charged) without periodic images or the compensating background contaminating the result?

**Background for non-experts.** A charged slab in a periodic cell is an ill-posed problem: QE spreads a uniform background of opposite charge over the whole volume, vacuum included, and the energy depends on the cell size without converging to anything. The **ESM** (Effective Screening Medium) replaces periodicity along $z$ by an explicit boundary condition: the Poisson equation is solved inside the cell and matched to an analytic solution outside. Three variants: `bc1` (vacuum on both sides, neutral slabs; the vacuum level is zero by construction), `bc2` (two metal plates: a capacitor, admits a field) and `bc3` (vacuum/metal: an electrode that receives the counter-charge). With `bc2`/`bc3` the distance to the electrode is no longer a convergence parameter but **physics**: it fixes the capacitance.

**Formulas.** In `qekit/modules/esm.py`.

Centring (`esm.centrar`): $z_i \leftarrow z_i - \tfrac{1}{2}(z_{\min}+z_{\max})$ (ESM measures $z$ from the cell centre).

Vacuum level (`esm.nivel_vacio`): average of $V_{\mathrm{tot}}(z)$ from the `.esm1` in the region $|z| > t/2 + m$, with $t$ the slab thickness and a margin $m$ that starts at `MARGEN_VACIO` = 2 Å and grows in 0.5 Å steps (up to `margen_max` = 8 Å) until the standard deviation of the potential drops below `tol` = 1e-3 eV; with `bc3` only the $z<0$ side.

$$\Phi = V_{\mathrm{vac}} - E_F$$

Capacitance (`esm.capacitancia`), linear fit $q = C' V + b$:
$$C = \frac{dq}{dV}\,\frac{1}{A}\cdot 1.602176634\times10^{3}\quad[\mu\mathrm{F/cm^2}],\qquad R^2 = 1 - \frac{\sum(q-\hat q)^2}{\sum(q-\bar q)^2}$$
- $q$ in e per cell, $V$ in V (eV/e), $A$ = cell area in Å² (`|(\mathbf a\times\mathbf b)_z|`); `E_A2_A_UF_CM2` = $1.602176634\times10^{3}$ converts e/(Å²·V) to µF/cm².

Linearity (`esm.linealidad`): $\max|P - \hat P| / (\max P - \min P) \le$ `tol` = 0.02.

Potential of zero charge (`esm.potencial_de_carga_cero`): linear interpolation of $\Phi(q)$ at $q = 0$.

Grand canonical (`esm.gran_canonico`, library only): $\Omega = E + q\,\Phi$.

**How Olla-DFT computes it.**
1. `esm.comprobar`: rejects `bc1` with charge (*"bc1 es vacío por los dos lados… la energía diverge"*) and cells not orthogonal in $z$; warns if vacuum < `VACIO_MINIMO` = 6 Å, if the slab was not centred and, with `bc2/bc3` and charge, that the vacuum is physics.
2. `esm.prepare` centres the slab, computes thickness, vacuum and area, and writes one `scf` per charge in `q00/`, `q01/`… (`inputgen.build_pw_input`, `conv_thr=1e-8`, `mv` smearing with `degauss=0.02`, grid $n_1\times n_2\times 1$, `tot_charge=q`) and inserts into `&SYSTEM`: `assume_isolated='esm'`, `esm_bc`, `esm_nfit=--nfit` (4), `esm_w=--esm-w` if ≠ 0, `esm_efield=--field` only with `bc2`. Writes `run.sh`.
3. `--run` or by hand: `pw.x` in each folder.
4. `esm.collect` reads from each folder the XML (`total_energy`, `fermi`) and the `<prefix>.esm1` (`esm.leer_esm1`: z (Å), charge (e/Å), $V_H$, $V_{\mathrm{loc}}$, $V_{\mathrm{tot}}$ in eV); `nivel_vacio` and $\Phi$.
5. `esm.report`: table $q, E, E_F, V_{\mathrm{vac}}, \Phi$; with `bc1` checks $|V_{\mathrm{vac}}| < 10^{-3}$ eV; with several charges, capacitance from $V_{\mathrm{vac}}(q)$ (cell voltage) and, if $\Phi(q)$ is linear, also from $\Phi(q)$ with the PZC.
6. `export` (`ESM.dat`, `ESM_perfil_qNN.dat`, `ESM.txt`) and `plot` (profiles and $q$ vs $\Phi$).

**Where each datum comes from.**

| Datum | Source | Detail |
|---|---|---|
| $V_{\mathrm{tot}}(z)$, charge$(z)$ | `<prefix>.esm1` written by `pw.x` with ESM | `esm.leer_esm1`, columns 0–4 |
| $E$, $E_F$ | `pw.x` XML | `qeout.read_xml` (Ha → eV) |
| Area $A$ | cell vectors $\mathbf a,\mathbf b$ | `esm.prepare` |
| µF/cm² factor | `esm.E_A2_A_UF_CM2` | $e/(10^{-8}\,\mathrm{cm})^2$ |
| Charges | `--charge` (list) | e per cell |
| Field | `--field` (Ry/a.u.) | `bc2` only |

**Limits and pitfalls.**
- *"Con bc2 o bc3 la capacitancia depende de la distancia al contraelectrodo: es una capacitancia DE ESTE MONTAJE, no una propiedad del material."*
- Energies with net charge are not comparable with each other: *"la energía de ESM incluye la interacción con la carga imagen del electrodo, que crece como q²."*
- If $\Phi(q)$ is not linear: *"Φ(q) = V_vac − E_F NO es una recta… no doy un potencial de carga cero sobre ella."*
- `gran_canonico` (Ω = E + qΦ) exists in the module but **no command uses it**; the "grand canonical" of the module title is not exposed in the CLI.
- The slab is centred automatically; if the user had already centred it at $c/2$ (ASE) the warning explains why it was re-centred.
- The calculation always uses smearing (`insulator=False`): intended for metals/electrodes.

**References.**
- M. Otani, O. Sugino, *Phys. Rev. B* 73, 115407 (2006) — ESM.
- N. Bonnet, T. Morishita, O. Sugino, M. Otani, *Phys. Rev. Lett.* 109, 266101 (2012) — constant potential with ESM.

---

### `olla-dft echem` — Computational hydrogen electrode: HER and OER

**What it answers.** What potential must be applied so that all the steps of hydrogen evolution (HER) or oxygen evolution (OER) become downhill, and how far is it from the equilibrium potential (overpotential)?

**Background for non-experts.** Computing a solvated proton is a very hard problem. The trick of the computational hydrogen electrode (CHE) is to notice that, at 0 V versus the standard hydrogen electrode and pH 0, the pair $\mathrm{H^+ + e^-}$ has the same free energy as $\tfrac12\mathrm{H_2(g)}$, which can be computed. Every step that releases a $(\mathrm{H^+ + e^-})$ is evaluated that way, and the potential $U$ and the pH enter afterwards as additive terms. The step with the largest $\Delta G$ is the "limiting" one: the potential that makes it exergonic is the limiting potential, and its distance to the equilibrium one is the overpotential. This is thermodynamics of intermediates: no kinetic barriers and no solvent.

**Formulas.** In `qekit/modules/echem.py`.

Dependence on $U$ and pH (`Echem.dG`):
$$\Delta G_i(U, \mathrm{pH}) = \Delta G_i(0,0) - eU - k_B T\ln 10\cdot\mathrm{pH} = \Delta G_i(0,0) - e\,U_{\mathrm{RHE}}$$
$$U_{\mathrm{RHE}} = U_{\mathrm{SHE}} + k_B T\ln 10\cdot\mathrm{pH}\quad(\text{`echem.u_rhe`; } 0.0592\,\mathrm{pH\ V\ at\ 298\ K})$$
- $k_B$ = `KB_EV` = $8.617333262\times10^{-5}$ eV/K; $T$ = `--temperature` (298.15 K); $U$ = `-U` in V **versus SHE** (at pH 0 it coincides with RHE); the pH term is exactly the SHE → RHE conversion, so on the RHE scale the $\Delta G$ do not depend on pH. One electron per step.

HER (`echem.her`):
$$\Delta G_{\mathrm{H^*}} = E_{\mathrm{ads}}(\mathrm{H}) + c_{\mathrm{H}},\qquad \text{steps: } (+\Delta G_{\mathrm{H^*}},\ -\Delta G_{\mathrm{H^*}})$$
- $E_{\mathrm{ads}}(\mathrm{H})$: `--her`, referred to $\tfrac12\mathrm{H_2}$ (eV); $c_{\mathrm{H}}$ = ZPE − TΔS = 0.24 eV by default (`CORRECCIONES`).

OER (`echem.oer`), with $G_X = E_{\mathrm{ads}}(X) + c_X$:
$$\Delta G_1 = G_{\mathrm{OH}},\quad \Delta G_2 = G_{\mathrm{O}} - G_{\mathrm{OH}},\quad \Delta G_3 = G_{\mathrm{OOH}} - G_{\mathrm{O}},\quad \Delta G_4 = 4.92\ \mathrm{eV} - (\Delta G_1+\Delta G_2+\Delta G_3)$$
- $c_{\mathrm{OH}} = 0.35$, $c_{\mathrm{O}} = 0.05$, $c_{\mathrm{OOH}} = 0.40$ eV by default; `DG_AGUA_TOTAL` = 4.92 eV (experimental, $2\mathrm{H_2O} \to \mathrm{O_2} + 2\mathrm{H_2}$).

Limiting potential and overpotential (`Echem.U_limitante`, `Echem.sobrepotencial`):
$$U_L = \max_i \Delta G_i(0,0)/e,\qquad \eta = U_L - U_{\mathrm{eq}},\quad U_{\mathrm{eq}}^{\mathrm{OER}} = 1.229\ \mathrm{V},\ U_{\mathrm{eq}}^{\mathrm{HER}} = 0$$
- $\eta$ is returned **with sign**: positive = at $U_{\mathrm{eq}}$ the limiting step is still uphill (with the profiles built here it never comes out negative; it only could with a `dG_total` different from the experimental one).

Scaling relation (`echem.escala_ooh_oh`, OER) and its limit (`echem.sobrepotencial_minimo_escala`):
$$\Delta_{\mathrm{sc}} = G_{\mathrm{OOH}} - G_{\mathrm{OH}}\ \text{(compared with `ESCALA_OOH_OH` = 3.2 ± 0.2 eV)},\qquad \eta_{\min} = \frac{\Delta_{\mathrm{sc}}}{2} - \frac{\Delta G_{\mathrm{total}}}{4} = 0.37\ \mathrm{V}$$
- If OOH* and OH* are separated by a fixed $\Delta_{\mathrm{sc}}$, steps 2 and 3 add up to $\Delta_{\mathrm{sc}}$ and the worse one cannot drop below $\Delta_{\mathrm{sc}}/2$ = 1.6 eV; against 4.92/4 = 1.23 V that leaves ~0.37 V.

Pourbaix-like grid (`echem.pourbaix`, library only): $\Delta G_{\lim}(U,\mathrm{pH}) = \max_i\Delta G_i(0,0) - eU - k_BT\ln10\cdot\mathrm{pH}$ over $U\in[-0.5,2]$ V and pH $\in[0,14]$.

**How Olla-DFT computes it.**
1. `_cmd_echem` requires exactly one of `--her E` or `--oer OH=..,O=..,OOH=..`; `--corrections X=eV` overrides the thermal corrections.
2. `echem.her` or `echem.oer` build the list of steps $(\text{name}, \Delta G_i)$; `oer` warns if $\Delta G_4 < 0$ and if tabulated corrections were used.
3. The user's `U` (vs SHE) and `pH` are set; `echem.report` also prints $U_{\mathrm{RHE}}$ (`Echem.U_rhe`) if pH ≠ 0, tabulates $\Delta G(0)$ and $\Delta G(U,\mathrm{pH})$, the limiting step, $U_L$ (vs RHE), signed $\eta$, the descriptor $\Delta G_{\mathrm{H^*}}$ (HER) or the scaling relation and the $\eta_{\min}$ it imposes (OER).
4. `export` writes `ECHEM.dat` and `ECHEM.txt`; `plot` draws the staircase diagram at $U = 0$, $U_{\mathrm{eq}}$ and $U_L$.

**Where each datum comes from.**

| Datum | Source | Detail |
|---|---|---|
| $E_{\mathrm{ads}}$ of H, OH, O, OOH | parameters `--her`, `--oer` | eV, referred to H₂O and ½H₂ (from `adsorb`) |
| ZPE − TΔS corrections | `echem.CORRECCIONES` or `--corrections` | H 0.24, OH 0.35, O 0.05, OOH 0.40 eV (Nørskov et al.) |
| Total $\Delta G$ of water | `echem.DG_AGUA_TOTAL` | 4.92 eV, experimental |
| $U_{\mathrm{eq}}$ | `echem.U_EQ_OER`, `U_EQ_HER` | 1.229 V, 0 V |
| $k_B$ | `echem.KB_EV` | $8.617333262\times10^{-5}$ eV/K (CODATA) |
| Universal $\Delta_{\mathrm{sc}}$ | `echem.ESCALA_OOH_OH` | 3.2 eV (Man et al. 2011) |

**Limits and pitfalls.**
- *"El CHE es termodinámica de intermedios: NO hay barreras cinéticas, ni disolvente explícito, ni doble capa."*
- `-U` is versus the **SHE** (CLI help: *"a pH 0 es el mismo que frente al RHE; el pH lo convierte"*); $U_L$ and $\eta$ are on the RHE scale. For the HER $U_L = |\Delta G_{\mathrm{H^*}}|$, so $\eta \ge 0$ always.
- Fourth step by difference: *"El cuarto paso sale NEGATIVO… o hay un error en las referencias, o tu superficie liga los intermedios muchísimo."*
- `pourbaix()` is not wired to any command: the "Pourbaix diagram" of the module title is not produced from the CLI.

**References.**
- J. K. Nørskov, J. Rossmeisl, A. Logadottir, L. Lindqvist, J. R. Kitchin, T. Bligaard, H. Jónsson, *J. Phys. Chem. B* 108, 17886 (2004) — CHE. DOI: 10.1021/jp047349j.
- J. K. Nørskov, T. Bligaard, A. Logadottir, J. R. Kitchin, J. G. Chen, S. Pandelov, U. Stimming, *J. Electrochem. Soc.* 152, J23 (2005) — HER volcano.
- I. C. Man et al., *ChemCatChem* 3, 1159 (2011) — OER scaling relation.

---

### `olla-dft adsorb` — Adsorption sites and adsorption energy

**What it answers.** On which inequivalent sites of a surface can a molecule sit, and how much does the system gain (or lose) by doing so on each?

**Background for non-experts.** A molecule on a metal sits on top of an atom (*top*), over the midpoint between two (*bridge*) or over the centre of a triangle of atoms (*hollow*; on fcc(111) there are two: with or without an atom underneath in the second layer). Many of those sites are copies by symmetry, so they are grouped by their "fingerprint": the sorted list of distances to their 24 nearest neighbours counting all layers. The adsorption energy is a subtraction of three total energies that only makes sense if the three calculations share cell, cutoffs, k grid and pseudos; that is why they are generated together.

**Formulas.** `thermochem.adsorcion` (called from `AdsorbRun.energias_ads`):
$$E_{\mathrm{ads}} = E(\text{slab}+\text{mol}) - E(\text{slab}) - n\,E(\text{mol})$$
- All in eV; $n$ = number of molecules (`n_mol`, 1). Negative = favourable.

Geometry after relaxation (`adsorb.collect`):
$$h = \min_{a\in\mathrm{ads}} z_a - \max_{s\in\mathrm{slab}} z_s,\qquad d_{\mathrm{contact}} = \min_{a,s}|\mathbf r_a - \mathbf r_s|$$

Site fingerprint (`adsorb._huella`): sorted distances to the $k$ = `N_VECINOS_HUELLA` = 24 nearest atoms (with periodic replicas); two sites are the same if $\max|\Delta d| <$ `TOL_HUELLA` = 0.05 Å.

**How Olla-DFT computes it.**
1. `adsorb.prepare` requires vacuum along $c$ (`kpoints.direcciones_con_vacio`) and loads the molecule (`cargar_molecula`: file or ASE G2 database).
2. `adsorb.sitios`: exposed layer = atoms within `TOL_CAPA` = 0.6 Å of the extreme $z$; *top* over each of them; *bridge* between pairs closer than `R_VECINO` = 3.6 Å; *hollow* at the centroids of the Delaunay triangulation (triangles with a side > 1.6·3.6 Å are discarded); they are brought into the cell and deduplicated by fingerprint; labelled `top1`, `bridge1`, `hollow1`…
3. With `--rotations N` and a polyatomic molecule, each site is repeated with rotations of $360k/N$ degrees around $z$.
4. `sweep.prepare_common` is resolved on the **union** slab + molecule (same pseudos and cutoffs for everything). `_losa/`, `_molecula/` (molecule centred in the **same** cell) and one folder per site are written (`adsorb.colocar`: atom `--anchor` at `--height` = 2.0 Å above the site), all `relax` unless `--fixed-ions`; `run.sh`. With `--dipole`, `dipole_correction=3` enters the **three** calculations (`inputgen.build_pw_input`: `tefield`, `dipfield`, `edir=3`, `emaxpos`/`eopreg` at the centre of the vacuum gap via `inputgen._region_vacio`, `eamp=0`); without ≥ 5 Å of vacuum it aborts.
5. `--run`/`--collect`: `adsorb.collect` reads the XMLs (`qeout.read_xml`), energies, convergence, height and contact.
6. `adsorb.report`: table sorted by $E_{\mathrm{ads}}$, best site, diagnosis by ranges (>0: does not bind; > −0.30 eV: weak physisorption; < −2 eV: probable reaction/dissociation or atomic chemisorption), difference with the second (< 50 meV: indistinguishable). `export`: `ADSORCION.dat/.txt`; `plot`: bars.

**Where each datum comes from.**

| Datum | Source | Detail |
|---|---|---|
| $E$(slab), $E$(mol), $E$(slab+mol) | `pw.x` XML in `_losa/`, `_molecula/`, `<site>/` | `total_energy` (Ha → eV) |
| Relaxed positions | XML (`atomic_positions`) | height and contact |
| Molecule | file or `ase.build.molecule` | `--mol` |
| Neighbour radius, tolerances | `adsorb.R_VECINO`, `TOL_CAPA`, `N_VECINOS_HUELLA`, `TOL_HUELLA` | 3.6 Å, 0.6 Å, 24, 0.05 Å |
| Initial height | `--height` | 2.0 Å |
| vdW correction | `--vdw` | passed to `inputgen.build_pw_input` |
| Dipole correction | `--dipole` | `dipole_correction=3` in slab, molecule and slab+molecule |

**Limits and pitfalls.**
- Without `--vdw`: *"AVISO: sin corrección de van der Waals. En fisisorción… la energía sale cerca de cero y la geometría desligada."*
- Without `--dipole` on the `top` face: *"Sugerencia: una molécula adsorbida en una sola cara deja la losa polar. Con --dipole se cancela el dipolo artificial a través del vacío."* The sawtooth is put into all three calculations on purpose: *"si la referencia se calcula sin corregir, la resta arrastra el error."*
- $E_{\mathrm{ads}} > 0$ with fixed ions: *"lo más probable es que la altura inicial… no sea la de equilibrio y estés midiendo la repulsión."*
- The reference is the molecule exactly as given: with `--mol H` the reference is the **atom**, not ½H₂ (the report warns about it for $|E_{\mathrm{ads}}| > 2$ eV).
- The isolated molecule is computed with the same k grid as the slab (deliberate consistency, not a separate box).
- Site enumeration is geometric: it does not detect sites over second layers or reconstructions.

**References.**
- B. Hammer, J. K. Nørskov, *Adv. Catal.* 45, 71 (2000) — adsorption on metal surfaces.
- S. Grimme, J. Antony, S. Ehrlich, H. Krieg, *J. Chem. Phys.* 132, 154104 (2010) — DFT-D3.

---

### `olla-dft surface` — Cutting an (hkl) slab with vacuum

**What it answers.** Given a crystal, what does the $(hkl)$ surface slab with $N$ layers and vacuum look like, is it symmetric, is it polar, and how much real vacuum is left?

**Background for non-experts.** A surface is simulated with a "slab": a few atomic layers parallel to the $(hkl)$ plane and, above them, enough vacuum so that the slab does not see its periodic copy. If the two faces are not the same (*polar* slab), an artificial dipole appears across the vacuum and shifts the work functions; QE corrects it with `dipfield`. The vacuum that matters is the one between atoms, not between cell edges.

**Formulas.** `builder.surface`:
$$t = z_{\max} - z_{\min},\qquad v_{\mathrm{real}} = c - t$$
- Symmetric: the sorted profile $z_i - \bar z$ coincides with its mirror within `tol` = 0.3 Å. Polar: composition of the top layer ≠ that of the bottom layer (atoms within `tol` of the extreme).

**How Olla-DFT computes it.**
1. `structure.conventional` → `ase.build.surface(base, miller, layers, vacuum=vacuum/2, periodic=True)` and `slab.center(vacuum=vacuum/2, axis=2)`.
2. `builder.surface` computes thickness, real vacuum, number of atomic planes (`_planos_z`, tolerance 0.3 Å), symmetry and polarity; with `--fix N` it marks the atoms of the $N$ lowest planes in two ways (`_fijar_capas`): the array `slab.arrays['qekit_fijo']` and an ASE `FixAtoms` constraint. `inputgen.fixed_atoms` reads either and writes `0 0 0` in the third column of `ATOMIC_POSITIONS`.
3. Warnings: > 1.5 atoms per plane (cell is a multiple of the minimal one), real vacuum < 10 Å, polar slab, < 4 layers, freezing all the planes.
4. `report_slab` and, with `-o`, `structure.convert` writes CIF/POSCAR/XYZ. If there are fixed atoms and the format does not keep them (`structure.conserva_fijos`: only POSCAR/CONTCAR/`.vasp` store them as *Selective dynamics*), the CLI warns and recommends `builder.FORMATO_CON_FIJOS` (POSCAR or `.vasp`) or `olla-dft gamma --fix`.

**Where each datum comes from.**

| Datum | Source | Detail |
|---|---|---|
| Slab | `ase.build.surface` on the conventional cell | `--miller`, `--layers` (6), `--vacuum` (15 Å) |
| Conventional cell | spglib via `structure.conventional` | reference for the hkl indices |
| Atomic planes | distinct $z$ heights (tol 0.3 Å) | `builder._planos_z` |

**Limits and pitfalls.**
- *"la losa es POLAR… Añade 'dipfield = .true.' y 'edir = 3' al input, o corta una losa simétrica."*
- `--fix` is lost when exporting to CIF or XYZ: *"el CIF no tiene dónde ponerlo, así que al volver a cargarlo se relajaría todo. Escribe la losa en POSCAR (o .vasp)…"*. Only POSCAR keeps the `FixAtoms` constraint, which `inputgen.fixed_atoms` translates into `0 0 0`.
- Polarity detection only compares compositions of the extreme layers: a slab with terminations of the same composition but different geometry is not flagged.
- Cutting on the conventional cell may give a surface cell larger than the minimal one (a warning is issued).

**References.**
- P. W. Tasker, *J. Phys. C* 12, 4977 (1979) — polar surfaces.
- ASE: A. H. Larsen et al., *J. Phys.: Condens. Matter* 29, 273002 (2017).

---

### `olla-dft defect` — Building a point defect

**What it answers.** What do the perfect supercell and the supercell with a vacancy, a substitution or an interstitial look like, and what is the formation-energy formula that will have to be filled in?

**Background for non-experts.** A point defect is modelled by repeating the primitive cell $n_1\times n_2\times n_3$ times and modifying one atom. The supercell must be large so that the defect does not interact with its periodic images. This command only builds the two structures and writes the formula with its terms; `eform` does the calculation.

**Formulas.** `builder.formation_energy_text` writes:
$$E_f = E(\text{defect}) - E(\text{perfect}) \pm \mu(\cdot)\ \ [+\,q(E_F + E_v) + E_{\mathrm{corr}}]$$
- vacancy: $+\mu(\text{species that leaves})$; substitution: $+\mu(\text{leaves}) - \mu(\text{enters})$; interstitial: $-\mu(\text{enters})$.

**How Olla-DFT computes it.**
1. `structure.primitive` → `repeat(supercell)` (default 2×2×2).
2. `builder.defect`: vacancy (`del d[site]`), substitution (`d[site].symbol = new`), interstitial (fractional position `--position` of the supercell; warns if it ends up < 1.0 Å from a neighbour, minimum-image distance).
3. Warning if the shortest side of the supercell < 10 Å.
4. `report_defect` and writing of `perfecto.cif` and `defecto.cif` in `--outdir`.

**Where each datum comes from.**

| Datum | Source | Detail |
|---|---|---|
| Primitive cell | spglib via `structure.primitive` | basis of the supercell |
| Site, species, position | `--site`, `--new-element`, `--position` | 0-based indices in the supercell |

**Limits and pitfalls.** *"la supercelda mide X Å en su lado más corto: el defecto se ve con sus imágenes periódicas. Para energías de formación conviene ≥ 10-12 Å."* It relaxes nothing and computes no energies; the `--site` index refers to the repeated supercell, not to the input crystal.

**References.** C. Freysoldt, B. Grabowski, T. Hickel, J. Neugebauer, G. Kresse, A. Janotti, C. G. Van de Walle, *Rev. Mod. Phys.* 86, 253 (2014).

---

### `olla-dft eform` — Formation energy of charged defects

**What it answers.** How much does it cost to form the defect in each charge state, how does it vary with the Fermi level, where are the charge-transition levels and what is the finite-size correction?

**Background for non-experts.** Forming a defect costs an energy that depends on three things: where the atoms come from or go to (chemical potential $\mu$, fixed by the synthesis conditions), where the electrons come from or go to (Fermi level $\varepsilon_F$, measured from the valence-band maximum) and an artefact: a charged periodic cell interacts with its own images and with the neutralising background that QE adds. That artefact is corrected with the electrostatic energy of a point charge in a lattice of image charges (Makov–Payne) screened by the dielectric constant, or with the Lany–Zunger version that includes a shape term. The point where two lines $E_f(q)$ cross is a transition level: the Fermi level at which the defect changes charge.

**Formulas.** In `qekit/modules/defects.py`.

Formation energy (`DefectRun.E_f`):
$$E_f[D^q](\varepsilon_F) = E[D^q] - E[\mathrm{perf}] - \sum_i n_i\mu_i + q\,(\varepsilon_{\mathrm{VBM}} + \varepsilon_F) + E_{\mathrm{corr}}(q) + q\,\Delta V$$
- $n_i$: atoms **added** of species $i$ (−1 for the one that leaves); $\mu_i$: `--mu EL=eV` (for an elemental crystal, $\mu = E[\mathrm{perf}]/N$ automatically, `asignar_mu_elemental`); $\varepsilon_{\mathrm{VBM}}$ = `highestOccupiedLevel` of the perfect supercell (eV); $\varepsilon_F \in [0, E_g]$; $\Delta V$: potential alignment (`--dv` or `--align`).

Madelung constant by Ewald summation (`defects.madelung_xi`, `constante_madelung`):
$$\xi = \sum_{\mathbf R\neq 0}\frac{\mathrm{erfc}(\eta R)}{R} + \frac{4\pi}{V}\sum_{\mathbf G\neq 0}\frac{e^{-G^2/4\eta^2}}{G^2} - \frac{2\eta}{\sqrt\pi} - \frac{\pi}{\eta^2 V},\qquad \alpha_M = -\xi\,L,\quad L = V^{1/3}$$
- $\eta = \sqrt\pi / V^{1/3}$; real- and reciprocal-space cutoffs set by `tol` = 1e-10. Gives $\alpha_M = 2.8372974$ for the simple cubic lattice.

Image correction (`defects.correccion_imagen`):
$$E_{\mathrm{MP}} = \frac{k_e\,q^2\,\alpha_M}{2\,\varepsilon\,L},\qquad E_{\mathrm{LZ}} = E_{\mathrm{MP}}\left[1 + c_{\mathrm{sh}}\left(1 - \frac{1}{\varepsilon}\right)\right]$$
- $k_e$ = `KE` = 14.399645 eV·Å; $\varepsilon$ = `--epsilon`; $c_{\mathrm{sh}}$ = `C_SHAPE` = −0.35 (single value; LZ give −0.369 sc, −0.343 fcc, −0.342 bcc). `--correction` ∈ {`ninguna`, `makov-payne`, `lany-zunger`}.

Alignment (`defects.alineamiento`): $\Delta V = f\,\langle \bar V_{\mathrm{def}}(z) - \bar V_{\mathrm{perf}}(z)\rangle$ averaged over the 25 % of the cell opposite to the point of largest $|\Delta V - \mathrm{median}|$, with its standard deviation; $f$ = `UNIDADES_POTENCIAL[unidades_cube]` converts the `pp.x` cubes (`plot_num=11`, Ry, `unidades_cube="Ry"` by default, $f$ = `RY_EV`) to eV; with `"eV"`, $f = 1$. The result (`dV`, `sigma`, `perfil`) is always in eV.

Transition levels (`defects.niveles_transicion`), one entry per pair of consecutive charges $a<b$ (sorted by $q$), with the flag `dentro` = $0 \le \varepsilon \le E_g$; it is **not** filtered by the lower envelope (for the observable levels cross it with `envolvente`):
$$\varepsilon(a/b) = \frac{E_f(a, 0) - E_f(b, 0)}{b - a}$$

**How Olla-DFT computes it.**
1. `defects.prepare`: requires `--epsilon` if there are charges ≠ 0 and correction ≠ `ninguna`; builds the cells with `builder.defect`; resolves pseudos on the union of species.
2. Parity: if `--insulator` and some charge state leaves an odd number of electrons (`defects.electrones` with the `z_valence` of the UPFs), it activates `nspin=2` in **all** states with `tot_magnetization` 1 (odd) or 0 (even).
3. Writes `_perfecto/` (scf) and `qm1/`, `qp0/`, `qp1/`… (`relax` unless `--fixed-ions`, `tot_charge=q`, same estimated `nbnd`) and `run.sh`.
4. `--run`/`--collect`: `defects.collect` reads energies, convergence, `homo` (VBM) and `lumo` (gap) of the perfect cell; `--mu`; `--align POT_DEF POT_PERF` or `--dv`.
5. `report`: table $q$, $E$, $E_{\mathrm{corr}}$, $E_f(\varepsilon_F=0)$, $E_f(\varepsilon_F=E_g)$; transition levels, flagging those outside the gap; lower envelope (`envolvente`) and stable charges across the gap. `export`: `FORMACION.dat` (table and $E_f(\varepsilon_F)$ at 51 points); `plot`: $E_f$ vs $\varepsilon_F$.

**Where each datum comes from.**

| Datum | Source | Detail |
|---|---|---|
| $E[D^q]$, $E[\mathrm{perf}]$ | `pw.x` XML | `total_energy` |
| $\varepsilon_{\mathrm{VBM}}$, $E_g$ | XML of the perfect supercell | `highestOccupiedLevel`, `lowestUnoccupiedLevel` |
| $\mu_i$ | `--mu` or $E[\mathrm{perf}]/N$ (elemental) | eV/atom |
| $\varepsilon$ | `--epsilon` | e.g. $\varepsilon_1(0)$ from `optics` |
| $\alpha_M$ | Ewald sum over the real cell | `defects.madelung_xi` |
| $k_e$, $c_{\mathrm{sh}}$ | `defects.KE`, `defects.C_SHAPE` | 14.399645 eV·Å, −0.35 |
| $\Delta V$ | two potential `.cube` files (`--align`, Ry → eV) or `--dv` (already in eV) | `defects.alineamiento`, `UNIDADES_POTENCIAL` |
| Electrons per cell | `z_valence` of the UPFs | `defects.electrones` |

**Limits and pitfalls.**
- Without `--epsilon`: *"la constante dieléctrica es lo que apantalla la interacción del defecto con sus imágenes; sin ella la corrección sale ε veces de más."*
- With `--correction ninguna`: *"SIN CORREGIR: las E_f de los estados cargados están sistemáticamente bajas, y el error crece con q²."*
- `--dv` is given directly in eV (not converted); `--align` assumes `pp.x` cubes in Ry and the report says so: *"entra en E_f como q·ΔV = … eV por unidad de carga (el potencial de pp.x viene en Ry y se pasó a eV)"*. If $\sigma_{\Delta V} > 0.3\,|\Delta V|$: *"el defecto todavía se nota en la zona 'lejana', o sea que la supercelda es pequeña."*
- The listed transition levels include crossings between states that are never the most stable; the report marks *"<< fuera del gap"* those outside $[0, E_g]$, but a level inside the gap between two states that are not on the envelope is not observable either.
- Without a VBM (metal, no empty bands): *"No pude leer el VBM… E_f de los estados cargados no está definida."*
- Missing $\mu$ in a compound: *"FALTA el potencial químico… las DIFERENCIAS entre cargas y los niveles de transición sí valen, el valor absoluto de E_f no."*
- The correction only removes the leading $\propto q^2/L$ term; side < 10 Å with charge: warning.

**References.**
- G. Makov, M. C. Payne, *Phys. Rev. B* 51, 4014 (1995).
- S. Lany, A. Zunger, *Phys. Rev. B* 78, 235104 (2008); *Modelling Simul. Mater. Sci. Eng.* 17, 084002 (2009).
- C. Freysoldt, J. Neugebauer, C. G. Van de Walle, *Phys. Rev. Lett.* 102, 016402 (2009).
- C. Freysoldt et al., *Rev. Mod. Phys.* 86, 253 (2014). DOI: 10.1103/RevModPhys.86.253.

---

### `olla-dft interface` — Heterostructures and lattice mismatch

**What it answers.** Which common supercell allows stacking two 2D materials (or two slabs) with the least possible strain, how large is that strain and what does the initial structure look like?

**Background for non-experts.** Two crystal lattices almost never fit. To put them in the same periodic cell one must look for integer multiples of the vectors of each that resemble each other and stretch one of the two. That strain is the number that decides whether the calculation describes the material or a stretched version of it: 1 % is tolerable, 8 % is already another material.

**Formulas.** In `qekit/modules/interface.py`.

Candidate supercells (`_celdas_candidatas`): $\mathbf A' = M\mathbf a$, $\mathbf B' = N\mathbf b$ with $M, N \in \mathbb Z^{2\times2}$, $|M_{ij}|,|N_{ij}| \le$ `--max-index` (4), $\det > 0$, grouped by determinant (the areas must match within $2\cdot$`tol`).

Strain (`_deformacion`):
$$\boldsymbol\epsilon = B'^{-1}A' - I,\qquad \epsilon_{\max} = \max_{ij}|\epsilon_{ij}| \le \texttt{--tol}\ (0.05)$$

Lagrange–Gauss reduction (`reducir_2d`) so as not to repeat the same lattice with different bases; tie-break by "simplicity" of $M, N$ (`_simplicidad`: sum of |entries|, maximum, negatives, non-zeros).

Initial separation (`separacion_vdw`): $d_0 = 0.85\,(r_1 + r_2)$ with van der Waals radii from `R_VDW` (Bondi; 2.0 Å if missing).

With `--strain both`: target cell $= (w A' + v B')/(w+v)$ with $w = n_1\,|\det \mathbf a|$, $v = n_2\,|\det\mathbf b|$.

**How Olla-DFT computes it.**
1. `interface.buscar`: enumerates, filters by atoms (`--max-atoms` 200) and strain, deduplicates by $(n_1, n_2, \text{reduced shape}, \epsilon_{\max})$, sorts by $(\epsilon_{\max}, N_{\mathrm{at}}, \text{simplicity})$ and returns the `--top` (10) best. `--list` only prints them.
2. `interface.emparejar` chooses `--index` and `construir`: `ase.build.make_supercell` for each material, the in-plane cell is taken to the target dragging fractional positions (`_supercelda_deformada`), material 2 is stacked at `--separation` (or $d_0$) above material 1, `--shift` is applied (fractions of the common cell), `--vacuum` (20 Å) is added and the cell is centred.
3. Warnings: $\epsilon_{\max} > 3\,\%$, vdW separation as a starting point, registry not optimised.
4. `export`: `<name>.cif` and `<name>.txt`.

**Where each datum comes from.**

| Datum | Source | Detail |
|---|---|---|
| In-plane vectors | cells of `file1`, `file2` (rows 0–1, columns 0–1) | `interface._plano` |
| vdW radii | table `interface.R_VDW` | Å; `R_VDW_DEFECTO` = 2.0 |
| Search limits | `--max-index`, `--tol`, `--max-atoms` | 4, 0.05, 200 |

**Limits and pitfalls.**
- The **largest component** $\max|\epsilon_{ij}|$ of the matrix is reported, not a norm or an average: *"una deformación de 0 % en una dirección y 6 % en la otra no es '3 %'."*
- *"La deformación es del X %. Por encima de ~3 % no se está modelando el material sino una versión estirada de él."*
- The separation is a starting point: *"con un funcional sin corrección de dispersión la distancia de equilibrio saldrá demasiado grande."*
- *"El REGISTRO… no está optimizado. Dos apilamientos distintos pueden diferir en decenas de meV por átomo."*
- $c$ is assumed to be the normal and the cell a slab; the actual strain with `--strain both` does not coincide with the reported $\boldsymbol\epsilon$ (which is that of taking B to A).

**References.**
- A. Bondi, *J. Phys. Chem.* 68, 441 (1964) — van der Waals radii.
- P. Lazić, *Comput. Phys. Commun.* 197, 324 (2015) — CellMatch, lattice matching.

---

### `olla-dft neb` — Reaction barriers with neb.x

**What it answers.** What is the minimum-energy path between reactant and product and how high is the activation barrier (forward and backward)?

**Background for non-experts.** Between two energy minima there is a "mountain pass": the transition state. The nudged elastic band (NEB) stretches a chain of images between reactant and product, joined by springs, and relaxes each image perpendicular to the path until the chain rests in the valley. The climbing image (CI) pushes the highest image up to the exact pass; without it the barrier is underestimated.

**Formulas.** In `qekit/modules/neb.py`, `neb.collect`:
$$E_a^{\rightarrow} = E_{\max} - E_1,\qquad E_a^{\leftarrow} = E_{\max} - E_N,\qquad \Delta E = E_N - E_1$$
- Energies in eV relative to the first image (column 2 of `<prefix>.dat`); if `neb.out` contains `activation energy (->)`/`(<-)`, those are used. Conversion to kJ/mol: × 96.485.

**How Olla-DFT computes it.**
1. `neb.comprobar_extremos`: same number and **order** of atoms, same cell (tol 1e-4), non-identical structures; if it fails, it aborts.
2. `neb.build_neb_input` writes `neb.in`: `&PATH` with `string_method='neb'`, `nstep_path=--nstep` (50), `ds=1`, `opt_scheme='broyden'`, `num_of_images=--images` (7), `k_max=0.3`, `k_min=0.2`, `CI_scheme='auto'` (or `'no-CI'` with `--no-ci`), `path_thr=--path-thr` (0.05 eV/Å); `pw.x` engine trimmed from `inputgen.build_pw_input` (no positions or cell); `FIRST_IMAGE`/`LAST_IMAGE` in Å with `0 0 0` on the `--fix` atoms; `CELL_PARAMETERS`.
3. The user runs `neb.x -inp neb.in > neb.out`.
4. `neb.collect --collect`: reads `<prefix>.dat` (s, E, F), `<prefix>.int` (interpolation), and from `*.out` the barriers, convergence (`convergence achieved`), iterations, `CI_scheme` and the images with *"scf convergence NOT achieved on image"*.
5. `report`: barriers, table per image, warning if the interpolated maximum falls more than 0.4 steps from any image; `export` (`NEB.dat`, `NEB.txt`); `plot`.

**Where each datum comes from.**

| Datum | Source | Detail |
|---|---|---|
| $s$, $E$, $F$ per image | `<prefix>.dat` from `neb.x` | `neb.collect` |
| Interpolated curve | `<prefix>.int` from `neb.x` | optional |
| Barriers, convergence, iterations, CI | `neb.out` (regex) | take precedence over the own computation |
| eV → kJ/mol | 96.485 in `neb.report` | — |

**Limits and pitfalls.**
- *"Esta barrera es ELECTRÓNICA, a 0 K y sin energía de punto cero."* Thermal corrections in `thermochem`.
- Without CI: *"esta barrera es una COTA INFERIOR."* Few images (< 5): warning.
- Images with unconverged scf: *"El scf NO convergió en la(s) imagen(es)…: por eso el perfil sale dentado."*
- The endpoints must be relaxed with the same parameters; the module does not check this.

**References.**
- G. Henkelman, B. P. Uberuaga, H. Jónsson, *J. Chem. Phys.* 113, 9901 (2000) — climbing-image NEB. DOI: 10.1063/1.1329672.
- G. Henkelman, H. Jónsson, *J. Chem. Phys.* 113, 9978 (2000) — improved tangent.

---

### `olla-dft amorphous` — Amorphous solid by melt-quench with an MLIP

**What it answers.** How to generate an amorphous structure of given composition and density, and what coordination and first-neighbour distances does it have?

**Background for non-experts.** A glass is not drawn: it is manufactured by heating the material until it melts and cooling it so fast that it has no time to crystallise. On the computer the quench is millions of times faster than in the laboratory, so the result is more disordered and somewhat less dense than the real one. Here the dynamics is done with a machine-learned interatomic potential (MACE by default), not with DFT, because thousands of steps are needed; the resulting structure is a starting point that must then be relaxed with `pw.x`.

**Formulas.** In `qekit/modules/amorphous.py`.

Edge of the cubic cell (`celda_para_densidad`) and density (`densidad_de`):
$$L = \left(\frac{\sum_i m_i\,u}{\rho}\right)^{1/3}\times 10^{8},\qquad \rho = \frac{\sum_i m_i\,u}{V}$$
- $m_i$ in amu; $u = 1.66053906660\times10^{-24}$ g; $\rho$ in g/cm³; $V$ in Å³ (× $10^{-24}$ cm³).

Quench rate (`Protocolo.velocidad_temple`):
$$\dot T = \frac{T_{\mathrm{melt}} - T_{\mathrm{final}}}{N_{\mathrm{quench}}\,\Delta t}$$
- By default $(3000 - 300)\,\mathrm{K}/(1000\times 1\ \mathrm{fs}) = 2.7\times10^{15}$ K/s.

Coordination (`coordinaciones`): $Z_{ab} = \frac{1}{N_a}\sum_{i\in a}\#\{j\in b: d_{ij} < 1.25\,(r_a^{\mathrm{cov}} + r_b^{\mathrm{cov}})\}$ with minimum image; mean first-neighbour distance with the same cutoff (`distancia_media`).

**How Olla-DFT computes it.**
1. `formula_a_simbolos` expands `SiO2` × `--units` (8).
2. `empaquetar` places atoms at random (seed `--seed`) rejecting distances < `--min-dist` × (sum of covalent radii), `FACTOR_MINIMO` = 0.75; up to 20000 attempts per atom; error if they do not fit.
3. `fundir_y_templar` (unless `--pack-only`): calculator `mlip.calculator(--model)`; Maxwell–Boltzmann velocities at `--melt` (3000 K); ASE `Langevin` with `friction=0.02` and `--dt` (1 fs); `--melt-steps` (500) at $T_{\mathrm{melt}}$; quench in 20 segments of $N_{\mathrm{quench}}/20$ steps lowering the thermostat temperature linearly to `--final` (300 K); `--anneal-steps` (200) at $T_{\mathrm{final}}$. $E$ and $T$ are recorded every 10 steps (`traza.dat`).
4. Warnings: final temperature $> 2.5\,T_{\mathrm{final}} + 200$ K (the thermostat did not follow the ramp) and $\dot T > 10^{13}$ K/s.
5. `report` (density, protocol, coordinations, distances, final $T$) and `export` (`amorfo.cif`, `AMORFO.dat`, `AMORFO.txt`).

**Where each datum comes from.**

| Datum | Source | Detail |
|---|---|---|
| Masses and covalent radii | `ase.data.atomic_masses`, `covalent_radii` | — |
| amu → g | local constant $1.66053906660\times10^{-24}$ | CODATA 2018 |
| Energies and forces | MLIP potential (`mlip.calculator`) | MACE-MP-0 small, CHGNet or M3GNet |
| Target density | `--density` | g/cm³ |
| Protocol | `--melt`, `--final`, `--melt-steps`, `--quench-steps`, `--anneal-steps`, `--dt` | K, steps, fs |

**Limits and pitfalls.**
- *"Esta estructura viene de un potencial aprendido, NO de DFT… relájala con 'olla-dft gen -p relax'… y compara varias realizaciones (--seed distintas)."*
- The default protocol is an **exploration** one: 2.7×10¹⁵ K/s, and the report warns about it (*"Velocidad de temple X K/s. Un vidrio de verdad se enfría a 1-100 K/s"*). The docstring and the `--quench-steps` help say so: 27 000 steps bring it down to 10¹⁴ K/s, ten times more to 10¹³ K/s, where the warning disappears.
- NVT dynamics at fixed volume: the final density is the imposed one, it is not relaxed.
- With `friction=0.02` and fast ramps the system may end up liquid: *"El sistema acabó a X K, no a los Y K pedidos."*
- Requires `torch` + the model package (not dependencies of Olla-DFT).

**References.**
- I. Batatia et al., *MACE-MP-0* (arXiv:2401.00096, 2023).
- ASE Langevin: A. H. Larsen et al., *J. Phys.: Condens. Matter* 29, 273002 (2017).

---

### `olla-dft mlip` — Pre-relaxation, volume scan and phonon screening with a machine-learned potential

**What it answers.** Before spending DFT: what is a nearly relaxed geometry, where approximately is the $E(V)$ minimum, and does the structure have imaginary frequencies?

**Background for non-experts.** A machine-learned interatomic potential (MLIP) gives energies and forces thousands of times cheaper than DFT. It does not replace `pw.x` — it is trained on PBE data from Materials Project and describes *another* energy surface — but it helps reach the DFT calculation with the geometry almost ready, bound the range of an equation of state, and detect before DFPT that a structure is not at a minimum.

**Formulas.** In `qekit/modules/mlip.py`.

Relaxation (`mlip.relax`): ASE BFGS until $f_{\max} <$ `--fmax` (0.01 eV/Å) or `--steps` (300), with `FrechetCellFilter` if the cell is relaxed. Pressure:
$$P = -\tfrac{1}{3}\,\mathrm{tr}\,\boldsymbol\sigma\times 160.21766208\ \ [\mathrm{GPa}]$$

Volume scan (`mlip.volume_scan`): 15 scales in $[1-s, 1+s]$, $s$ = `--span` (0.10); parabola $E = aV^2 + bV + c$:
$$V_0 = -\frac{b}{2a},\qquad B_0 \approx 2aV_0\times160.21766208\ \mathrm{GPa},\qquad \text{scale} = (V_0/V)^{1/3}$$

Finite-difference phonons (`mlip.phonon_check`, `frequencies`): Hessian $H_{i\alpha,j\beta} = -\partial F_{j\beta}/\partial u_{i\alpha}$, central differences with $\delta$ = 0.01 Å in a `--supercell` (2×2×2), symmetrised; dynamical matrix $D = H/\sqrt{m_im_j}$;
$$\omega = \mathrm{sign}(\lambda)\sqrt{|\lambda|}\times 521.4708\ \mathrm{cm^{-1}}$$
- $\lambda$: eigenvalues of $D$ in eV/(Å²·amu); imaginary if $\omega < -5$ cm⁻¹.

**How Olla-DFT computes it.**
1. `mlip.calculator` loads MACE (`mace_mp(model=--size, default_dtype='float64')`), CHGNet or M3GNet; if the package is missing it explains what to install.
2. `relax`: initial/final forces and pressure, maximum displacement, volume change; warnings if it does not converge or if some atom moved > 0.5 Å. Writes the structure (`relajado_mlip.cif`) and `MLIP_PROCEDENCIA.json` (`write_provenance`) so that `audit` knows it is not DFT.
3. `scan`: `report_scan` suggests `olla-dft eos --scale X --span 0.04`; warns if the minimum falls outside the range.
4. `phonons`: `report_phonon`; exit code 1 if there are imaginary modes.

**Where each datum comes from.**

| Datum | Source | Detail |
|---|---|---|
| $E$, $F$, $\sigma$ | MLIP calculator | `mace_mp`, `CHGNetCalculator`, `PESCalculator` |
| eV/Å³ → GPa | 160.21766208 | local constant |
| $\sqrt{\mathrm{eV/(Å^2\,amu)}}$ → cm⁻¹ | 521.4708 | `CONV` |
| Masses | `atoms.get_masses()` (ASE) | amu |

**Limits and pitfalls.**
- *"ESTO NO ES EL RESULTADO FINAL. El modelo está entrenado con datos PBE… no mezcles sus energías con las de QE."* Example from the report: Si, MACE 5.464 Å vs LDA 5.402 Å.
- `phonon_check` diagonalises the **full** dynamical matrix of the supercell: the Γ modes of the primitive cell come out plus those of the q points the supercell folds onto Γ (the docstring states this). It is not a dispersion.
- The $B_0$ of the scan comes from a parabola: *"sirve para saber el orden de magnitud, no para reportarlo."*
- Without `torch`/`mace-torch`: *"para usar 'mace' hace falta instalar 'mace-torch'… Ocupa algo más de 1 GB."*

**References.**
- I. Batatia, D. P. Kovács, G. N. C. Simm, C. Ortner, G. Csányi, *NeurIPS* 35 (2022) — MACE.
- B. Deng et al., *Nat. Mach. Intell.* 5, 1031 (2023) — CHGNet.
- C. Chen, S. P. Ong, *Nat. Comput. Sci.* 2, 718 (2022) — M3GNet.

---

### `olla-dft audit` and `olla-dft db` — Comparability between calculations and local index

**What it answers.** Can the total energies of this set of calculations be subtracted? And `db`: which calculations do I have, with which parameters and what came out?

**Background for non-experts.** Two QE total energies can only be subtracted if they come from the same "recipe": same functional, same pseudopotentials, same cutoffs and same treatment of occupations. Otherwise the difference is a perfectly well-formed number without meaning, and QE does not warn. The audit computes a fingerprint with those parameters and groups: more than one group = not comparable. The k grid is treated separately as a warning, comparing the **density** of k points, which is what is comparable between different cells.

**Formulas.** `audit.kdensity`:
$$\rho_k = \frac{n_1 n_2 n_3}{(2\pi)^3 / V}\quad[\text{points}/\text{Å}^{-3}]$$

Fingerprint (`qeout.QEResult.fingerprint` + `origen`): (origin, functional, {element: UPF}, `ecutwfc`, `ecutrho`, `smearing`, `degauss`, `occupations`, `nspin`).

**Implemented rules.**

| Rule | Where | Effect |
|---|---|---|
| DFT vs MLIP origin enters the fingerprint | `audit.audit` (reads `MLIP_PROCEDENCIA.json` via `mlip.read_provenance`) | different groups: NOT COMPARABLE |
| Functional, pseudos, ecutwfc, ecutrho, smearing, degauss, occupations, nspin | `_campos`/`ETIQUETAS` | the differing ones are listed |
| Unconverged SCF | only `scf/relax/vc-relax/md/vc-md` with `converged=False` | "NO CONVERGIERON — sus energías no sirven" |
| `nscf`/`bands` | by calculation type | "Sin energía utilizable" |
| Disparate k density | $\max\rho_k/\min\rho_k > 2$ | WARNING, not an incompatibility |
| Folder without its own XML but with children | `audit.collect` | the children are audited (a sweep) |

**How Olla-DFT computes it.**
1. `audit.collect(paths)`: for each folder reads the MLIP mark and the XML (`qeout.read_xml`).
2. `audit.audit`: groups by fingerprint, lists differences, unconverged ones and those without energy.
3. `audit.report`; exit code 1 if not comparable. `--index` registers them in `olla-dft.db`.
4. `db folder/…` indexes (`audit.index`, `INSERT OR REPLACE` by absolute path); `db --query "SELECT …"` (SELECT only); `db --formula/--calculation/--gap-min/--gap-max` (`audit.search`); `db --export` (JSON); with no arguments, `audit.summary`.

**Where each datum comes from.** Everything from the `pw.x` XML (`qeout.read_xml`): functional, `pseudo_files`, cutoffs, smearing, occupations, `nspin`, energy (Ha → eV), volume, pressure, maximum force, `homo/lumo` → gap, magnetisation, convergence, SCF steps, `nk` (k points used), `nbnd`, BFGS steps, wall time; plus `MLIP_PROCEDENCIA.json` if present. Columns of the `calculos` table in `audit.ESQUEMA`.

**Limits and pitfalls.**
- The fingerprint includes neither the k grid nor the cell: *"un bulk y una losa necesitan mallas distintas por construcción."*
- It compares the **names** of the UPFs, not their content: two different files with the same name pass.
- `hull` and `thermo.from_runs` rely on this audit and refuse to mix origins.
- `db --query` only accepts `SELECT`; old databases are migrated by adding `nk`, `nbnd`, `n_bfgs` (`_migrar`).

**References.** Quantum ESPRESSO manual (`qes` XML schema); K. Lejaeghere et al., *Science* 351, aad3000 (2016) — why pseudos and cutoffs fix the energy reference.

---

### `olla-dft hull` — Formation energies and convex hull

**What it answers.** Is each phase stable against decomposing into the others, and how much energy per atom is it above the convex hull?

**Background for non-experts.** The formation energy per atom is plotted against composition. The lowest curve that envelops all points from below (convex hull) joins the stable phases; any phase above it gains energy by decomposing into the two (or three) hull phases surrounding it, and that vertical distance is $E_{\mathrm{hull}}$. It is energy at 0 K without entropy: a phase 25 meV/atom above is sometimes synthesised anyway.

**Formulas.** In `qekit/modules/thermo.py`.
$$E_f = \frac{E(\text{compound}) - \sum_i n_i\,\mu_i}{N},\qquad \mu_i = \min_{\text{pure phases of } i}\frac{E}{N}$$
$$E_{\mathrm{hull}} = E_f - E_{\mathrm{hull\ line}}(\mathbf x)$$
- Binary (`_casco`): lower envelope by monotone chain over $x$ and linear interpolation. Ternary or higher: `scipy.spatial.ConvexHull` in $(x_1,\dots,x_{n-1}, E_f)$, keeping facets whose normal points downward in energy (`eq[-2] < 0`), and $E_{\mathrm{hull\ line}}$ is obtained by barycentric coordinates inside the facet (`Delaunay.find_simplex`).

**How Olla-DFT computes it.**
1. `audit.collect` + `audit.audit`; if not comparable, it prints the audit and refuses unless `--force`.
2. `thermo.from_runs`: discards `nscf/bands`, no-energy or unconverged runs; refuses to mix DFT and MLIP; formula with `ase.Atoms`; elemental references = lowest energy per atom of the pure phases (warning if any is missing).
3. `_casco`; `report` with metastability threshold `--threshold` (0.025 eV/atom): ESTABLE / metaestable / inestable / fuera del dominio.
4. `export` (`CASCO_CONVEXO.dat`); `plot` only for binaries.

**Where each datum comes from.** Total energies and symbols from the XML of each folder (`qeout.read_xml`); element order from `--elements` or alphabetical.

**Limits and pitfalls.**
- *"Esto es energía a 0 K, sin punto cero ni entropía."*
- Without elemental references: *"hay que calcular cada elemento puro en su fase estable, con los mismos parámetros."*
- `--force` builds the hull with non-comparable calculations at the user's own risk.
- A pure element with several phases: the lowest is the reference; the others come out with $E_f > 0$.
- The plot is only for binaries.

**References.** S. P. Ong, L. Wang, B. Kang, G. Ceder, *Chem. Mater.* 20, 1798 (2008); W. Sun et al., *Sci. Adv.* 2, e1600225 (2016) — metastability scale.

---

### `olla-dft doctor` — pw.x convergence diagnostics

**What it answers.** Is this calculation usable and, if the SCF did not converge, is it because of charge sloshing (mix less) or slowness (mix more or more steps)?

**Background for non-experts.** The self-consistent cycle mixes the new density with the old one. If it mixes too much, the charge "sloshes" from one side of the cell to the other (oscillation, typical in slabs and metals) and the error goes up and down; if it mixes too little, the error always goes down but slowly. The two remedies are opposite, so the module looks at the **shape** of the `estimated scf accuracy` curve.

**Implemented rules** (`diagnose._clasificar`, only if not converged):

| Condition | Diagnosis | Advice |
|---|---|---|
| < 8 iterations | `pocos_datos` | raise `electron_maxstep` to ≥ 100 |
| (≥ 6 points and > 25 % of rises after the first 2 iterations) **or** the error grows > 5× in one iteration | `oscilacion` | `mixing_beta = max(0.05, β/3)`, `mixing_mode='local-TF'`, `mixing_ndim=12` |
| dropped < 3 orders of magnitude in total | `estancada` | check `starting_magnetization`, smearing, distances |
| otherwise, with β ≥ 0.6 | `lenta` | `electron_maxstep = 300` (do not raise β) |
| otherwise, with β < 0.6 | `lenta` | `mixing_beta = min(0.7, max(1.75β, 0.3))`, `electron_maxstep = 300` |

Problems from the XML (`diagnose.diagnose`): unconverged SCF; residual force > 0.05 eV/Å; $|P| > 1$ GPa in `scf/relax/vc-relax`; `Error in routine`. Relaxation: warning if the energy rose in more than $N/3$ steps.

**How Olla-DFT computes it.**
1. `qeout.find_xml` + `read_xml` (convergence, steps, error, forces, pressure, magnetisation, timings).
2. `diagnose.find_stdout` looks for the file containing `Program PWSCF`; `read_scf_history` splits the stdout into SCF cycles with `_ciclos_scf` (each `iteration #  1` opens one; in a `relax` there is one per ionic step), stores `n_ciclos` and extracts **from the last cycle only** `estimated scf accuracy`, `total energy`, `beta`, `convergence has been achieved` / `convergence NOT achieved`; `read_trajectory` reads the `!    total energy`, `Total force`, `P=` lines from the whole file.
3. `report` and `plot` (SCF accuracy on a log scale and energy per ionic step). Exit code 1 if there are problems. `--system` delegates to `health.check` (installation).

**Where each datum comes from.** XML (`converged`, `n_scf_steps`, `scf_error`, `max_force` in eV/Å, `pressure` in GPa, `wall_time`) and `pw.x` stdout (regex `_RE_ACC`, `_RE_ETOT`, `_RE_ITER`, `_RE_FORCE`, `_RE_PRESS`, `_RE_WARN`, `_RE_MAXSTEP`). $\beta$ defaults to 0.4 if not found.

**Limits and pitfalls.** In a `relax` only the last SCF cycle is diagnosed (the report says so: *"en el último de N ciclos SCF (uno por paso iónico; se diagnostica solo el último)"*); an intermediate cycle that oscillated is not seen. The thresholds (0.05 eV/Å, 1 GPa) are fixed. It does not detect symmetry or pseudopotential problems.

**References.** D. D. Johnson, *Phys. Rev. B* 38, 12807 (1988) — Broyden mixing; G. Kresse, J. Furthmüller, *Phys. Rev. B* 54, 11169 (1996) — charge sloshing and `local-TF`.

---

### `olla-dft crosscheck` — The same quantity by two independent routes

**What it answers.** Do two physically independent routes to the same quantity agree? If not, something is wrong in one of them.

**Background for non-experts.** Comparing against the literature detects errors in one module, but not a shared systematic bias. Computing $B_0$ from the equation of state and from the elastic constants, or the band gap and the Tauc gap, are routes that share no code: if they agree, it is hard for both to be wrong in the same way.

**Implemented checks** (`crosscheck.run`; relative deviation $|b-a|/|a|$, or absolute if $a = 0$):

| # | Quantity | Route A | Route B | Tolerance | Data |
|---|---|---|---|---|---|
| 1 | $B_0$ | `EOS.txt` (line with `B0` and `GPa`) | $B_{\mathrm{Hill}}$ from `ELASTIC_C.dat` (`elastic.moduli`) | 5 % | both files |
| 2 | $v_L[100]$, $v_T[100]$ | $\sqrt{C_{11}/\rho}$, $\sqrt{C_{44}/\rho}$ (`derived.cubic_directional`) | LA/TA slope at Γ from `FONONES_BANDAS.dat` | 10 % | Cij, bands, masses, volume |
| 3 | $\Theta_D$ | sound velocities (`derived.debye_from_velocity`) | second moment of `FONONES_DOS.dat` | 30 % (different definitions) | Cij, DOS, N |
| 4 | optical gap | `--gap-bandas` | `--gap-tauc` | 6 % | parameters |
| 5 | $C_v$ at 1500 K | $3Nk_B$ (Dulong–Petit) | $k_B\int x^2 e^x/(e^x-1)^2\,g(\omega)\,d\omega$ with $g$ normalised to $3N$ (`_cv_alta_T`) | 3 % | DOS |
| 6 | number of modes | $3N$ | $\int g(\omega)\,d\omega$ | 5 % | DOS |
| 7 | $\kappa_L$ | `KAPPA.dat` at ~300 K | Slack model from Cij (`derived.slack`) | 60 % | KAPPA, Cij |
| 8 | Berry phase | `BERRY.dat` (column 3 at charge 0) | $-2\sum_n (\bar r_n\cdot b)/2\pi$ from `WANNIER_centros.dat`, same branch mod 2 | 0.05 | both, cell |
| 9 | work function | `ESM.dat` (Φ at $q = 0$) | `WF.dat` (`Phi_eV`) | 5 % | both |
| 10 | $B_0$ (third route) | `EOS.txt` | $-\tfrac{1}{3}\,dP/d\epsilon$ from `STRAIN.dat` (kbar → GPa × 0.1) | 10 % | STRAIN (hydrostatic) |

Constants: `KB_EV` = $8.617333262\times10^{-5}$ eV/K; cm⁻¹ → eV: $1.239841984\times10^{-4}$.

**How Olla-DFT computes it.** `crosscheck._cargar` searches recursively for the result files in the project folder; with `-f structure` it takes masses, volume, N and cell; `run` executes every check for which data exist; `report` marks OK/FALLA with the diagnosis of what to look at first. Exit code 1 if any fails.

**Limits and pitfalls.** *"Un cruce que falla NO dice cuál de los dos caminos está mal."* Check 3 compares different definitions of $\Theta_D$ (*"coincidir al 1 % sería sospechoso"*); check 10 is only valid if the sweep was hydrostatic; check 2 blames the q grid before the Cij. Checks 8–10 swallow any exception silently (`except Exception: pass`).

**References.** R. Hill, *Proc. Phys. Soc. A* 65, 349 (1952); G. A. Slack, *Solid State Phys.* 34, 1 (1979); R. D. King-Smith, D. Vanderbilt, *Phys. Rev. B* 47, 1651 (1993).

---

### `olla-dft selftest` — Validation against known physics

**What it answers.** Does Olla-DFT reproduce measured, published or exact values, and not just what it says about itself?

**Background for non-experts.** Unit tests compare the code with itself. Here each test computes a quantity with a known answer (Ewald constants, Sackur–Tetrode entropy, $T_c$ of aluminium, topological invariants…) and checks it against that reference and its source. `--quick` (default) runs those that do not need `pw.x`; `--full` adds those that do; `--mlip` the one that needs MACE.

**Tests and references** (`selftest.PRUEBAS`; relative deviation, or absolute if the reference is 0):

| Key | Quantity | Reference | Tol. | Source (as stated in the code) | Function tested |
|---|---|---|---|---|---|
| `madelung` | $\alpha_M$ simple cubic | 2.8372974 | 1e-5 | classical Ewald value | `defects.constante_madelung` |
| `lorenz` | $L/L_0$ of the free-electron gas | 1.0 | 12 % | Sommerfeld limit | `transport.compute`, `lorenz` |
| `npw` | plane waves of Si at 30 Ry | 725 | 6 % | what `pw.x` reports (V = 39.5 Å³) | `cost.n_ondas_planas` |
| `sackur` | $S_{\mathrm{trans}}$ of N₂ at 298 K | 150.4 J/(mol·K) | 1 % | Sackur–Tetrode, NIST-JANAF | `thermochem.S_traslacional` |
| `allen_dynes` | $T_c$ of Al (λ=0.44, ω_log=270 K, µ*=0.12) | 1.18 K | 12 % | Allen–Dynes 1975, exp. | `elph.allen_dynes` |
| `allen_dynes_mu` | $T_c(0.10)/T_c(0.12)$ | 1.56 | 5 % | exponential dependence on µ* | `elph.allen_dynes` |
| `born2d` | $Y_{2D}$ with C11=352, C12=60 N/m | 341.8 N/m | 1 % | $Y = C_{11} - C_{12}^2/C_{11}$ (graphene DFT) | `elastic.modulos_2d` |
| `gap_invariante` | ΔE_v of a material with itself | 0 eV | 1e-9 | exact identity | `align.alinear` |
| `ewald_escala` | $\lvert\alpha(3) - \alpha(30)\rvert$ | 0 | 1e-6 | scale invariance | `defects.constante_madelung` |
| `chern_qwz` | $C$ of the Qi–Wu–Zhang model (m=−1) | −1 | 1e-10 | PRB 74, 085308 (2006) | `topology.invariants_from_vectors` |
| `umklapp` (`--mlip`) | exponent $n$ in $\kappa\propto T^{-n}$ of Si | 1.0 | 25 % | Umklapp law above $\Theta_D$ | `kappa.*` with MACE |
| `her_pt` | $\Delta G_{\mathrm{H^*}}$ with $E_{\mathrm{ads}} = -0.33$ eV | −0.09 eV | 5 % | Nørskov 2005, Pt(111) | `echem.her` |
| `oer_ruo2` | η with ΔG(OH,O,OOH)=(0.77, 2.16, 3.87) | 0.48 V | 10 % | Man et al. 2011 | `echem.oer` |
| `escala_oer` | ΔG(OOH) − ΔG(OH) of the RuO₂ profile | 3.2 eV | 10 % | universal scaling relation | `echem.oer` + `echem.escala_ooh_oh` |
| `escala_eta_min` | $\eta_{\min}$ = Δ/2 − ΔG_total/4 | 0.37 V | 2 % | Man et al. 2011 | `echem.sobrepotencial_minimo_escala` |
| `fonon_si` (`--full`) | optical ω(Γ) of Si | 520 cm⁻¹ | 10 % | exp. Raman 520.7 cm⁻¹ | `phonons.*` with `ph.x` |
| `wannier_si` (`--full`) | Si–Si Wannier centre | 1.17563 Å | 2 % | $\sqrt3\,a/8$ with a = 5.43 Å | `wannier.*` |
| `condensador` (`--full`) | slope of $1/C$ vs $d$ / $(1/\varepsilon_0)$ | 1.0 | 6 % | parallel-plate capacitor electrostatics | `esm.*` bc3 Al(111) |
| `born_si` (`--full`) | $Z^*$ of Si | 0 e | 0.05 | acoustic sum rule | `berry.*` |
| `gamma_al` (`--full`) | γ of Al(111) | 1.10 J/m² | 25 % | Vitos 1998 (1.20), exp. 1.14 | `surfen.*` |
| `bulk_si` (`--full`) | $B$ of Si by strain | 95 GPa | 15 % | LDA 93–97 (Nielsen & Martin 1985), exp. 98 | `strain.*` |
| `sitio_h_al` (`--full`) | $E_{\mathrm{ads}}$(top) − $E_{\mathrm{ads}}$(hollow), H/Al(111) | 5.6 eV | 60 % | hollow < bridge < top ordering | `adsorb.*` |

**How Olla-DFT computes it.** `selftest.ejecutar` filters by `--only`, `--full`, `--mlip`; creates a temporary folder (`--keep` to preserve it); runs each `fn(ctx)` and times it; `report` lists value, reference, deviation, tolerance and source. Exit code 1 if any fails or errors. `--list` prints the table without running anything.

**Limits and pitfalls.** *"Las que salen MAL no siempre son un fallo del código: una tolerancia ajustada, un pseudopotencial distinto o un cutoff bajo también las mueven."* The `--full` tests depend on the pseudos in `--pseudo-dir` and on `pw.x`/`ph.x` working.

**Note on `qekit/modules/uncertainty.py`.** It has no command of its own. It offers `propagate(f, values, sigmas)` — propagation in quadrature with central derivatives, $\sigma_f^2 = \sum_i (\partial f/\partial x_i)^2\sigma_i^2$, relative step $10^{-6}$, independent inputs — and `weighted_mean` — weighted mean with $w_i = 1/\sigma_i^2$ and $\sigma = (\sum w_i)^{-1/2}$. No module in this part calls it; only `validation`/`results` check that declared uncertainties are finite and non-negative.

**References.** P. B. Allen, R. C. Dynes, *Phys. Rev. B* 12, 905 (1975); X.-L. Qi, Y.-S. Wu, S.-C. Zhang, *Phys. Rev. B* 74, 085308 (2006); L. Vitos et al., *Surf. Sci.* 411, 186 (1998); O. H. Nielsen, R. M. Martin, *Phys. Rev. B* 32, 3792 (1985).

---

### `olla-dft suggest` — Parameters from your own history

**What it answers.** Based on the calculations that already converged with these elements, which `ecutwfc`, dual, k density and `electron_maxstep` should be used?

**Background for non-experts.** With a few dozen calculations there is no point in training anything: similar calculations are looked up (share elements, similar size) and what worked for them is examined, always stating how many cases back each number.

**Implemented rules** (`recommend.similares`, `recommend.sugerir`):

| Rule | Detail |
|---|---|
| Similarity | only calculations with `convergido`; score = Jaccard of elements $\lvert A\cap B\rvert/\lvert A\cup B\rvert$; × 0.5 if $N_{\mathrm{at}}$ differs by more than a factor of 2 |
| `ecutwfc` | **maximum** among the similar ones (not the mean), with range |
| dual | maximum of `ecutrho/ecutwfc` |
| k density | median of `kdensity` (points/Å⁻³) |
| `electron_maxstep = 300` | if the median of `n_scf` > 40 |
| `mixing_beta = 0.3` + `local-TF` | if the structure is a slab (vacuum along $c$ > 8 Å), general rule, 0 cases |
| Confidence | high ≥ 8 cases, medium ≥ 3, low < 3 |
| No history | refers to the cutoffs of the UPF itself / SSSP |

**How Olla-DFT computes it.** `_cmd_suggest` loads the structure, reads `SELECT * FROM calculos` from `--db` (`olla-dft.db`), detects whether it is a slab and calls `recommend.sugerir`; `report` prints value, number of cases and reason.

**Where each datum comes from.** `calculos` table of `olla-dft.db` (`audit.index`): `formula`, `natoms`, `ecutwfc`, `ecutrho`, `kdensity`, `n_scf`, `convergido`.

**Limits and pitfalls.** *"No sustituyen a una prueba de convergencia: 'olla-dft converge' sigue siendo la forma de saberlo de verdad."* With "low" confidence the report marks *"UN SOLO CASO: tómalo como indicio"* for 1 case and *"SOLO n CASOS"* for 2. It does not invent cutoffs without history.

**References.** G. Prandini, A. Marrazzo, I. E. Castelli, N. Mounet, N. Marzari, *npj Comput. Mater.* 4, 72 (2018) — SSSP.

---

### `olla-dft pseudos` — Choosing pseudopotentials with criteria

**What it answers.** Of the UPFs available for each element, which are usable for the task (optics, spin-orbit, XANES, DFT+U, phonons) and which one is advisable?

**Background for non-experts.** A folder usually holds several pseudopotentials per element, from different families and functionals. Picking the first one alphabetically fails silently: a scalar-relativistic pseudo with `lspinorb` gives a zero splitting, an ultrasoft one with `epsilon.x` gives a whole wrong spectrum, and mixing functionals between elements invalidates the total energy. The selector applies hard requirements (which discard) and preferences (which rank) and explains every decision.

**Implemented rules** (`pseudos.TAREAS`, `pseudos.evaluar`):

| Task | Hard requirement | Preference |
|---|---|---|
| `optics` | type ∈ {NC} | — |
| `soc` | relativistic = `full` (except elements with Z < 19: note and −0.5 points) | — |
| `xanes` | UPF with `PP_GIPAW` sections | — |
| `hubbard` | — | +0.15 × `z_valence` (semicore) |
| `fonones` | — | +2.0 if type ∈ {NC, US} |
| `general` | — | — |
| all | functional equal to `--functional` (aliases PBE/`SLA PW PBX PBC`, PZ/LDA, PBEsol, BLYP, revPBE) | +max(0, (90 − ecutwfc)/30); −0.5 without a declared cutoff; +1.0 US/PAW with `--cheap`; +0.3 with GIPAW; +0.2 if `full` |

Final order: non-discarded first, then descending points, then name. Coherence across elements (`pseudos.coherencia`): warning if functionals are mixed, if NC is mixed with US/PAW (the ultrasoft dual rules) and if the suggested cutoffs differ by more than 2.5×.

**How Olla-DFT computes it.**
1. `pseudos.candidatos`: `pseudo.find_for_element` (`.UPF` files whose name starts with the symbol) and `pseudos.leer` (type, functional normalised by `_funcional`/`NOMBRE_CORTO`, relativistic, `z_valence`, suggested cutoffs, GIPAW, size).
2. `pseudos.evaluar` and `elegir`; `report` with table and discarded ones; `report_coherencia` if there is more than one element; prints the `--pseudo EL=file` line for reuse.
3. The same selector is used by `sweep.prepare_common` in every command (`pseudo.resolve` → `_elegir` with the task) and `_coherencia_de_funcional` re-selects to unify the functional (preference PBE > PBEsol > revPBE > PZ > BLYP).

**Where each datum comes from.** UPF header (first 20–30 kB): `pseudo_type`, `functional`, `relativistic`, `z_valence`, `wfc_cutoff`/`rho_cutoff` (or their v1 equivalents), presence of `PP_GIPAW`. `Z_SOC` = 19 in `pseudos.py`.

**Limits and pitfalls.** *"Esto es una recomendación, no una verdad… hay que converger el cutoff con 'olla-dft converge'."* Type/functional are inferred by regex from the header: a UPF without those fields shows as `?` and is not discarded. The suggested cutoffs declared by the UPF are a starting point, not a convergence.

**References.** M. J. van Setten et al., *Comput. Phys. Commun.* 226, 39 (2018) — PseudoDojo; A. Dal Corso, *Comput. Mater. Sci.* 95, 337 (2014) — pslibrary; G. Prandini et al., *npj Comput. Mater.* 4, 72 (2018) — SSSP.
