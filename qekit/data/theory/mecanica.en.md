## Mechanics, vibrations, temperature and transport

This part documents the physics Olla-DFT implements in the commands that go from the total energy to mechanical, vibrational, thermal and transport properties: from convergence tests (`converge`, `tune`) and the equation of state (`eos`), through elastic constants (`elastic`, `derived`), strain sweeps (`strain`), surfaces and layered materials (`gamma`, `layers`, `xrd`, `exfoliate`), phonons and everything derived from them (`phonons`, `qha`, `thermochem`, `kappa`, `elph`), molecular dynamics (`md`), diffusive and ballistic electronic transport (`transport`, `ballistic`) and the cost estimator (`cost`). Every section was written by reading the code in `qekit/modules/*.py` and `qekit/cli.py`, and it only lists the formulas the code really executes, with the constants and defaults exactly as written. Whenever a docstring promises something the code does not do, it is said in "Limits and pitfalls". A note on file names: the module `qekit/modules/thermo.py` does NOT contain the harmonic thermodynamics (that lives in `phonons.thermodynamics`) but the convex hull of formation energies used by the `hull` command, which is documented elsewhere.

---

### `olla-dft converge` — Convergence of cutoffs and k-mesh

**What it answers.** From which `ecutwfc`, `ecutrho` or k-point mesh does the total energy stop changing by more than a threshold (1 meV/atom by default)? It is the first question for any new system.

**Background for non-experts.** A plane-wave calculation describes the electrons with a sum of waves; `ecutwfc` sets how many waves are included (the "resolution" of the wavefunctions), `ecutrho` the resolution of the charge density, and the k-mesh how many points of the Brillouin zone are sampled. With too few waves or points the result is coarse; with too many, the calculation costs more without gaining anything. A convergence test repeats the same calculation while raising the parameter and looks at when the energy "flattens", like adjusting the zoom of a microscope until the image stops changing.

The criterion Olla-DFT uses has an important subtlety: it compares each point against the DENSEST one in the series (the last), not against its previous neighbour. Two adjacent points may look alike by chance in the middle of a curve that has not flattened yet; comparing them with each other is the usual mistake.

**Formulas.** Per-atom difference with respect to the densest point (`converge.ConvergenceRun.per_atom_diffs`):

$$
\Delta E_i = \frac{|E_i - E_{\mathrm{ref}}|}{N_{\mathrm{at}}} \times 1000
$$

- $E_i$: total energy of point $i$, in eV per cell (read from the XML and converted from Hartree with $27.211386245988$ eV).
- $E_{\mathrm{ref}}$: energy of the last point that finished (the densest).
- $N_{\mathrm{at}}$: atoms in the cell.
- $\Delta E_i$: in meV/atom.

Convergence index (`converge.ConvergenceRun.converged_index`): the first $i$ such that every $\Delta E_j$ with $j \ge i$ satisfies $\Delta E_j \le$ threshold (failed points are ignored). k-mesh from a spacing (`kpoints.kgrid_from_spacing`):

$$
n_i = \left\lceil \frac{|\mathbf{b}_i|}{k_{\mathrm{spacing}}} \right\rceil, \qquad |\mathbf{b}_i| \text{ including the factor } 2\pi
$$

- $\mathbf{b}_i$: reciprocal vectors in Å⁻¹; $k_{\mathrm{spacing}}$ in Å⁻¹ (configuration `kspacing`, default 0.20). Directions with ≥ 8 Å of vacuum get a single point.

**How Olla-DFT computes it.**
1. `qekit/cli.py: _cmd_converge` loads the structure and calls `qekit/modules/converge.py: prepare`.
2. `sweep.prepare_common` resolves pseudopotentials and cutoffs (`pseudo.recommend_cutoffs`: the maximum declared by the UPF files; if none declares any, `ecutwfc` from the configuration (60 Ry) and `dual` (8); `ecutrho` never below $4\,\mathrm{ecutwfc}$).
3. Default series: `ecutwfc` = 30, 40, …, 100 Ry with `ecutrho = dual × ecutwfc`; `ecutrho` = 4, 6, 8, 10, 12 × ecutwfc; `kmesh` = meshes for the spacings 0.40, 0.30, 0.25, 0.20, 0.15, 0.12 Å⁻¹ (without repeats). `--values` replaces the series (for `kmesh` it accepts `8x8x8` or spacings).
4. One `pw.in` (`calculation='scf'`, `conv_thr = 1e-8`, `tstress`/`tprnfor` on) per point via `sweep.write_scf_job`, plus `run.sh` and `run.py`.
5. With `--run`, `runner.run_all` executes `pw.x`; with `--collect`, `converge.collect` reads `out/*.xml` (`qeout.read_xml`, tag `<total_energy><etot>`).
6. `converge.report` prints the table, the convergence point and the recommendation (`--ecutwfc N` or the mesh); `converge.export` writes `CONVERGENCIA.dat` and `.txt`; `converge.plot` draws $|\Delta E|$ on a log scale with the threshold band.

**Where each datum comes from.**

| Datum | Source | Detail |
|---|---|---|
| Total energy | pw.x XML (`output/total_energy/etot`, Hartree) | `qeout.read_xml` → eV |
| scf convergence | XML (`convergence_info/scf_conv/convergence_achieved`) | an unconverged point counts as failed under `--run` |
| Threshold | `--threshold` parameter | 1.0 meV/atom by default |
| Base cutoffs | UPF headers or `olla-dft config` | `pseudo.recommend_cutoffs` |
| Fixed k-mesh | `sweep.default_grid` | configuration `kspacing` (0.20 Å⁻¹) |
| Ry ↔ eV | `qeout.RY_EV` | 13.605693122994 eV |

**Limits and pitfalls.** It only looks at the total energy; the report warns: "convergence depends on the property: the total energy converges before stresses or phonons". If only the last point passes, it says: "Only the last point is below … there is no margin to be sure it has already flattened there". If none passes: "does NOT converge within … Extend the series towards denser values". The `energies` field of the dataclass is commented as "eV per cell", but the table is printed in Ry (divided by `RY_EV`): not a bug, just a display conversion. With `--collect` the inputs are not rewritten (`sweep.set_write_inputs(False)`), so the report describes what actually ran.

**References.** Quantum ESPRESSO manual (`pw.x`, variables `ecutwfc`, `ecutrho`, `K_POINTS`). Monkhorst and Pack, *Phys. Rev. B* 13, 5188 (1976), DOI 10.1103/PhysRevB.13.5188.

---

### `olla-dft tune` — Adaptive convergence recommendation

**What it answers.** Given an already generated `CONVERGENCIA.dat`, is the series converged, and if not, which value should be tried next?

**Background for non-experts.** It is pure post-processing of the `converge` table: it applies the same criterion ("from this point on, the whole tail stays within the threshold") and, when it is not met, proposes the next value with a sensible step instead of leaving the user to guess.

**Formulas.** Criterion (`tuning.analyze`): minimum index $i$ with $|\Delta E_j| \le$ threshold for all $j \ge i$. States: `ready` (such $i$ exists and is not the last), `confirm` (only the last passes), `extend` (none). Next value (`tuning._next_value`):

$$
v_{\mathrm{next}} = v_{\mathrm{last}} + \max\!\left(\mathrm{median}\{v_{k+1}-v_k > 0\},\; 0.10\,|v_{\mathrm{last}}|\right)
$$

- With fewer than two values, or no positive steps: $v_{\mathrm{next}} = 1.25\,v_{\mathrm{last}}$ (or $v_{\mathrm{last}}+1$ if it is not positive).

**How Olla-DFT computes it.**
1. `qekit/cli.py: _cmd_tune` → `qekit/modules/tuning.py: read` reads the numeric rows (column 1 value, 2 energy in Ry, 3 $\Delta E$ in meV/atom; comments and NaN are skipped).
2. `tuning.analyze` applies the criterion and picks the state and the recommended value.
3. `tuning.report` prints it; with `-o`, `tuning.export` writes a JSON (`CONVERGENCIA_RECOMENDACION.json` by default).

**Where each datum comes from.**

| Datum | Source | Detail |
|---|---|---|
| Value, E, ΔE | `CONVERGENCIA.dat` (from `olla-dft converge`) | `tuning.read`; ΔE taken in absolute value |
| Threshold | `--threshold` | 1.0 meV/atom if omitted; must be > 0 |

**Limits and pitfalls.** It runs nothing and reads no QE outputs: only the table. It uses the ΔE column as written, which `converge` computed against the densest point of THAT series; if points are added later, the table must be regenerated. The report reminds: "The energy property may converge before forces, phonons or tensors".

**References.** None specific; it is the same logic as `converge`.

---

### `olla-dft eos` — E–V equation of state and bulk modulus

**What it answers.** What are the equilibrium volume $V_0$, minimum energy $E_0$, bulk modulus $B_0$ and its derivative $B_0'$ of the crystal? And, if cubic, the lattice parameter $a_0$.

**Background for non-experts.** The cell is compressed and stretched a little around the starting size, the energy is computed at each volume, and a valley-shaped curve is fitted. The bottom of the valley is the equilibrium volume; the "stiffness" of the valley (its curvature) is the bulk modulus, which measures how hard you must squeeze to reduce the volume. Olla-DFT fits three different equations; if all three agree, the fit is reliable, and if they disagree, usually the range is too narrow or there are noisy points.

**Formulas.** Third-order Birch–Murnaghan (`eos.birch_murnaghan`), with $\eta = (V_0/V)^{2/3}$:

$$
E(V) = E_0 + \frac{9 V_0 B_0}{16}\left[(\eta-1)^3 B_0' + (\eta-1)^2 (6 - 4\eta)\right]
$$

Murnaghan (`eos.murnaghan`):

$$
E(V) = E_0 + \frac{B_0 V}{B_0'}\left[\frac{(V_0/V)^{B_0'}}{B_0'-1} + 1\right] - \frac{B_0 V_0}{B_0'-1}
$$

Vinet (`eos.vinet`), with $x = (V/V_0)^{1/3}$ and $\xi = \tfrac{3}{2}(B_0'-1)$:

$$
E(V) = E_0 + \frac{9 B_0 V_0}{\xi^2}\left[1 + \left(\xi(1-x) - 1\right) e^{\xi(1-x)}\right]
$$

- $V$, $V_0$: volumes in Å³; $E$, $E_0$: eV; $B_0$: eV/Å³ inside the fit, converted to GPa with `EV_A3_GPA = 160.21766208`; $B_0'$: dimensionless.
- Fit seed (`eos.fit`): parabola $E = aV^2+bV+c$ → $V_0 = -b/2a$, $B_0 = 2aV_0$, $B_0' = 4$.
- RMSE: $\sqrt{\langle (E - E_{\mathrm{fit}})^2 \rangle}/N_{\mathrm{at}}$, in eV/atom (printed in meV/atom).
- Cubic lattice parameter (`eos.fit`, field `EOSFit.a0`, with $V_{\mathrm{conv}}/V_{\mathrm{prim}}$ measured in `prepare`): $a_0 = (V_0 \cdot V_{\mathrm{conv}}/V_{\mathrm{prim}})^{1/3}$.

**How Olla-DFT computes it.**
1. `qekit/cli.py: _cmd_eos` → `qekit/modules/eos.py: prepare`. Requires `--npoints` ≥ 5 (default 9); `--span` 0.10 (±10 % in VOLUME); `--scale` (linear centring factor) 1.0.
2. Equally spaced volume factors $f \in [c^3(1-s),\, c^3(1+s)]$; linear factor $f^{1/3}$ applied to the cell with `set_cell(..., scale_atoms=True)`.
3. Cubic or not is decided by asking spglib (`structure.symmetry_dataset`, space group ≥ 195) and $V_{\mathrm{conv}}/V_{\mathrm{prim}}$ is stored via `structure.conventional`.
4. One `scf` (or `relax` with `--relax-ions`) per volume, all with the SAME k-mesh (`sweep.default_grid`), written by `sweep.write_scf_job` into `V_<factor>/pw.in`.
5. `--run` executes `pw.x`; `--collect`/`eos.collect` reads `etot` from the XML.
6. `eos.fit_all` fits the three equations with `scipy.optimize.curve_fit` (`maxfev=20000`); the fit is rejected if $V_0 \notin (0.6 V_{\min}, 1.4 V_{\max})$ or $B_0 \le 0$.
7. `eos.report` prints the table, the three fits, the Birch–Murnaghan result and the spread between equations; `eos.export` writes `EOS.dat` and `EOS.txt`; `eos.plot` draws $E - E_0$ with residuals.

**Where each datum comes from.**

| Datum | Source | Detail |
|---|---|---|
| Total energy per volume | pw.x XML (`etot`) | `qeout.read_xml` |
| Volume of each point | scaled cell (ASE) | $|\det(\mathbf{a})|$ in Å³ |
| Cubic symmetry and conventional cell | spglib via `structure` | `symmetry_dataset`, `conventional` |
| eV/Å³ → GPa | `eos.EV_A3_GPA` | 160.21766208 |
| Non-linear fit | `scipy.optimize.curve_fit` | library |

**Limits and pitfalls.** It does not relax the cell shape: only isotropic scaling (for non-cubic crystals $c/a$ stays fixed unless `--relax-ions` is used, which relaxes positions, not the cell). It warns: "V₀ falls OUTSIDE the computed range. Re-run the sweep centred on that volume" and, if the three equations differ by more than 5 %: "usually indicates missing points or a very narrow volume range". `--relax-ions` uses `calculation='relax'`, so those inputs carry `forc_conv_thr = 1e-4`.

**References.** F. Birch, *Phys. Rev.* 71, 809 (1947), DOI 10.1103/PhysRev.71.809. F. D. Murnaghan, *Proc. Natl. Acad. Sci. USA* 30, 244 (1944). P. Vinet, J. Ferrante, J. R. Smith and J. H. Rose, *J. Phys. C* 19, L467 (1986).

---

### `olla-dft elastic` — Elastic constants by stress–strain

**What it answers.** What are the elastic constants $C_{ij}$ of the crystal (or of the sheet, with `--2d`), the bulk, shear and Young moduli, the Poisson ratio, and is the structure mechanically stable?

**Background for non-experts.** A slightly deformed solid responds with a stress (force per unit area) proportional to the strain: this is the generalised Hooke's law, and the proportionality constants are the $C_{ij}$. Olla-DFT deforms the cell by ±1 % (and ±0.5 %) along each of the six independent directions (three stretches and three shears), asks `pw.x` for the stress tensor in each, and fits a straight line. Since `pw.x` gives all six stresses at once, every strain provides six equations, far fewer runs than the energy method.

In a sheet (graphene, MoS₂) stretching the vacuum direction makes no sense: only the two in-plane directions and the in-plane shear are strained, and the constants are given in N/m by multiplying by the cell height, so the vacuum cancels.

**Formulas.** Applied strain (`elastic.strain_matrix`): $\mathbf{a}' = \mathbf{a}(\mathbf{I}+\boldsymbol{\varepsilon})$ with $\varepsilon_{ii}=\delta$ for normal components and $\varepsilon_{ij}=\varepsilon_{ji}=\delta/2$ for shears (Voigt convention, $\varepsilon_4 = 2\varepsilon_{23}$). Fit (`elastic.fit`), with the sign inverted because the tensor `pw.x` writes is the opposite of the elasticity one:

$$
C_{ij} = -\frac{\partial\,\sigma^{\mathrm{pw}}_i}{\partial \varepsilon_j}\Big|_{\text{least squares}}, \qquad \sigma^{\mathrm{pw}}_i \to \sigma^{\mathrm{pw}}_i - \sigma^{\mathrm{pw}}_i(\text{ref})
$$

Voigt, Reuss and Hill averages (`elastic.moduli`), with $S = C^{-1}$:

$$
B_V = \frac{(C_{11}+C_{22}+C_{33}) + 2(C_{12}+C_{23}+C_{13})}{9}, \quad
G_V = \frac{(C_{11}+C_{22}+C_{33}) - (C_{12}+C_{23}+C_{13}) + 3(C_{44}+C_{55}+C_{66})}{15}
$$

$$
B_R = \frac{1}{(S_{11}+S_{22}+S_{33}) + 2(S_{12}+S_{23}+S_{13})}, \quad
G_R = \frac{15}{4(S_{11}+S_{22}+S_{33}) - 4(S_{12}+S_{23}+S_{13}) + 3(S_{44}+S_{55}+S_{66})}
$$

$$
B_H = \tfrac{1}{2}(B_V+B_R),\quad G_H = \tfrac{1}{2}(G_V+G_R),\quad
E = \frac{9 B_H G_H}{3B_H + G_H},\quad \nu = \frac{3B_H - 2G_H}{2(3B_H+G_H)},\quad
A^U = 5\frac{G_V}{G_R} + \frac{B_V}{B_R} - 6
$$

- $C_{ij}$, $B$, $G$, $E$ in GPa; $\nu$ and $A^U$ dimensionless; Pugh ratio $B_H/G_H$ (ductility threshold 1.75).
- Stability (generalised Born): all eigenvalues of $\tfrac{1}{2}(C+C^T)$ positive.

Sheet (`elastic.constantes_2d`, `modulos_2d`, `born_2d`): $C^{2D}_{ij} = C_{ij}\,c\times 0.1$ (GPa·Å → N/m), with $c$ the cell height;

$$
Y_x = \frac{C_{11}C_{22}-C_{12}^2}{C_{22}},\quad \nu_x = \frac{C_{12}}{C_{22}},\quad
K = \frac{C_{11}+C_{22}+2C_{12}}{4},\quad G = C_{66};\qquad
C_{11}>0,\; C_{66}>0,\; C_{11}C_{22}-C_{12}^2>0
$$

**How Olla-DFT computes it.**
1. `qekit/cli.py: _cmd_elastic` → `qekit/modules/elastic.py: prepare`. In 3D the structure is ALWAYS taken to spglib's standardised primitive cell (`structure.primitive`) so that the Cartesian axes coincide with the crystal-physical ones; in `--2d` it is not (it requires vacuum along $c$ via `kpoints.direcciones_con_vacio`).
2. Crystal family from the space-group number (`elastic.crystal_family`: ≥195 cubic, ≥168 hexagonal, ≥143 trigonal, ≥75 tetragonal, ≥16 orthorhombic).
3. Strains: `--delta` 0.010, `--npoints` 4 (even) → ±δ/2, ±δ. Components: the 6 Voigt ones, or (1, 2, 6) in 2D. Plus an undeformed reference cell.
4. `--ion-mode auto` (default): `scf` (fixed ions) for ε1–ε3 and `relax` for ε4–ε6; `relax`: all relaxed; `fixed`: all fixed. Deformed cells use `conv_thr = 1e-9`.
5. `pw.x` with `tstress = .true.`; `elastic.collect` reads `<stress>` from the XML (Ha/bohr³ → GPa with `qeout.HA_BOHR3_GPA = 29421.026`).
6. `elastic.fit` fits column by column with `np.polyfit(..., 1)`; `elastic.symmetrize` averages the equivalents of the family (cubic, hexagonal with $C_{66}=(C_{11}-C_{12})/2$, partial tetragonal); `elastic.moduli` computes VRH and Born.
7. `elastic.report`/`_report_2d`, `elastic.export` (`ELASTIC_C.dat`, `ELASTIC.txt`), `elastic.plot` (σ–ε lines).

**Where each datum comes from.**

| Datum | Source | Detail |
|---|---|---|
| Stress tensor | pw.x XML (`output/stress`, Ha/bohr³, Fortran order) | `qeout.read_xml`; requires `tstress=.true.` (always set) |
| Space group and family | spglib (`structure.symmetry_dataset`) | `elastic.crystal_family` |
| Cell height (2D) | `|a_3|` of the input cell | `ElasticRun.altura` |
| GPa·Å → N/m | `elastic.GPA_A_NM` | 0.1 |
| Assumed thickness (2D) | `--thickness` | only for the GPa equivalent, "a convention, not a measurement" |

**Limits and pitfalls.** Symmetrisation only covers cubic, hexagonal and (partially) tetragonal; trigonal, orthorhombic and monoclinic/triclinic are left with the symmetrised matrix $\tfrac{1}{2}(C+C^T)$ and nothing else. The Born criterion is the general one (eigenvalues), not the family-specific inequalities. It warns if the residual stress of the reference cell exceeds 0.5 GPa: "it is high. Relax the cell with vc-relax before computing the elastic constants". In 2D it warns that with `--ion-mode auto` the identity $C_{66}=(C_{11}-C_{12})/2$ stops holding even if the sheet is isotropic. The GPa equivalent of a sheet depends on the chosen thickness: "This thickness is a CONVENTION, not a measurement". With fewer than 3 stresses read, no fit is done.

**References.** R. Hill, *Proc. Phys. Soc. A* 65, 349 (1952), DOI 10.1088/0370-1298/65/5/307. S. I. Ranganathan and M. Ostoja-Starzewski, *Phys. Rev. Lett.* 101, 055504 (2008), DOI 10.1103/PhysRevLett.101.055504 ($A^U$ index). F. Mouhat and F.-X. Coudert, *Phys. Rev. B* 90, 224104 (2014), DOI 10.1103/PhysRevB.90.224104 (Born criteria). S. F. Pugh, *Philos. Mag.* 45, 823 (1954).

---

### `olla-dft strain` — Strain sweep: gap, energy and magnetic moment

**What it answers.** How do the energy, band gap, pressure and magnetic moment change when the cell is strained (biaxial, uniaxial, hydrostatic or shear)? What is the deformation potential $dE_{\mathrm{gap}}/d\varepsilon$ and at which strain does the gap close?

**Background for non-experts.** Stretching or compressing a crystal changes the interatomic distances and with them the electronic structure: the gap can open, close or change type, and a magnetic material can lose its moment. The deformation potential is the slope of that response, and it is what gets compared with pressure experiments or with sheets on substrates that stretch them. Olla-DFT ALWAYS applies each strain to the original cell (not to the previous point's cell, which would accumulate error) and relaxes the internal positions at each point.

**Formulas.** Strain (`strain.matriz`): $\mathbf{a}' = \mathbf{a}_0(\mathbf{I}+\boldsymbol{\varepsilon})$ with the Voigt components of each mode (`strain.MODOS`: biaxial (xx, yy), uniaxial-a/b/c, hydrostatic (xx, yy, zz), xy shear with $\varepsilon_{xy}=\varepsilon_{yx}=\varepsilon/2$). Energy minimum by a local parabola (`strain.minimo`, up to 3 points on each side of the sampled minimum): $\varepsilon^* = -b/2a$. Deformation potential (`strain.potencial_deformacion`):

$$
E_{\mathrm{gap}}(\varepsilon) \approx m\,\varepsilon + b, \qquad R^2 = 1 - \frac{\sum (y - \hat y)^2}{\sum (y-\bar y)^2}
$$

- $m$: in eV per unit strain (fraction, not per cent). Gap $= E_{\mathrm{LUMO}} - E_{\mathrm{HOMO}}$ from the XML.
- Gap closing (`strain.cierre_de_gap`): linear interpolation of the strain at which the gap crosses 0.02 eV.

2D biaxial modulus (`strain.modulo_biaxial`, biaxial mode only, points with $|\varepsilon|\le 0.03$):

$$
Y_{2D} = \frac{1}{A_0}\frac{d^2E}{d\varepsilon^2} \times 16.021766 \;\; [\mathrm{N/m}], \qquad \frac{d^2E}{d\varepsilon^2} = 2a
$$

- $A_0 = |\mathbf{a}_1\times\mathbf{a}_2|$ in Å²; $E$ in eV; 1 eV/Å² = 16.021766 N/m. It is the combination $C_{11}+2C_{12}+C_{22}$, NOT the Young modulus.

**How Olla-DFT computes it.**
1. `qekit/cli.py: _cmd_strain` → `qekit/modules/strain.py: prepare`. `--range MIN:MAX:N` in PER CENT (default `-5:5:11`; it rejects $|\varepsilon| > 30$ %, requires $N \ge 3$ and adds $\varepsilon=0$ if missing).
2. Estimates `nbnd` with empty bands (`inputgen._estimate_nbnd`: $\lceil 1.25\,N_{\mathrm{occ}} + 4\rceil$, ×1.2+2 if `nspin=2`) so that a LUMO exists.
3. Calculation type: `relax` (default), `scf` with `--fixed-ions`, `vc-relax` with `cell_dofree` (`z`, `shape` or `2Dxy` depending on the mode) with `--relax-perp`.
4. One input per strain (`sweep.write_scf_job`, accepts `--nspin/--mag`, `--hubbard`, `--vdw`).
5. `strain.collect` reads from each XML `etot`, `highestOccupiedLevel`, `lowestUnoccupiedLevel`, `stress` (pressure = trace/3), `magnetization/total` and `convergence_achieved`.
6. `strain.report` prints the table, the minimum, the deformation potential, the gap closing, the moment and the biaxial modulus (if there is vacuum along $c$); `strain.export` (`STRAIN.dat`, `.txt`); `strain.plot` (two panels).

**Where each datum comes from.**

| Datum | Source | Detail |
|---|---|---|
| Energy, HOMO, LUMO | pw.x XML (`etot`, `highestOccupiedLevel`, `lowestUnoccupiedLevel`) | `qeout.read_xml` |
| Pressure | XML (`stress`, trace/3, QE sign) | `QEResult.pressure` in GPa |
| Magnetic moment | XML (`magnetization/total`) | μ_B per cell |
| Convergence | XML (`convergence_achieved`) | rows flagged `<< SIN CONVERGER` |
| Reference area and volume | input cell | `StrainRun.area0`, `volume0` |
| eV/Å² → N/m | literal constant in `modulo_biaxial` | 16.021766 |

**Limits and pitfalls.** If the HOMO exists but the LUMO does not (no empty bands) the gap column stays empty and it warns: "No gap in the table: the calculations have no empty bands". If the minimum is not at $\varepsilon=0$ (|ε*| > 0.3 %): "The starting structure was not relaxed". With $R^2 < 0.9$: "The gap does not respond linearly in this range". Biaxial without vacuum: "if it is bulk material, perhaps you wanted 'hidrostatica'". `--relax-perp` with hydrostatic strain is rejected. Unconverged points DO enter the table (they are read from the XML even if the runner marks them failed) but it warns that they "are NOT comparable with the rest". The "gap" is that of the scf k-mesh, not the fundamental gap from a band path.

**References.** J. Bardeen and W. Shockley, *Phys. Rev.* 80, 72 (1950), DOI 10.1103/PhysRev.80.72 (deformation potentials). C. G. Van de Walle, *Phys. Rev. B* 39, 1871 (1989).

---

### `olla-dft gamma` — Surface energy and the Fiorentini–Methfessel fit

**What it answers.** How much energy per unit area does it cost to create the (hkl) surface of a crystal, $\gamma$ in J/m², and how does it converge with slab thickness?

**Background for non-experts.** Cutting a crystal leaves atoms with fewer neighbours: that costs energy, and the surface energy is that cost per unit area. It is computed with a "slab" (a few atomic layers with vacuum above and below) by subtracting what the same atoms would be worth inside the crystal. The problem is that the bulk energy comes from ANOTHER calculation, with another k-mesh, and any residual error per atom is multiplied by the number of atoms: $\gamma$ does not converge, it drifts. The way out is to fit a straight line $E_{\mathrm{slab}}(N)$ over several thicknesses: the intercept gives $2\gamma A$ and the slope a bulk energy consistent with the slabs themselves.

**Formulas.** Direct (`surfen.GammaRun.gamma_directo`):

$$
\gamma_{\mathrm{dir}}(N) = \frac{E_{\mathrm{slab}}(N) - N_{\mathrm{at}}\,E_{\mathrm{bulk}}}{2A}
$$

Fit (`surfen.ajustar`), least squares over the pairs $(N_{\mathrm{at}}, E_{\mathrm{slab}})$:

$$
E_{\mathrm{slab}}(N_{\mathrm{at}}) = 2\gamma A + N_{\mathrm{at}}\,E_{\mathrm{bulk}}^{\mathrm{fit}}
$$

- $E_{\mathrm{slab}}$: eV per slab cell; $N_{\mathrm{at}}$: atoms in the slab; $E_{\mathrm{bulk}}$: eV/atom from the separate conventional-cell calculation; $A = |\mathbf{a}_1\times\mathbf{a}_2|$ in Å² (one face); the 2 stands for the two faces (`GammaRun.caras` is always 2).
- $\gamma$ in eV/Å² → J/m² with `EV_A2_A_J_M2 = 16.021766`. Cleavage energy $= 2\gamma$.

**How Olla-DFT computes it.**
1. `qekit/cli.py: _cmd_gamma` → `qekit/modules/surfen.py: prepare`. `--miller` (default `1 0 0`), `--layers` 3,4,5,6 (at least two, minimum 3 layers), `--vacuum` 20 Å.
2. `builder.surface` cuts each slab from the CONVENTIONAL cell with `ase.build.surface`, centres the vacuum, detects whether it is symmetric ($z$ profile equal to its mirror, tol 0.3 Å) and polar (composition of the top layer ≠ bottom) and emits warnings.
3. Unless `--no-reduce` or `--fix`, `surfen.reducir_losa` replaces the slab by its spglib primitive if the $c$ axis does not change (same $\gamma$, fewer atoms).
4. k-mesh fixed with the smallest slab and reused for all (`sweep.default_grid`); the bulk (conventional cell) gets its own mesh. `scf` or `relax` (`--relax`) calculations, with options `--vdw`, `--dipole` (`dipole_correction=3`), `--nspin/--mag`.
5. `surfen.collect` reads `etot` and `convergence_achieved` from each XML and fits (`surfen.ajustar`).
6. `surfen.report` prints the direct-γ table, the drift, the fit with $R^2$ and the difference $E_{\mathrm{bulk}}^{\mathrm{fit}} - E_{\mathrm{bulk}}$; `surfen.export` (`GAMMA.dat`, `GAMMA.txt`); `surfen.plot`.

**Where each datum comes from.**

| Datum | Source | Detail |
|---|---|---|
| $E_{\mathrm{slab}}(N)$ | XML of each `capasNN/` (`etot`) | `qeout.read_xml` |
| $E_{\mathrm{bulk}}$ | XML of `_bulto/` (`etot`) / atoms of the conventional cell | skipped with `--no-bulk` |
| Area $A$ | cell of the thinnest slab | `np.cross(a1, a2)` |
| Symmetry / polarity | `builder.surface` (geometry) | tolerance 0.3 Å |
| eV/Å² → J/m² | `surfen.EV_A2_A_J_M2` | 16.021766 |

**Limits and pitfalls.** It is the LINEAR Fiorentini–Methfessel fit; Boettger's incremental scheme (which takes $E_{\mathrm{bulk}}$ from the difference between consecutive slabs) is not implemented, and the docstring says so. Non-symmetric slab: "γ is the AVERAGE of its two faces, not that of one". Polar: "use --dipole". Without `--relax`: "Unrelaxed: γ comes out high. Surface relaxation lowers it by 5 to 20 %". If the drift between slabs exceeds 0.05 J/m²: "It does not converge … It is the residual error of E_bulk multiplied by the number of atoms, not physics. The good value is the fitted one". $R^2 < 0.999$: "either some point lacks convergence, or the thin slabs do not yet have a bulk interior". With `--fix` the cell is not reduced (constraints refer to specific atoms).

**References.** V. Fiorentini and M. Methfessel, *J. Phys.: Condens. Matter* 8, 6525 (1996), DOI 10.1088/0953-8984/8/36/005. J. C. Boettger, *Phys. Rev. B* 49, 16798 (1994), DOI 10.1103/PhysRevB.49.16798.

---

### `olla-dft layers` — Layer detection by connectivity

**What it answers.** Is the structure layered? How many layers per cell, along which axis are they stacked, what are the basal spacing $d$ and the interlayer gap, and where would the (00l) basal peak fall in a diffractogram?

**Background for non-experts.** The geometry is not eyeballed; bonds are: two atoms are bonded if their distance does not exceed the sum of covalent radii plus a tolerance. The bond network is built respecting periodicity and the connected pieces are separated. A piece that repeats in exactly two directions is a layer; in three, a 3D framework; in one, a chain; in none, a molecule. The dimensionality is read from the "closure vectors": walking the bonds while assigning each atom a cell displacement relative to a root atom, every bond that "does not fit" contributes an integer vector, and the rank of the set of those vectors is the number of periodic directions.

**Formulas.** Bond criterion (`layers.bonds`, with `ase.neighborlist.neighbor_list` and per-atom radii $r_i + \mathrm{tol}/2$):

$$
d_{ij} \le r^{\mathrm{cov}}_i + r^{\mathrm{cov}}_j + \mathrm{tol}, \qquad \mathrm{tol} = 0.45\ \text{Å (default)}
$$

Dimensionality (`layers._components_and_dim`): $\dim = \operatorname{rank}\{\mathbf{d}\}$ with $\mathbf{d} = \mathbf{o}_a + \mathbf{S}_{ab} - \mathbf{o}_b \ne 0$. Spacings (`layers.analyze`):

$$
d_{\mathrm{basal}} = \frac{P}{n_{\mathrm{layers}}}, \qquad P = |\mathbf{a}_{\mathrm{stack}}\cdot\hat{\mathbf{n}}|, \qquad
\mathrm{gap} = \min_k\left(z^{\mathrm{bot}}_{k+1} - z^{\mathrm{top}}_k\right)
$$

Basal reflection (`layers.report`), Bragg: $2\theta = 2\arcsin\!\left(\lambda/(2 d_{\mathrm{basal}}/l)\right)$ for $l = 1, 2, 3$.

- $\hat{\mathbf{n}}$: unit normal to the plane of the two non-stacking cell vectors; $z^{\mathrm{top/bot}}_k$: centre ± thickness/2 of layer $k$ (no van der Waals radii); $\lambda$ in Å.

**How Olla-DFT computes it.**
1. `qekit/cli.py: _cmd_layers` → `qekit/core/layers.py: analyze` (`--tol` 0.45 Å).
2. Bonds with ASE; connected components and rank with `np.linalg.matrix_rank` over the closure vectors.
3. Stacking axis: the fractional direction outside the plane spanned by the closure vectors (SVD), taking the cell vector with the largest out-of-plane component.
4. Each layer is rebuilt contiguous (BFS with Cartesian displacements) to measure centre and thickness without cell cuts.
5. `layers.report` prints layers, $d$, gap, period and the 00l reflections with the λ from `--wavelength` (default CuKα = 1.54184 Å, `xrd.wavelength_value`), labelled with the actual radiation name (`xrd.wavelength_name`: "Cu Kα", "Mo Kα1", or "λ dada" if a number was given).
6. With `--slab FILE`, `layers.make_slab` isolates the first layer: it unwraps it along the stacking axis (minimum image in fractional coordinates relative to the first atom, so a layer crossing the cell boundary is not split), replaces the stacking vector by the normal with height thickness + `--vacuum` (20 Å), centres it and `structure.convert` writes it.

**Where each datum comes from.**

| Datum | Source | Detail |
|---|---|---|
| Covalent radii | `ase.data.covalent_radii` | library |
| Periodic neighbour list | `ase.neighborlist.neighbor_list("ijS")` | library |
| Wavelength | `xrd.WAVELENGTHS` or a value in Å | CuKa 1.54184, MoKa 0.71073, CoKa 1.79026, … |
| Tolerance | `--tol` | 0.45 Å (`layers.DEFAULT_TOL`) |

**Limits and pitfalls.** No QE calculation: it is pure geometry. If there are no 2D components: "No layers detected … you can try a smaller --tol". The stacking axis and the normal are computed from the FIRST layer; layers with different orientations will not be noticed. `make_slab` only unwraps along the stacking axis (not in-plane), which is all that affects centring and thickness.

**References.** M. Ashton, J. Paul, S. B. Sinnott and R. G. Hennig, *Phys. Rev. Lett.* 118, 106101 (2017), DOI 10.1103/PhysRevLett.118.106101 (topological dimensionality criterion). W. H. Bragg and W. L. Bragg, *Proc. R. Soc. Lond. A* 88, 428 (1913).

---

### `olla-dft xrd` — Simulated powder diffractogram

**What it answers.** Where do the powder X-ray diffraction peaks of this structure appear, with which relative intensity and which hkl indices? Does it look like the measured diffractogram?

**Background for non-experts.** A crystal diffracts X-rays in directions fixed by Bragg's law: each family of planes with spacing $d$ produces a peak at an angle $2\theta$. The intensity depends on how the waves scattered by each atom in the cell interfere (structure factor), on how much each element scatters (atomic scattering factor, which decays with angle), and on geometric factors of the experiment (Lorentz–polarisation). Small crystallites broaden the peaks (Scherrer). Olla-DFT computes all of that and overlays, if given, an experimental diffractogram, to see at a glance whether the structural model is the right one.

**Formulas.** (`xrd.compute`) For each $hkl$ with $|\mathbf{g}| = 1/d$ inside the accessible sphere:

$$
\sin\theta = \frac{\lambda |\mathbf{g}|}{2}, \qquad s^2 = \left(\frac{\sin\theta}{\lambda}\right)^2 = \frac{|\mathbf{g}|^2}{4}
$$

$$
f(s) = Z - 41.78214\, s^2 \sum_{i=1}^{4} a_i\, e^{-b_i s^2}, \qquad f \to f\, e^{-B_{\mathrm{iso}} s^2}
$$

$$
F(hkl) = \sum_j f_j(s)\, e^{2\pi i\,(hkl)\cdot\mathbf{r}_j}, \qquad
I \propto |F|^2 \cdot \frac{1+\cos^2 2\theta}{\sin^2\theta\,\cos\theta}
$$

- $\lambda$: Å; $\mathbf{g} = (hkl)\,\mathbf{B}$ with $\mathbf{B} = (\mathbf{A}^{-1})^T$ without $2\pi$; $Z$: atomic number; $a_i, b_i$: analytical coefficients (data file taken from pymatgen, values from the *International Tables*); $B_{\mathrm{iso}}$: global temperature factor in Å² (`--biso`, 0 by default); $\mathbf{r}_j$: fractional positions.
- Multiplicity: emerges by enumerating ALL hkl and merging those coinciding in $2\theta$ (tolerance 0.02°). Intensities normalised to 100; peaks < 0.1 are dropped.

Profile (`xrd.broaden`), pseudo-Voigt with $\eta = 0.5$ and width $w$ (FWHM in ° 2θ):

$$
y(x) = \sum_p I_p\left[(1-\eta)\, e^{-\frac{(x-x_p)^2}{2\sigma^2}} + \eta\,\frac{1}{1+\left(\frac{x-x_p}{w/2}\right)^2}\right], \quad \sigma = \frac{w}{2\sqrt{2\ln 2}}, \quad
w_{\mathrm{Scherrer}} = \frac{K\lambda}{L\cos\theta},\; K = 0.9
$$

- $L$: crystallite size (`--size` in nm, converted to Å); without `--size`, $w$ = `--fwhm` (0.15°).

**How Olla-DFT computes it.**
1. `qekit/cli.py: _cmd_xrd` → `qekit/modules/xrd.py: compute`. With `--basis conventional` (default) the cell is standardised to the conventional one (`structure.conventional`) so the hkl match the PDF cards; `input` uses the cell as given.
2. Enumerates hkl in the box $|h_i| \le \lceil g_{\max}/|\mathbf{b}_i|\rceil + 1$, filters $g_{\min} \le |\mathbf{g}| \le g_{\max}$ (range `--tt-min` 5°, `--tt-max` 70°).
3. Factors $f_j(s)$ from `qekit/data/atomic_scattering_params.json` (`xrd.scattering_params`), vectorised structure factor, LP, removal of extinct reflections ($I < 10^{-8} I_{\max}$), merging by 2θ and the "most readable" hkl label (`Peak.label`, Friedel-oriented).
4. `xrd.broaden` generates the continuous profile (step 0.02°).
5. `xrd.read_experimental` reads `--exp` (two columns, ≥ 10 rows; subtracts the minimum and normalises to 100).
6. `xrd.report` (12 strongest peaks), `xrd.export` (`XRD.dat`, `XRD_HKL.dat`), `xrd.plot` (experimental offset +105 above the simulated), `--suite` (exchange JSON).

**Where each datum comes from.**

| Datum | Source | Detail |
|---|---|---|
| Coefficients $a_i, b_i$ | `qekit/data/atomic_scattering_params.json` | from pymatgen (MIT), *International Tables* |
| $Z$ | `ase.data.atomic_numbers` | library |
| Wavelength | `xrd.WAVELENGTHS` | CuKa 1.54184 Å, CuKa1 1.54056, MoKa 0.71073, CoKa 1.79026, FeKa 1.93735, CrKa 2.29100, AgKa 0.56087 |
| Scherrer constant | `xrd.SCHERRER_K` | 0.9 |
| Conventional cell | spglib (`structure.conventional`) | if it fails, the input cell is used |

**Limits and pitfalls.** There is no Rietveld refinement nor R factor: the comparison with `--exp` is purely visual (overlay). No absorption factor, preferred orientation, anomalous-dispersion correction or Kα1/Kα2 doublet (a single λ). The temperature factor is a single isotropic $B$ for all atoms. Scherrer broadening ignores strain broadening. The $f(s)$ formula is pymatgen's parametrisation (the same as the *International Tables* in the form $Z - 41.78214 s^2\sum a_i e^{-b_i s^2}$), valid for X-rays. With `--basis input` on a primitive cell, "the hkl are NOT those of the PDF card".

**References.** *International Tables for Crystallography*, Vol. C (scattering factors). P. Scherrer, *Nachr. Ges. Wiss. Göttingen* 2, 98 (1918). S. P. Ong et al., *Comput. Mater. Sci.* 68, 314 (2013), DOI 10.1016/j.commatsci.2012.10.028 (pymatgen, origin of the coefficients). B. E. Warren, *X-ray Diffraction*, Dover (1990).

---

### `olla-dft exfoliate` — Exfoliation energy

**What it answers.** How much does it cost to separate one layer from the layered crystal, in J/m² (and meV/Å², meV/atom)? Is it exfoliable?

**Background for non-experts.** The energy of the crystal (per layer) is compared with that of an isolated monolayer in vacuum. The difference per unit area is the exfoliation energy; typical layered materials lie between 0.2 and 0.6 J/m² (graphite ≈ 0.35 J/m² experimental). Interlayer cohesion is mostly van der Waals dispersion, which LDA and PBE describe poorly: without a dispersion correction the number is not comparable with experiment, and the module says so.

**Formulas.** (`exfoliate.report_result`)

$$
E_{\mathrm{exf}} = \frac{E_{\mathrm{mono}} - E_{\mathrm{bulk}}/N_{\mathrm{layers}}}{A}
$$

- $E_{\mathrm{mono}}$, $E_{\mathrm{bulk}}$: total energies in eV; $N_{\mathrm{layers}}$: layers per bulk cell detected by `layers.analyze`; $A = |\mathbf{a}_i\times\mathbf{a}_j|$ of the two non-stacking vectors, in Å². Conversion with `EV_A2_TO_J_M2 = 16.02176634`.

**How Olla-DFT computes it.**
1. `qekit/cli.py: _cmd_exfoliate` → `qekit/modules/exfoliate.py: prepare`: `layers.analyze(atoms, tol)` (`--tol` 0.45 Å); without layers, usage error.
2. `layers.make_slab` builds the monolayer (first layer) with `--vacuum` 20 Å.
3. Bulk k-mesh from `sweep.default_grid`; the monolayer's is the same in-plane and 1 along the stacking axis.
4. Two `scf` runs (`bulk/pw.in`, `monocapa/pw.in`; `relax` for the monolayer with `--relax-slab`), both with the same `--vdw` (`grimme-d2`, `grimme-d3`, `DFT-D`, `ts-vdw`, `xdm`, `mbd`).
5. `exfoliate.collect` reads `etot` from both XML files; `exfoliate.report_result` prints the result.

**Where each datum comes from.**

| Datum | Source | Detail |
|---|---|---|
| $E_{\mathrm{bulk}}$, $E_{\mathrm{mono}}$ | pw.x XML (`etot`) | `qeout.read_xml` |
| $N_{\mathrm{layers}}$, stacking axis | `layers.analyze` | bond connectivity |
| Area $A$ | bulk cell | two non-stacking vectors |
| eV/Å² → J/m² | `exfoliate.EV_A2_TO_J_M2` | 16.02176634 |

**Limits and pitfalls.** It assumes all layers in the cell are equivalent (divides $E_{\mathrm{bulk}}$ by $N_{\mathrm{layers}}$ and uses only the first one). Without `--vdw`: "WITHOUT van der Waals correction: PBE barely binds the layers and LDA binds by error cancellation". With pseudopotentials that look like LDA (name containing `pz`, `lda`, `pw92`) plus Grimme: "combining them counts dispersion twice". If negative: "Almost always means the vdW correction is missing or some calculation is not well converged". There is no bulk relaxation. It writes no `.dat` files (only the on-screen report).

**References.** J. H. Jung, C.-H. Park and J. Ihm, *Nano Lett.* 18, 2759 (2018), DOI 10.1021/acs.nanolett.7b04201 (exfoliation vs. interlayer binding energy). S. Grimme, *J. Comput. Chem.* 27, 1787 (2006), DOI 10.1002/jcc.20495 (D2); S. Grimme et al., *J. Chem. Phys.* 132, 154104 (2010), DOI 10.1063/1.3382344 (D3).

---

### `olla-dft phonons` — DFPT phonons: dispersion, DOS, thermodynamics, IR, Raman and electronic temperature

**What it answers.** What are the vibrational frequencies of the crystal (at Γ or across the Brillouin zone), is it dynamically stable (no imaginary frequencies), what are its harmonic zero-point energy, free energy, entropy and heat capacity, which modes are IR and Raman active (`--raman`) and, with `--tscan`, does a soft mode stabilise as the electronic temperature rises?

**Background for non-experts.** The atoms of a crystal vibrate around their equilibrium positions as if joined by springs. Density-functional perturbation theory (DFPT, what `ph.x` does) computes the stiffness of those springs (the force constants) from the electronic density, without displacing atoms by hand. From there come the frequencies of all vibrational waves (phonons). An "imaginary" frequency (negative in the output) means the structure is not at a minimum: either it was not relaxed well or it is unstable. With the phonon density of states the harmonic thermodynamics is computed: even at 0 K the atoms vibrate (zero-point energy) and on heating more modes get populated (entropy, heat capacity). The Γ modes are the ones an infrared or Raman spectrometer sees.

**Formulas.** Harmonic thermodynamics per cell from the DOS $g(\omega)$ (`phonons.thermodynamics`), with $\epsilon = \hbar\omega$ in eV, $x = \epsilon/k_BT$ (capped at 500), $n = 1/(e^x - 1)$, and $g$ renormalised to $\int g\,d\omega = 3N_{\mathrm{at}}$ (only $\omega > 1$ cm⁻¹):

$$
E_{\mathrm{ZPE}} = \int \tfrac{1}{2}\epsilon\, g\, d\omega, \qquad
F(T) = E_{\mathrm{ZPE}} + k_B T \int \ln\!\left(1 - e^{-x}\right) g\, d\omega
$$

$$
U(T) = \int \left(\tfrac{1}{2} + n\right)\epsilon\, g\, d\omega, \qquad
C_v(T) = k_B \int x^2 e^{x} n^2\, g\, d\omega, \qquad S(T) = \frac{U - F}{T}
$$

- $k_B$ = `KB_EV` = 8.617333262e-5 eV/K; cm⁻¹ → eV with `CM1_TO_EV` = 1.239841984e-4; cm⁻¹ → THz with `CM1_TO_THZ` = 0.0299792458. Trapezoidal integrals. $T$ = 0…1000 K in steps of 10.

Stokes Raman spectrum (`phonons.raman_spectrum`) from the activity $A$ (Å⁴/amu) of `dynmat.x`:

$$
I(\omega) \propto \frac{(\omega_L - \omega)^4}{\omega}\,[n(\omega,T)+1]\,A(\omega), \qquad \omega_L = \frac{10^7}{\lambda_{\mathrm{laser}}[\mathrm{nm}]}\ \mathrm{cm^{-1}}
$$

convolved with Lorentzians of FWHM 5 cm⁻¹ at $T$ = 300 K; $\omega \le 1$ cm⁻¹ is excluded. Electronic temperature (`tphonons.degauss_de_T`):

$$
\mathrm{degauss} = k_B T, \qquad k_B = 6.333621\times10^{-6}\ \mathrm{Ry/K}, \qquad \text{smearing = fermi-dirac}
$$

Stabilisation temperature (`tphonons.temperatura_de_estabilizacion`): linear interpolation of the $T$ at which the softest mode (minimum of the frequencies with $|\omega| > 10$ cm⁻¹) crosses from negative to non-negative.

**How Olla-DFT computes it.**
1. `qekit/cli.py: _cmd_phonons` → `qekit/modules/phonons.py: prepare`. The structure is taken to the standardised primitive cell (`structure.primitive`). It writes `scf.in` (`conv_thr = 1e-12`, k-mesh from `kspacing`), `ph.in` (`tr2_ph = 1e-14`, `fildyn='dyn'`; `epsil=.true.` if `--insulator` or `--raman`; `lraman=.true.` and `trans=.true.` with `--raman`; `ldisp` with `nq` mesh = `--qgrid` or `kgrid_from_spacing(atoms, 0.6)`).
2. Mesh mode: `q2r.in` (`zasr='simple'`, `flfrc='fuerzas.fc'`), `matdyn_band.in` (seekpath path via `kpoints.get_kpath`, 30 points per segment, `q_in_band_form`, `q_in_cryst_coord`), `matdyn_dos.in` (`dos=.true.`, 12×12×12 mesh, `fldos='fonones.dos'`). Γ mode (`--gamma` or `--raman`): `dynmat.in` (`asr='simple'`).
3. `--raman` requires norm-conserving pseudopotentials (`p["type"] == "NC"`), otherwise a usage error.
4. `--run`: `runner.run_all` runs `pw.x`; `phonons.run_chain` executes `ph.x` → (`dynmat.x` | `q2r.x` → `matdyn.x` ×2), skipping steps whose `.out` already says `JOB DONE`.
5. `phonons.collect`: at Γ it reads the table `# mode [cm-1] [THz] IR [Raman depol]` from `dynmat.out` (`read_dynmat_table`); on a mesh it reads `bandas.freq` (`_read_flfrq`, `&plot` format, q in Cartesian 2π/alat) and `fonones.dos`.
6. `phonons.report` / `report_gamma_activities` (mutual exclusion rule, depolarisation 0.75), `phonons.thermodynamics`, `phonons.export` (`FONONES_GAMMA.dat` or `FONONES_BANDAS.dat`, `FONONES_DOS.dat`, `FONONES_TERMO.dat`), and `phonons.plot` only if `phonons.has_dispersion(run)` (there are `band_freqs` and `qdist`, i.e. it is not a Γ-only run).
7. `--tscan T1,T2,...` → `qekit/cli.py: _cmd_phonons_tscan` → `qekit/modules/tphonons.py: prepare`: one full chain per temperature in `T00300/` etc., with `insulator=False`, `smearing='fermi-dirac'` and `degauss = k_B T`; `tphonons.collect`, `report` (table of imaginary modes per T, monotonicity, $T_{\mathrm{stab}}$), `export` (`FONONES_T.dat`), `plot`.

**Where each datum comes from.**

| Datum | Source | Detail |
|---|---|---|
| Γ frequencies, IR, Raman, depol | `dynmat.out` (`# mode` table) | `phonons.read_dynmat_table`; IR in (D/Å)²/amu, Raman in Å⁴/amu |
| Dispersion | `bandas.freq` from `matdyn.x` | `phonons._read_flfrq` |
| Phonon DOS | `fonones.dos` from `matdyn.x` | `np.loadtxt`, states/cm⁻¹ |
| High-symmetry path | seekpath (`kpoints.get_kpath`) | labels and discontinuities |
| $k_B$, cm⁻¹→eV, cm⁻¹→THz | `phonons.KB_EV`, `CM1_TO_EV`, `CM1_TO_THZ` | CODATA |
| $k_B$ in Ry/K | `tphonons.KB_RY` | 6.333621e-6 |
| Imaginary threshold | literal −5 cm⁻¹ (`phonons.report`), `tphonons.UMBRAL_IMAGINARIO` = 10 | numerical noise below that |

**Limits and pitfalls.** It is harmonic: no thermal expansion or anharmonicity (for that, `qha`). The `prepare` docstring defaults to `insulator=True`, but the CLI passes `args.insulator`, which is `False` unless `--insulator`: by default the scf uses smearing and `epsil` is NOT enabled (no LO–TO splitting). It warns: "there are imaginary (negative) frequencies. Either the structure is not relaxed, or it is unstable at Γ" and "the structure must be relaxed (vc-relax) with these same cutoffs". The absolute scale of the `matdyn` DOS depends on the mesh; it is renormalised to $3N$ in the thermodynamics. With `--raman` (which forces Γ mode even without `--gamma`) no dispersion is drawn: the CLI decides with `has_dispersion(run)`, not with the flag. The thermodynamics ignores $\omega \le 1$ cm⁻¹. In `--tscan`, "this is ELECTRONIC temperature. The ions stay still"; if the number of imaginary modes does not decrease monotonically: "Usually means the k-mesh is not converged". Only `smearing='fermi-dirac'` corresponds to a real temperature.

**References.** S. Baroni, S. de Gironcoli, A. Dal Corso and P. Giannozzi, *Rev. Mod. Phys.* 73, 515 (2001), DOI 10.1103/RevModPhys.73.515 (DFPT). M. Lazzeri and F. Mauri, *Phys. Rev. Lett.* 90, 036401 (2003), DOI 10.1103/PhysRevLett.90.036401 (Raman via 2n+1). D. Porezag and M. R. Pederson, *Phys. Rev. B* 54, 7830 (1996) (Raman intensities). A. A. Maradudin et al., *Theory of Lattice Dynamics in the Harmonic Approximation*, Academic Press (1971).

---

### `olla-dft qha` — Quasi-harmonic approximation

**What it answers.** How does the crystal expand with temperature ($V(T)$, $\alpha(T)$, $a(T)$), what is the Grüneisen parameter, $C_p$ versus $C_v$ and $B(T)$?

**Background for non-experts.** The harmonic approximation gives no expansion: if the frequencies do not depend on volume, the minimum of the free energy does not move. The QHA keeps the harmonic modes but lets their frequencies change with VOLUME. At each temperature the static energy $E(V)$ and the vibrational free energy $F_{\mathrm{vib}}(V,T)$ are added, and the volume minimising the sum is found. On heating, the vibrational free energy favours larger volumes (modes soften) and the minimum shifts: that is thermal expansion.

**Formulas.** (`qha.f_vib`, `qha.cv_modos`, `qha.run`) For each volume $V_i$ with its modes $\omega_k$ (cm⁻¹, only $\omega > 1$), $\epsilon_k = \hbar\omega_k$:

$$
F(V_i,T) = E(V_i) + \frac{1}{N_{\mathrm{cells}}}\left[\sum_k \tfrac{1}{2}\epsilon_k + k_B T \sum_k \ln\!\left(1 - e^{-\epsilon_k/k_BT}\right)\right]
$$

Minimum by a local parabola (up to 2 points on each side of the sampled minimum): $V(T) = -b/2a$ (clipped to the range), $B(T) = 2a\,V(T)\times 160.21766208$ GPa.

$$
\alpha(T) = \frac{1}{V}\frac{dV}{dT}\ (\texttt{np.gradient}), \qquad
C_v = k_B\sum_k x_k^2 \frac{e^{x_k}}{(e^{x_k}-1)^2}\ \text{(interpolated at } V(T)), \qquad
C_p = C_v + \alpha^2 B V T
$$

$$
\gamma = -\frac{d\ln\langle\omega\rangle}{d\ln V}\ \text{(straight line over } \ln V\text{, } \ln\bar\omega\text{)}, \qquad
a(T) = \begin{cases} (V \cdot V_{\mathrm{conv}}/V_{\mathrm{prim}})^{1/3} & \text{with } \texttt{--structure} \\ V_{\mathrm{prim}}^{1/3} & \text{without it} \end{cases}\ (\texttt{--cubic})
$$

- $E$ in eV, $V$ in Å³, $C_v$, $C_p$ in meV/K per cell, $\alpha$ in K⁻¹, $B$ in GPa; `KB_EV` = 8.617333262e-5, `CM1_EV` = 1.239841984e-4.

**How Olla-DFT computes it.**
1. `qekit/cli.py: _cmd_qha` reads a TABLE (`data`): columns $V$ (Å³), $E$ (eV), $\omega_1, \omega_2, \ldots$ (cm⁻¹) per volume; values ≤ −1000 are discarded as padding.
2. With `--structure CIF`, `qha.factor_convencional` counts how many primitive cells fit in the conventional one ($N_{\mathrm{conv}}/N_{\mathrm{prim}}$ via spglib: 4 for fcc/diamond, 2 for bcc, 1 for simple cubic) and `qha.es_cubico` (space group ≥ 195) turns on `--cubic` automatically.
3. `qekit/modules/qha.py: run` with $T$ = 0 … `--tmax` (1000) in steps of `--dt` (5), `--natoms` (1), `--cells` (primitive cells per supercell of the modes, 1), `--cubic`, `factor_conv`.
4. Warnings if there are < 4 volumes or frequencies < −5 cm⁻¹ at some volume.
5. `qha.report` (at `--temp` 300 K; $a(T)$ labelled "lattice parameter (conventional cell)" or "V_prim^(1/3) (NOT the conventional lattice parameter)" according to `QHAResult.a_convencional`), `qha.export` (`QHA.dat`), `qha.plot` (V, α, $C_v$/$C_p$).

**Where each datum comes from.**

| Datum | Source | Detail |
|---|---|---|
| $E(V)$ and $\omega_k(V)$ | user table (`data`) | from `eos` + `phonons` (or `mlip phonons`) per volume; Olla-DFT does not generate it |
| $V_{\mathrm{conv}}/V_{\mathrm{prim}}$ and cubicity | `--structure` via spglib | `qha.factor_convencional`, `qha.es_cubico` |
| Constants | `qha.KB_EV`, `CM1_EV`, `EV_A3_GPA` | CODATA; 160.21766208 |

**Limits and pitfalls.** It launches no QE calculation: it receives the table. It uses discrete frequencies (one set of modes per volume), not a DOS: the thermodynamics is done on the list it is given, so the quality depends on the supercell/mesh of those modes. The Grüneisen parameter is an average over the mean frequency, not mode by mode. The QHA "holds up to ~half the melting temperature". Without `--structure`, $a(T)$ is only $V_{\mathrm{prim}}^{1/3}$ and the report warns: "In fcc, bcc or diamond that is NOT the conventional lattice parameter (they differ by 4^(1/3) or 2^(1/3)). Pass the structure with --structure". With a single temperature, $\alpha$ is NaN and a warning is issued.

**References.** A. Togo, L. Chaput, I. Tanaka and G. Hug, *Phys. Rev. B* 81, 174301 (2010), DOI 10.1103/PhysRevB.81.174301. G. Grimvall, *Thermophysical Properties of Materials*, North-Holland (1999). P. Pavone et al., *Phys. Rev. B* 48, 3156 (1993) (Si, negative expansion).

---

### `olla-dft derived` — Debye, sound velocities and Slack from the $C_{ij}$

**What it answers.** From an already computed elastic matrix: what are the density, the sound velocities, the elastic Debye temperature, the Poisson ratio, an approximate Grüneisen parameter and an estimate of the lattice thermal conductivity?

**Background for non-experts.** The stiffness of a solid sets how fast sound waves travel through it, and that speed sets the fastest possible vibration; the Debye temperature is that maximum frequency expressed in kelvin. With it and an anharmonicity parameter (Grüneisen), the Slack model estimates how much heat the lattice conducts. Everything is post-processing: it costs no new calculation.

**Formulas.** (`derived.density`, `sound_velocities`, `debye_from_velocity`, `poisson_ratio`, `gruneisen_from_poisson`, `slack`, `cubic_directional`)

$$
\rho = \frac{\sum_i m_i}{V}, \qquad
v_l = \sqrt{\frac{B + \tfrac{4}{3}G}{\rho}}, \qquad v_t = \sqrt{\frac{G}{\rho}}, \qquad
v_m = \left[\frac{1}{3}\left(\frac{2}{v_t^3} + \frac{1}{v_l^3}\right)\right]^{-1/3}
$$

$$
\Theta_D = \frac{\hbar}{k_B}\left(6\pi^2 n\right)^{1/3} v_m, \qquad
\nu = \frac{3B - 2G}{2(3B+G)}, \qquad
\gamma = \frac{3(1+\nu)}{2(2-3\nu)}
$$

$$
\kappa_L = A\,\frac{\bar M\,\Theta_D^3\,\delta}{\gamma^2\, n_{\mathrm{at}}^{2/3}\, T}, \qquad
A = \frac{3.1\times10^{-6}}{1 - 0.514/\gamma + 0.228/\gamma^2}, \qquad
v_L^{[100]} = \sqrt{C_{11}/\rho},\ v_T^{[100]} = \sqrt{C_{44}/\rho}
$$

- $B$, $G$: Hill averages in GPa (×10⁹ to Pa); $\rho$ in kg/m³ (masses in amu × 1.66053906660e-27, $V$ in Å³ × 1e-30); $n$: atoms per m³; $\hbar$ = 1.054571817e-34 J·s, $k_B$ = 1.380649e-23 J/K; $\bar M$: mean mass in amu; $\delta = (V/n_{\mathrm{at}})^{1/3}$ in Å; $T$ in K (`--temp`, 300); $\kappa_L$ in W/(m·K).

**How Olla-DFT computes it.**
1. `qekit/cli.py: _cmd_derived` loads the structure (masses and volume) and `--cij` (`ELASTIC_C.dat` from `elastic`, 6×6 matrix).
2. `elastic.moduli` → $B_H$, $G_H$; `qekit/modules/derived.py: analyze` computes everything above.
3. `derived.cubic_directional` prints $v_L$, $v_T$ along [100] only if the structure is cubic according to spglib (`elastic.crystal_family`) or the tensor has cubic form (`derived.is_cubic_tensor`: $C_{11}=C_{22}=C_{33}$, $C_{12}=C_{13}=C_{23}$, $C_{44}=C_{55}=C_{66}$ and zeros elsewhere, with a 5 % or 2 GPa tolerance).
4. `derived.report`; `derived.export` writes `DERIVED.dat`.

**Where each datum comes from.**

| Datum | Source | Detail |
|---|---|---|
| $C_{ij}$ | `ELASTIC_C.dat` (`--cij`) | `np.loadtxt`, GPa |
| Masses and volume | structure (ASE `get_masses`, `get_volume`) | amu, Å³ |
| $\hbar$, $k_B$, amu | `derived.HBAR`, `KB`, `AMU` | CODATA 2018 |
| Slack prefactor | literal 3.1e-6 and correction $(1 - 0.514/\gamma + 0.228/\gamma^2)$ | Slack / Julian |

**Limits and pitfalls.** The $\Theta_D$ is the ELASTIC one (low-temperature acoustic limit): "The one from the phonon DOS uses the whole spectrum and gives another number; they are not the same quantity" (`derived.debye_from_dos`, $\Theta_D = (\hbar/k_B)\sqrt{5\langle\omega^2\rangle/3}$, exists and is used by `crosscheck`, not by this command). The Grüneisen parameter "comes from an empirical correlation with the Poisson ratio" (Belomestnykh) and Slack "is an order-of-magnitude estimate". Negative Poisson: auxetic-material warning. The Slack κ is labelled with the temperature actually used (`Termoelastico.T`, key `kappa_Slack_<T>K` in `DERIVED.dat`). If $G \le 0$ there are no velocities.

**References.** O. L. Anderson, *J. Phys. Chem. Solids* 24, 909 (1963), DOI 10.1016/0022-3697(63)90067-2 (elastic $\Theta_D$). G. A. Slack, *Solid State Phys.* 34, 1 (1979); D. T. Morelli and G. A. Slack, in *High Thermal Conductivity Materials*, Springer (2006) (prefactor with Julian's correction). V. N. Belomestnykh and E. P. Tesleva, *Tech. Phys.* 49, 1098 (2004) (Grüneisen–Poisson).

---

### `olla-dft thermochem` — ZPE, entropy and free energy

**What it answers.** How much must be added to a DFT energy (at 0 K, without vibrations) to obtain a free energy $G(T,p)$ comparable with experiment, for a solid, an adsorbate, an ideal gas or a transition state?

**Background for non-experts.** A DFT energy is electronic and at 0 K. What is measured is a free energy at the laboratory temperature and pressure. Between the two there are three terms: the zero-point energy (modes vibrate even at 0 K), the enthalpy correction (modes get populated on heating) and the entropic term $-TS$, which for a gas molecule includes the translational (Sackur–Tetrode) and rotational (rigid rotor) entropies and can amount to about 1 eV at 500 K. Forgetting it can flip the sign of an adsorption energy.

**Formulas.** (`thermochem.zpe`, `H_vib`, `S_vib`, `Cv_vib`, `S_traslacional`, `S_rotacional`, `corregir`) With $\epsilon_k = h c\,\tilde\nu_k$, $x_k = \epsilon_k/k_BT$ (capped at 500):

$$
E_{\mathrm{ZPE}} = \tfrac{1}{2}\sum_k \epsilon_k, \quad
H_{\mathrm{vib}} = \sum_k \frac{\epsilon_k}{e^{x_k}-1}, \quad
S_{\mathrm{vib}} = k_B\sum_k\left[\frac{x_k}{e^{x_k}-1} - \ln\!\left(1-e^{-x_k}\right)\right], \quad
C_v = k_B\sum_k \frac{x_k^2 e^{x_k}}{(e^{x_k}-1)^2}
$$

$$
S_{\mathrm{trans}} = k_B\left[\ln\!\left(\frac{V}{\Lambda^3}\right) + \tfrac{5}{2}\right], \quad V = \frac{k_BT}{p}, \quad \Lambda = \frac{h}{\sqrt{2\pi m k_B T}}
$$

$$
S_{\mathrm{rot}}^{\mathrm{linear}} = k_B\left[\ln\frac{T}{\sigma\Theta_r} + 1\right], \qquad
S_{\mathrm{rot}}^{\mathrm{non\,linear}} = k_B\left[\tfrac{1}{2}\ln\frac{\pi T^3}{\sigma^2\Theta_A\Theta_B\Theta_C} + \tfrac{3}{2}\right], \qquad \Theta_i = \frac{\hbar^2}{2 I_i k_B}
$$

$$
G - E_{\mathrm{DFT}} = E_{\mathrm{ZPE}} + H_{\mathrm{corr}} - TS, \qquad
H_{\mathrm{corr}}^{\mathrm{gas}} = H_{\mathrm{vib}} + \left(\tfrac{3}{2} + n_{\mathrm{rot}} + 1\right)k_BT, \qquad S_{\mathrm{elec}} = k_B\ln(\text{multiplicity})
$$

- $\tilde\nu$ in cm⁻¹ (`C_CM` = 2.99792458e10 cm/s, `H_EVS` = 4.135667696e-15 eV·s); $m$: molecular mass (amu → kg); $p$ in Pa (`--pressure` in bar); $\sigma$: symmetry number (`--symmetry`, 1); $I_i$: principal moments of inertia (amu·Å² → kg·m²); linear if $I_1 < 10^{-3} I_3$; $n_{\mathrm{rot}}$ = 1 (linear), 1.5 (non-linear), 0 (atom). Adsorption energy (`thermochem.adsorcion`): $E_{\mathrm{ads}} = E_{\mathrm{slab+ads}} - E_{\mathrm{slab}} - nE_{\mathrm{gas}}$; $G_{\mathrm{ads}} = E_{\mathrm{ads}} + G^{\mathrm{corr}}_{\mathrm{ads}} - n\,G^{\mathrm{corr}}_{\mathrm{gas}}$.

**How Olla-DFT computes it.**
1. `qekit/cli.py: _cmd_thermochem` reads the frequencies (`_leer_frecuencias`: a one- or multi-column file —last column—, or an inline list) and, for gas, the structure with `ase.io.read`.
2. `qekit/modules/thermochem.py: limpiar_frecuencias`: separates imaginary ones ($\tilde\nu < -1$), discards $|\tilde\nu| \le 1$ (residual translations/rotations), raises soft modes to the `--floor` (e.g. 100 cm⁻¹) and emits warnings according to `--phase` (`solido`, `adsorbato`, `gas`, `transicion`).
3. `thermochem.corregir` sums the terms at `--temp` (298.15 K) and `--pressure` (1 bar), with `--multiplicity`.
4. `thermochem.report` (with `G(T)` if `--energy` is given); with `-o`, it writes `TERMOQUIMICA.txt`.

**Where each datum comes from.**

| Datum | Source | Detail |
|---|---|---|
| Frequencies | file (`FONONES_GAMMA.dat`, last column) or list | `_leer_frecuencias` |
| Masses and geometry (gas) | `--structure` via ASE | moments of inertia |
| $h$, $k_B$, $c$, amu, $\hbar$ | `thermochem.H_EVS`, `KB_EV`, `C_CM`, `AMU_KG`, `HBAR_JS`, `KB_J` | CODATA |
| Soft-mode floor | `--floor` (`PISO_BLANDO` = 100 only as a reference) | without `--floor` nothing is raised |

**Limits and pitfalls.** Verified in the tests against NIST (H₂O, N₂, CH₄ to 0.5 %). Transition state without an imaginary mode or with more than one: explicit warning. Imaginary modes at a minimum: "the structure is NOT a minimum … They are excluded from the sums". Raised modes: "it is a CORRECTION, not a calculation: say so if you publish". Gas with a number of modes ≠ $3N-6$ (or $3N-5$): "they are counted twice with the translational and rotational terms". In the gas phase no anharmonicity or internal rotors are included; the solid carries no $pV$ term. `adsorcion` applies no vibrational corrections to the clean slab (assumes they do not change). The CLI `thermochem` command does not expose `adsorcion` (used from `adsorb`/API).

**References.** D. A. McQuarrie, *Statistical Mechanics*, University Science Books (2000). C. J. Cramer, *Essentials of Computational Chemistry*, Wiley (2004), ch. 10. O. Sackur, *Ann. Phys.* 36, 958 (1911); H. Tetrode, *Ann. Phys.* 38, 434 (1912).

---

### `olla-dft md` — Analysis of a molecular-dynamics trajectory

**What it answers.** From a `pw.x` output with `calculation='md'`: what structure does the system have (g(r), coordination numbers), do the atoms diffuse (MSD and diffusion coefficient $D$) and what is its vibrational spectrum (VDOS) including temperature and anharmonicity?

**Background for non-experts.** A molecular dynamics run is a "movie" of the atoms moving. Three functions summarise the movie: the radial distribution function $g(r)$ says how many neighbours there are at each distance (its first peak is the bond length, its area the coordination number); the mean-square displacement (MSD) says whether atoms move away from where they were (if it grows linearly, they diffuse; if it flattens, they only vibrate); and the Fourier transform of the velocity autocorrelation gives the frequencies at which they vibrate. Before trusting any of them, the initial equilibration stretch must be discarded and the trajectory must be long enough.

**Formulas.** (`dynamics.rdf`, `coordinacion`, `msd`, `difusion`, `vdos`) With minimum-image distances in fractional coordinates:

$$
g(r) = \frac{h(r)}{N_{\mathrm{steps}}\,\frac{N(N-1)}{2}\,\frac{4\pi r^2\,\Delta r}{V}}, \qquad
g_{AB}(r) = \frac{h_{AB}(r)}{N_{\mathrm{steps}}\,N_{\mathrm{pairs}}\,\frac{4\pi r^2\Delta r}{V}}, \qquad
n_{\mathrm{coord}} = \int_0^{r_{\min}} 4\pi r^2 \rho\, g(r)\, dr
$$

$$
\mathrm{MSD}(\tau) = \left\langle |\mathbf{r}_i(t+\tau) - \mathbf{r}_i(t)|^2 \right\rangle_{i,t}, \qquad
\mathrm{MSD} = 6 D \tau + b \;\Rightarrow\; D = \frac{m}{6}\times 10^{-1}\ [\mathrm{cm^2/s}]
$$

$$
C(t) = \frac{\sum_{i,\alpha}\langle v_{i\alpha}(0)v_{i\alpha}(t)\rangle}{C(0)}\ (\text{via FFT}), \qquad
\mathrm{VDOS}(\tilde\nu) = \left|\mathcal{F}\{C(t)\,w_{\mathrm{Hann}}(t)\}\right|, \quad \tilde\nu = \frac{f[\mathrm{fs^{-1}}]\times 10^{15}}{c[\mathrm{cm/s}]}
$$

- $r_{\max}$ = half the shortest cell edge (or `--rmax` if smaller); `--bins` 200; $N_{\mathrm{pairs}} = N_AN_B$ or $N_A(N_A-1)/2$; $\rho = N/V$; $r_{\min}$: first local minimum with $g<1$ after the first maximum; the $D$ fit uses only the 20–80 % stretch of the lags (lags up to $n/2$), slope $m$ in Å²/fs; velocities by `np.gradient` of the UNWRAPPED positions (no periodic jumps), not mass-weighted.

**How Olla-DFT computes it.**
1. `qekit/cli.py: _cmd_md` → `qekit/modules/dynamics.py: leer_md`: reads from `pw.out` (text) the `ATOMIC_POSITIONS` blocks (alat, bohr, angstrom or crystal → Å), the cell from `a(i) = (...)`·alat or `CELL_PARAMETERS`, `temperature = … K`, `!    total energy` and `Time step = … femto-seconds`; discards `--skip` steps.
2. `dynamics.analizar`: `rdf`, `coordinacion` per pair, `desdoblar` + `msd` (total and per species), `difusion`, `vdos` (≥ 8 steps), temperature drift between halves.
3. `dynamics.report`, `dynamics.export` (`MD_RDF.dat`, `MD_MSD.dat`, `MD_VDOS.dat`, `MD.txt`), `dynamics.plot` (three panels).

**Where each datum comes from.**

| Datum | Source | Detail |
|---|---|---|
| Positions per step | `pw.out` (`ATOMIC_POSITIONS`) | `dynamics._leer_marcos`; units detected |
| Cell | `pw.out` (`a(1..3) = (...)` × alat, or `CELL_PARAMETERS`) | assumed constant |
| Time step | `pw.out` (`Time step = … a.u., X femto-seconds`) | 1 fs if absent |
| Temperature and energy | `pw.out` (`temperature =`, `!    total energy`) | K; Ry → eV with 13.605693122994 |
| bohr → Å | `dynamics.BOHR_A` | 0.529177210903 |

**Limits and pitfalls.** Constant cell: not usable for `vc-md`. If $g(r)$ is empty up to the cutoff: "the cell is too small to extract structure from it: build a supercell". Less than 2 ps: "Good for looking at the structure, not for a diffusion coefficient". $R^2 < 0.95$: "the MSD is NOT linear: no diffusion, or not enough time". Temperature drift > 15 %: "it is still equilibrating: discard more steps with --skip". The coordination number is integrated up to the FIRST MINIMUM, "a convention, not a measurement". The VDOS is not mass-weighted (it is not the phonon DOS) and takes the modulus of the spectrum, not the real part; its resolution is $1/(N_{\mathrm{steps}}\,\Delta t)$. `KB_RY` in `dynamics.py` is unused.

**References.** M. P. Allen and D. J. Tildesley, *Computer Simulation of Liquids*, Oxford (2017). A. Einstein, *Ann. Phys.* 17, 549 (1905). J.-P. Hansen and I. R. McDonald, *Theory of Simple Liquids*, Academic Press (2013).

---

### `olla-dft kappa` — Lattice thermal conductivity (fc3 + BTE with phono3py)

**What it answers.** How much heat does the crystal lattice conduct, $\kappa_L(T)$ in W/(m·K), which exponent does it follow with temperature and which phonon mean free paths carry the heat (to know whether nanostructuring helps)?

**Background for non-experts.** In a perfectly harmonic crystal a phonon would travel forever and the conductivity would be infinite. What makes it finite is that a phonon can split into two or two can merge into one: that is allowed by the cubic term of the energy (the third-order force constants, fc3). They are obtained by displacing two atoms at a time in a supercell and computing the forces; with them, the phonon Boltzmann equation (in the relaxation-time approximation, RTA) gives $\kappa$. It is expensive because the number of configurations grows fast with the supercell. Olla-DFT allows the forces to be computed with `pw.x` (the real calculation) or with a learned potential (MACE, etc.) for exploration.

**Formulas.** They are solved by phono3py (`kappa.resolver`); Olla-DFT post-processes:

$$
\kappa_L^{\alpha\beta} = \frac{1}{NV}\sum_\lambda C_\lambda\, v_\lambda^\alpha v_\lambda^\beta\, \tau_\lambda, \qquad \tau_\lambda = \frac{1}{2\Gamma_\lambda}, \qquad \Lambda_\lambda = |\mathbf{v}_\lambda|\,\tau_\lambda
$$

$$
\bar\kappa = \frac{\kappa_{xx}+\kappa_{yy}+\kappa_{zz}}{3}, \qquad
\kappa \propto T^{-n}\ (n \text{ by a straight line in } \ln\kappa\text{–}\ln T,\ T \ge 200\ \mathrm{K}), \qquad
\kappa_{\mathrm{cum}}(\Lambda) = \frac{\sum_{\lambda:\Lambda_\lambda<\Lambda} w_\lambda C_\lambda \tfrac{|\mathbf{v}_\lambda|^2}{3}\tau_\lambda}{\sum_\lambda w_\lambda C_\lambda \tfrac{|\mathbf{v}_\lambda|^2}{3}\tau_\lambda}
$$

- $\Gamma_\lambda$: linewidth (THz) from phono3py; $\mathbf{v}_\lambda$: group velocity (THz·Å); $C_\lambda$: modal heat capacity; $w_\lambda$: q-point weight; $\Lambda$ in Å (reported in nm). Modes with $\Gamma = 0$ (acoustic at Γ) are discarded.

**How Olla-DFT computes it.**
1. `qekit/cli.py: _cmd_kappa` → `qekit/modules/kappa.py: preparar`: `Phono3py(..., supercell_matrix=--dim (2x2x2), phonon_supercell_matrix=--dim-fc2, primitive_matrix="auto", symprec=1e-5)` and `generate_displacements(distance=--distance 0.03 Å)`.
2. `kappa.configuraciones` converts the displaced supercells to ASE (fc3 and, if any, fc2).
3. Forces: (a) `--model mace|chgnet|m3gnet` → `kappa.fuerzas_mlip`; (b) without `--model` → `kappa.escribir_inputs` writes one `scf` per configuration in `fc3/dNNNN/pw.in` (and `fc2/`), `conv_thr = 1e-10`, mesh from `--kspacing` 0.35 Å⁻¹, `occupations='fixed'` unless `--metal` (smearing), plus `correr.sh`; it refuses above 150 configurations without `--force`; (c) `--collect` → `kappa.leer_fuerzas` reads `<forces>` from each XML (Ha/bohr → eV/Å) and requires ALL of them.
4. `kappa.resolver`: `produce_fc3`, `produce_fc2`, symmetrisation, `mesh_numbers = --mesh (13)`, `init_phph_interaction`, `run_thermal_conductivity(temperatures=--temps 100:800:8, is_isotope=--isotopes, boundary_mfp=--grain µm ×1e4 Å or 1e6)`.
5. `kappa.recoger` stores κ (Voigt 6), Γ, velocities, $C_\lambda$, weights; `kappa.report`, `export` (`KAPPA.dat`, `KAPPA_recorrido.dat`, `KAPPA.txt`), `plot` (κ(T) log-log with a $T^{-1}$ guide; cumulative vs Λ).

**Where each datum comes from.**

| Datum | Source | Detail |
|---|---|---|
| Forces | pw.x XML (`output/forces`) or MLIP potential | `qeout.read_xml` / `mlip.calculator` |
| fc2, fc3, Γ, v, C, κ | `phono3py` library | `Phono3py.thermal_conductivity` (RTA) |
| Mean-free-path grid | `kappa.RECORRIDOS` | `np.logspace(0, 7, 141)` Å |
| Isotopes | phono3py (natural abundances) | `--isotopes` |

**Limits and pitfalls.** "It is RTA, not the exact solution of the Boltzmann equation. RTA underestimates κ (≈10-15 % in silicon, much more in graphene or diamond)". "Only three-phonon scattering is included". Without `--isotopes`: "Natural silicon conducts ~10 % less than isotopically pure". With a learned potential: "the absolute value may be far off: with small MACE-MP silicon gives ~51 W/mK at 300 K where experiment is ~140". Supercell ≤ 8 cells: "κ must converge in the supercell size AND the q-mesh at the same time". By default the fc2/fc3 scf runs use `occupations='fixed'` ("the right choice for insulators"); for a metal `--metal` must be passed, or the scf runs will not converge. There is no option to write/read phono3py's `fc2.hdf5`/`fc3.hdf5`: every `--collect` rebuilds everything.

**References.** A. Togo, L. Chaput and I. Tanaka, *Phys. Rev. B* 91, 094306 (2015), DOI 10.1103/PhysRevB.91.094306 (phono3py). J. M. Ziman, *Electrons and Phonons*, Oxford (1960). L. Lindsay, D. A. Broido and T. L. Reinecke, *Phys. Rev. B* 87, 165201 (2013) (RTA vs. exact solution).

---

### `olla-dft elph` — Electron–phonon coupling: λ, ω_log, Tc and τ

**What it answers.** How strongly do electrons couple to phonons ($\lambda$), what is the Allen–Dynes superconducting critical temperature with its strong-coupling corrections, and what is the phonon-limited relaxation time $\tau(T)$ that the CRTA transport lacks?

**Background for non-experts.** Electrons moving through a metal collide with the lattice vibrations; how much they collide is measured by $\lambda$, a dimensionless number. In a conventional superconductor that same coupling is what pairs the electrons, and the Allen–Dynes formula turns $\lambda$ and a typical phonon frequency ($\omega_{\log}$) into a critical temperature. The same $\lambda$ gives the time between collisions at high temperature, $\tau$. `ph.x` computes the coupling for several numerical broadenings; the good value is the one at the "plateau", where it stops depending on the broadening.

**Formulas.** (`elph.lambda_de_a2F`, `omega_log_de_a2F`, `omega_2`, `factores_correccion`, `allen_dynes`, `tau_elph`)

$$
\lambda = 2\int \frac{\alpha^2F(\omega)}{\omega}\,d\omega, \qquad
\omega_{\log} = \exp\!\left[\frac{2}{\lambda}\int \ln\omega\,\frac{\alpha^2F(\omega)}{\omega}\,d\omega\right], \qquad
\bar\omega_2 = \left[\frac{2}{\lambda}\int \omega\,\alpha^2F(\omega)\,d\omega\right]^{1/2}
$$

$$
T_c = f_1 f_2\,\frac{\omega_{\log}}{1.2}\exp\!\left[\frac{-1.04(1+\lambda)}{\lambda - \mu^*(1+0.62\lambda)}\right], \qquad
f_1 = \left[1 + \left(\frac{\lambda}{\Lambda_1}\right)^{3/2}\right]^{1/3}, \quad \Lambda_1 = 2.46(1+3.8\mu^*)
$$

$$
f_2 = 1 + \frac{(r-1)\lambda^2}{\lambda^2 + \Lambda_2^2}, \quad r = \frac{\bar\omega_2}{\omega_{\log}}, \quad \Lambda_2 = 1.82(1+6.3\mu^*)\,r, \qquad
\frac{1}{\tau} = \frac{2\pi\lambda k_B T}{\hbar}
$$

- $\omega$ in THz in `a2F.dos*`; $\omega_{\log}$, $\bar\omega_2$ in K (`THZ_K` = 47.9924 K/THz); $\mu^*$ = 0.10, 0.13, 0.16 (a range, not computed); $T_c$ = 0 if the denominator ≤ 0; $\hbar$ = `HBAR_EVS` = 6.582119569e-16 eV·s, $k_B$ = 8.617333262e-5 eV/K; $\tau$ in s. Plateau (`elph.plato`): longest stretch of ≥ 3 consecutive λ that do not differ by more than 5 % from the first; its midpoint.

**How Olla-DFT computes it.**
1. Preparation (`qekit/cli.py: _cmd_elph` without `--collect` → `qekit/modules/elph.py: prepare`): writes `1_scf.in` (mesh `--kgrid` or `kspacing`), `2_nscf.in` with `la2F = .true.` and mesh `--kgrid-nscf` (default $q_i\cdot\max(2, \lceil 2k_i/q_i\rceil)$, a multiple of the q-mesh), and `3_ph.in` with `electron_phonon='interpolated'`, `el_ph_sigma = --sigma (0.005 Ry)`, `el_ph_nsigma = --nsigma (10)`, `fildvscf='dvscf'`, `tr2_ph = 1e-12`, `ldisp` with `--qgrid` (2x2x2). Smearing `methfessel-paxton`, `--degauss` 0.02 Ry, `conv_thr = 1e-10`.
2. The user runs `pw.x` ×2, `ph.x` and, optionally, `lambda.x` (Olla-DFT has `elph.build_lambda_input` but the CLI never writes it).
3. `--collect`: `elph.leer_elph_ph` reads from `ph.out` the broadenings (`Gaussian Broadening: X Ry`) and `DOS = … states/spin/Ry`; `elph.leer_lambda_out` reads `lambda.dat` (columns σ, λ, ∫α²F, ⟨log ω⟩, N(E_F)) or the text of `lambda.out`, takes $\mu^*$ from `lambda.in` (last numeric line) and fills the per-broadening $T_c$ column: from the final `lambda omega_log T_c` table of `lambda.out` if present and of matching size, and otherwise computes it row by row with `allen_dynes(λ_i, ω_log,i, μ*, correcciones=False)`; `ElPhRun.Tc_fuente` records which of the two was done and the report prints it. `elph.leer_a2F` reads `a2F.dos*` (or `A2F.dat`) and from it λ, $\omega_{\log}$ (if missing) and, always, $\bar\omega_2$ (`elph.omega_2`, computed in the CLI), which is what enables the $f_2$ factor in the summary $T_c$ values.
4. `elph.plato`, `elph.report` (table per broadening, regime, $T_c$ for three $\mu^*$, τ at 100/300/500/800 K), `elph.export` (`ELPH.dat`, `A2F.dat`, `ELPH.txt`), `elph.plot`.

**Where each datum comes from.**

| Datum | Source | Detail |
|---|---|---|
| Broadenings, N(E_F) | `ph.out` (`Gaussian Broadening`, `DOS =`) | `elph.leer_elph_ph` |
| λ, ⟨log ω⟩ per broadening | `lambda.dat` / `lambda.out` from `lambda.x` | `elph.leer_lambda_out` |
| $\alpha^2F(\omega)$ | `a2F.dos*` (or `A2F.dat`) | `elph.leer_a2F`, column 1 THz, column 2 α²F |
| THz → K | `elph.THZ_K` | 47.9924 |
| $\mu^*$ | literal (0.10, 0.13, 0.16) | empirical |
| $T_{\mathrm{Debye}}$ | `--debye` | only to mark the validity regime of τ |

**Limits and pitfalls.** The λ values in `ph.out` are NOT read (only σ and N(E_F)); λ comes from `lambda.x` or from the $\alpha^2F$. The "Tc(K)" column of the per-broadening table is the one from `lambda.x` (Allen–Dynes WITHOUT $f_1 f_2$, with the $\mu^*$ of `lambda.in`) or the same thing recomputed by Olla-DFT; the corrected $T_c$ values for the three $\mu^*$ are those of the "Critical temperature" block, and without `a2F.dos*` there is no $\bar\omega_2$ and $f_2 = 1$. No plateau: "the k-mesh is insufficient … Any lambda reported from here is arbitrary". $\mu^*$ "is empirical (0.10-0.16) and is NOT computed here". τ "holds ABOVE the Debye temperature; below it overestimates the scattering". `lambda.x` with a coarse q-mesh leaves $\omega_{\log}$ as NaN: it is recomputed from the $\alpha^2F$ if available. τ is NOT injected automatically into `transport`: the module docstring says so explicitly and gives the sequence (`transport --collect` → `elph --collect` → $\sigma(T) = [\sigma/\tau](T)\cdot\tau(T)$ by hand on the columns of `TRANSPORTE.dat`).

**References.** P. B. Allen and R. C. Dynes, *Phys. Rev. B* 12, 905 (1975), DOI 10.1103/PhysRevB.12.905. W. L. McMillan, *Phys. Rev.* 167, 331 (1968), DOI 10.1103/PhysRev.167.331. P. B. Allen, *Phys. Rev. B* 3, 305 (1971) (high-T τ). G. Grimvall, *The Electron–Phonon Interaction in Metals*, North-Holland (1981).

---

### `olla-dft transport` — Electronic transport in CRTA: Seebeck, σ/τ, κ_e/τ, Lorenz and spin

**What it answers.** From the bands on a dense mesh: what are the Seebeck coefficient $S$, the conductivity per relaxation time $\sigma/\tau$, the electronic thermal conductivity $\kappa_e/\tau$, the power factor $S^2\sigma/\tau$ and the carrier concentration as functions of chemical potential and temperature? Does Wiedemann–Franz hold? How is transport shared between the two spin channels?

**Background for non-experts.** An electron in a band moves with velocity $v = (1/\hbar)\,dE/dk$. At a given temperature only the states within a few $k_BT$ of the chemical potential take part in transport (the $-\partial f/\partial E$ "window"). Summing velocity times velocity over that window gives the conductivity; weighting additionally by $(E-\mu)$ gives the Seebeck coefficient, which measures how much voltage appears per degree of temperature difference. The constant relaxation-time approximation (CRTA) assumes all electrons collide at the same rate $1/\tau$: then $\tau$ cancels in $S$ and in the Lorenz number (real predictions) but not in σ or κ_e, which are reported divided by τ.

**Formulas.** (`transport._fd_derivative`, `transport.compute`, `lorenz`, `cancelacion`, `TransporteEspin`) With $x = (E-\mu)/k_BT$, $-\partial f/\partial E = \mathrm{sech}^2(x/2)/(4k_BT)$, weights $w_k = 1/N_k$, $V$ the cell volume:

$$
\mathbf{v}_{n\mathbf{k}} = \frac{1}{\hbar}\nabla_{\mathbf{k}}E_{n\mathbf{k}}\ (\text{periodic finite differences, } \texttt{np.gradient}), \qquad
\mathbf{s}_m = \sum_{n\mathbf{k}} w_k\, \mathbf{v}\otimes\mathbf{v}\,(E-\mu)^m\left(-\frac{\partial f}{\partial E}\right)
$$

$$
\frac{\boldsymbol\sigma}{\tau} = \frac{e}{V}\mathbf{s}_0, \qquad
\mathbf{S} = -\frac{1}{T}\,\mathbf{s}_1\mathbf{s}_0^{-1}, \qquad
\frac{\boldsymbol\kappa_e}{\tau} = \frac{e}{VT}\mathbf{s}_2 - \mathbf{S}\mathbf{S}\,\frac{\boldsymbol\sigma}{\tau}\,T, \qquad
\mathrm{PF} = \bar S^2\,\bar\sigma/\tau
$$

$$
n = \frac{N_{\mathrm{elec}} - 2\sum_{n\mathbf{k}} w_k f(E_{n\mathbf{k}})}{V}, \qquad
L = \frac{\bar\kappa_e}{\bar\sigma T}, \qquad L_0 = 2.44\times10^{-8}\ \mathrm{W\,\Omega/K^2}, \qquad
c = \frac{|\bar\kappa_e|}{|\bar\kappa_e + \bar S^2\bar\sigma T|}
$$

$$
\sigma_{\mathrm{tot}} = \sigma_\uparrow + \sigma_\downarrow, \qquad
S_{\mathrm{tot}} = \frac{S_\uparrow\sigma_\uparrow + S_\downarrow\sigma_\downarrow}{\sigma_\uparrow+\sigma_\downarrow}, \qquad
P = \frac{\sigma_\uparrow-\sigma_\downarrow}{\sigma_\uparrow+\sigma_\downarrow}, \qquad S_{\mathrm{spin}} = S_\uparrow - S_\downarrow
$$

- $e$ = 1.602176634e-19 C; $\hbar$ = 6.582119569e-16 eV·s; $k_B$ = 8.617333262e-5 eV/K; σ/τ in S/(m·s); $S$ in V/K (µV/K in the report); κ_e/τ in W/(m·K·s); $n$ in cm⁻³ (positive = holes); bar = trace/3. The "cancellation" $c$ measures which fraction survives the subtraction $\kappa^0 - S^2\sigma T$.

**How Olla-DFT computes it.**
1. `qekit/cli.py: _cmd_transport` → `qekit/modules/transport.py: prepare`: standardised primitive cell; `scf.in` (mesh from `--kspacing` or the configured one) and `nscf.in` with `K_POINTS automatic` `--grid` (16x16x16), `nosym=.true.` (full mesh), `nbnd = 2 × estimated nbnd`; `--metal` turns off `occupations='fixed'`; `--nspin 2` and `--mag EL=value` (which implies `nspin=2`) write scf and nscf with spin polarisation, required for `--spin-resolved`.
2. `--run` executes scf and nscf; `--collect` → `transport.load` reads the first `out/*.xml`, rebuilds the grid from the fractional coordinates (rejects anything that is not a complete grid) and differentiates $E(\mathbf{k})$ with `np.gradient` on the periodically wrapped mesh, converting to Cartesian with $\mathbf{B}^{-T}$.
3. `transport.compute` over $T$ = `--temperatures` (300) and 201 values of µ in $E_F \pm$ `--mu-span` (1 eV).
4. `transport.report` (best p- and n-type $S$, maximum PF), `report_lorenz`, with `--spin-resolved` it loads `spin=1` and `report_espin`; `transport.export` (`TRANSPORTE.dat`), `transport.plot`.

**Where each datum comes from.**

| Datum | Source | Detail |
|---|---|---|
| Eigenvalues, k, weights | nscf XML (`ks_energies`, Hartree → eV) | `qeout.read_xml`; `weights` replaced by $1/N_k$ |
| $E_F$, $N_{\mathrm{elec}}$, volume, cell | XML (`fermi_energy`, `nelec`, `cell`) | `qeout.read_xml` |
| Constants | `transport.E_CHARGE`, `HBAR_EVS`, `KB_EV`, `L0_SOMMERFELD` | CODATA; $L_0$ = 2.44e-8 |

**Limits and pitfalls.** CRTA only: "To give σ in S/m you need a τ that comes from a fit to a measurement or from an electron–phonon calculation — Olla-DFT does not invent it". It does NOT compute ZT (that would need κ_lattice and τ) nor couple the τ from `elph` automatically. It does not interpolate bands (unlike BoltzTraP): with a mesh < 24 per side or < 12000 points it warns "INSUFFICIENT … sigma comes out as isolated spikes". The Lorenz number inside the gap suffers catastrophic cancellation: "DO NOT TRUST THIS NUMBER … X % survives"; only points with $c > 0.10$ are summarised. `--spin-resolved` on an XML with `nspin = 1` is rejected with the instruction to re-prepare with `--nspin 2 --mag EL=0.7 --run`. Two-current model: "Valid as long as spin-flip scattering is slow … it stops holding near [the Curie temperature]".

**References.** G. K. H. Madsen and D. J. Singh, *Comput. Phys. Commun.* 175, 67 (2006), DOI 10.1016/j.cpc.2006.03.007 (BoltzTraP, same CRTA formulation). N. W. Ashcroft and N. D. Mermin, *Solid State Physics*, ch. 13 (Wiedemann–Franz). N. F. Mott, *Proc. R. Soc. A* 153, 699 (1936) (two-current model).

---

### `olla-dft ballistic` — Landauer conductance with `pwcond.x`

**What it answers.** How many conduction channels does an electrode have at each energy (complex bands) and what transmission probability $T(E)$ does a nanocontact or a molecule between electrodes have? What is its conductance in units of $G_0$?

**Background for non-experts.** In a macroscopic crystal the electron collides many times (diffusive transport, `transport`). In a few-atom contact it crosses in one go: there is no conductivity, there is conductance, given by the Landauer formula: $G = G_0 T(E_F)$, with $T$ the probability of passing summed over all open "lanes". Since $T$ cannot exceed the number of lanes, the conductance comes out quantised in steps of $G_0 = 2e^2/h$, and seeing those steps is the sign that the calculation is right.

**Formulas.** (`ballistic.G0`, `CondRun`)

$$
G = G_0\,T(E_F), \qquad G_0 = \frac{2e^2}{h} = 7.748091729\times10^{-5}\ \mathrm{S}, \qquad R = \frac{1}{G} = \frac{12.906\ \mathrm{k\Omega}}{T(E_F)}, \qquad T(E) \le N_{\mathrm{channels}}(E)
$$

- $T(E_F)$: transmission at the energy closest to $E - E_F = 0$ in the window. Region limits (`ballistic.longitud_z`): $\mathrm{bdl} = |\mathbf{a}_3|_{\mathrm{electrode}}/a$, $\mathrm{bds} = |\mathbf{a}_3|_{\mathrm{scatterer}}/a$ with $a = |\mathbf{a}_1|$ (alat units): the boundary of each region is the end of ITS cell, not the height of the last atom.

**How Olla-DFT computes it.**
1. `qekit/cli.py: _cmd_ballistic` → `qekit/modules/ballistic.py: prepare`: `comprobar_geometria` requires $\mathbf{a}_3 \parallel z$, $\mathbf{a}_{1,2} \perp z$ and, with `--scatterer`, the same in-plane cell.
2. Writes `scf_electrodo.in` (and `scf_dispersor.in`) with `insulator=False` and prefixes `electr`/`disper`, and `cond.in` (`&inputcond`: `ikind` = `--ikind` (only 0 or 1) or 1 if there is a scatterer / 0 if not; `ikind=1` without `--scatterer` and `ikind=2` are rejected with an explanation; `energy0 = --emax` (3), `denergy = -(emax-emin)/(n-1)` (`--emin` −3, `--points` 61), `ewind = 1`, `epsproj = 1e-3`, `nz1 = --nz1` (3), `bdl` = `longitud_z(electrode)`, `bds` = `longitud_z(scatterer)`, one k-point (0, 0, 1)).
3. The user runs `pw.x` and `pwcond.x`; `--collect` → `ballistic.collect` reads `trans*.dat` (E, T) or the `T_tot` lines of `cond*.out`, and `Nchannels of the left tip` per energy (maximum over k); `ikind` from the `.out` or from `cond.in`.
4. `ballistic._avisar`, `report`, `export` (`BALISTICO.dat`, `.txt`), `plot` (T(E) and channel steps).

**Where each datum comes from.**

| Datum | Source | Detail |
|---|---|---|
| $T(E)$ | `trans.dat` from `pwcond.x` (or `T_tot` in `cond.out`) | `ballistic.collect` |
| Open channels | `cond.out` (`Nchannels of the left tip`) | maximum per energy |
| $G_0$ | `ballistic.G0` | 7.748091729e-5 S (CODATA) |
| bdl/bds limits | cell length along $z$ (`longitud_z`) | in alat units |

**Limits and pitfalls.** `ikind=0` "is NOT the conductance. It is the number of open channels, which bounds the conductance from above". If $T > N$: "That is impossible: T <= N by construction. Check that the bdl/bds limits…". Negative transmissions: "the calculation did not converge or … the geometry of the regions is badly cut". Different electrodes (`ikind=2` of `pwcond.x`) are not supported: `--ikind` only accepts 0 and 1, and asking for 2 gives "not implemented … the third scf and 'prefixr' and 'bdr' in cond.in must be written by hand". A single transverse k-point by default (0, 0, 1). There is no `--run`: it is always run by hand.

**References.** R. Landauer, *IBM J. Res. Dev.* 1, 223 (1957); M. Büttiker, *Phys. Rev. Lett.* 57, 1761 (1986), DOI 10.1103/PhysRevLett.57.1761. H. J. Choi and J. Ihm, *Phys. Rev. B* 59, 2267 (1999), DOI 10.1103/PhysRevB.59.2267; A. Smogunov, A. Dal Corso and E. Tosatti, *Phys. Rev. B* 70, 045417 (2004) (`pwcond.x`).

---

### `olla-dft cost` — Cost estimator calibrated with your history

**What it answers.** How long will this sweep take on THIS machine, and with what uncertainty? (`cost` shows the model; `--estimate` on any sweep applies it.)

**Background for non-experts.** The time of a plane-wave calculation scales in a known way with the number of k-points, plane waves, bands and iterations. What is not known in advance is the proportionality constant of each machine. Olla-DFT takes the shape from the physics and fits the scale with the calculations the user has already indexed in `olla-dft db` (with their wall times), and measures how wrong it is by leaving one system out and predicting it with the others.

**Formulas.** (`cost.n_ondas_planas`, `trabajo`, `iteraciones`, `_ajusta`, `estimar`)

$$
N_{\mathrm{PW}} = \frac{V\,E_{\mathrm{cut}}^{3/2}}{6\pi^2}\ (\text{bohr}^3,\ \mathrm{Ry}), \qquad
w_1 = n_k\, s\, N_{\mathrm{PW}}\, n_{\mathrm{b}}, \qquad
w_2 = n_k\, s\, N_{\mathrm{PW}}\, n_{\mathrm{b}}^2
$$

$$
t = t_0 + \left(C_1 w_1 + C_2 w_2\right)\, n_{\mathrm{scf}}\, n_{\mathrm{ion}}, \qquad
\text{NNLS fit with weights } t^{-1/2},\qquad
[t/\mathrm{disp},\ t\cdot\mathrm{disp}],\ \mathrm{disp} = e^{\sigma(\ln(\mathrm{pred}/\mathrm{real}))}
$$

- $V$: volume (Å³ → bohr³ with `A3_BOHR3`); $n_k$: irreducible k-points (spglib `get_ir_reciprocal_mesh`, or the real `number of k points` from a `pw.out`); $s$: 1, 2 or 4 (`nspin`, non-collinear); $n_{\mathrm{b}}$: `nbnd` from the input or $\max(4, 2N_{\mathrm{at}})$; $n_{\mathrm{scf}}$: median of the history (14 by default); $n_{\mathrm{ion}}$: median per type (`relax` 8, `vc-relax` 12 by default). The three-coefficient fit only if there are ≥ 8 calculations and $w_1^{\max}/w_1^{\min} \ge 5$; otherwise $C_1$ = geometric median of $t/w_1$.

**How Olla-DFT computes it.**
1. `qekit/cli.py: _cmd_cost` → `qekit/modules/cost.py: calibrar(--db olla-dft.db)`: `cost.historial` queries the `calculos` table (natoms, ecutwfc, kgrid, nspin, volume, n_scf, nk, nbnd, n_bfgs, wall_s, calculation).
2. `cost._prepara` builds $w_1 n_{it}$, $w_2 n_{it}$ and $t$; `cost._ajusta` (`scipy.optimize.nnls`, or `lstsq` clipped to ≥ 0).
3. Out-of-sample validation per system (`_clave_sistema`: natoms, ecutwfc, nk, calculation, nspin) if there are ≥ 4 systems and ≥ 8 remaining calculations: bias and dispersion of $\ln(\mathrm{pred}/\mathrm{real})$; otherwise the fit residual.
4. `cost.report_modelo` prints $t_0$, $C_1$, $C_2$, iterations, accuracy and warnings. In a sweep, `cli._run_or_explain` calls `cost.estimar_barrido` (reads each `pw.in` with `descriptores_de_input`; reuses the real $n_k$ of a point already run or from the history with the same formula and mesh) and `cost.report`.

**Where each datum comes from.**

| Datum | Source | Detail |
|---|---|---|
| Timing history | `olla-dft.db` (SQLite, table `calculos`, `wall_s`) | from `olla-dft db` |
| $N_{\mathrm{PW}}$ | volume and `ecutwfc` from the `pw.in` | `cost.n_ondas_planas`; verified against QE in the tests |
| Irreducible $n_k$ | spglib or `pw.out` (`number of k points`) | `cost.k_irreducibles`, `nk_de_salida` |
| Fit weight | `cost.EXP_PESO` | 0.5 (chosen over 63 real calculations) |

**Limits and pitfalls.** "This tool tells ten minutes from six hours … It is not a stopwatch". Without history: "Nothing to calibrate with: the calculation database is empty or stores no times". Poorly varied history (`extrapola_bien` requires ≥ 8 calculations and a range ≥ 5): "predicting a system of another size can be off by a factor of two". spglib and `pw.x` do not always see the same symmetry: "that is where a factor of two or three goes". It does not model MPI parallelism (the time with `-j N` is simply total/N), nor the cost of `ph.x`, nor memory.

**References.** M. C. Payne et al., *Rev. Mod. Phys.* 64, 1045 (1992), DOI 10.1103/RevModPhys.64.1045 (scaling of plane-wave methods). C. L. Lawson and R. J. Hanson, *Solving Least Squares Problems*, SIAM (1995) (NNLS).
