## Electronic structure

This chapter documents the physics that Olla-DFT actually implements — not what the Quantum ESPRESSO manual promises — in the commands that prepare an electronic-structure calculation (`gen`, `kpath`, `info`, `prim`, `conv`, `supercell`) and in those that read its results (`bands`, `gap`, `dos`, `plot`, `effmass`, `fermi`, `unfold`, `wannier`, `topology`, `berry`, `hubbard`, `align`). Each section states which question the command answers, which formulas the code uses (citing the Python function that contains them), which QE file or physical constant every number comes from, and in which cases the result is not valid. Constants and defaults were read from the source code of version 0.35.0 (`qekit/config.py`, `qekit/core/qeout.py`, and each module).

---

### `olla-dft gen` — generate the pw.x and post-processing inputs

**What it answers.** It translates a structure (CIF, POSCAR, XYZ with cell, pw.x input or output) into a consistent set of input files for pw.x, dos.x, projwfc.x and bands.x, choosing for you the cutoffs, the k-point mesh, the high-symmetry path, the number of bands and the treatment of occupations.

**Background for non-experts.** A plane-wave DFT calculation needs four numerical decisions before it starts: (1) how many plane waves to use to describe the orbitals (the *cutoff* `ecutwfc`, a maximum kinetic energy in Rydberg: the higher, the finer the description and the more expensive the run), (2) how many for the density (`ecutrho`, which for ultrasoft or PAW pseudopotentials must be considerably larger), (3) how many k-points to sample in the Brillouin zone (the "resolution" with which one integrates over the infinite crystal), and (4) how to distribute the electrons among bands near the Fermi level (in an insulator the occupation is fixed; in a metal it is "smeared" so that the sum over k converges). Olla-DFT takes (1) and (2) from the header of the pseudopotential itself, (3) from a spacing in reciprocal space (the same idea as VASP's `KSPACING`) and (4) from whether the user declares the system to be an insulator.

A pseudopotential is a "smoothed version" of an atom: it replaces the nucleus and the core electrons by an effective potential, and only the valence electrons are solved explicitly. Every UPF file carries recommended cutoffs and the number of valence electrons, and Olla-DFT reads them. For the band structure one also needs a *high-symmetry path*: a route of straight segments between special points of the Brillouin zone (Γ, X, L…). That path is decided by the seekpath library with the convention of Hinuma et al., and it refers to a standardized primitive cell, which is why `gen` uses that cell whenever the preset includes bands.

**Formulas.**

k-point mesh from a spacing (`qekit/core/kpoints.py: kgrid_from_spacing`):

$$
n_i = \max\!\left(1,\ \left\lceil \frac{|\mathbf{b}_i|}{\Delta k} \right\rceil\right), \qquad \mathbf{b}_i = 2\pi\,(\mathbf{A}^{-1})^{\mathsf T}_{i}
$$

- $n_i$: number of mesh points along reciprocal vector $i$ (dimensionless).
- $\mathbf{b}_i$: reciprocal-lattice vector **including the factor $2\pi$**, in Å⁻¹; $\mathbf{A}$ is the cell matrix with the vectors as rows (Å).
- $\Delta k$: requested spacing in Å⁻¹. Defaults: `kspacing = 0.20` (scf) and `kspacing_nscf = 0.12` (nscf/DOS), read from `qekit/config.py: DEFAULTS`. The `--klevel` levels are `coarse 0.30`, `medium 0.20`, `fine 0.15`, `very-fine 0.10` and `gamma` (Γ only).

If along some axis the widest gap between atoms exceeds `VACIO_MINIMO = 8.0` Å (`kpoints.direcciones_con_vacio`), that $n_i$ is forced to 1. The mesh is written without shift (`0 0 0`, Γ-centred); if it is $1\times1\times1$ the card is `K_POINTS gamma`.

Cell thickness and vacuum gap (`qekit/modules/inputgen.py: espesor_celda`, `hueco_vacio`):

$$
d_i = \frac{V}{|\mathbf{a}_j \times \mathbf{a}_k|}, \qquad h_{\text{Å}} = d_i \cdot \max_m \left(f^{(i)}_{m+1} - f^{(i)}_m\right)
$$

- $d_i$: height of the cell along the normal to the plane of the other two vectors (Å); $V$ is the volume (Å³).
- $h_{\text{Å}}$: the largest gap between sorted fractional coordinates $f^{(i)}$ along axis $i$ (including the gap that wraps around the periodic boundary), converted to Å.

Recommended cutoffs (`qekit/core/pseudo.py: recommend_cutoffs`):

$$
E_{\text{wfc}} = \max_s E^{\text{UPF}}_{\text{wfc},s}, \qquad
E_{\rho} = \max\!\left(\max_s E^{\text{UPF}}_{\rho,s},\ 4\,E_{\text{wfc}}\right)
$$

- $E^{\text{UPF}}_{\text{wfc},s}$, $E^{\text{UPF}}_{\rho,s}$: suggested cutoffs in the UPF header of species $s$ (Ry), read by `pseudo.suggested_cutoffs` (attributes `wfc_cutoff`/`rho_cutoff` in UPF v2, or the text "Suggested minimum cutoff for wavefunctions/charge density" in UPF v1). Values $\le 1$ are ignored.
- If no UPF declares a cutoff, `ecutwfc = 60.0` Ry and `ecutrho = dual × ecutwfc` with `dual = 8` (config) are used. The floor $4E_{\text{wfc}}$ is the physical minimum for plane waves.

Estimated number of bands for nscf/bands (`inputgen._estimate_nbnd`):

$$
n_{\text{bnd}} = \left\lceil 1.25\cdot\frac{N_{\text{el}}}{2} + 4 \right\rceil, \qquad N_{\text{el}} = \sum_{\text{atoms}} Z^{\text{UPF}}_{\text{val}}
$$

With `--nspin 2` it is enlarged to $\lfloor 1.2\,n_{\text{bnd}}\rfloor + 2$. If any UPF does not declare `z_valence`, `nbnd` is not written and pw.x uses its own default.

MD time step (`inputgen.build_pw_input`): $\mathrm{dt}_{\text{Ry}} = \mathrm{dt}_{\text{fs}} / 0.048378$, because pw.x asks for `dt` in Rydberg atomic units (`_FS_POR_UA = 4.8378e-2` fs).

Dipole correction (`inputgen._region_vacio`): the maximum of the sawtooth is placed at the centre of the vacuum gap, `emaxpos = centre`, and its decrease occupies `eopreg = clip(gap/3, 0.02, 0.20)` (fractions of the axis). $h_{\text{Å}} \ge 5$ Å is required.

Estimated cost of a hybrid (`inputgen.generate`), measured on 2-atom silicon: $\text{factor} \approx 3 + 2.6\,n_q$, with $n_q = n_{q1}n_{q2}n_{q3}$ the exact-exchange mesh. This is an empirical rule, not a formula.

**How Olla-DFT computes it.**

1. `qekit/cli.py: _cmd_gen` reads the structure with `qekit/core/structure.py: load` (ASE; for `POSCAR/CONTCAR` it forces `format="vasp"`, and if the file carries several images it keeps the last one). It combines `--klevel`/`--kspacing`, `--mag` (which switches on `nspin=2`), `--hubbard`, `--soc`, `--functional`, `--exx-grid` and the MD options into a `GenOptions`.
2. `qekit/modules/inputgen.py: generate` decides the working cell: if the preset is `bands` or `all`, it calls `kpoints.get_kpath` (seekpath) and uses the **standardized primitive cell** it returns; with `--primitive` it uses `structure.primitive` (spglib); otherwise the cell as given.
3. `qekit/core/pseudo.py: resolve` looks for one UPF per element in `pseudo_dir` (a file whose name starts with the symbol followed by a non-alphabetic character, extension `.upf` case-insensitive). It honours `--pseudo El=file`, and `_coherencia_de_funcional` re-selects whatever is needed so that all pseudopotentials share the same functional (preference `PBE > PBESOL > REVPBE > PZ > BLYP`).
4. `pseudo.recommend_cutoffs` sets `ecutwfc`/`ecutrho`; the user can override them with `--ecutwfc`/`--ecutrho`.
5. `kpoints.kgrid_from_spacing` produces the scf and nscf meshes; `_estimate_nbnd` the number of bands.
6. `inputgen.build_pw_input` writes `&CONTROL` (with `tprnfor`, `tstress`, `outdir='./out'`), `&SYSTEM` (`ibrav=0`, cutoffs, occupations, spin, SOC, Hubbard, `tot_charge`, dipole, hybrid, `nosym`/`noinv`), `&ELECTRONS` (`conv_thr=1e-8`, `mixing_beta=0.4`, `electron_maxstep=200`), `&IONS`/`&CELL` according to the preset (BFGS for relax, Verlet for MD with the requested thermostat, `press_conv_thr=0.05`), and the cards `ATOMIC_SPECIES` (masses from `ase.data.atomic_masses`), `ATOMIC_POSITIONS crystal`, `CELL_PARAMETERS angstrom`, the `HUBBARD` card if `--hubbard-style card`, and `K_POINTS`.
7. For bands, `kpoints.kpath_card` writes `K_POINTS crystal_b` with `band_points` points per segment (20 by default) and `KPATH.txt` with the labels; `build_bandsx_input` writes `bands_pp.in` (`lsym=.true.`). For DOS, `build_dos_input` and `build_projwfc_input` write `dos.in` and `projwfc.in` with `DeltaE = 0.02` eV.
8. `build_run_script` and `build_run_python_script` write `run.sh` (with `set -e -o pipefail` and `mpirun -np $NP`) and `run.py` in the order `pw.x → dos.x/projwfc.x → pw.x (bands) → bands.x`.

**Where each datum comes from.**

| Datum | Source | Detail |
|---|---|---|
| Cell and positions | user file (CIF/POSCAR/…) | `structure.load` via `ase.io.read` |
| Primitive cell and k-path | seekpath library | `kpoints.get_kpath` with `symprec = 1e-4` Å (`structure.SYMPREC`) |
| Suggested cutoffs | UPF header (`wfc_cutoff`, `rho_cutoff`) | `pseudo.suggested_cutoffs`; reads only the first 20 000 characters |
| Valence electrons | `z_valence` in the UPF | `pseudo.z_valence`, used in `_estimate_nbnd` |
| Pseudopotential type and relativity | `pseudo_type`, `relativistic` in the UPF | `pseudo.pseudo_type`, `pseudo.relativistic` |
| `ecutwfc`, `dual`, `kspacing`, `degauss`, `smearing`, `nproc` | `~/.config/qekit/config.ini` or `config.DEFAULTS` | 60 Ry, 8, 0.20 Å⁻¹, 0.01 Ry, `cold`, 4 |
| Atomic masses | `ase.data.atomic_masses` | `ATOMIC_SPECIES` card |
| Hybrid parameters | table `inputgen.HIBRIDOS` | HSE: `exx_fraction 0.25`, `screening_parameter 0.106` bohr⁻¹; PBE0 0.25; B3LYP 0.20; Gau-PBE 0.24 |
| Hubbard orbital (card) | atomic number (`inputgen._orbital_hubbard`) | 3d (Z 21–30), 4d (39–48), 5d (72–80), 4f (57–71), 5f (89–103), `2p` otherwise |
| fs → Rydberg a.u. conversion | constant `inputgen._FS_POR_UA` | 4.8378e-2 fs |

**Limits and pitfalls.**

- The "automatic" cutoffs are those **suggested by the UPF**, not a convergence: the report says "(automático)". If the UPF does not declare them, the code falls back to 60 Ry / 480 Ry with no warning beyond the report.
- `--soc` writes `noncolin=.true.` and `lspinorb=.true.` only if every pseudopotential declares `relativistic="full"`: `inputgen.generate` calls `sweep.check_soc_pseudos` before writing anything and, if any of them is scalar-relativistic or does not declare it, it stops with "el acoplamiento espín-órbita necesita pseudopotenciales TOTALMENTE RELATIVISTAS (relativistic='full'), y estos no lo son". The reason, quoted in the error itself: with scalar pseudopotentials lspinorb "devuelve un desdoblamiento espín-órbita de cero que parece un resultado válido". `--soc` and `--nspin 2` are rejected together.
- `--hubbard-style legacy` (the default) writes `lda_plus_u` and `Hubbard_U(i)`, a syntax removed in QE ≥ 7.1; for those versions `--hubbard-style card` is needed, which writes the `HUBBARD (ortho-atomic)` card. The orbital in the card is deduced only from the atomic number.
- With hybrids, `nqx` **must divide** the k mesh; otherwise the command stops with "la malla de intercambio exacto tiene que DIVIDIR la de k". The report also warns that with `1x1x1` "el resultado va a salir claramente sobrestimado" and that pw.x cannot do a `calculation='bands'` with EXX.
- Without `--mag`, `--nspin 2` starts with zero magnetization and the report warns: "sin magnetización inicial el cálculo suele converger a la solución no magnética".
- The dipole correction requires a vacuum gap ≥ 5 Å; otherwise the error reads "la corrección dipolar necesita vacío en la dirección …".
- The scf mesh of a `bands` preset is computed on the seekpath primitive cell, which may not be the one the user supplied; the report says so ("AVISO: se usó la celda primitiva estandarizada").
- For MD, `nosym` is forced, and a warning is issued for fewer than 20 atoms or `dt > 2` fs.
- `tot_charge` is compensated by a uniform background; the report reminds that the energy of a charged cell is not comparable to the neutral one.
- The mesh is always uniform and Γ-centred (no shift), as the `kpoints.py` docstring states: it is not a shifted Monkhorst-Pack mesh, and with even $n$ it contains Γ where an MP mesh would not.

**References.**

- Y. Hinuma, G. Pizzi, Y. Kumagai, F. Oba, I. Tanaka, *Comput. Mater. Sci.* **128**, 140 (2017) — seekpath k-path convention. DOI 10.1016/j.commatsci.2016.10.015.
- A. Togo, I. Tanaka, "Spglib: a software library for crystal symmetry search", arXiv:1808.01590 (2018).
- N. Marzari, D. Vanderbilt, A. De Vita, M. C. Payne, *Phys. Rev. Lett.* **82**, 3296 (1999) — "cold" smearing. DOI 10.1103/PhysRevLett.82.3296.
- J. Heyd, G. E. Scuseria, M. Ernzerhof, *J. Chem. Phys.* **118**, 8207 (2003) — HSE.
- L. Bengtsson, *Phys. Rev. B* **59**, 12301 (1999) — dipole correction for slabs.
- G. Prandini et al., *npj Comput. Mater.* **4**, 72 (2018) — the SSSP pseudopotential library and cutoffs.

---

### `olla-dft kpath` — high-symmetry path

**What it answers.** Which standard path through the Brillouin zone (and in which cell it is expressed) must be used to draw the band structure of this structure.

**Background for non-experts.** The Brillouin zone is the "unit cell" of wave-vector space; bands $E(\mathbf k)$ are drawn along a route through its highest-symmetry points, because that is where bands touch, cross or have their extrema. Which points and in which order depends on the space group and the cell shape, and there are several incompatible conventions in the literature. Olla-DFT delegates the choice to seekpath (the convention of Hinuma et al., the same used by Materials Cloud), which also **standardizes the cell**: the coordinates of the special points are only valid in that primitive cell, not necessarily in the one in the user's CIF.

**Formulas.** No formulas of its own: the command calls `seekpath.get_path` and transcribes its result. The only arithmetic is the "same cell" criterion (`qekit/core/kpoints.py: get_kpath`):

$$
\text{cell\_changed} = \neg\left[N_{\text{prim}} = N_{\text{in}} \ \wedge\ \max_{ij} |A^{\text{prim}}_{ij} - A^{\text{in}}_{ij}| \le 10^{-5}\ \text{Å}\right]
$$

- $A^{\text{prim}}$, $A^{\text{in}}$: cell matrices (Å) of the seekpath primitive and of the input; $N$: number of atoms.

**How Olla-DFT computes it.**

1. `qekit/cli.py: _cmd_kpath` loads the structure with `structure.load`.
2. `kpoints.get_kpath` converts to the spglib tuple (`structure.to_spglib_cell`), calls `seekpath.get_path(..., symprec=1e-4)` and rebuilds the primitive cell (`structure.from_spglib_cell`) from `primitive_lattice`, `primitive_positions`, `primitive_types`.
3. `kpoints.kpath_text` prints the space group (`spacegroup_international`, `spacegroup_number`), the compacted path (`Γ — X — U | K — Γ …`), the fractional coordinates of every point (`point_coords`) with the "pretty" labels of `pretty_label` (GAMMA→Γ, DELTA_0→Δ0), and a warning if the cell changed.

**Where each datum comes from.**

| Datum | Source | Detail |
|---|---|---|
| Path and coordinates | seekpath library (`get_path`) | keys `path`, `point_coords` |
| Space group | seekpath (spglib underneath) | `spacegroup_international`, `spacegroup_number` |
| Symmetry tolerance | constant `structure.SYMPREC` | 1e-4 Å |
| Primitive cell | seekpath | `primitive_lattice/positions/types` |

**Limits and pitfalls.**

- The coordinates are **in the standardized primitive cell**. If `cell_changed` is true, the text warns: "el k-path está referido a la celda primitiva estandarizada, que difiere de la celda de entrada. Usa esa celda primitiva en el cálculo de bandas". Using those coordinates in the original cell gives a wrong path with no error whatsoever.
- With `symprec = 1e-4` Å, a relaxed structure with numerical noise can lose symmetry and receive the path of a lower space group; there is no command-line option to change the tolerance.
- seekpath ignores magnetism and spin-orbit coupling when choosing the space group.

**References.**

- Y. Hinuma, G. Pizzi, Y. Kumagai, F. Oba, I. Tanaka, *Comput. Mater. Sci.* **128**, 140 (2017). DOI 10.1016/j.commatsci.2016.10.015.
- W. Setyawan, S. Curtarolo, *Comput. Mater. Sci.* **49**, 299 (2010) — the alternative convention that seekpath does **not** use.

---

### `olla-dft info` — structure and symmetry

**What it answers.** What the structure file contains: formula, lattice parameters, volume, space group, point group, Wyckoff positions and how many atoms the primitive cell would have.

**Background for non-experts.** Before computing anything it is worth knowing whether the cell is the smallest possible one (the primitive) or a larger one (conventional or supercell), because the cost of the calculation grows with the number of atoms, and whether the structure has the symmetry one believes. spglib compares every atom with the candidate symmetry operations within a tolerance and returns the space group in international notation (for example `Fd-3m`, No. 227 for silicon), the Hall symbol, the point group and the Wyckoff letter of each site (a label for the kind of symmetry position each atom occupies).

**Formulas.** No formulas beyond the cell geometry that ASE computes (`atoms.cell.cellpar()` returns $a, b, c, \alpha, \beta, \gamma$ and `atoms.cell.volume` the volume $V = |\det \mathbf A|$ in Å³).

**How Olla-DFT computes it.**

1. `qekit/cli.py: _cmd_info` → `structure.load`.
2. `qekit/core/structure.py: info_text` calls `symmetry_dataset` (`spglib.get_symmetry_dataset` with `symprec = 1e-4`) and `primitive` (`spglib.standardize_cell(to_primitive=True)`) to count the atoms of the primitive cell.
3. It prints formula, composition, number of atoms, volume, lattice parameters, space group (`international`, `number`), Hall symbol, point group, atoms in the primitive cell, Wyckoff positions (sorted set of `ds.wyckoffs`) and the cell vectors.

**Where each datum comes from.**

| Datum | Source | Detail |
|---|---|---|
| Space group, Hall, point group, Wyckoff | spglib library | `structure.symmetry_dataset` |
| Lattice parameters and volume | ASE (`Cell.cellpar`, `Cell.volume`) | a, b, c in Å; angles in degrees |
| Atoms in the primitive cell | spglib `standardize_cell` | `structure.primitive` |
| Tolerance | `structure.SYMPREC` | 1e-4 Å |

**Limits and pitfalls.**

- If spglib cannot determine the symmetry, the command fails with `RuntimeError("spglib no pudo determinar la simetría de la estructura")`.
- The tolerance is fixed (1e-4 Å); there is no `--symprec`.
- Only the distinct Wyckoff letters are listed, not how many atoms sit in each.

**References.**

- A. Togo, I. Tanaka, "Spglib: a software library for crystal symmetry search", arXiv:1808.01590 (2018).
- International Tables for Crystallography, Vol. A (IUCr) — notation of space groups and Wyckoff positions.

---

### `olla-dft prim` — standardized primitive cell

**What it answers.** Which is the smallest cell that reproduces the crystal by translation, written in spglib's standard orientation.

**Background for non-experts.** Many CIFs come in the conventional cell (the one that makes the symmetry visible, e.g. the 8-atom cube of silicon), but the calculation only needs the primitive cell (2 atoms in silicon). Reducing it saves a factor equal to the ratio of atom counts in cost, without changing the physics. "Standardized" means that spglib reorients it and expresses it with the choice of vectors fixed by the International Tables convention, so that two equivalent inputs give the same output.

**Formulas.** No arithmetic of its own: it is `spglib.standardize_cell(cell, to_primitive=True, symprec=1e-4)`.

**How Olla-DFT computes it.**

1. `qekit/cli.py: _cmd_prim` → `structure.load`.
2. `structure.primitive` → spglib → `structure.from_spglib_cell` (an `Atoms` with `pbc=True`).
3. `structure.convert` writes the result according to the extension of `-o` (default `primitive.cif`): CIF, POSCAR/`.vasp` (with `direct=True, sort=True`) or any format ASE can infer.

**Where each datum comes from.**

| Datum | Source | Detail |
|---|---|---|
| Primitive cell | spglib `standardize_cell(to_primitive=True)` | `structure.primitive` |
| Tolerance | `structure.SYMPREC` | 1e-4 Å |
| Output format | file extension | `structure.convert` |

**Limits and pitfalls.**

- spglib's primitive cell is **not** necessarily the same as the seekpath one used by `gen -p bands` (seekpath applies its own additional standardization); for bands, let `gen` choose.
- When writing POSCAR the atoms are **reordered** by species (`sort=True`); if the user relied on a specific order (for example for a `--displace`), it is lost.
- If spglib fails: `RuntimeError("spglib no pudo estandarizar la celda")`.

**References.**

- A. Togo, I. Tanaka, arXiv:1808.01590 (2018).

---

### `olla-dft conv` — standardized conventional cell

**What it answers.** Which is the conventional cell (the one showing the full symmetry of the crystal system) of the structure.

**Background for non-experts.** It is the inverse operation to `prim`: start from any cell and obtain the "textbook" cell (cubic for silicon, hexagonal for graphite), useful for building surfaces, supercells or for comparison with diffraction data, even though it has more atoms than strictly needed for the calculation.

**Formulas.** No arithmetic of its own: `spglib.standardize_cell(cell, to_primitive=False, symprec=1e-4)`.

**How Olla-DFT computes it.**

1. `qekit/cli.py: _cmd_conv` → `structure.load`.
2. `structure.conventional` → spglib → `from_spglib_cell`.
3. `structure.convert` writes `-o` (default `conventional.cif`).

**Where each datum comes from.**

| Datum | Source | Detail |
|---|---|---|
| Conventional cell | spglib `standardize_cell(to_primitive=False)` | `structure.conventional` |
| Tolerance | `structure.SYMPREC` | 1e-4 Å |

**Limits and pitfalls.**

- Same as `prim`: fixed tolerance, atom reordering in POSCAR, error if spglib cannot standardize.
- For a low-symmetry (triclinic) cell the "conventional" cell coincides with the primitive one and the command changes nothing.

**References.**

- A. Togo, I. Tanaka, arXiv:1808.01590 (2018).

---

### `olla-dft supercell` — build a supercell

**What it answers.** The structure repeated $n_x \times n_y \times n_z$ times along its three cell vectors.

**Background for non-experts.** A supercell is several cells glued together and treated as a single periodic unit. It is needed to place a defect, a dopant or an adsorbed molecule at low concentration, for molecular dynamics, or to compute phonons by finite displacements. The price is that the Brillouin zone shrinks by the same factor and the bands "fold" (see `unfold`).

**Formulas.** ASE's `atoms.repeat((nx, ny, nz))`: the new cell is $\mathbf a'_i = n_i \mathbf a_i$ and every atom is copied into the $n_x n_y n_z$ translations $\sum_i m_i \mathbf a_i$ with $0 \le m_i < n_i$.

**How Olla-DFT computes it.**

1. `qekit/cli.py: _cmd_supercell` → `structure.load`.
2. `structure.supercell` checks that the three factors are ≥ 1 (otherwise `ErrorDeUso("los factores de la supercelda deben ser >= 1")`) and calls `Atoms.repeat`.
3. `structure.convert` writes `-o` (default `supercell.cif`).

**Where each datum comes from.**

| Datum | Source | Detail |
|---|---|---|
| Factors $n_x, n_y, n_z$ | user parameters (positional) | integers ≥ 1 |
| Repetition | ASE `Atoms.repeat` | diagonal multiples only |

**Limits and pitfalls.**

- Only **diagonal** supercells (multiples of each vector); general matrices $\mathbf A' = \mathbf M \mathbf a$ such as those `unfold` can recognize cannot be built.
- No symmetry reduction and no check that the supercell is "reasonable" (e.g. cubic).

**References.**

- A. H. Larsen et al., "The atomic simulation environment—a Python library for working with atoms", *J. Phys.: Condens. Matter* **29**, 273002 (2017). DOI 10.1088/1361-648X/aa680e.

---
### `olla-dft bands` — band structure and band gap

**What it answers.** How the energy of every electronic state varies along the high-symmetry path, whether the system has a gap or is metallic, where the valence-band maximum (VBM) and the conduction-band minimum (CBM) are, whether the gap is direct or indirect, and — with `--fat` — which atomic orbital each band "is made of".

**Background for non-experts.** In a crystal electrons do not have discrete energies but *bands*: for every wave vector $\mathbf k$ (a "direction and wavelength" of the electron) there is a list of allowed energies $\varepsilon_n(\mathbf k)$. Drawing them along a path in the Brillouin zone gives the classic "spaghetti" plot. The *band gap* is the distance between the highest occupied band and the lowest empty one. If the maximum of one and the minimum of the other are at the same $\mathbf k$ the gap is *direct* (the material absorbs and emits light efficiently); otherwise it is *indirect*. If some band crosses the Fermi level (the energy up to which states are filled), it is a metal and there is no gap.

*Fat bands* answer a different question: what weight each atomic orbital (nickel $d$, oxygen $p$) has in each state. projwfc.x projects every wavefunction onto atomic orbitals and writes the weights; Olla-DFT draws them as dots whose size is proportional to the weight on top of the bands.

**Formulas.**

Unit conversion of the pw.x XML (`qekit/core/qeout.py: read_xml`):

$$
E_{\text{eV}} = E_{\text{Ha}} \cdot 27.211386245988, \qquad
\mathbf{k}_{\text{Å}^{-1}} = \mathbf{k}_{2\pi/a} \cdot \frac{2\pi}{a_{\text{bohr}}\cdot 0.529177210903}, \qquad
\mathbf{k}_{\text{frac}} = \mathbf{k}_{\text{Å}^{-1}}\, \mathbf{B}^{-1}
$$

- $a_{\text{bohr}}$: `alat` from the XML (bohr). $\mathbf B$: matrix with the reciprocal vectors $\mathbf b_i = 2\pi(\mathbf A^{-1})^{\mathsf T}_i$ as rows (Å⁻¹).

Cumulative distance on the x axis (`qekit/modules/bands.py: _build_kdist`):

$$
x_0 = 0, \qquad x_i = x_{i-1} + \begin{cases} 0 & i \in \text{breaks} \\ |\mathbf{k}_i - \mathbf{k}_{i-1}| & \text{otherwise} \end{cases}
$$

- `breaks`: indices where two special points appear consecutively (a `U|K` discontinuity of the path), detected by `_detect_breaks`.

Metal / insulator classification (`bands.analyze_gap`), with `CROSS_TOL = 1e-6` eV and reference $E_{\text{ref}}$:

$$
\text{crosses}_n = \left[\min_{\mathbf k}\varepsilon_n < E_{\text{ref}} - \delta\right] \wedge \left[\max_{\mathbf k}\varepsilon_n > E_{\text{ref}} + \delta\right]
$$

If any band crosses, it is a metal. Otherwise $n_v = \max\{n : \max_{\mathbf k}\varepsilon_n \le E_{\text{ref}}+\delta\}$ and $n_c = \min\{n : \min_{\mathbf k}\varepsilon_n > E_{\text{ref}}-\delta\}$ (forcing $n_c \ge n_v+1$), and

$$
E_{\text{VBM}} = \max_{\mathbf k}\varepsilon_{n_v}(\mathbf k), \quad
E_{\text{CBM}} = \min_{\mathbf k}\varepsilon_{n_c}(\mathbf k), \quad
E_g = E_{\text{CBM}} - E_{\text{VBM}}, \quad
E_g^{\text{dir}} = \min_{\mathbf k}\left[\varepsilon_{n_c}(\mathbf k) - \varepsilon_{n_v}(\mathbf k)\right]
$$

The gap is direct if $\arg\max\varepsilon_{n_v} = \arg\min\varepsilon_{n_c}$ (same k-point index, not same coordinate).

Reference $E_{\text{ref}}$ (in this order): `<fermi_energy>` from the XML; if absent, `<highestOccupiedLevel>`; if also absent and `nspin = 1`, the midpoint $\tfrac12[\max_{\mathbf k}\varepsilon_{n_{occ}-1} + \min_{\mathbf k}\varepsilon_{n_{occ}}]$ with $n_{occ} = \mathrm{round}(N_{\text{el}}/2)$; as a last resort, the median of all energies.

Energy zero in the plot and in the exported data (`bands.reference_energy`): `--ref auto` uses the VBM when there is a gap and $E_F$ for a metal; `fermi`, `vbm`, `none` as their names say.

Weight of a selector in the fat bands (`bands.peso_de`):

$$
w_{n}(\mathbf k) = \sum_{i \in \text{selector}} |\langle \phi_i | \psi_{n\mathbf k}\rangle|^2
$$

- $|\langle \phi_i | \psi_{n\mathbf k}\rangle|^2$: the coefficients `psi = 0.498*[# 1] + …` from the text output of projwfc.x for atomic state $i$, read by `bands.leer_proyecciones`. **Not normalized**: the part missing up to 1 is the fraction of the wavefunction that falls in no atomic sphere, which `report_fat` quantifies as $1 - \langle\sum_i w_i\rangle$.

**How Olla-DFT computes it.**

1. `qekit/cli.py: _cmd_bands` → `bands.load(path, prefix)`.
2. `qeout.find_xml` locates the XML (`./out/*.xml`, `./*.xml` or `*.save/data-file-schema.xml`, checking that it contains "espresso"); `qeout.read_xml` reads `<atomic_structure>` (cell in bohr), `<band_structure>` (`nbnd`, `nelec`, `lsda`, `noncolin`, `fermi_energy`, `highestOccupiedLevel`, `lowestUnoccupiedLevel`), and every `<ks_energies>` (`k_point` with `weight`, `eigenvalues`, `occupations`). With `lsda`, the list of eigenvalues per k is split into two halves (up/down).
3. Labels come from `KPATH.txt` (`qeout.read_kpath_labels`) or, if it does not exist, from the `K_POINTS crystal_b` card of `bands.in` with `! G` comments (`qeout.read_crystal_b_card`). `qeout.match_labels_to_kpoints` assigns them to indices with tolerance `1e-3` in fractional coordinates, always moving forward and tolerating reciprocal-lattice translations.
4. `bands.analyze_gap` per spin channel; `bands.gap_report` prints the summary and the reminder that "los funcionales GGA/LDA subestiman el gap sistemáticamente (típicamente 30–50 %)".
5. With `--fat`: `bands.leer_proyecciones` reads `projwfc.out` (or `proj.out`, `projwfc_bands.out`) from the same bands calculation; `comprobar_compatibilidad` requires the same number of k-points; `peso_de` sums the states that match the selector (`Ni`, `Ni-d`, `d`, `atomo:3`).
6. `bands.export` writes `BAND.dat` (or `BAND_up.dat`/`BAND_dw.dat`), `KLABELS.dat` and `BAND_GAP.txt`; `bands.plot` draws with matplotlib (bands in ink, spin ↓ dashed, VBM as a circle and CBM as a square, fat bands as `scatter` with `s = w · fat_scale`, `fat_scale = 55`).

**Where each datum comes from.**

| Datum | Source | Detail |
|---|---|---|
| Eigenvalues, k-points, weights | `prefix.xml` from pw.x (`<ks_energies>`) | `qeout.read_xml`; Ha → eV |
| Fermi energy | `<fermi_energy>` in `prefix.xml` | only if the scf used smearing |
| HOMO / LUMO | `<highestOccupiedLevel>` / `<lowestUnoccupiedLevel>` | fixed occupations |
| Number of electrons, bands, spin | `<nelec>`, `<nbnd>`, `<lsda>`, `<noncolin>` | `nbnd` is recomputed from the length of the eigenvalue list |
| Cell and `alat` | `<atomic_structure alat=…>` and `<cell>` | bohr → Å with 0.529177210903 |
| High-symmetry labels | `KPATH.txt` from `olla-dft gen` or `bands.in` | matching tolerance 1e-3 |
| Orbital weights (fat bands) | `projwfc.out` (`psi = …` blocks) | `bands.leer_proyecciones`; `state #` gives atom, element and $l$ |
| Hartree in eV, bohr in Å | constants `qeout.HARTREE_EV`, `qeout.BOHR_ANG` | 27.211386245988 eV; 0.529177210903 Å (CODATA 2018) |

**Limits and pitfalls.**

- In a `bands` calculation with fixed occupations the XML may not carry `<fermi_energy>`; the reference is then `<highestOccupiedLevel>`, which QE inherits from the scf. If the scf used smearing, $E_F$ may fall in the middle of the gap or inside a flat band; a band touching $E_F \pm 10^{-6}$ eV at a single point is classified as a metal.
- If `nbnd` only covers the occupied bands, the report says "No hay bandas de conducción en el cálculo (aumenta nbnd para obtener el gap)".
- With `nspin = 2` each channel is analysed separately with the **same** reference; the report does not compute the global gap between different channels (e.g. VBM up and CBM down).
- The plot and `--ref auto` always use the analysis of channel 0 (spin up) to decide the zero and mark the extrema.
- "Direct" is decided by comparing k-point **indices**; two symmetry-equivalent k-points at different indices count as indirect.
- Fat bands: if projwfc.x was run on the DOS nscf rather than on the bands, `comprobar_compatibilidad` stops: "las bandas tienen N puntos k y las proyecciones M. No son del mismo cálculo". If more than 10 % of the mean weight falls outside atomic spheres, `report_fat` warns: "De media, un X % de cada función de onda NO cae dentro de ninguna esfera atómica".
- Weights of states with $l>3$ are labelled `l4`, `l5`… and cannot be selected by letter.
- SOC (`noncolin`) calculations are read as a single channel; projwfc weights with $j$ (`p_j1.5`) are grouped only by the orbital letter.

**References.**

- P. Giannozzi et al., *J. Phys.: Condens. Matter* **29**, 465901 (2017) — Quantum ESPRESSO (XML format, projwfc.x). DOI 10.1088/1361-648X/aa8f79.
- J. P. Perdew, M. Levy, *Phys. Rev. Lett.* **51**, 1884 (1983) and L. J. Sham, M. Schlüter, *Phys. Rev. Lett.* **51**, 1888 (1983) — why DFT underestimates the gap.
- CODATA 2018, E. Tiesinga et al., *Rev. Mod. Phys.* **93**, 025010 (2021) — constants.

---

### `olla-dft gap` — band-gap report only

**What it answers.** The same as the analysis part of `bands` — metal or not, VBM, CBM, fundamental gap and minimum direct gap per spin channel — without exporting data or plotting.

**Background for non-experts.** It is the most frequent question asked of a band calculation ("how big is the gap?") decoupled from the figure. It works equally on an `scf`, an `nscf` mesh or a `bands` path: it reads any pw.x XML with eigenvalues. On a mesh, the gap obtained is that of the sampled points, which can be larger than the true one if the extremum falls between points.

**Formulas.** Exactly those of `bands.analyze_gap` described under `olla-dft bands` (classification with `CROSS_TOL = 1e-6` eV, $E_g = E_{\text{CBM}} - E_{\text{VBM}}$, $E_g^{\text{dir}} = \min_{\mathbf k}[\varepsilon_{n_c} - \varepsilon_{n_v}]$).

**How Olla-DFT computes it.**

1. `qekit/cli.py: _cmd_gap` → `bands.load(path, prefix)` (reads the XML and, if present, `KPATH.txt` or `bands.in` to label the points).
2. `bands.gap_report` loops over `range(nspin)` calling `analyze_gap`, prints the XML path, calculation type, `nbnd`, `nk`, `nelec`, $E_F$ if present, and the result per channel.

**Where each datum comes from.**

| Datum | Source | Detail |
|---|---|---|
| Eigenvalues and k | `prefix.xml` from pw.x | `qeout.read_xml` |
| Reference | `<fermi_energy>` → `<highestOccupiedLevel>` → electron count → median | `bands.analyze_gap` |
| Label of the k-point of the extremum | `KPATH.txt` / `bands.in`, or the fractional coordinates | `bands._label_for` |

**Limits and pitfalls.**

- On a symmetry-reduced scf/nscf mesh, "direct" can only be detected if VBM and CBM fall at the **same index** of the irreducible-point list.
- An XML without `<output>` (unfinished calculation) yields `FaltanDatos("… no contiene una sección <output>")`.
- `gap_report` calls `analyze_gap` twice per channel (once for the report and once to decide whether to print the GGA reminder); this is only cost, it does not change the result.

**References.**

- The same as for `olla-dft bands`.

---

### `olla-dft dos` — total and projected density of states

**What it answers.** How many electronic states there are per unit energy (DOS), how they are distributed among elements and orbitals (PDOS), how large the DOS is at the Fermi level, and — with `--dband` — the centre, width and filling of a projected band (the "d-band centre" of catalysis).

**Background for non-experts.** The DOS is the histogram of the energies of all states in the Brillouin zone: where it is high there are many states, where it is zero there is a gap. dos.x computes it from the eigenvalues of the nscf (with the tetrahedron method that `gen -p dos` requests); projwfc.x decomposes every state into atomic orbitals and gives a PDOS per atom and orbital. Summing the files by element and by orbital letter ($s, p, d, f$) gives the chemical decomposition that gets published.

The *d-band centre* is the mean energy of the $d$ PDOS of a transition metal relative to the Fermi level. It is an empirical descriptor: the closer to the Fermi level (less negative), the more strongly the surface adsorbs (Hammer–Nørskov model).

**Formulas.**

File columns (`qekit/modules/dos.py: read_dos_file`, `read_pdos_file`): from `<prefix>.dos` one takes $E$, DOS (one or two columns depending on spin) and the integrated DOS; from `pdos_atm#N(El)_wfc#M(l)` one takes the `ldos` column (already summed over $m$), or `ldosup`/`ldosdw` with spin, according to the number of columns $1 + n_s(1 + (2l+1))$.

Aggregated PDOS (`dos.load`): $\rho_{\text{El},l}(E) = \sum_{\text{atoms } a \in \text{El}} \sum_{\text{wfc with } l} \text{ldos}_{a,l}(E)$; if dos.x and projwfc.x use different energy meshes, the projections are interpolated linearly (`np.interp`, zero outside the range) onto the total-DOS mesh. Without `<prefix>.dos`, the total DOS is defined as $\sum_{\text{El},l}\rho_{\text{El},l}$.

DOS at the Fermi level (`dos.report`): $\rho(E_F) = \sum_s \rho_s(E_{i^*})$ with $i^* = \arg\min_i |E_i - E_F|$; it is called "compatible with a gap" when $\rho(E_F) < 10^{-3}$ states/eV.

Moments of a projected band (`dos.momentos`), with $e = E - E_F$ and $\rho$ the PDOS of the selector (per spin channel, integrated with the trapezoidal rule `np.trapezoid`):

$$
N = \int \rho(e)\,de, \qquad
\varepsilon_c = \frac{1}{N}\int e\,\rho(e)\,de, \qquad
W = \sqrt{\frac{1}{N}\int (e-\varepsilon_c)^2\rho(e)\,de}, \qquad
f = \frac{1}{N}\int_{e \le 0}\rho(e)\,de
$$

- $N$: integrated states (dimensionless); $\varepsilon_c$: centre (eV relative to $E_F$); $W$: rms width (eV); $f$: filling (fraction).
- With two channels, the reported value is the average of each quantity weighted by $N_s$; the "exchange splitting" is $\varepsilon_c^{\uparrow} - \varepsilon_c^{\downarrow}$.
- Relative tail: $\max(\rho)$ over the last $\max(3, n/50)$ points divided by the global $\max(\rho)$; if it exceeds 0.05 the band is cut off at the top.

**How Olla-DFT computes it.**

1. `qekit/cli.py: _cmd_dos` → `dos.load(path, prefix)`.
2. `dos.load` tries to read the XML with `qeout.read_xml` to get $E_F$ and the prefix; it looks for `<prefix>.dos` (or `*.dos`) and all `*pdos_atm#*`; it parses the name with the expression `pdos_atm#(\d+)\(([A-Za-z]+)\)_wfc#(\d+)\(([A-Za-z])…\)` to obtain element and orbital letter (also with SOC's `p_j1.5`).
3. If the XML gives no $E_F$, it takes it from the `EFermi = …` comment in the header of the `.dos`.
4. It orders the projections by element (order of appearance) and orbital $s,p,d,f$; `by_element` sums orbitals.
5. `dos.report` prints the energy range, $E_F$, the origin of the zero (`reference_energy`: Fermi unless `--ref none`), channels, projections and $\rho(E_F)$.
6. `dos.export` writes `DOS.dat` (E, DOS[_up/_dw], DOS_integrada) and `PDOS.dat` (per element_orbital and totals per element); `dos.plot` / `dos.draw` draw the total with a fill, the PDOS by `--mode orbital|element|total`, spin ↓ mirrored downwards.
7. With `--dband El[-orb]`: `dos.momentos` requires $E_F$ (otherwise `ErrorDeUso("no se encontró la energía de Fermi…")`) and the key (El, orb); `report_momentos` prints centre, width, filling and warnings.

**Where each datum comes from.**

| Datum | Source | Detail |
|---|---|---|
| Total, integrated DOS | `<prefix>.dos` from dos.x | columns E, dos[, dosup, dosdw], int |
| PDOS per atom/orbital | `<prefix>.pdos.pdos_atm#N(El)_wfc#M(l)` from projwfc.x | column `ldos` (or `ldosup`, `ldosdw`) |
| Fermi energy | `<fermi_energy>` in `prefix.xml`; else `EFermi` in the `.dos` header | `dos.load` |
| Energy mesh | that of dos.x (`DeltaE = 0.02` eV in `gen`) | PDOS are interpolated onto it |
| "Gap" threshold | constant in `dos.report` | 1e-3 states/eV |
| Tail threshold | constant in `dos.momentos`/`report_momentos` | 5 % of the peak |

**Limits and pitfalls.**

- A `pdos_atm#` file with a different number of points from the first one is skipped, and `dos.load` records it in `DOSData.avisos`; the report prints it: "se han SALTADO N archivo(s) de projwfc.x cuya malla de energía no coincide con la del primero (… puntos), así que la PDOS está incompleta… Casi siempre es que hay archivos de dos corridas de projwfc.x mezclados en la misma carpeta". The exported data still lack that orbital: move the old files away and reload.
- The total DOS defined as a sum of PDOS (when `.dos` is missing) omits the part of the wavefunctions that falls outside atomic spheres: it will be smaller than the real one.
- `momentos` integrates the **whole** available range except for `--dband-emax`; if the PDOS has not decayed, it warns: "al final del rango todavía queda un X % del pico de PDOS. La banda está CORTADA por arriba, así que el centro sale más bajo de lo que debería". The text recommends "Vuelve a correr projwfc.x con un Emax mayor".
- The d-band centre "es una correlación empírica dentro de una misma familia de metales, no una ley" (report text).
- With SOC, the orbital letter is taken from `p_j1.5` → `p`; the $j$ components are summed.
- `--ref vbm` does not exist for the DOS: `reference_energy` only distinguishes `none` from the rest (always Fermi).

**References.**

- P. E. Blöchl, O. Jepsen, O. K. Andersen, *Phys. Rev. B* **49**, 16223 (1994) — tetrahedron method (dos.x, `tetrahedra_opt`).
- B. Hammer, J. K. Nørskov, *Surf. Sci.* **343**, 211 (1995); *Adv. Catal.* **45**, 71 (2000) — d-band centre model.
- P. Giannozzi et al., *J. Phys.: Condens. Matter* **29**, 465901 (2017) — projwfc.x.

---

### `olla-dft plot` — combined bands + DOS figure

**What it answers.** It produces the standard figure of an electronic-structure paper: bands on the left and the DOS rotated on the right, sharing the energy axis and the same zero.

**Background for non-experts.** Bands say *where* in k-space the states are; the DOS says *how many* there are at each energy. Placing them side by side with the same zero (the VBM if there is a gap, the Fermi level if metallic) lets one read at a glance which orbitals form each band.

**Formulas.** None of its own: `qekit/modules/combined.py: plot` only draws; the zero is taken from the band analysis (`bands.reference_energy`) and applied to both panels.

**How Olla-DFT computes it.**

1. `qekit/cli.py: _cmd_plot` loads `bands.load` and `dos.load` on the same folder.
2. It prints `bands.gap_report`.
3. `combined.plot` creates two axes with ratio `ratio = 2.6`, calls the band-drawing logic and `dos.draw(vertical=True)` with the same energy shift.

**Where each datum comes from.**

| Datum | Source | Detail |
|---|---|---|
| Bands and gap | `prefix.xml` (see `bands`) | `bands.load` |
| DOS/PDOS | `.dos` and `pdos_atm#` (see `dos`) | `dos.load` |
| Energy zero | `bands.reference_energy(bs, ref)` | the same for both panels |

**Limits and pitfalls.**

- Bands and DOS usually come from different calculations (path vs. mesh) with **the same scf**; if they come from different scfs, their $E_F$ do not coincide and the right panel is shifted without any warning.
- The zero of the DOS is that of the bands (VBM in `auto` with a gap), whereas `olla-dft dos` on its own would use the Fermi level.

**References.**

- Those of `bands` and `dos`.

---
### `olla-dft effmass` — effective mass by parabolic fit

**What it answers.** How much an electron at the bottom of the conduction band and a hole at the top of the valence band "weigh": the effective mass $m^*/m_e$ in each direction, which governs mobilities, effective densities of states and exciton levels.

**Background for non-experts.** Near an extremum a band looks like a parabola, just like the kinetic energy of a free particle $E = \hbar^2 k^2/2m$. The curvature of that parabola defines a mass: a strongly curved ("open") band corresponds to a light, fast carrier; a flat band to a heavy one. The mass may depend on direction (in silicon the electron has a longitudinal mass of ~0.92 and a transverse one of ~0.19), so parabolas must be fitted along specific straight lines in k-space.

Olla-DFT works in two stages. First it fits on the bands you already have (fast, but with few points and only along the path directions). Then it writes a dedicated `bands` calculation with very fine lines crossing the VBM and the CBM in three directions (for a valley away from Γ: the radial Γ→k₀ direction, "longitudinal", and two perpendicular ones, "transverse"; at Γ: [100], [110], [111]) and, when it finishes, fits one parabola per line.

**Formulas.**

Quadratic fit and mass (`qekit/modules/effmass.py: from_bands`, `collect_fine`, `_mass_from_quadratic`):

$$
E(k) \approx a\,k^2 + b\,k + c, \qquad \frac{m^*}{m_e} = \frac{\hbar^2/m_e}{2a}, \qquad \frac{\hbar^2}{m_e} = 7.6199682\ \text{eV·Å}^2
$$

- $k$: signed distance to the extremum along the line (Å⁻¹); $a$ in eV·Å²; the fit is `np.polyfit(x, y, 2)`.
- The sign is kept: $a<0$ (downward curvature) gives $m^*<0$, which is the report's convention for a hole.
- The linear term $b$ is fitted but **not** used: the extremum is assumed at $k=0$.

Fit quality (`effmass._r2`): $R^2 = 1 - \sum(y - \hat y)^2 / \sum (y - \bar y)^2$.

Valley directions (`effmass.valley_directions`): if $|\mathbf k_0| < 10^{-6}$ Å⁻¹ (extremum at Γ), $\{[100], [110]/\sqrt2, [111]/\sqrt3\}$; otherwise $\hat e_1 = \mathbf k_0/|\mathbf k_0|$ and two perpendiculars built with cross products.

Points of the fine line (`effmass.prepare`): $\mathbf k_j = \mathbf k_0 + t_j\,\hat e$, $t_j \in [-h, h]$ with `half_width = 0.06` Å⁻¹ and `npts = 21` (forced odd), converted to fractional coordinates with $\mathbf k_{\text{frac}} = \mathbf k\,\mathbf B^{-1}$.

Identification of the valence band in the fine calculation (`collect_fine`): $n_v = \mathrm{round}(N_{\text{el}}/2) - 1$ (0-based index) if `nspin = 1`; otherwise the last band whose maximum lies below $E_F$ (or the HOMO, or the median).

**How Olla-DFT computes it.**

1. `qekit/cli.py: _cmd_effmass` requires `--bands-dir` (a folder with a finished bands calculation) unless `--collect`. It loads the structure and `bands.load(bands_dir)`.
2. `effmass.from_bands`: `bands.analyze_gap` gives VBM/CBM; for holes it takes the VBM band and those lying within `DEGEN_TOL = 0.05` eV below it at that k; for electrons the CBM band and those degenerate above it. For each band it relocates the extremum (`argmax`/`argmin`), delimits the segment without crossing special points or discontinuities (`_segment_bounds`), collects the points with $|k - k_0| \le$ `--window` (a half-width; by default `WINDOW_DEFAULT = PARABOLIC_MAX/2 = 0.06` Å⁻¹, widening up to `--min-points` = 7) on both sides if the extremum is interior or on one side if it sits on a special point (`_collect_window`), and fits.
3. It warns if there are fewer than 5 points ("solo N puntos: el ajuste no es confiable; haz el cálculo dedicado (effmass sin --collect y luego --collect)") or if the **total fitted span** ($k_{\max} - k_{\min}$, `MassFit.window`) exceeds `PARABOLIC_MAX = 0.12` Å⁻¹ ("tramo ajustado de X Å⁻¹ (límite parabólico 0.12): el camino no tiene puntos más finos; haz el cálculo dedicado"). With the default half-width a fit centred on the extremum spans exactly 0.12, so the warning fires only when the window had to be widened for lack of points.
4. `effmass.prepare` reduces to the primitive cell (`structure.primitive`), resolves pseudopotentials and cutoffs with `sweep.prepare_common`, writes `masa.in` (`calculation='bands'`, `K_POINTS crystal` with the 6 lines, `nbnd` equal to that of the bands calculation, `occupations='fixed'` because `insulator=True`) and `scf.in` (mesh `sweep.default_grid`), and saves `masa_meta.json` with the description of every line.
5. With `--run`, `runner.run_all` launches pw.x on `scf.in` and `masa.in`; with `--collect`, `effmass.collect_fine` reads `out/*.xml`, slices the k list into chunks of `npts`, computes $t_j = \pm|\mathbf k_j - \mathbf k_c|$ and fits all bands degenerate with the extremum.
6. `effmass.report` prints the table (carrier, band, $m^*/m_e$, $R^2$, points, Δk, direction) and `export` writes `MASA_EFECTIVA.dat`.

**Where each datum comes from.**

| Datum | Source | Detail |
|---|---|---|
| Eigenvalues and Cartesian k | `prefix.xml` from pw.x | `qeout.read_xml` (previous bands and fine calculation) |
| VBM, CBM and their k | `bands.analyze_gap` | see `olla-dft bands` |
| $\hbar^2/m_e$ | constant `effmass.HBAR2_OVER_ME` | 7.6199682 eV·Å² |
| Number of electrons | `<nelec>` from the XML | to identify the valence in `collect_fine` |
| Window, minimum points, half-width, points per line | user parameters | `--window` (default `effmass.WINDOW_DEFAULT` = 0.06 Å⁻¹), `--min-points 7`, `--half-width 0.06`, `--points 21` |
| Parabolic limit | `effmass.PARABOLIC_MAX` | 0.12 Å⁻¹ of total span (slack `_TOL_VENTANA = 1e-6`) |
| Degeneracy tolerance | `effmass.DEGEN_TOL` | 0.05 eV |
| Description of the lines | `masa_meta.json` written by `prepare` | `effmass.load_meta` |

**Limits and pitfalls.**

- There is no "quick fit only" mode: with `--bands-dir` the command always fits on the path and then prepares the fine calculation in `--outdir` (or runs it with `--run`); with `--collect` it reads the fine one. This is what the module docstring describes.
- For a metal the command stops: "El sistema es metálico: no hay un extremo de banda aislado que ajustar".
- The report warns that the calculation **does not include spin-orbit coupling**: "cerca de Γ hay un triplete degenerado, no el par hueco pesado / hueco ligero del modelo de Luttinger".
- An $R^2$ of 1.0000 with 3 or 4 points "no dice nada — una parábola pasa exacta por tres puntos cualesquiera" (report text).
- `--window` is a **half-width**: the fitted span is up to twice it, and it is that span that is compared with `PARABOLIC_MAX`. Raising `--window` above 0.06 Å⁻¹ triggers the non-parabolic warning even if the fit looks "good" in $R^2$.
- `--collect` uses the **first** XML in `out/`; with several prefixes it may read the wrong one.
- The fine calculation is always written with `occupations='fixed'` and without spin; for magnetic systems `masa.in` must be edited by hand.
- The transverse lines for a valley away from Γ are chosen with an arbitrary cross product: they are not necessarily crystallographic axes, and in a non-ellipsoidal valley the two transverse masses may differ.

**References.**

- N. W. Ashcroft, N. D. Mermin, *Solid State Physics* (1976), ch. 12 — definition of effective mass.
- J. M. Luttinger, W. Kohn, *Phys. Rev.* **97**, 869 (1955) — k·p model of degenerate valence bands.
- Silicon reference values: M. Cardona, F. H. Pollak, *Phys. Rev.* **142**, 530 (1966).

---

### `olla-dft fermi` — Fermi surface in BXSF format

**What it answers.** Which bands cross the Fermi level and what the surface $\varepsilon_n(\mathbf k) = E_F$ of each one looks like, written to a BXSF file that XCrySDen or FermiSurfer render in 3D.

**Background for non-experts.** In a metal the states fill up to an energy $E_F$; the set of k-points whose energy is exactly $E_F$ forms a surface in the Brillouin zone, the Fermi surface. Its shape determines conductivity, quantum oscillations and many instabilities (charge-density waves, superconductivity). To draw it one needs $\varepsilon_n(\mathbf k)$ on a **complete and uniform** mesh of the Brillouin zone, which is what `olla-dft transport` produces (an nscf with `nosym`, `noinv`).

**Formulas.**

Bands crossing $E_F$ (`qekit/modules/transport.py: crossing_bands`), with `tol = 1e-6` eV:

$$
\min_{\mathbf k}\varepsilon_n(\mathbf k) < E_F - \delta \quad\wedge\quad \max_{\mathbf k}\varepsilon_n(\mathbf k) > E_F + \delta
$$

Mesh reconstruction (`transport.load`): fractional coordinates are brought to $[0,1)$ with $f \leftarrow f - \lfloor f + 10^{-6}\rfloor$, rounded to 6 decimals, and $n_i$ is the number of distinct values along each axis; $n_1 n_2 n_3 = N_k$ is required.

BXSF grid (`transport.export_bxsf`): $(n_i + 1)$ points per axis are written, repeating the first plane at the end (`np.pad(..., mode="wrap")`), in C order (last index fastest), with the reciprocal vectors $\mathbf b_i = 2\pi(\mathbf A^{-1})^{\mathsf T}_i$ in Å⁻¹ and the energies in eV.

**How Olla-DFT computes it.**

1. `qekit/cli.py: _cmd_fermi` looks for `out/*.xml` inside `--outdir` (default `transporte`).
2. `transport.load` reads the XML with `qeout.read_xml`; it rejects an `scf`-type XML ("es de un cálculo SCF, no del nscf de malla densa"); it reconstructs the mesh and reorders the energies with `np.lexsort`. It also computes band velocities by finite differences (unused here) and warns if the mesh is smaller than 24×24×24 or has fewer than 12 000 points.
3. `transport.crossing_bands` lists the metallic bands; if none, it prints "Ninguna banda cruza E_F: el sistema no es metálico y no tiene superficie de Fermi".
4. `transport.export_bxsf` writes `superficie_fermi.bxsf` with `Fermi Energy`, the grid and one `BAND:` block per band.

**Where each datum comes from.**

| Datum | Source | Detail |
|---|---|---|
| Eigenvalues on the mesh | `prefix.xml` of the `olla-dft transport` nscf | `qeout.read_xml`; spin channel 0 only |
| Fermi energy | `<fermi_energy>` from the XML | `run.fermi`; without it, `ErrorDeUso("no hay nivel de Fermi…")` |
| Cell | `<atomic_structure>` from the XML | reciprocal vectors with $2\pi$ |
| Crossing tolerance | argument `tol` of `crossing_bands` | 1e-6 eV |

**Limits and pitfalls.**

- It only works on the `olla-dft transport` folder (same full-mesh nscf); a band path or a symmetry-reduced mesh fails with "los N puntos k no forman una malla uniforme".
- Only spin channel 0 is exported (`transport.load(spin=0)`); a ferromagnetic metal would need two files and the command does not produce them.
- $E_F$ is taken as-is from the nscf XML; with `occupations='fixed'` it does not exist and the command fails.
- The Fermi level of a dense nscf is not recomputed: it is the one inherited from the scf (coarser mesh).
- The reciprocal vectors are written in Å⁻¹ with the factor $2\pi$; the viewer must interpret them in those units.

**References.**

- A. Kokalj, *Comput. Mater. Sci.* **28**, 155 (2003) — XCrySDen and the BXSF format. DOI 10.1016/S0927-0256(03)00104-6.
- M. Kawamura, *Comput. Phys. Commun.* **239**, 197 (2019) — FermiSurfer. DOI 10.1016/j.cpc.2019.01.017.

---

### `olla-dft unfold` — band unfolding of a supercell

**What it answers.** What fraction of every supercell state "belongs" to each k-point of the primitive cell: the spectral weight that lets one see the band of the original material (and how much a defect, a dopant or disorder blurs it) from a supercell calculation.

**Background for non-experts.** A supercell of $N$ primitive cells has a Brillouin zone $N$ times smaller, so its bands come out *folded*: where the primitive cell had one band, there are $N$ branches piled up. Every supercell state is a sum of plane waves $e^{i(\mathbf K + \mathbf G)\cdot\mathbf r}$, and every plane wave has a well-defined wave vector. Asking "how much of this state lives at point $\mathbf k$ of the primitive cell?" has an exact answer: the sum of $|C(\mathbf G)|^2$ over the plane waves whose $\mathbf K + \mathbf G$ coincides with $\mathbf k$ modulo the primitive reciprocal lattice. If the supercell is perfect, every state has weight 1 at a single $\mathbf k$ and the primitive band is recovered; if there is a defect, the weight spreads and the band looks blurred. That blurring is the physical result.

**Formulas.**

Supercell matrix (`qekit/modules/unfold.py: matriz_supercelda`): $\mathbf M = \mathbf A_{\text{sc}}\,\mathbf a_{\text{prim}}^{-1}$, rounded to integers; accepted if $\max|\mathbf M - \mathrm{round}(\mathbf M)| \le 10^{-3}$. If it fails because of orientation, `_m_por_metricas` searches for an integer $\mathbf M$ such that $\mathbf G_{\text{sc}} = \mathbf M\,\mathbf G_p\,\mathbf M^{\mathsf T}$ with $\mathbf G = \mathbf X\mathbf X^{\mathsf T}$ the metric tensor (rotation-invariant), row by row among integer vectors of the right length. The primitive cell is then **re-derived** as $\mathbf a = \mathbf M^{-1}\mathbf A_{\text{sc}}$ so that both share axes. $N = |\det\mathbf M|$.

Coordinates: since $\mathbf B_{\text{sc}} = \mathbf M^{-\mathsf T}\mathbf b_{\text{prim}}$, a vector with coordinates $\mathbf c_{\text{sc}}$ in the supercell reciprocal basis has coordinates $\mathbf c_p = \mathbf c_{\text{sc}}\mathbf M^{-\mathsf T}$ in the primitive one, and $\mathbf k_{\text{prim}} = \mathbf k_{\text{sc}}\mathbf M^{-\mathsf T}$ (`desdoblar`).

Spectral weight (`unfold.pesos_de_k`), with $\mathbf m_0 = \mathbf k_{\text{prim}}\mathbf M^{\mathsf T} - \mathbf k_{\text{sc}}$ (must be integer to `TOL_ENTERO = 1e-4`; otherwise the weight is 0 because that $\mathbf k$ does not fold onto this $\mathbf K$):

$$
P_{n}(\mathbf k) = \frac{\sum_{\mathbf G \in S}\ \sum_{\sigma}|C_{n\sigma}(\mathbf G)|^2}{\sum_{\mathbf G}\ \sum_{\sigma}|C_{n\sigma}(\mathbf G)|^2}, \qquad
S = \left\{\mathbf G : (\mathbf G - \mathbf m_0)\,\mathbf M^{-\mathsf T} \in \mathbb{Z}^3\right\}
$$

- $C_{n\sigma}(\mathbf G)$: plane-wave coefficients of band $n$ (spinor component $\sigma$ if `npol = 2`) read from `wfc<N>.dat`; $\mathbf G$ given by its Miller indices in the supercell reciprocal basis.
- The denominator normalizes in case the coefficients are not normalized; $P_n \in [0,1]$.

Distance on the x axis (`unfold._distancias`): sum of $|\Delta\mathbf k|$ with $\mathbf k = \mathbf k_{\text{frac}}\,\mathbf b_{\text{prim}}$; a jump larger than 5 times the median of the non-zero steps counts as zero (branch change).

**How Olla-DFT computes it.**

1. `qekit/cli.py: _cmd_unfold` loads the primitive structure and calls `unfold.desdoblar(path, primitive_cell, bandas=range(--bands), spin=--spin)`; `--spin` is `up` (default) or `dw`.
2. `desdoblar` reads the XML (`qeout.read_xml`), locates the `.save` folder (`_carpeta_save`) and the wavefunction files of the requested channel (`qekit/core/wfc.py: buscar_wfc(save, spin)`, sorted by k number): if `wfc.es_lsda` detects `wfcup*`/`wfcdw*`, it returns only the `wfc{up|dw}<N>.dat` of that channel; otherwise the `wfc<N>.dat` of a spin-unpolarized run. Without them: "El cálculo no guardó las funciones de onda: eso pasa con disk_io='nowf' o 'low'" (or, with lsda, "falta el canal '…'").
3. `matriz_supercelda` obtains $\mathbf M$; the k-points of the calculation are converted to primitive coordinates.
4. For every k-point, `wfc.leer_wfc` reads the unformatted Fortran file: record 1 (`ik`, `xk`, `ispin`, `gamma_only`, `scalef`), record 2 (`ngw`, `igwx`, `npol`, `nbnd`), record 3 (`b1,b2,b3`), record 4 (Miller indices) and one record per band with `npol·igwx` complex numbers (only the requested bands are materialized).
5. `pesos_de_k` computes $P_n(\mathbf k)$; the energies come from the XML, from the same spin channel as the wavefunctions (`res.eigenvalues[0]` for `up`, `[1]` for `dw`).
6. `unfold.report` prints $N$, $\mathbf M$, the weight distribution (mean, fraction > 0.9, fraction < 0.1) and warnings; `export` writes `UNFOLD.dat` (distance, $E - E_F$, weight) and `UNFOLD.txt`; `plot` draws a `scatter` with size $= 60\,P$ for weights > 0.005.

**Where each datum comes from.**

| Datum | Source | Detail |
|---|---|---|
| Coefficients $C(\mathbf G)$, Miller indices, `npol` | `out/<prefix>.save/wfc<N>.dat` from pw.x | `wfc.leer_wfc` (sequential Fortran format, little-endian) |
| Eigenvalues, fractional k, $E_F$, supercell cell | `prefix.xml` | `qeout.read_xml` |
| Spin channel | `--spin up|dw` (user) | `wfc.buscar_wfc`, `wfc.es_lsda`; in a spin-unpolarized run it changes nothing |
| Primitive cell | user file | re-derived as $\mathbf M^{-1}\mathbf A_{\text{sc}}$ |
| Integer tolerance | `unfold.TOL_ENTERO` | 1e-4 (and 1e-3 to accept $\mathbf M$) |

**Limits and pitfalls.**

- `disk_io='medium'` or `'high'` is needed in the supercell bands calculation; `olla-dft gen` does not set it by default.
- **One spin channel per run** is unfolded. For an `lsda` calculation the report warns: "el cálculo es de espín polarizado (lsda) y aquí solo se ha desdoblado el canal 'up' (wfcup<N>.dat y sus energías). El otro canal no se mezcla ni se suma: para verlo repite el desdoblamiento con --spin dw". The two channels are never combined in a single figure.
- If the supercell is relaxed and the primitive is not (or vice versa), $\mathbf M$ is not integer: "la celda de la supercelda no es un múltiplo entero de la primitiva (error …)".
- If almost all weights are 1, the report warns: "la supercelda parece PERFECTA (sin defecto ni desorden). En ese caso el desdoblamiento reproduce exactamente las bandas primitivas — que es la comprobación de que funciona, pero no un resultado nuevo".
- The k-points of the calculation are interpreted as supercell k-points and converted to the primitive; no primitive path is generated and no check is made that the supercell k-points are the correct folds of the desired path.
- Ultrasoft/PAW wavefunctions are handled only through their plane-wave part (the augmentation term $S$ is not included); the weight is that of the smooth part.

**References.**

- V. Popescu, A. Zunger, *Phys. Rev. B* **85**, 085201 (2012) — unfolding spectral weight. DOI 10.1103/PhysRevB.85.085201.
- P. V. C. Medeiros, S. Stafström, J. Björk, *Phys. Rev. B* **89**, 041407(R) (2014) — plane-wave unfolding (BandUP). DOI 10.1103/PhysRevB.89.041407.
- W. Ku, T. Berlijn, C.-C. Lee, *Phys. Rev. Lett.* **104**, 216401 (2010).

---
### `olla-dft wannier` — Wannier functions and band interpolation

**What it answers.** From a DFT calculation on a coarse k-point mesh it builds a small model $H_{mn}(\mathbf R)$ in a basis of localized functions (Wannier functions) with which the bands can be evaluated at **any** k-point without re-running pw.x; it also gives where each function is centred, how far it extends (its spread $\Omega$), and how closely the interpolated band matches the DFT one.

**Background for non-experts.** Bloch states $\psi_{n\mathbf k}$ are delocalized over the whole crystal. Their Fourier transform in k gives functions $|\mathbf R n\rangle$ localized around a cell $\mathbf R$: the Wannier functions. In that basis the Hamiltonian is a small matrix $H(\mathbf R)$ that decays with $|\mathbf R|$, and transforming back to k gives the band at any point (a "Fourier interpolation" that is exact at the starting points). The difficulty is that every $\psi_{n\mathbf k}$ is defined only up to a phase (and, with degenerate bands, up to a unitary rotation among them): that freedom is called the *gauge*. With an arbitrary gauge the Wannier functions are not localized and the interpolation is garbage. Marzari and Vanderbilt proposed choosing the gauge that minimizes the total spread $\Omega$ (the sum of the squared "widths"); a good starting point is to project onto trial atomic orbitals and orthonormalize.

When the bands of interest cross others (metals, conduction bands) there is no isolated group to transform: one must choose at every k a subspace of $J$ states that "connects smoothly" with that of its neighbours. That is the *disentanglement* of Souza, Marzari and Vanderbilt, with an *outer* window (where one may choose from) and optionally a *frozen* one (states that must be kept exactly). Olla-DFT implements both in Python, using only the overlaps and projections computed by `pw2wannier90.x` (shipped with QE), without needing wannier90; if the user has wannier90, it also reads its `seedname_hr.dat`.

**Formulas.**

Full mesh and $\mathbf b$ vectors (`qekit/modules/wannier.py: malla_completa`, `capas_b`, `residuo_completitud`): $\mathbf k_{ijk} = (i/n_1, j/n_2, k/n_3)$ in QE's order (last index fastest). Neighbour shells $\mathbf b = (h_1/n_1, h_2/n_2, h_3/n_3)\,\mathbf B$ are added by distance until, by least squares over the 6 independent components,

$$
\sum_{\mathbf b} w_{\mathbf b}\, b_\alpha b_\beta = \delta_{\alpha\beta}, \qquad \text{residual } = \left\|\textstyle\sum_{\mathbf b} w_{\mathbf b}\,\mathbf b\otimes\mathbf b - \mathbf 1\right\|_\infty < 10^{-5}
$$

- $w_{\mathbf b}$: weight of each shell (Å²); shells that add no rank (SVD) or with $|w| < 10^{-8}$ are discarded.

Projection gauge (`wannier.gauge_proyeccion`), with $A_{mn}(\mathbf k) = \langle\psi_{m\mathbf k}|g_n\rangle$ from the `.amn` and the SVD $A = u\,s\,v^\dagger$:

$$
U(\mathbf k) = A\,(A^\dagger A)^{-1/2} = u\,v^\dagger
$$

- $U$: $N_b\times J$ matrix with orthonormal columns (Löwdin); the smallest singular value $s_{\min}$ is reported (warning if $< 0.2$).

Gauge-invariant spread and disentanglement (`wannier.omega_I`, `gauge_desenredo`), with $M^{\mathbf k,\mathbf b}_{mn} = \langle u_{m\mathbf k}|u_{n,\mathbf k+\mathbf b}\rangle$ from the `.mmn`:

$$
\Omega_I = \frac{1}{N_k}\sum_{\mathbf k}\sum_{\mathbf b} w_{\mathbf b}\left[J - \left\|U^\dagger(\mathbf k)\,M^{\mathbf k,\mathbf b}\,U(\mathbf k+\mathbf b)\right\|_F^2\right]
$$

$$
Z(\mathbf k) = \sum_{\mathbf b} w_{\mathbf b}\, M^{\mathbf k,\mathbf b}\,U(\mathbf k+\mathbf b)\,U^\dagger(\mathbf k+\mathbf b)\,M^{\mathbf k,\mathbf b\,\dagger}
$$

- At every iteration $Z$ is mixed with the previous one ($Z \leftarrow \mu Z_{\text{new}} + (1-\mu)Z_{\text{old}}$, $\mu = 0.5$ initially, halved if $\Omega_I$ goes up), restricted to the bands in the outer window, the frozen ones are projected out ($Q Z Q$ with $Q = 1 - P_{\text{frozen}}$) and the $J - N_{\text{frozen}}$ eigenvectors of largest eigenvalue are taken. At most 200 steps, tolerance $10^{-10}$ Å². At the end it re-projects onto the trial orbitals inside the subspace ($U \leftarrow U\,\mathrm{polar}(U^\dagger A)$) to obtain a smooth starting gauge.

Real-space Hamiltonian and interpolation (`wannier.hamiltoniano_k`, `a_reales`, `interpolar`, `celda_wigner_seitz`):

$$
H(\mathbf k) = U^\dagger(\mathbf k)\,\mathrm{diag}\big(\varepsilon_n(\mathbf k)\big)\,U(\mathbf k), \qquad
H(\mathbf R) = \frac{1}{N_k}\sum_{\mathbf k} e^{-2\pi i\,\mathbf k\cdot\mathbf R}\,H(\mathbf k), \qquad
H^{\text{int}}(\mathbf k) = \sum_{\mathbf R}\frac{e^{2\pi i\,\mathbf k\cdot\mathbf R}}{\deg(\mathbf R)}\,H(\mathbf R)
$$

- $\mathbf k$ and $\mathbf R$ in fractional coordinates (hence the explicit $2\pi$). $\mathbf R$ runs over the vectors of the Wigner-Seitz cell of the $n_1\times n_2\times n_3$ superlattice; $\deg(\mathbf R)$ is the number of equidistant images (tolerance $10^{-5}$ Å²). The bands are the eigenvalues of $\tfrac12(H^{\text{int}} + H^{\text{int}\dagger})$.

Centres and spread (`wannier.dispersion`), equations 31 and 34–36 of Marzari-Vanderbilt, with $M^W = U^\dagger(\mathbf k) M^{\mathbf k,\mathbf b} U(\mathbf k+\mathbf b)$ and $\phi_n = \operatorname{Im}\ln M^W_{nn}$:

$$
\bar{\mathbf r}_n = -\frac{1}{N_k}\sum_{\mathbf k,\mathbf b} w_{\mathbf b}\,\mathbf b\,\phi_n, \qquad
\Omega_n = \frac{1}{N_k}\sum_{\mathbf k,\mathbf b} w_{\mathbf b}\left[\left(1-|M^W_{nn}|^2\right) + \phi_n^2\right] - |\bar{\mathbf r}_n|^2
$$

$$
\Omega_I = \frac{1}{N_k}\sum_{\mathbf k,\mathbf b} w_{\mathbf b}\Big[J - \sum_{mn}|M^W_{mn}|^2\Big], \quad
\Omega_{OD} = \frac{1}{N_k}\sum_{\mathbf k,\mathbf b} w_{\mathbf b}\sum_{m\ne n}|M^W_{mn}|^2, \quad
\Omega_D = \frac{1}{N_k}\sum_{\mathbf k,\mathbf b} w_{\mathbf b}\sum_n\left(\phi_n + \mathbf b\cdot\bar{\mathbf r}_n\right)^2
$$

- $\bar{\mathbf r}_n$ in Å, $\Omega$ in Å²; $\Omega = \sum_n\Omega_n = \Omega_I + \Omega_D + \Omega_{OD}$ (the report prints the sum as a check).

Minimization (`wannier._gradiente`, `_rotar`, `minimizar`), eqs. 52–57 of Marzari-Vanderbilt, with $R_{mn} = M_{mn}M^*_{nn}$, $T_{mn} = (M_{mn}/M_{nn})\,q_n$, $q_n = \phi_n + \mathbf b\cdot\bar{\mathbf r}_n$, $\mathcal A(B) = (B - B^\dagger)/2$, $\mathcal S(B) = (B + B^\dagger)/2i$:

$$
G(\mathbf k) = -\frac{4}{N_k}\sum_{\mathbf b} w_{\mathbf b}\left[\mathcal A(R^{\mathbf k,\mathbf b}) - \mathcal S(T^{\mathbf k,\mathbf b})\right], \qquad
U(\mathbf k) \leftarrow U(\mathbf k)\,\exp\!\left(-\Delta t\,G(\mathbf k)\right), \qquad
\Delta t_0 = \frac{\alpha}{4\sum_{\mathbf b} w_{\mathbf b}},\ \alpha = 2
$$

- If the step raises $\Omega$ it is halved up to 12 times; at most 500 steps (`--iterations`); stop when $|\Delta\Omega| < 10^{-10}$. It checks that $\Omega_I$ does not change (`deriva_I`).

Interpolated DOS (`wannier.dos_interpolada`): $\rho(E) = \frac{1}{N_k\,\sigma\sqrt{2\pi}}\sum_{\mathbf k,n} e^{-(E-\varepsilon_n(\mathbf k))^2/2\sigma^2}$ on an $N^3$ mesh (`--dos N`), $\sigma$ = `--sigma` 0.05 eV; it integrates to $J$ states per cell, **without** the spin factor 2. The header of `WANNIER_dos.dat` declares it via `wannier.DOS_UNIDADES`: "estados/eV/celda, sin factor de espín: integra a num_wann (x2 para comparar con dos.x sin espín)".

**How Olla-DFT computes it.**

1. *Prepare* (`qekit/cli.py: _cmd_wannier` → `wannier.prepare`): it translates `--projections` (`Si:sp3`, `O:p;Ti:d`, `f=0.25,0.25,0.25:s`, or `auto` = $s$ and $p$ on every atom) into $(l, m_r)$ orbitals from the `ORBITALES` table (wannier90 convention); it writes `1_scf.in`, `2_nscf.in` (full `--grid` mesh, 4×4×4 by default, `K_POINTS crystal`, `nosym`, `noinv`, `conv_thr 1e-10`, `nbnd = --bands` or $J$ + excluded), the `<prefix>.nnkp` (`escribir_nnkp`: real and reciprocal lattice, k-points, projections, `nnkpts` neighbours with their $\mathbf G$, `exclude_bands`), `3_pw2wan.in` (`write_amn`, `write_mmn`, `write_unk=.false.`), `<prefix>.win` (in case one prefers wannier90) and `4_bands.in` (DFT bands along the seekpath path, 30 points per segment, `outdir='./out_bandas'`).
2. *Run* (`--run` → `wannier.correr`): `pw.x` (scf), `pw.x` (nscf), `pw2wannier90.x`, `pw.x` (bands), in that order, stopping at the first failure.
3. *Collect* (`--collect` → `wannier.collect`): reads `.eig` (`leer_eig`), `.amn` (`leer_amn`), `.mmn` (`leer_mmn`, with $m$ running fastest → `reshape(order="F")`); recomputes shells and neighbours from the `.nnkp` (`_leer_nnkp`); if there are more bands than functions or `--window`/`--frozen` were given, `gauge_desenredo`; otherwise `gauge_proyeccion`. `dispersion` before and after `minimizar` (unless `--no-minimize`). `celda_wigner_seitz`, `hamiltoniano_k`, `a_reales`; it checks that `interpolar` reproduces the mesh (`error_malla` < `TOL_EXACTA = 1e-6` eV) and, if DFT bands exist (`out_bandas`, `--dft-bands`), compares at points that were not in the mesh. As a negative control, it repeats the interpolation with $U = 1$ (`E_sin_gauge`).
4. If a wannier90 `*_hr.dat` exists (other than Olla-DFT's own `WANNIER_hr.dat`), `leer_hr` uses it directly and skips the localization.
5. `wannier.report` prints mesh, neighbours and residual, windows, $\Omega_I$, the decomposed $\Omega$, centres with assignment to atom or bond (`asignar`, bond window 0.5–3.2 Å), decay of $H(\mathbf R)$, exactness on the mesh and error against DFT; `export` writes `WANNIER_hr.dat` (wannier90 format), `WANNIER_centros.dat`, `WANNIER_bandas.dat`, `WANNIER.txt` and optionally `WANNIER_dos.dat`; `plot` draws Wannier over DFT and the $\Omega$ trace.

**Where each datum comes from.**

| Datum | Source | Detail |
|---|---|---|
| Energies $\varepsilon_n(\mathbf k)$ | `seedname.eig` from pw2wannier90.x | absolute eV; `wannier.leer_eig` |
| Projections $A_{mn}(\mathbf k)$ | `seedname.amn` | `wannier.leer_amn` |
| Overlaps $M^{\mathbf k,\mathbf b}_{mn}$ | `seedname.mmn` | `wannier.leer_mmn` |
| Cell, mesh, excluded bands | `seedname.nnkp` (written by Olla-DFT) | `wannier._leer_nnkp` |
| External $H(\mathbf R)$ | `seedname_hr.dat` from wannier90 | `wannier.leer_hr` |
| DFT validation bands | `out_bandas/*.xml` from step 4 | `qeout.read_xml`, channel 0 |
| Trial orbitals $(l, m_r)$ | table `wannier.ORBITALES` | Table 3.1/3.2 of the wannier90 manual |
| Tolerances | `TOL_COMPLETITUD 1e-5`, `TOL_PESO 1e-8`, `TOL_EXACTA 1e-6` eV | module constants |
| High-symmetry path | seekpath (`wannier.camino_denso`) | 30 points per segment (`--points`) |

**Limits and pitfalls.**

- The `--window` and `--frozen` windows are compared with the **absolute** energies of the `.eig` (not relative to $E_F$), as in wannier90.
- With disentanglement and no frozen window, nothing has to be reproduced exactly; the report warns: "Sin ventana congelada no hay ninguna banda que la interpolación tenga que reproducir exactamente… Si quieres que la valencia salga exacta, pásala en --frozen".
- The `auto` projections ($s$ and $p$ per atom) fail for transition metals (the $d$ are missing) and for strongly covalent bonds; the report always says so.
- If $H(\mathbf R)$ at the edge of the superlattice exceeds 5 % of $H(0)$: "H(R) apenas ha decaído al borde de la superred: la base no está localizada".
- Without `--dft-bands` (or `out_bandas`) the report warns: "No has comparado con bandas de DFT. Que la interpolación reproduzca la malla es trivial".
- Only spin channel 0 of the DFT bands is read; the workflow is not designed for `nspin = 2` or SOC (pw2wannier90 supports them, but `prepare` does not write `nspin`).
- The minimization is gradient descent with line search, not wannier90's conjugate gradient: it may need more steps and may stop in a local minimum.
- Very anisotropic meshes may admit no shells satisfying completeness: "no encuentro un conjunto de capas de vecinos que cumpla la condición de completitud con esta malla".
- `--collect` without the structure as first argument fails: "para analizar hace falta la estructura".
- The interpolated DOS carries no spin factor 2 (the file header says so and asks to multiply by 2 to compare with spin-unpolarized dos.x) and is only valid within the energy range covered by the Wannier functions.

**References.**

- N. Marzari, D. Vanderbilt, *Phys. Rev. B* **56**, 12847 (1997) — maximally localized Wannier functions. DOI 10.1103/PhysRevB.56.12847.
- I. Souza, N. Marzari, D. Vanderbilt, *Phys. Rev. B* **65**, 035109 (2001) — disentanglement. DOI 10.1103/PhysRevB.65.035109.
- N. Marzari, A. A. Mostofi, J. R. Yates, I. Souza, D. Vanderbilt, *Rev. Mod. Phys.* **84**, 1419 (2012) — review. DOI 10.1103/RevModPhys.84.1419.
- G. Pizzi et al., *J. Phys.: Condens. Matter* **32**, 165902 (2020) — Wannier90 v3 (`.nnkp`, `.amn`, `.mmn`, `_hr.dat` formats). DOI 10.1088/1361-648X/ab51ff.
- P.-O. Löwdin, *J. Chem. Phys.* **18**, 365 (1950) — symmetric orthonormalization.

---

### `olla-dft topology` — Chern number and Wilson loops

**What it answers.** Whether the occupied subspace of a Wannier model, on a two-dimensional section of the Brillouin zone, has a non-zero Chern number (an integer topological invariant), and how the hybrid Wannier centres (Wilson loops) evolve across that section.

**Background for non-experts.** Besides their energies, bands have a "geometry": when traversing a closed loop in k-space, the occupied states accumulate a phase (Berry phase) that does not depend on how the phases of each state are chosen. Summing that phase over a whole 2D section of the Brillouin zone yields an integer, the Chern number, which does not change under smooth deformations of the system: it is *topological*. A non-zero Chern number implies dissipationless edge currents (quantum anomalous Hall effect). The *Wilson loop* is the "slice by slice" version: for every $k_2$ one computes the product of overlaps along $k_1$; the phases of its eigenvalues are the positions (modulo 1) of the hybrid Wannier functions, and their "pumping" as $k_2$ varies is another way of seeing the Chern number.

**Formulas.**

Mesh and states (`qekit/modules/topology.py: kmesh`, `analyze`): $\mathbf k_{ij}$ with $k_a = i/n_1$, $k_b = j/n_2$ and the third coordinate fixed at `--fixed` (mod 1), in the `--plane` (`xy`, `xz`, `yz`); the eigenvectors $|u_n(\mathbf k)\rangle$ come from `wannier.interpolar(..., vectores=True)`.

Unitary links and discrete Berry curvature (`topology._unitary_overlap`, `invariants_from_vectors`), with $V(\mathbf k)$ the $N_w\times N_{\text{occ}}$ matrix of occupied eigenvectors:

$$
O_\mu(\mathbf k) = V^\dagger(\mathbf k)\,V(\mathbf k+\hat\mu), \qquad
Q_\mu = u\,v^\dagger \ \text{(unitary part of } O_\mu = u\,s\,v^\dagger), \qquad
U_\mu(\mathbf k) = \frac{\det Q_\mu(\mathbf k)}{|\det Q_\mu(\mathbf k)|}
$$

$$
F_{12}(\mathbf k) = \arg\!\left[U_1(\mathbf k)\,U_2(\mathbf k+\hat 1)\,U_1^*(\mathbf k+\hat 2)\,U_2^*(\mathbf k)\right], \qquad
C = \frac{1}{2\pi}\sum_{\mathbf k} F_{12}(\mathbf k)
$$

- $\hat\mu$: mesh step in direction $\mu$ (periodic). $F_{12} \in (-\pi, \pi]$ per plaquette (rad); $C$ is rounded to the nearest integer and the residual $|C - \mathrm{round}(C)|$ is reported.
- The smallest singular value of all $O_\mu$ (`min_overlap`) is also reported; if $< 10^{-6}$, a warning about a too-coarse mesh.

Wilson loops (`invariants_from_vectors`):

$$
W(k_2) = \prod_{i=0}^{n_1-1} Q_1(k_1^{(i)}, k_2), \qquad
x_n(k_2) = \frac{\arg\lambda_n\!\left[W(k_2)\right]}{2\pi} \bmod 1
$$

- $x_n$: sorted hybrid Wannier centres, in fractions of the lattice vector along direction 1.

Section gaps: $E_g^{\text{dir}} = \min_{\mathbf k}[\varepsilon_{N_{\text{occ}}+1} - \varepsilon_{N_{\text{occ}}}]$, $E_g^{\text{ind}} = \min_{\mathbf k}\varepsilon_{N_{\text{occ}}+1} - \max_{\mathbf k}\varepsilon_{N_{\text{occ}}}$. $E_g^{\text{dir}} > $ `--gap-tol` (1e-8 eV) is required.

**How Olla-DFT computes it.**

1. `qekit/cli.py: _cmd_topology` requires exactly one of `--occupied N` or `--fermi EV`, and a `--grid` of at least 3×3 (40×40 by default).
2. `topology.resolve_model` accepts a `*_hr.dat` or a folder containing `WANNIER_hr.dat` (or a single `*_hr.dat`; with several, error "indica el archivo exacto").
3. `wannier.leer_hr` reads $H(\mathbf R)$, $\mathbf R$ and degeneracies; `wannier.interpolar` diagonalizes on the section mesh.
4. With `--fermi`, it counts the states with $\varepsilon < E_F$ at every k; if the count varies, error: "el nivel de Fermi corta bandas… El sistema es metálico en esta sección y el Chern de 'las ocupadas' no está definido".
5. `invariants_from_vectors` computes curvature, Chern and Wilson loops.
6. `topology.report` prints gaps, discrete and integer Chern, residual and minimum overlap; `export` writes `TOPOLOGY_curvature.dat` (flux per plaquette), `TOPOLOGY_wilson.dat` (centres vs. $k_2$) and `TOPOLOGY.txt`; `plot` draws the flux map and the centres.

**Where each datum comes from.**

| Datum | Source | Detail |
|---|---|---|
| $H(\mathbf R)$, $\mathbf R$, $\deg(\mathbf R)$ | `WANNIER_hr.dat` (Olla-DFT) or `seedname_hr.dat` (wannier90) | `wannier.leer_hr` |
| Eigenvectors on the mesh | `wannier.interpolar` | fractional coordinates, phase $e^{2\pi i\mathbf k\cdot\mathbf R}$ |
| Occupation | `--occupied` or `--fermi` (user) | never guessed |
| Gap tolerance | `--gap-tol` | 1e-8 eV |

**Limits and pitfalls.**

- 2D sections only: for a 3D material `--fixed` must be scanned by hand; the Chern number of a section is the invariant of a 2D Chern insulator or of a slice.
- "La señal cambia al invertir la orientación del plano" (report text): the sign of $C$ depends on the `(a, b)` order of the chosen plane.
- $\mathbb Z_2$ is not computed: "no se asigna un Z2 automático sin comprobar simetría de reversión temporal". With time-reversal symmetry the Chern number is always 0; the exported Wilson loops let one read $\mathbb Z_2$ by eye, but the code does not do it.
- If the direct gap closes on the mesh: "el subespacio ocupado no está aislado… El número de Chern no está definido".
- If the discrete Chern does not close to an integer within $10^{-6}$: "refina la malla y revisa la localización del modelo Wannier".
- The result inherits every defect of the Wannier model (bad projections, undecayed $H(\mathbf R)$).

**References.**

- T. Fukui, Y. Hatsugai, H. Suzuki, *J. Phys. Soc. Jpn.* **74**, 1674 (2005) — discrete Chern number on a mesh. DOI 10.1143/JPSJ.74.1674.
- R. Yu, X. L. Qi, A. Bernevig, Z. Fang, X. Dai, *Phys. Rev. B* **84**, 075119 (2011) — Wilson loops and hybrid centres.
- A. A. Soluyanov, D. Vanderbilt, *Phys. Rev. B* **83**, 235401 (2011) — hybrid Wannier centres and invariants.
- D. Vanderbilt, *Berry Phases in Electronic Structure Theory* (Cambridge, 2018).
- X.-L. Qi, Y.-S. Wu, S.-C. Zhang, *Phys. Rev. B* **74**, 085308 (2006) — test model used in `tests/test_topology.py`.

---
### `olla-dft berry` — Berry-phase polarization, Born charges

**What it answers.** How much the electric polarization of an insulating crystal changes when going from a reference structure (usually the centrosymmetric one) to the polar one — the spontaneous polarization of a ferroelectric — and how much effective charge "moves" when an atom is displaced (Born effective charge $Z^*$).

**Background for non-experts.** The polarization of a periodic solid **cannot** be computed as the dipole moment of the cell: that number depends on where the cell boundaries are cut. King-Smith and Vanderbilt showed that what is well defined is a geometric phase (Berry phase) accumulated by the occupied states when traversing the Brillouin zone along one direction. That phase is defined modulo $2\pi$, so the polarization is defined modulo a "quantum" $e\mathbf R/\Omega$: only **differences** between two structures connected by a path are measurable, exactly as in experiment (one measures the charge that flows while the structure changes, not $P$). pw.x computes that phase with `lberry = .true.` on "strings" of k-points parallel to a reciprocal vector; Olla-DFT prepares the strings correctly, follows the branch of the phase along the path, and checks the ionic part against its exact formula.

**Formulas.**

k strings (`qekit/modules/berry.py: cuerdas`): for every point $(i/n_\perp^{(1)}, j/n_\perp^{(2)})$ of the perpendicular mesh (`--kperp` 6×6), `nppstr` points (9 by default) along $\mathbf b_{\text{gdir}}$ with coordinate $l/(n_{\text{pp}}-1)$, $l = 0,\dots,n_{\text{pp}}-1$: the last point is the first one plus $\mathbf G$.

Ionic phase (`berry.fase_ionica`), in QE's units (the quantum is `MOD_TOT` = 2 if all valences are even, 1 if any is odd; `berry.modulo_de`):

$$
\varphi_{\text{ion}} = \sum_a \left[Z_a f_a^{(g)}\right]_{\bmod\, m_a}\Big|_{\bmod\, m}, \qquad m_a = \begin{cases}1 & Z_a \text{ odd}\\ 2 & Z_a \text{ even}\end{cases}
$$

- $Z_a$: valence charge of the pseudopotential of atom $a$ (electrons); $f_a^{(g)}$: fractional coordinate along `gdir`. The per-ion folding and the final folding reproduce what pw.x does; folding to $[-m/2, m/2)$ uses Fortran's `NINT` (`berry._nint`, half rounds away from zero), so that half a quantum comes out as $-1$ as in QE.

Electronic phase from Wannier centres (`berry.desde_wannier`), as an independent check:

$$
\varphi_{\text{el}} = -f_s\sum_n \bar r_n^{(g)}, \qquad f_s = 2
$$

- $\bar r_n^{(g)}$: fractional coordinate of Wannier centre $n$ along `gdir`; $f_s$ is the spin factor. The total phase is $\varphi_{\text{el}} + \varphi_{\text{ion}}$ folded.

Polarization and quantum (`berry.polarizacion`, `berry.cuanto`):

$$
P_g = \varphi\,\frac{|\mathbf R_g|}{\Omega}, \qquad
\Delta P_{\text{quantum}} = m\,\frac{|\mathbf R_g|}{\Omega}, \qquad
1\ e/\text{Å}^2 = 16.02176634\ \text{C/m}^2
$$

- $\varphi$: total phase in QE units (dimensionless, quantum $m$); $\mathbf R_g$: lattice vector `gdir` (Å); $\Omega$: volume (Å³). $P_g$ is the **projection** of $P\Omega/e$ onto $\mathbf R_g$, not the modulus of $\mathbf P$.

Branch tracking (`berry.desenrollar`): $\tilde\varphi_0 = \varphi_0$, $\tilde\varphi_i = \varphi_i + m\cdot\mathrm{round}\big((\tilde\varphi_{i-1} - \varphi_i)/m\big)$; a warning is issued if any jump $|\tilde\varphi_i - \tilde\varphi_{i-1}| > 0.25\,m$ (`FRACCION_SOSPECHOSA`).

Born effective charge (`berry.analizar`), with $\mathbf u$ the total displacement (Å) and $\mathbf B_g$ the `gdir` reciprocal vector (with $2\pi$):

$$
Z^*_{g} = \frac{2\pi\,\dfrac{d\tilde\varphi}{d\lambda}}{\mathbf u\cdot\mathbf B_g}
$$

- $d\tilde\varphi/d\lambda$: slope of the linear fit of the tracked phase versus $\lambda \in [0,1]$ (`np.polyfit` degree 1 if more than 2 points; finite difference otherwise). It is the $Z^*_{g,\hat u}$ component of the tensor. If $\mathbf u\perp\mathbf B_g$, it is not computed.

Adiabatic path (`berry._interpolar_estructuras`): positions interpolated in fractional coordinates by minimum image, $f(\lambda) = f_a + \lambda\,[(f_b - f_a) - \mathrm{round}(f_b - f_a)]$, and the cell interpolated linearly.

**How Olla-DFT computes it.**

1. `qekit/cli.py: _cmd_berry` loads the polar structure, optionally `--reference` (centrosymmetric) or `--displace ATOM:dx,dy,dz` (Å, atom 1-based), and `--kperp`.
2. `berry.prepare` builds the list of structures (`--nlambda` 5 values of $\lambda$; a single point if there is no path), resolves pseudopotentials and cutoffs (`sweep.prepare_common(insulator=True)`) and, in every `pNN/`, writes `1_scf.in` and `2_berry.in` (`calculation='nscf'`, `occupations='fixed'`, `conv_thr 1e-10`, `nosym`, `noinv`, with `lberry`, `gdir` and `nppstr` inserted into `&CONTROL`), plus `correr.sh`/`correr.py`.
3. `--run` → `berry.correr`: pw.x on scf and berry at every point, skipping those that already contain `JOB DONE` unless `--redo`.
4. `--collect` → `berry.collect`: `leer_berry` extracts from `2_berry.out` `Ionic Phase`, `Electronic Phase`, `TOTAL PHASE`, `MOD_TOT`, `P = … (mod …) (e/Omega).bohr`, `direction of vector`, `Number of k-points per string`, `Number of different strings`; `valencias_de` reads the table "atomic species valence mass pseudopotential" from `1_scf.out`.
5. `berry.analizar`: unwraps the phases, converts to C/m², computes $\Delta P$ and, if the path is a displacement, $Z^*$; `comprobar_ionica` compares pw.x's ionic phase with $\sum Z_a f_a$ (warning if they differ by more than $10^{-4}$).
6. `berry.report` prints the table $\lambda$ / ionic / electronic / total / tracked / $P$; `export` writes `BERRY.dat` and `BERRY.txt`; `plot` draws $P(\lambda)$, the folded values from pw.x and a band one quantum wide.

**Where each datum comes from.**

| Datum | Source | Detail |
|---|---|---|
| Ionic, electronic, total phases, `MOD_TOT` | `pNN/2_berry.out` from pw.x (`lberry`) | `berry.leer_berry`, regular expression over the text |
| Valences $Z_a$ | `atomic species / valence` table in `1_scf.out` | `berry.valencias_de` |
| Cell, volume, $\mathbf R_g$, $\mathbf B_g$ | user structure (last one of the path) | `berry.cuanto`, `berry.analizar` |
| e/Å² → C/m² conversion | constant `berry.E_A2_A_C_M2` | 16.02176634 |
| Suspicious-jump threshold | `berry.FRACCION_SOSPECHOSA` | 0.25 of the quantum |
| Wannier centres (check) | `olla-dft wannier` | `berry.desde_wannier`, API/tests only |

**Limits and pitfalls.**

- **Insulators** only: the nscf is written with `occupations='fixed'`; in a metal the phase is undefined.
- A single point is useless: "Un solo punto. P está definida módulo el cuanto, así que este número por sí solo no significa nada".
- If a step moves the phase by more than 25 % of the quantum: "El seguimiento de la rama supone que el paso es pequeño; con saltos así, elegir la imagen más cercana es una apuesta. Sube --nlambda". If $|\Delta P| > 0.9$ quanta: "Comprueba con más puntos que no es un salto de rama disfrazado".
- Only **one component** (`--gdir`) is computed; for the vector $\mathbf P$ three runs are needed.
- If pw.x stops with "Wrong k-strings", `nosym`/`noinv` were almost certainly missing; Olla-DFT forces them, but a hand-edited input may lose them.
- In the figure, the markers "lo que escribe pw.x (plegado)" come from `berry.polarizacion_plegada`: $P = \varphi_{\text{tot}}/m \cdot \Delta P_{\text{quantum}}$ with the same `MOD_TOT` that `analizar` uses, so they coincide with the tracked branch at $\lambda = 0$ and differ from it only by integer multiples of the quantum.
- `desde_wannier` (check against Wannier centres) is not wired to the CLI; it is only used from Python or in the tests.
- No correction is applied for spin polarization or SOC (the spin factor is a fixed 2 in `desde_wannier`; pw.x handles it internally in `lberry`).

**References.**

- R. D. King-Smith, D. Vanderbilt, *Phys. Rev. B* **47**, 1651 (1993) — modern theory of polarization. DOI 10.1103/PhysRevB.47.1651.
- R. Resta, *Rev. Mod. Phys.* **66**, 899 (1994). DOI 10.1103/RevModPhys.66.899.
- N. A. Spaldin, "A beginner's guide to the modern theory of polarization", *J. Solid State Chem.* **195**, 2 (2012). DOI 10.1016/j.jssc.2012.05.010.
- D. Vanderbilt, *Berry Phases in Electronic Structure Theory* (Cambridge, 2018).

---

### `olla-dft hubbard` — Hubbard U by linear response (hp.x)

**What it answers.** How large the DFT+U parameter $U$ is for the localized ($d$ or $f$) orbitals of your system, computed by linear response with `hp.x` instead of copied from a paper, and — with `--cycle` — its self-consistent value.

**Background for non-experts.** Semilocal functionals (LDA, GGA) let an electron "see itself" (self-interaction), which over-delocalizes $d$ and $f$ orbitals and turns insulating oxides such as NiO into metals. DFT+U adds a penalty $U$ to fractional occupation of those orbitals. The value of $U$ is not a property of the element but of the system and of the *projection scheme* with which the occupations are counted. Cococcioni and de Gironcoli obtain it by measuring how the orbital occupation responds to a small perturbation of the potential: the "bare" response $\chi_0$ (without letting the rest of the system readjust) and the full one $\chi$. Their difference is the spurious curvature that $U$ must cancel. `hp.x` does that calculation with perturbation theory (DFPT) on a mesh of $\mathbf q$ vectors equivalent to a supercell. Since the $U$ obtained depends on the $U$ used in the starting scf, one must iterate until it stabilizes.

**Formulas.**

Linear response (computed by `hp.x`, not by Olla-DFT; `qekit/modules/hubbard.py`, docstring):

$$
U_I = \left(\chi_0^{-1} - \chi^{-1}\right)_{II}
$$

- $\chi_0$, $\chi$: response matrices of the occupations $n_I$ of Hubbard site $I$ to the perturbation $\alpha_J$ of the potential on site $J$, without and with self-consistent readjustment (eV⁻¹). $U_I$ in eV.

Self-consistency cycle (`hubbard.ciclo`):

$$
U^{(k+1)}_s = (1 - \mu)\,U^{(k)}_s + \mu\,U^{\text{hp}}_s\!\left[U^{(k)}\right], \qquad
\text{converged if } k \ge 1 \ \wedge\ \max_s\left|U^{\text{hp}}_s - U^{(k)}_s\right| < \text{tol}
$$

- $\mu$ = `--mixing` (1.0 by default), tol = `--tol` (0.05 eV), at most `--max-iter` = 6 rounds; $U^{(0)}_s$ = `U_SEMILLA` = $10^{-8}$ eV. The $U$ reported per species is the mean over its sites (`HubbardRun.U`).

**How Olla-DFT computes it.**

1. `qekit/cli.py: _cmd_hubbard` loads the structure; `--species` or, by default, `hubbard.elementos_hubbard` (those in the `ORBITAL_HUBBARD` table: 3d Sc–Zn, 4d Y–Cd, 5d Hf–Hg, 4f La–Lu); `--qgrid` 2×2×2; `--hubbard-style legacy|card` (the same selector as `gen`).
2. `hubbard.prepare` writes `scf.in` with `inputgen.build_pw_input(hubbard={s: U_seed}, hubbard_style=…, conv_thr=1e-15)`. With `legacy` (default, QE ≤ 7.0): `lda_plus_u = .true.`, `Hubbard_U(i) = 1e-8` and `U_projection_type = 'ortho-atomic'` inserted into `&SYSTEM` (`_fijar_proyeccion`). With `card` (QE ≥ 7.1): a `HUBBARD (<projection>)` card with `U El-orb 1e-8` and no `U_projection_type`, which is an error in those versions. `--projection` accepts `atomic`, `ortho-atomic`, `norm-atomic`, `wannier`, `pseudo`. And `hp.in` (`build_hp_input`: `nq1..3`, `conv_thr_chi = 1e-8`, `iverbosity = 2`). Fixed occupations unless `--metal`; `--nspin 2` and `--mag` are passed through.
3. `--cycle` → `hubbard.ciclo`: per iteration it creates `iterNN/`, runs pw.x (`runner.run_all`) and hp.x (`run_hp`, searched next to pw.x), reads `*.Hubbard_parameters.dat` (`collect` → `leer_parametros`, section "Hubbard U parameters", columns site/type/label/spin/new type/new label/U) and mixes.
4. `--collect` → `hubbard.collect` reads the first `*.Hubbard_parameters.dat` in the folder; `--intersite` adds `leer_v` (section "Hubbard V parameters", table atom 1 / atom 2 / distance in bohr / V) and writes `HUBBARD.card` with `tarjeta_hubbard` (`U El-orb value` and `V El-orb El-orb i j value`, with hp.x supercell indices and threshold `--v-threshold` 0.01 eV).
5. `hubbard.report` prints the table of $U$ per site, the cycle history and the warnings; `export` writes `HUBBARD_U.dat` and `HUBBARD_U.txt`, and suggests the line `olla-dft gen … --hubbard El=U`.

**Where each datum comes from.**

| Datum | Source | Detail |
|---|---|---|
| $U$ per site | `<prefix>.Hubbard_parameters.dat` from hp.x | `hubbard.leer_parametros` |
| Intersite $V$, neighbour supercell | same output, section "Hubbard V parameters" | `hubbard.leer_v` |
| Corrected orbital | table `hubbard.ORBITAL_HUBBARD` | per element; `3d` if not in the table (`2p` for the second atom of a V) |
| Seed $U$ | `hubbard.U_SEMILLA` | 1e-8 eV |
| scf `conv_thr`, `conv_thr_chi` | constants in `prepare` / `build_hp_input` | 1e-15 Ry, 1e-8 |
| $\mathbf q$ mesh | `--qgrid` | 2×2×2 (8 cells) |

**Limits and pitfalls.**

- Olla-DFT **does not compute** $\chi$ or $U$: it reads them from hp.x. Without hp.x compiled (`make hp`) the command fails: "no se encontró hp.x junto a pw.x".
- By default the scf uses the `lda_plus_u`/`Hubbard_U(i)` syntax (QE ≤ 7.0); with QE ≥ 7.1 you must request `--hubbard-style card`. The `tarjeta_hubbard` docstring warns that the card "está probado contra la sintaxis documentada, no contra una corrida de QE 7.1, porque el QE de esta máquina es 6.6".
- A single round gives "U de PRIMERA ITERACIÓN. Depende del U que llevaba el scf de partida".
- With `nq = 1×1×1`: "la perturbación ve sus propias imágenes periódicas y el U sale mal. Usa al menos 2x2x2".
- The $U$ "solo vale con la MISMA proyección"; the report repeats it in every output.
- If the cycle does not converge in `--max-iter`: "Se hicieron N vueltas sin bajar de tol eV… si el número oscila arriba y abajo, baja --mixing a 0.5; si baja despacio pero siempre en el mismo sentido, sube --max-iter"; the command returns exit code 1.
- The HUBBARD-card orbital for an element outside the table is `3d` (or `2p` as the second atom of a $V$), which may be wrong (e.g. `O-2p` is fine, `S` would get `2p`).
- The indices of the $V$ pairs are in hp.x's **supercell** numbering; the card copies them verbatim, as QE requires.

**References.**

- M. Cococcioni, S. de Gironcoli, *Phys. Rev. B* **71**, 035105 (2005) — U by linear response. DOI 10.1103/PhysRevB.71.035105.
- I. Timrov, N. Marzari, M. Cococcioni, *Phys. Rev. B* **98**, 085127 (2018) — hp.x, DFPT for U. DOI 10.1103/PhysRevB.98.085127.
- I. Timrov, N. Marzari, M. Cococcioni, *Phys. Rev. B* **103**, 045141 (2021) — self-consistent U and V, ortho-atomic. DOI 10.1103/PhysRevB.103.045141.
- V. L. Campo Jr., M. Cococcioni, *J. Phys.: Condens. Matter* **22**, 055602 (2010) — DFT+U+V.
- S. L. Dudarev et al., *Phys. Rev. B* **57**, 1505 (1998) — simplified DFT+U formulation used by QE.

---

### `olla-dft align` — band alignment between two materials

**What it answers.** Where the valence band (and the conduction band) of one material sits relative to that of the other when they are brought into contact: the *offsets* $\Delta E_v$ and $\Delta E_c$ and the heterojunction type (I nested, II staggered, III broken).

**Background for non-experts.** Every periodic calculation fixes the zero of its potential arbitrarily (the $G = 0$ term of the Hartree potential), so directly subtracting the VBMs of two different calculations gives a meaningless number. There are two ways to put them on a common scale. In **vacuum mode**, each material is computed as a slab with vacuum and its VBM is measured relative to the vacuum level of its own calculation (its ionization potential); the offset is the difference. It ignores the charge transferred on forming the contact. In **interface mode** (Van de Walle and Martin) both bulks and also the interface are computed, and the macroscopically averaged electrostatic potential on each side of the interface serves as a bridge between the two scales: it is the only term that knows about the contact dipole.

**Formulas.**

Offsets (`qekit/modules/align.py: alinear`), with $E_v^{A}$ the VBM and $V^{A}_{\text{ref}}$ the reference of calculation $A$:

$$
\Delta E_v = \left(E_v^{A} - V_{\text{ref}}^{A}\right) - \left(E_v^{B} - V_{\text{ref}}^{B}\right) + \Delta\bar V, \qquad
\Delta E_c = \left(E_c^{A} - V_{\text{ref}}^{A}\right) - \left(E_c^{B} - V_{\text{ref}}^{B}\right) + \Delta\bar V
$$

- Vacuum mode: $V_{\text{ref}}$ = vacuum level (maximum of the planar potential, mean over a 20 % window around it; `fields.work_function`), $\Delta\bar V = 0$.
- Interface mode: $V_{\text{ref}}$ = mean electrostatic potential of the bulk cell ($\langle V\rangle$ of the planar average), and $\Delta\bar V = \bar V_A - \bar V_B$ measured at the interface (`align.puente_interfaz`).
- Everything in eV; the pp.x potential (`plot_num = 11`, $V_{\text{bare}} + V_H$) comes in Ry and is multiplied by `RY_EV = 13.605693122994`.

Interface bridge (`align.puente_interfaz`): planar average $\bar V(z)$ of the cube, periodic moving macroscopic average with window $w$ (`fields.macroscopic_average`; $w$ = `--window` or $L/8$), and

$$
\bar V_A = \langle \bar{\bar V}\rangle_{z \in [L/8,\, L/4]}, \qquad \bar V_B = \langle \bar{\bar V}\rangle_{z \in [5L/8,\, 3L/4]}, \qquad \Delta\bar V = \bar V_A - \bar V_B
$$

- Material $A$ is assumed to occupy the first half of the interface cell and $B$ the second.

Alignment type (`align.alinear`), on $B$'s scale (VBM of $B$ at 0): $v_A = \Delta E_v$, $c_A = E_g^{B} + \Delta E_c$, $v_B = 0$, $c_B = E_g^{B}$:

- `=` if $|\Delta E_v| < 0.05$ and $|\Delta E_c| < 0.05$ eV (`TOL_ALINEADOS`);
- I if one gap contains the other ($v_A \le v_B \wedge c_A \ge c_B$, or the reverse);
- III if $c_A \le v_B$ or $c_B \le v_A$;
- II in every other case.

**How Olla-DFT computes it.**

1. `qekit/cli.py: _cmd_align` receives folders `a` and `b`, `--interface FOLDER` (switches on interface mode), `--axis` (c by default), `--window`, `--names`.
2. `align.leer_lado` reads the XML (`qeout.read_xml`): VBM = `<highestOccupiedLevel>`, CBM = `<lowestUnoccupiedLevel>`, $E_F$; without HOMO it fails: "no da un VBM. En un metal no hay banda de valencia que alinear; y si es un aislante, al cálculo le faltan bandas vacías (nbnd) o no usó occupations='fixed'". Without LUMO the side is flagged `es_metal` and only $\Delta E_v$ is given.
3. `align._potencial` reuses `potencial.cube` or runs `pp.x` (`fields.run_pp` with `plot_num = 11`, `output_format = 6`) and reads it with `fields.read_cube`.
4. Vacuum mode: `fields.work_function` gives `v_vacuum` and the plateau flatness; interface mode: `fields.planar_average` and its mean.
5. With `--interface`, `align.puente_interfaz` computes $\Delta\bar V$.
6. `align.alinear`, `report` (table VBM/CBM/gap relative to the reference, offsets, type and which material each carrier goes to in type II), `export` (`ALINEAMIENTO.dat`, `.txt`) and `plot` (box diagram). The box positions come from `align.posiciones_en_escala_de_b` — $v_A = \Delta E_v$, $c_A = E_g^{B} + \Delta E_c$, $v_B = 0$, $c_B = E_g^{B}$ — the same convention with which `alinear` classifies the type, so that report, export and figure cannot disagree.

**Where each datum comes from.**

| Datum | Source | Detail |
|---|---|---|
| VBM, CBM, $E_F$ | `<highestOccupiedLevel>`, `<lowestUnoccupiedLevel>`, `<fermi_energy>` in `prefix.xml` | `align.leer_lado`; requires fixed occupations |
| Electrostatic potential | `potencial.cube` from pp.x (`plot_num = 11`) | `fields.read_cube`; Ry → eV with 13.605693122994 |
| Vacuum level and flatness | maximum of the planar average, 20 % window | `fields.work_function` |
| Macroscopic window | `--window` or $L/8$ | `align.puente_interfaz` |
| "Aligned" threshold | `align.TOL_ALINEADOS` | 0.05 eV |
| Flatness threshold | constant in `alinear` | 0.05 eV |

**Limits and pitfalls.**

- Vacuum mode: the report always warns: "son las dos superficies AISLADAS. Al ponerlas en contacto se transfiere carga y aparece un dipolo de interfaz que desplaza el offset, típicamente entre 0.1 y 0.5 eV".
- If the vacuum plateau varies by more than 0.05 eV: "O falta vacío, o la losa tiene dipolo neto: usa --dipole al generarla. El nivel de vacío es la referencia de todo esto, así que ese error entra entero en el offset".
- Interface mode assumes $A$ lies in the first half of the cell and $B$ in the second, and uses two fixed windows ($[L/8, L/4]$ and $[5L/8, 3L/4]$); an asymmetric interface or layers of different thickness give a wrong bridge without warning.
- VBM/CBM are read from `highestOccupiedLevel`/`lowestUnoccupiedLevel`, which depend on the scf k-mesh; no band analysis is performed.
- If $A$ has no CBM (a metal or no empty bands), the figure draws its gap as a box of height $E_g^{A}$ (or 1 eV if there is no gap either) above $v_A$: a visual filler, not a datum.
- The offsets carry the systematic error of the functional; `TIPOS["="]` reminds that "con funcionales semilocales el error frente al experimento es de varias décimas".

**References.**

- C. G. Van de Walle, R. M. Martin, *Phys. Rev. B* **35**, 8154 (1987) — alignment via macroscopic potential. DOI 10.1103/PhysRevB.35.8154.
- A. Baldereschi, S. Baroni, R. Resta, *Phys. Rev. Lett.* **61**, 734 (1988) — macroscopic average. DOI 10.1103/PhysRevLett.61.734.
- L. Kleinman, *Phys. Rev. B* **24**, 7412 (1981) — the arbitrary zero of the potential in periodic calculations.
- J. Tersoff, *Phys. Rev. B* **30**, 4874 (1984) — alignment and interface dipoles.
