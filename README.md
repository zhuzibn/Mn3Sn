# Mn3Sn spin dynamics

MATLAB code for a two-dimensional, spatially resolved classical-spin model of Mn$_3$Sn, associated with the three publications listed below. The repository constructs a masked kagome-like array and integrates fixed-length LLGS dynamics using an NVIDIA GPU.

**Status:** the checked-in configuration is a small, undriven relaxation example. It is not a complete reproduction package for any of the three papers. The October 2026 audit found numerical and interaction-normalization problems; see [plan.md](plan.md) for evidence, priorities, and verification gates.

## Publications and model scope

1. Zhengde Xu, Xue Zhang, Yixiao Qiao, Gengchiau Liang, Shuyuan Shi, Zhifeng Zhu, Deterministic spin-orbit torque switching including the interplay between spin polarization and kagome plane in Mn3Sn, Physical Review B **109**, 134433 (2024). [DOI: 10.1103/PhysRevB.109.134433](https://doi.org/10.1103/PhysRevB.109.134433).
2. Zhengde Xu, Yue Zhou, Xue Zhang, Yixiao Qiao, Zhuo Xu, Dingfu Shao, Zhifeng Zhu, Spin-orbit torque-driven chiral domain wall motion in Mn3Sn, Applied Physics Letters **127**, 072402 (2025). [DOI: 10.1063/5.0275467](https://doi.org/10.1063/5.0275467).
3. Boyu Zhao, Zhengde Xu, Xue Zhang, Zhenhang Kong, Shouzhong Peng, Shuyuan Shi, Zhifeng Zhu, Investigation of sub-configurations reveals stable spin–orbit torque switching polarity in polycrystalline Mn3Sn, Applied Physics Letters **127**, 022407 (2025). [DOI: 10.1063/5.0273330](https://doi.org/10.1063/5.0273330).

The 2024 switching paper uses a three-sublattice macrospin model. The 2025 domain-wall paper extends it to a spatial atomistic nanoribbon and reports a different exchange coefficient. The 2025 sub-configuration paper evolves three Mn moments while scanning rotated orientations; it does not simulate coupled grains. These distinctions matter when choosing parameters and comparing results. Detailed source locators are recorded in `plan.md`.

## Requirements and running

- MATLAB with Parallel Computing Toolbox and a MATLAB-supported NVIDIA GPU.
- Statistics and Machine Learning Toolbox if `thermalenable=1`, because the code calls `normrnd`.
- Run from this repository's root. Scripts share the MATLAB workspace; `main.m` starts with `clear all`, selects GPU device 1, and closes figures.
- GNU Octave and a CPU-only simulation path are not supported by the current implementation.

From MATLAB:

```matlab
run('main.m')
```

From a shell with MATLAB on `PATH`:

```bash
matlab -batch "run('main.m')"
```

**Output:** `main.m` saves the entire workspace into `finalnew.mat` in the current directory, overwriting a file of that name. Preserve existing results before running. There is no implemented output-path option, configuration object, or restart reader. Parameters are edited directly in `main.m`; setting workspace variables before running it does not override its defaults.

The optional visualization block is disabled with `if(0)`. Enabling it requires an external `threedplott` function and a machine-specific path that are not included here. Plotting is not required for integration.

## Files and data flow

| File | Role |
|---|---|
| `main.m` | Controls, physical parameters, GPU selection, time allocation, integration, workspace save |
| `constantfile.m` | Constants; notably `hbar` is in eV s |
| `systemgeneration.m` | Initial spin construction; calls `distrib.m` |
| `distrib.m` | Occupation mask, local anisotropy axes, exchange/DMI scaling masks, optional distribution I/O |
| `integrate_llg.m` | Interaction arrays, field evaluation, GPU kernel dispatch, edge pinning, trajectory assembly |
| `atomgpurk4.m` | Four-stage per-spin update with SOT/STT and end-step normalization; fields are frozen across stages |
| `atomgpu.m` | Legacy Heun update without torque inputs |
| `atomgpupc4.m` | Legacy predictor-corrector update; currently broken and unverified |
| `dipole_calc.m` | Disabled, direct, and macrocell dipole branches |

Although `systemgeneration()` and `integrate_llg()` resemble function calls, these files are scripts that depend on shared workspace variables. The three `atomgpu*` files are functions.

## Annotated notebook

The [Mn3Sn.pdf](Mn3Sn.pdf) notebook uses the same local format as the Alpha-Fe2O3 notebook. The HTML and transcript retain 23 pages after removal of original page 23; References is now page 23.

| File | Role |
|---|---|
| [docs/Mn3Sn-notebook.html](docs/Mn3Sn-notebook.html) | Interactive notebook with original color pages, transcription, equations, page navigation, image zoom, and annotation fields |
| [docs/Mn3Sn-notebook-transcription.md](docs/Mn3Sn-notebook-transcription.md) | Editable page-by-page transcription, stable item IDs, switching tables, mask arrays, and extraction provenance |
| `docs/Mn3Sn-pages/page-01.jpg` through `page-22.jpg`, plus `page-24.jpg` | Full-page images rendered at 160 dpi; landscape pages retain their orientation. References keeps its original image filename and annotation key. |

Open the HTML file directly in a modern browser. It works offline and renders equations with native MathML. Keep the page-image directory beside the HTML and the source PDF one directory above it. Author-supplied meaning and code-correspondence notes are saved in browser storage when available; **Export annotations** downloads a JSON copy. These fields do not change the PDF or Markdown file.

The extraction preserves historical questions, differing coordinate conventions, incomplete references, and uncertain handwriting. Notebook source labels [1]–[5] refer to its own page-23 reference list. The transcription does not establish physical correctness or agreement between the notebook masks and the current MATLAB implementation.

## Interactive structure and interactions

Open either page directly in a modern browser; both work offline with their local companion files:

| Page | What it shows |
|---|---|
| [Cell and array](docs/unit-cell.html) | Notebook-based triangular schematic, outlined repeating motifs, selectable Mn/empty entries, finite-patch repetition, exact array correspondence, initial spins and easy axes |
| [Structure to interactions](docs/structure-exchange.html) | Selected receiver and neighbor reads linked across the schematic, array, coefficient table, and MATLAB source; exchange, signed DMI reads, onsite anisotropy, and both boundary flags |

The geometry has arbitrary display spacing. The 4-row × 2-column outline represents a repeating code pattern, **not a verified crystallographic unit cell**. Physical lattice vectors, Sn positions, stacking, code-to-crystal labels, and physical DMI-vector assignments remain unresolved. The drawing follows notebook pages 19–20; spin arrows come from the current initializer and do not depict relaxed states.

The pages preserve the current receiver-dependent interaction weights and edge cuts. In the 5×6 open example there are 20 occupied entries, 30 exchange connections, and 12 pairs with unequal reverse weights. `bc=0` shows legacy wrapped reads with retained masks, not validated periodic topology. The table distinguishes active reads, coefficient cuts, open edges, empty sources, and empty receivers. Its coefficients for excluded reads are raw algebraic values, not physical contributions.

Select entries in either diagram or array; choose a neighbor row to highlight its path and source expressions. Drag to pan, use the zoom slider or `+`/`−`, and press `0` to reset the focused diagram. Cell repetition regenerates a finite patch with its termination rules. Links between the two pages retain the selected array and entry.

The exchange explorer also shows [all unequal exchange pairs](docs/structure-exchange.html#unequal-exchange) together: numbered orange links correspond to a table with both source-to-receiver weights in units of `Jintra`. Node colors and labels show the receiving-site factor `hex_constant`, separately from the easy-axis class colors in the detailed view. Select a numbered link or table button to inspect that pair's individual read and source expression below. The overview follows the array and boundary controls and appears only for exchange. In the 5×6 open patch, all 12 highlighted pairs have weight `Jintra` into one endpoint and `Jintra/2` into the other; the remaining 18 connections have equal reverse weights. These are implemented field-sum coefficients, not validated reciprocal bond energies.

The shared `docs/mn3sn-data.js` snapshot comes from a CPU-only MATLAB export/check harness. It includes source excerpts and SHA-256 hashes, masks, axes, initial spins, raw coefficients, and independently measured impulse responses. After relevant source changes, refresh it from the repository root:

```bash
matlab -batch "addpath('docs'); export_visualization_data('docs/mn3sn-data.js')"
node docs/check_visualization_data.cjs
node --check docs/mn3sn-model.js
node --check docs/mn3sn-view.js
```

The exporter executes reviewed setup, shift, and static-field excerpts, with wall initialization and distribution I/O disabled. It does not call `main.m`, launch integration, select a GPU, or write simulation results. The command overwrites only the generated viewer snapshot. It checks ten sizes (5×5 plus all combinations of 5/9/13 rows and 6/10/14 columns), both boundary flags, and every occupied Cartesian spin impulse. It currently requires `DMIenable=1`, positive nonzero exchange/DMI/anisotropy coefficients, and the existing `musRE=3*mub` mapping; changes to these assumptions require reviewing the export scope. Node checks the display helpers against the measured operators and rejects stale source hashes. These tools need MATLAB and Node respectively; viewing the pages needs only a browser.

Validation on 2026-10-21: MATLAB R2024b Update 1 static checks passed; Node matched 367,700 exchange/DMI operator entries and all source hashes. Chromium 151 checks exercised offline loading, selections, cross-page links, controls, and desktop/390 px layouts. These are visualization/source-consistency checks, not scientific validation of the existing interaction model. No dynamics or publication benchmark was run. Exact procedures and remaining evidence gaps are recorded in [plan.md](plan.md#visualization-task--2026-10-21).

## Current conventions and defaults

Spin components are arrays indexed `(W, L, time)`. Array row/column indices are not a documented crystallographic coordinate system. The spin initialization and easy axes lie in the spin-space $yz$ plane; the default polarization is $+x$, perpendicular to that plane.

For the notebook's 5×6 example with open boundaries (`bc=1`), the code and [Mn3Sn.pdf, pages 19–20](Mn3Sn.pdf#page=19) are structurally equivalent in their occupation pattern and nearest-neighbor connectivity: 20 occupied sites, 10 empty slots, and 30 undirected connections, including the same removal of site `(3,6)`. This equivalence concerns the lattice connectivity; exchange weights and physical crystallographic coordinates remain separate issues. See the [A02 follow-up in plan.md](plan.md#a02-follow-up-page-19-structure-comparison-2026-10-21) for the recorded verification.

Ignoring the scalar prefactor, the implemented anisotropy matches the methodology in [Mn3Sn.pdf, page 4](Mn3Sn.pdf#page=4): onsite uniaxial easy-axis anisotropy with field proportional to $(\mathbf e_i\cdot\mathbf s_i)\mathbf e_i$. The revised page explicitly writes $(\sin\theta\cos\varphi,\sin\theta\sin\varphi,\cos\theta)$ with $\varphi=\pi/2$ and $\theta=\pi/6,3\pi/2,5\pi/6$, preserving the same three $yz$ easy-axis vectors. Their directions and site assignment agree with pages 4 and 21; a MATLAB check of the page-21 5×6 example matched all 20 occupied sites and 10 empty slots. The prefactor difference remains a separate normalization issue. See the [A04 follow-up in plan.md](plan.md#a04-follow-up-anisotropy-correspondence-to-pages-4-and-21-2026-10-21) for source locators and verification.

`atomtype_==1` denotes an occupied Mn spin site, and `atomtype_==3` denotes an empty slot. The `RE`/`TM` names and comments are inherited from ferrimagnetic code: occupied Mn sites use the `RE` coefficients. They do not describe two physical Mn species. Empty slots must remain zero vectors; unit-norm checks apply only to occupied sites.

| Quantity | Checked-in value | Meaning or caveat |
|---|---|---|
| Array dimensions | `natomW=5`, `natomL=5` | 25 slots, 16 occupied spins, 9 empty slots |
| Boundaries | `bc=1` | Open neighbor shifts; `bc=0` requests periodic shifts but interaction masks still remove seam bonds |
| Solver | `rk4=1` | Default path exercised by the audit; full coupled RK4 order is not established |
| Exchange | `Jintra=2.8083e-21 J` | Approximately 17.53 meV before site-dependent scaling |
| Anisotropy | `Ksim=0.196 meV` | Field uses `Ksim/mu`; mapping to the published energy coefficient needs resolution |
| DMI | `Dsim=0.833 meV` | Enabled, with site-dependent scaling and directed stencil signs |
| Moment | `musRE=3*mub` | Per occupied Mn site |
| Damping | `alp=0.003` | Dimensionless |
| Spin Hall angle | `thetaSH=0.06` | Dimensionless |
| Length/thickness | `d=tz=0.4 nm` | Not the 30 nm thickness of the switching papers |
| Volume and magnetization | `volume=0.3778 nm^3`, `msRE=musRE/volume` | Switching papers instead use `Ms=6*moment/Vcell` |
| Drives | `jcSOT=jcSTT=0`, `Hext=[0,0,0]` | No time-dependent pulse schedule |
| Thermal noise | Off; `T=100 K` | Temperature has no effect while thermal noise is off |
| Dipoles / edge pinning / wall initialization | All off | Optional branches require separate validation |
| Timestep | `tstep=5e-16 s` | 0.5 fs; switching papers report 5 fs |
| Output blocks | 10 blocks, `gpusave=1 ps` | Nominal runtime 10 ps |
| Sampling | `savetstep=200` | Nominal 100 fs sampling; block boundaries change actual spacing |

Field sums are flux densities in tesla despite variable names such as `Hext`, `hex_x`, and `hhx`. Moments are J/T, link energies are J, and `gamRE` is rad/(s T). `hbar` in eV s is multiplied by `ele` when calculating gamma; in the torque formula, the numerical conversion from eV to J cancels the elementary-charge denominator. Do not insert another charge factor without deriving the units.

## Output and known limitations

`mmx`, `mmy`, and `mmz` contain the sampled spin components. At the defaults their shapes are `5 x 5 x 100`. `atomtype_`, initial arrays, constants, intermediate fields, and other workspace variables are also saved. No versioned result schema or explicit code/hardware provenance is provided.

The saved `t` vector is not a reliable physical clock. Each block allocates 2000 states but advances only 1999 steps. The default total is therefore 19,990 updates (9.995 ps), while the initial state is labeled 0.5 fs. Saved timestamps also drift relative to the block states. Do not infer switching times or domain-wall velocities from the current labels without reconstructing the actual update count.

The most consequential audit findings are:

- Coupled exchange, DMI, and anisotropy fields are evaluated once per step and reused at every RK4 stage.
- Site-dependent exchange weights are not reciprocal on some default-lattice bonds, preventing identification with a conventional symmetric pair Hamiltonian.
- Published anisotropy, exchange, magnetization, and thickness conventions have not been reconciled with the code.
- Heun generates NaNs at empty slots. Predictor-corrector startup passes an obsolete argument list. Torque flags and alternative solver guards are incomplete.
- Dipole mode 1 returns different variable names from those consumed by integration; mode 2 is explicitly marked suspect; mode 3 uses the $y$ spin component for the $z$ moment.
- `dwcalc=1` initializes occupied sites with collinear domain states and initializes empty slots as spins. It does not construct the published noncollinear chiral wall.
- `loadstartm` defines a filename but never loads it. Loading a fixed distribution occurs after derived masks are built.
- No octupole observable, orientation scan, wall-position measurement, pulse protocol, accepted baseline, automated tests, or CI is included.

## Audit verification

The [2026-11-02 three-spin comparison](docs/three-spin-comparison.md) reduces the array to three occupied sites and checks it against `Mn3Sn_3atoms`. Source-based field and LLGS checks pass after explicit parameter matching. Across 24 Figure 5 cases, a diagnostic coupled RK4 variant reproduces the corrected reference and shows fourth-order timestep convergence. Frozen fields remain inaccurate: even with matched parameters at 0.5 fs, the 275 Oe case ends near 90° instead of 30°. The report includes plots, 342 full-length diagnostic runs, short MATLAB cross-checks, and reproduction commands under `verification/`. Production source is unchanged; this does not resolve the larger-lattice audit findings.

On 2026-10-14, MATLAB R2024b Update 1 and an NVIDIA GeForce RTX 3060 Laptop GPU were available. A temporary two-block GPU probe used six allocated states per block and saved every two indices. Its six saved frames were finite, occupied spins had zero measured norm error, and empty slots remained zero. It confirmed 10 actual updates against 12 nominal updates. CPU probes reproduced the Heun empty-slot failure, obsolete RK4 startup call, and missing dipole-mode-1 field outputs. An effective exchange matrix reconstructed from the current stencils had 24 directed mismatches (12 unordered bonds), with maximum normalized asymmetry 0.5. MATLAB Code Analyzer returned 18 messages across the nine source files; these are diagnostics, not a physics-validation result.

The production-length simulation, convergence order, stochastic dynamics, energy derivatives, and published results were not validated. The audit wrote no simulation results and changed no MATLAB source. [plan.md](plan.md) records the staged improvement and reproduction work.

## License

See [LICENSE](LICENSE), which contains the GNU General Public License, version 3.
