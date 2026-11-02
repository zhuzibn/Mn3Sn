# Mn3Sn audit and improvement plan

The historical three-spin comparison is recorded in [the diagnostic plan below](#three-spin-comparison--2026-11-02) and the [results report](docs/three-spin-comparison.md). The subsequent [frozen-field repair](docs/frozen-field-repair.md) is followed by the approved reference-parameter/protocol alignment in the final section and [its report](docs/reference-protocols.md). Earlier audit statements describe the checkout at their recorded stage.

Audit date: 2026-10-14. Base revision: `2307530`, branch `main`; working tree was clean before documentation work. Status: audit and documentation complete; implementation stages below are proposed and have not been approved or executed. No MATLAB source or physical defaults were changed.

## Goal and Success Criteria

The requested outcome is a code audit and root `AGENTS.md`, `README.md`, and improvement plan. Acceptance: describe actual behavior and commands, provide source-located findings with evidence levels, separate the three publication scopes, and define ordered improvement gates. This task does not implement solver fixes, change the material model, install dependencies, or reproduce the papers.

Future implementation should yield a traceable Mn3Sn model with verified topology, energy/field consistency, LLGS dynamics, physical timestamps, safe results, and separately validated publication protocols. Preserve the legacy configuration as a documented comparison; do not silently relabel corrected results as the old baseline.

## Model Classification and Equation Family

Class 4: material-specific spatial atomistic model intended to represent a Mn kagome layer. This is a classification of the code's intended scope, not certification of a crystallographic construction. The code has independent spins at occupied sites and three local easy-axis classes. The published spatial model resolves Mn moments [2, p. 072402-2 and Fig. 1]; the switching models evolve three sublattice moments [1, Sec. II; 3, p. 022407-3].

Equation family: fixed-length LLGS, including precession, Gilbert damping, and SOT; the code also contains optional STT and noise branches. State: $(s_{ix},s_{iy},s_{iz})$ on occupied sites, with $|\mathbf s_i|=1$ after each kernel update; empty slots have $\mathbf s_i=0$. Inputs are workspace parameters, occupation/mask arrays, initial spins, and numerical controls; outputs are sampled component arrays and a nominal time vector. Electronic structure and a quantitative Hall-response model are outside this implementation.

## Existing Behavior and Reference Implementations

All nine MATLAB files were read. There were no existing project instructions, README, plan, tests, benchmarks, changelog, or CI files. The default calculation has 16 occupied sites in a 5×5 mask, no current or external field, no thermal noise, and no domain wall. It is neither the coherent three-sublattice switching protocol nor a long spatial nanoribbon protocol.

### Findings ranked by consequence

Severity reflects impact if the affected path is used. “Confirmed” means source inspection plus a focused executed probe where stated; “scientific gap” means a published convention or required physical mapping has not been reconciled.

| ID / priority | Finding and source locator | Evidence and consequence | Proposed action |
|---|---|---|---|
| A01 / P0 | Frozen coupled fields: `integrate_llg.m:191–247`; `atomgpurk4.m:7,98,190,282` | Confirmed statically. Every RK4 stage receives the same fields, although exchange, DMI, and anisotropy depend on all spins. It is RK4 for a local frozen-field problem, not demonstrated fourth-order integration of the coupled equations. Heun has the same field-freezing limitation. | Recompute conservative fields from the complete stage state; verify coupled convergence before claiming order. |
| A02 / P0 | Nonreciprocal exchange: `distrib.m:48–57`; `integrate_llg.m:19–51,191–199` | Executed effective-weight reconstruction found 24 directed mismatches / 12 unordered bonds on the default active lattice, maximum asymmetry 0.5 after dividing by `Jintra`. Equal Mn moments make this incompatible with the usual symmetric pair-energy derivative. Site-wise half weighting cannot simply be explained by counting both bond directions. | Resolve the intended bond normalization with sources; represent and verify symmetric physical bonds. |
| A03 / P0 | Time contract: `main.m:120–128`; `integrate_llg.m:132–143,291–306` | Confirmed by two-block probe: 10 updates versus 12 nominal. Each block allocates states rather than the required steps-plus-one; initial state is labeled `tstep`. At defaults, 19,990 updates give 9.995 ps instead of 10 ps. Decimated times drift at block boundaries. | Count updates explicitly and derive saved times from those counts; define initial/final-state policy. |
| A04 / P0 | Anisotropy normalization: `main.m:61`; `integrate_llg.m:203–205` | Scientific gap. Code field is $K_{sim}(\mathbf e_i\cdot\mathbf s_i)\mathbf e_i/\mu_i$, implying energy $-K_{sim}(\mathbf e_i\cdot\mathbf s_i)^2/2$. Published single-site energy is written without this half factor [1, Sec. II; 2, p. 072402-2; 3, p. 022407-3]. | Resolve scalar/vector $K$ notation and compare a derived field with energy finite differences; do not double the code coefficient without approval. |
| A05 / P0 | SOT scale and publication protocol: `main.m:56,72–75,84–101,120–128` | Confirmed mismatch against [1, Sec. II and Fig. 2; 3, p. 022407-3]. Code uses $M_s=\mu/V$, thickness 0.4 nm, and 0.5 fs; switching papers use $6\mu/V$, 30 nm, and 5 fs. For equal nonzero current and other factors, the code torque amplitude would be 450 times the paper formula: $(6\mu/V)(30\,nm)/[(\mu/V)(0.4\,nm)]$. No pulse schedule exists. Current defaults are zero, so this does not drive the audit smoke. | Separate paper-specific configurations after units and coefficients are agreed. |
| A06 / P1 | Periodic mode incomplete: `integrate_llg.m:19–51,177–189`; `distrib.m:21–25,49–57,66` | Confirmed statically: exchange masks remove boundary bonds independent of `bc`; the occupation and DMI masks also contain termination rules. `circshift` alone does not make this a periodic kagome supercell. | Define a compatible periodic cell and separately test all exchange/DMI seams and open terminations. |
| A07 / P1 | Empty-slot handling differs by solver: `atomgpu.m:50–53`; `atomgpurk4.m:370–373`; `integrate_llg.m:248–259` | CPU Heun zero-spin probe returned nonfinite components. All kernels normalize zero spins without a guard, but only RK4 dispatch resets empty outputs to zero. The short default RK4 run was finite after that reset. | Exclude empty slots from physical kernels or guard normalization consistently; reject unvalidated solvers first. |
| A08 / P1 | Predictor-corrector: `integrate_llg.m:237–238`; `atomgpupc4.m:27–35,46–54,83–85,89–96` | Obsolete nine-input RK4 startup call reproduced an argument failure. Additional static defects: AB4 predictor uses `-3` rather than `-9` for oldest derivative; historical and predicted damping terms reuse the current spin; historical fields are not retained. | Disable selection clearly; restoration is a separate numerical task, not a signature-only repair. |
| A09 / P1 | Torque flags/guards: `main.m:95–116,139–141`; `integrate_llg.m:71–72,257–259` | FL amplitudes are derived from already DL-gated amplitudes; `chi`/`chiSTT` are zero in both branches. Non-RK4 guard covers SOT flags only, while Heun/PC accept no STT inputs. FL-only and STT-only selections can fail to represent the selected drive. | Verify independent flags, signed amplitudes, and dispatch rejection for both drives; scope STT separately from cited SOT papers. |
| A10 / P1 | Domain-wall initialization: `systemgeneration.m:6–30` | Occupied sites all have `atomtype_==1`, so legacy branch assigns the same nearly collinear spin on each side; empty slots enter the alternate species branch and acquire nonzero spins. This is not the six-state noncollinear CDW in [2, Fig. 1]. | Build a sourced, sublattice-resolved wall and relax it before applying SOT; preserve vacancies. |
| A11 / P1 | Dipole branches: `dipole_calc.m:8–34,35–50,77–79`; `integrate_llg.m:108–123,229–231` | Mode 1 probe created `cpuhdipole*` but no consumed `hdipole*_`. Mode 2 is explicitly marked suspect and uses a different coordinate assignment. Mode 3 calculates `mui_sz` from `mmytmp`; vacancy dummy moments also corrupt macrocell centers. Dipole distances use a square index grid rather than a sourced kagome coordinate mapping. | Keep off for the current scoped model; reject or separately repair each mode with physical coordinates and dipole reference tests. |
| A12 / P1 | Restart/fixed mask: `main.m:15–18`; `distrib.m:125–130` | `loadstartm` never loads `startmname`. Fixed occupation loads after anisotropy and DMI masks are derived, without validation; old RE/TM masks are incompatible with the Mn encoding. | Validate loaded Mn masks first, derive all dependent arrays afterward; implement restart only with agreed state/time schema. |
| A13 / P1 | Edge pinning: `integrate_llg.m:264–271` | Right slice `natomL-fixededgeL:end` pins one more column than the left slice. Both assignments overwrite vacancy zeros and impose a single spin on all sublattices. | Keep off; if requested, pin valid Mn sites with sourced sublattice-resolved states and symmetric widths. |
| A14 / P1 | Thermal contract: `main.m:31,118`; `integrate_llg.m:216–222` | Noise is sampled once per outer step, then held across frozen-field stages; seed uses `rng('shuffle')`. No stochastic interpretation, invariant measure, convergence, or thermal benchmark is specified. | Maintain noise-off deterministic work first; separately approve and validate stochastic LLGS and reproducible seeds. |
| A15 / P2 | Performance/output: `integrate_llg.m:132–134,233,248–255,274–289`; `main.m:174` | CPU block arrays and field calculations plus repeated GPU transfers and scalar vacancy loops prevent a fully GPU-resident pipeline. Multistep history is maintained for RK4 too. Full-workspace save overwrites a fixed filename; no observables or provenance schema exists. | After correctness gates, profile, keep working states/fields on the GPU, sample deliberately, and save explicit results safely. |

### A02 follow-up: page-19 structure comparison (2026-10-21)

**Structure equivalence: verified for the notebook's 5×6 open-boundary example.** The code's occupation pattern and nearest-neighbor connectivity match [Mn3Sn.pdf, page 19](Mn3Sn.pdf#page=19), including removal of site `(3,6)`. [Page 20's explicit occupation table](Mn3Sn.pdf#page=20) confirms the same mask: 20 occupied sites and 10 empty slots. The resulting graph has 30 undirected connections and zero connectivity mismatches in the recorded probe below. Source mapping: `distrib.m:11–32` defines occupation; `integrate_llg.m:19–51,177–199` defines the exchange masks, open-boundary shifts, and neighbor reads. This establishes equivalence of the illustrated site pattern and connectivity, without claiming verified physical distances, lattice vectors, or crystallographic site labels.

**Exchange-weight status: open; recorded for a future fix.** Structural equivalence does not validate the exchange weights; the A02 discrepancy below remains separate from the confirmed notebook/code topology match.

Historical context: [A02-Nonreciprocal_exchange.md](A02-Nonreciprocal_exchange.md) summarizes Zhengde's discussion and the March–April 2024 weekly-report evidence, including system sizes, exchange scaling, reported relaxation outcomes, and unresolved timestep/DMI effects. These records explain the reported motivation for the correction without validating its physical necessity or resolving A02.

Published-code comparison (2026-10-29): the same note now records an isolated verification of the PRB 2024 `code/fig5` three-sublattice model. Its active exchange is reciprocal, all six zero-drive states are stable minima, and the sampled Figure 5 destinations are reproduced. A replay matches a fresh short MATLAB probe and the archived default trajectory. The comparison also identifies an incorrect half-step fourth RK stage, unsigned angle postprocessing, and a small time-label drift; these do not authorize source repairs or resolve the spatial A02 finding. Evidence, commands, limitations, and audit artifacts are linked in the note.

The receiving-site scaling in `distrib.m:48–57`, applied in `integrate_llg.m:191–199`, makes some exchange bonds nonreciprocal. For example, with $J=\texttt{Jintra}$, the bond `(1,1)` ↔ `(2,2)` contributes with weight $J$ at `(1,1)` and $J/2$ at `(2,2)`. The 5×6 probe found **12 unordered nonreciprocal bond pairs** (24 directed mismatches), with maximum asymmetry $0.5J$.

Verification used MATLAB R2024b Update 1 on the CPU, with `natomW=5`, `natomL=6`, `bc=1`, and fixed-distribution loading/saving disabled. The probe ran `systemgeneration.m`, evaluated the exchange-mask setup and exchange-field expressions from `integrate_llg.m`, and reconstructed the coupling matrix by exciting each occupied site's x component in turn, using `Jintra=1` and unit moments to expose normalized weights. It found 20 occupied sites, 30 undirected bonds, and zero connectivity mismatches against the page-19 neighbor rules. No time integration or simulation-result writes were performed.

**Evidence limit:** page 19 does not establish easy-axis angles, initial spin directions, or DMI signs. Easy-axis correspondence to notebook pages 4 and 21 is established separately in the A04 follow-up below; initial spin and DMI mappings remain unresolved by this structure comparison.

**Future fix gate:** resolve the intended bond normalization from sources, derive reciprocal bond weights from the agreed Hamiltonian, and verify both reciprocity and energy-gradient consistency on the 5×6 example and the default 5×5 patch. Do not change the site factors solely to force a passing symmetry check. No MATLAB source changes are authorized by this note.

### A04 follow-up: anisotropy correspondence to pages 4 and 21 (2026-10-21)

**Conclusion: ignoring the scalar prefactor, the code and notebook use the same anisotropy methodology.** [Mn3Sn.pdf, page 4](Mn3Sn.pdf#page=4) and `integrate_llg.m:203–205` describe onsite uniaxial easy-axis anisotropy: the energy depends on $-(\mathbf e_i\cdot\mathbf s_i)^2$, up to an additive constant and a positive coefficient, and the field is proportional to $(\mathbf e_i\cdot\mathbf s_i)\mathbf e_i$. The easy-axis directions and their site assignment also agree with notebook pages 4 and 21. This comparison sets aside the prefactor difference at the user's request; A04's normalization issue remains a separate finding, and no coefficient or MATLAB source is changed.

The revised page 4 explicitly writes the component formula $(\sin\theta\cos\varphi,\sin\theta\sin\varphi,\cos\theta)$ and annotates all three Configuration II axes with $\varphi=\pi/2$, with $\theta=\pi/6,3\pi/2,5\pi/6$ for $\mathbf k_a,\mathbf k_b,\mathbf k_c$, respectively. Thus $\mathbf e_i=(0,\sin\theta_i,\cos\theta_i)$ as before; the notebook's full-turn in-plane value $3\pi/2$ gives $(0,-1,0)$. This coordinate angle is distinct from the spin-to-easy-axis angle in the energy derivation. The vectors equal the axes in `distrib.m:84–115`, where $9\pi/6=3\pi/2$. For the 5×6 example, the row/column assignments match [page 21](Mn3Sn.pdf#page=21) exactly (`X` means an empty slot):

```text
3π/2   X    3π/2   X    3π/2   X
5π/6  π/6   5π/6  π/6   5π/6  π/6
 X    3π/2   X    3π/2   X     X
π/6   5π/6  π/6   5π/6  π/6   5π/6
3π/2   X    3π/2   X    3π/2   X
```

**Verification:** directly inspected the original page-4 and page-21 images and their transcription items `P04-ANISOTROPY`, `P04-CONFIGURATION-II-AXES`, and `P21-ANISOTROPY-EXCHANGE-MASKS`. A CPU-only `matlab -batch` probe set `natomW=5`, `natomL=6`, and both fixed-distribution I/O flags to zero, then ran `constantfile.m` and `distrib.m`. It compared the occupation mask and `coeffs_2`, `coeffs_3` against the table above: all 20 occupied sites and 10 empty slots matched, with component differences below $10^{-14}$. It also evaluated the actual anisotropy expressions from `integrate_llg.m:203–205` for occupied spins $(0.3,0.4,\sqrt{0.75})$, confirming the same projection/axis dependence after accounting for the scalar prefactor (zero measured residual). MATLAB exited successfully. No GPU, time integration, or simulation-result writes were involved. The default 5×5 configuration was not substituted for the notebook's 5×6 example.

The field comparison uses the per-atom convention in `presentations/equations/sections/hamiltonian.tex`, Eq. `eq:atomistic-effective-field`, $\mathbf B_i=-\mu_i^{-1}\partial E/\partial\mathbf s_i$. Agreement here is with the local notebook's anisotropy construction; it does not resolve the mapping to published crystallographic site labels or validate other interaction terms.

**Page-4 documentation refresh (2026-10-21):** visually compared the revised PDF page against the previous image, refreshed `docs/Mn3Sn-pages/page-04.jpg`, and synchronized the Markdown/HTML transcription, source hash, and README. `pdftoppm -r 160 -jpeg -jpegopt quality=90 Mn3Sn.pdf /tmp/mn3sn-page4-rendered/page` produced byte-identical images for the other 23 pages. `/home/zhuzibn/.venvs/webpage-capture/bin/python /tmp/check_mn3sn_page4_refresh.py` checked the three spherical-coordinate vectors (component errors below $10^{-14}$), PDF hash/size, MathML, angle table, page links, navigation, and 1600 px/390 px layouts in Chromium 151.0.7922.34; final checks passed with no page errors or network requests. Chromium required approved escalation after the sandbox blocked launch; an initial check counted one extra standalone angle entry and was corrected. The MATLAB checks above are prior evidence and were not rerun for this documentation refresh. The Class 4 fixed-length LLGS model, source files, and remaining physical-reference gaps are unchanged.

### Reference implementations and equation conventions inspected

Read-only ASTRO reference: `/mnt/c/Users/zzf-m/OneDrive/code_softwares/spin_dynamics/atomistic_model/2D/atomistic/atomistic_v2/{AGENTS.md,README.md,field_calc.m,integrate_llg.m}`. It demonstrates explicit solver rejection, separate deterministic fields, stage evaluation, step-count handling, and a documented result schema. Its Gd-FeCo geometry and material values are not Mn3Sn evidence. Its accepted baselines are not baselines for this repository. The macrospin implementation was not inspected because this audit concerns the spatial checkout rather than a cross-model code migration.

Equation files consulted under `/mnt/c/Users/zzf-m/OneDrive/presentations/equations/`: `spintronic_methods.tex`, `sections/hamiltonian.tex` (labels `eq:atomistic-hamiltonian`, `eq:atomistic-effective-field`, `eq:atomistic-pair-fields`), `sections/llg_equation.tex` (`eq:llg-gilbert`, `eq:ll-with-flt`), `sections/atomistic_boundaries.tex`, and the Mn3Sn subsection of `sections/material_parameters.tex` (`tab:mn3sn-parameters`). They informed derivative units, pair counting, LL/Gilbert conversion, and the separation of model-specific parameter sets. The parent `presentations/AGENTS.md` and `README.md` were absent. The equation notes are conventions/navigation, not qualifying evidence for physical topology.

The paper index and relevant switching topic under `_knowledge` were inspected; all three indexed summaries were read after user authorization, and the model pages of their PDFs were directly inspected. A targeted presentation-knowledge catalog/index search returned no match for Mn3Sn or these article numbers. No network access or supplemental-PDF inspection was performed.

## State Variables, Geometry, and Boundaries

Current storage is `(W,L)` with an occupation mask generated by row modulo 4 and column parity. `natom=natomW*natomL` counts slots, not Mn atoms. Local spin-space axes are $(0,\sin\theta,\cos\theta)$, with $\theta=\pi/6,5\pi/6,9\pi/6$. Default spin angles are 121°, 241°, and 361° under another parity mapping. These are software facts from `distrib.m` and `systemgeneration.m`; the exact relationship to the published site numbering and crystallographic coordinates remains unresolved.

### Claim-to-implementation structure table

| Physical detail | Current value or proposed rule | Builder/interaction mapping | Qualifying reference and exact locator | Evidence status |
|---|---|---|---|---|
| Magnetic layer and species | Mn spin sites in one modeled kagome layer; Sn does not carry a modeled spin | `atomtype_==1`; mask and three easy-axis classes | [1, Sec. II, Fig. 1(a)]; [2, p. 072402-2, Fig. 1(a)]; [3, p. 022407-3] | General model supported; 5×6 occupation mask verified against notebook page 20 (A02 follow-up); exact mapping to published crystallographic sites remains unresolved |
| Lattice vectors, basis coordinates, physical axis transform | Code has only `d=0.4 nm` and array-index geometry | `distrib.m`; square coordinates in dipole branches | Figures in [1–3] show topology but do not substantiate these numerical coordinate rules | **REFERENCE MISSING** for current metric construction; blocks physical length/velocity and dipole mapping |
| Sublattice/site map and easy axes | Three parity-dependent axes in the $yz$ spin plane | `coeffs_2`, `coeffs_3`; initial spin parity rules | [1, Sec. II and Fig. 1(a)], [2, Fig. 1(a)], [3, Fig. 1] identify nearest-Sn easy axes | Easy-axis directions and 5×6 site assignment verified against notebook pages 4 and 21 (A04 follow-up); **REFERENCE MISSING** for exact parity-to-publication labels and rotation |
| Layer stacking/interlayer bonds | No interlayer bonds or layer coordinates | One 2D spin array | [1, Sec. II] restricts interactions to same-layer nearest neighbors; [2, Fig. 1(a)] shows the modeled plane | Supported model omission for [1]; no claim of complete bulk structure |
| Exchange shell | Eight possible index shifts selected by parity masks, intended same-layer nearest neighbors | `Asim*` arrays in `integrate_llg.m` | [1, Sec. II]; [2, Fig. 1(a)]; [3, p. 022407-3] | 5×6 open-boundary connectivity verified against notebook page 19 (A02 follow-up); **REFERENCE MISSING** for physical distances, crystallographic shell assignments, periodic seam mapping, and interaction normalization |
| Exchange coefficient/sign | Positive `Jintra≈17.53 meV`, then 0.5 or 1 site scaling | `hex_constant.*Asim*` | [1, Sec. II] and [3, p. 022407-3]: 17.53 meV; [2, p. 072402-2]: 9 meV | Published model-specific coefficients supported; **REFERENCE MISSING** for conversion to each code bond and the edge/interior weights |
| DMI direction/sign/multiplicity | $B_x^{DMI}=0$; signed $yz$ differences, multiplied by 0.5 or 1 | `hdmi_constant*`, `hdmi_y`, `hdmi_z` | [1, Sec. II]; [2, p. 072402-2] states positive D favors clockwise sublattice order | **REFERENCE MISSING** for exact directed-bond table, orientation convention, edge normalization, and paper-sum multiplicity |
| Open termination / periodic cell | Default open shifts; terminal-column sites removed; periodic option partly truncated | `distrib.m:14–29`; stencil boundary masks | [2, Fig. 1 and ribbon discussion] supports a spatial ribbon, not these exact edge rules | **REFERENCE MISSING** for physical termination and compatible periodic supercell; blocks topology replacement |
| Array ordering and empty-slot code | `(W,L,time)`, occupied `1`, empty `3` | All script arrays | Software bookkeeping; no scientific citation required | Preserve unless explicit migration is approved |

Development sequence: verify a sourced, compatible periodic bulk cell, then separately validate open or mixed boundaries for the finite-ribbon target. This proposed verification sequence does not replace the current open default or establish the paper's exact boundary conditions. Validate every seam bond, retained edge bond, and finite-size response; edge pinning is a distinct physical constraint.

## Hamiltonian

The published models include positive-sign antiferromagnetic exchange, a positive-sign DMI energy, easy-axis anisotropy, and Zeeman energy [1, Sec. II; 2, p. 072402-2; 3, p. 022407-3]. Their displayed pair sums do not by themselves settle the numerical mapping to this code's site-scaled stencil. Preserve the reported coefficients as source values until pair multiplicity is resolved.

A proposed unambiguous unique-bond specification, subject to that reconciliation, is

$$
E=\sum_{\langle ij\rangle}J_{ij}\mathbf s_i\cdot\mathbf s_j
+\sum_{\langle ij\rangle}\mathbf D_{ij}\cdot(\mathbf s_i\times\mathbf s_j)
-\sum_i K_i(\mathbf s_i\cdot\mathbf e_i)^2
-\sum_i\mu_i\mathbf s_i\cdot\mathbf B_{ext}.
$$

Each physical bond is counted once, $J_{ji}=J_{ij}$ and $\mathbf D_{ji}=-\mathbf D_{ij}$. In an ordered-pair implementation both pair sums instead carry $1/2$. This is a bookkeeping convention; **TO BE RESOLVED**: derive the effective $J_{ij}$ and $\mathbf D_{ij}$ from the papers' reported sums before adoption. The general energy/field relation is supported by [1–3]; exact code shell assignment remains source-gated above. Energy is total J, $K_i$ is J/site, and $J_{ij}$ and DMI components are J/bond.

The current code has no total-energy evaluator. Because of A02, a conventional pair energy reproducing all current exchange fields cannot be assumed. Its anisotropy field corresponds to $K_i=K_{sim}/2$ under the above scalar convention. This exposes the discrepancy rather than choosing a correction.

Dipole energy: not applicable to the initial approved-scope proposal because dipoles remain off. Restoration requires separate physical coordinates and validation. SOT, STT, and thermal noise are not Hamiltonian terms.

## Parameters and Evidence

Values below record current behavior and source comparisons, not a new adopted parameter set.

| Parameter | Current value and units | Published comparison / evidence | Status and model entry |
|---|---|---|---|
| `ele`, `mu_0`, `mub`, `hbar`, `kb` | $1.602\times10^{-19}$ C; $1.25663706143592\times10^{-6}$ T m/A; $9.274\times10^{-24}$ J/T; $6.58211951440\times10^{-16}$ eV s; $1.38064852\times10^{-23}$ J/K | `constantfile.m`; physical constants, not fitted material inputs | Preserve rounded legacy values; consistency matters more than incidental constant updates |
| `Jintra` | $2.8083\times10^{-21}$ J, ≈17.53 meV | [1, Sec. II; 3, p. 022407-3]; [2, p. 072402-2] instead reports 9 meV | Bond normalization **TO BE RESOLVED**; exchange fields |
| `Ksim`, `Dsim` | 0.196 and 0.833 meV, converted to J | [1, Sec. II; 2, p. 072402-2; 3, p. 022407-3] | Values supported; field prefactors and stencil assignments unresolved |
| `musRE`, `gRE`, `gamRE` | $3\mu_B$, 2, $g\mu_B/(\hbar_{eV}e)$ ≈$1.759\times10^{11}$ rad/(s T) | Moment: [1–3, model sections]; exact g value not substantiated in inspected model pages | Moment supported; g provenance **REFERENCE MISSING**; occupied-site LLGS |
| `alp`, `thetaSH` | 0.003, 0.06 | [1–3, model sections] | Supported model inputs; damping and SOT |
| `d`, `tz` | Both 0.4 nm | [1, Sec. II; 3, p. 022407-3] give film thickness 30 nm; [2] has a distinct ribbon model | 0.4 nm metric/thickness **REFERENCE MISSING**; separate site spacing from magnetic thickness |
| `volume`, `msRE` | 0.3778 nm³; $3\mu_B/V≈7.364\times10^4$ A/m | [1, Sec. II; 3, p. 022407-3]: $V_{cell}=0.3778$ nm³, $M_s=6\mu/V_{cell}≈4.418\times10^5$ A/m | Volume supported as unit-cell volume; current per-moment use requires reconciliation |
| `jcSOT`, `Hext`, `psjSHE` | 0 A/m², zero T, $(1,0,0)$ | Configuration II in [1, Fig. 1(c); 3, Fig. 2(a)] has polarization normal to kagome plane | Current undriven case is a software example; reproduction needs explicit current, field, and pulse schedule |
| SOT gates and `chi` | DL on, FL off, ratio 0 | [1–3] primary model equations contain damping-like SOT | Preserve DL-only default; FL ratio source **REFERENCE MISSING** for an extension |
| STT parameters | Current 0, gates off, polarization $(0,0,1)$, efficiency 0.8, ratio 0 | No STT reproduction target in these three SOT model equations | Not applicable to initial deterministic SOT improvement; provenance **REFERENCE MISSING** if enabled |
| Thermal settings | Off; 100 K; shuffled RNG; noise amplitude in `integrate_llg.m:219` | Current implementation only; no qualifying stochastic-method reference examined for this branch | **REFERENCE MISSING** for a validated stochastic scheme; keep off in initial gates |
| Numerical/output controls | RK4 mode 1; 0.5 fs; 1 ps blocks ×10; stride 200 | [1, p. 134433-3; 3, p. 022407-3] report 5 fs; [1, Fig. 2] uses ns pulses | Existing settings retained as legacy observations; not source-matched numerical protocol |
| Other controls | Open boundary; DMI on; no wall, dipole, pinning, restart, or distribution I/O | Code defaults; source correspondence described above | Disabled modes are not validated capabilities |
| Legacy TM/RE values | `Jfefe=-2.835e-21 J`, `Jfegd=1.09e-21 J`, `gTM=2.2`, `musTM=2.217*mub` | Inherited ferrimagnet comments | Not applicable to current occupied Mn mapping; remove only in approved cleanup |

Physical geometry, reference baselines, digitized figure data, and exact tolerances for paper agreement are **TO BE RESOLVED**. Do not derive wall velocity from `d` until its physical coordinate interpretation is sourced.

## Effective Fields and Units Check

For occupied sites use

$$
\mathbf B_i=-\frac{1}{\mu_i}\frac{\partial E}{\partial\mathbf s_i}.
$$

This follows the effective-field definitions of [1–3, model sections], interpreting their moment-based derivative as flux density. With the proposed positive-sign unique-bond DMI convention, differentiation gives

$$
\mathbf B_i^{ex}=-\frac{1}{\mu_i}\sum_jJ_{ij}\mathbf s_j,\qquad
\mathbf B_i^{DMI}=\frac{1}{\mu_i}\sum_j\mathbf D_{ij}\times\mathbf s_j,\qquad
\mathbf B_i^{an}=\frac{2K_i}{\mu_i}(\mathbf s_i\cdot\mathbf e_i)\mathbf e_i.
$$

The DMI sign differs algebraically from the negative-DMI-energy convention in `sections/hamiltonian.tex`; changing energy sign requires changing the vector convention, not importing the ASTRO field expression unchanged. The paper-to-stencil orientation is unresolved.

$[J/\mu]=J/(J/T)=T$, and $[\gamma B]=s^{-1}$. External input to this flux-density RHS is $\mathbf B_{ext}=\mu_0\mathbf H_{ext}$ when a source supplies A/m. For a quoted Oe field, 100 Oe corresponds to 0.01 T in this conversion. Code variable `Hext` is already summed with tesla fields; passing a raw Oe or A/m value would be wrong.

Torque field amplitude in SI is

$$
B_D=\frac{\hbar_{SI}\theta_{SH}J_c}{2eM_st_{mag}}.
$$

Its dimensions are $(J\,s)(A/m^2)/[C(A/m)m]=T$. Since `hbar` is numerically in eV s, $\hbar_{SI}=e\hbar_{eV}$ cancels the explicit $e$ in this numerical expression. The actual scale discrepancy is the magnetization/thickness mapping (A05), not an automatically missing elementary charge.

## Dynamical Equation and Nonconservative Torques

The current RK4 kernel is algebraically consistent with the explicit LL form below for occupied unit spins. Here $q$ denotes SOT or STT, each with its own polarization; current `scalgpu=gamma/(1+alpha^2)` supplies the prefactor:

$$
\dot{\mathbf s}=\frac{\gamma}{1+\alpha^2}\left[
-\mathbf s\times\mathbf B-\alpha\mathbf s\times(\mathbf s\times\mathbf B)
+\sum_q\left\{-B_D^q\mathbf s\times(\mathbf s\times\boldsymbol\sigma_q)
+\alpha B_D^q\mathbf s\times\boldsymbol\sigma_q
+B_F^q\mathbf s\times\boldsymbol\sigma_q
+\alpha B_F^q\mathbf s\times(\mathbf s\times\boldsymbol\sigma_q)\right\}\right].
$$

DLT and Gilbert signs follow [1–3, model equations]; the full FLT extension and its damping-dependent terms follow the local `eq:ll-with-flt` algebra. No FLT material value is adopted. Do not remove these cross terms when separating the RHS. Verify the residual of the original implicit Gilbert equation at finite damping.

Thermal fields are generated outside conservative fields in the current loop but then included in `hh*`. A deterministic RK4 correction does not select a stochastic interpretation. Stochastic LLGS, noise statistics, and field-like/STT extensions require their own scientific scope and evidence. LLB is not applicable: the modeled occupied spins have fixed length.

## Initial Conditions, Solver, and Outputs

Current default initial state: three parity-dependent directions separated by 120° in the spin-space $yz$ plane, offset by 121°, with zero vacancies. The lattice and local-axis correspondence must be resolved before assigning a paper state label. `dwcalc` currently produces a legacy collinear pattern rather than the published wall. `loadstartm` is inert.

Current solver: normalized frozen-field RK4, fixed timestep 0.5 fs, no adaptive tolerance, no validated coupled order. Stochastic noise is off; `rng('shuffle')` still runs and the console misleadingly reports a random distribution although the active Mn mask is deterministic.

Proposed output contract: save an explicit configuration, mask/site map, initial state, sampled occupied-spin or masked component arrays, update-derived times, units, code revision, solver and hardware information, seed/noise policy, and the requested observables into a chosen candidate path. Include the initial state at $t=0$ and the terminal state at the declared runtime; do not duplicate a block seam. Exact schema and compatibility with downstream workspace readers require approval. No such downstream reader is tracked here.

Publication observables: switching destinations and octupole angles [1, Fig. 2 and Table I]; sub-configuration final angles/polarities [3, Figs. 2–5]; spatial cluster-angle profile, wall position, velocity, and local oscillation [2, Figs. 1–4]. Do not assume that averaging all occupied spins defines the published octupole: derive and cite the site/cluster transformation and normalization before implementing it.

## Implementation Mapping

First preserve a legacy candidate record and define its limitations. Implementation then follows the gates below; do not begin a later stage until the preceding applicable gate passes. Each stage should be a small reviewable change. Paper reproduction and optional physics are later scopes, not incidental additions to solver work.

| Stage | Implementation | Verification gate |
|---|---|---|
| 0. Evidence and legacy record | Resolve target model, structure/normalization gaps, and safe run procedure; record defaults and a deterministic candidate without overwriting user data | All structure rules used by subsequent work have qualifying source locators; legacy input/output/environment record exists and is explicitly labeled as containing known defects |
| 1. Structure | Make Mn occupation, sublattice map, axes, and physical coordinates explicit; validate loaded masks before derived arrays | No duplicate physical sites; vacancies stay empty; mappings match sourced cell; compatible periodic replication and separately sourced open termination pass |
| 2. Interaction topology | Derive reciprocal exchange and antisymmetric DMI bonds with shell, units, and multiplicity; resolve scaling factors | Every active bond has the appropriate reverse coefficient, no open wraparound, periodic seams match bulk, shell counts and source mapping pass; current 12 nonreciprocal bonds are explained and removed only under the approved contract |
| 3. Static energy/fields | Separate conservative-field evaluation from orchestration; implement the agreed energy/field convention | Each term passes analytical and energy-gradient checks; no unapproved factor/sign/parameter change; energy and units agree with selected source |
| 4. RHS and drives | Separate LLGS RHS from field construction and torque amplitudes; define pulses and independently gated flags; reject unsupported paths | Implicit Gilbert residual passes; all torque-flag combinations and signs are exercised; invalid solver/configuration fails before GPU integration |
| 5. Integrator and clock | Recompute fields on complete coupled stage states; correct step/state allocations, vacancy handling, and saved times | Coupled fourth-order convergence in an appropriate smooth deterministic case; finite unit occupied spins and zero vacancies; identical saved states/times for different block sizes; exact terminal update count |
| 6. Simulations and output | Explicit safe persistence and sourced octupole/cluster observables; separate switching, angular scan, and ribbon protocols | Save/load round trip and provenance pass; approved paper-specific benchmark checks pass with declared tolerances and coordinates |
| 7. Performance | Profile; keep states/fields on GPU, vectorize masks, reduce histories and transfers without changing physics | Compare against stage-6 candidate within justified hardware-aware tolerances; measured speed/memory benefit; no output or stochastic convention drift |
| 8. Optional branches | Separately scoped thermal dynamics, dipoles, restart, pinning, or legacy solver restoration | Branch-specific scientific reference, focused benchmark, vacancy/boundary interaction checks, and documented support status; not required to complete initial deterministic SOT work |

Affected existing files: `distrib.m` and `systemgeneration.m` for structure/initial states; `integrate_llg.m` for topology, fields, stages, time, and GPU residency; `atomgpurk4.m` for the RHS boundary; `main.m` for parameters, dispatch, pulse/output policy; optional kernels and `dipole_calc.m` only under explicit scope. Start with minimal functions for a verified structure/topology, conservative fields, RHS, and run/output orchestration. Do not introduce a general simulation framework or copy ASTRO wholesale.

Architectural boundary: `verified geometry -> verified interactions -> conservative fields + explicit drive amplitudes -> LLGS RHS -> integrator -> observables/save`. The RHS receives fields and torque inputs; it must not construct geometry, interpret material tables, or derive coefficients. Keep deterministic conservative fields, nonconservative drives, noise, and persistence separable.

## Verification and Acceptance Criteria

### Checks executed during this audit

The initial command `matlab -batch "fprintf('MATLAB %s\n',version); fprintf('CWD %s\n',pwd); fprintf('GPU count %d\n',gpuDeviceCount);"` reported MATLAB R2024b Update 1 and one GPU. The audit probe identified an NVIDIA GeForce RTX 3060 Laptop GPU, compute capability 8.6, toolkit 12.2, graphics driver 596.08. A subsequent sandbox launch failed with a WSL socket error and was rerun with approval outside the sandbox; this was a launch issue, not a numerical failure.

A single inline MATLAB batch probe evaluated the prefix of `main.m` ending before `integrate_llg(); toc`, then overrode only numerical probe controls and invoked `integrate_llg`. It did not execute the `save` line or modify source. Exact probe settings: `gpusteps=6`, `gpurun_number=2`, `savetstep=2`, `final_m_savestep=3`, `gpusave=6*tstep`, `runtime=2*gpusave`, and `t=(1:12)*tstep`; all physical inputs retained default values. These checks and results were:

| Procedure | Observed result | What it establishes |
|---|---|---|
| Count `atomtype_==1` and `==3` | 16 occupied / 9 empty | Actual default occupancy |
| GPU integration, then finite/norm/mask checks over all saved frames | Six frames; finite components; max occupied norm error 0; max empty norm 0 | Short default-path operation after vacancy reset, not convergence or published agreement |
| Count loop updates vs allocated block steps | 10 actual / 12 nominal; first saved label $5\times10^{-16}$ s | A03 reproduced |
| `atomgpu(0,0,0,1,alp,tstep,0,0,0)` on CPU | Nonfinite components | A07 reproduced |
| `atomgpurk4(0,0,1,1,alp,tstep,0,0,0)` on CPU | Argument failure, as in PC startup dispatch | A08 reproduced; predictor-corrector trajectory was not run |
| Set initial `mm*tmp`, `dipolemode=1`, clear consumed field names, run `dipole_calc` | `hdipolex_` absent; `cpuhdipolex` produced | A11 output-interface failure reproduced |
| Build effective exchange matrix from each `circshift` source index and `hex_constant.*Asim*/Jintra`, restricting both endpoints to occupied slots; compare with transpose | Max asymmetry 0.5; 24 directed mismatches | A02 reproduced independently of spin evolution |
| `checkcode(file,'-id')` for every root `.m` file | Counts: 0 for all three kernels, constants, and initializer; 2 dipole; 2 distribution; 1 integration; 13 main; total 18 | Code Analyzer executed; zero diagnostics in a kernel does not validate its physics or numerics |

Direct PDF checks used `pdftotext -layout` for the three authorized indexed PDFs and visual inspection of rendered model pages: PRB pp. 134433-2 and -3; APL p. 072402-2; APL p. 022407-3. Temporary extraction/render files were placed in `/tmp`. No new benchmark/test files were added and no simulation-result MAT file was saved.

Not run: full default-length simulation; 5 fs paper runs; RK4 convergence; complete energy-gradient tests; DMI reciprocity; thermal ensemble checks; dipole modes 2/3; wall propagation; 360-orientation scan; complete paper reproduction. Confirmed software defects do not by themselves quantify changes to published outcomes.

### Proposed focused checks

- Structure/topology: hand-enumerate one sourced cell and a small replicated cell. Require exact occupation, sublattice, neighbor, shell, and seam agreement. Compare open and periodic modes explicitly; do not infer true periodic support from `bc=0` alone.
- Fields: compare analytical isolated-bond/site results and centered finite-difference energy gradients over an appropriate dimensionless perturbation sweep (for example $10^{-4}$ to $10^{-7}$). Proposed double-precision relative target $10^{-6}$ at a converged perturbation; use an absolute scale near zero. Final tolerances must follow measured conditioning and cannot mask a factor-of-two or sign error.
- LLGS: substitute the explicit RHS into the implicit Gilbert equation for nonzero damping, fields, DL and FL torques. Proposed normalized residual $10^{-12}$ on a CPU double-precision analytical probe; record actual GPU differences before setting a GPU tolerance. Verify tangent RHS and signs under current/polarization reversal.
- Integration: use a genuinely coupled two-spin or sourced small-cell case with smooth drive and several timestep halvings. Expected asymptotic error ratio is near 16 for fourth-order global accuracy; fit an order near 4 over at least three intervals before rounding/noise dominates. A fixed-field single spin alone cannot catch A01. Also check energy conservation at zero damping/no drive and decay with damping/no drive, using step-refinement rather than arbitrary exact equality.
- Masks/clock: proposed occupied norm tolerance $10^{-12}$ for double-precision short deterministic runs, exact zero vacancies, finite values, exact integer update counts, and times reconstructed from those counts. Include a block-boundary case, multiple block sizes, initial time zero, and terminal time; equality of configuration and integration sequence should make deterministic block partitioning irrelevant.
- Paper-specific comparisons after earlier gates: [1, Fig. 2/Table I] six initial states and current/field sign destinations; [3, Fig. 2] representative 55°/175°/295° and 30°/150°/270° states before a full 360-angle sweep; [2, Figs. 1–4] relaxed wall structure, chirality/current direction, velocity vs current and oscillation onset. Treat 545.3 m/s at $9\times10^{10}$ A/m² and approximately 770 nm wall extent as published targets [2], not achieved code results. Fix measurement definitions and tolerances from the available figure/source resolution before execution.
- Spatial studies: converge ribbon size and edge effects with central and edge-resolved observables. A periodic cell and an open finite ribbon require different checks; ordinary boundaries are not absorbing boundaries.

There is no tracked accepted baseline. New tests/benchmark commands should be added and documented only in their approved implementation stages. Protect user result files and use separately named candidate directories; never run the present whole-workspace `save` during read-only probes.

## References

[1] Zhengde Xu, Xue Zhang, Yixiao Qiao, Gengchiau Liang, Shuyuan Shi, Zhifeng Zhu, Deterministic spin-orbit torque switching including the interplay between spin polarization and kagome plane in Mn3Sn, Physical Review B **109**, 134433 (2024). [DOI: 10.1103/PhysRevB.109.134433](https://doi.org/10.1103/PhysRevB.109.134433). Consulted reviewed summary and indexed main PDF; directly inspected model/protocol pages 134433-2 and -3. Paper ID: `AFM/2024-Deterministic spin-orbit torque switching including the interplay between spin-PRB-zhengde Xu`. Summary path relative to the papers root: `AFM/2024-Deterministic spin-orbit torque switching including the interplay between spin-PRB-zhengde Xu/2024-Deterministic spin-orbit torque switching including the interplay between spin-PRB-zhengde Xu-summary.md`. Supplemental material was not inspected.

[2] Zhengde Xu, Yue Zhou, Xue Zhang, Yixiao Qiao, Zhuo Xu, Dingfu Shao, Zhifeng Zhu, Spin-orbit torque-driven chiral domain wall motion in Mn3Sn, Applied Physics Letters **127**, 072402 (2025). [DOI: 10.1063/5.0275467](https://doi.org/10.1063/5.0275467). Consulted reviewed summary and indexed main PDF; directly inspected p. 072402-2 (model and Fig. 1), with text extraction for result locations. Paper ID: `AFM/2025-Spin-orbit torque-driven chiral domain wall motion in Mn3Sn-APL-Zhengde Xu`. Summary path: `AFM/2025-Spin-orbit torque-driven chiral domain wall motion in Mn3Sn-APL-Zhengde Xu/2025-Spin-orbit torque-driven chiral domain wall motion in Mn3Sn-APL-Zhengde Xu-summary.md`. The reported field symbol and Zeeman normalization use mixed moment/magnetization notation; conversion to per-site SI energy still requires reconciliation.

[3] Boyu Zhao, Zhengde Xu, Xue Zhang, Zhenhang Kong, Shouzhong Peng, Shuyuan Shi, Zhifeng Zhu, Investigation of sub-configurations reveals stable spin–orbit torque switching polarity in polycrystalline Mn3Sn, Applied Physics Letters **127**, 022407 (2025). [DOI: 10.1063/5.0273330](https://doi.org/10.1063/5.0273330). Consulted reviewed summary and indexed main PDF; directly inspected p. 022407-3 (Hamiltonian, LLGS, parameters, Fig. 2). Paper ID: `AFM/2025-Investigation of sub-configurations reveals stable spin–orbit-APL-boyu zhao`. Summary path: `AFM/2025-Investigation of sub-configurations reveals stable spin–orbit-APL-boyu zhao/2025-Investigation of sub-configurations reveals stable spin–orbit-APL-boyu zhao-summary.md`.

Local paper root: `/mnt/c/Users/zzf-m/OneDrive/papers`. The reviewed notes contain PDF cross-checks and extraction limitations; exact equation/parameter comparisons in this audit additionally used the source model pages. The switching topic `/mnt/c/Users/zzf-m/OneDrive/papers/_knowledge/topics/mn3sn-spin-orbit-torque-switching.md` supplied navigation and model-scope comparison, not physical construction data. Repository source and equation notes describe implementation/conventions and do not replace the qualifying publications.

## Assumptions, Missing References, Risks, and Approval Questions

The audit assumes the user wants this checkout documented as it exists and a staged improvement proposal. Association with all three papers is user-supplied provenance; this revision has no tracked evidence proving it is the exact code used to generate their figures.

**REFERENCE MISSING:** exact lattice vectors/basis/metric; parity-to-crystallographic and sublattice map; shell counts/distances/multiplicities; source-to-code exchange scaling; directed DMI vectors and edge factors; periodic compatibility and physical termination; exact g-factor provenance; optional STT/FLT coefficients; validated stochastic method. Search scope: active source, the three reviewed summaries and authorized main PDFs, relevant `_knowledge` index/topic, targeted presentation catalog/index, and local equation notes. Main papers provide schematic/parameter evidence but not a complete numerical mapping for this builder. No supplemental files, external datasets, or online search were used. These gaps block the corresponding physical construction or extension stage, while pure software investigations can continue under their own scope.

**TO BE RESOLVED:** anisotropy factor, paper pair sums, per-site/volume Zeeman notation in the domain-wall paper, magnetic thickness versus site spacing, octupole/cluster definition, downstream result compatibility, accepted baseline, and quantitative reproduction tolerances. Do not treat schematic topology or approximate coefficient similarity (17.53/2≈8.765 meV versus 9 meV) as evidence of an approved normalization.

Important risks: correcting nonreciprocal fields or normalization changes trajectories; true coupled RK4 changes the historical algorithm; replacing time labels changes reported times/velocities; the small finite patch can contain substantial edge effects; stochastic and dipole changes introduce independent scientific decisions. Separate preservation of legacy observations from acceptance of scientifically corrected behavior.

Before source implementation, review this plan and select the target scope:

1. Start with deterministic infrastructure/topology/LLGS correctness and safe output, or prioritize a particular paper reproduction after those gates?
2. Supply or authorize access to the exact structure, bond normalization, boundary, and simulation inputs needed to resolve the source-gated rows; the main PDFs inspected here are insufficient for all of them.
3. Approve a concrete model contract and the affected stages before changing physics, solver behavior, output compatibility, or defaults. The current task authorizes the audit and these three documentation files only.

---

## Visualization task — 2026-10-21

**Status: implemented and verified on 2026-10-21.** The user first approved appending this section while retaining the audit, then explicitly approved its visualization scope. The earlier solver/model repair stages remain separate and unapproved. This task creates documentation views of existing behavior and a visibly schematic cell illustration; it does not resolve the audit's physical-model gaps.

### Goal and Success Criteria

Create two linked, offline interactive pages: `docs/unit-cell.html` and `docs/structure-exchange.html`. The first explains the illustrated repeating motif and finite array; the second connects a selected site to its exchange/DMI reads, anisotropy axis, array indices, coefficients, and source excerpts. Use the existing notebook as historical visual context. Success requires agreement with current MATLAB arrays and field expressions, functional controls, readable desktop/mobile views, and visible evidence limits.

Non-goals: change MATLAB simulation behavior, repair exchange reciprocity or periodic boundaries, build a crystallographically validated Mn3Sn structure, simulate dynamics, reproduce a publication, install dependencies, or overwrite result files. A schematic drawing must never become simulation geometry.

### Model Classification and Equation Family

Class 4 describes the repository's intended material-specific spatial atomistic scope; its crystallographic mapping is unresolved. The viewer is a software/documentation companion to fixed-length LLGS. Each occupied Mn entry has an independent unit spin; three code-defined easy-axis classes do not establish publication sublattice IDs. Empty entries remain zero and do not carry spin arrows. The reviewed knowledge-base topic [V1, Scope and evidence base; Comparison] distinguishes three-sublattice switching models from experimental systems; it does not supply a crystal builder for this spatial checkout. Electronic structure and Hall-response prediction are outside scope.

### Existing Behavior and Reference Implementations

Inspected `AGENTS.md`, `README.md`, the complete existing audit plan, `main.m`, `distrib.m`, `systemgeneration.m`, and the setup/shift/field portion of `integrate_llg.m`. The untracked notebook, documentation, and audit existed before this task and must be preserved. `main` is the current branch; no tracked source changes were present at preflight.

The default 5×5 mask is documented as 16 occupied sites. The earlier 5×6 audit reports 20 sites, 30 undirected exchange connections, and 12 unequal reverse-weight pairs; these are prior results to recheck, not tests executed in this task. Notebook transcription items `P19-ATOMISTIC-5X6`, `P19-NEIGHBOR-DIRECTION-SKETCHES`, `P20-OCCUPATION`, and `P21-ANISOTROPY-EXCHANGE-MASKS` provide historical comparison [V2].

Read-only reference inspection: ASTRO `AGENTS.md` and `README.md` establish that its species, axes, supported solvers, and baselines differ. Fe2O3 `docs/unit-cell.html` and `docs/structure-exchange.html` supply interaction-design examples only. No geometry, material data, or masks will be copied from either model.

Equation sources inspected under `/mnt/c/Users/zzf-m/OneDrive/presentations/equations`: `AGENTS.md`, `README.md`, `spintronic_methods.tex`, `sections/hamiltonian.tex` (`eq:atomistic-hamiltonian`, `eq:atomistic-effective-field`, `eq:atomistic-pair-fields`), `sections/llg_equation.tex` (`eq:llg-gilbert`, `eq:ll-with-flt`), `sections/atomistic_boundaries.tex`, and the Mn3Sn parameter table in `sections/material_parameters.tex`. These govern explanatory units, pair counting, and separation of boundaries and torque conventions. They are not independent physical structure evidence. The negative-DMI-energy convention in the equation notes must not be silently assigned to the current Mn3Sn stencil.

### State Variables, Geometry, and Boundaries

Array identity is one-based `(row W, column L)`, with `1=occupied Mn` and `3=empty slot`. The unit-cell page will explicitly call its repeat outline an **illustrated motif**, not a verified crystallographic cell. Provide a notebook-derived schematic and an exact array view, both with selectable indices. Illustration coordinates have arbitrary display units; spin-space axes are shown separately from array axes. Physical lattice vectors, lengths, Sn positions, stacking, and crystallographic directions remain unavailable. Do not fabricate them.

| Physical detail | Adopted value or rule | Builder or interaction mapping | Qualifying reference and exact locator | Evidence status |
|---|---|---|---|---|
| Material and magnetic-order context | Mn3Sn, noncollinear three-sublattice order as literature context | Introductory explanation only | [V1, Scope and evidence base; Comparison], underlying publications [1, Sec. II] and [3, model section] in the audit | Context only; not proof of array-to-crystal mapping |
| Metric lattice vectors, fractional basis, Sn coordinates, stacking | No numerical physical construction adopted | Excluded from viewer data | No qualifying structure data established in this task | **REFERENCE MISSING**; blocks a metric or complete chemical-cell reconstruction |
| Physical site/sublattice mapping and nearest-Sn axes | Retain code classes with neutral labels; no publication site IDs | `distrib.m` parity classes and axis components | [V1] does not identify the code parity map | **REFERENCE MISSING** for physical mapping; code values shown only as implementation facts |
| Repeat motif and site drawing | Historical sketch, arbitrary spacing; array index display is exact | [V2] pages 19–20 plus occupancy from `distrib.m` | Notebook is historical evidence, not qualifying crystallography | Schematic only; keep separate from model input |
| Neighbor shells, physical distance, coordination | Show actual nonzero field reads, not claimed physical shells | Eight exchange shifts and three DMI masks in `integrate_llg.m` | No qualifying code-to-shell map established | **REFERENCE MISSING** for physical assignment; exact read counts are software diagnostics |
| Exchange coefficient and pair multiplicity | Display `Jintra` and receiver-specific weights separately | `hex_constant .* Asim*` | Audit references [1–3] concern distinct published models; no conversion to every code bond established | **REFERENCE MISSING** for physical normalization; preserve unequal reverse weights |
| Directed DMI vectors | Display signed component reads only | `hdmi_constant*`, `hdmi_y`, `hdmi_z` | No qualifying stencil-to-vector map established | **REFERENCE MISSING**; no physical DMI-vector or chirality claims |
| Boundary termination and periodic registry | Show current open mode and diagnostic `bc=0` behavior | Existing mask cuts plus zeroed or wrapped shifts | No qualifying physical termination/periodic-cell mapping established | **REFERENCE MISSING** for physical geometry; label `bc=0` as legacy periodic shifts with retained cuts |
| Occupancy, storage order, display IDs | Preserve exact indices and distinguish rendered sites from repeated motifs | Current scripts | Software bookkeeping; scientific citation not required | Verify against MATLAB output |

Repeat controls enlarge the illustrative pattern or regenerate the selected finite array. They must not claim that tiling a terminated mask constructs a periodic crystal. Use 5×5 and notebook 5×6 presets, plus a larger patch containing every row-modulo-four/column-parity class and interior sites. Record omitted edge reads, vacancy reads, and wrapped reads distinctly. There is no production-boundary transition in this documentation task.

### Hamiltonian

No new Hamiltonian is implemented. The audit above records a proposed unique-bond energy and unresolved publication normalization. Because receiver-dependent exchange weights are not symmetric, do not present that proposed energy as a generator of all current fields. Explain the actual directed sums and compare both directions of a connection. Anisotropy's existing field corresponds algebraically to $E_i=-K_{sim}(\mathbf e_i\cdot\mathbf s_i)^2/2$ under $\mathbf B_i=-\partial E_i/\partial\mathbf s_i\,/\mu_i$; label this a derivation of current behavior, not a resolved published coefficient. DMI energy mapping remains unresolved. Torques and thermal noise are not bond energies.

### Parameters and Evidence

The following are source-code snapshot values, not newly adopted material parameters. Refresh from source and record source hashes when generating viewer data.

| Display quantity | Value/units | Provenance and use |
|---|---|---|
| Exchange | `Jintra=2.8083e-21 J` | `main.m`; multiply by `Asim*/Jintra` and receiver `hex_constant` in the table |
| Receiver scale | 0.5 or 1 | `distrib.m`; expose both directed values without symmetrizing |
| DMI | `Dsim=0.833 meV`, converted with `ele=1.602e-19 C` | `main.m`, `constantfile.m`; show signed terms and receiver scale |
| Anisotropy | `Ksim=0.196 meV`, converted with the same constant | `main.m`; show local axis and current field prefactor |
| Moment | `musRE=3*mub`, `mub=9.274e-24 J/T` | Current occupied-site field denominator; no rare-earth interpretation |
| Easy axes | $(0,\sin\theta,\cos\theta)$; $\theta=30°,150°,270°$ | `distrib.m`; code classes, not verified crystallographic labels |
| Initial spins | Angles 121°, 241°, 361° at azimuth 90° | `systemgeneration.m`; optional arrows explicitly labeled initial state |
| Geometry scale | Arbitrary display units | No adoption of `d=0.4 nm` as a verified Mn–Mn distance |
| Boundaries and size | Default `bc=1`, 5×5; additional diagnostic presets | `main.m`; viewer choices do not change production inputs |

Damping, torque efficiencies, time steps, thermal parameters, and dipole parameters do not enter these static views. Their existing values remain documented in the audit; no viewer dynamics or new physical controls are proposed.

### Effective Fields and Units Check

For an occupied receiver $i$, expose the exact current exchange operation as

$$
\mathbf B_i^{ex}=-\frac{h_i}{\mu_i}\sum_q A_{i,q}\mathbf s_{j(i,q)},
\qquad h_i=\texttt{hex\_constant}(i).
$$

Here $q$ labels an implemented shift, not a verified crystallographic shell. The table reports $h_i A_{i,q}$ in J, and the field coefficient $-h_i A_{i,q}/\mu_i$ in T. For DMI display the signed scalar reads in $B_y$ and $B_z$, including masks and receiver scaling; $B_x=0$. For anisotropy display $\mathbf B_i^{an}=K_{sim}(\mathbf e_i\cdot\mathbf s_i)\mathbf e_i/\mu_i$. These are source-expression descriptions. $J/(J/T)=T$; no extra $\mu_0$ enters these formulas. Empty-slot denominators are bookkeeping, not magnetic moments.

### Dynamical Equation and Nonconservative Torques

Not applicable to implementation: no time integration or RHS is added. Existing fixed-length LLGS and its damping-dependent SOT/STT cross terms remain as documented above and in `eq:ll-with-flt`. Optional arrows show only the current initializer. Neither arrow animation nor static selection will be represented as computed dynamics. Drives and noise are outside the structure/interaction explorer.

### Initial Conditions, Solver, and Outputs

Use `dwcalc=0` initializer data for the optional spin overlay; show easy axes with a different glyph. Presets are finite static arrays, not equilibrium or simulated states. Solver, seed, timestep convergence, and trajectory outputs are not applicable. The new outputs are local HTML and small local JS/CSS/data companions. They must open through `file://` without a server, CDN, framework, or installation. Existing notebook assets are reused by links rather than rewritten. Export/check artifacts go in a temporary directory; `finalnew.mat` and `atomtypee.mat` must remain untouched.

### Implementation Mapping

Use the smallest shared display-data module needed to keep both pages consistent. A reproducible MATLAB export/check harness reads existing construction and field expressions without running the integration loop or save commands. It records source/configuration provenance and exports only static display data. Its exact extraction boundaries must be reviewed before execution. No source refactor is needed.

| Stage | Implementation | Verification gate |
|---|---|---|
| V0. Scope | Approve this appended plan and the schematic evidence boundary | Explicit approval; unresolved physical construction excluded |
| V1. Structure/data | Add static export/check harness and `docs/unit-cell.html` with local companions; selectable sites, repeated schematic motif, exact array view, labels, zoom/pan/reset, optional initial-spin/axis overlays | Match every occupied/empty entry and axis class to existing MATLAB output on 5×5, 5×6, and larger patch; distinguish motif marker counts from simulated sites; inspect notebook image before redrawing |
| V2. Interaction explorer | Add `docs/structure-exchange.html`; synchronize selected site/read across schematic, array, table, and source excerpt; exchange, DMI, and onsite anisotropy selections | Compare every displayed coefficient and endpoint to independent MATLAB impulse probes for both boundary flags, including edges/corners/vacancies/seams; preserve and expose reverse-weight discrepancies |
| V3. Static-field explanation | Add signs, units, directional weights, and current field expressions | Basis-spin responses agree with source expressions; source excerpts/hash records match; no unsupported shell or physical-vector claim |
| V4. RHS and integrator | Not applicable: visualization-only task | Confirm production MATLAB files unchanged |
| V5. Browser and documentation | Cross-link pages and notebook; update README with scope and refresh/check procedure; maintain this task's actual verification record | Exercise all controls, inspect desktop/narrow renders, check console/assets/offline loading; review changed and untracked files |

Do not advance until the preceding applicable gate passes. Boundary: existing source -> static provenance-tagged data -> display state -> render. No visualization module enters the LLGS computation. If a later scientific implementation is requested, the earlier audit's verified-fields/torques -> RHS -> integrator boundary still applies.

### Verification and Acceptance Criteria

1. Run a CPU-only static MATLAB harness in an isolated temporary output directory. Disable fixed-distribution load/save and wall initialization. Execute only reviewed setup/field snippets; never run `main.m` or the integration/save loop. Compare masks, axes, initial spins, neighbor endpoints, and raw coefficients on all presets and both boundary flags.
2. Independently excite each occupied source spin component to reconstruct exchange and DMI field operators. Require exact integer indices/signs/masks; compare numerical coefficients to normalized double-precision tolerance $10^{-12}$, with an absolute zero check after normalization. A failure must be investigated rather than hidden by increasing tolerance. Anisotropy is checked with component basis inputs. These tests validate the display against code, not physical correctness.
3. Recheck the prior 5×6 occupation/connectivity and exchange-asymmetry observations. Include all inequivalent row/parity classes and bulk, edge, corner, and wrapped-read cases; count only occupied endpoints as magnetic neighbors. Do not require reciprocity to pass when faithfully displaying a known nonreciprocal operator.
4. Check embedded source excerpts and hashes against the source files and document the refresh procedure. Run JS syntax checks with available Node. Use an available browser/automation environment without installing packages; exercise selection, linked highlighting, boundary/term selection, size/repetition, padding, arrows, pan/zoom, and reset. Inspect rendered desktop and narrow views and local asset loads. Report unavailable browser checks explicitly.
5. Run `git diff --stat` and `git diff --check`, include untracked additions in review, and confirm no production-source or result changes. No dynamics, GPU smoke, solver convergence, physical energy-gradient certification, or publication benchmark is needed for these static views.

#### Executed verification and deliverables

Created `docs/unit-cell.html` and `docs/structure-exchange.html`, with shared `mn3sn-view.css`, `mn3sn-view.js`, `mn3sn-model.js`, and the generated `mn3sn-data.js`. Added `docs/export_visualization_data.m` and `docs/check_visualization_data.cjs`. Updated README with links, scope, refresh commands, controls, and validation limits; added the four generated/cache ignore patterns required by global instructions. Existing notebook files, audit content, production MATLAB source, and simulation results were preserved. No branch, commit, dependency installation, or network access was needed.

| Check / exact command or procedure | Result |
|---|---|
| `matlab -batch "addpath('docs'); export_visualization_data('docs/mn3sn-data.js')"` | Passed on MATLAB 24.2.0.2740171 (R2024b) Update 1, CPU only. Ten sizes: 5×5 and each combination of W=5/9/13 with L=6/10/14; both boundary flags; every occupied x/y/z impulse checked against unmodified static exchange, DMI and anisotropy expressions. Normalized differences below $10^{-12}$; expected zero responses exactly zero. No GPU or integration invoked. |
| `node docs/check_visualization_data.cjs` | Passed on Node v20.20.1: all source hashes/excerpts, occupation/axes/spin invariants, 367,700 operator entries, reverse weights, and connection counts. Coverage included 170 active wrapped reads, 246 empty-source reads, 8,064 empty-receiver reads, 9,970 coefficient cuts, and 302 open-edge reads across configurations and terms. |
| `node --check docs/mn3sn-model.js` and `node --check docs/mn3sn-view.js` | Passed syntax checks. |
| 5×6 static reference | 20 occupied sites, 30 exchange connections, 12 unequal reverse pairs; notebook page-20 mask matched exactly. |
| `/home/zhuzibn/.venvs/webpage-capture/bin/python /tmp/check_mn3sn_visualizations.py` | Passed in existing Playwright environment with Chromium 151.0.7922.34. Both pages loaded via `file://`; checked site/read selection, cross-page selection retention, presets, repetition, boundary/term switching, empty-slot visibility, arrows/axes/labels, zoom, drag, keyboard/reset, and 1440 px / 390 px layouts. No page/console/asset errors, HTTP assets, or document-wide horizontal overflow. Temporary script and screenshots are audit artifacts under `/tmp`, not required viewer dependencies. |
| Visual inspection | Inspected original notebook page images 19–21 before drawing, then browser screenshots of the desktop cell view, both mobile pages, and a wrapped-DMI case. Narrow layouts use an adjusted SVG view box. |
| Independent read-only verifier | Independently reconstructed default/5×6/larger-patch occupation and exchange counts; reviewed extraction boundaries and evidence labels; ran the Node comparison. Found and resolved the DMI-refresh guard and fractional URL-index issues. Final focused review reported no remaining material issue. |
| Git review | Tracked-source diff remained empty; checked new/untracked files separately because pre-existing documentation is untracked. Whitespace checks passed. |

Transient failures: the first exporter run rejected the indented `Dsim` assignment before writing output; whitespace-tolerant extraction fixed it. Chromium initially hit sandbox process restrictions and succeeded after approved escalation. An initial browser cross-page assertion ran before page initialization; the check now waits for the destination selection. Final relevant checks passed after the fixes. The exporter explicitly rejects a changed `DMIenable` setting or nonpositive normalization coefficient so that unsupported refresh assumptions cannot silently produce a new snapshot.

Not run: dynamics, GPU simulations, coupled-solver convergence, physical Hamiltonian gradient validation, or publication reproduction. The snapshots validate software correspondence only; all physical evidence gaps in the structure table remain unresolved.

### References

[V1] Reviewed topic record: `/mnt/c/Users/zzf-m/OneDrive/papers/_knowledge/topics/mn3sn-spin-orbit-torque-switching.md`, **Mn3Sn spin-orbit-torque switching: geometry, polycrystallinity, and speed**, read completely in this task. Relevant locators: Scope and evidence base; Comparison; From ideal configuration II to polycrystalline operation; Evidence strength and limitations. Underlying model references are the full journal entries [1] and [3] in the audit's References section. Only the topic record and index were newly consulted; linked paper summaries/PDFs outside `_knowledge` were not reopened. The topic does not establish spatial lattice construction.

[V2] Repository notebook transcription: `docs/Mn3Sn-notebook-transcription.md`, items `P19-ATOMISTIC-5X6`, `P19-NEIGHBOR-DIRECTION-SKETCHES`, `P20-OCCUPATION`, and `P21-ANISOTROPY-EXCHANGE-MASKS`; original local images `docs/Mn3Sn-pages/page-19.jpg` through `page-21.jpg`. Transcription inspected in planning; all three original images inspected before drawing. Historical illustration evidence only, not a qualifying crystallographic reference.

[V3] Current implementation: `main.m`, `constantfile.m`, `distrib.m`, `systemgeneration.m`, and `integrate_llg.m`. Code establishes snapshot values, indices, masks, and algebra, not material correctness. Exact excerpts and hashes will accompany exported display data.

Scientific journal entries [1–3] above are retained as the prior audit's bibliography and evidence record; prior PDF inspection is not represented as new inspection. Targeted searches of `papers/_knowledge/papers.jsonl` and `presentation-knowledge/{presentations.jsonl,presentation-catalog.md}` found relevant paper records and no matching presentation record. No network access was used.

### Assumptions, Missing References, Risks, and Approval Questions

The proposed deliverable interprets “unit cell and interactions” as a notebook-based schematic plus an exact current-code explorer. A metric crystallographic unit cell is a distinct scope requiring qualifying lattice/basis evidence. The physical gaps listed in the structure table remain blocking for those details; user approval does not establish missing structure evidence.

Search scope: active code/audit/notebook transcription; processed paper index and reviewed switching topic; targeted presentation index/catalog; equation notes; Fe2O3 display examples. Linked paper summaries/source PDFs outside the pre-authorized `_knowledge` tree and online sources were not accessed. Request specific source access only if the user instead selects a physical crystallographic reconstruction.

Main risks are a convincing schematic being mistaken for a validated crystal, a receiving-site weight being mislabeled a reciprocal bond energy, code DMI signs being mistaken for sourced physical vectors, and copied data becoming stale. Explicit labels, directional tables, source hashes, and independent static probes address these risks without claiming the underlying model is repaired.

**Approved decision:** the user approved V0–V5 on 2026-10-21. Both pages, the narrow static export/check harness, and relevant README/plan updates are complete. The schematic/code distinction, unresolved mappings, and current interaction defects remain explicit; no approval of the earlier solver/physics repair stages is implied.

## Unequal exchange pair visualization — 2026-10-27

User request: visualize every unequal exchange pair identified in the notebook's 5×6 open patch. This is a documentation-view extension for the existing Class 4, fixed-length LLGS model; no simulation source, physical coefficients, or geometry construction changes are in scope.

Reuse the existing MATLAB snapshot and schematic coordinates. Add an overview to `docs/structure-exchange.html` that numbers all unequal unordered exchange pairs, identifies both endpoints, and reports each directed weight separately in units of `Jintra`. Select a bond or table row to inspect its receiver/source in the existing detailed view. The overview follows the existing array and boundary controls, and is hidden for non-exchange terms. Preserve the distinction between an implemented field coefficient and a reciprocal bond energy.

Verification gates: compare each displayed pair and both weights against the independent MATLAB impulse operator already stored in the snapshot, including absent reverse contributions; require 12 unequal pairs out of 30 connections for the 5×6 open case. Verify browser selection, term/boundary/preset changes, keyboard interaction, and desktop/mobile layout. Inspect a rendered overview. Update README with the interaction and evidence scope. No new dependencies, network access, MATLAB runs, or scientific-reference claims are needed.

Completed: added the numbered overview and directional table to the existing interaction page, a pair aggregation helper, and independent impulse comparisons to the existing Node check. Pair selection updates the receiver, selected read, array, and source expressions. README documents the overview. The 5×6 open patch contains 12 unequal pairs with weights 1 and 0.5 in units of `Jintra`, plus 18 equal pairs. The schematic retains the missing `(3,6)` slot.

Validation on 2026-10-27:

- `node docs/check_visualization_data.cjs`: passed source hashes/excerpts and 367,700 existing operator comparisons across ten sizes and both boundary flags; additionally checked every overview endpoint and directed weight against the stored MATLAB impulse operator, pair uniqueness, and unequal counts.
- `node --check docs/mn3sn-model.js` and `node --check docs/mn3sn-view.js`: passed.
- `/home/zhuzibn/.venvs/webpage-capture/bin/python /tmp/check_mn3sn_unequal_pairs.py`: passed with Chromium 151.0.7922.34. Verified all 12 table selections, map/keyboard selection, focus retention, directional labels/weights against the independent impulse data, four presets with both boundaries, term switching, empty receivers, the cell page, and 1440/390 px layouts. No browser errors, network requests, or document-wide horizontal overflow. Chromium required sandbox escalation for the full browser check; the initial sandboxed launch failed before opening the page.
- Inspected desktop and mobile screenshots. Preview files `/tmp/mn3sn-unequal-exchange-overview.png`, `/tmp/mn3sn-unequal-exchange-mobile.png`, and the standalone map `/tmp/mn3sn-unequal-exchange-pairs.svg` are temporary inspection artifacts; the maintained deliverable is the HTML viewer and its local companions.
- No new MATLAB/GPU simulation, physical-model validation, source coefficient change, dependency installation, or network access. The source snapshot and its MATLAB provenance are unchanged. The prior equation-note convention in `sections/hamiltonian.tex` informed the distinction between directed field coefficients and reciprocal pair energies; this task adds no crystallographic or publication-agreement claim.

## Three-spin comparison — 2026-11-02

### Goal and Success Criteria

The user authorized the proposed comparison of this spatial code reduced to three occupied sites against `Mn3Sn_3atoms`, followed by matched Figure 5 trajectories and timestep refinement. Preserve production source and historical results. Success means measured term-level correspondence, explicit accounting for parameter and integration differences, and a reproducible report distinguishing the current algorithm from a diagnostic coupled RK4 variant. This extends the audit; it does not approve production repairs or larger-lattice changes.

### Model Classification and Equation Family

Class 5 for the diagnostic: three noncollinear, fixed-length Mn sublattice vectors with LLGS [1, Sec. II]. The spatial code supplies a three-site software realization; its larger-lattice Class 4 crystallographic limitations remain unresolved. Electronic structure and transport are outside scope.

### Existing Behavior and Reference Implementations

Inspect the current `distrib.m`, `systemgeneration.m`, `integrate_llg.m`, `atomgpurk4.m`, `main.m`, and constants. Compare with `Mn3Sn_3atoms/main.m`, `field_eta.m`, `LLG_solver.m`, corrected `rk4_4llg.m`, and its existing verification tools/report. Source inspection identifies half the reference anisotropy field, approximately 450 times its SOT amplitude at equal current, slightly different exchange and gamma, and fields frozen across RK stages. These differences must be measured separately. The ASTRO README/instructions were consulted only for model boundaries; no ferrimagnetic geometry or parameters are imported.

### State Variables, Geometry, and Boundaries

Use `natomW=natomL=2`, `bc=1`, with distribution I/O, domain walls, dipoles, pinning, and noise disabled. Source inspection gives three occupied entries and one empty entry. Map reference spins `(a,b,c)` to array entries `((2,2),(1,1),(2,1))` by their easy-axis angles `(30,270,150)` degrees. Verify this mapping and all directed interaction coefficients by executing the original setup and Cartesian impulses. The empty `(1,2)` entry must stay exactly zero.

| Detail | Evidence and implementation boundary |
|---|---|
| Three Mn sublattice moments; yz spin plane and +x polarization | [1, Sec. II and Fig. 1], checked against the reference implementation |
| 2×2 storage and index permutation | Software bookkeeping from current source, not a physical lattice-constant claim |
| Three reciprocal exchange pairs and cyclic DMI | Expected algebraic reduction of the source; must pass measured operator checks |
| Physical coordinates, stacking, finite-sample surface interpretation | Not required for the three-vector comparison; existing spatial mapping remains REFERENCE MISSING and is not implemented |

### Hamiltonian

For the reference triangle, use the implemented single-counted pair convention:

$$
E=J\sum_{(i,j)=(a,b),(b,c),(c,a)}\mathbf m_i\cdot\mathbf m_j
+D\hat x\cdot\sum_{(i,j)=(a,b),(b,c),(c,a)}\mathbf m_i\times\mathbf m_j
-K\sum_i(\mathbf e_i\cdot\mathbf m_i)^2-\mu\sum_i\mathbf B_{\rm ext}\cdot\mathbf m_i.
$$

The terms are those of [1, Sec. II]; the exact pair-counting convention is established by the reference code, not independently resolved by the paper's general pair sums. The native spatial anisotropy corresponds instead to coefficient `Ksim/2` in this energy. SOT is nonconservative.

### Parameters and Evidence

Reference values: `mu=3*mub`, exchange field `100.812 T`, `K=0.196 meV`, `D=0.833 meV`, `alpha=0.003`, `thetaSH=0.06`, `Ms=441850 A/m`, thickness `30 nm`, and `gamma=1.760859644e11 rad/(s T)`. Material scales and the 5 fs timestep are from [1, Sec. II]; preserve exact reference-code constants for equality checks and record their rounding relative to the paper's `J=17.53 meV` and `Ms=6*mu/Vcell`. Read native values independently from source. No thermal noise, STT, FLT, dipoles, or interlayer coupling is introduced.

### Effective Fields and Units Check

Use `B_i=-(1/mu)*dE/dm_i` in tesla. For cyclic `(i,j,k)`, exchange is `-(J/mu)*(m_j+m_k)`, DMI is `(D/mu)*cross(xhat,m_j-m_k)`, and anisotropy is `(2*K/mu)*dot(e_i,m_i)*e_i`. Check each term against source impulses and finite-difference energy derivatives. Compare native and explicitly matched coefficients separately. Oe-to-tesla conversion is `1e-4`; the stored `hbar` is numerically hbar/e in SI torque formulas, so do not divide by charge again.

### Dynamical Equation and Nonconservative Torques

With `p=(1,0,0)` and signed `b=(hbar/e)*thetaSH*Jc/(2*Ms*thickness)`, compare

$$
\dot{\mathbf m}=\frac{\gamma}{1+\alpha^2}\left[-\mathbf m\times\mathbf B-\alpha\mathbf m\times(\mathbf m\times\mathbf B)-b\mathbf m\times(\mathbf m\times\mathbf p)+\alpha b\mathbf m\times\mathbf p\right].
$$

This is the explicit form of [1, Sec. II], retaining Gilbert cross terms as in `presentations/equations/sections/llg_equation.tex`, `eq:ll-with-flt` and `eq:llg-gilbert`. Verify the RHS against both original implementations and the implicit Gilbert residual.

### Initial Conditions, Solver, and Outputs

Use reference first-spin seeds `30,90,150,210,270,330` degrees, with successive spins offset by 120 degrees. Figure 5(a) uses seed 210 degrees, current `+5e10 A/m²`, and fields `175,200,225,250,275,300,325,350,375,400,450,500 Oe`. Figure 5(b,c) uses all six seeds, `1000 Oe`, and both current signs [1, Fig. 5]. Apply field/current from 1 to 4 ns in a 5 ns run. Use common physical sample times, explicit update counts, and signed `atan2d(mean_y,mean_z)`; the perfect initial 120-degree state's zero mean has no defined angle.

Compare native parameters/frozen fields, matched parameters/frozen fields, and matched parameters/coupled RK4 with the corrected reference. The coupled variant is diagnostic only. Start at 5, 2.5, and 1.25 fs; include native 0.5 fs or further refinement if needed to characterize failure. Record preprojection norm error as well as normalized norms. Store outputs in a fresh directory under `.sisyphus/`; never run either production `main.m` or overwrite result files. A compiled replay may accelerate sweeps only after short MATLAB comparisons pass.

### Implementation Mapping

| Stage | Implementation | Verification gate |
|---|---|---|
| Mapping and topology | Isolated setup excerpts and impulse measurements | Exactly three occupied sites, correct axes/permutation, three reciprocal pairs, antisymmetric directed DMI |
| Static fields | Native and matched coefficients; explicit triangle energy | Matched term fields agree within 1e-11 T; scaled finite-difference gradient residual below 1e-7 |
| RHS and one-step updates | Original kernels, explicit RHS, optional replay | Relative RHS/Gilbert residual below 1e-12; replay spin discrepancy below 1e-10 over short deterministic checks |
| Integrators | Compare frozen and recomputed stage fields | Record coupled timestep convergence and preprojection errors; do not require the legacy method to pass fourth order |
| Figure cases and report | All listed fields/seeds/signs; plots and source hashes | Report endpoints, transient errors and convergence without tuning parameters to force agreement |

Do not begin a dependent stage until the preceding implementation gate passes. A verified defect in native behavior is a diagnostic result, not permission to repair production. Fields, torque/RHS, integration, and sweep orchestration stay separate. New crystallographic illustrations are not applicable: no new lattice is constructed and the existing viewers remain unchanged. A read-only independent verification may audit the checks, following the spin-dynamics skill.

### Verification and Acceptance Criteria

Require finite occupied spins, postprojection norm error below 1e-12, exact-zero empty entries, aligned physical times, and protected source/result hashes. For matched coupled dynamics, target mean-spin differences below 1e-5 under timestep halving and final-angle agreement within 0.1 degree. Compare with the reference's recorded 30-degree destinations at 175–275 Oe, 90 degrees at 300–500 Oe, and 90 degrees at 1000 Oe for both current signs. Report the sampled transition bracket rather than claim a precise critical field. A failing native/frozen trajectory remains visible; qualitative endpoint agreement alone does not validate its integration.

### References

[1] Zhengde Xu, Xue Zhang, Yixiao Qiao, Gengchiau Liang, Shuyuan Shi, Zhifeng Zhu, Deterministic spin-orbit torque switching including the interplay between spin polarization and kagome plane in Mn3Sn, Physical Review B **109**, 134433 (2024). DOI: https://doi.org/10.1103/PhysRevB.109.134433. Evidence already inspected before the environment change: reviewed paper record `AFM/2024-Deterministic spin-orbit torque switching including the interplay between spin-PRB-zhengde Xu`, source PDF model page 134433-2 and Fig. 5 page 134433-5, plus text extraction of pp. 2–5. The presentation record `notes/ee237/week15a-afm/slides.md` was checked for context; the paper supplies quantitative inputs. No further access to the denied papers directory is required.

### Assumptions, Missing References, Risks, and Approval Questions

- The user's “proceed your suggestion” and “continue” authorize the isolated diagnostic comparison described above. Production changes and larger-lattice repairs remain outside scope.
- Existing historical sections are retained. No exact digitized published curves are available; distinguish numerical agreement with the reference implementation from qualitative agreement with the published plot.
- Physical crystallographic mappings remain REFERENCE MISSING, as recorded in the original audit. They do not block the algebraic three-vector comparison and must not be inferred from it.
- Execution completed on 2026-11-02 at the user's request. Fresh source snapshots and outputs are under `.sisyphus/three-spin-20261102-02/`, with maintained diagnostic harnesses under `verification/`. Production source and existing results remain protected. Results and rerun commands are in [the comparison report](docs/three-spin-comparison.md).

### Completed comparison

- **Mapping/static gates passed.** The 2×2 source setup gives exactly the proposed three-spin mapping, three reciprocal exchange pairs, and cyclic antisymmetric DMI coefficients. Matched fields differ from `field_eta.m` by at most $2.84\times10^{-14}$ T; scaled energy-gradient residuals are below $2.23\times10^{-11}$. Native/reference anisotropy differs by two, SOT by 449.997115592, and exchange/gamma by the small source rounding differences tabulated in the report.
- **RHS/replay gates passed.** Explicit/source RHS and implicit Gilbert residuals are below $3.6\times10^{-16}$. Twelve startup and eight driven 5 ps MATLAB comparisons have maximum spin discrepancy $1.75\times10^{-13}$. A scalar GPU dispatch/reset check agrees with the CPU and preserves empty output slots. This is not full production GPU integration.
- **Coupled convergence and Figure 5 gates passed.** 336 full 5 ns runs cover all 24 physical cases and four variants; six additional frozen-field refinements bring the full-length total to 342. At 5/2.5/1.25 fs, both coupled variants give 30° for 175–275 Oe, 90° for 300–500 Oe, and 90° for all six seeds/current signs at 1000 Oe. Maximum coupled mean-component differences under halving are $2.14\times10^{-9}$ and $1.34\times10^{-10}$, approximately sixteenfold reduction. The sampled transition bracket is 275 < H ≤ 300 Oe.
- **Frozen-field accuracy failed, as a retained diagnostic result.** Matched frozen fields at 0.5 fs give 89.9338° at 275 Oe rather than 30°. Even at 0.0625 fs, the 300 Oe matched frozen transient differs from the reference by 0.00108219 in a mean component. Native/frozen trajectories are also unconverged in the tested range and combine coefficient and integration differences. Near-machine-precision norms do not establish trajectory accuracy.
- **Scope and artifacts.** All outputs use update-derived physical times, exact 5 ns endpoints, and 0.25 ps sampling. Full-length calculations use the validated compiled replay; both production entry scripts remain unexecuted. MATLAB R2024b Update 1, GCC 11.4.0, Python 3.10.12 standard library, and the RTX 3060 Laptop GPU scalar probe are recorded in the report. No dependency installation, network access, new physical lattice, production repair, or publication-curve digitization occurred. The source-hash preservation check and independently reviewed plots complete the diagnostic; larger-lattice reference gaps remain open.


## Frozen-field repair and Figures 2, 4(a), 5 — 2026-11-02

### Goal and Success Criteria

The user's explicit request authorizes fixing the frozen-field issue and checking Figures 2, 4 and 5 against `Mn3Sn_3atoms`, excluding the ferromagnet. Preserve historical results and reference sources. Pass requires production stage-state field evaluation, coupled fourth-order convergence, exact-zero vacancies, unchanged LLGS coefficients, and measured agreement with the reference figure cases. This supersedes the earlier audit-only restriction for this specific solver repair.

### Model Classification and Equation Family

Class 4 spatial fixed-length LLGS; the verification reduction is the already defined Class 5 three-vector noncollinear Mn model [1, Sec. II]. No electronic-structure calculation or new crystal construction.

### Existing Behavior and Reference Implementations

`integrate_llg.m` builds fields once, and `atomgpurk4.m` advances each spin with those fields frozen. The earlier comparison above establishes the failure and reference correspondence. Consulted the reference repository's settings and corrected-RK4 report, existing static/replay harnesses, and ASTRO instructions/README for scope only. Preserve the legacy scalar kernel for historical diagnostics and unsupported legacy solver paths. The active RK4 path will use a scalar derivative kernel containing its unchanged first-stage algebra.

### State Variables, Geometry, and Boundaries

Arrays remain `(W,L,time)`; occupied sites 1, empty sites 3. Keep the source topology and coefficients exactly. Three-site bookkeeping uses `[4,1,2]` in column-major 2×2 storage, as verified above.

| Physical detail | Adopted rule / mapping | Evidence and status |
|---|---|---|
| Three Mn vectors and yz easy axes | Existing reduction and permutation above | [1, Sec. II]; existing measured static correspondence |
| Spatial coordinates, shells and periodic seams | No changes to `distrib.m` or interaction masks | REFERENCE MISSING for physical construction, as in the audit; blocks physical topology claims, not this solver repair |
| Boundary masks | Execute the same source at every stage | Software-preservation choice; no new physical claim |

### Hamiltonian

Unchanged. The reference triangle energy and native factor-of-two anisotropy distinction are given in the preceding comparison. Do not repair nonreciprocal larger-lattice weights as part of this change.

### Parameters and Evidence

Production parameters remain unchanged. Figure checks explicitly use the preceding reference parameter table and `Mn3Sn_3atoms/FIG2_FIG4_SETTINGS.md`: alpha 0.003, exchange field 100.812 T, K 0.196 meV with field factor 2, D 0.833 meV, moment 3 mu_B, Ms 441850 A/m, thickness 30 nm, thetaSH 0.06, gamma 1.760859644e11 rad/(s T). These are the reference-code realization of [1, Sec. II]; they are not new production defaults.

### Effective Fields and Units Check

Unchanged exchange, DMI and anisotropy expressions in tesla are evaluated from each complete unnormalized stage state. Keep external fields and torque amplitudes constant within a step. Noise, when selected, retains one draw per outer step; no stochastic convergence claim. Preserve the moment factors and energy-gradient checks of the preceding comparison.

### Dynamical Equation and Nonconservative Torques

Use precisely the existing scalar LLGS RHS, including both SOT/STT damping-dependent cross terms and field-like terms, matching `sections/llg_equation.tex` equations `eq:llg-gilbert` and `eq:ll-with-flt`. Extract the first-stage derivative without normalization or timestep factors; integration owns those factors. Validate all four torque amplitudes against the implicit Gilbert equation on unit spins.

### Initial Conditions, Solver, and Outputs

Classical RK4 stage states are m, m+dt*k1/2, m+dt*k2/2, m+dt*k3; recompute fields at every stage, combine with weights 1,2,2,1, normalize only the final occupied spins. Preserve empty entries at every stage. Preserve current chunking/time-label behavior and document its existing limitation; figure checks use exact update-derived physical times and aligned pulse edges.

Figures: 24 low-field combinations (six seeds, both current signs, both 100 Oe field signs), including Fig. 4(a)'s seed 210; the 24 Figure 5 cases already listed above. Current magnitude 5e10 A/m², drives 1–4 ns, total 5 ns, dt 5/2.5/1.25 fs, output interval 0.25 ps. No ferromagnetic calculation. Fresh `.sisyphus/frozen-field-*` directories protect production outputs. Full sweeps may use the existing compiled replay only after validation against the changed production MATLAB loop, including startup and driven states.

### Implementation Mapping

| Stage | Implementation | Verification gate |
|---|---|---|
| Structure/topology | Unchanged source and static checker | Three-site mapping, impulse fields and energy gradients pass |
| RHS | New scalar `atom_llgs_rhs.m`, unchanged original stage-one algebra | All torque signs/cross terms and implicit Gilbert residual <1e-12 |
| Integration | Four field/RHS passes in `integrate_llg.m`; occupied normalization | Actual production-loop replay <1e-10; vacancy and norm checks; short larger-array run; fourth-order convergence |
| Figures | Extend isolated benchmark tooling, reuse verified C++ replay | All 48 destinations within 0.1 degree; mean trajectory halving errors <1e-5; Fig. 4(a) crossing shift <1 ps |
| Documentation/viewer provenance | README, report, plan, source-hash snapshot | Commands and limitations accurate; viewer source/impulse check passes |

Proceed only after each dependent gate passes. Preserve existing static excerpts so the visualization and historical checks continue to operate. No geometry viewer redesign is needed. A read-only independent verifier may review the production patch and evidence, as recommended by the skill.

### Verification and Acceptance Criteria

Use MATLAB R2024b and the existing GPU where available; record actual environment. Check finite states, occupied norms <1e-12, exact zero empty slots, preprojection errors, and order near four on a coupled deterministic case. Replay the actual production loop, not merely another corrected formula. Full figure checks require the source-matched compiled pilot to pass first. Protect source/result hashes in the reference checkout and all historical MAT files. No digitized paper curves are available, so quantitative claims concern reference-code agreement; published-figure agreement remains qualitative.

### References

[1] Zhengde Xu, Xue Zhang, Yixiao Qiao, Gengchiau Liang, Shuyuan Shi, Zhifeng Zhu, Deterministic spin-orbit torque switching including the interplay between spin polarization and kagome plane in Mn3Sn, Physical Review B **109**, 134433 (2024). DOI: https://doi.org/10.1103/PhysRevB.109.134433. Source locators and prior reviewed-PDF evidence are recorded in the preceding comparison. This task consults that existing evidence plus the reference settings/report and executable sources. A fresh targeted `_knowledge` filename search returned no accessible matching record; no new physical construction is inferred.

### Assumptions, Missing References, Risks, and Approval Questions

The explicit fix request authorizes this narrow repair and benchmark scope; no further approval is needed for it. Production defaults cannot reproduce these figures without the explicitly matched parameters. Larger-lattice nonreciprocity, physical geometry, legacy solver defects, output-time drift, thermal/dipole/fixed-edge validation remain separate limitations. No new dependencies, network access, reference-repository writes, commits, or result overwrites are authorized or needed.


### Completed repair and verification

Implemented the four complete-lattice field/RHS stages in `integrate_llg.m`, using `atom_llgs_rhs.m` with the original derivative algebra. Static fields, production default parameters, source line endings, and historical kernels remain unchanged. The visualization snapshot was refreshed and its 367,700 operator comparisons and source hashes pass.

The actual production loop passes startup and eight driven GPU replay checks (maximum spin discrepancy 7.38e-14), exact-zero vacancy/finite/norm checks, both legacy boundary flags on short 5×5 runs, all-torque Gilbert residual 3.20e-16, and a generic noncoplanar convergence test (orders 3.8016 and 3.9136). CPU convergence changes only transfer/scalar-expansion mechanics in the isolated harness. These checks do not validate the physical periodic seams or optional paths.

All 192 full diagnostic trajectories completed. All 48 physical cases match the expected destinations at 5/2.5/1.25 fs. Maximum mean-component differences under halving are 2.13797e-9 and 1.33636e-10 (about sixteenfold); Figure 4(a) crosses mean-z zero at 1.399001533608 ns; Figure 5's sampled transition bracket remains 275 < H ≤ 300 Oe. Full runs use the validated C++ diagnostic; short source checks execute actual MATLAB/GPU integration. The final-angle error is at most 6.63e-7 degree and every-update projected norm error at most 2.22e-16.

An additional stricter equal-dt per-spin comparison, introduced during harness construction, misses 1e-9 at 275 Oe: the discrepancy is 2.93e-9, while the mean-component difference is only 7.77e-12. Retain this as `strict_spin_match_pass=false`; the planned mean-observable and convergence gates pass. The compiled modes differ only in field arithmetic ordering. No Hamiltonian, parameter, or model choice was adjusted to force agreement.

Results, limitations, plots and complete commands are in [docs/frozen-field-repair.md](docs/frozen-field-repair.md). Raw artifacts use `.sisyphus/frozen-field-20261102-01/`. Independent read-only review found no blocking defect in the default deterministic repair or harnesses; fixed-edge stage constraints remain outside the validated scope and are explicitly documented. Reference sources/results and historical MAT files are protected by hashes. No network, dependency installation, commits or pushes occurred.


## Reference parameter alignment and explicit protocols — 2026-11-02

### Goal and Success Criteria

The user approved the preceding recommendation: adopt the three-atom conventions, explicit figure protocols and an undriven relaxation mode, while deferring unequal spatial exchange weights. Preserve the prior uncommitted RK4 repair and all existing results. Success requires native three-site parameters/fields to equal the reference without material overrides, correct finite pulse timing and saved-state times, verified figure presets, and unchanged spatial masks.

### Model Classification and Equation Family

Class 4 fixed-length spatial LLGS; figure presets use the established Class 5 three-vector reduction [1, Sec. II]. The equation, torque cross terms and normalization remain unchanged. Electronic structure is not modeled.

### Existing Behavior and Reference Implementations

Current source and the completed repair above are the starting point. `Mn3Sn_3atoms/main.m`, `field_eta.m`, `FIG2_FIG4_SETTINGS.md`, and `docs/rk4-verification.md` provide exact implementation values/protocols. Reuse the already inspected equation references and static/dynamical verification. No reference repository edits. Record the pre-task source under `.sisyphus/parameter-alignment-before/`.

### State Variables, Geometry, and Boundaries

| Physical detail | Adopted mapping | Evidence / status |
|---|---|---|
| Spatial lattice and interactions | Keep `distrib.m` and all exchange/DMI masks byte-for-byte unchanged | Existing audited implementation; physical geometry/shell mapping remains REFERENCE MISSING; explicitly deferred |
| Figure state | 2×2 open array, occupied reference order `[4,1,2]` | Software representation verified by earlier impulses; three-spin model [1, Sec. II] |
| Relaxation state | Existing 5×5 open array and 121-degree seed | Existing software example, no new physical-size claim |
| Film thickness vs spacing | `tz=30 nm`; retain `d=0.4 nm` | Film thickness [1, Sec. II]; d is unchanged legacy geometry, not newly certified |

### Hamiltonian

Adopt onsite energy `-Ksim*sum((e_i dot m_i)^2)` with Ksim=0.196 meV, giving field `2*Ksim/mu*(e_i dot m_i)*e_i`. This replaces the previous half-strength field while preserving the energy parameter's physical meaning [1, Sec. II; exact normalization confirmed by reference `hk6`]. Triangle exchange, DMI and Zeeman energy remain as specified above. Do not infer a reciprocal Hamiltonian for the unchanged larger-array weights.

### Parameters and Evidence

Use J=100.812*(3*mub) J, gamma=gam=1.760859644e11 rad/(s T), msRE=441850 A/m and tz=30e-9 m, matching reference source. Keep moment 3*mub, K=0.196 meV, D=0.833 meV, alpha=0.003, thetaSH=0.06, polarization +x, no STT/FLT/noise/dipoles. The exact exchange field is a code-matching value (17.50806 meV); the paper quotes nominal 17.53 meV [1, Sec. II]. Ms approximates 6*mu/Vcell; do not use one atomic moment per crystallographic cell as the film magnetization. Keep T=100 K inactive and other legacy inactive parameters unchanged.

### Effective Fields and Units Check

`B=-dE/dm / mu` gives tesla, with J/T moments and J energies. Derive the anisotropy factor two explicitly above. The torque amplitude remains `(hbar/e)*thetaSH*Jc/(2*Ms*tz)` in SI; the stored eV-second hbar already supplies the numerical charge cancellation. No extra charge factor. Check both the energy gradient and identical native/reference coefficients.

### Dynamical Equation and Nonconservative Torques

Unchanged explicit LLGS from the prior specification and `sections/llg_equation.tex`, `eq:ll-with-flt`. Multiply applied field and all existing torque amplitudes by a common per-step pulse gate for the selected figure protocol. Keep that gate fixed across stages; align pulse edges with step boundaries. No change to thermal-noise interpretation or other optional paths.

### Initial Conditions, Solver, and Outputs

Add a small `mn3sn_protocol(name,value)` function for explicit fixed presets. Default `relaxation`: 5×5, 121-degree seed, zero drive, 10 ps, 0.5 fs, 100 fs samples, 1 ps chunks. Figure presets: 2×2, 5 ns, 5 fs, 0.25 ps samples, 10 ps chunks, on at 1 ns/off at 4 ns. Fig2a/b/c use seed270 and signs (+,+),(+,-),(-,-) with 100 Oe; fig2d accepts the six seeds at (+,+); fig4a uses seed210; fig5a accepts the twelve tested fields with seed210; fig5b/c accept the six seeds at 1000 Oe and opposite current signs [1, Figs. 2,4,5]. A preset selects one trajectory, not an automatic full sweep. The ferromagnet stays excluded.

Correct the chunk accounting needed for exact pulses: each chunk stores its initial state plus all requested updates, retain each sampled boundary once, include t=0 and the final state, and derive output times from completed updates. This necessary timing repair changes default output from 100 to 101 frames and actual duration from 9.995 to 10 ps. Keep the legacy MAT save path and schema; tests never execute that save. Manual preset selection remains in main.m because it clears the workspace.

### Implementation Mapping

| Stage | Implementation | Verification gate |
|---|---|---|
| Structure | No mask/topology change; configurable initial angle only | Preserve distrib hash and normalized exchange/DMI impulse matrices |
| Material/fields | Main parameter constants and factor-two anisotropy field | Native/reference coefficients equal; independent finite-difference anisotropy and all static gates pass |
| Protocols/time | Named presets, integer pulse gate, correct chunk states/sample indices | All presets and seeds checked; exact t=0/final, no duplicate boundaries, chunk-invariance and pulse-edge replays |
| Regression | Adapt existing checks to K as energy coefficient; use native parameters | Actual source GPU replays plus inherited figure diagnostic convergence; no material overrides masking a mismatch |
| Views/docs | Refresh factor-two display and source hashes; README/AGENTS/changelog | Node static checks, MATLAB export, accurate commands and explicit spatial limitation |

Keep field evaluation, scalar RHS and orchestration separate. Advance only after each dependent gate passes. A read-only verifier may review the changes per the spin-dynamics skill.

### Verification and Acceptance Criteria

Reuse source-based triangle impulses/energy-gradient checks (<1e-11 T and <1e-7 scaled residual), all-torque Gilbert checks (<1e-12), GPU/source replay (<1e-10 spin difference), norms (<1e-12), exact vacancies and coupled convergence. Add deterministic multi-chunk/pulse tests against a single chunk and reference gate scheduling. Validate every named preset/choice, and audit all 48 figure cases against native coefficients and the existing full-trajectory diagnostic. Reuse completed figure evidence when equivalent parameters/equations are demonstrated; do not claim full 5 ns GPU runs. Snapshot/hash-check existing results and the reference checkout. No installs or network.

### References

[1] Zhengde Xu, Xue Zhang, Yixiao Qiao, Gengchiau Liang, Shuyuan Shi, Zhifeng Zhu, Deterministic spin-orbit torque switching including the interplay between spin polarization and kagome plane in Mn3Sn, Physical Review B **109**, 134433 (2024). DOI: https://doi.org/10.1103/PhysRevB.109.134433. Existing reviewed-paper/source evidence and exact local record are recorded above; this step uses the verified reference source and figure settings. No new physical geometry is introduced.

### Assumptions, Missing References, Risks, and Approval Questions

The user's latest message approves these changes and explicitly postpones spatial exchange-weight repair. Larger-lattice reciprocity, geometric evidence, dipoles, noise and fixed edges remain unvalidated. Presets use original ideal 120-degree seeds; startup ringing is retained. Full GPU figure runs are potentially slow; the verified compiled diagnostic remains the regression sweep backend. Do not repeat approval for the authorized changes; update the user about the necessary clock correction. No commits or pushes.


### Alignment results and implementation checkpoints

The approved native parameters now equal the reference's stored values for mu, J, D, Kfield=2*Ksim, gamma, alpha, Ms, thickness and SOT. Actual-source static fields and independent energy derivatives pass (2.84e-14 T and <2.23e-11 scaled residual). The 2K factor is explicit in the field; Ksim retains its 0.196 meV energy meaning.

`mn3sn_protocol.m` exposes one relaxation case and named Fig. 2/4(a)/5 cases. Actual main initialization passed all 35 valid preset/choice combinations and rejected unsupported choices, including the ferromagnet. Eight production GPU pulse fixtures cover all drive-sign combinations and edges inside/on chunk boundaries: independent-reference discrepancy <=1.05e-15 and single-vs-multiple-chunk difference zero. Irregular sampling, unique boundaries, t=0/final samples and exact update counts pass. Default output is now 101 states at 0–10 ps; figure presets produce 20,001 states at 0–5 ns.

The native-parameter production convergence orders remain 3.8016/3.9136, and startup GPU replay discrepancies are <=4.43e-14. All 192 full compiled figure trajectories are bit-for-bit identical to the preceding explicitly matched figure runs. All destinations and refinement gates pass; the recorded strict auxiliary per-spin comparison exception at 275 Oe is unchanged. See [docs/reference-protocols.md](docs/reference-protocols.md) for commands, results and limits.

The source/display impulse checks pass across ten array sizes and both boundary flags, with 367,700 exchange/DMI operator comparisons. Browser inspection confirmed the updated 2K/mu display. Exact source comparisons preserve distrib.m, exchange masks, shifts, exchange field and DMI expressions; no spatial-weight repair was made. Independent read-only review found no blocking defect in the physics, presets, pulse indexing or chunk tests. All outputs use `.sisyphus/parameter-alignment-20261102-01/`; historical results and reference source remain protected.

Final driven GPU replays passed (maximum spin discrepancy 7.38e-14), and the source/result preservation check passed for all 30 recorded paths. Whitespace checks passed for tracked diffs and new text files. No production output or historical result was overwritten; no commit or push was made.
