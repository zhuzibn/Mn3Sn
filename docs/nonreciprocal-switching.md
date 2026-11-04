# Nonreciprocal switching comparison

**The existing nonreciprocal method reproduces the reference switching curves for the balanced sublattice mean on all three spatial sizes.** Each size matches all 48 figure cases, including the 275–300 Oe sampled transition and 1.399 ns z reversal. The raw site-weighted magnetization differs because the finite patches contain excess b sites. This establishes numerical equivalence in the tested coherent subspace and perturbed-start recovery, not physical validation of the asymmetric bonds.

This diagnostic tests the existing receiver-weighted exchange **and DMI**, with uniform anisotropy, against the same PRB 2024 Figures 2, 4(a), and 5 used in the [unequal-anisotropy comparison](anisotropy-switching.md). Production physics and defaults are unchanged. The current branch remains `experiment/unequal-anisotropy-prb2024`; no branch, commit or push was created for this follow-up.

## Model and reason to expect coherent agreement

The input graph and directed coefficient matrices are the preserved MATLAB-source exports. The three spatial graphs have 16, 130 and 336 occupied sites, with a,b,c populations (5,6,5), (42,46,42) and (110,116,110). Each site is integrated independently; the solver does not impose sublattice coherence.

For dimensionless exchange matrix A and DMI matrix C, the sums into the three sublattices are exactly:

| Receiving class | Exchange sums into a,b,c | DMI sums into a,b,c |
|---|---|---|
| a | 0, 1, 1 | 0, −1, 1 |
| b | 1, 0, 1 | 1, 0, −1 |
| c | 1, 1, 0 | −1, 1, 0 |

These sums hold at every occupied site, including edges, and equal the three-spin control. Consequently, when spins are identical within each sublattice, the internal fields equal the control fields. With uniform K, damping, moment, external field and SOT, the complete driven LLGS equations preserve that coherent subspace. This is an exact algebraic property of these code graphs, independently checked from the exported coefficients. It is not an independent justification of a physical nonreciprocal interaction.

The fields are

$$
\mathbf B_i=-\frac J\mu\sum_j A_{ij}\mathbf m_j+\frac D\mu\sum_j C_{ij}(0,m_{jz},-m_{jy})+\frac{2K}{\mu}(\mathbf e_i\cdot\mathbf m_i)\mathbf e_i+\mathbf B_{\rm ext}.
$$

The existing LLGS RHS retains the damping-dependent SOT cross term. The existing RK4 kernel recomputes all fields on each unnormalized coupled stage, then normalizes occupied spins after combination. No solver or physical coefficient was changed.

The spatial graphs contain respectively 12, 44 and 76 unordered pairs with unequal reverse exchange weights; the same counts have nonantisymmetric DMI weights. The maximum normalized defect is 0.5. With equal physical moments, the full operator is not the gradient of the usual reciprocal pair Hamiltonian. The legacy output column `energy_over_mu_per_spin` is therefore **not used as a physical energy** in this comparison.

## Observables and interpretation

The paper-compatible three-sublattice observable is the unnormalized balanced mean,

$$
\mathbf m_{\rm bal}=\frac{\overline{\mathbf m}_a+\overline{\mathbf m}_b+\overline{\mathbf m}_c}{3}.
$$

The equal-moment finite patch has site-weighted mean $\mathbf m_{\rm site}=N^{-1}\sum_i\mathbf m_i$. For populations $(n,n+\delta,n)$,

$$
\mathbf m_{\rm site}=\frac{3n}{N}\mathbf m_{\rm bal}+\frac{\delta}{N}\overline{\mathbf m}_b.
$$

The excess-b fractions are 0.0625, 0.030769 and 0.017857. They exceed the reference weak-moment magnitude, about 0.00346, so matching balanced curves does not establish matching raw finite-patch magnetization. Both means and all three sublattice means are saved, together with edge/interior means and maximum within-sublattice vector spread. Neither mean is normalized before comparison.

Angles are atan2(mean_y,mean_z), with circular comparison modulo 360 degrees. At exactly compensated starts the angle is undefined. The optional plotting function omits the first 0.05 ns from angle panels and masks in-plane magnitudes below 1e-10; component panels retain startup. No smoothing is used.

## Protocol and verification scope

Use the [existing verified reference protocol](reference-protocols.md): $J/\mu=100.812$ T, $D=0.833$ meV, $K=0.196$ meV, $\mu=3\mu_B$, $\alpha=0.003$, $\gamma=1.760859644\times10^{11}$ rad/(s T), $M_s=441850$ A/m, thickness 30 nm, spin Hall angle 0.06 and polarization +x. The implemented exchange is 17.50806 meV; retain its documented difference from the paper's nominal 17.53 meV. No fitting or parameter tuning is performed.

Each trajectory lasts 5 ns, with 5 fs steps and external field/current active during [1,4) ns. The 48-case matrix per model contains 24 low-field sign/seed combinations, twelve fields from 175 to 500 Oe and twelve +1000 Oe cases spanning six seeds and both current signs. Figure 4(a) reuses the +current,+100 Oe,210-degree first-spin seed. The field sweep samples 175:25:400,450,500 Oe. Figure 4(b)'s ferromagnet is excluded. Noise, dipoles, pinning, STT and field-like SOT remain disabled.

Initial states are ideal 120-degree seeds, with first-spin angles 30:60:330 degrees and exactly zero x components in the diagnostic. Startup relaxation is retained; no pre-relaxed state is silently substituted.

Eight representative protocols per model have 2.5 fs refinements. Four protocols per spatial size additionally start with the existing deterministic site-dependent perturbation of amplitude 0.001, at both timesteps. These perturbations are injected at t=0, allowing 1 ns of relaxation before the pulse. They test **perturbed-start recovery**, not an ensemble of thermal fluctuations or general stability to disturbances injected during switching.

All full trajectories use the unchanged C++ diagnostic. Fresh MATLAB/GPU execution is unavailable in this environment. Validation replays eight preserved current-mode undriven GPU fixtures and twenty preserved driven GPU fixtures against a newly compiled binary. The latter contain sixteen spatial reciprocal-candidate cases and four current-mode three-spin cases; they do not constitute fresh nonreciprocal spatial driven GPU runs. Kernel SHA-256 and production-source hashes match the prior validated run. A separate independent Python implementation, using vector triple-product identities, checks sixteen short nonreciprocal driven trajectories with both pulse edges and noncoplanar initial states. Reciprocal candidate energy-gradient checks invoked by shared preparation are not evidence of a nonreciprocal energy.

Long-run checks require finite states, occupied norm error below 1e-12, exact clocks/update counts, and full-spin refinement errors below 1e-5. Source fixtures retain exact-zero vacancies. Endpoint agreement uses a 0.1-degree circular tolerance. Full transient mean differences and retained final z are assessed separately. No digitized paper curves are fitted: quantitative claims compare the numerical reference, while agreement with published figures remains qualitative.

## Results and artifacts

**248 full trajectories completed:** 192 primary runs (48 cases across four models), 32 additional pilot refinements and 24 perturbed-start/refinement runs. All numerical gates passed. An independent verifier recomputed all 192 primary and 24 perturbed endpoints directly from raw final spins, all twelve perturbed full-spin refinement pairs, and the dense Fig. 4 traces.

| Model | Low-field matches | Field-sweep matches | High-field matches | Maximum transient balanced-mean component difference from numerical reference |
|---|---:|---:|---:|---:|
| Full-strength three-spin control | 24/24 | 12/12 | 12/12 | 7.75e-12 |
| 16 spins | 24/24 | 12/12 | 12/12 | 7.71e-12 |
| 130 spins | 24/24 | 12/12 | 12/12 | 7.71e-12 |
| 336 spins | 24/24 | 12/12 | 12/12 | 7.72e-12 |

Figure 2 retains the four current/field destination combinations from all six starts: (+J,+H) → 30 degrees; (+J,−H) → 210 degrees; (−J,+H) → 150 degrees; (−J,−H) → 330 degrees. Figure 5(a) ends at 30 degrees for 175–275 Oe and 90 degrees for 300–500 Oe on every size. Figures 5(b,c) end at 90 degrees for all six starts and either current sign at +1000 Oe. These reproduce the reference's current-dependent low-field selection and sampled transition, both absent in the unequal-anisotropy spatial candidate. The largest endpoint error across primary and perturbed runs is 6.63e-7 degrees.

For Figure 4(a), every balanced mean changes from approximately −0.00299732242 to +0.00299732243. The negative-to-positive zero crossing is **1.3990015336083 ns**, consistent across all sizes to numerical roundoff. The patch's physical site mean has the opposite sign in these finite, unbalanced graphs:

| Occupied spins | Balanced z at 0.9 ns | Balanced z at 5 ns | Site-weighted z at 0.9 ns | Site-weighted z at 5 ns |
|---:|---:|---:|---:|---:|
| 16 | −0.00299732 | +0.00299732 | +0.05112862 | −0.05112862 |
| 130 | −0.00299732 | +0.00299732 | +0.02364929 | −0.02364930 |
| 336 | −0.00299732 | +0.00299732 | +0.01246723 | −0.01246723 |

The JSON `balanced_crossing_ns` and `site_crossing_ns` fields scan specifically for a negative-to-positive crossing after 1 ns. A null site crossing does not mean no reversal: the site mean above reverses in the opposite direction. Do not identify its curve with the paper's three-sublattice mean.

Uniform-start maximum within-sublattice spread is 3.03e-14, so these trajectories remain coherent without an explicit coherence constraint. The perturbed-start cases reach a transient spread as large as 0.00691 and finish with spread below 2.59e-14 across both timesteps. At 5 fs, their final states agree with the corresponding unperturbed states within 7.64e-14 per component. This recovery result retains the t=0 perturbation limitation described above.

| Numerical check | Measured maximum |
|---|---:|
| Eight archived current-mode undriven GPU replays | 3.34e-16 per component |
| Twenty archived driven GPU replays | 4.05e-16 per component |
| Sixteen independent driven CPU equation checks | 3.34e-16 per component |
| 32 pilot full-spin comparisons, 5 vs 2.5 fs | 1.54e-8 |
| Twelve perturbed full-spin timestep comparisons | 1.92e-7 |
| Postprojection occupied norm error | 2.23e-16 |
| Preprojection norm error, sweep trajectories | 1.16e-12 |
| Saved means versus selected raw full-spin snapshots | 1.37e-15 per component |

Source/result preservation passed. Fresh environment: GNU/Linux, Python 3.10.12 standard library and GNU C++ 11.4.0; exact compiler command and versions are saved in `compiler.json` and `provenance.json`. The archived source fixtures were produced with MATLAB R2024b Update 1 and the RTX 3060 Laptop GPU, as documented in the earlier reports. They are not fresh hardware measurements from this run.

Data and logs are saved outside OneDrive:

```text
D:\Git_control\spin_dynamics\simulation_results\Mn3Sn\nonreciprocal-switching-20261104-01
/mnt/d/Git_control/spin_dynamics/simulation_results/Mn3Sn/nonreciprocal-switching-20261104-01
```

| Artifact | Content |
|---|---|
| `figure-results.json`, `.csv`, `figure-summary.json` | Per-case and grouped results |
| `pilot-checks.json`, `coherent-reduction.json` | Refinement and exact sublattice-sum checks |
| `validation-provenance.json`, `driven-gates.json` | Archived fixture origins, hashes and replay errors |
| `independent-short-replay.json`, `independent-verification.json` | Independent equation and output checks |
| `*-observables.csv` | Means/spread every 0.25 ps |
| `*.bin`, `*-final.bin`, `*-snapshot-times.csv` | Full spins and exact sample times |
| `*-run.json`, `compiler.json`, `provenance.json` | Commands, parameters and source/compiler provenance |
| `register-matlab.sh` | User-requested WSL MATLAB registration procedure |

Full-spin snapshots occur every 5 ps, plus denser samples around pulse boundaries, and include t=0 and 5 ns. Reference datasets, `finalnew.mat` and `atomtypee.mat` are protected from writes.

## Reproduction and MATLAB registration

Run from the repository root and choose a fresh external directory. This runner uses standard-library Python and the existing GNU C++ compiler; no packages were installed.

```bash
python3 verification/spatial_relaxation/nonreciprocal_figures.py OUTPUT prepare
matlab -wait -batch "addpath('verification/spatial_relaxation'); prepare_driven('.sisyphus/spatial-relaxation-20261102-01','WINDOWS_OUTPUT','current')"
python3 verification/spatial_relaxation/nonreciprocal_figures.py OUTPUT check
python3 verification/spatial_relaxation/nonreciprocal_figures.py OUTPUT pilot
python3 verification/spatial_relaxation/nonreciprocal_figures.py OUTPUT pilot-check
python3 verification/spatial_relaxation/nonreciprocal_figures.py OUTPUT sweep
python3 verification/spatial_relaxation/nonreciprocal_figures.py OUTPUT analyze
python3 verification/spatial_relaxation/nonreciprocal_figures.py OUTPUT preserve
```

For the exact fallback used here, replace the MATLAB preparation and `check` commands with `python3 verification/spatial_relaxation/nonreciprocal_figures.py OUTPUT archived-check`. This explicitly depends on the preserved datasets and refuses a changed kernel/source hash.

`prepare_driven.m` now accepts the optional third argument `'current'`, preserving original receiver masks on spatial patches. Its default reciprocal-candidate behavior is retained. The new selection is reviewed but **not executed in MATLAB here**.

Optional saved-result plotting:

```matlab
addpath('verification/spatial_relaxation');
plot_nonreciprocal('WINDOWS_OUTPUT',fullfile('WINDOWS_OUTPUT','figures'));
```

The plotting function has **not been executed or visually verified**, and no new MATLAB PNG is claimed. The supplied MATLAB location, `D:\Program Files\MATLAB\R2024b\bin`, is absent from this execution environment even outside the sandbox. Windows cmd/PowerShell are also absent. The user's MATLAB README instructions were supplied in chat after the file itself proved inaccessible. Two wrapper-registration attempts created `~/.local/bin/matlab` and the `~/.bashrc` PATH entry during their commands, but the wrapper was absent in the following tool session. The same-process launch failed with exit code 127 because `matlab.exe` was missing. Persistent host registration therefore remains unverified. Run the saved `register-matlab.sh` in local WSL to install the wrapper, add PATH without duplicating its line, and test `matlab -wait -batch "disp(version)"`.

## Sources and limitations

Zhengde Xu, Xue Zhang, Yixiao Qiao, Gengchiau Liang, Shuyuan Shi, Zhifeng Zhu, Deterministic spin-orbit torque switching including the interplay between spin polarization and kagome plane in Mn3Sn, Physical Review B **109**, 134433 (2024). DOI: [10.1103/PhysRevB.109.134433](https://doi.org/10.1103/PhysRevB.109.134433). This task reuses the repository's prior reviewed-source evidence (Sec. II and pp. 134433-2 through -5) and verified reference dataset. It does not claim a fresh PDF reading; the current permission policy denies the papers tree.

Consulted `presentations/equations/spintronic_methods.tex`, `sections/material_parameters.tex`, `sections/hamiltonian.tex`, `sections/llg_equation.tex` and `sections/atomistic_boundaries.tex` for units, field differentiation and the explicit Gilbert cross terms. The reviewed `presentation-knowledge/notes/ee237/week15a-afm/slides.md`, especially slides 25–27, supplies qualitative context. ASTRO AGENTS/README supplied workflow context only.

Physical nonreciprocal equal-moment bonds, crystallographic mapping, shell assignments and a finite-device interpretation remain **REFERENCE MISSING**. This is a Class 3 fixed-length LLGS graph diagnostic. Its numerical agreement cannot establish a Class 4 material-specific spatial model or electronic/Hall predictions.
