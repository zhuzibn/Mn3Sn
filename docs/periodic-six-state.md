# Periodic boundaries recover six stable relaxation states

**Production update (2026-11-09):** this tested periodic model is now the default in the root source. The investigation below records the earlier isolated run and its preserved inputs; see [production promotion and verification](periodic-production.md) for the current behavior and commands.

**All six reference orientations survive on the tested periodic arrays with reciprocal J/2,D/2 bonds and uniform K.** The 48-, 108- and 300-spin cells relax to mean-spin angles 30°, 90°, 150°, 210°, 270° and 330°. The six final configurations are distinct, have equal energies to numerical precision, and recover after small spatial perturbations. On the largest cell, all twelve +/-5° orientation-offset starts also return to their corresponding state.

This is a deterministic, zero-drive test of the existing mathematical interior stencil continued periodically in both directions. It is not a crystallographic validation or an exhaustive search for metastable textures. At the time of this isolated experiment, production source and defaults were unchanged; they have since been promoted as described above. The [approved plan](../plan.md#proposed-periodic-six-state-relaxation-test--2026-11-09) records the scope and evidence limitations.

**Driven follow-up:** the [PRB 2024 switching comparison](periodic-prb2024.md) also passes for this same periodic reciprocal model with uniform K: all 48 reference cases per size, including Figures 2, 4(a), and 5, across 248 total trajectories.

## Results

| Array | Spins in a,b,c | Distinct final states | Final mean angles | Energy spread (meV/spin) |
|---|---|---:|---|---:|
| 8x8 | 16,16,16 | 6 | 30°,90°,150°,210°,270°,330° | 3.02e-14 |
| 12x12 | 36,36,36 | 6 | 30°,90°,150°,210°,270°,330° | 3.20e-14 |
| 20x20 | 100,100,100 | 6 | 30°,90°,150°,210°,270°,330° | 5.15e-14 |

The final energy is approximately **-9.573769075686 meV/spin**, matching the three-spin control. Across all 78 trajectories, the largest independently recomputed target-angle error is below 2.29e-10 degrees, terminal transverse field is 1.03e-12 T, and within-sublattice final spread is 2.86e-14. Distinctness and recovery use the full occupied-spin arrays, not rounded mean angles.

| Check | Measured maximum | Acceptance threshold |
|---|---:|---:|
| Perturbed/restarted/offset final state versus corresponding ideal result | 2.79e-12 per component | 1e-6 |
| 5 versus 2.5 fs, all common samples of six largest-cell perturbed runs | 1.38e-7 per component | 1e-5 |
| Occupied norm error across saved frames | 3.34e-16 | 1e-12 |
| Norm error before final RK projection | 1.15e-12 | Recorded; trajectory convergence checked separately |
| Largest sampled energy increase | 3.74e-13 meV/spin | 1e-9 |
| Recomputed mean versus stored mean / site versus balanced mean | 1.29e-15 per component | 1e-12 |

The evidence supports six locally attracting states in the tested cells and perturbation range. It does not establish that only six stable states exist, a global basin partition, an energy barrier, or finite-temperature retention.

The saved MATLAB figures are outside the source repository:

- [Relaxation angles, residual torques and spatial uniformity](../../simulation_results/Mn3Sn/periodic-six-state-20261109-01/periodic-relaxation.png).
- [Six final sublattice-spin configurations](../../simulation_results/Mn3Sn/periodic-six-state-20261109-01/periodic-six-states.png).

The trajectory figure shows the six spatially perturbed starts for each size. Mean-angle plots omit the first 0.05 ns and undefined in-plane means below 1e-10. Torque and spread retain the complete saved clock; values below 1e-15 are floored only for log plotting. All samples are 5 ps apart. No smoothing or normalization of mean spins is applied. The spin figure shows actual sublattice averages on schematic triangle positions, with the small mean multiplied by 180 for visibility.

## Why this differs from the open patches

The earlier [uniform-K reciprocal open-patch test](six-state-relaxation.md) had excess b sites. For the ideal 120° texture, its extra anisotropy contribution was

$$
E_{\mathrm{extra}}=-K\Delta N\sin^2\varphi,
$$

favoring two orientations. Complete periodic repeats have $N_a=N_b=N_c=N_0$, so this explicit imbalance term vanishes. The remaining ideal onsite sum is $-3KN_0/2$, independent of the common orientation. Canting under the existing exchange, DMI and onsite anisotropy then produces the six relaxed destinations observed numerically; no explicit sixth-order anisotropy or boundary-dependent K was added.

For every site, the periodic graph has two neighbors in each other sublattice. Its dimensionless coupling sums into a,b,c are exactly

$$
A_{\mathrm{sub}}=\begin{pmatrix}0&1&1\\1&0&1\\1&1&0\end{pmatrix},\qquad
C_{\mathrm{sub}}=\begin{pmatrix}0&-1&1\\1&0&-1\\-1&1&0\end{pmatrix}.
$$

Thus a spatially coherent state has the full-strength three-spin field at every site. This identity explains agreement of the coherent solutions. Independent spatial perturbations and recovery runs additionally test attraction back to that subspace. The observed maximum initial/trajectory class spread is about 0.00237, while final spreads are below 2.86e-14.

Equal populations also make the ordinary site mean equal to the balanced mean:

$$
\mathbf m_{\mathrm{site}}=\frac1N\sum_i\mathbf m_i
=\frac13(\overline{\mathbf m}_a+\overline{\mathbf m}_b+\overline{\mathbf m}_c).
$$

Angles use $\operatorname{atan2}(m_y,m_z)$, measured from +z toward +y. The ideal 120° seeds are compensated, so their initial mean angle is undefined. The six labels refer to the relaxed mean orientations.

## What was changed for the diagnostic

Before promotion, simply changing `bc=1` to `bc=0` was insufficient: `distrib.m` deleted terminal-column sites, `integrate_llg.m` cut exchange coefficients at the edges, and DMI followed a different mask path. Those historical cuts could leave exchange and DMI with different wrapped topology. The current periodic source has removed those cuts for bc=0.

The archived preparation helper applied the following **in-memory diagnostic overrides**, while retaining the production dynamics. The current [prepare_periodic.m](../verification/spatial_relaxation/prepare_periodic.m) executes production source directly because those periodic rules are now implemented there:

1. Remove the two terminal-column occupation exceptions from the repeating parity pattern.
2. Restore the terminal DMI selector and remove all fourteen exchange edge-cut assignments.
3. Use complete four-row/two-column repeats and periodic shifts in both directions.
4. Use uniform 1/2 exchange and DMI factors, while keeping uniform K unchanged.

The exported bond matrices satisfy $A_{ij}=A_{ji}$ and $C_{ij}=-C_{ji}$, including every seam. The isolated ordered-pair energy is

$$
E=\frac J2\sum_{ij}A_{ij}\mathbf m_i\cdot\mathbf m_j
-\frac D2\sum_{ij}C_{ij}\hat{\mathbf x}\cdot(\mathbf m_i\times\mathbf m_j)
-K\sum_i(\mathbf e_i\cdot\mathbf m_i)^2,
$$

with $A_{ij}=1/2$, $C_{ij}=\pm1/2$ on each directed bond. The outer 1/2 factors remove double counting, separately from the chosen half-strength bonds. The derivative field is $\mathbf B_i=-\mu^{-1}\partial E/\partial\mathbf m_i$ in tesla.

Parameters remain J/mu=100.812 T, D=0.833 meV, K=0.196 meV, mu=3 mu_B, alpha=0.003 and gamma=1.760859644e11 rad/(s T). No field, current, thermal noise, dipoles or pinned spins are enabled. This tests the same reciprocal mathematical normalization as the earlier open-array comparison; it does not resolve its physical exchange-shell assignment. L>=4 is required by production and this diagnostic: a two-column periodic array would collapse distinct image bonds onto the same source and require explicit multiplicity handling.

Inspect the verified data in the [periodic cell and interaction companion](periodic-cell.html). It shows cell repetition, class populations, source/receiver indices, neighbor offsets, image translations, exchange/DMI signs and the onsite term. Positions are schematic; physical lattice vectors, stacking, shell distances and crystal registry remain **REFERENCE MISSING**.

## Execution and verification

MATLAB R2024b Update 1, Parallel Computing Toolbox and an NVIDIA GeForce RTX 3060 Laptop GPU (compute capability 8.6) were used for fresh field and integration checks. Environment details and toolbox versions are in `environment.json`. Long trajectories use the unchanged `spatial.cpp` diagnostic compiled with `g++ -O3 -std=c++17 -Wall -Wextra`. Its `current` selector preserves the supplied coefficients; here the supplied matrices are explicitly reciprocal. This is not the old nonreciprocal spatial model.

Before the long sweep:

- Independently enumerate the wrapped stencil and compare every coefficient and occupied site with MATLAB impulse exports on 4x4, 8x8, 12x12 and 20x20 cells.
- Check equal populations, four distinct neighbors, two neighbors per other class, both seams, reciprocal exchange, antisymmetric DMI, matching supports, translation invariance and repeated-cell energy.
- Check termwise energy gradients, direct stencil versus matrix fields and the implicit Gilbert identity.
- Replay two fresh MATLAB/GPU trajectories per size plus two triangle controls: ten fixtures, 20 updates each, across two chunks. Check clocks, shapes, norms and exactly zero vacancies.

| Gate | Maximum discrepancy |
|---|---:|
| Direct source stencil versus matrix field | 5.69e-14 T |
| Independent termwise energy derivative | 1.59e-8 T |
| Implicit Gilbert residual divided by gamma | 3.36e-14 T |
| Energy per spin under cell repetition | 4.80e-14 T in E/mu units |
| C++ versus fresh MATLAB/GPU trajectory | 2.23e-16 per component |

The sweep contains 78 trajectories: 36 periodic primary runs (three sizes, six seeds, ideal/perturbed), six triangle controls, six largest-cell refinements, 18 perturbed restarts and twelve largest-cell +/-5° offsets. Primary/offset/refinement runs last 5 ns; restarts last 2 ns. Base step is 5 fs and refinement step is 2.5 fs. All states are independently evolved with coupled RK4; projection occurs after combining stages.

The first analysis attempt rejected 5000.000000000001 ps because it compared floating-point times with exact equality. The check now allows 1e-9 ps roundoff while still requiring the exact integer update count and each saved sample's clock. Analysis was rerun on the original completed trajectories; no numerical parameters or trajectories changed. The initial and final analysis sources/hashes are retained separately.

An independent read-only verifier checked all 67,278 saved configurations (11,998,026 spin vectors), initial-state construction, sample shapes/clocks, endpoint angles/fields/energies, configuration clustering, recovery and full-trajectory refinement. Unique-bond energies using accurate summation are degenerate within 3.56e-15 meV/spin; the slightly larger spread in the table reflects ordinary scalar accumulation. Intermediate energy monotonicity was checked from saved scalar energies, and preprojection errors remain integrator-recorded. Both MATLAB PNGs were visually inspected; 54 line series were reconciled with saved data. The companion passed 48 browser interaction cases and desktop/mobile layout inspection without console errors.

MATLAB Code Analyzer records 43 preparation notices, chiefly assignments/variables consumed through evaluated script source, plus growth and scalar-length advisories; the plotting function has two scalar-length advisories. All execution checks pass, but static analysis is not warning-free. The notice details are retained in `code-analyzer.json`. Python/JavaScript syntax, protected source/result hashes, whitespace and tracked line-ending checks pass.

## Files and rerun commands

Results directory: `D:\Git_control\spin_dynamics\simulation_results\Mn3Sn\periodic-six-state-20261109-01` (WSL: `/mnt/d/Git_control/spin_dynamics/simulation_results/Mn3Sn/periodic-six-state-20261109-01`).

Key artifacts are `gpu-gates.json`, `verified-gates.json`, `commands.json`, `refinement.json`, `periodic-summary.json`, `periodic-checks.json`, `periodic-results.mat`, `plot-checks.json`, per-run binary/CSV/JSON files, `.model` inputs, static MATLAB files, source snapshots, protected-file hashes and MATLAB logs. `periodic-results.mat` contains the summary, acceptance results and plotting checks; full trajectories remain in the documented binary format, one sample of N consecutive xyz triples at a time, little-endian float64. Their paired CSVs supply the exact saved times. Angles and means are also exported as `*-angles.csv`.

From the repository root, choose a new empty result directory (the commands below deliberately use a new name):

```bash
mkdir /mnt/d/Git_control/spin_dynamics/simulation_results/Mn3Sn/periodic-six-state-new
matlab -wait -batch "cd('D:/Git_control/spin_dynamics/Mn3Sn'); addpath('verification/spatial_relaxation'); prepare_periodic('D:/Git_control/spin_dynamics/simulation_results/Mn3Sn/periodic-six-state-new');"
python3 verification/spatial_relaxation/periodic_six_states.py /mnt/d/Git_control/spin_dynamics/simulation_results/Mn3Sn/periodic-six-state-new gates
python3 verification/spatial_relaxation/periodic_six_states.py /mnt/d/Git_control/spin_dynamics/simulation_results/Mn3Sn/periodic-six-state-new run
matlab -wait -batch "cd('D:/Git_control/spin_dynamics/Mn3Sn'); addpath('verification/spatial_relaxation'); plot_periodic('D:/Git_control/spin_dynamics/simulation_results/Mn3Sn/periodic-six-state-new');"
```

Python uses only its standard library. `analyze` can replace `run` to recompute diagnostics without integration. The sweep saves source hashes and refuses to overwrite partial trajectories; completed jobs can be reused only with unchanged runner sources and matching commands. Long runs begin with all six largest-cell refinement pairs and proceed to the remaining runs only after that numerical gate passes.

To refresh the local companion's data after successful preparation and gates:

```python
import json, hashlib
from pathlib import Path
run = Path('/mnt/d/Git_control/spin_dynamics/simulation_results/Mn3Sn/periodic-six-state-new')
data = json.loads((run/'periodic-view-data.json').read_text())
data['source_hashes'] = {p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in Path('.').glob('*.m')}
Path('docs/periodic-data.js').write_text(
    '/* Verified production graph; refresh command in periodic-six-state.md. */\n'
    'window.PERIODIC_DATA = '+json.dumps(data, separators=(',', ':'))+';\n')
```

Scientific source context remains the earlier reviewed Xu et al. PRB **109**, 134433 (2024) comparison documented in [reference-protocols.md](reference-protocols.md). This task consulted the local `_knowledge/topics/mn3sn-spin-orbit-torque-switching.md` synthesis and `presentation-knowledge/notes/ee237/week15a-afm/slides.md`, slides 25–27, for scope. Equation conventions were checked in `spintronic_methods.tex` and its material-parameter, Hamiltonian, LLGS and atomistic-boundary sections. No new source-PDF reading, physical lattice construction, literature reproduction or finite-temperature claim is made.
