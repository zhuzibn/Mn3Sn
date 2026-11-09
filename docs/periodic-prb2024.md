# Periodic reciprocal model: PRB 2024 switching comparison

**Production update (2026-11-09):** this tested periodic model is now the default in the root source. The investigation below records the earlier isolated run and its preserved inputs; see [production promotion and verification](periodic-production.md) for the current behavior and commands.

**Yes: the reciprocal periodic model reproduces the verified switching behavior of Figures 2, 4(a), and 5.** All three sizes match all 48 primary cases, and complete ordinary magnetization trajectories agree with the numerical reference within 7.72e-12 per component. All 248 control, primary, refinement and perturbed-start trajectories pass. This isolated comparison preceded promotion of the periodic model; current production now implements its verified periodic rules directly.

The comparison targets Figures 2, 4(a), and 5 of Zhengde Xu, Xue Zhang, Yixiao Qiao, Gengchiau Liang, Shuyuan Shi and Zhifeng Zhu, “Deterministic spin-orbit torque switching including the interplay between spin polarization and kagome plane in Mn3Sn,” Physical Review B **109**, 134433 (2024), [DOI: 10.1103/PhysRevB.109.134433](https://doi.org/10.1103/PhysRevB.109.134433). Quantitative comparisons use the preserved, previously verified three-spin numerical reference. Published-figure agreement is qualitative; no curves were digitized or fitted. Figure 4(b)'s ferromagnetic control and other untested paper results are outside scope.

## Model and reason to expect agreement

Use the same 48-, 108- and 300-spin periodic graphs verified in the [six-state relaxation test](periodic-six-state.md): reciprocal J/2,D/2 bonds, uniform K, equal a,b,c populations, and independent spin dynamics at every occupied site. No fields, coefficients, solver source or topology were changed for this comparison. The existing kernel already supports uniform field and damping-like SOT pulses.

Each site's summed exchange and signed DMI couplings into a,b,c equal those of the full-strength three-spin reference:

$$
A_{\rm sub}=\begin{pmatrix}0&1&1\\1&0&1\\1&1&0\end{pmatrix},\qquad
C_{\rm sub}=\begin{pmatrix}0&-1&1\\1&0&-1\\-1&1&0\end{pmatrix}.
$$

When the spins are uniform within each class, every site's internal field therefore equals its corresponding reference spin's field. A uniform applied field and uniform SOT add identical terms to these identical spins. The complete driven LLGS equation and every coupled RK stage preserve that subspace. This is an algebraic reduction of the tested code graph, not merely an inference from six stable equilibrium states.

Equal populations additionally give

$$
\mathbf m_{\rm site}=\frac1N\sum_i\mathbf m_i
=\frac13(\overline{\mathbf m}_a+\overline{\mathbf m}_b+\overline{\mathbf m}_c)
=\mathbf m_{\rm bal}.
$$

Thus the ordinary magnetization can match the reference directly. This avoids the excess-b contribution that made the raw magnetization of earlier nonreciprocal open patches differ from their balanced observable.

The historical script name `nonreciprocal_figures.py` does not describe the supplied periodic model. Its `current` mode preserves the input A,C coefficients, which are explicitly reciprocal/antisymmetric here. It does not impose nonreciprocity. Periodic outputs contain zero-valued `edge_*` columns because no degree-two sites exist; those are empty-category placeholders, not a physical edge magnetization, and are excluded from interpretation.

## Protocol and numerical checks

Retain [the verified reference settings](reference-protocols.md): J/mu=100.812 T, D=0.833 meV, K=0.196 meV, mu=3 mu_B, alpha=0.003, gamma=1.760859644e11 rad/(s T), spin Hall angle 0.06, M_s=441850 A/m and film thickness 30 nm. The SOT polarization is +x, with the existing damping-dependent Gilbert-to-LL cross term. At current +5e10 A/m2, its field-equivalent amplitude is 7.448364280185585e-5 T. The implemented exchange is about 17.50806 meV, retaining the reference's documented difference from the paper's rounded 17.53 meV.

Each run lasts 5 ns; current and field are active during [1,4) ns and held constant across each RK update. Initial first-spin angles are 30:60:330 degrees, with the other classes offset by 120 and 240 degrees. These are ideal 120-degree starts, retaining startup relaxation. Fields point along y. Thermal noise, dipoles, pinning, field-like SOT and STT remain disabled.

The 48 primary cases per model comprise 24 +/-100 Oe cases spanning both current signs and six starts, twelve fields 175:25:400,450,500 Oe at positive current, and twelve +1000 Oe cases spanning both signs and six starts. Figure 4(a) uses positive current, +100 Oe and first-spin seed 210 degrees. Means are unnormalized; angles are atan2(mean_y,mean_z), undefined at compensation.

The three-spin control plus three spatial sizes give 192 primary runs. Eight representative protocols per model are refined from 5 to 2.5 fs, adding 32 unique runs. Four selected protocols per spatial size additionally receive a deterministic site-dependent perturbation of amplitude 0.001 at both timesteps, adding 24 runs: **248 unique trajectories** in total. Perturbations occur at t=0, leaving 1 ns before pulse activation. They test recovery from perturbed starts, not disturbance injection during switching or a thermal ensemble.

The unchanged C++ diagnostic performs the full-length integrations. Exact coherent-field sums, model hashes and seeds were rechecked. Periodic source fields and undriven dynamics were freshly checked against MATLAB/GPU in the preceding relaxation task. Twenty archived driven MATLAB/GPU fixtures were replayed with the same kernel and passed within 4.05e-16/component. Those fixtures test the shared drive kernel on the earlier spatial candidate and three-spin control; they are not fresh periodic driven GPU executions. MATLAB environment provenance for the periodic fixtures is R2024b Update 1, RTX 3060 Laptop GPU, compute capability 8.6.

Pilot refinement must pass before the full sweep. Acceptance uses full-spin refinement error <1e-5/component, projected norm error <1e-12, exact integer update/pulse counts, aligned sample clocks, finite occupied states and endpoint error <0.1 degree. The analyzer also compares complete mean trajectories, the Figure 4 reversal and the sampled Figure 5 transition. Energy monotonicity is not required under nonconservative SOT. Mean sampling is 0.25 ps; full spins are saved every 5 ps with additional pulse-edge samples, including t=0 and t=5 ns.

## Results

| Model | Low-field cases | Field sweep | High-field cases | Maximum full-transient site-mean difference from numerical reference |
|---|---:|---:|---:|---:|
| Three-spin control | 24/24 | 12/12 | 12/12 | 7.75e-12 |
| 48-spin periodic cell | 24/24 | 12/12 | 12/12 | 7.71e-12 |
| 108-spin periodic cell | 24/24 | 12/12 | 12/12 | 7.72e-12 |
| 300-spin periodic cell | 24/24 | 12/12 | 12/12 | 7.70e-12 |

**Figure 2: current/field selection.** All six starts select the same reference destination for each sign combination:

| Current | Field along y | Final mean angle |
|---|---:|---:|
| +5e10 A/m2 | +100 Oe | 30° |
| +5e10 A/m2 | -100 Oe | 210° (-150°) |
| -5e10 A/m2 | +100 Oe | 150° |
| -5e10 A/m2 | -100 Oe | 330° (-30°) |

**Figure 4(a): reversal.** The ordinary site mean starts at $m_z\simeq-0.00299732242$ before the pulse and ends at $+0.00299732243$. Its negative-to-positive crossing occurs at **1.3990015336083 ns**, about **0.399 ns after pulse activation**, on every periodic size. The positive final sign is retained after pulse removal. Site and balanced means agree throughout within 1.46e-15 per component.

**Figure 5(a): field dependence.** Every tested 175–275 Oe case ends at 30°, and every tested 300–500 Oe case ends at 90°. The measured discrete transition bracket is **275 < H <= 300 Oe**; this is not a continuously resolved critical field.

**Figures 5(b,c): high field.** At +1000 Oe, all six starts end at 90° for both current signs, as in the reference.

The maximum endpoint error over primary and perturbed sweep cases is below 6.63e-7 degrees. Uniform-start maximum within-class spread is 2.28e-14, so the computed trajectories remain coherent without imposing a coherence constraint.

| Numerical check | Measured maximum |
|---|---:|
| 32 pilot full-spin comparisons, 5 versus 2.5 fs | 1.54e-8 per component |
| Twelve perturbed-start full-spin refinement comparisons | 1.39e-7 per component |
| Perturbed final state versus corresponding unperturbed state, both timesteps | 7.36e-14 per component |
| Occupied norm error | 2.23e-16 |
| Norm error before projection, primary/perturbed sweep | 9.43e-13 |
| Stored means versus selected raw snapshots, analyzer | 1.17e-15 per component |

The independent read-only audit covered all 248 unique runs, including all raw saved spin vectors, endpoints, exact seeds/perturbations, pulse/update counts, sample clocks and full-spin timestep comparisons. It also compared all 3,840,192 dense primary samples against the preserved numerical reference. Its ordinary-site-mean comparisons independently agree with the analyzer's balanced-mean conclusions.

The t=0 perturbations have already decayed to within-class spread below 1.54e-14 by pulse activation. Their success therefore demonstrates perturbed-start recovery, not robustness to a disturbance introduced during the drive.

Saved-data figures:

- [Figures 2, 4(a), field sweep and spatial spread](../../simulation_results/Mn3Sn/periodic-prb2024-20261109-01/periodic-prb2024.png).
- [Figures 5(b,c): both current signs and six starts](../../simulation_results/Mn3Sn/periodic-prb2024-20261109-01/periodic-prb2024-high-field.png).

Both figures overlay the three-spin reference where indicated. The periodic curves overlap the reference to numerical precision. Figure 4(a)'s site and balanced panels coincide because populations are equal. This agreement contrasts with the [open reciprocal Kz_i/4 trial](anisotropy-switching.md), which retained six undriven states but failed the driven low-field comparison.

MATLAB plot checks cover 51 series with 20,001 samples per time trace. Both exported PNGs were visually inspected for trace overlays, labels and clipping.

## Files and reproduction

External result directory:

```text
D:\Git_control\spin_dynamics\simulation_results\Mn3Sn\periodic-prb2024-20261109-01
/mnt/d/Git_control/spin_dynamics/simulation_results/Mn3Sn/periodic-prb2024-20261109-01
```

Existing execution functions are [driven_figures.py](../verification/spatial_relaxation/driven_figures.py) and [nonreciprocal_figures.py](../verification/spatial_relaxation/nonreciprocal_figures.py). Their source is unchanged. `invocation.json` records the module configuration and exact calls; every `*-run.json` records its executable command, parameters, elapsed time and kernel hash. `driven-gates.json` records the actual archived fixture replay commands and periodic validation provenance. Model copies and seed files are retained.

The MATLAB plotting copy `plot_periodic_prb.m` is stored in the external result directory, along with its source hash and configuration substitutions in `plot-provenance.json`. It adapts the existing saved-data plotting routine to the periodic sizes and adds high-field panels; it never reruns integration. Plots omit t<0.05 ns only from angle traces and mask in-plane means below 1e-10, unwrap angles continuously and apply no smoothing. Angles differing by 360 degrees represent the same orientation.

The following records the exact pre-promotion analysis invocation. Its source-preservation check intentionally rejects a changed checkout. Run it from the archived pre-promotion source snapshot to repeat that historical audit; do not replace its hashes with current ones. Current production checks are documented in [periodic-production.md](periodic-production.md).

```python
from pathlib import Path
import sys
sys.path.insert(0, 'verification/spatial_relaxation')
import driven_figures as d
import nonreciprocal_figures as n
out = Path('/mnt/d/Git_control/spin_dynamics/simulation_results/Mn3Sn/periodic-prb2024-20261109-01')
d.SIZES = ['2x2', '8x8', '12x12', '20x20']
d.verify_sources(out)
n.analyze(out)
d.preserve(out)
```

To repeat the same historical experiment from that archived source snapshot with the preserved validated inputs on this machine, copy only those inputs to a fresh directory before calling the existing functions:

```python
from pathlib import Path
import shutil, sys
sys.path.insert(0, 'verification/spatial_relaxation')
import driven_figures as d
import nonreciprocal_figures as n
old = Path('/mnt/d/Git_control/spin_dynamics/simulation_results/Mn3Sn/periodic-prb2024-20261109-01')
out = old.with_name('periodic-prb2024-new')
out.mkdir()  # Fail rather than overwrite an existing run.
names = ['spatial', 'source-hashes.json', 'compiler.json', 'protected.json', 'driven-gates.json']
files = [old/name for name in names] + list(old.glob('*.model')) + list(old.glob('*-seed*.bin'))
for path in files:
    shutil.copy2(path, out/path.name)
d.SIZES = ['2x2', '8x8', '12x12', '20x20']
d.verify_sources(out)
n.pilot(out)
n.pilot_check(out)
n.sweep(out)
n.analyze(out)
d.preserve(out)
```

This reuses validated model inputs and archived-fixture evidence; it does not claim a fresh source/GPU validation. Model/kernel changes require re-preparation and new gates. The executable is local-machine data, not a portable environment definition. Python uses only its standard library; no new dependencies were installed.

Use MATLAB to regenerate the saved-result plots:

```bash
matlab -wait -batch "out='D:/Git_control/spin_dynamics/simulation_results/Mn3Sn/periodic-prb2024-20261109-01'; addpath(out); plot_periodic_prb(out,out);"
```

Paper scope and parameter evidence are reused from the earlier reviewed-source comparison; this task also consulted the local `_knowledge/topics/mn3sn-spin-orbit-torque-switching.md` synthesis. The mathematical periodic continuation still has unresolved crystallographic mapping and physical exchange-shell normalization. Matching these coherent switching protocols does not validate finite-device edges, domain walls, arbitrary spatial instabilities, thermal retention or electronic Hall responses.
