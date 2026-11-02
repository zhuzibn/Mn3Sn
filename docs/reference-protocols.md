# Reference parameters and explicit figure protocols

This update adopts the approved `Mn3Sn_3atoms` material conventions in the
production source, adds explicit relaxation/figure settings, and fixes the
time accounting required to apply their pulses correctly. **Spatial exchange
receiver weights, neighbor masks and DMI expressions are unchanged.** Their
larger-array Hamiltonian correspondence remains deferred.

The target is Zhengde Xu, Xue Zhang, Yixiao Qiao, Gengchiau Liang, Shuyuan Shi,
Zhifeng Zhu, “Deterministic spin-orbit torque switching including the interplay
between spin polarization and kagome plane in Mn3Sn,” Physical Review B **109**,
134433 (2024), [DOI](https://doi.org/10.1103/PhysRevB.109.134433), Sec. II,
Fig. 2/Table I and Figs. 4–5. Exact numerical conventions and protocol settings
come from the verified local reference implementation and its
`FIG2_FIG4_SETTINGS.md` and `docs/rk4-verification.md`. The paper evidence
record, equation references and implementation stages are in [plan.md](../plan.md).

## Material changes

| Quantity | Before alignment | Current source |
|---|---:|---:|
| Anisotropy field coefficient | Ksim/μ = 1.128574509 T | 2 Ksim/μ = 2.257149019 T |
| Ksim | 0.196 meV | 0.196 meV, explicitly an energy coefficient |
| Exchange J/μ before spatial weights | 100.938106534 T | 100.812 T |
| Gamma, rad/(s T) | 1.75901203865e11 | 1.760859644e11 |
| Magnetization entering torque, A/m | 73642.1387 | 441850 |
| Film thickness tz | 0.4 nm | 30 nm |
| Legacy spatial spacing d | 0.4 nm | 0.4 nm, unchanged |
| Damping-like SOT at +5e10 A/m² | 0.0335174244 T | 0.0000744836428 T |

The anisotropy convention is now

$$
E_{\mathrm{ani}}=-K_{\mathrm{sim}}\sum_i(\mathbf e_i\cdot\mathbf m_i)^2,
\qquad
\mathbf B_{\mathrm{ani},i}=-\frac{1}{\mu_i}\frac{\partial E}{\partial\mathbf m_i}
=\frac{2K_{\mathrm{sim}}}{\mu_i}(\mathbf e_i\cdot\mathbf m_i)\mathbf e_i.
$$

The factor two is in `integrate_llg.m`, rather than hidden in a doubled energy
parameter. `gamRE` now uses `gam` from `constantfile.m`. The exact exchange
coefficient corresponds to 17.508061573 meV using the shared constants; the
paper's nominal 17.53 meV and the reference's rounded 100.812 T are slightly
different. This update deliberately matches the executable reference.

The magnetization 441850 A/m approximates six Mn moments per crystallographic
cell volume, instead of one moment per cell. Film thickness is independent of
the array-spacing variable. The SOT expression remains
`hbar/2*thetaSH*jcSOT/(msRE*tz)` with the stored eV-second `hbar`; there is no
extra elementary-charge division. The amplitude at equal current is reduced
by approximately 450. Default relaxation still has zero current and field.

Moment 3 μB, DMI 0.833 meV, damping 0.003, spin Hall angle 0.06, polarization
+x, the LLGS derivative and its damping cross terms are unchanged. Thermal
noise, dipoles, STT and field-like SOT remain off. Inactive T=100 K is retained.

## Selecting a calculation

Edit the `protocol=mn3sn_protocol(...)` assignment near the top of `main.m`.
The [README protocol table](../README.md#relaxation-and-paper-protocols)
lists every name, initial seed and field/current sign. `mn3sn_protocol`
returns one trajectory's settings; it does not execute a sweep.

- `relaxation`: existing 5×5 open array, 121° first-spin seed, zero drive,
  10 ps, 0.5 fs steps, 100 fs samples.
- `fig2a`, `fig2b`, `fig2c`, `fig2d`, `fig4a`, `fig5a`, `fig5b`, `fig5c`:
  three occupied spins in the verified 2×2 open reduction, 5 ns, 5 fs steps,
  0.25 ps samples, synchronized field/current pulses during 1–4 ns.
- `fig2d`, `fig5b`, `fig5c` take an optional seed from 30:60:330° (default
  210°). `fig5a` takes a field from 175:25:400, 450 or 500 Oe (default 250 Oe).
- Fig. 4(b)'s ferromagnet is excluded. The original ideal 120-degree seeds
  retain startup ringing; the calculation does not silently pre-relax them.

For example, `mn3sn_protocol('fig5a',300)` selects the positive-current
300 Oe case, and `mn3sn_protocol('fig5c',330)` selects the negative-current
1000 Oe case with first-spin seed 330°. `main.m` clears the workspace, so
setting a workspace variable before calling it is not an override mechanism.

The existing save behavior remains: `main.m` overwrites `finalnew.mat` in its
working directory. Use a separate copy of the root `.m` files for new runs
and edit its main assignment; the verification commands below do not execute
this save. No restart reader or output-path API was added.

## Physical clock and output compatibility

Each chunk now allocates its initial state **plus** the requested updates.
Shared chunk boundaries are sampled once. `saved_steps` contains zero-based
update indices; `t=saved_steps*tstep` is in seconds. The final state is included
even if the sampling stride does not divide the update count. The pulse gate
is on-inclusive/off-exclusive, evaluated from completed updates and held
constant through the four RK stages.

| Configuration | Updates | Saved states | Final time |
|---|---:|---:|---:|
| Previous default | 19,990 | 100, with incorrect labels | 9.995 ps |
| Current relaxation | 20,000 | 101 | 10 ps |
| Figure presets | 1,000,000 | 20,001 | 5 ns |

Downstream scripts must not hard-code the previous 100-frame default shape.
Use the saved `t` and array dimensions. Existing MAT files were not rewritten
and retain their historical conventions. The full-workspace MAT schema is
otherwise unchanged; the selected `protocol` is included in that workspace.

## Verification and reproduction

Use a fresh output directory under `.sisyphus`. From the repository root:

```bash
python3 verification/three_spin.py prepare .sisyphus/parameter-check-new
```

In MATLAB, from the repository root:

```matlab
addpath('verification');
out = fullfile(pwd,'.sisyphus','parameter-check-new');
check_three_spin(out);
check_protocols(out);
```

Then Bash:

```bash
python3 verification/three_spin.py probes .sisyphus/parameter-check-new
```

MATLAB:

```matlab
check_three_spin_dynamics(out);
check_coupled_rk4(out);
```

Bash after those gates:

```bash
python3 verification/coupled_figures.py sweep .sisyphus/parameter-check-new
python3 verification/coupled_figures.py analyze .sisyphus/parameter-check-new
python3 verification/coupled_figures.py replays .sisyphus/parameter-check-new
```

MATLAB:

```matlab
check_coupled_rk4(out,'driven');
plot_coupled_figures(out);
```

Finally Bash:

```bash
python3 verification/three_spin.py preserve .sisyphus/parameter-check-new
```

The static check independently reads current material assignments and requires
native/reference equality before its field and finite-difference energy checks.
`check_protocols` evaluates the actual `main.m` initialization without running
its save, checking all 35 preset/choice combinations. Its GPU integration
fixtures shorten the time grid, place pulse edges inside and on chunk
boundaries, exercise all four drive-sign combinations, and compare single vs
multiple chunks and irregular sampling against an independent coupled RHS.
The remaining production-loop checks use the native material values.

Compiled mode 2 now reads `native.params`; mode 3 retains `reference.params`.
All 48 figure cases are checked at 5/2.5/1.25 fs, with an additional 5 fs
reference trajectory per case. Full sweeps use the validated C++ diagnostic;
the actual MATLAB/GPU production loop is checked on short trajectories.
No full 5 ns GPU production run is claimed. The earlier
[figure report and plots](frozen-field-repair.md) explain the observables,
acceptance tolerances, and retained strict per-spin comparison exception.

The static viewer now shows 2K/μ and the corresponding −K energy. Refresh it
with the [README commands](../README.md#interactive-structure-and-interactions).
The source/impulse checks retain the unequal spatial bonds; they do not
silently make those interactions reciprocal.

## Results — 2026-11-02

Fresh artifacts are under `.sisyphus/parameter-alignment-20261102-01/`, with
the pre-task source preserved in `.sisyphus/parameter-alignment-before/`.
MATLAB R2024b Update 1 and the NVIDIA GeForce RTX 3060 Laptop GPU were used;
the compiled sweeps use the existing GCC 11.4.0 diagnostic and Python 3
standard library. No dependencies or reference-repository files were changed.

| Check | Measured result |
|---|---|
| Native/reference μ, J, D, anisotropy coefficient, gamma, alpha, Ms, thickness, SOT | Identical stored numerical values |
| Static field comparison | Maximum discrepancy 2.84e-14 T |
| Scaled energy-gradient residual | Below 2.23e-11 |
| Actual main initialization | All 35 preset/choice combinations passed |
| Eight pulse/chunk tests | Reference spin discrepancy ≤1.05e-15; zero measured single/multi-chunk difference |
| Irregular sampling and final sample | Exact match to the corresponding unsampled states/times |
| Coupled convergence fixture | Orders 3.8016 and 3.9136 |
| Production GPU startup replays | Maximum spin discrepancy ≤4.43e-14 |
| Eight production GPU driven replays | Maximum spin discrepancy ≤7.38e-14 |
| Full figure regression | 192 trajectories, 48 physical cases; all expected destinations retained |
| Figure timestep-halving mean errors | 2.13797e-9 and 1.33636e-10 |
| Every-update normalized spin error | ≤2.22e-16 |
| Static viewer/source checks | 367,700 exchange/DMI operator entries; all ten sizes and both boundary flags passed |

All 192 full trajectory files are **bit-for-bit identical** to the earlier
explicitly matched runs, confirmed by SHA-256 in `prior-figure-comparison.json`.
The existing [Figure 2, 4(a), and 5 plots](frozen-field-repair.md#measured-results--2026-11-02)
therefore remain valid; no new rendering is needed. Fig. 4(a)'s zero crossing
remains 1.399001533608 ns, and the sampled Fig. 5 transition bracket remains
275 < H ≤ 300 Oe. The auxiliary equal-dt per-spin 1e-9 test still misses at
275 Oe (2.93e-9), while the mean discrepancy is only 7.77e-12; its failure
flag remains recorded and is not counted as a pass.

Browser inspection confirmed the live anisotropy table and formula show
2K/μ = 2.2571 T and energy −K(e·s)². Byte comparisons confirmed unchanged
`distrib.m`, exchange coefficient masks, neighbor shifts, and exchange/DMI
field expressions. The spatial-weight issue remains open as requested;
these tests do not establish physical relaxation of the larger array.

The source/result preservation check passed for all 30 recorded paths.
`git diff --check` and whitespace checks for new text files passed. MATLAB
Code Analyzer reports no messages for the new protocol function or updated
initializer; script/harness diagnostics remain and are not treated as a
scientific validation result. No commit or push was made.
