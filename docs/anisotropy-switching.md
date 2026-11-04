# Unequal-anisotropy switching comparison

**The coordination-weighted candidate does not reproduce Figures 2, 4(a), and 5 of PRB 109, 134433 (2024) as a set.** All three spatial sizes fail the low-field destination map and sustained z reversal. Their field sweep lacks the reference transition on the sampled grid. They reproduce the high-field final angle, which is only partial agreement.

Branch: `experiment/unequal-anisotropy-prb2024`, based on `0af802c`. This is an isolated deterministic diagnostic; production MATLAB source and defaults are unchanged. The [undriven six-state result](anisotropy-six-state.md) remains valid, but preserving six minima does not preserve driven switching.

## Figure results

The comparison uses identical current, field, pulse timing, moments and damping. Angles below are modulo 360° and refer to the balanced sublattice mean at 5 ns. The full-strength triangle is the numerical control associated with the paper; the half-strength triangle isolates the candidate's coefficient scaling.

| Model | Low-field destinations matching the reference | +1000 Oe destinations matching the reference | Sampled field-sweep behavior |
|---|---:|---:|---|
| Full-strength three-spin reference | 24/24 | 12/12 | 30° to 90° between 275 and 300 Oe |
| Half-strength three-spin candidate | 4/24 | 12/12 | 30° to 90° between 400 and 425 Oe after local refinement |
| 16-spin candidate | 0/24 | 12/12 | 90° at every sampled field from 175 to 500 Oe |
| 130-spin candidate | 0/24 | 12/12 | 90° at every sampled field from 175 to 500 Oe |
| 336-spin candidate | 0/24 | 12/12 | 90° at every sampled field from 175 to 500 Oe |

**Figure 2:** each spatial patch gives the following result from all six initial states:

| Current sign | Field sign, magnitude 100 Oe | Reference destination | Spatial candidate destination |
|---|---|---:|---:|
| + | + | 30° | 90° |
| + | − | 210° | 270° |
| − | + | 150° | 90° |
| − | − | 330° | 270° |

The candidate loses the reference's current-dependent destination selection. It retains dependence on the field sign.

**Figure 4(a):** the reference balanced z component changes from approximately −0.002997 to +0.002997. The spatial candidates settle near +y with negligible z. The site-weighted mean also differs substantially because the patches have excess b sites.

| Occupied spins | Balanced z at 0.9 ns | Balanced z at 5 ns | Site-weighted z at 0.9 ns | Maximum within-sublattice vector spread during this trajectory |
|---:|---:|---:|---:|---:|
| 16 | −0.00299732 | −2.32×10⁻¹² | +0.0511286 | 0.03442 |
| 130 | −0.00299732 | +6.71×10⁻¹² | +0.0236493 | 0.04106 |
| 336 | −0.00299732 | +1.13×10⁻¹¹ | +0.0124672 | 0.04227 |

Transient balanced-z zero crossings occur at 1.08284, 1.12448 and 1.17131 ns, compared with the reference's 1.39900 ns. These are **not successful reversal times**, because the candidates do not retain the target positive-z final state. The raw output records the full-spin snapshots bracketing each crossing. The half-strength triangle has positive final z but negative final y, ending near 330° rather than 30°; checking z alone would also misclassify that case.

**Figure 5(a):** each spatial patch matches seven of the twelve individual reference endpoints, because the seven higher-field reference cases also end at 90°. This does not reproduce the figure's transition: every candidate point is already at 90°. No conclusion is drawn about an untested lower-field transition. The half-strength triangle's 400–425 Oe bracket is specific to the stated pulse and observation protocol, not a universal critical field.

**Figures 5(b,c):** all six starts reach 90° for either current sign at +1000 Oe on each spatial patch. This is final-state agreement; the transient waveforms remain different.

## Model and observables

The candidate retains the source-derived open graph and DMI signs, with reciprocal $J/2,D/2$ bonds and $K_i=Kz_i/4$. Here $z_i$ counts actual occupied exchange neighbors. Each site retains the same moment $3\mu_B$ and uniform external-field/SOT inputs. The anisotropy field is $2K_i(\mathbf e_i\cdot\mathbf m_i)\mathbf e_i/\mu$.

The diagnostic integrates the existing fixed-length LLGS convention, including the damping-dependent SOT cross term, with coupled unnormalized RK stages and normalization after the combined step. The target remains the mathematical graph model; physical crystallographic mapping and the coordination-dependent edge anisotropy law are **REFERENCE MISSING**.

The saved means are

$$
\mathbf m_{\rm bal}=\frac{\overline{\mathbf m}_a+\overline{\mathbf m}_b+\overline{\mathbf m}_c}{3},
\qquad
\mathbf m_{\rm site}=\frac{1}{N}\sum_i\mathbf m_i.
$$

Neither mean is normalized. Sublattice, edge and interior means and maximum within-sublattice spread are also saved. The three-spin control has no interior sites; those columns explicitly contain NaN.

For these patches, the a,b,c populations are $(5,6,5)$, $(42,46,42)$ and $(110,116,110)$. Writing them as $(n,n+\delta,n)$ gives the exact identity

$$
\mathbf m_{\rm site}=\frac{3n}{N}\mathbf m_{\rm bal}+\frac{\delta}{N}\overline{\mathbf m}_b.
$$

The excess fractions $\delta/N$ are 0.0625, 0.030769 and 0.017857, respectively. These exceed the reference weak moment's magnitude, approximately 0.00346. The Zeeman energy couples to the actual site sum, so plotting the balanced mean cannot remove the excess-b contribution from the dynamics. Together with the degree-dependent internal fields, this provides a plausible explanation for the field-aligned destinations. It is an inference consistent with the result, not a separately isolated causal proof.

## Protocol and numerical evidence

All runs use $J/\mu=100.812$ T before bond weighting, $D=0.833$ meV, $K=0.196$ meV, $\alpha=0.003$, $\gamma=1.760859644\times10^{11}$ rad/(s T), $M_s=441850$ A/m, 30 nm film thickness, spin Hall angle 0.06 and polarization +x. The slight difference between the executable exchange coefficient and the paper's nominal 17.53 meV is retained from the [approved reference alignment](reference-protocols.md).

Ideal 120° seeds evolve for 5 ns; field and current are applied during $[1,4)$ ns. The drive gate is held fixed across RK stages. The 48-case matrix comprises 24 low-field sign/seed combinations, twelve fields from 175 to 500 Oe, and twelve +1000 Oe cases spanning six seeds and both current signs. Figure 4(a) reuses the +current, +100 Oe, 210° first-spin seed. Noise, dipoles, pinning, field-like SOT and STT are disabled. Initial x components are exactly zero in the diagnostic, rather than the source initializer's trigonometric roundoff; full-strength control comparisons verify this has negligible effect here.

**300 full trajectories completed:** 240 primary runs (48 cases × five models), 32 pilot timestep refinements, 24 perturbed-start/refinement runs, and four additional half-triangle transition runs. The 20 short MATLAB/GPU fixtures and short zero-drive regressions are additional checks, not included in 300.

| Check | Maximum discrepancy or result |
|---|---:|
| Source static fields versus independently assembled fields | 4.27×10⁻¹⁴ T |
| Candidate-plus-Zeeman finite-difference energy gradient | 2.13×10⁻⁹ T |
| Explicit versus implicit Gilbert residual divided by gamma | 4.44×10⁻¹⁴ T |
| Twenty actual-source MATLAB/GPU driven replays | 4.05×10⁻¹⁶ per spin component |
| Old/new zero-drive modes | Binary and scalar outputs byte-identical for twelve size/mode combinations |
| 32 pilot comparisons, 5 versus 2.5 fs | 1.13×10⁻⁷ per spin component; gate 10⁻⁵ |
| Twelve perturbed-start timestep comparisons | 2.31×10⁻⁷ per spin component; gate 10⁻⁵ |
| Perturbed versus unperturbed final states | 9.29×10⁻¹⁴ per component |
| Half-triangle transition refinement | 9.61×10⁻¹⁰ per component |
| All full-strength control means versus preserved reference | 7.75×10⁻¹² per component; gate 10⁻¹⁰ |
| Maximum occupied norm error after normalization | 2.23×10⁻¹⁶ |
| Maximum norm error before normalization | 2.00×10⁻¹² |
| Saved means recomputed from sampled full-spin states | 1.10×10⁻¹⁵ per component |

All update counts, pulse-active update counts and saved clocks pass. MATLAB fixtures verify exact-zero vacancies and pulse boundaries across two integration blocks. Long trajectories use the source-checked compiled diagnostic; they are not full-length production GPU runs. The independent read-only verifier recomputed all 240 primary endpoints from binary spins, all twelve perturbed timestep comparisons, and the summary classifications.

No parameter tuning, noise model, physical geometry repair or production promotion was performed. No digitized paper curves were fitted. Qualitative publication comparison and numerical reference-code comparison remain distinct. Three sizes establish the observed trend, not a thermodynamic limit.

## Output location and figures

At the user's request, new run data, MATLAB logs and generated figures are stored **outside OneDrive**:

```text
D:\Git_control\spin_dynamics\simulation_results\Mn3Sn\unequal-anisotropy-switching-20261103-02
/mnt/d/Git_control/spin_dynamics/simulation_results/Mn3Sn/unequal-anisotropy-switching-20261103-02
```

Each run's `relocation.json` records its earlier move out of OneDrive. The workspace-level `D:\Git_control\spin_dynamics\relocation-20261103.json` records the earlier project/result consolidation and SHA-256 hashes. The [workspace README](../../README.md) maps those recorded paths to the current root after the later move into `D:\Git_control\spin_dynamics`. Historical reference datasets remain at their recorded paths. The preliminary `-01` setup-only folder was also relocated. Original command records retain the paths actually used when the computations ran.

| External file | Content |
|---|---|
| `figures/anisotropy-switching-fig2.png` | Polarity and six-start angle comparison for all three spatial sizes |
| `figures/anisotropy-switching-fig4a.png` | Balanced and site-weighted component traces |
| `figures/anisotropy-switching-fig5.png` | Final-angle field sweep and both high-field polarities |
| `figures/anisotropy-switching-triangle.png` | Full-strength and half-strength triangle transients |
| `figure-results.json`, `figure-results.csv`, `figure-summary.json` | Per-case values, numerical differences and group outcomes |
| `pilot-checks.json`, `triangle-transition.json` | Timestep and sampled-transition evidence |
| `gpu-gates.json`, `driven-gates.json`, `independent-*.json` | Source/GPU and independent verification |
| `*-observables.csv` | Dense means every 0.25 ps, sublattice/edge/interior data and spread |
| `*.bin`, `*-snapshot-times.csv`, `*-final.bin` | Full-spin snapshots, exact times and final states |
| `*-run.json`, `compiler.json`, `provenance.json`, `environment.json` | Commands, configuration, compiler/source hashes and environment |

Full-spin snapshots are sampled every 5 ps, with denser snapshots near pulse edges; analysis records the snapshots bracketing each measured reversal. Plots unwrap angles on each continuous valid interval, omit the first 0.05 ns from angle traces and mask in-plane mean magnitudes below $10^{-10}$. Component plots retain startup. No smoothing or unit normalization is applied, and axes retain the full angular trajectories.

All four PNG exports were visually inspected. `plot-checks.json` confirms 163 plotted series and 20,001 samples per trajectory. Unwrapped endpoints such as −270° represent the same orientation as 90°; endpoint classifications use circular angle differences.

## Reproduction

Run from the repository root. The prepared source matrices and reference results are the existing datasets documented in [the spatial relaxation report](spatial-relaxation.md), [six-state report](anisotropy-six-state.md), and [reference-protocol report](reference-protocols.md). Use a fresh external output directory. The new runner rejects output paths inside OneDrive.

```bash
python3 verification/spatial_relaxation/driven_figures.py /mnt/d/Git_control/spin_dynamics/simulation_results/Mn3Sn/unequal-anisotropy-new prepare
matlab -wait -batch "addpath('verification/spatial_relaxation'); prepare_driven('.sisyphus/spatial-relaxation-20261102-01','D:/Git_control/spin_dynamics/simulation_results/Mn3Sn/unequal-anisotropy-new')"
python3 verification/spatial_relaxation/driven_figures.py /mnt/d/Git_control/spin_dynamics/simulation_results/Mn3Sn/unequal-anisotropy-new check
python3 verification/spatial_relaxation/driven_figures.py /mnt/d/Git_control/spin_dynamics/simulation_results/Mn3Sn/unequal-anisotropy-new pilot
python3 verification/spatial_relaxation/driven_figures.py /mnt/d/Git_control/spin_dynamics/simulation_results/Mn3Sn/unequal-anisotropy-new pilot-check
python3 verification/spatial_relaxation/driven_figures.py /mnt/d/Git_control/spin_dynamics/simulation_results/Mn3Sn/unequal-anisotropy-new sweep
python3 verification/spatial_relaxation/driven_figures.py /mnt/d/Git_control/spin_dynamics/simulation_results/Mn3Sn/unequal-anisotropy-new transition
python3 verification/spatial_relaxation/driven_figures.py /mnt/d/Git_control/spin_dynamics/simulation_results/Mn3Sn/unequal-anisotropy-new analyze
matlab -wait -batch "addpath('verification/spatial_relaxation'); out='D:/Git_control/spin_dynamics/simulation_results/Mn3Sn/unequal-anisotropy-new'; mkdir(fullfile(out,'figures')); plot_driven(out,fullfile(out,'figures'))"
python3 verification/spatial_relaxation/driven_figures.py /mnt/d/Git_control/spin_dynamics/simulation_results/Mn3Sn/unequal-anisotropy-new preserve
```

Environment: MATLAB R2024b Update 1, NVIDIA GeForce RTX 3060 Laptop GPU, GNU C++ and standard-library Python; exact versions and commands are in the external artifacts. No dependencies were installed. The MATLAB metadata-aggregation issue was fixed and all fixtures rerun; a plotting legend-format error was corrected before successful export. Neither issue changed the physical simulation.

## Sources and interpretation

Zhengde Xu, Xue Zhang, Yixiao Qiao, Gengchiau Liang, Shuyuan Shi, Zhifeng Zhu, *Deterministic spin-orbit torque switching including the interplay between spin polarization and kagome plane in Mn3Sn*, Physical Review B **109**, 134433 (2024), DOI: [10.1103/PhysRevB.109.134433](https://doi.org/10.1103/PhysRevB.109.134433). Read the user-supplied local reviewed summary and fulltext; visually checked PDF pp. 134433-2, -3 and -5 and checked Table I on p. 134433-4. The paper supplies the target model and figure protocols, not the candidate edge law.

The equation conventions were checked against `presentations/equations/spintronic_methods.tex`, `sections/hamiltonian.tex`, `sections/llg_equation.tex` and `sections/atomistic_boundaries.tex`. The reviewed `presentation-knowledge/notes/ee237/week15a-afm/slides.md`, slides 25–27, supplies qualitative context. ASTRO AGENTS/README supplied workflow context only.

Keep this branch as evidence about the candidate's limitations. Restoring equilibrium degeneracy is insufficient to justify replacing the production spatial interactions with this coefficient rule. A subsequent model investigation would need to isolate finite-population and boundary effects before claiming paper-equivalent spatial switching.
