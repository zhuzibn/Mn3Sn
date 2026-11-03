# A02 — Nonreciprocal exchange: historical relaxation evidence

Evidence summary prepared on 2026-10-29 from [discussion_zhengde.md][D] and Zhengde Xu's weekly-report presentations. This note records the historical reason for the middle/interior exchange correction and the limits of the available results. The current implementation audit remains in [plan.md, A02](plan.md#a02-follow-up-page-19-structure-comparison-2026-10-21).

**The strongest documented result is that halving the middle/interior exchange produced the expected relaxation in both 3×3 and 15×15 samples. However, a 2:1 edge-to-interior ratio alone did not guarantee success: doubling edge exchange was subsequently reported as unsuccessful, and relaxation depended strongly on timestep.**

## Evidence conventions

- **Reported correct** means the slide or discussion calls the result correct. It does not mean independently validated physics or a demonstrated converged trajectory.
- **Reported unsuccessful** may mean an unexpected final state, slow relaxation, or continued oscillation. These outcomes are distinguished below.
- $A_0$ denotes the baseline exchange coefficient for the particular historical run. The relevant slides do not consistently state its absolute value or establish identical baseline parameters between reports. Do not assign the current code's approximately 17.53 meV to every historical run.
- Sizes such as 3×3 denote simulation-array dimensions, not the number of Mn atoms. The discussion identifies the 3×3 example as five Mn atoms forming two triangles sharing one atom. April 23, slide 2 identifies the 2×2 example as three atoms.
- `Aex_mid/2`, `Aex_side*2`, and similar expressions preserve the slide terminology. They describe the reported site/region scaling; they do not establish a reciprocal bond-energy convention or the exact historical masks.
- Slide numbers below are positions in the presentation, starting at 1. The historical presentation summary did not rerun those spatial simulations. The separate Figure 5 code verification below includes new numerical checks.

## System size, exchange setting, and relaxation outcome

| System | Exchange setting | Recorded relaxation outcome | Evidence |
|---|---|---|---|
| **2×2: three Mn atoms** | No special middle-site correction indicated; absolute $A_{\rm ex}$ unspecified | **Successful:** relaxation completed within 10 ns at $\Delta t=2\times10^{-15}$ s. This supports the recollection that the isolated three-atom system worked, but the slide does not specify a complete parameter set. | [April 23, slide 2][P5] |
| **3×3: five Mn atoms, two triangles sharing one atom** | Same exchange strength at middle and edge sites | **Incorrect according to Zhengde's recollection:** the shared atom needed its exchange halved. The discussion supplies this comparison, rather than a fully specified pair of slide results. | [Discussion][D] |
| **9×9** | Not specified on the slides; tested original code and version 3.1 | **Incorrect relative to the six-state expectation:** only two stable states were reported across the tested initial angles. The slides do not independently establish the exchange mask used. | [March 28, slides 1–5][P1] |
| **3×3: five Mn atoms** | **Middle atom $A_0/2$**, others unchanged in the described modification | **Explicitly reported correct:** initial 40° → final 30°. The figure shows five spins. | [April 2, slide 4][P2] |
| **15×15** | **Interior atoms $A_0/2$**, others unchanged in the described modification | **Explicitly reported correct**, with initial/final spatial spin plots. | [April 2, slide 5][P2] |
| **3×3** | **Edge atoms $2A_0$**, middle unchanged in the described modification | **Reported unsuccessful:** “Aex_side*2不行.” Tests at $\Delta t=2\times10^{-15}$ and $5\times10^{-15}$ s show very slow relaxation. This does not establish that the eventual equilibrium would be wrong. | [April 23, slide 2][P5]; corresponding comparison in [April 12, slide 1][P4] |
| **15×15, with DMI added** | **Edge atoms $2A_0$** | **Unresolved/unsuccessful relaxation:** continued oscillation, “仍然振荡.” | [April 9, slide 5][P3] |
| **5×5 and 15×15** | **Interior exchange and DMI both halved** | **Strong timestep dependence:** slow approach at $5\times10^{-16}$ s; 15×15 tests at $2\times10^{-16}$ and $10^{-16}$ s approach the apparent target faster. No quantitative convergence criterion is supplied. | [April 23, slides 4–6][P5] |
| **2×2, 3×3 and 15×15** | Exchange/DMI masks **not restated** | **Later encouraging comparison:** plotted $m_y(t)$ curves at the labelled $5\times10^{-16}$ s nearly overlap across these sizes and reach the same plateau. This observation cannot be assigned confidently to a particular exchange prescription. | [April 28, slides 3–5][P6] |

## Additional recorded observations

### Timestep and damping

- **3×3 with middle exchange halved:** April 23, slide 3 compares $\Delta t=10^{-15}$, $5\times10^{-16}$, and $2\times10^{-16}$ s and explicitly states “Tstep影响很大.” It also includes a case labelled `Aex_mid/2 Aex_constant*2`; this is a separate parameter change, not a repetition of the same absolute exchange setting. [Source][P5]
- **Larger samples:** April 23, slides 4–6 state “样本越大，就需要减小tstep.” Some comparisons also change DMI scaling or multiply the overall exchange and DMI constants by two, so they do not isolate system size alone. [Source][P5]
- **Later timestep scans:** April 28, slides 3–5 show timestep-dependent relaxation rates for 2×2, 3×3, and 15×15, alongside the nearly overlapping size-comparison curves. The plots are evidence of the displayed trajectories, not a reported convergence tolerance or validated integration order. [Source][P6]
- **Damping:** April 2, slide 6 compares $\alpha=0.003$, $0.03$, and $0.3$ in a 3×3 system and reports a strong effect on relaxation. April 23, slides 7–8 separately report slower relaxation at larger damping in the three-sublattice macrospin model. These are distinct model comparisons. [April 2][P2], [April 23][P5]

### Geometry and periodic boundaries

- **Adding complete triangles at boundaries did not fix the result.** April 2, slides 2–3 show an initial angle of 40° and label the final result “Final not 30,” both when two triangles were added below and when one was added at each end. [Source][P2]
- **An inverted-triangle sample was labelled correct, but the criterion is ambiguous.** April 2, slide 8 says the geometry ensures equal populations of the three atom types and calls the result correct. Its annotations nevertheless read **initial 40° → final 40°**, so it does not demonstrate the same 40° → 30° relaxation as slide 4. [Source][P2]
- **Periodic boundaries were not an established cure.** April 2, slide 1 describes periodicity as giving correct results, while slide 7 reports slower relaxation. April 9, slide 4 records continued oscillation after adding DMI in a periodic case and suggests that the DMI implementation might be wrong. The exact size and complete settings are not supplied on those slides. [April 2][P2], [April 9][P3]

## Interpretation and limits

**Halving the middle and doubling the edges are different parameter sets.** Their respective `(middle, edge)` coefficients are

$$
(A_0/2,\ A_0)
\qquad\text{and}\qquad
(A_0,\ 2A_0).
$$

They share the same edge-to-middle ratio, but the latter doubles the absolute exchange throughout relative to the former, assuming the same $A_0$ and site partition. Thus the successful and unsuccessful reports cannot be reduced to the ratio alone. The reports also do not establish identical DMI, damping, timestep, initial conditions, or code versions across all comparisons.

The discussion supports the historical explanation that equal middle/edge exchange failed in expanded samples and that the shared-site correction was introduced to recover the expected states. The slides directly support particular successful middle-halved examples, while also recording failures and numerical sensitivities. They do not establish that this correction is physically necessary.

The records do **not** establish:

1. A controlled, sufficiently converged comparison of uniform versus site-scaled exchange with every other setting identical.
2. Recovery of **all six** expected stable states for the successful examples; the 40° → 30° example establishes only that particular reported relaxation.
3. A common absolute $A_{\rm ex}$, complete exchange/DMI masks, or exact archived code version for every displayed run.
4. That a slow trajectory has the wrong eventual equilibrium, or that continued oscillation uniquely identifies an exchange error.
5. Reciprocity or energy-gradient consistency of the historical interaction coefficients.

**Relation to A02:** the existing [A02 audit](plan.md#a02-follow-up-page-19-structure-comparison-2026-10-21) separately records unequal reverse exchange weights in the current code. The historical reports help explain why site-dependent scaling was introduced; they do not resolve that audit finding. No source-code change or physical-model correction is made or authorized by this note.

## Published PRB Figure 5 code verification — 2026-10-29

**Result:** the published `code/fig5` model reproduces the sampled Figure 5 final-state outcomes and relaxes to all six expected zero-drive states. Its active exchange is reciprocal and does not contain the middle/edge weighting implicated in A02. However, the archived integrator has an incorrect fourth RK stage, and the supplied angle-processing script loses the sign of the octupole angle. Correct equilibrium destinations do not certify the transient integration accuracy.

### Scope, model, and paper protocol

The user identified this source folder:

```text
C:\Users\zzf-m\OneDrive\projects\published\2024.04-SOT Mn3Sn-zhengde-PRB\code\fig5
```

This is a **three-sublattice noncollinear macrospin model with fixed-length LLGS dynamics**, not a spatial finite-size lattice. The script stores six vectors under legacy Fe/Gd names, but `Hex_1=0` decouples the two triplets. The reported observable uses only the first three vectors. The second triplet has no active anisotropy or DMI and zero SOT polarization, so it does not affect the audited first-triplet trajectories. No boundary/interior coefficient or shared-site correction appears in the active triplet.

The source PDF of [Xu et al., Physical Review B 109, 134433 (2024)][PRB] was read directly, including Figure 5 on printed page 134433-5 and the methodology/protocol on pages 134433-2–3. Figure 5(a) uses $J_c=+5\times10^{10}\,\mathrm{A/m^2}$ and fields of 175, 200, 225, 250, 275, 300, 325, 350, 375, 400, 450, and 500 Oe. Figure 5(b,c) uses 1000 Oe, the six initial states, and opposite current signs. The drive is applied from 1 to 4 ns within a 5 ns run, with 1 ns relaxation before and after.

The supplied `main.m` selects **one 250 Oe case**, with $J_c=+5\times10^{10}\,\mathrm{A/m^2}$, first-spin initial angle 210°, $\Delta t=5$ fs, $\alpha=0.003$, $\theta_{\rm SH}=0.06$, thickness 30 nm, and $M_s=441850\,\mathrm{A/m}$. The folder does not contain a complete Figure 5 sweep/plot driver. With the source constants, `Hex_2=100.812` T corresponds to **17.50806 meV**, about 0.125% below the nominal paper value 17.53 meV; this audit retains the implemented value. The anisotropy and DMI fields are $2K/\mu=2.2571490$ T and $D/\mu=4.7964417$ T, with $K=0.196$ meV, $D=0.833$ meV, and $\mu=3\mu_B$.

### Executed verification and reproduction results

The original published folder was read-only. Its 14 files, including `final.mat` and `result.txt`, were SHA-256 hashed before and after the audit and remained unchanged. A temporary copy under `/tmp/fig5-audit-rd2emqsx` held the checks and outputs.

For efficient sweeps, a C++ replay transcribed the active deterministic three-spin equations, including all damping-dependent torque cross terms and the legacy fourth-stage error. Its reduction to three spins uses the implemented zero inter-triplet coupling. It was checked against both a fresh short execution of the original MATLAB routines and the existing full-length saved trajectory; it was not assumed equivalent solely from inspection.

| Check | Executed result |
|---|---|
| Environment | MATLAB R2024b Update 1, `24.2.0.2740171`, Windows `PCWIN64`; C++ compiled with `g++ -O3 -std=c++17` under WSL. CPU checks; no GPU required. |
| Original MATLAB probe | 2000 allocated states / 1999 updates over a nominal 10 ps, with original field and integration routines. Overrides: drive starts at 2 ps and lasts 6 ps; save every state; suppress text output by making its interval exceed the run. Maximum spin-component difference from the matching replay: **$4.72\times10^{-14}$**. |
| Full archived default trajectory | Compared all **20000 saved frames** of `final.mat` to the 5 ns legacy replay. Maximum first-triplet component difference: **$9.35\times10^{-10}$**. Archived maximum $|\|\mathbf m_i\|^2-1|$ was $1.56\times10^{-15}$. |
| Figure 5(a), 175–275 Oe | All five sampled fields ended at **30°** after drive removal. At 250 Oe, the angle immediately before removal was **30.8404323°**, and the final angle was **29.999999993°**. |
| Figure 5(a), 300–500 Oe | All seven sampled fields ended at **90°**. This brackets the sampled transition between 275 and 300 Oe, consistent with the paper's statement that the change occurs above 275 Oe. It does not locate an exact critical field. |
| Figure 5(b,c), 1000 Oe | All six relaxed initial states ended at **90° for either current sign**. Immediately before removal, the angles were approximately **83.831109°** for positive current and **96.168891°** for negative current; both relaxed to 90° afterward. |
| Diagnostic full-step fourth RK stage | Repeated every sampled Figure 5 case with the fourth-stage state evaluated at the full step. **All final-state classifications were unchanged.** This variant existed only in the audit replay; published code was not repaired. |
| Sweep size | 80 recorded cases: 24 Figure 5(a) runs, 24 high-field/current-sign/initial-state runs, 24 perturbed zero-drive relaxation runs, and 8 additional timestep runs. Two default comparisons and the short replay probe were also executed. |

The source's first-spin initialization angle is **not generally the relaxed octupole angle**. For nominal first-spin angles 30°, 90°, 150°, 210°, 270°, and 330°, the measured pre-pulse mean-spin directions are respectively 30°, −30°, −90°, −150°, 150°, and 90°. The six-state sweep was classified by these measured directions. At the exact initial 120° arrangement the mean spin is zero up to rounding, so its angle is initially undefined.

The Figure 5 comparison establishes the sampled final-state map and visually consistent trajectories, not pixel-level agreement with digitized published data. The full 5 ns original MATLAB script was not rerun; the full-length check uses its existing archive plus the independently executed short MATLAB probe.

### Does it relax correctly?

**Yes for the tested zero-drive equilibrium structure.** Twelve initial conditions used first-spin angles of $30°,90°,150°,210°,270°,330°$, each perturbed by $\pm10°$, with the other two spins initially offset by 120° and 240°. Current and external field were zero. Both legacy and full-step-stage methods reached the six expected mean-spin directions:

$$
\phi_{\rm oct}=30°,\ 90°,\ 150°,\ 210°,\ 270°,\ 330°.
$$

The mean-spin magnitude was **0.00346100981244** in every final state. The maximum legacy final angular residual was $7.35\times10^{-11}\,\mathrm{rad/ns}$; final-angle deviations from the six directions were below $10^{-9}$ degrees. These are measured numerical values, not a claim of corresponding physical precision.

Independent checks used the **original `field_eta.m`**, rather than only the replay:

- The reconstructed 9×9 conservative-field matrix had maximum antisymmetric entry **$1.11\times10^{-16}$ T** for unit-spin inputs. The three exchange pairs have equal forward/reverse coefficient `Hex_2`; the A02 nonreciprocity is absent here.
- A central finite-difference energy-gradient check, with component perturbation $10^{-5}$, agreed with the original field to **$6.63\times10^{-10}$ T**.
- All six final states had the same energy, $E/\mu=-165.378240880411$ T. Their constrained tangent-space energy Hessians were positive definite; the smallest eigenvalue was approximately **$2.65906\times10^{-4}$ T**. Thus these states are local minima of the implemented energy, not merely stationary points.

The energy used for these checks was derived from the active source fields, with each pair counted once:

$$
\frac{E}{\mu}=J_B\sum_{(i,j)=(1,2),(2,3),(3,1)}\mathbf m_i\cdot\mathbf m_j
+D_B\sum_{(i,j)=(1,2),(2,3),(3,1)}\hat{\mathbf x}\cdot(\mathbf m_i\times\mathbf m_j)
-\frac{K_B}{2}\sum_i(\mathbf e_i\cdot\mathbf m_i)^2.
$$

Here $J_B=100.812$ T, $D_B=4.79644166487$ T, $K_B=2.25714901876$ T, and $\mathbf e_i=(0,\sin\theta_i,\cos\theta_i)$ with $\theta_i=(30°,270°,150°)$. This establishes consistency of this code's conservative fields with its derived energy. It does not independently resolve every pair-counting convention in the paper's compact Hamiltonian notation or validate a spatial lattice extension.

### Numerical and output defects that remain

1. **Fourth RK stage uses a half step.** In `rk4_4llg.m`, lines 158–164 construct the fourth-stage states as `m + k3*ts1/2`, rather than `m + k3*ts1`. Fields are recomputed at every stage, unlike the frozen-field spatial code, but the stage location is wrong. The routine is not classical fourth-order RK4.

   | Method and timestep | Maximum spin-component difference from full-step-stage RK4 at 1.25 fs, default 250 Oe trajectory |
   |---|---:|
   | Legacy, 5 fs | $5.15\times10^{-3}$ |
   | Legacy, 2.5 fs | $2.23\times10^{-3}$ |
   | Legacy, 1.25 fs | $1.06\times10^{-3}$ |
   | Full-step fourth stage, 5 fs | $1.63\times10^{-8}$ |
   | Full-step fourth stage, 2.5 fs | $9.60\times10^{-10}$ |

   Comparisons use matched physical sample times and exclude the unmatched final update. The default legacy driven-angle discrepancy reaches **0.241°**. Final destinations are robust in the tested cases, while the transient error decreases much more slowly than fourth order.

2. **Energy relaxation is not strictly monotonic with the legacy integrator.** For the zero-drive 40° initial case, the largest per-update increase in $E/\mu$ was **$8.63\times10^{-7}$ T** at 5 fs, versus approximately $1.42\times10^{-13}$ T for the full-step-stage comparison. The maximum pre-normalization norm error in the perturbed relaxation tests was about **$7.97\times10^{-8}$** for legacy integration and $3.87\times10^{-13}$ for the full-step-stage comparison. Post-step unit normalization alone hides this numerical error.

3. **`theta_mtot.m` loses the angle sign.** Lines 8–10 replace the mean $y$ component by $\sqrt{\bar m_x^2+\bar m_y^2}$ before `atan2`. This folds angles into 0–180°: the archived −150° initial relaxed state is displayed as +150°. The audit uses the signed in-plane angle `atan2(mean_y,mean_z)`, matching the signed-angle presentation in Figure 5.

4. **Saved time labels are slightly inconsistent with update counts.** `totstep=runtime/tstep`, `linspace(0,runtime,totstep)`, and `while ct1<totstep` yield 999999 updates, not 1000000. The actual last integrated time is 4.999995 ns. The last decimated archived state is at actual 4.999750 ns but is labelled 4.999754999755 ns, a drift of about 5 fs. This is small for Figure 5 but remains an output-contract defect.

5. **Default output is not a complete Figure 5 reproduction.** `main.m` runs only the 250 Oe case and overwrites `final.mat`; `result.txt` is appended. Sweeps and signed-angle postprocessing were supplied by the isolated audit harness. No production run was launched in the published folder.

### Artifacts, commands, and evidence limits

- [Verification plots, PNG](docs/fig5-verification/verification.png) and [PDF](docs/fig5-verification/verification.pdf): Figure 5 sweeps, six-state relaxation, energy comparison, and timestep errors.
- [80-case numerical summary](docs/fig5-verification/summary.csv).
- [Original-source hashes](docs/fig5-verification/source-hashes.json), [original-field matrix](docs/fig5-verification/source-field-matrix.csv), and [static-check log](docs/fig5-verification/static-checks.txt).
- [Audit harness archive](docs/fig5-verification/audit-harness.zip): source snapshots, replay, run manifest, original-MATLAB probe, plotting/static checks, and rerun instructions. Extract into an empty directory; do not execute `main.m` in the published folder.

Representative executed commands from the temporary audit directory:

```bash
g++ -O3 -std=c++17 replay.cpp -o replay
./replay legacy 5e-15 5 250 5e10 210 1 3 default_legacy.csv
./replay rk4 5e-15 5 250 5e10 210 1 3 default_rk4.csv
matlab -batch "run('//wsl.localhost/Ubuntu-22.04/tmp/fig5-audit-rd2emqsx/audit_probe.m')"
matlab -batch "run('//wsl.localhost/Ubuntu-22.04/tmp/fig5-audit-rd2emqsx/audit_static_and_plot.m')"
```

The common equation reference consulted was `presentations/equations/spintronic_methods.tex` and `sections/llg_equation.tex`, especially `eq:llg-gilbert` and `eq:ll-with-flt`; these supplied the tesla convention and Gilbert-to-explicit torque cross-term check. The paper's indexed summary supplied navigation, and the source PDF supplied the Figure 5 protocol. No network access, dependency installation, source repair, or modification of the published archive was performed. Thermal noise, interlayer coupling, spatial extensions, other paper figures, and a fully converged critical-field boundary were outside this check.

**Implication for A02:** the working reciprocal three-sublattice relaxation provides a useful small-system comparison. It supports Zhengde's recollection that the isolated three-spin model worked. It does not establish whether uniform reciprocal exchange in a correctly implemented extended finite lattice would fail or succeed; that requires a separate spatial-code audit and controlled comparison.

## Source chronology and provenance

The search covered the slide text of 95 dated presentations from 2024–2026, with visual inspection of the key figures. The source folder is:

```text
C:\Users\zzf-m\OneDrive\COSMAL\本科生研究生\研究生\archived\2021\徐正德\weekly report
```

Zhengde recalled the first discussion as **2024.2.3**. No presentation with that date was found in the folder. The earliest explicit Mn3Sn troubleshooting mention found in the searched material is [March 19, slide 1][P0]; the explicit middle-exchange-halving result appears in [April 2, slides 4–5][P2]. These are dates of available records, not proof of when the first simulation or conversation occurred.

| Reference | Source | Relevant content |
|---|---|---|
| D | [discussion_zhengde.md][D] | Retrospective account of uniform exchange, the three-/five-atom distinction, and shared-site halving |
| P0 | [2024.3.19.pptx][P0], slide 1 | Mn3Sn atomistic-model troubleshooting underway |
| P1 | [2024.3.28.pptx][P1], slides 1–5 | Two stable states in 9×9; original code and version 3.1 |
| P2 | [2024.4.2.pptx][P2], slides 1–8 | Geometry tests, middle-halved 3×3 and 15×15, damping and periodicity |
| P3 | [2024.4.9.pptx][P3], slides 1–5 | Middle-halved versus edge-doubled exchange; slow relaxation and DMI oscillations |
| P4 | [2024.4.12.pptx][P4], slide 1 | 3×3 edge-doubled and 2×2 timestep comparison |
| P5 | [2024.4.23.pptx][P5], slides 2–8 | Three-atom success, edge-doubled failure, timestep/size/DMI/damping comparisons |
| P6 | [2024.4.28.pptx][P6], slides 3–5 | Later timestep and system-size comparisons; exchange masks not restated |

External presentation and paper links use absolute Windows paths because this checkout is on D:. In WSL, replace `C:/Users/zzf-m/OneDrive/` with `/mnt/c/Users/zzf-m/OneDrive/`. The original discussion and presentations were not modified. The historical presentation review involved no new spatial simulations; the separate Figure 5 verification above documents its own numerical checks.

PRB source: Zhengde Xu, Xue Zhang, Yixiao Qiao, Gengchiau Liang, Shuyuan Shi, Zhifeng Zhu, Deterministic spin-orbit torque switching including the interplay between spin polarization and kagome plane in Mn3Sn, Physical Review B **109**, 134433 (2024). DOI: 10.1103/PhysRevB.109.134433. [Local source PDF][PRB]; methodology and protocol pp. 134433-2–3, Figure 5 and discussion p. 134433-5.

[D]: discussion_zhengde.md
[P0]: <C:/Users/zzf-m/OneDrive/COSMAL/本科生研究生/研究生/archived/2021/徐正德/weekly report/2024.3.19.pptx>
[P1]: <C:/Users/zzf-m/OneDrive/COSMAL/本科生研究生/研究生/archived/2021/徐正德/weekly report/2024.3.28.pptx>
[P2]: <C:/Users/zzf-m/OneDrive/COSMAL/本科生研究生/研究生/archived/2021/徐正德/weekly report/2024.4.2.pptx>
[P3]: <C:/Users/zzf-m/OneDrive/COSMAL/本科生研究生/研究生/archived/2021/徐正德/weekly report/2024.4.9.pptx>
[P4]: <C:/Users/zzf-m/OneDrive/COSMAL/本科生研究生/研究生/archived/2021/徐正德/weekly report/2024.4.12.pptx>
[P5]: <C:/Users/zzf-m/OneDrive/COSMAL/本科生研究生/研究生/archived/2021/徐正德/weekly report/2024.4.23.pptx>
[P6]: <C:/Users/zzf-m/OneDrive/COSMAL/本科生研究生/研究生/archived/2021/徐正德/weekly report/2024.4.28.pptx>
[PRB]: <C:/Users/zzf-m/OneDrive/papers/AFM/2024-Deterministic spin-orbit torque switching including the interplay between spin-PRB-zhengde Xu.pdf>
