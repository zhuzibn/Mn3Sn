# Deterministic benchmark validation — 2026-11-19

This report records the initial benchmark implementation. The later
[configuration/result migration](configuration-results.md) replaces `main.m`
text extraction with the callable production API while preserving these accepted
artifacts, comparator and thresholds. Use [benchmark/README.md](../benchmark/README.md)
for current commands and missing-reference behavior.

The new benchmark freezes the existing periodic Mn3Sn software behavior and
compares future full saved trajectories and GPU timings. Its runner, comparison
and failure tests are under `benchmark/`; usage and the result contract are in
[benchmark/README.md](../benchmark/README.md). The user approved the scope in
`plan.md` before implementation. No production `.m` file was changed.

## Adopted scope

Class 3 simplified spatial atomistic code graph with fixed-length LLGS. The
benchmark preserves periodic half-strength exchange/DMI, uniform onsite
anisotropy, local easy axes, Mn parameters, explicit damping-dependent torque
terms and coupled RK4 stages. It adds software orchestration, not a material model
or crystallographic construction. Physical lattice and exchange-shell mapping
remain **REFERENCE MISSING**; no physical construction is implemented here.

Quick/current modes use 4x4/8x8 slots (12/48 spins). Each runs production relaxation
for 10 ps at 0.5 fs, plus a nonuniform 60-update, 5 fs driven stress test with a
pulse active on updates 13 through 43. The latter tests nonuniform coupled dynamics
and pulse edges inside 20-update chunks, not physical long-time switching.
Vacancies stay zero, occupied spins have unit length, and ordinary/class means
exclude vacancies. The retained relaxation gate has 20,000 active updates with
zero drive amplitude. No stochastic path, dipole, pinning or open-edge case is used.

In this initial implementation, initialization came from guarded extraction of `main.m` before its
visualization/integration/save tail. Production initialization and integration
execute in a fresh external source snapshot. Checks fail if expected anchors or
supported option assignments change. The runner has no duplicated material-default
block. The nonuniform initial state is generated once by a deterministic algebraic
formula and saved as actual arrays; candidates load those frozen arrays only after
checking their current configuration, mask, axes, classes and generated seed.

## Environment and provenance

- Production source revision: `8bacbeb47a2c5bf285945c08261a9cadf1485583` on `main`.
  Benchmark files and documentation were uncommitted during validation; actual
  source hashes and snapshots identify the tested code.
- Windows MATLAB `24.2.0.2740171 (R2024b) Update 1`, Parallel Computing Toolbox 24.2.
- NVIDIA GeForce RTX 3060 Laptop GPU, compute capability 8.6, 6 GiB VRAM.
- MATLAB-reported driver `610.60`, toolkit `12.2`, Windows 11 Home Chinese edition.
- Commands used the existing WSL `matlab -wait -batch` launcher. No dependencies
  were installed. Toolbox inventories are retained in each `metadata.json`.
- Results root: `D:/Git_control/spin_dynamics/simulation_results/Mn3Sn/benchmark/20261119-v1`
  (WSL: `/mnt/d/Git_control/spin_dynamics/simulation_results/Mn3Sn/benchmark/20261119-v1`).

All eleven production `.m` files were byte-for-byte equal to their baseline
snapshots and metadata SHA-256 hashes, independently checked by a read-only
reviewer. Root `finalnew.mat` and `atomtypee.mat` were absent before the work and
remained absent; historical data were not modified.

## Commands and staged checks

From the repository root, initialization was checked without integration:

```matlab
addpath('benchmark'); root=pwd;
for mode={'quick','current'}
    for name={'relaxation','driven_nonuniform'}
        r=benchmark_case(mode{1},name{1},root,false);
        fprintf('INIT PASS %s %s %d spins\n',mode{1},name{1},nnz(r.atomtype==1));
    end
end
```

All four combinations passed mask, dimension, class-count and initial-angle
checks. The initial baselines were then generated (no manifest existed yet):

```matlab
addpath('benchmark');
out='D:/Git_control/spin_dynamics/simulation_results/Mn3Sn/benchmark/20261119-v1';
run_baseline_benchmark('quick',fullfile(out,'quick-baseline'));
run_baseline_benchmark('current',fullfile(out,'current-baseline-retry'));
```

Both modes/cases passed all saved-state, vacancy, norm, update-count and clock
invariants. Three warm repeats per case agreed exactly with the first trajectory.

The comparator and safe-path checks were run after those gates passed:

```matlab
check_benchmark(fullfile(out,'quick-baseline'),fullfile(out,'comparator-tests-final'));
assert(compare_to_baseline(fullfile(out,'current-baseline-retry'), ...
    fullfile(out,'current-baseline-retry')));
```

The final check passed 21 candidate mutations: intermediate rotation, NaN, Inf,
vacancy contamination, shape, mask, axes, classes, time, saved steps, pulse count,
configuration, initial state, result schema, stored summary, missing field,
missing case, missing case inventory, incomplete status, mode and suite schema.
An invalid baseline was also rejected. The intermediate rotation preserved the
spin norm and final state; the full-trajectory gate detected a maximum component
change of approximately `7.05e-5`, above the `1e-10` tolerance.

The tests also verified that a source change can pass, invalid warm timings and a
complex-valued median suppress speedup reporting, existing output is rejected
without altering baseline metadata, and in-repository output is rejected before
creation. MATLAB `checkcode` found no syntax errors. Its remaining diagnostics
concern variables created by the guarded script `eval` and minor style suggestions;
these are recorded in `comparator-tests-final/checks.json`.

Fresh candidate procedure (each runner invocation includes one first run and
three warm repeats for each case):

```matlab
run_baseline_benchmark('quick',fullfile(out,'quick-candidate'),fullfile(out,'quick-baseline'));
[passed,report]=compare_to_baseline(fullfile(out,'quick-candidate'),fullfile(out,'quick-baseline'));
f=fopen(fullfile(out,'quick-comparison.json'),'w');
fprintf(f,'%s',jsonencode(report,PrettyPrint=true));fclose(f);assert(passed);
run_baseline_benchmark('current',fullfile(out,'current-candidate'),fullfile(out,'current-baseline-retry'));
[passed,report]=compare_to_baseline(fullfile(out,'current-candidate'),fullfile(out,'current-baseline-retry'));
f=fopen(fullfile(out,'current-comparison.json'),'w');
fprintf(f,'%s',jsonencode(report,PrettyPrint=true));fclose(f);assert(passed);
```

Both fresh candidates passed: every saved spin component, occupied mean and
class mean matched its numerical reference exactly (maximum difference zero).
The reports are `quick-comparison.json` and `current-comparison.json` under the
external results root. All timed repeats also matched exactly.

The accepted manifest selects `quick-baseline/` and `current-candidate/`. The
latter retains its original run-directory name but is deliberately promoted as
the current baseline after verification, because the earlier two current runs
had overlapping timing windows. Its frozen inputs are included in its result
files; future runs do not need to follow the provenance chain to an older run.
The final quick/current candidates ran sequentially in one MATLAB batch; Windows
process inspection confirmed only that batch's launcher and worker were active.
Earlier reference directories and comparison reports are preserved.

After registering the manifest, a CPU-only MATLAB smoke check loaded both
registered directories, validated both MAT cases, called
`benchmark_case(mode,name,pwd,false,saved_result)` to verify frozen-input
compatibility, and passed self-comparisons for both accepted suites. A separate
SHA-256 check verified all six manifest artifacts (four MAT files and two metadata
files). Documentation links and whitespace checks, including new untracked files,
passed; tracked production source remained unchanged.

## Accepted baseline timings and numerical results

| Mode / case | First run (s) | Warm median (s) | Warm range (s) | Maximum occupied norm error |
|---|---:|---:|---:|---:|
| Quick relaxation | 101.481973 | 61.620746 | 58.953476–67.411696 | 1.11e-16 |
| Quick driven nonuniform | 0.252757 | 0.185428 | 0.184337–0.187670 | 2.22e-16 |
| Current relaxation | 100.651682 | 62.960307 | 60.847234–63.458042 | 1.11e-16 |
| Current driven nonuniform | 0.268517 | 0.192756 | 0.187531–0.208832 | 2.22e-16 |

Both baselines have exactly zero vacant spins at all saved times. PNG diagnostics
for quick driven and current relaxation were visually inspected; masks, initial
and final arrays, ordinary/class means and norm errors render correctly. Array
indices are schematic and are not labeled as crystallographic coordinates.

The two preliminary current runs overlapped (see below), so their timings are
not accepted performance references. The table uses the later isolated current
candidate. The quick candidate also ran without another benchmark worker, yet
its relaxation median was 110.098784 s versus the original baseline's 61.620746 s
with unchanged production source. No code speedup or slowdown is inferred from
this environmental/runtime variation; reported ratios are diagnostic only.

Timing includes the unchanged integration implementation, including GPU transfers
and in-memory saved states. First-run times may include compilation/warm-up. Timings
are evidence of this environment, not universal performance thresholds. Setup and
save/plot times are separate in metadata; the comparator does not gate correctness
on timing or source-hash equality.

## Interruptions and limitations

Several MATLAB launch attempts failed before execution with WSL
`UtilBindVsockAnyPort: socket failed 1`; retries succeeded. The WSL launcher for the first current-mode attempt (`current-baseline/`)
exited with code 143 while its directory still showed a running suite. WSL
`nvidia-smi` reported GPU access blocked by the operating system; a subsequent
Windows MATLAB `gpuDevice(1)` probe succeeded. However, independent review later
found that the original Windows child continued and completed successfully. Its
10:21:32–10:29:06 run overlapped `current-baseline-retry/` (10:24:19–10:31:56).
Both outputs are preserved, and their deterministic trajectories remain valid;
their timing measurements are confounded and excluded from accepted performance
references. The later sequential candidate run supplies an isolated current
reference. A wrapper exit alone is not sufficient evidence that Windows MATLAB
has stopped; check Windows process start times and output completion before retry.
The cause of the launcher interruption was not established.

No simulation or failure is represented as successful without output evidence.
Raw baselines are external; cloning Git alone is insufficient to run comparisons.
The manifest records their locations and artifact hashes for transfer/archival.
No CI, production output-path refactor, CPU solver fallback or new dependencies
were added. Cross-GPU/driver/MATLAB tolerances were not calibrated. No full 5 ns
switching run or broader scientific suite was rerun because production source is
unchanged. These software regressions supplement, not replace, field-gradient,
convergence, seam/topology and published-protocol checks for future physics edits.

## Consulted sources

- Current `main.m`, `mn3sn_protocol.m`, initialization/integration scripts,
  `verification/check_protocols.m`, `check_periodic_production.m`,
  `check_coupled_rk4.m`, [periodic production report](periodic-production.md), and
  [protocol report](reference-protocols.md) define unchanged software behavior.
- ASTRO's benchmark README, runner, comparator and baseline metadata in
  `/mnt/d/Git_control/spin_dynamics/atomistic_v2/benchmark/`; its
  [GitHub README](https://github.com/zhuzibn/ASTRO/blob/master/benchmark/README.md)
  was checked. Mn3Sn's comparator explicitly gates intermediate saved states.
- Equations workspace under
  `/mnt/c/Users/zzf-m/OneDrive/presentations/equations`: `AGENTS.md`, `README.md`,
  `spintronic_methods.tex`, `sections/hamiltonian.tex`
  (`eq:atomistic-hamiltonian`) and `sections/llg_equation.tex`
  (`eq:llg-gilbert`, `eq:ll-with-flt`). These informed preservation of flux-density
  units, pair-counting and damping cross terms; no equations were changed.
- Local paper `papers.jsonl`/`paper-catalog.md` were searched for PRB 109, 134433.
  The indexed summary resolves outside the preauthorized `_knowledge` root and
  was not read. No publication claim or physical mapping was derived from metadata.
  No new literature evidence or scientific agreement is claimed by this task.
