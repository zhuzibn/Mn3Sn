# Configuration, production results and protected output

CQ-01–04 implementation, approved 2026-11-19 against
`01683afd458ddef90a27a973d7e1c6b546224955`. This migration preserves the existing
fixed-length LLGS calculation and accepted deterministic benchmark. It adds one
public configuration path, an explicit production result and protected output.
The user authorized implementation after reviewing `plan.md`; no commit, push,
baseline promotion or other-repository edit is part of this work.

## Interface and compatibility

```matlab
cfg=mn3sn_default_config('fig5a',300); % same preset/value semantics as mn3sn_protocol
cfg.protocol.natomW=8;
cfg.protocol.natomL=8;
cfg.rng_seed=7;                       % [] retains the default shuffle policy
state=mn3sn_initialize(cfg);          % CPU only; no output or GPU selection
result=mn3sn_run(cfg,'D:/Git_control/spin_dynamics/simulation_results/Mn3Sn/new-run');
```

`main.m` calls the relaxation defaults and `mn3sn_run(cfg)` without clearing the
workspace. Omitted output allocates a fresh `run-*` directory in the sibling
`simulation_results/Mn3Sn` root. Explicit `''` or `[]` returns in memory; a relative
path resolves from the caller's original directory. Caller cwd, path and CPU RNG
state are restored on success and failure. Reusing the current selected GPU does
not reset it. Requesting another device changes selection; restoration of GPU
selection is not promised.

The public configuration has one nested `protocol`, containing name, dimensions,
seed, current, field, pulse edges, runtime, timestep, chunk duration and sampling
stride. Material/control inputs are top-level fields in
[`mn3sn_default_config.m`](../mn3sn_default_config.m). Derived coefficients and
counts belong to `state.config` / `result.config`; they are not editable duplicate
inputs. `protocol.name` records preset origin; explicit field overrides are the
actual input. Preset names must match exactly, including case. Nonempty preset
values must be full, real, finite double scalars; `mn3sn_protocol` and the default
factory reject other types with `Mn3Sn:ProtocolValue` before seed assignment or
field conversion. Omitted/empty values retain their preset defaults.

[`mn3sn_validate_config.m`](../mn3sn_validate_config.m) rejects missing/unknown
fields, malformed numeric inputs, nonunit polarizations, invalid dimensions and
inconsistent time grids. Physical/numerical inputs use full real doubles; flags
also accept logicals. Initialization additionally checks derived coefficients
for overflow and validates the actual initial state. Runtime must contain an
integer number of chunks, and a chunk an integer number of updates, using the
existing relative grid tolerance of 1e-10. Require at least two updates per chunk:
the unchanged integrator reads the last three allocated states. Sampling stride
is a positive integer. Pulse edges retain their original rounded update indices,
with `pulse_end=Inf` allowed; effective indices are recorded in `result.solver`.

| Selection | Supported behavior |
|---|---|
| `rk4=1` | Coupled deterministic RK4; only supported solver |
| `bc=0` | Periodic W>=4 divisible by four; even L>=4 |
| `bc=1` | Retained legacy open code graph, W,L>=2; not a certified finite-device model |
| `DMIenable`, `SOT_DLT` | Validated 0/1; effective DMI/DL amplitude follows the switch |
| `initial_state=[]` | Original generated seed |
| `initial_state` supplied | Finite `(W,L,3)` doubles, unit occupied norms within 1e-12, exactly zero vacancies; no silent normalization or restart interpretation |
| `rng_seed=[]` / integer | Shuffle / twister seed in [0,2^32-1]; actual states recorded |
| Heun, predictor-corrector, FL SOT, STT, thermal, dipoles, pinning, distribution I/O, restart flag, wall initialization | Rejected before GPU integration; legacy source is retained, not repaired |

Material/drive overrides remain numerical inputs, not a claim of scientific
validation for every value. The default material numbers, protocol values and
effective benchmark configuration are exactly preserved. `T=100 K`, STT parameters
and inactive TM coefficients remain recorded legacy quantities; disabled paths
do not acquire physical support by appearing in the schema.

## Production result schema 1

One variable named `result` is saved in `result.mat` (`-v7.3`).

| Fields | Meaning |
|---|---|
| `schema_name`, `schema_version` | `Mn3Sn production result`, 1; distinct from the benchmark contract |
| `requested_config`, `config` | Validated public inputs and exact effective legacy inputs/constants/counts |
| `atomtype`, `axes`, `classes`, `occupied_indices` | Site mask, onsite axes, a/b/c initialization classes and MATLAB column-major occupied indexing |
| `generated_initial`, `initial` | Generated and actual starting `(W,L,3)` arrays |
| `mmx`, `mmy`, `mmz` | Raw double `(W,L,time)` trajectories |
| `t`, `saved_steps`, `actual_updates`, `pulse_active_updates`, `solver` | Physical clock, counts, effective pulse gates and RK4/normalization convention |
| `invariants`, `units`, `array_order`, `site_meaning`, `limitations` | Raw norm/vacancy checks and interpretation |
| `output_root`, `setup_seconds`, `integration_seconds` | Actual destination and separate setup/integration timing |
| `provenance` | UTC timestamp, MATLAB/toolboxes/OS/GPU, RNG policy/states, Git query status and source identity |

Provenance includes names, exact bytes and SHA-256 hashes of every root `.m` file,
including all executed local helpers. A changed/untracked source is identifiable
without a commit. Git HEAD/branch/status are reported only if Git's top-level
directory equals this source root. A snapshot beneath another repository reports
that distinction, rather than borrowing the ancestor's HEAD. Source bytes are
preserved for replay; this is not a restart loader or an automatic environment
reconstruction feature.

Before saving, the runner verifies finite arrays, shape, exact initial frame,
occupied norm error <=1e-12, exactly zero vacancies, saved steps, update counts and
clock alignment within 8 eps(runtime). It holds pulse gates over RK stages and
normalizes only the combined update, as before.

Compatibility changes are deliberate: there is no whole-workspace save or
`finalnew.mat` overwrite. Old `mmx/mmy/mmz`, `t` and `saved_steps` map to fields of
`result`; `atomtype_` maps to `result.atomtype`; parameters/protocol map to
`result.config`; initial component arrays map to `result.initial(:,:,1:3)`.
Intermediate solver fields are not exported. Historical MAT files are untouched.

## Output and baseline ownership

`mn3sn_prepare_output` canonicalizes the requested path and atomically claims a
fresh directory. It rejects an existing file/directory, source-tree or OneDrive
location, either registered baseline tree (including descendants and resolvable
aliases), and an explicitly selected benchmark input tree. Failed validation
happens before GPU selection/integration. A later execution failure can leave its
newly claimed directory; inspect it and use another fresh path for retry.

Ordinary saved runs need the tracked baseline manifest to establish protection;
they do not need the external baseline MAT files. Ordinary benchmark candidates
need both manifest and input artifacts. Missing references fail before output
creation. Explicit empty third-argument benchmark bootstrap is retained only for
separately authorized baseline generation. No baseline was generated or promoted
in this task. A source-only snapshot without a manifest can run in memory.

The benchmark uses the public initializer and in-memory runner, then projects
results into its existing schema 1. The comparator, validator, manifest, exact
input/configuration/clock checks, full-trajectory thresholds and four timed repeats
are unchanged. The accepted current-size reference is intentionally stored in the
historically named `current-candidate` directory; its manifest identity governs.

## Model and evidence boundary

This is **Class 3: a simplified spatial atomistic code graph, fixed-length LLGS**.
Every occupied Mn site evolves independently; a/b/c are initialization classes,
not chemical species. Occupation is 1 for Mn and 3 for vacancy. Defaults remain
8×8 slots/48 spins, periodic J/2 and D/2 bonds, uniform onsite K, yz-plane seeds
and axes, and +x SOT polarization. RE coefficients apply to all occupied sites.

The preserved fields use B in T, moments in J/T, energies in J and gamma in
rad/(s T), with $\mathbf B_i=-(1/\mu_i)\partial E/\partial\mathbf m_i$.
Onsite energy is $-K(\mathbf e_i\cdot\mathbf m_i)^2$, giving the existing factor
$2K/\mu$. The periodic exchange sum counts directed pairs with the energy's
one-half factor. All LLGS damping-dependent torque cross terms, unnormalized
coupled RK stages, final projection, pulse gates and sampled clocks are unchanged.
The effective default coefficients/units remain those in the
[README](../README.md#current-conventions-and-defaults) and
[reference alignment report](reference-protocols.md).

Immutable implementation files: `constantfile.m`, `distrib.m`,
`systemgeneration.m`, `integrate_llg.m`, `atom_llgs_rhs.m`, `dipole_calc.m` and the
three legacy `atomgpu*` kernels. SHA-256 checks cover these, comparison code,
manifest, six accepted external artifacts and archived reference files.

Physical lattice vectors/basis, Sn/stacking, crystallographic site labels,
exchange-shell assignment/normalization and physical DMI/edge registry remain
**REFERENCE MISSING**. Verification below certifies the implemented code graph,
not those physical mappings, finite devices or electronic/AHE predictions.

Planning consulted the shared [CQ catalog][catalog] and [adoption matrix][matrix],
then current source and baselines. The matrix was not treated as authoritative or
edited. ASTRO configuration/result code supplied software examples only. Equation
references under `/mnt/c/Users/zzf-m/OneDrive/presentations/equations` were
`AGENTS.md`, `README.md`, `spintronic_methods.tex`, `references.tex`,
`sections/hamiltonian.tex` (`eq:atomistic-hamiltonian`,
`eq:atomistic-effective-field`, `eq:atomistic-pair-fields`),
`sections/llg_equation.tex` (`eq:llg-gilbert`, `eq:ll-with-flt`) and the Mn3Sn section
of `sections/material_parameters.tex`. They informed preservation of pair, unit
and Gilbert conventions; they do not establish this graph's physical geometry.

The local paper index identified Xu et al., PRB 109, 134433 (2024),
[DOI 10.1103/PhysRevB.109.134433](https://doi.org/10.1103/PhysRevB.109.134433),
and the two APL records in README. At planning time, their indexed full-text/summary
paths resolved outside the then-authorized `_knowledge` root and were not opened.
The subsequently updated instructions register `D:/Git_control/paper-knowledge-base`
as the processed knowledge base; linked OneDrive artifacts retain a separate access
boundary. No new physics claim required a further literature lookup here. Reviewed EE237
presentation records `week11-sot-publish` (summary and extracted slide 19) and
`week15a-afm` (summary and extracted slides 25–27), under
`/mnt/c/Users/zzf-m/OneDrive/presentation-knowledge/notes/ee237/`, supplied secondary
torque context/source discovery only. No source PDF, online paper, new literature
parameter or physical construction was adopted.

## Validation record

Validation completed in the fresh external directory
`D:/Git_control/spin_dynamics/simulation_results/Mn3Sn/cq01-04-20261119-implementation`
(called `task` below). It contains the approved plan, original source/init/short
trajectory capture, logs, candidate suites and final reports. Environment:
Windows MATLAB R2024b Update 1, Parallel Computing Toolbox 24.2, NVIDIA GeForce
RTX 3060 Laptop GPU, invoked through the WSL MATLAB wrapper. Python checks use
the standard library only. No dependency installation or network access was used.

The accepted artifacts were available locally: all six manifest hashes matched
before and after. The archived three-spin reference, MATLAB/toolboxes, GPU and
existing g++ compiler were available. No baseline or prerequisite was substituted.
`before/source/` and `before/capture.mat` preserve the old implementation separately
from the new defaults; `capture_before.m` records its exact harness. That one-time
capture used 35 preset combinations and four short trajectories: 2×2 and 5×5 open
driven controls, 8×8 periodic driven, and a short periodic relaxation with stride 5.

| Gate | Observed result / artifact under task |
|---|---|
| Old/new initialization | All 35 presets and four controls match exactly; masks, axes, generated spins and complete effective configuration |
| Configuration | 60 invalid-input cases rejected; valid coefficient overrides and caller isolation pass (`interface-and-fixtures-final.log`) |
| Production interface | Four captured trajectories have maximum difference 0; exact clock/configuration/mapping; save/load, wrapper, RNG/source provenance, same-GPU array preservation, DMI/DL off-switch equivalence and protected paths pass (`production-tests-final/checks.json`) |
| Accepted quick/current suites | All four cases pass with trajectory, occupied/class mean and time differences 0 (`quick-comparison-final.json`, `current-comparison.json`) |
| Raw benchmark invariants | Relaxation: 20,000 updates/101 frames, norm error 1.11e-16; driven: 60 updates/61 frames/31 active updates, norm error 2.22e-16; vacancies exactly zero; all timed repeats agree exactly |
| Comparator/protection | Original 21 corruptions rejected; invalid baseline, missing artifacts/manifest, protected subtrees/aliases, source change and invalid timing cases pass (`comparator-tests/checks.json`) |
| Fresh static fixtures | 2×2 open and 4×4/8×8/12×12/20×20 periodic; maximum source-field discrepancy 5.69e-14 T; independent gradient discrepancy <=1.59e-8 T (<1e-6 T gate); ten C++/GPU comparisons <=2.23e-16 (`focused/verified-gates.json`) |
| Three-spin fields | Source discrepancy <=2.85e-14 T; scaled gradient residual <=2.23e-11 (`focused/static.json`) |
| Protocol/periodic dynamics | 35 presets and eight pulse/chunk cases pass; reference discrepancy <=9.16e-16, chunk difference 0. Twelve nonuniform driven cases across 48/108/300 spins pass, maximum discrepancy 4.00e-15; five invalid periodic sizes rejected (`focused/protocols.json`, `periodic-production.json`) |
| Coupled RK4 | Gilbert relative residual 3.20e-16; orders 3.80159 and 3.91363 within [3.7,4.3]; preprojection norm errors decrease from 2.12e-6 to 1.82e-11 with refinement; CPU/GPU difference 0; three driven replay errors <=4.43e-14 (`focused/production-gates.json`, `production-initial.json`) |
| Viewer data | All numeric data unchanged; only main/protocol source identity and the new default-factory entry differ. Node passes ten open sizes/183,850 operator entries (`viewer-comparison.json`) |
| Source analysis/preservation | 17 changed/new MATLAB files checked, no `SYNER`; 133 nonfatal diagnostics recorded in `checkcode.json`. All 32 protected entries match, excluding intentionally changed old main/protocol; no legacy result file created. Final quick/current root source snapshots match the working source byte-for-byte. |

The analyzer diagnostics concern script/eval aliases invisible to static analysis,
existing exporter/check code, and style/performance suggestions. The wrapper test's
return and input are used through `eval`. They were reviewed; no blanket suppression
or unrelated formatting cleanup was applied.

Timing remains informational. Final warm medians (baseline / candidate seconds):
quick relaxation 61.620746 / 94.906629; quick driven .185429 / .166475; current
relaxation 62.960307 / 55.058529; current driven .192756 / .189709. Suites ran
sequentially without another agent GPU workload. Timing varied substantially;
there is no performance acceptance target or speedup claim. An earlier quick
candidate also matched exactly but preceded the same-device reuse correction;
`quick-candidate-final` and `current-candidate` contain the final runner, SHA-256
`176a88b3db0b75f3e1f4ab13eabb69d273dac6d6efde6090014553628ee175ae`.

One intermediate interface test failed because a source snapshot inherited its
parent Git repository. Root identity checking fixed this; the final no-checkout
test passes. Independent review also identified repeated indexed GPU selection
as a reset; reuse and the retained-array test cover that issue. Earlier logs are
retained, with final artifacts identified above. No physical source or accepted
artifact changed to make a check pass.

The browser spot-check was unavailable: Chrome MCP could not find its debugging
port, and a separate headless launch failed at the WSL socket bridge. No rendered
browser pass is claimed. Viewer scripts/layout and displayed field-expression
sources are unchanged; the refreshed source hashes/excerpts and all numeric data
passed the MATLAB/Node checks. Multi-GPU selection, an unavailable default device,
other MATLAB/GPU environments, full 5 ns switching sweeps and rejected optional
physics were not tested. Existing physical-mapping gaps remain unresolved.

## Commands and rerun procedure

These are the executed procedures. Use new output child names for a repeat; keep
the original `before` capture and accepted baselines read-only. MATLAB commands
were invoked as `matlab -wait -logfile 'D:/.../logfile.log' -batch "..."` from the
repository root; Python/Node commands below also run there. All final numerical
checks used the same root `.m` bytes as the two final benchmark snapshots.

```matlab
task='D:/Git_control/spin_dynamics/simulation_results/Mn3Sn/cq01-04-20261119-implementation';
addpath('benchmark');
manifest=jsondecode(fileread('benchmark/baseline_manifest.json'));
% Executed sequentially; quick was repeated after the device-reuse correction.
run_baseline_benchmark('quick',fullfile(task,'quick-candidate-final'));
assert(compare_to_baseline(fullfile(task,'quick-candidate-final'),fullfile(pwd,manifest.quick.directory)));
run_baseline_benchmark('current',fullfile(task,'current-candidate'));
assert(compare_to_baseline(fullfile(task,'current-candidate'),fullfile(pwd,manifest.current.directory)));
check_benchmark(fullfile(pwd,manifest.quick.directory),fullfile(task,'comparator-tests'));
addpath('verification');
check_configuration(fullfile(task,'before'));
check_production_interface(fullfile(task,'production-tests-final'),fullfile(task,'before'));
```

Before preparing fixtures, create `task/focused/native` from all current root
`.m` files and copy `.sisyphus/parameter-alignment-20261102-01/reference` unchanged
to `task/focused/reference`, using the
[existing preparation recipe](periodic-production.md#repeat-the-checks). No prior
simulation outputs are copied. Then execute these steps in order:

```matlab
addpath('verification/spatial_relaxation');
prepare_periodic(fullfile(task,'focused'));
```

```bash
python3 verification/spatial_relaxation/periodic_six_states.py ../simulation_results/Mn3Sn/cq01-04-20261119-implementation/focused gates
```

```matlab
root=pwd; out=fullfile(task,'focused');
addpath('verification'); check_three_spin(out);
addpath(root); check_protocols(out); check_periodic_production(out);
```

```bash
python3 verification/three_spin.py probes ../simulation_results/Mn3Sn/cq01-04-20261119-implementation/focused
```

```matlab
addpath('verification'); check_coupled_rk4(fullfile(task,'focused'));
addpath('docs'); export_visualization_data('docs/mn3sn-data.js');
eval(fileread(fullfile(task,'run_checkcode.m'))); % retained exact 17-file analyzer harness
```

```bash
node docs/check_visualization_data.cjs
git diff --stat
git diff --check
```

A Node VM loaded `before/mn3sn-data.js` and the refreshed data, serialized both
to plain objects, and deep-compared every key except `sources`/`hashes`; only
`main`, `mn3sn_protocol` and `mn3sn_default_config` may differ in source identity.
Python standard-library SHA-256 compared the immutable entries in
`protected-before.json` and every root `.m` byte in both final source snapshots.
Documentation links and line endings were checked separately. Existing ignore
rules already cover the generated/cache directories; `.gitignore` is unchanged.

An independent read-only verifier reviewed the approved scope, source/tests,
protected hashes and documentation. Its findings were resolved without concurrent
writers or changes to physics. The shared matrix owner can use this report as
evidence for CQ-01–04 adoption; no shared matrix update is included here.

## CQ-02 follow-up — 2026-11-20

The user approved the scoped plan against HEAD
`e477b6f1e508098997262a06bd600c9760ce9d21`. The shared matrix still described the
pre-migration production interface; current source already had the configuration
factory, validator and callable runner. This follow-up closes two narrower gaps.

`mn3sn_protocol` previously accepted generic numeric preset values. For example,
`mn3sn_default_config('fig5a',uint16(300))` returned an integer zero field because
of arithmetic in the input type. The production validator then rejected that
configuration before integration. Single and sparse values similarly produced
configurations rejected downstream. A single early guard now enforces the
existing full, real, finite double scalar contract, reporting
`Mn3Sn:ProtocolValue`. Valid double values, selection errors and omitted/empty
defaults retain their original behavior; no automatic type conversion is added.

`check_configuration()` now runs current-interface checks without external data
or GPU calls. It reports `Historical initialization comparison: not requested`.
The explicit `check_configuration(before_root)` form additionally checks the
original 35 preset records and four open/short initialization controls. Missing,
empty or structurally incomplete captures fail; they cannot turn into a skipped
comparison. The historical argument remains separate from production inputs.

The only production-source edit is the preset guard in `mn3sn_protocol.m`.
`verification/check_configuration.m` retains its prior checks and adds preset
coverage. README documents the types and two check modes; the changelog records
the fix. `docs/mn3sn-data.js` retains all numerical data while updating the
protocol source/hash through the existing exporter for the
[open cell](unit-cell.html) and [interaction viewer](structure-exchange.html).
The default factory, validator, initializer, runner, physics kernels, benchmark
and output guards retain their source bytes. Valid-input numerics and defaults
are unchanged. Existing ignore rules already cover generated
and cache directories; no instruction or ignore-file change is needed.

### Follow-up evidence and commands

All new artifacts are in
`D:/Git_control/spin_dynamics/simulation_results/Mn3Sn/cq02-20261120-6o2898i6`
(`task` below). `approved-plan.md`, `source-before/`, `factory-before.mat`,
`check_configuration-before.m`, `mn3sn-data-before.js` and
`protected-before.json` preserve the approved scope and pre-change evidence.
The historical fixture remains read-only at
`D:/Git_control/spin_dynamics/simulation_results/Mn3Sn/cq01-04-20261119-implementation/before`.
No capture or accepted baseline was generated, replaced or promoted.

Environment: Windows MATLAB `24.2.0.2740171 (R2024b) Update 1`, Parallel Computing
Toolbox 24.2, NVIDIA GeForce RTX 3060 Laptop GPU (compute capability 8.6, graphics
driver 610.60, MATLAB toolkit 12.2, 6 GiB). The GPU preflight reported no running
GPU processes; only this sequential benchmark was launched. Node and Python's
standard library supplied data/source checks without dependency installation.

MATLAB commands use the installed WSL launcher with Windows paths:

```bash
matlab -wait -logfile 'D:/Git_control/spin_dynamics/simulation_results/Mn3Sn/cq02-20261120-6o2898i6/before.log' -batch "run('D:/Git_control/spin_dynamics/simulation_results/Mn3Sn/cq02-20261120-6o2898i6/capture_before.m')"
matlab -wait -logfile 'D:/Git_control/spin_dynamics/simulation_results/Mn3Sn/cq02-20261120-6o2898i6/preset-patch.log' -batch "run('D:/Git_control/spin_dynamics/simulation_results/Mn3Sn/cq02-20261120-6o2898i6/check_preset_patch.m')"
matlab -wait -logfile 'D:/Git_control/spin_dynamics/simulation_results/Mn3Sn/cq02-20261120-6o2898i6/interface.log' -batch "run('D:/Git_control/spin_dynamics/simulation_results/Mn3Sn/cq02-20261120-6o2898i6/check_interface.m')"
matlab -wait -logfile 'D:/Git_control/spin_dynamics/simulation_results/Mn3Sn/cq02-20261120-6o2898i6/regression.log' -batch "run('D:/Git_control/spin_dynamics/simulation_results/Mn3Sn/cq02-20261120-6o2898i6/check_regression.m')"
```

These retained harnesses record the exact calls and fixtures. `capture_before.m`
ran before the source edit and must not be rerun as a replacement capture. The
interface harness first resets the MATLAB path and runs a source-only copy with
neither `capture.mat` nor a benchmark directory. It then checks four invalid
historical arguments and the valid original capture. The regression harness
compares CPU initialization with all four accepted quick/current inputs, exports
viewer data externally, runs Code Analyzer on the two edited MATLAB files, then
runs the unchanged quick GPU suite and asserts its full-trajectory comparator.
Use fresh child paths for future output-producing checks.

| Check | Observed result / artifact |
|---|---|
| Pre-change configuration and reproduction | Existing 35 historical presets, four initialization controls and 60 invalid configurations pass; six single/integer/sparse constructor-versus-validator mismatches reproduced (`before.log`) |
| Valid factory preservation | All 35 complete configurations equal the saved pre-change outputs exactly; 80 focused typed rejections pass (`preset-patch.log`) |
| Standalone current checks | 35 valid presets, 132 rejected preset calls and 60 rejected configurations pass, including effective overrides and workspace isolation; no historical comparison claimed (`interface.log`) |
| Explicit historical checks | Same current checks pass plus exact comparison of 35 records/four controls; nonexistent, malformed, numeric-empty and empty-text fixture arguments all fail (`interface.log`) |
| Accepted CPU inputs | Complete effective configuration, initial spins and mapping match all four quick/current cases (`regression.log`) |
| Viewer export and consistency | Fresh MATLAB impulse checks pass for ten open sizes; Node deep comparison finds only `sources.mn3sn_protocol` and `hashes.mn3sn_protocol` changed; current-source/excerpt and 183,850 operator-entry checks pass (`viewer-check.json`) |
| Code Analyzer | No syntax errors in either changed MATLAB file; protocol has no diagnostics, configuration check has three style/performance diagnostics (`checkcode.json`) |
| Quick GPU trajectory regression | Both cases pass; every trajectory component, occupied/class mean and time matches the accepted baseline exactly. All three warm repeats of each case have zero trajectory difference. Norm errors are 1.11e-16 (relaxation) and 2.22e-16 (driven); vacancies remain exactly zero (`quick-comparison.json`, `quick/metadata.json`) |
| Source and baseline preservation | All 27 protected files other than the intentionally edited protocol retain their hashes, including the manifest and all six accepted artifacts. Every root MATLAB file matches the final benchmark snapshot byte-for-byte (`preservation-after.json`) |

The Node comparison loaded old/new `window.MN3SN_DATA` in separate VM contexts,
converted them to plain JSON objects, removed only the two protocol source/hash
entries, and required deep equality before copying the generated file into the
repository. The new excerpt and hash were also checked against the actual source.
Then `node docs/check_visualization_data.cjs` verified the installed snapshot.
No browser/layout change was made. Code Analyzer's diagnostics are the existing
compact `try,` style plus the same style in the new helper and one cell-array
construction suggestion; no diagnostics were suppressed.

All four MATLAB harness processes exited 0. The quick suite used the existing
20,000-update/101-frame relaxation and 60-update/61-frame driven cases, including
the latter's 31 active pulse updates. The unchanged comparator enforces full
trajectory/mean error <=1e-10, occupied norm error <=1e-12, exact vacancy zeros,
configuration/mapping/count compatibility and clock roundoff limits. The observed
trajectory and clock errors were zero. Current-size inputs were checked on CPU;
current-size trajectories were not rerun.

Timing is informational: warm medians were 97.338478 s for relaxation and
0.252224 s for the driven case, versus accepted 61.620746 s and 0.185429 s. No
performance improvement is claimed and no tolerance or repeat count was changed.
The final protocol SHA-256 is
`618adcdd017655a342ae3048c4e9aa83c89eb13ad3b3e7c66d04b6ea317085cb`;
the configuration check is
`3ba25e404408a2268ee5324365c51294fac44c17f866675420ac8142cf53b58f`.
Local documentation links, preserved LF line endings, `git diff --stat`,
`git diff --check` and the untracked-file list were reviewed. The completed task
has been removed from the active plan; its approved specification remains in the
external artifact directory named above and its outcome is recorded here.

The model remains Class 3, fixed-length LLGS on the existing spatial code graph.
No material value, field/torque expression, damping cross term, boundary rule,
initialization, time grid or solver algorithm changes. The equation files and
labels listed in the earlier Model and evidence boundary section were consulted
to preserve those conventions. The local processed
[Mn3Sn synthesis](/mnt/d/Git_control/paper-knowledge-base/topics/mn3sn-spin-orbit-torque-switching.md)
was read for model-scope context; the index located the PRB 109, 134433 (2024)
record and the two APL records cited in README. Indexed per-paper summaries,
extractions and PDFs remain outside the authorized processed KB root and were
not opened. The reviewed EE237 Week 15a summary/slides 25–27 and Week 11
summary/slide 19 supplied secondary context only. No source deck, new literature
parameter or online evidence was adopted.

Physical lattice/basis/stacking, crystallographic site assignment, exchange-shell
normalization and physical DMI/edge registry remain **REFERENCE MISSING**. These
existing limitations do not affect the software-only guard. Full switching
sweeps, current-size timing, multi-GPU behavior, browser rendering and a new
physics/convergence audit are outside this follow-up. No commit, push, dependency
installation, network access or shared-matrix edit is part of the task.

[catalog]: /home/zhuzibn/.agents/skills/spin-dynamics-development/references/code-quality.md
[matrix]: /mnt/d/Git_control/spin_dynamics/docs/repository-quality-matrix.md
