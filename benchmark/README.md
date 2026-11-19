# Mn3Sn deterministic benchmark

Run this benchmark before and after a code change to compare saved trajectories
and synchronized GPU runtime. It adapts the [ASTRO benchmark workflow](https://github.com/zhuzibn/ASTRO/tree/master/benchmark)
to Mn3Sn's periodic production model. Production MATLAB files are unchanged.
The baseline records existing behavior; it does not establish physical correctness.

## Run and compare

Run from the repository root with MATLAB, Parallel Computing Toolbox and a
supported NVIDIA GPU. No Python packages or archived verification fixtures are
required. The examples below use MATLAB paths; the installed WSL launcher invokes
Windows MATLAB, so use `D:/...` rather than `/mnt/d/...` inside MATLAB expressions.

```matlab
addpath('benchmark');
manifest=jsondecode(fileread('benchmark/baseline_manifest.json'));
baseline=fullfile(pwd,manifest.quick.directory);
candidate='D:/Git_control/spin_dynamics/simulation_results/Mn3Sn/benchmark/my-change-quick';
run_baseline_benchmark('quick',candidate);
[passed,report]=compare_to_baseline(candidate,baseline);
assert(passed);
```

Use `'current'` and `manifest.current.directory` for the production-size suite.
Use a new candidate directory for every invocation. The runner refuses any
existing output directory, a location inside this repository, or a OneDrive path.
It copies the production source before execution, records its hashes, and runs in
that isolated directory. It never calls `main.m`'s final save or visualization.

A shell equivalent is:

```bash
matlab -wait -batch "addpath('benchmark'); run_baseline_benchmark('quick','D:/Git_control/spin_dynamics/simulation_results/Mn3Sn/benchmark/my-change-quick');"
```

Without a third argument, the runner loads frozen inputs from the mode's accepted
baseline in `baseline_manifest.json`. To use a baseline stored elsewhere:

```matlab
run_baseline_benchmark('current',candidate,'D:/archive/mn3sn-current-baseline');
```

The current mask, axes, class mapping, generated initial state and effective
configuration must match those saved in the selected baseline. For the nonuniform
case, the actual initial arrays are then loaded from that baseline. A configuration
change is an explicit incompatibility, not a valid performance comparison.

## Cases and timing

| Mode | Array slots | Occupied spins |
|---|---|---|
| `quick` | 4x4 | 12 |
| `current` | 8x8 | 48 |

Both modes execute these deterministic, periodic, unpinned, dipole-free coupled
RK4 cases with thermal noise disabled:

| Case | Initial state and drive | Clock |
|---|---|---|
| `relaxation` | Production 121-degree seed; zero current and field | 10 ps, 0.5 fs updates, 1 ps chunks, every 200 updates: 20,000 updates / 101 saved states |
| `driven_nonuniform` | Frozen normalized nonuniform spins; +5e10 A/m2, +0.01 T along y, SOT polarization +x | 5 fs updates, 60 updates, 20-update chunks, every update saved; pulse gate [13,44): 31 active updates / 61 states |

The second case is a software stress test, not a 5 ns switching protocol. The
relaxation gate spans all 20,000 updates even though its drive amplitudes are zero.
Initial states and shared chunk boundaries are saved once; the final state is
included and `t=saved_steps*tstep`. Means exclude vacancies and are not normalized.
Classes a/b/c correspond to seed/seed+120/seed+240 degrees in initialization.

Each case runs once for first-run timing and three more times in fresh function
workspaces. `wait(gpuDevice)` brackets the integration timer. The median and range
of the three warm runs are reported; all repeats must agree within 1e-10 per spin
component. Integration timing includes production field work, GPU transfers and
in-memory trajectory storage. Initialization and disk/plot time are separate.
It is not isolated GPU kernel throughput. The small quick case still performs
20,000 relaxation updates, so it is not an instant smoke test.

Timing is informational. A speedup is reported only for matching recorded MATLAB,
toolbox, OS and GPU environments and valid timing data. Background load and warm-up
still affect runtime; inspect the range before interpreting a small difference.
No fixed performance threshold determines correctness.
Run one benchmark suite at a time and avoid other GPU workloads. Environment
matching cannot detect competing work. If the WSL launcher exits unexpectedly,
check the Windows MATLAB process and the suite log before retrying: a Windows
child can continue after its launcher exits. Never reuse its output directory.

## Result contract and comparison

Each directory contains:

- `relaxation.mat`, `driven_nonuniform.mat`: one `result` structure per case;
- corresponding PNGs: occupation, initial/final z, occupied and class means,
  and norm-error histogram (array indices, not crystallographic coordinates);
- `metadata.json`: completion status, command, environment, Git state, source and
  result SHA-256 hashes, timing and norm summaries;
- `matlab.log` and `source/`: execution log and source snapshot.

Schema version 1 includes configuration, mask, axes, classes, generated and actual
initial states, `mmx/mmy/mmz`, `saved_steps/t`, actual and pulse-active update counts,
observables, timing and provenance. Raw arrays remain double precision. Vacancies
remain exactly zero. Per-case provenance is captured during execution; the final
**suite** completion status is in `metadata.json`.

The read-only comparator requires both cases and completed, compatible suites.
It validates both baseline and candidate, recomputes observables from raw arrays,
and rejects inconsistent stored summaries. It checks exact configuration, inputs,
mask/class/axis mapping, saved steps and update counts. Every saved spin must be
finite, every occupied norm within 1e-12 of one, and every vacancy exactly zero.
Full saved trajectory components and occupied/class means must differ by at most
1e-10; clocks must agree within roundoff. These tolerances were checked on the
recorded environment, not calibrated across different GPUs or MATLAB releases.

`[passed,report]=compare_to_baseline(candidate,baseline)` returns a structured
report; missing/malformed data returns `passed=false` with a diagnostic. Use
`assert(passed)` to make a batch command fail. Git/source changes are recorded and
allowed: the purpose is to compare later implementations. Invalid timing suppresses
speedup reporting without turning a correct trajectory into a failure.

## Baseline ownership and portability

The tracked [manifest](baseline_manifest.json) identifies accepted external MAT
artifacts and their SHA-256 hashes. Raw data, figures and logs stay outside this
repository under its result-storage policy. A Git clone alone does not include
them. Transfer the two MAT files and `metadata.json` for each mode to the recorded
relative location, or supply an explicit baseline directory to the runner and
comparator. Preserve the complete directory when archiving provenance and plots.
Verify transferred artifact hashes against the manifest.

Ordinary runs never update the manifest or overwrite a baseline. Creating a new
baseline requires explicit authorization and review. For that action only, use a
fresh external output path and an explicit empty third argument to generate the
initial fixtures without loading the old baseline:

```matlab
run_baseline_benchmark('quick','D:/.../new-baseline-quick','');
```

Without an existing manifest, the initial bootstrap also works with two arguments.
The bootstrap nonuniform state uses a deterministic trigonometric construction;
its actual arrays are saved for later runs. Validate a fresh candidate, inspect the
comparison and source changes, then deliberately update the manifest. Never loosen
tolerances or replace golden data merely to make a failing change pass.

## Check the benchmark itself

```matlab
addpath('benchmark');
check_benchmark(baseline,'D:/Git_control/spin_dynamics/simulation_results/Mn3Sn/benchmark/new-comparator-tests');
```

This requires a fresh test directory and MATLAB, but no GPU simulation. It tests
21 deliberate corruptions, including an intermediate norm-preserving rotation
with an unchanged final state, invalid baseline data, incomplete suites, stale
summaries, incompatible clocks/configuration and missing files. It also checks
source-change acceptance, invalid timing, output protection and MATLAB syntax.

The source extraction intentionally checks exact `main.m` initialization anchors.
A structural refactor may require a reviewed harness update; it must not silently
fall back to duplicated parameters. See [validation evidence](../docs/benchmark.md)
for commands, environments, timings, interruptions and limitations.

No coverage is claimed for full-length switching, open boundaries, thermal noise,
dipoles, pinning, restart, alternate solvers, physical crystallographic topology,
energy conservation or new material parameters. Existing scientific verification
suites remain necessary when changing those aspects.
