# Mn3Sn repository instructions

## Scope and source of truth

This repository implements spatial classical-spin dynamics for Mn3Sn in MATLAB. Read `README.md` for actual entry points, defaults, requirements, outputs, and solver limitations. Read `plan.md` for pending work and implementation scope. Completed findings and validation live in `docs/`; the historical audit and completed plans are in `docs/audit-history.md`.

The current source defines implemented behavior. The three papers cited in `README.md` define their own published models and protocols; do not assume that this checkout reproduces them. Distinguish spatial atomistic calculations from three-sublattice switching models and isolated orientation scans.

Use the installed `spin-dynamics-development` skill for physics, field, torque, solver, and scientific-validation work. For substantial implementation, update the concrete root `plan.md` and obtain approval of the relevant model and implementation scope before editing source. An audit or documentation request alone does not approve those changes.

## Repository workflow

- Check `git status --short --branch` before edits. Preserve unrelated changes.
- Keep changes limited to the requested behavior and preserve existing source line endings.
- Do not create branches, commits, or pushes without task authorization.
- Do not modify other spin-dynamics repositories or knowledge bases as a side effect.
- Do not install dependencies or access the network without authorization.
- Keep durable instructions here, user-facing behavior in `README.md`, pending scope and acceptance checks in `plan.md`, and completed findings/validation in `docs/`. Remove completed task details from `plan.md` after recording them in the relevant report.
- Update `CHANGELOG.md` for meaningful behavior changes or significant fixes, following the global changelog policy. Group related edits into one entry.
- Check documentation impact, `git diff --stat`, and `git diff --check` before completion. Include new untracked files in the review.

## Entry points and execution

- `main.m` calls `mn3sn_default_config` and `mn3sn_run`; numerical defaults and GPU index 1 are preserved. Configure runs through `cfg`, with one authoritative `cfg.protocol` struct. Do not reintroduce workspace overrides or parse `main.m` for parameters.
- `mn3sn_initialize` supplies CPU initialization and derived coefficients; `mn3sn_run` supplies function-local state to the unchanged scripts. Bind helpers to their source tree and restore caller cwd/path/RNG. Reuse an already-selected matching GPU without resetting its arrays.
- MATLAB with Parallel Computing Toolbox and a supported NVIDIA GPU is required. Thermal noise is rejected by the supported API.
- Do not assume Octave compatibility or a CPU simulation backend.
- Production command: `matlab -batch "run('main.m')"`. It creates a fresh external directory and saves only `result` in `result.mat`. `mn3sn_run(cfg,out)` selects a fresh destination; explicit `''` or `[]` is in-memory only. There is no overwrite option, restart loader or CI workflow.
- Preserve the versioned production schema, raw arrays/invariants, requested/effective settings, mapping, actual pulse indices and source bytes/hashes. Git HEAD alone is insufficient provenance; snapshots must not inherit a parent repository's identity. See `docs/configuration-results.md` for the interface and checks, `docs/reference-protocols.md` for clocks and `docs/frozen-field-repair.md` for integrator evidence.

## Scientific conventions to preserve

- Each occupied site carries a fixed-length spin; arrays use `(W,L,time)`.
- `atomtype_==1` means occupied Mn; `atomtype_==3` means empty. Empty slots remain exactly zero and are excluded from norm checks and physical averages.
- Occupied sites currently use `musRE`, `gamRE`, and the `RE` torque amplitudes. These legacy names do not establish rare-earth/transition-metal composition in Mn3Sn.
- Initial spins and easy axes lie in the spin-space $yz$ plane, with default SOT polarization $+x$. Do not equate array indices with crystallographic coordinates without a sourced mapping.
- `bc=0` is the production default: complete periodic repeats, uniform half-strength exchange/DMI bonds and uniform K. Require W>=4 divisible by four and even L>=4. All protocol defaults use 8x8 slots (48 spins). `bc=1` retains the legacy open weights and terminal cuts; the 2x2 open array remains a verification control.
- Effective fields entering LLGS are in tesla; gamma is rad/(s T), moments are J/T, and interaction coefficients are energies in J. `H` variable names do not imply A/m.
- `constantfile.m` stores `hbar` in eV s. Derive conversions before changing gamma or torque prefactors.
- `Ksim` is the onsite energy coefficient in `-Ksim*(e dot m)^2`; its field uses `2*Ksim/mu`. Keep `tz=30 nm` film thickness separate from `d=0.4 nm` legacy spatial spacing. Native Mn parameters match `Mn3Sn_3atoms`. Periodic bonds use Jintra/2 and Dsim/2; open receiver weights remain legacy. Physical lattice mapping and exchange-shell normalization remain unresolved.
- Do not silently alter interaction signs, pair counting, site scaling, anisotropy factors, moment/volume factors, thickness, damping cross terms, timestep, drive timing, noise generation, or normalization.
- Source every physical lattice, basis, neighbor, DMI orientation, or parameter change. Label unsupported details `REFERENCE MISSING` and block that portion of implementation.
- Reconcile published Hamiltonians with their field derivatives explicitly. Do not transplant the ASTRO ferrimagnetic DMI geometry or material parameters into this model.
- The explicit LLGS kernel already includes damping-dependent SOT/STT cross terms. Preserve them during any approved refactor and verify against the implicit Gilbert equation.

## Solver and optional-path policy

`rk4=1` is the default coupled RK4 path. Recompute fields from the complete unnormalized lattice at each stage through `integrate_llg.m`, evaluate the explicit derivative with `atom_llgs_rhs.m`, and normalize occupied spins only after combining stages. `atomgpurk4.m` is retained as a historical frozen-field kernel, not the active RK4 step. Keep the production GPU replay and coupled-convergence checks in `verification/check_coupled_rk4.m` passing when modifying stages. Validation is restricted to deterministic, unpinned, dipole-free cases; periodic seam/topology checks are in `verification/spatial_relaxation/prepare_periodic.m` and `periodic_six_states.py`, with driven production checks in `verification/check_periodic_production.m`. Their code-graph verification does not certify crystallographic topology.

Heun (`rk4=0`) and predictor-corrector (`rk4=2`) are legacy paths with known defects and are rejected by configuration validation. Do not restore support without scoped tests.

The public API also rejects FL SOT, STT, thermal noise, dipoles, fixed edges, distribution I/O, `loadstartm` and `dwcalc`. Their source remains unvalidated. Any restoration must exercise the selected branch, its interaction with empty slots, flags and output contract. Supported DMI/DL SOT on/off settings must reach the actual coefficients. Keep stochastic conventions separate from deterministic RK4 fixes.

## Benchmark ownership

- `benchmark/run_baseline_benchmark.m` runs quick (4x4) or current (8x8) periodic suites; `benchmark/compare_to_baseline.m` compares complete saved trajectories and reports timing separately.
- Read `benchmark/README.md` before running or modifying the benchmark. Use a fresh external output directory. The runner loads frozen inputs from the manifest's baseline by default; an explicit third argument selects another baseline directory.
- Run timing suites sequentially without competing GPU work. After an unexpected WSL launcher exit, check Windows MATLAB processes and output completion before retrying.
- Treat `benchmark/baseline_manifest.json` and its external baseline artifacts as accepted regression data. Do not regenerate or promote baselines as a side effect of code changes; obtain explicit authorization for replacement.
- Raw baseline data are external and are not included in a Git clone. Missing baseline artifacts/manifest fail ordinary candidate runs; they are not permission to bootstrap replacements. Only an explicitly authorized baseline-generation run may use an empty third argument.
- Preserve full-trajectory comparison, exact configuration/clock checks, zero vacancies, occupied norms and safe-path guards. Use `benchmark/check_benchmark.m` for changes to the benchmark itself. Existing scientific checks remain required for physics or solver changes.

## Validation and result protection

- Save new simulation data, MATLAB logs and generated simulation figures outside OneDrive. Default to `D:\Git_control\spin_dynamics\simulation_results\Mn3Sn` (`/mnt/d/Git_control/spin_dynamics/simulation_results/Mn3Sn` in WSL), with a fresh directory per run. Keep source and written reports in this repository. Preserve historical datasets unless their relocation is explicitly requested.
- Use the narrowest checks appropriate to the requested work. Deterministic regression baselines are registered in `benchmark/baseline_manifest.json`; see `benchmark/README.md` and `docs/benchmark.md`.
- Preserve historical results before establishing a baseline. A record of legacy behavior is not proof of physical correctness.
- Production, the benchmark runner, `check_benchmark` and `check_production_interface` use `mn3sn_prepare_output` to claim fresh directories, reject source/OneDrive/existing destinations and protect registered/selected baseline subtrees before integration. A missing manifest blocks ordinary persistence; a source copy can run in memory. Preserve `finalnew.mat`, `atomtypee.mat` and all historical results. Other scientific/study harnesses write into prepared fixture directories and still require an explicitly reviewed fresh external destination.
- Check shapes, finite occupied spins, occupied norms, exact zero empty slots, sublattice mapping, actual update counts, and alignment of saved states with physical times.
- Each block includes its initial state plus every update. Save shared boundaries once, include t=0 and the final state, and derive `t` from `saved_steps*tstep`. Pulse gates use zero-based completed updates, on-inclusive/off-exclusive, held constant across RK stages. Keep `verification/check_protocols.m` passing for changes to protocols or time accounting.
- For field changes, check reciprocal exchange weights, antisymmetric DMI bonds, term-level analytical limits, and energy-gradient consistency on small deterministic lattices.
- For integration changes, check coupled-system timestep convergence and norm error before projection; a normalized trajectory alone does not establish solver accuracy.
- Verify periodic seams separately from open-edge truncation. Finite domain-wall studies require edge and size checks, a sourced chiral initial profile, and physical coordinates.
- Reproduction claims require the target paper's geometry, parameter normalization, pulses, and observables. Do not claim AHE or electronic-structure predictions from the classical spins alone.
- Record exact commands or harness procedure, MATLAB/toolbox/GPU environment, overrides, results, and omissions. Never invent a pass, tolerance, or published agreement.
- When source work introduces generated outputs, review appropriate ignore rules. Do not change `.gitignore` during a read-only audit or documentation-only task unless requested.

Use `$...$` and `$$...$$` for Markdown mathematics intended for Typora.
