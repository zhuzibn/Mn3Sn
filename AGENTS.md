# Mn3Sn repository instructions

## Scope and source of truth

This repository implements spatial classical-spin dynamics for Mn3Sn in MATLAB. Read `README.md` for actual entry points, defaults, requirements, outputs, and solver limitations. Read `plan.md` for the audit, scientific evidence, and proposed improvements.

The current source defines implemented behavior. The three papers cited in `README.md` define their own published models and protocols; do not assume that this checkout reproduces them. Distinguish spatial atomistic calculations from three-sublattice switching models and isolated orientation scans.

Use the installed `spin-dynamics-development` skill for physics, field, torque, solver, and scientific-validation work. For substantial implementation, update the concrete root `plan.md` and obtain approval of the relevant model and implementation scope before editing source. An audit or documentation request alone does not approve those changes.

## Repository workflow

- Check `git status --short --branch` before edits. Preserve unrelated changes.
- Keep changes limited to the requested behavior and preserve existing source line endings.
- Do not create branches, commits, or pushes without task authorization.
- Do not modify other spin-dynamics repositories or knowledge bases as a side effect.
- Do not install dependencies or access the network without authorization.
- Keep durable instructions here, user-facing behavior in `README.md`, and findings and implementation checkpoints in `plan.md`.
- Update `CHANGELOG.md` for meaningful behavior changes or significant fixes, following the global changelog policy. Group related edits into one entry.
- Check documentation impact, `git diff --stat`, and `git diff --check` before completion. Include new untracked files in the review.

## Entry points and execution

- Run from the repository root. `main.m`, `constantfile.m`, `systemgeneration.m`, `distrib.m`, `integrate_llg.m`, and `dipole_calc.m` share script workspace state.
- `main.m` clears the workspace, selects GPU device 1, runs integration, and saves the complete workspace to `finalnew.mat`.
- MATLAB with Parallel Computing Toolbox and a supported NVIDIA GPU is required. `normrnd` additionally requires Statistics and Machine Learning Toolbox when thermal noise is enabled.
- Do not assume Octave compatibility or a CPU simulation backend.
- Production command: `matlab -batch "run('main.m')"`. It overwrites `finalnew.mat`; preserve existing data and obtain authorization for overwriting user results.
- There is currently no production output-path override, restart loader, or CI workflow. Select `mn3sn_protocol` inside `main.m` for relaxation or a three-spin figure case. See `docs/reference-protocols.md` for parameter and clock checks and `docs/frozen-field-repair.md` for the earlier integrator evidence.
- The disabled visualization block uses an external `threedplott` function and a machine-specific path. Do not enable or port it unless requested.

## Scientific conventions to preserve

- Each occupied site carries a fixed-length spin; arrays use `(W,L,time)`.
- `atomtype_==1` means occupied Mn; `atomtype_==3` means empty. Empty slots remain exactly zero and are excluded from norm checks and physical averages.
- Occupied sites currently use `musRE`, `gamRE`, and the `RE` torque amplitudes. These legacy names do not establish rare-earth/transition-metal composition in Mn3Sn.
- Initial spins and easy axes lie in the spin-space $yz$ plane, with default SOT polarization $+x$. Do not equate array indices with crystallographic coordinates without a sourced mapping.
- `bc=1` requests open boundaries; `bc=0` requests periodic shifts. Existing coefficient masks prevent treating the latter as verified periodic topology.
- Effective fields entering LLGS are in tesla; gamma is rad/(s T), moments are J/T, and interaction coefficients are energies in J. `H` variable names do not imply A/m.
- `constantfile.m` stores `hbar` in eV s. Derive conversions before changing gamma or torque prefactors.
- `Ksim` is the onsite energy coefficient in `-Ksim*(e dot m)^2`; its field uses `2*Ksim/mu`. Keep `tz=30 nm` film thickness separate from `d=0.4 nm` legacy spatial spacing. Native Mn parameters match `Mn3Sn_3atoms`; spatial receiver weights remain deliberately unresolved.
- Do not silently alter interaction signs, pair counting, site scaling, anisotropy factors, moment/volume factors, thickness, damping cross terms, timestep, drive timing, noise generation, or normalization.
- Source every physical lattice, basis, neighbor, DMI orientation, or parameter change. Label unsupported details `REFERENCE MISSING` and block that portion of implementation.
- Reconcile published Hamiltonians with their field derivatives explicitly. Do not transplant the ASTRO ferrimagnetic DMI geometry or material parameters into this model.
- The explicit LLGS kernel already includes damping-dependent SOT/STT cross terms. Preserve them during any approved refactor and verify against the implicit Gilbert equation.

## Solver and optional-path policy

`rk4=1` is the default coupled RK4 path. Recompute fields from the complete unnormalized lattice at each stage through `integrate_llg.m`, evaluate the explicit derivative with `atom_llgs_rhs.m`, and normalize occupied spins only after combining stages. `atomgpurk4.m` is retained as a historical frozen-field kernel, not the active RK4 step. Keep the production GPU replay and coupled-convergence checks in `verification/check_coupled_rk4.m` passing when modifying stages. Validation is restricted to deterministic, unpinned, dipole-free cases; larger-lattice topology is not certified by the three-site figure comparison.

Heun (`rk4=0`) and predictor-corrector (`rk4=2`) are legacy paths with known defects. Do not use them for production or restore support without scoped tests. Early rejection is a proposed improvement, not behavior already implemented.

Treat thermal noise, dipole modes 1–3, fixed edges, loaded distributions, `loadstartm`, and `dwcalc` as unvalidated options. Any work on them must exercise the selected branch, its interaction with empty slots, and its output contract. Keep stochastic conventions separate from deterministic RK4 fixes.

## Validation and result protection

- Save new simulation data, MATLAB logs and generated simulation figures outside OneDrive. Default to `D:\Git_control\spin_dynamics\simulation_results\Mn3Sn` (`/mnt/d/Git_control/spin_dynamics/simulation_results/Mn3Sn` in WSL), with a fresh directory per run. Keep source and written reports in this repository. For unmodified `main.m`, use an isolated source copy in the external run directory because its save path is fixed. Preserve historical datasets unless their relocation is explicitly requested.
- Use the narrowest checks appropriate to the requested work; there is no accepted local regression baseline yet.
- Preserve historical results before establishing a baseline. A record of legacy behavior is not proof of physical correctness.
- Candidate runs must avoid overwriting `finalnew.mat`, `atomtypee.mat`, and external result files. The current script has no built-in safe output selection, so use an explicitly reviewed temporary harness or isolated copy for audit runs.
- Check shapes, finite occupied spins, occupied norms, exact zero empty slots, sublattice mapping, actual update counts, and alignment of saved states with physical times.
- Each block includes its initial state plus every update. Save shared boundaries once, include t=0 and the final state, and derive `t` from `saved_steps*tstep`. Pulse gates use zero-based completed updates, on-inclusive/off-exclusive, held constant across RK stages. Keep `verification/check_protocols.m` passing for changes to protocols or time accounting.
- For field changes, check reciprocal exchange weights, antisymmetric DMI bonds, term-level analytical limits, and energy-gradient consistency on small deterministic lattices.
- For integration changes, check coupled-system timestep convergence and norm error before projection; a normalized trajectory alone does not establish solver accuracy.
- Verify periodic seams separately from open-edge truncation. Finite domain-wall studies require edge and size checks, a sourced chiral initial profile, and physical coordinates.
- Reproduction claims require the target paper's geometry, parameter normalization, pulses, and observables. Do not claim AHE or electronic-structure predictions from the classical spins alone.
- Record exact commands or harness procedure, MATLAB/toolbox/GPU environment, overrides, results, and omissions. Never invent a pass, tolerance, or published agreement.
- When source work introduces generated outputs, review appropriate ignore rules. Do not change `.gitignore` during a read-only audit or documentation-only task unless requested.

Use `$...$` and `$$...$$` for Markdown mathematics intended for Typora.
