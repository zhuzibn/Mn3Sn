# Project Guidelines & Documentation

All contributors (including AI Agents) must follow these protocols. Record meaningful changes concisely under the relevant sections.

---

## Section 1: Source Code Changes (Changelog)
Update `CHANGELOG.md` for meaningful behavior changes or significant fixes,
following the global changelog policy. Group related edits into one entry.

Include meaningful features, breaking changes, and setup or migration impacts.
Skip typo fixes, formatting, behavior-preserving internal refactors, test-only
changes, and routine documentation cleanup unless project conventions require them.

Append one concise entry per logical change to the bottom of this section.
Use the fields below where useful; omit inapplicable fields:

### [vX.X.X] | YYYY-MM-DD
- **Type:** (feat | fix | change)
- **Description:** Resulting behavior and user impact.
- **Motivation:** Why was this change necessary?
- **Reference:** Optional issue, commit, or detailed project documentation.
- **Test Strategy:** How was this verified?

### 2026-11-02
- **Type:** fix
- **Description:** RK4 now recomputes coupled fields from every full-lattice stage state, preserving the existing LLGS torque coefficients and final occupied-spin normalization.
- **Motivation:** Frozen fields could produce inaccurate transients and incorrect switching destinations even at small timesteps.
- **Reference:** [Repair and Figures 2, 4(a), 5 verification](docs/frozen-field-repair.md). Reference parameters are explicitly matched; the ferromagnet is excluded.
- **Test Strategy:** Static field/energy checks, production GPU replay, coupled timestep convergence, and the documented three-spin figure suite.

### 2026-11-02 | Reference parameters and figure protocols
- **Type:** change
- **Description:** Align Mn anisotropy, exchange, gamma, film magnetization and thickness with the verified three-spin reference; add named relaxation and Fig. 2/4(a)/5 presets. Spatial receiver weights remain unchanged.
- **Motivation:** The previous anisotropy field was half-strength and SOT was approximately 450 times the reference at equal current. Exact pulse protocols also require removing the missing update and duplicate/drifting sample conventions at chunk boundaries.
- **Compatibility:** Default output now has 101 frames spanning exactly 0–10 ps; figure presets have 20,001 frames spanning 0–5 ns. Existing MAT files retain their old clocks. Select the preset inside main.m; its existing output file behavior is unchanged.
- **Reference:** [Parameter and protocol verification](docs/reference-protocols.md).
- **Test Strategy:** Native/reference coefficient and field-gradient checks, all preset choices, production GPU pulse/chunk comparisons, figure regression suite and source-mask preservation.

### 2026-11-02 | Continuous README angle trajectories
- **Type:** fix
- **Description:** Figures 2 and 5 now unwrap angles through ±180° and fit all trajectories within shared panel limits. README captions explain equivalent orientations modulo 360°.
- **Motivation:** Signed-angle wrapping and deliberate line breaks made complete switching paths look truncated.
- **Test Strategy:** Re-export saved trajectories with MATLAB; check plotted samples against raw angles modulo 360°, continuity, series counts and axis coverage, then inspect the PNGs. No dynamics rerun.

### 2026-11-03 | Simulation results outside OneDrive
- **Type:** change
- **Description:** New simulation data, MATLAB logs and generated simulation figures use fresh directories under `D:\simulation_results\Mn3Sn`. The driven anisotropy diagnostic rejects new output paths inside OneDrive; unmodified production scripts require an external source copy because their save path remains fixed.
- **Migration:** Relocated this investigation's two generated run folders after verifying every file's SHA-256. Historical reference datasets retain their existing locations.
- **Reference:** [Unequal-anisotropy switching comparison](docs/anisotropy-switching.md).

### 2026-11-03 | Relocation to D:\spin_dynamics
- **Type:** change
- **Description:** Move the checkout and existing local datasets to `D:\spin_dynamics\Mn3Sn`, and the two external run folders to `D:\spin_dynamics\simulation_results\Mn3Sn`. Update current output instructions and cross-drive reference links.
- **Test Strategy:** Verify every relocated file with SHA-256 and preserve the Git branch, HEAD and pre-existing working-tree changes before updating path documentation. No numerical source changes or simulations.

### 2026-11-03 | Relocation into Git_control
- **Type:** change
- **Description:** Move the checkout to `D:\Git_control\spin_dynamics\Mn3Sn` and its external results under `D:\Git_control\spin_dynamics\simulation_results\Mn3Sn`. Update current output instructions and reproduction commands; retain historical provenance paths.
- **Test Strategy:** Compare branch, HEAD, status and binary diffs with the pre-move snapshot before path edits; check documentation diffs. No numerical source changes or simulations.


### 2026-11-09 | Periodic production model and README figures
- **Type:** change
- **Description:** Make complete periodic boundaries the default for relaxation and figure presets, using 8×8 slots (48 spins), reciprocal Jintra/2 exchange, antisymmetric Dsim/2 DMI and uniform anisotropy. Restore periodic seam bonds and reject dimensions that cannot tile the motif. README figures now show the verified periodic simulations.
- **Compatibility:** Preset array shapes change to 8×8; time grids, pulse timing, material coefficients, RK4 stages and output filenames are unchanged. The legacy open branch remains available with bc=1. The default 10 ps relaxation is still a short smoke example.
- **Reference:** [Production verification](docs/periodic-production.md), [periodic PRB comparison](docs/periodic-prb2024.md).
- **Test Strategy:** Exact model/fixture equivalence to the approved diagnostic; seam, gradient, fresh MATLAB/GPU drive, protocol, coupled-convergence and viewer/source checks. Physical crystallographic mapping remains unresolved.

---

## Section 2: Error & Solution Log
For significant resolved errors, add a concise note here when useful. Skip
routine transient errors and avoid duplicating Section 1. Keep lengthy
investigations in issues or project documentation and link them here.

### [ERROR-ID] | Short Descriptive Title
- **Context:** Where/when did it happen?
- **Symptoms:** Brief description of the failure.
- **Root Cause:** Brief explanation of the cause.
- **Resolution:** Summary of the fix.
- **Prevention:** How to avoid recurrence, when useful.

---

## Section 3: Workflow for Agents
1. **Analyze** the task or error.
2. **Execute** code changes.
3. **Verify** with tests.
4. Update `CHANGELOG.md` for meaningful behavior changes or significant fixes,
   following the global changelog policy. Group related edits into one entry.
   Use Section 1 for changes and Section 2 for useful significant-error notes.
