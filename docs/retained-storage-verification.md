# Output guards and retained storage — 2026-11-20

Revision 1 was approved with “proceed”. After measuring the GPU-resident variant,
the user selected: **“Keep the smaller buffers and previous CPU/GPU execution
placement; record GPU residency as a remaining CQ-03 limitation.”** The selected
implementation retains two working state planes and the requested output frames.
It preserves the existing CPU field/stage arithmetic and GPU RHS evaluation.
The guard now resolves real existing ancestors, including Windows junctions,
before claiming an output directory. Relative paths use MATLAB's current folder.

All selected full-suite and focused gates passed. Results below distinguish
the discarded GPU experiment from the user-selected implementation.

## Unchanged model and contract

The adopted Class 3 masked-array fixed-length LLGS model, parameters, units,
initialization and supported paths remain those in
[the configuration report](configuration-results.md),
[protocol conventions](reference-protocols.md) and
[coupled-RK4 specification](frozen-field-repair.md).
Default 8-by-8 periodic slots contain 48 occupied spins (label 1); vacancies
(label 3) remain exactly zero. Field/neighbor expressions, axes, half-strength
periodic exchange/DMI weights, explicit Gilbert torque cross terms, RK stages
and final occupied-spin normalization are unchanged. No scalar default changed.

The initial/final states remain included once, with
`saved_steps=unique([0:savetstep:actual_updates actual_updates])` and
`t=saved_steps*tstep`. Pulses retain zero-based completed-update gates, inclusive
on/exclusive off, held across each step's stages. Schema 1 and output ownership
are unchanged. Unsupported Heun/PC4, thermal, dipole, pinning, restart, FL SOT
and STT remain rejected. No checkpoint or new observable was added.

Physical crystallographic coordinates, shell normalization and legacy open
reciprocity remain unresolved. The reviewed local
`paper-knowledge-base/topics/mn3sn-spin-orbit-torque-switching.md` (“Scope and
evidence base”, “Comparison”) and EE237 week15a summary (slides 25–28) supply
context, not the missing spatial mapping. Shared equation sources
`spintronic_methods.tex`, `symbols.tex`, `sections/llg_equation.tex`
(`eq:llg-gilbert`, `eq:ll-with-flt`) informed unchanged conventions. No physical
source was substituted, and no paper PDF or equation source was modified.

## Implementation and measured scope

Only `integrate_llg.m` and `mn3sn_prepare_output.m` changed among captured root
MATLAB files. Rolling storage applies to RK4; unsupported legacy branches retain
their historical storage. The two-plane current state is sampled after each
update only when its global index is requested; output cadence cannot change
stage evaluation. The final remainder and chunk endpoints are retained correctly.

The old path guard admitted a destination through an actual Windows junction
into a disposable source fixture. No actual source/baseline tree was used as a
write target. The new test rejects ten protected direct/relative/case/junction
paths covering source, registered/selected reference trees and a simulated
OneDrive directory. Fresh nested output still succeeds, and existing output fails.

At the default configuration, instrumented allocations before stepping were:

| Array storage | Before | Selected implementation |
|---|---:|---:|
| Three working component buffers | 3,073,536 bytes | 3,072 bytes |
| Retained host trajectories | 155,136 bytes | 155,136 bytes |

The measurement counts actual double-array elements (8-by-8, 2,001 dense states
versus two working states). It excludes other field/stage buffers, process peak
memory and allocator overhead. On a 24-update/three-chunk control, instrumented
explicit transfers remain 288 host uploads and 288 GPU gathers before and after.
Implicit arrayfun argument transfers are not included in that inventory.
Thus selective history retention is implemented, while **GPU residency remains
a CQ-03 limitation by explicit user choice**.

The discarded variant passed four short old/new controls with zero error, but
its full quick relaxation took 358.860247 s first, then 273.820663 s and
170.835304 s in two completed warm runs. The user stopped it before completion.
Its incomplete `quick-candidate` is marked aborted; no full pass or median is
claimed. The selected placement's complete timing results are informational,
not a controlled cross-run speedup claim.

## Evidence, commands and current results

Source snapshot: `main`, `ef36eac5ad2e5d76d1fd750d424bf0463ad98c5b` plus this
working-tree patch. Windows MATLAB R2024b Update 1, PCT 24.2, NVIDIA RTX 3060
Laptop GPU, driver 610.60, CUDA toolkit 12.2. No competing GPU suites ran.

External task root:
`D:/simulation_results/Mn3Sn/cq-followup-20261120-v1/`.
Source captures, approved revision and harnesses:
`D:/Git_control/spin_dynamics/.sisyphus/cq-followup-20261120-v1/`.
The pre-migration interface fixture at
`../simulation_results/Mn3Sn/cq01-04-20261119-implementation/before/`
was read only.

- Before quick regression: both accepted cases match exactly; three warm
  repeats per case agree. Four new multi-chunk/stride/open/periodic captures
  are in `dense-edges.mat`.
- Guard/contract checks: ten protected alias cases, 21 deliberate comparator
  mutations, missing/invalid baseline and timing/source policies pass.
- Configuration: 35 presets, 132 invalid preset calls, 60 invalid configurations
  and all historical initialization records pass.
- Selected rolling storage: all four old/new controls have maximum error zero;
  dense-stride subsets match exactly. The production-interface check passes
  four historical trajectories, save/load, main wrapper, provenance and guards.
- Selected quick and current suites both pass with zero trajectory/mean/class
  differences and exact clocks. Each case completed three repeat runs with zero
  repeat differences. Relaxation warm medians were 49.458076 s (quick) and
  59.353176 s (current); timing remains informational.
- Five periodic preparations and ten independent C replays pass: maximum field
  discrepancy `5.68e-14 T`, replay component discrepancy `2.220446049250313e-16`.
- Static field/energy checks pass (source discrepancy `2.84e-14 T`, scaled
  gradient residual `2.22e-11`). Viewer export/Node checks pass all 183,850
  operator entries across ten open sizes; deep comparison proves only the
  integrator source excerpt/hash changed, with all numerical viewer data equal.
- Final focused checks pass: 35 presets/eight pulse cases (maximum independent
  reference error `9.16e-16`, chunk difference zero), 12 nonuniform driven
  periodic cases (maximum `4e-15`), implicit-Gilbert torque residual `3.2e-16`,
  coupled orders 3.80159 and 3.91363, CPU/GPU discrepancy zero, and three coupled
  production replays (maximum `4.43e-14`). Existing gates were not relaxed.

The exact selected MATLAB harness is `finish_mn_selected.m`, invoked with:

```bash
matlab -wait -logfile 'D:/simulation_results/Mn3Sn/cq-followup-20261120-v1/selected-final.log' -batch "run('D:/Git_control/spin_dynamics/.sisyphus/cq-followup-20261120-v1/finish_mn_selected.m')"
```

It runs retention and production-interface tests, allocation/explicit-transfer
probes, quick/current suites with asserted full-trajectory comparisons, and
`prepare_periodic` on fresh `focused/native` and unchanged reference fixtures.
The subsequent independent preparation/check commands were:

```bash
python3 verification/spatial_relaxation/periodic_six_states.py /mnt/d/simulation_results/Mn3Sn/cq-followup-20261120-v1/focused gates
python3 verification/three_spin.py probes /mnt/d/simulation_results/Mn3Sn/cq-followup-20261120-v1/focused
```

```matlab
addpath('verification','docs');
out='D:/simulation_results/Mn3Sn/cq-followup-20261120-v1/focused';
check_three_spin(out); % CPU static/reference checks before C probes
export_visualization_data('docs/mn3sn-data.js');
check_protocols(out); check_periodic_production(out); check_coupled_rk4(out);
```

`finish_mn_focused.m` executed the last three checks with log
`focused-final.log` (exit 0). `static-viewer.log` records CPU checks/export;
`node docs/check_visualization_data.cjs` and a Node deep comparison of old/new
viewer data passed. See the original
[configuration recipe](configuration-results.md#commands-and-rerun-procedure)
for the fixture preparation details. Logs and MAT/JSON evidence remain in the
task root. Accepted references, their manifests and numerical helper sources
remain immutable. Read-only independent reviews covered each storage stage and
the real-parent guard; they are separate from executed test evidence.
