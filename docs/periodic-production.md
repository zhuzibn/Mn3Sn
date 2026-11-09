# Periodic production model

Production now implements the previously verified periodic graph directly. `main.m` selects `bc=0`; every `mn3sn_protocol` preset uses 8×8 slots with 48 spins, sixteen in each easy-axis class. Periodic dimensions require W>=4 divisible by four and even L>=4. The retained `bc=1` branch preserves its old terminal cuts and receiver weights.

The change restores both periodic seams, retains all sites in the complete repeating motif and uses uniform exchange/DMI factors of one half. Thus each included bond has exchange `Jintra/2` and directed DMI magnitude `Dsim/2`; uniform onsite `Ksim` is unchanged. Every spin evolves independently under the existing coupled RK4 and explicit LLGS derivative. This is the approved mathematical code graph, with physical crystallographic mapping and exchange-shell normalization still unresolved.

All material constants, seed conventions, pulse gates and time grids are unchanged. The default relaxation remains a short 10 ps example at 0.5 fs, saving 101 frames with shape `8 x 8 x 101`; it does not establish complete relaxation. Figure presets run 5 ns at 5 fs and save `8 x 8 x 20001`. Historical MAT files and `main.m`'s fixed `finalnew.mat` output behavior are unchanged. Run source copies in fresh external directories to protect results.

## Production verification

Fresh results are stored in `D:/Git_control/spin_dynamics/simulation_results/Mn3Sn/periodic-production-20261109-01`. MATLAB R2024b Update 1 and an NVIDIA GeForce RTX 3060 Laptop GPU (compute capability 8.6) were used. Python uses its standard library, and the existing C++ diagnostic was rebuilt with g++.

| Check | Result |
|---|---|
| Current-source models versus approved periodic diagnostic | All five model files byte-identical: 2×2 open control, 4×4/8×8/12×12/20×20 periodic |
| Fresh undriven GPU fixtures versus approved diagnostic | All ten trajectories byte-identical |
| Seam topology, class populations, translation and repeated-cell energy | Passed independent wrapped-stencil checks |
| Source field versus matrix field | Maximum 5.69e-14 T |
| Termwise energy-gradient discrepancy | Maximum 1.59e-8 T, below 1e-6 T tolerance |
| Ten C++/production GPU fixture comparisons | Maximum 2.23e-16 per component |
| All protocol choices | 35 presets passed; each has 48 occupied spins |
| Pulse signs, on/off edges, chunks and sampling | Eight cases passed; reference error 9.16e-16/component; chunk error zero |
| Nonuniform driven periodic GPU trajectories | Twelve cases across 48/108/300 spins and all current/field sign combinations; maximum error 4.00e-15/component |
| Invalid periodic dimensions | Five shapes rejected with `Mn3Sn:PeriodicSize` |
| Coupled RK4 regression | Observed orders 3.8016 and 3.9136; CPU/GPU discrepancy zero |
| Retained three-spin production GPU probes | Three driven probes passed, maximum 4.43e-14/component |
| Open-view compatibility | All ten old open snapshots and their measured operators unchanged; Node checked 183,850 operator entries |

The nonuniform check starts from fixed-seed random unit spins, assembles fields independently from verified A/C matrices, integrates the explicit vector LLGS on the CPU and compares every spin at every saved step with the production GPU loop. It spans three chunks and includes pulse-on/off steps inside chunks. All 60 updates and 31 active updates are checked; maximum norm error is 2.23e-16 and vacancies remain exactly zero.

The existing [78-run relaxation study](periodic-six-state.md) and [248-run switching study](periodic-prb2024.md) were not rerun. Exact operator/parameter equivalence, identical undriven fixtures and fresh driven checks establish the connection to production without changing those preserved datasets. Quantitative switching agreement is with the verified numerical reference; published-figure agreement remains qualitative.

An independent read-only audit confirmed the model/fixture equivalence, unchanged open snapshots, source hashes and README image provenance. The new nonuniform driven gate stores scalar discrepancies rather than full trajectories: that audit reviewed its independent implementation and recorded results; raw-state checks cover the fresh undriven binaries and preceding long sweeps. No root result files were overwritten.

## Repeat the checks

From the repository root, create a fresh external directory and source snapshots. The archived reference initialization is reused unchanged; this copies no prior simulation output.

```python
from pathlib import Path
import shutil
root = Path.cwd()
out = root.parent/'simulation_results/Mn3Sn/periodic-production-new'
out.mkdir()
(out/'native').mkdir()
for p in root.glob('*.m'):
    shutil.copy2(p, out/'native'/p.name)
shutil.copytree(root/'.sisyphus/parameter-alignment-20261102-01/reference', out/'reference')
```

Run these commands sequentially from the repository root; use matching Windows/WSL paths:

```bash
matlab -wait -batch "addpath('verification/spatial_relaxation'); prepare_periodic('D:/Git_control/spin_dynamics/simulation_results/Mn3Sn/periodic-production-new');"
python3 verification/spatial_relaxation/periodic_six_states.py ../simulation_results/Mn3Sn/periodic-production-new gates
matlab -wait -batch "out='D:/Git_control/spin_dynamics/simulation_results/Mn3Sn/periodic-production-new'; root=pwd; addpath('verification'); check_three_spin(out); addpath(root); check_protocols(out); check_periodic_production(out);"
python3 verification/three_spin.py probes ../simulation_results/Mn3Sn/periodic-production-new
matlab -wait -batch "addpath('verification'); check_coupled_rk4('D:/Git_control/spin_dynamics/simulation_results/Mn3Sn/periodic-production-new');"
```

`prepare_periodic.m` now executes current production source without text overrides. The 2×2 fixture explicitly sets bc=1 before initialization. The existing open preparation helpers likewise set bc before constructing sites. `check_protocols.m` checks every periodic spin against its class trajectory; `check_periodic_production.m` additionally tests arbitrary nonuniform states.

## README figures and viewers

The two README PNGs are byte-identical copies of `periodic-prb2024.png` and `periodic-prb2024-high-field.png` from the preceding verified switching run. Raw trajectories, the MATLAB plotting source, plot checks and original exports remain in the external run directory documented in [the comparison report](periodic-prb2024.md). No plotting transformations or simulation parameters changed during copying. The original PNGs were visually inspected; the combined figure contains size comparisons and the high-field figure shows the 300-spin cell.

The [periodic companion](periodic-cell.html) now identifies the production model and defaults to 8×8. Its data include current source hashes and match the fresh verified export. The [open cell](unit-cell.html) and [open interaction](structure-exchange.html) views explicitly show only bc=1. Their snapshots were regenerated from current source; the obsolete incomplete-periodic selector and probes were removed.

Refresh open data and check it with:

```bash
matlab -wait -batch "addpath('docs'); export_visualization_data('docs/mn3sn-data.js');"
node docs/check_visualization_data.cjs
```

Refresh periodic data only after the preparation and gates above pass, using the JSON-to-JavaScript command in the [periodic relaxation report](periodic-six-state.md#files-and-rerun-commands). Browser checks exercised 72 periodic class/seam/term cases, selected-neighbor highlighting, size controls and zoom/reset, with no horizontal overflow at 390 px. The mobile rendering was inspected. A desktop screenshot request stalled and was terminated; no desktop screenshot pass is claimed for this promotion.

The open interaction view passed twelve size/term cases and the open cell view loaded a larger repeated patch; desktop DOM checks found no horizontal overflow or console errors. Initial periodic browser-test attempts used the wrong event type for term controls; the corrected tests passed without changing the UI.

## Evidence and limits

The approved specification is in [plan.md](../plan.md#production-periodic-boundaries--2026-11-09). It reuses the prior PRB source review, local Mn3Sn synthesis and EE237 presentation record; no new paper or physical structure was introduced. Equation conventions were checked against `spintronic_methods.tex` and its Hamiltonian, material-parameter, boundary and LLGS sections. Thermal noise, dipoles, pinning, domain-wall initialization, finite-device edges and physical crystallographic claims remain outside this verification.
