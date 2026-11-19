# Mn3Sn pending work

This file tracks unresolved work and the scope/checks needed before implementation.
Current behavior belongs in [README.md](README.md); completed investigations and
validation belong in the reports under `docs/`. Listing an item here does not
approve its implementation. No implementation task is currently in progress.

## Remaining work

These items summarize existing limitations, not new feature commitments. Select
and scope the relevant item before starting substantial source changes.

| Area | Unresolved scope | Required check before completion |
|---|---|---|
| Physical geometry and interaction mapping | Source the lattice/basis, crystallographic labels, exchange-shell normalization, directed DMI and physical edge termination. Periodic code-graph verification does not establish this mapping. | Reconcile the sourced Hamiltonian with bonds and field derivatives; verify topology, seams and energy gradients. Unsupported physical details remain **REFERENCE MISSING**. |
| Legacy open-boundary model | Resolve the physical justification of receiver-dependent exchange/DMI weights. Periodic production is already reciprocal; the legacy open branch remains a separate issue. | Agree the intended open model before changing weights; check reciprocal bonds, field/energy consistency and edge/size dependence. See [A02 evidence](A02-Nonreciprocal_exchange.md). |
| Solver selection and torque flags | Heun and predictor-corrector remain defective; early rejection is still proposed. Independent FL/DL and STT flag behavior and guards need scoped validation. | Exercise selected/rejected paths, vacancies and flag combinations; preserve the explicit Gilbert damping cross terms. Restoring another solver requires separate convergence checks. |
| Optional paths | Thermal noise, dipoles, pinning, loaded distributions, restart and domain-wall initialization remain unvalidated. | Scope each requested path separately, including its physical conventions, occupied-spin/vacancy invariants and output contract. Domain walls require sourced sublattice states and physical coordinates. |
| Results and regression workflow | Fixed-name full-workspace output can overwrite results; no production output-path override, versioned result schema, accepted regression baseline or CI exists. | Define a compatible output/provenance contract and safe paths; preserve historical data and establish a deterministic baseline from the existing verification tools. |
| Performance | Repeated transfers and workspace storage remain candidates for profiling. | Profile an agreed representative case first; any optimization must preserve fields, trajectories, clocks and vacancy handling. |

Additional publication studies (orientation scans, domain-wall motion or octupole
observables) require a separate requested scope and acceptance criteria. Existing
switching comparisons do not establish those capabilities.

## Completed work: reports only

The visualization, frozen-field RK4 repair, parameter alignment, physical clock
repair and periodic production promotion are complete within their documented
scope. They are not pending tasks.

- [Structure viewers and usage](README.md#interactive-structure-and-interactions)
- [Coupled RK4 repair and figure checks](docs/frozen-field-repair.md)
- [Reference parameters, protocols and time accounting](docs/reference-protocols.md)
- [Periodic six-state relaxation](docs/periodic-six-state.md), [switching comparison](docs/periodic-prb2024.md), and [production verification](docs/periodic-production.md)
- [Model variants and earlier open-boundary studies](README.md#model-variants-exchange-anisotropy-and-boundaries)
- [Archived audit, source locators and completed task plans](docs/audit-history.md)

## Keeping this plan short

For the next approved task, add only its goal, scope, implementation steps and
acceptance checks. Once complete, put the findings and validation in the relevant
`docs/` report and remove the detailed task plan from this file. Keep durable
agent rules in `AGENTS.md` and user-facing behavior in `README.md`.
