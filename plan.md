# Mn3Sn pending work

CQ-01–04 adoption is complete within the approved deterministic scope. The
[implementation and validation report](docs/configuration-results.md) records the
configuration/result interface, unchanged accepted baselines, exact comparisons,
scientific checks, provenance and limitations. Browser rendering and multi-GPU
behavior were not verified. No commit or push was made. The approved detailed
plan is preserved with the external validation artifacts named in that report.

## Pending backlog

These items summarize existing limitations, not new feature commitments. Select
and scope the relevant item before starting substantial source changes.

| Area | Unresolved scope | Required check before completion |
|---|---|---|
| Physical geometry and interaction mapping | Source the lattice/basis, crystallographic labels, exchange-shell normalization, directed DMI and physical edge termination. Periodic code-graph verification does not establish this mapping. | Reconcile the sourced Hamiltonian with bonds and field derivatives; verify topology, seams and energy gradients. Unsupported physical details remain **REFERENCE MISSING**. |
| Legacy open-boundary model | Resolve the physical justification of receiver-dependent exchange/DMI weights. Periodic production is already reciprocal; the legacy open branch remains a separate issue. | Agree the intended open model before changing weights; check reciprocal bonds, field/energy consistency and edge/size dependence. See [A02 evidence](A02-Nonreciprocal_exchange.md). |
| Solver and torque restoration | Heun and predictor-corrector remain defective. The public API now rejects them; independent FL/DL and STT restoration remains separate work. | Restoring support requires selected-path, vacancy, flag and convergence checks; preserve explicit Gilbert damping cross terms. |
| Optional paths | Thermal noise, dipoles, pinning, loaded distributions, restart and domain-wall initialization remain unvalidated. | Scope each requested path separately, including its physical conventions, occupied-spin/vacancy invariants and output contract. Domain walls require sourced sublattice states and physical coordinates. |
| CI | No CI exists. The local validated configuration/output API is implemented. | Scope CI separately; do not infer that local GPU checks are portable automated coverage. |
| Performance | Repeated transfers and workspace storage remain candidates for profiling. | Profile an agreed representative case first; any optimization must preserve fields, trajectories, clocks and vacancy handling. |

Additional publication studies (orientation scans, domain-wall motion or octupole
observables) require a separate requested scope and acceptance criteria. Existing
switching comparisons do not establish those capabilities.

## Completed work: reports only

The following work is complete within each report's documented scope. It is not pending work.

- [CQ-01–04 configuration, production results and protected output](docs/configuration-results.md)
- [Deterministic regression and timing benchmark](docs/benchmark.md)
- [Structure viewers and usage](README.md#interactive-structure-and-interactions)
- [Coupled RK4 repair and figure checks](docs/frozen-field-repair.md)
- [Reference parameters, protocols and time accounting](docs/reference-protocols.md)
- [Periodic six-state relaxation](docs/periodic-six-state.md), [switching comparison](docs/periodic-prb2024.md), and [production verification](docs/periodic-production.md)
- [Model variants and earlier open-boundary studies](README.md#model-variants-exchange-anisotropy-and-boundaries)
- [Archived audit, source locators and completed task plans](docs/audit-history.md)
