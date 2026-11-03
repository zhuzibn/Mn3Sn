"""Check distinct spin states, angle agreement, and reference degeneracy in saved runs."""
import array
import json
from pathlib import Path
import sys

p = Path(sys.argv[1]).resolve()
rows = json.loads((p/'six-state-summary.json').read_text())
assert len(rows) == 126
assert all(r['stationary'] and r['norm_error'] < 1e-12 for r in rows)


def read(path):
    a = array.array('d')
    a.frombytes(path.read_bytes())
    return a


def distance(a, b):
    assert len(a) == len(b)
    return max(abs(x-y) for x, y in zip(a, b))


ref = sorted([r for r in rows if r['N'] == 3], key=lambda r: r['initial_phi_deg'])
ref_error = max(abs((r['final_phi_deg']-r['initial_phi_deg']+180) % 360-180) for r in ref)
ref_spread = max(r['energy_meV_per_spin'] for r in ref)-min(r['energy_meV_per_spin'] for r in ref)
assert ref_error < .1 and ref_spread < 1e-10
results = []
for size in ['5x5', '5x6', '9x10', '13x14', '21x22']:
    for mode in ['current', 'reciprocal']:
        selected = sorted([r for r in rows if r['size'] == size and r['mode'] == mode
                           and r['perturb'] == .001], key=lambda r: r['initial_phi_deg'])
        clusters = []
        max_seed_difference = 0
        for r in selected:
            phi = r['initial_phi_deg']
            state = read(p/f'{size}-{mode}-phi{phi}-p0.001-final.bin')
            exact = read(p/f'{size}-{mode}-phi{phi}-p0-final.bin')
            max_seed_difference = max(max_seed_difference, distance(state, exact))
            if not any(distance(state, old) < 1e-6 for old in clusters):
                clusters.append(state)
        # These assert the measured outcome, including failure of six-state retention.
        assert len(clusters) == (6 if mode == 'current' else 2)
        assert max_seed_difference < 1e-6
        results.append(dict(size=size, mode=mode, distinct_states=len(clusters),
                            initial_phi_deg=[r['initial_phi_deg'] for r in selected],
                            final_phi_deg=[r['final_phi_deg'] for r in selected],
                            whole_array_phi_deg=[r['site_mean_phi_deg'] for r in selected],
                            counts_abc=selected[0]['counts'],
                            ideal_vs_perturbed_max_difference=max_seed_difference,
                            energy_meV_per_spin=[r['energy_meV_per_spin'] for r in selected]
                            if mode == 'reciprocal' else None))
report = dict(reference_angle_error_deg=ref_error, reference_energy_spread_meV_per_spin=ref_spread,
              max_torque_T=max(r['max_torque_T'] for r in rows),
              max_norm_error=max(r['norm_error'] for r in rows), spatial=results)
(p/'six-state-checks.json').write_text(json.dumps(report, indent=2))
print(json.dumps(report, indent=2))
