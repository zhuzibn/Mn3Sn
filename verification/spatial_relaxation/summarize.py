from pathlib import Path
import array
import csv
import json
import sys

p = Path(sys.argv[1]).resolve()


def data(f):
    a = array.array('d')
    a.frombytes(f.read_bytes())
    return a


rows = []
for model in sorted(p.glob('*.model')):
    size = model.stem
    n = int(model.read_text().split()[0])
    for mode in (['current'] if size == '2x2' else ['current', 'reciprocal']):
        tag = f'{size}-{mode}'
        meta = json.loads((p / f'{tag}-5fs-5ns.json').read_text())
        trace = list(csv.DictReader((p / f'{tag}-5fs-5ns.csv').open()))
        a = data(p / f'{tag}-5fs-5ns.bin')
        b = data(p / f'{tag}-2.5fs-1ns.bin')
        assert len(a) == 51 * 3 * n and len(b) == 11 * 3 * n
        err = max(abs(x-y) for x, y in zip(a, b))
        assert err < 1e-5
        eq = data(p / f'{tag}-5fs-5ns-final.bin')
        recovery = data(p / f'{tag}-recovery-final.bin')
        pert = data(p / f'{tag}-perturbed-start-final.bin')
        recerr = max(abs(x-y) for x, y in zip(eq, recovery))
        perterr = max(abs(x-y) for x, y in zip(eq, pert))
        assert recerr < 1e-5
        # A perturbed nonequilibrium seed may select a different stable basin.
        pertmeta = json.loads((p / f'{tag}-perturbed-start.json').read_text())
        recmeta = json.loads((p / f'{tag}-recovery.json').read_text())
        assert max(pertmeta['max_torque_T'], recmeta['max_torque_T']) < 1e-6
        settled = next(float(r['time_ps']) for i, r in enumerate(trace)
                       if all(float(s['max_torque_T']) < 1e-6 for s in trace[i:]))
        assert meta['max_torque_T'] < 1e-6
        assert meta['drift_since_sample'] < 1e-5 and meta['norm_error'] < 1e-12
        energies = [float(r['energy_over_mu_per_spin']) for r in trace]
        rise = max(y-x for x, y in zip(energies, energies[1:]))
        if mode == 'reciprocal':
            assert rise < 1e-10
        rows.append(dict(size=size, N=n, mode=mode, settled_ps=round(settled),
                         torque_T=meta['max_torque_T'], drift_100ps=meta['drift_since_sample'],
                         refinement_error=err, recovery_error=recerr,
                         perturbed_start_difference=perterr,
                         perturbed_start_torque_T=pertmeta['max_torque_T'],
                         preprojection=meta['preprojection_error'],
                         max_sampled_energy_rise=rise))
(p / 'summary.json').write_text(json.dumps(rows, indent=2))
with (p / 'summary.csv').open('w') as f:
    w = csv.DictWriter(f, fieldnames=rows[0].keys())
    w.writeheader()
    w.writerows(rows)
short = []
for size in ['5x5', '21x22']:
    for mode in ['current', 'reciprocal']:
        a = data(p / f'{size}-{mode}-5fs-10ps.bin')
        b = data(p / f'{size}-{mode}-2.5fs-10ps.bin')
        c = data(p / f'{size}-{mode}-0.5fs-10ps.bin')
        assert len(a) == len(b) == len(c)
        e1 = max(abs(x-y) for x, y in zip(a, c))
        e2 = max(abs(x-y) for x, y in zip(b, c))
        short.append(dict(size=size, mode=mode, error_5fs=e1, error_2_5fs=e2, ratio=e1/e2))
        assert e1 < 1e-5
(p / 'short-convergence.json').write_text(json.dumps(short, indent=2))
print(json.dumps(dict(summary=rows, short=short), indent=2))
