"""Six-orientation tests using the previously source-validated spatial diagnostic."""
import argparse
import array
from concurrent.futures import ThreadPoolExecutor
import csv
import hashlib
import json
import math
from pathlib import Path
import shutil
import subprocess

parser = argparse.ArgumentParser()
parser.add_argument('source', type=Path)
parser.add_argument('output', type=Path)
parser.add_argument('phase', choices=['reference', 'spatial', 'analyze'])
args = parser.parse_args()
source, out = args.source.resolve(), args.output.resolve()
out.mkdir(exist_ok=True)


def read(path):
    a = array.array('d')
    a.frombytes(path.read_bytes())
    return a


def model_data(path):
    lines = path.read_text().splitlines()
    n = int(lines[0].split()[0])
    axes = [list(map(float, line.split()))[:3] for line in lines[1:n+1]]
    angles = [round(math.degrees(math.atan2(v[1], v[2]))) % 360 for v in axes]
    classes = [[30, 270, 150].index(a) for a in angles]
    return n, classes


if args.phase != 'analyze':
    assert (source/'verified-gates.json').exists(), 'Run the source replay gates first.'
    exe = source/'spatial'
    models = [source/'2x2.model'] if args.phase == 'reference' else sorted(
        p for p in source.glob('*.model') if p.stem != '2x2')
    jobs = []
    for model in models:
        shutil.copyfile(model, out/model.name)
        n, classes = model_data(model)
        for phi in range(30, 360, 60):
            seed = (60-phi) % 360
            init = array.array('d')
            for c in classes:
                theta = math.radians(seed+120*c)
                init.extend([0, math.sin(theta), math.cos(theta)])
            initial = out/f'{model.stem}-phi{phi}-initial.bin'
            assert not initial.exists(), 'Use a fresh output directory.'
            initial.write_bytes(init.tobytes())
            for mode in (['current'] if n == 3 else ['current', 'reciprocal']):
                for perturb in ([0] if n == 3 else [0, .001]):
                    tag = f'{model.stem}-{mode}-phi{phi}-p{perturb:g}'
                    jobs.append([str(exe), str(out/model.name), mode, '5', '5000', '5',
                                 str(out/tag), str(initial), str(perturb)])
    (out/f'{args.phase}-commands.json').write_text(json.dumps(jobs, indent=2))
    (out/f'{args.phase}-provenance.json').write_text(json.dumps({
        'executable_sha256': hashlib.sha256(exe.read_bytes()).hexdigest(),
        'source_hashes': json.loads((source/'source-hashes.json').read_text()),
        'seed': 'theta_a = 60 - phi; theta_b = theta_a + 120; theta_c = theta_a + 240',
        'perturbation': 'epsilon*[sin(1.7*i),cos(2.3*i),sin(.9*i)], then normalize'}, indent=2))

    def run(command):
        r = subprocess.run(command, capture_output=True, text=True, check=True)
        print(r.stdout.strip(), flush=True)

    with ThreadPoolExecutor(max_workers=3) as pool:
        list(pool.map(run, jobs))

rows = []
for path in sorted(out.glob('*-phi*-p*.json')):
    meta = json.loads(path.read_text())
    tag = path.stem
    size, mode, phi, perturb = tag.split('-')
    phi, perturb = int(phi[3:]), float(perturb[1:])
    n, classes = model_data(out/f'{size}.model')
    trajectory = read(out/f'{tag}.bin')
    trace = list(csv.DictReader((out/f'{tag}.csv').open()))
    assert len(trajectory) == len(trace)*3*n
    counts = [classes.count(c) for c in range(3)]
    records = []
    final_subs = None
    final_spin = None
    for frame, scalar in enumerate(trace):
        spin = [trajectory[frame*3*n+3*i:frame*3*n+3*i+3] for i in range(n)]
        subs = [[sum(spin[i][k] for i in range(n) if classes[i] == c)/counts[c]
                 for k in range(3)] for c in range(3)]
        balanced = [sum(s[k] for s in subs)/3 for k in range(3)]
        raw = [sum(s[k] for s in spin)/n for k in range(3)]
        def angle(v):
            return math.degrees(math.atan2(v[1], v[2])) % 360 if math.hypot(v[1], v[2]) > 1e-10 else float('nan')
        # Coherent AFM phase from the a-sublattice spin, mapped to reference labels.
        order = (60-math.degrees(math.atan2(subs[0][1], subs[0][2]))) % 360
        records.append([float(scalar['time_ps'])/1000, angle(balanced), order, angle(raw),
                        math.sqrt(sum(x*x for x in balanced)), float(scalar['max_torque_T'])])
        final_subs, final_spin = subs, spin
    with (out/f'{tag}-angles.csv').open('w') as f:
        w = csv.writer(f)
        w.writerow(['time_ns', 'balanced_phi_deg', 'order_phi_deg', 'site_mean_phi_deg',
                    'balanced_magnitude', 'max_torque_T'])
        w.writerows(records)
    row = dict(size=size, N=n, mode=mode, initial_phi_deg=phi, perturb=perturb,
               counts=counts, final_phi_deg=records[-1][1], order_phi_deg=records[-1][2],
               site_mean_phi_deg=records[-1][3], balanced_magnitude=records[-1][4],
               max_torque_T=meta['max_torque_T'], norm_error=meta['norm_error'],
               energy_meV_per_spin=meta['energy_over_mu_per_spin']*(3*9.274e-24)/(1.602e-22),
               sublattice_spins=final_subs)
    assert row['norm_error'] < 1e-12
    # Report stationarity rather than dropping failed cases or forcing six destinations.
    row['stationary'] = row['max_torque_T'] < 1e-6
    rows.append(row)
(out/'six-state-summary.json').write_text(json.dumps(rows, indent=2))
print('Analyzed', len(rows), 'trajectories.', flush=True)
