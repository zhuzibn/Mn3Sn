"""Isolated K_i=K*z_i/4 trial. Python standard library; no production edits.

Arguments: source_operator_folder previous_six_state_folder fresh_output check|run|analyze
Run prepare_anisotropy.m before check. Each phase must pass before the next.
"""
import argparse
import array
from concurrent.futures import ThreadPoolExecutor
import csv
import hashlib
import itertools
import json
import math
from pathlib import Path
import random
import shutil
import subprocess


def read(path):
    data = array.array('d')
    data.frombytes(path.read_bytes())
    return data


def dot(a, b):
    return sum(x*y for x, y in zip(a, b))


def cross(a, b):
    return [a[1]*b[2]-a[2]*b[1], a[2]*b[0]-a[0]*b[2], a[0]*b[1]-a[1]*b[0]]


def distance(a, b):
    assert len(a) == len(b)
    return max(abs(x-y) for x, y in zip(a, b))


def angle(v):
    return math.degrees(math.atan2(v[1], v[2])) % 360 if math.hypot(v[1], v[2]) > 1e-10 else math.nan


def angle_error(a, b):
    return abs((a-b+180) % 360-180)


def model(path):
    lines = path.read_text().splitlines()
    header = lines[0].split()
    n = int(header[0])
    J, D, K = map(float, header[3:6])  # J/mu, D/mu, 2K/mu in tesla
    axes = [list(map(float, line.split()))[:3] for line in lines[1:n+1]]
    classes = [[30, 270, 150].index(round(angle(v)) % 360) for v in axes]
    bonds = {}
    for line in lines[n+2:]:
        i, j, a, c = map(float, line.split())
        bonds[int(i), int(j)] = (J/2 if a else 0, math.copysign(D/2, c) if c else 0)
    degree = [sum(i == p and a != 0 for (p, _), (a, _) in bonds.items()) for i in range(n)]
    ki = [K*z/4 for z in degree]

    def energy(s):
        terms = [-.5*ki[i]*dot(s[i], axes[i])**2 for i in range(n)]
        terms += [a*dot(s[i], s[j])-c*cross(s[i], s[j])[0]
                  for (i, j), (a, c) in bonds.items() if i < j]
        return math.fsum(terms)

    def field(s):
        h = [[ki[i]*dot(s[i], axes[i])*x for x in axes[i]] for i in range(n)]
        for (i, j), (a, c) in bonds.items():
            for k in range(3):
                h[i][k] -= a*s[j][k]
            h[i][1] += c*s[j][2]
            h[i][2] -= c*s[j][1]
        return h

    return n, axes, classes, bonds, degree, energy, field


parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('source', type=Path)
parser.add_argument('reference', type=Path)
parser.add_argument('output', type=Path)
parser.add_argument('phase', choices=['check', 'run', 'analyze'])
args = parser.parse_args()
source, reference, out = args.source.resolve(), args.reference.resolve(), args.output.resolve()
here = Path(__file__).resolve().parent
root = here.parent.parent
phis = list(range(30, 360, 60))
conversion = 3*9.274e-24/1.602e-22  # E/mu in T -> meV
models = sorted(source.glob('*.model'))
assert len(models) == 6 and out.is_dir()
hashes = json.loads((source/'source-hashes.json').read_text())
for name, expected in hashes.items():
    assert hashlib.sha256((root/name).read_bytes()).hexdigest() == expected, name


def run(command):
    result = subprocess.run(command, check=True, capture_output=True, text=True)
    print(result.stdout.strip(), flush=True)


if args.phase == 'check':
    assert not (out/'candidate-gates.json').exists(), 'Use a fresh trial directory.'
    subprocess.run(['g++', '-O3', '-std=c++17', '-Wall', '-Wextra', str(here/'spatial.cpp'),
                    '-o', str(out/'spatial')], check=True)
    ref_n, _, ref_classes, _, _, _, _ = model(reference/'2x2.model')
    assert ref_n == 3
    refs = {}
    for phi in phis:
        flat = read(reference/f'2x2-current-phi{phi}-p0-final.bin')
        refs[phi] = {ref_classes[i]: list(flat[3*i:3*i+3]) for i in range(3)}
    gates = []
    for path in models:
        shutil.copyfile(path, out/path.name)
        n, axes, classes, bonds, degree, energy, field = model(path)
        assert set(degree) <= {2, 4}
        for (i, j), (a, c) in bonds.items():
            assert bonds[j, i] == (a, -c) and a > 0 and c != 0
        neighbors = [{j for p, j in bonds if p == i} for i in range(n)]
        triangles = [(i, j, k) for i in range(n) for j, k in itertools.combinations(sorted(neighbors[i]), 2)
                     if i < j < k and k in neighbors[j]]
        for i, j in bonds:
            assert sum(i in t and j in t for t in triangles) == 1
        for i in range(n):
            assert sum(i in t for t in triangles) == degree[i]/2
        patterns = set()
        for t in triangles:
            assert sorted(classes[i] for i in t) == [0, 1, 2]
            a, b, c = sorted(t, key=lambda i: classes[i])
            patterns.add(tuple(math.copysign(1, bonds[p, q][1]) for p, q in [(a, b), (b, c), (c, a)]))
        assert len(patterns) == 1
        weighted_counts = [sum(degree[i]/4 for i in range(n) if classes[i] == c) for c in range(3)]
        assert len(set(weighted_counts)) == 1
        rng = random.Random(20261102)
        s = [[rng.uniform(-1, 1) for _ in range(3)] for _ in range(n)]
        h = field(s)
        gradient_error = 0
        for i in sorted({0, n//3, n//2, n-1}):
            for k in range(3):
                old = s[i][k]
                s[i][k] = old+1e-5
                ep = energy(s)
                s[i][k] = old-1e-5
                em = energy(s)
                s[i][k] = old
                gradient_error = max(gradient_error, abs((ep-em)/2e-5+h[i][k]))
        assert gradient_error < 1e-6
        energies, torque = [], 0
        for phi in phis:
            s = [refs[phi][c] for c in classes]
            energies.append(energy(s)*conversion/n)
            torque = max(torque, max(math.sqrt(dot(v, v)) for v in map(cross, s, field(s))))
        spread = max(energies)-min(energies)
        assert torque < 1e-9 and spread < 1e-10
        replay = []
        for trial in range(2):
            gpu = out/f'{path.stem}-gpu{trial}.bin'
            tag = out/f'{path.stem}-replay{trial}'
            run([str(out/'spatial'), str(path), 'triangle', '5', '.1', '.005', str(tag), str(gpu)])
            actual, predicted = read(gpu), read(tag.with_suffix('.bin'))
            assert len(actual) == len(predicted) == 21*3*n
            error = distance(actual, predicted)
            assert error < 1e-10
            replay.append(error)
        # Unchanged modes must reproduce the previous compiled diagnostic exactly.
        for mode in ['current', 'reciprocal']:
            tag = out/f'{path.stem}-old-{mode}'
            old_tag = out/f'{path.stem}-baseline-{mode}'
            tail = [str(path), mode, '5', '.1', '.005']
            run([str(source/'spatial'), *tail, str(old_tag)])
            run([str(out/'spatial'), *tail, str(tag)])
            assert tag.with_suffix('.bin').read_bytes() == old_tag.with_suffix('.bin').read_bytes()
            assert tag.with_suffix('.csv').read_bytes() == old_tag.with_suffix('.csv').read_bytes()
        gates.append(dict(size=path.stem, N=n, triangles=len(triangles), degree_counts={str(z): degree.count(z) for z in set(degree)},
                          weighted_counts=weighted_counts, gradient_error_T=gradient_error,
                          reference_torque_T=torque, reference_energy_spread_meV_per_spin=spread, replay_errors=replay))
    (out/'candidate-gates.json').write_text(json.dumps(gates, indent=2))
    (out/'diagnostic-hashes.json').write_text(json.dumps({str(p.relative_to(root)): hashlib.sha256(p.read_bytes()).hexdigest()
                                                       for p in here.glob('*') if p.is_file()}, indent=2))
    print('All candidate gates passed.', flush=True)

if args.phase == 'run':
    assert len(json.loads((out/'candidate-gates.json').read_text())) == 6
    assert not (out/'commands.json').exists(), 'Use a fresh directory; existing runs are preserved.'
    jobs, followups = [], []
    for path in models:
        n, _, classes, _, _, _, _ = model(path)
        if n == 3:
            continue
        for phi in phis:
            initial = out/f'{path.stem}-phi{phi}-initial.bin'
            state = array.array('d')
            for c in classes:
                theta = math.radians(60-phi+120*c)
                state.extend([0, math.sin(theta), math.cos(theta)])
            initial.write_bytes(state.tobytes())
            command = [str(out/'spatial'), str(out/path.name), 'triangle']
            for perturb in ['0', '0.001']:
                tag = f'{path.stem}-triangle-phi{phi}-p{perturb}'
                jobs.append(command+['5', '5000', '5', str(out/tag), str(initial), perturb])
            tag = f'{path.stem}-triangle-phi{phi}'
            followups.append(command+['5', '2000', '5', str(out/f'{tag}-recovery'), str(out/f'{tag}-p0-final.bin'), '.001'])
            if path.stem == '21x22':
                followups.append(command+['2.5', '2000', '5', str(out/f'{tag}-refined'), str(initial), '.001'])
    (out/'commands.json').write_text(json.dumps(jobs+followups, indent=2))
    for batch in [jobs, followups]:
        with ThreadPoolExecutor(max_workers=3) as pool:
            list(pool.map(run, batch))

if args.phase == 'analyze':
    rows, groups = [], []
    max_recovery = max_refinement = max_increase = max_preprojection = 0
    _, _, reference_classes, _, _, _, _ = model(reference/'2x2.model')
    max_reference_error = max_sublattice_spread = 0
    for path in models:
        n, _, classes, _, _, energy, _ = model(path)
        if n == 3:
            continue
        counts = [classes.count(c) for c in range(3)]
        finals, energies, angle_errors, start_errors = [], [], [], []
        for phi in phis:
            states = []
            for perturb in ['0', '0.001']:
                tag = f'{path.stem}-triangle-phi{phi}-p{perturb}'
                meta = json.loads((out/f'{tag}.json').read_text())
                scalars = list(csv.DictReader((out/f'{tag}.csv').open()))
                data = read(out/f'{tag}.bin')
                assert len(scalars) == 1001 and len(data) == 1001*n*3 and all(map(math.isfinite, data))
                assert meta['steps'] == 1000000 and meta['norm_error'] < 1e-12 and meta['max_torque_T'] < 1e-6
                assert meta['drift_since_sample'] < 1e-5
                max_preprojection = max(max_preprojection, meta['preprojection_error'])
                sampled_energy = [float(r['energy_over_mu_per_spin']) for r in scalars]
                increase = max(b-a for a, b in zip(sampled_energy, sampled_energy[1:]))
                max_increase = max(max_increase, increase)
                assert increase < 1e-10
                angles = []
                for frame, scalar in enumerate(scalars):
                    assert abs(float(scalar['time_ps'])-frame*5) < 1e-8
                    spin = [data[frame*3*n+3*i:frame*3*n+3*i+3] for i in range(n)]
                    subs = [[math.fsum(spin[i][k] for i in range(n) if classes[i] == c)/counts[c] for k in range(3)] for c in range(3)]
                    balanced = [math.fsum(s[k] for s in subs)/3 for k in range(3)]
                    raw = [math.fsum(s[k] for s in spin)/n for k in range(3)]
                    angles.append([frame*.005, angle(balanced), angle(raw), float(scalar['max_torque_T']), math.sqrt(dot(balanced, balanced))])
                with (out/f'{tag}-angles.csv').open('w') as f:
                    writer = csv.writer(f)
                    writer.writerow(['time_ns', 'balanced_phi_deg', 'site_mean_phi_deg', 'max_torque_T', 'balanced_magnitude'])
                    writer.writerows(angles)
                state = read(out/f'{tag}-final.bin')
                assert distance(state, data[-3*n:]) == 0
                states.append(state)
                e = energy(spin)*conversion/n
                err = angle_error(angle(balanced), phi)
                assert err < .1
                ref = read(reference/f'2x2-current-phi{phi}-p0-final.bin')
                by_class = {reference_classes[i]: ref[3*i:3*i+3] for i in range(3)}
                reference_error = max(abs(spin[i][k]-by_class[classes[i]][k]) for i in range(n) for k in range(3))
                spatial_spread = max(abs(spin[i][k]-subs[classes[i]][k]) for i in range(n) for k in range(3))
                max_reference_error = max(max_reference_error, reference_error)
                max_sublattice_spread = max(max_sublattice_spread, spatial_spread)
                assert reference_error < 1e-5 and spatial_spread < 1e-5
                rows.append(dict(size=path.stem, initial_phi_deg=phi, perturb=float(perturb), counts=counts,
                                 final_phi_deg=angle(balanced), site_mean_phi_deg=angle(raw), sublattice_spins=subs,
                                 energy_meV_per_spin=e, angle_error_deg=err, reference_spin_error=reference_error,
                                 sublattice_spread=spatial_spread, **meta))
                energies.append(e)
                angle_errors.append(err)
            start_errors.append(distance(*states))
            assert start_errors[-1] < 1e-5
            finals.append(states[1])
            prefix = f'{path.stem}-triangle-phi{phi}'
            recovery_meta = json.loads((out/f'{prefix}-recovery.json').read_text())
            assert recovery_meta['max_torque_T'] < 1e-6 and recovery_meta['norm_error'] < 1e-12
            assert recovery_meta['steps'] == 400000 and abs(recovery_meta['time_ps']-2000) < 1e-8
            recovery = distance(states[0], read(out/f'{prefix}-recovery-final.bin'))
            max_recovery = max(max_recovery, recovery)
            assert recovery < 1e-5
            if path.stem == '21x22':
                refined = read(out/f'{prefix}-refined.bin')
                original = read(out/f'{prefix}-p0.001.bin')[:len(refined)]
                assert len(refined) == 401*n*3
                refined_meta = json.loads((out/f'{prefix}-refined.json').read_text())
                assert refined_meta['steps'] == 800000 and refined_meta['norm_error'] < 1e-12
                assert refined_meta['max_torque_T'] < 1e-6
                refinement = distance(original, refined)
                max_refinement = max(max_refinement, refinement)
                assert refinement < 1e-5
        minimum_separation = min(distance(a, b) for a, b in itertools.combinations(finals, 2))
        assert minimum_separation > 1e-6
        spread = max(energies)-min(energies)
        assert spread < 1e-9
        groups.append(dict(size=path.stem, N=n, distinct_states=6, angle_error_deg=max(angle_errors),
                           energy_spread_meV_per_spin=spread, minimum_state_separation=minimum_separation,
                           ideal_perturbed_error=max(start_errors)))
    assert len(rows) == 60
    report = dict(main_trajectories=60, recovery_trajectories=30, refinement_trajectories=6, groups=groups,
                  max_torque_T=max(r['max_torque_T'] for r in rows), max_norm_error=max(r['norm_error'] for r in rows),
                  max_preprojection_error=max_preprojection, max_recovery_error=max_recovery,
                  max_refinement_error=max_refinement, max_sampled_energy_increase_T=max_increase,
                  max_reference_spin_error=max_reference_error, max_sublattice_spread=max_sublattice_spread)
    (out/'six-state-summary.json').write_text(json.dumps(rows, indent=2))
    (out/'trial-checks.json').write_text(json.dumps(report, indent=2))
    print(json.dumps(report, indent=2))
