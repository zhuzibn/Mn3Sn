"""Periodic-stencil checks and six-state runs; Python standard library and existing C++ kernel."""
import argparse
import array
from concurrent.futures import ThreadPoolExecutor
import csv
import hashlib
import json
import math
from pathlib import Path
import random
import subprocess
import time

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent.parent
SIZES = ('8x8', '12x12', '20x20')
PHIS = tuple(range(30, 360, 60))
MEV = 3*9.274e-24/1.602e-22  # Convert E/mu in tesla to meV, using constantfile.m.


def write_json(path, value):
    path.write_text(json.dumps(value, indent=2, allow_nan=False)+'\n')


def read_binary(path):
    result = array.array('d')
    result.frombytes(path.read_bytes())
    return result


def model_data(path):
    lines = path.read_text().splitlines()
    h = lines[0].split()
    n, w, length = map(int, h[:3])
    sites = [list(map(float, line.split())) for line in lines[1:n+1]]
    bonds = {}
    for line in lines[n+2:]:
        i, j, a, c = line.split()
        key = int(i), int(j)
        assert key not in bonds
        bonds[key] = (float(a), float(c))
    assert len(bonds) == int(lines[n+1])
    axes = [v[:3] for v in sites]
    classes = [1 if e[1] < -.5 else (2 if e[2] < 0 else 0) for e in axes]
    return dict(n=n, w=w, length=length, parameters=list(map(float, h[3:])),
                axes=axes, classes=classes, bonds=bonds)


def periodic_graph(w, length):
    """Independent directed stencil, in source/receiver offsets, not physical coordinates."""
    assert w % 4 == 0 and length % 2 == 0 and length >= 4
    sites = []
    for c in range(1, length+1):
        for r in range(1, w+1):
            if (r % 4 == 1 and c % 2 == 0) or (r % 4 == 3 and c % 2 == 1):
                continue
            cls = 1 if r % 2 else (0 if (r % 4 == 2) == (c % 2 == 0) else 2)
            sites.append(dict(row=r, col=c, cls=cls))
    lookup = {(s['row'], s['col']): i for i, s in enumerate(sites)}
    # (source row offset, source column offset, directed C sign).
    offsets = {0: [(0, -1, 1), (-1, -1, -1), (0, 1, 1), (1, 0, -1)],
               1: [(-1, 0, 1), (-1, 1, -1), (1, 1, 1), (1, 0, -1)],
               2: [(-1, 0, 1), (0, 1, -1), (1, -1, 1), (0, -1, -1)]}
    reads, bonds = [], {}
    for i, site in enumerate(sites):
        for dr, dc, sign in offsets[site['cls']]:
            rr, cc = site['row']+dr, site['col']+dc
            j = lookup[((rr-1) % w+1, (cc-1) % length+1)]
            assert i != j and (i, j) not in bonds
            bonds[i, j] = (.5, .5*sign)
            reads.append(dict(i=i, j=j, dr=dr, dc=dc, image_row=(rr-1)//w,
                              image_col=(cc-1)//length, A=.5, C=.5*sign))
    return sites, bonds, reads


def energy(m, model, term=None):
    jcoef, dcoef, kcoef = model['parameters'][:3]
    value = 0.
    if term in (None, 'anisotropy'):
        value -= .5*kcoef*sum(sum(x*y for x, y in zip(s, e))**2
                             for s, e in zip(m, model['axes']))
    for (i, j), (a, c) in model['bonds'].items():
        if term in (None, 'exchange'):
            value += .5*jcoef*a*sum(x*y for x, y in zip(m[i], m[j]))
        if term in (None, 'dmi'):
            value -= .5*dcoef*c*(m[i][1]*m[j][2]-m[i][2]*m[j][1])
    return value


def field(m, model, term=None):
    jcoef, dcoef, kcoef = model['parameters'][:3]
    result = [[0., 0., 0.] for _ in m]
    if term in (None, 'anisotropy'):
        for i, e in enumerate(model['axes']):
            dot = sum(x*y for x, y in zip(m[i], e))
            result[i] = [kcoef*dot*x for x in e]
    for (i, j), (a, c) in model['bonds'].items():
        if term in (None, 'exchange'):
            for k in range(3):
                result[i][k] -= jcoef*a*m[j][k]
        if term in (None, 'dmi'):
            result[i][1] += dcoef*c*m[j][2]
            result[i][2] -= dcoef*c*m[j][1]
    return result


def execute(command):
    started = time.monotonic()
    result = subprocess.run(command, capture_output=True, text=True)
    Path(command[6]+'.log').write_text(result.stdout+result.stderr)
    if result.returncode:
        raise RuntimeError(f'{command[6]} failed: {result.stderr}')
    print(f'{Path(command[6]).name}: {time.monotonic()-started:.1f}s', flush=True)


def gates(out):
    gpu = json.loads((out/'gpu-gates.json').read_text())
    assert len(gpu) == 5
    subprocess.run(['g++', '-O3', '-std=c++17', '-Wall', '-Wextra',
                    str(HERE/'spatial.cpp'), '-o', str(out/'spatial')], check=True)
    records, snapshots = [], []
    base_model = model_data(out/'4x4.model')
    base_sites, _, _ = periodic_graph(4, 4)
    rng = random.Random(20261109)
    base_state = {(s['row'], s['col']): [rng.uniform(-1, 1) for _ in range(3)] for s in base_sites}
    reference_energy = energy(list(base_state.values()), base_model)/12
    commands = []
    for tag in ('2x2', '4x4', *SIZES):
        model = model_data(out/f'{tag}.model')
        n = model['n']
        gradient = {}
        if tag != '2x2':
            sites, bonds, reads = periodic_graph(model['w'], model['length'])
            assert set(bonds) == set(model['bonds'])
            assert max(abs(x-y) for pair in bonds for x, y in zip(bonds[pair], model['bonds'][pair])) < 1e-15
            exported = json.loads((out/f'{tag}-sites.json').read_text())
            assert [s['row'] for s in sites] == exported['rows']
            assert [s['col'] for s in sites] == exported['cols']
            assert [s['cls'] for s in sites] == exported['classes'] == model['classes']
            counts = [model['classes'].count(c) for c in range(3)]
            assert counts == [model['w']*model['length']//4]*3
            for (i, j), (a, c) in bonds.items():
                assert bonds[j, i] == (a, -c)
            for i in range(n):
                assert sum(k[0] == i for k in bonds) == 4
                assert [sum(a for (ii, j), (a, _) in bonds.items() if ii == i and sites[j]['cls'] == c)
                        for c in range(3)] == [[0, 1, 1], [1, 0, 1], [1, 1, 0]][sites[i]['cls']]
            lookup = {(s['row'], s['col']): i for i, s in enumerate(sites)}
            for dr, dc in ((4, 0), (0, 2)):
                perm = [lookup[((s['row']-1+dr) % model['w']+1,
                                (s['col']-1+dc) % model['length']+1)] for s in sites]
                assert all(bonds[perm[i], perm[j]] == value for (i, j), value in bonds.items())
            repeated = [base_state[((s['row']-1) % 4+1, (s['col']-1) % 4+1)] for s in sites]
            repeat_error = abs(energy(repeated, model)/n-reference_energy)
            assert repeat_error < 1e-11
            snapshots.append(dict(tag=tag, W=model['w'], L=model['length'], sites=sites, reads=reads))
        else:
            repeat_error = 0
        m = [[rng.uniform(-1, 1) for _ in range(3)] for _ in range(n)]
        # Include seam/corner sites as well as interior sites in the derivative check.
        for term in ('exchange', 'dmi', 'anisotropy'):
            h = field(m, model, term)
            error = 0.
            for i in sorted({0, 1, n//3, n//2, n-2, n-1}):
                for k in range(3):
                    old = m[i][k]
                    eps = 1e-5
                    m[i][k] = old+eps
                    plus = energy(m, model, term)
                    m[i][k] = old-eps
                    minus = energy(m, model, term)
                    m[i][k] = old
                    error = max(error, abs((plus-minus)/(2*eps)+h[i][k]))
            assert error < 1e-6, (tag, term, error)
            gradient[term] = error
        replay = []
        for trial in range(2):
            fixture = out/f'{tag}-gpu{trial}.bin'
            command = [str(out/'spatial'), str(out/f'{tag}.model'), 'current', '5', '.1', '.005',
                       str(out/f'{tag}-replay{trial}'), str(fixture)]
            commands.append(command)
            execute(command)
            actual, expected = read_binary(Path(command[6]+'.bin')), read_binary(fixture)
            assert len(actual) == len(expected) == 21*n*3
            error = max(abs(x-y) for x, y in zip(actual, expected))
            assert error < 1e-10
            replay.append(error)
        records.append(dict(size=tag, gradient_error_T=gradient, replay_errors=replay,
                            repeated_energy_error_T_per_spin=repeat_error))
    write_json(out/'gate-commands.json', commands)
    write_json(out/'verified-gates.json', records)
    write_json(out/'periodic-view-data.json', dict(snapshots=snapshots,
               provenance='prepare_periodic.m impulse matrices checked against periodic_six_states.py directed stencil',
               parameters=base_model['parameters'], physical_geometry='schematic; crystallographic mapping unresolved'))
    print(json.dumps(records, indent=2))


def seed_file(out, tag, phi, offset=0):
    model = model_data(out/f'{tag}.model')
    seed = array.array('d')
    for cls in model['classes']:
        theta = math.radians(60-phi+120*cls+offset)
        seed.extend((0., math.sin(theta), math.cos(theta)))
    path = out/f'{tag}-phi{phi}-offset{offset:g}-initial.bin'
    if path.exists():
        assert path.read_bytes() == seed.tobytes()
    else:
        path.write_bytes(seed.tobytes())
    return path


def make_job(out, tag, phi, kind, dt=5, duration=5000, perturb=0., offset=0):
    label = f'{tag}-phi{phi}-{kind}'
    initial = out/f'{tag}-phi{phi}-ideal-final.bin' if kind == 'recovery' else seed_file(out, tag, phi, offset)
    command = [str(out/'spatial'), str(out/f'{tag}.model'), 'current', str(dt), str(duration), '5',
               str(out/label), str(initial), str(perturb)]
    return dict(tag=tag, phi=phi, kind=kind, dt_fs=dt, duration_ps=duration,
                perturb=perturb, offset_deg=offset, label=label, command=command)


def refinement_check(out):
    results = []
    for phi in PHIS:
        a = read_binary(out/f'20x20-phi{phi}-perturbed.bin')
        b = read_binary(out/f'20x20-phi{phi}-refined.bin')
        assert len(a) == len(b) == 1001*300*3
        error = max(abs(x-y) for x, y in zip(a, b))
        assert error < 1e-5, (phi, error)
        results.append(dict(phi=phi, max_component_error=error))
    write_json(out/'refinement.json', results)
    return results


def sweep(out):
    assert (out/'verified-gates.json').exists()
    jobs = []
    # Pass the complete largest-cell refinement gate before expanding the sweep.
    for phi in PHIS:
        jobs.append(make_job(out, '20x20', phi, 'perturbed', perturb=.001))
        jobs.append(make_job(out, '20x20', phi, 'refined', dt=2.5, perturb=.001))
    pilot_count = len(jobs)
    for tag in ('2x2', *SIZES):
        for phi in PHIS:
            jobs.append(make_job(out, tag, phi, 'ideal'))
            if tag in SIZES and tag != '20x20':
                jobs.append(make_job(out, tag, phi, 'perturbed', perturb=.001))
    primary_count = len(jobs)
    for tag in SIZES:
        for phi in PHIS:
            jobs.append(make_job(out, tag, phi, 'recovery', duration=2000, perturb=.001))
    for phi in PHIS:
        for offset in (-5, 5):
            jobs.append(make_job(out, '20x20', phi, f'offset{offset:+d}', offset=offset))
    manifest = out/'commands.json'
    if manifest.exists():
        assert json.loads(manifest.read_text()) == jobs
    else:
        write_json(manifest, jobs)
    hashes = {str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest()
              for p in [HERE/'spatial.cpp', HERE/'prepare_periodic.m', Path(__file__)]}
    if (out/'run-source-hashes.json').exists():
        assert json.loads((out/'run-source-hashes.json').read_text()) == hashes
    else:
        write_json(out/'run-source-hashes.json', hashes)
    for start, end in ((0, pilot_count), (pilot_count, primary_count), (primary_count, len(jobs))):
        batch = []
        for job in jobs[start:end]:
            completed = out/(job['label']+'.json')
            if completed.exists():
                meta = json.loads(completed.read_text())
                assert meta['steps'] == round(job['duration_ps']*1000/job['dt_fs'])
            else:
                assert not (out/(job['label']+'.bin')).exists(), 'Partial run: preserve it and choose a fresh directory.'
                batch.append(job['command'])
        with ThreadPoolExecutor(max_workers=3) as pool:
            list(pool.map(execute, batch))
        if start == 0:
            refinement_check(out)
    analyze(out)


def analyze(out):
    jobs = json.loads((out/'commands.json').read_text())
    summary = []
    for job in jobs:
        model = model_data(out/f"{job['tag']}.model")
        n = model['n']
        classes = model['classes']
        counts = [classes.count(c) for c in range(3)]
        raw = read_binary(out/(job['label']+'.bin'))
        final = read_binary(out/(job['label']+'-final.bin'))
        meta = json.loads((out/(job['label']+'.json')).read_text())
        with (out/(job['label']+'.csv')).open() as f:
            scalars = list(csv.DictReader(f))
        assert len(raw) == len(scalars)*n*3 == (round(job['duration_ps']/5)+1)*n*3
        assert list(raw[-n*3:]) == list(final) and meta['steps'] == round(job['duration_ps']*1000/job['dt_fs'])
        assert abs(meta['time_ps']-job['duration_ps']) < 1e-9 and all(math.isfinite(v) for v in raw)
        norm_error = mean_error = max_spread = 0.
        observations = []
        for step, scalar in enumerate(scalars):
            assert abs(float(scalar['time_ps'])-step*5) < 1e-9
            state = raw[step*n*3:(step+1)*n*3]
            means = [[0., 0., 0.] for _ in range(3)]
            for i, c in enumerate(classes):
                spin = state[3*i:3*i+3]
                norm_error = max(norm_error, abs(math.sqrt(sum(x*x for x in spin))-1))
                for k in range(3):
                    means[c][k] += spin[k]/counts[c]
            mean = [sum(v[k] for v in means)/3 for k in range(3)]
            site = [sum(state[k::3])/n for k in range(3)]
            spread = max(math.sqrt(sum((state[3*i+k]-means[c][k])**2 for k in range(3)))
                         for i, c in enumerate(classes))
            max_spread = max(max_spread, spread)
            mean_error = max(mean_error, max(abs(a-b) for a, b in zip(mean, site)),
                             max(abs(site[k]-float(scalar[f'mean_{axis}'])) for k, axis in enumerate('xyz')))
            magnitude = math.hypot(mean[1], mean[2])
            phi = math.degrees(math.atan2(mean[1], mean[2])) % 360 if magnitude >= 1e-10 else math.nan
            observations.append([step*.005, phi, magnitude, spread, float(scalar['max_torque_T']),
                                 float(scalar['energy_over_mu_per_spin'])*MEV, *mean])
        assert norm_error < 1e-12 and mean_error < 1e-12
        with (out/(job['label']+'-angles.csv')).open('w', newline='') as f:
            writer = csv.writer(f)
            writer.writerow(['time_ns', 'phi_deg', 'inplane_mean', 'sublattice_spread', 'max_torque_T',
                             'energy_meV_per_spin', 'mean_x', 'mean_y', 'mean_z'])
            writer.writerows(observations)
        energy_rise = max(b[5]-a[5] for a, b in zip(observations, observations[1:]))
        direct_energy = energy([list(final[i:i+3]) for i in range(0, len(final), 3)], model)/n*MEV
        assert abs(direct_energy-observations[-1][5]) < 1e-11
        angle = observations[-1][1]
        error = abs((angle-job['phi']+180) % 360-180) if math.isfinite(angle) else 360.
        numerical_pass = energy_rise <= 1e-9 and meta['norm_error'] < 1e-12
        summary.append({k: v for k, v in job.items() if k != 'command'} | dict(
            N=n, counts=counts, final_phi_deg=angle if math.isfinite(angle) else None,
            angle_error_deg=error, energy_meV_per_spin=direct_energy, max_energy_rise_meV=energy_rise,
            max_torque_T=meta['max_torque_T'], norm_error=norm_error, mean_error=mean_error,
            preprojection_error=meta['preprojection_error'], max_spread=max_spread,
            final_spread=spread, sublattice_spins=means, numerical_pass=numerical_pass,
            stationary=meta['max_torque_T'] < 1e-7, target_retained=error < .1))
    write_json(out/'periodic-summary.json', summary)
    checks = []
    for tag in SIZES:
        rows = [r for r in summary if r['tag'] == tag]
        ideals = [read_binary(out/f'{tag}-phi{phi}-ideal-final.bin') for phi in PHIS]
        clusters = []
        for state in ideals:
            if not any(max(abs(a-b) for a, b in zip(state, other)) < 1e-6 for other in clusters):
                clusters.append(state)
        distances = {}
        for kind in ('perturbed', 'recovery', 'offset-5', 'offset+5', 'refined'):
            selected = [r for r in rows if r['kind'] == kind]
            if selected:
                distances[kind] = max(max(abs(a-b) for a, b in zip(
                    read_binary(out/(r['label']+'-final.bin')), ideals[PHIS.index(r['phi'])])) for r in selected)
        energies = [r['energy_meV_per_spin'] for r in rows if r['kind'] == 'ideal']
        spread = max(energies)-min(energies)
        passed = (len(clusters) == 6 and spread < 1e-9 and max(distances.values()) < 1e-6
                  and all(r['numerical_pass'] and r['stationary'] and r['target_retained'] for r in rows))
        checks.append(dict(size=tag, N=rows[0]['N'], distinct_states=len(clusters),
                           energy_spread_meV_per_spin=spread, recovery_errors=distances,
                           six_state_pass=passed))
    reference = [r for r in summary if r['tag'] == '2x2']
    reference_energy_spread = max(r['energy_meV_per_spin'] for r in reference)-min(r['energy_meV_per_spin'] for r in reference)
    reference_pass = (len(reference) == 6 and reference_energy_spread < 1e-9
                      and all(r['numerical_pass'] and r['target_retained'] and r['stationary'] for r in reference))
    refinement = refinement_check(out)
    report = dict(trajectories=len(summary), sizes=checks, reference_pass=reference_pass,
                  reference_energy_spread_meV_per_spin=reference_energy_spread,
                  all_pass=reference_pass and all(r['six_state_pass'] for r in checks),
                  max_torque_T=max(r['max_torque_T'] for r in summary),
                  max_norm_error=max(r['norm_error'] for r in summary),
                  max_preprojection_error=max(r['preprojection_error'] for r in summary),
                  max_refinement_error=max(r['max_component_error'] for r in refinement),
                  max_energy_rise_meV_per_spin=max(r['max_energy_rise_meV'] for r in summary),
                  max_mean_error=max(r['mean_error'] for r in summary))
    write_json(out/'analysis-source-hash.json', {str(Path(__file__).relative_to(ROOT)):
               hashlib.sha256(Path(__file__).read_bytes()).hexdigest()})
    write_json(out/'periodic-checks.json', report)
    print(json.dumps(report, indent=2), flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('output', type=Path)
    parser.add_argument('phase', choices=('gates', 'run', 'analyze'))
    args = parser.parse_args()
    out = args.output.resolve()
    assert out.is_dir() and 'onedrive' not in str(out).lower()
    {'gates': gates, 'run': sweep, 'analyze': analyze}[args.phase](out)


if __name__ == '__main__':
    main()
