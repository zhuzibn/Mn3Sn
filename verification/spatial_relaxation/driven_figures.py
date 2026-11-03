"""Driven candidate comparison, using standard-library Python and source-checked C++.

Usage: python3 verification/spatial_relaxation/driven_figures.py OUTPUT PHASE
PHASE: prepare, check, pilot, pilot-check, sweep, transition, analyze, preserve.
Run prepare_driven.m between prepare and check. OUTPUT must be a fresh folder outside OneDrive.
"""
import argparse
import array
from concurrent.futures import ThreadPoolExecutor, as_completed
import csv
import hashlib
import itertools
import json
import math
from pathlib import Path
import random
import shutil
import subprocess
import sys
import time

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent.parent
SOURCE = ROOT / '.sisyphus/spatial-relaxation-20261102-01'
UND = ROOT / '.sisyphus/anisotropy-trial-20261102-01'
REF = ROOT / '.sisyphus/frozen-field-20261102-01'
SIZES = ['2x2', '5x5', '13x14', '21x22']
B = 7.448364280185585e-5
sys.path.insert(0, str(HERE.parent))
from coupled_figures import cases


def digest(path):
    h = hashlib.sha256()
    with path.open('rb') as f:
        for block in iter(lambda: f.read(1024*1024), b''):
            h.update(block)
    return h.hexdigest()


def writejson(path, data):
    path.write_text(json.dumps(data, indent=2, allow_nan=False) + '\n')


def readbin(path):
    a = array.array('d')
    a.frombytes(path.read_bytes())
    assert sys.byteorder == 'little'
    return a


def distance(a, b):
    assert len(a) == len(b)
    return max(abs(x-y) for x, y in zip(a, b))


def graph(path):
    lines = path.read_text().splitlines()
    header = list(map(float, lines[0].split()))
    n = int(header[0])
    axes = [list(map(float, s.split()))[:3] for s in lines[1:n+1]]
    classes = [1 if e[1] < -.5 else (2 if e[2] < 0 else 0) for e in axes]
    bonds = {(int(i), int(j)): (a, c) for i, j, a, c in
             (map(float, s.split()) for s in lines[n+2:])}
    degree = [sum(p == i and a != 0 for (p, _), (a, _) in bonds.items()) for i in range(n)]
    return n, header, axes, classes, bonds, degree


def dot(a, b):
    return sum(x*y for x, y in zip(a, b))


def cross(a, b):
    return [a[1]*b[2]-a[2]*b[1], a[2]*b[0]-a[0]*b[2], a[0]*b[1]-a[1]*b[0]]


def verify_sources(out):
    for name, expected in json.loads((out/'source-hashes.json').read_text()).items():
        assert digest(ROOT/name) == expected, name
    compiler = json.loads((out/'compiler.json').read_text())
    assert digest(HERE/'spatial.cpp') == compiler['source_sha256'], 'Reprepare after kernel changes.'
    assert digest(out/'spatial') == compiler['binary_sha256']


def prepare(out):
    assert not any(part.lower() == 'onedrive' for part in out.parts), 'Save simulation results outside OneDrive.'
    assert not out.exists(), 'Use a fresh output directory.'
    out.mkdir(parents=True)
    shutil.copy2(SOURCE/'source-hashes.json', out/'source-hashes.json')
    command = ['g++', '-O3', '-std=c++17', '-Wall', '-Wextra', str(HERE/'spatial.cpp'), '-o', str(out/'spatial')]
    subprocess.run(command, check=True)
    writejson(out/'compiler.json', dict(command=command, source_sha256=digest(HERE/'spatial.cpp'),
              binary_sha256=digest(out/'spatial'), compiler=subprocess.check_output(['g++', '--version'], text=True)))
    protected = {str(ROOT/p): digest(ROOT/p) if (ROOT/p).exists() else None
                 for p in ['finalnew.mat', 'atomtypee.mat']}
    writejson(out/'protected.json', protected)
    for size in SIZES:
        shutil.copy2(SOURCE/f'{size}.model', out/f'{size}.model')
        n, _, _, classes, _, _ = graph(out/f'{size}.model')
        for seed in range(30, 360, 60):
            state = array.array('d')
            for c in classes:
                theta = math.radians(seed+120*c)
                state.extend([0, math.sin(theta), math.cos(theta)])
            (out/f'{size}-seed{seed}.bin').write_bytes(state.tobytes())
    writejson(out/'provenance.json', dict(revision=subprocess.check_output(['git', 'rev-parse', 'HEAD'], text=True).strip(),
              branch=subprocess.check_output(['git', 'branch', '--show-current'], text=True).strip(),
              source=str(SOURCE), undriven=str(UND), reference=str(REF), python=sys.version,
              models={s: digest(out/f'{s}.model') for s in SIZES}))
    verify_sources(out)
    print(f'Prepared {out}', flush=True)


def check(out):
    assert not (out/'driven-gates.json').exists()
    gates = []
    for size in SIZES:
        n, header, axes, classes, bonds, degree = graph(out/f'{size}.model')
        J, D, K = header[3:6]
        candidate = {(i, j): (J/2 if a else 0, math.copysign(D/2, c) if c else 0)
                     for (i, j), (a, c) in bonds.items()}
        assert set(degree) <= {2, 4}
        for (i, j), (a, c) in candidate.items():
            assert candidate[j, i] == (a, -c)
        neighbors = [{j for p, j in bonds if p == i} for i in range(n)]
        triangles = [(i, j, k) for i in range(n) for j, k in itertools.combinations(sorted(neighbors[i]), 2)
                     if i < j < k and k in neighbors[j]]
        assert all(sum(i in t and j in t for t in triangles) == 1 for i, j in bonds)
        assert all(sorted(classes[i] for i in t) == [0, 1, 2] for t in triangles)
        assert all(sum(i in t for t in triangles) == degree[i]/2 for i in range(n))
        weights = [sum(degree[i]/4 for i in range(n) if classes[i] == c) for c in range(3)]
        assert len(set(weights)) == 1
        rng = random.Random(20261103)
        state = [[rng.uniform(-1, 1) for _ in range(3)] for _ in range(n)]

        def energy(s):
            return math.fsum([-.5*K*degree[i]/4*dot(s[i], axes[i])**2-.1*s[i][1] for i in range(n)] +
                             [a*dot(s[i], s[j])-c*cross(s[i], s[j])[0]
                              for (i, j), (a, c) in candidate.items() if i < j])

        h = [[K*degree[i]/4*dot(state[i], axes[i])*v for v in axes[i]] for i in range(n)]
        for i in range(n):
            h[i][1] += .1
        for (i, j), (a, c) in candidate.items():
            for k in range(3):
                h[i][k] -= a*state[j][k]
            h[i][1] += c*state[j][2]
            h[i][2] -= c*state[j][1]
        error = 0
        for i in sorted({0, n//3, n//2, n-1}):
            for k in range(3):
                old = state[i][k]
                state[i][k] = old+1e-5
                ep = energy(state)
                state[i][k] = old-1e-5
                em = energy(state)
                state[i][k] = old
                error = max(error, abs((ep-em)/2e-5+h[i][k]))
        assert error < 1e-6
        for mode in ['current', 'reciprocal', 'triangle']:
            for executable, suffix in [(UND/'spatial', 'old'), (out/'spatial', 'new')]:
                subprocess.run([str(executable), str(out/f'{size}.model'), mode, '5', '.1', '.005',
                                str(out/f'{size}-{mode}-{suffix}')], check=True, capture_output=True)
            for ext in ['.bin', '.csv']:
                assert (out/f'{size}-{mode}-old{ext}').read_bytes() == (out/f'{size}-{mode}-new{ext}').read_bytes()
        gates.append(dict(size=size, N=n, populations=[classes.count(c) for c in range(3)],
                          triangles=len(triangles), energy_gradient_error_T=error, zero_drive_regression=True))
    replays = []
    for r in json.loads((out/'gpu-gates.json').read_text()):
        tag = out/(r['name']+'-replay')
        cmd = [str(out/'spatial'), str(out/(r['size']+'.model')), r['mode'], '5', '.1', '.005', str(tag),
               str(out/(r['name']+'.bin')), '0', str(r['field_y_T']), str(r['sot_T']), '.025', '.075', '.005']
        subprocess.run(cmd, check=True, capture_output=True)
        expected, actual = readbin(out/(r['name']+'.bin')), readbin(Path(str(tag)+'.bin'))
        assert len(expected) == 21*3*r['N']
        error = distance(expected, actual)
        assert error < 1e-10
        meta = json.loads(Path(str(tag)+'.json').read_text())
        assert meta['steps'] == 20 and meta['active_updates'] == 10
        replays.append(dict(name=r['name'], max_spin_error=error, command=cmd))
    writejson(out/'driven-gates.json', dict(pass_=True, static=gates, replays=replays))
    print(f'PASS static, undriven regression and {len(replays)} MATLAB/GPU driven replays', flush=True)


def spec(size, case, dt=5, mode='triangle', perturb=0):
    label, seed, field, current, expected = case
    name = f'{size}-{mode}-{label}-dt{dt}-p{perturb}'
    return dict(name=name, size=size, mode=mode, case=label, seed=seed, field_Oe=field,
                current=current, expected=expected, dt_fs=dt, perturb=perturb)


def execute(out, s):
    done = out/(s['name']+'-run.json')
    if done.exists():
        old = json.loads(done.read_text())
        assert old['spec'] == s
        return old
    prefix = out/s['name']
    assert not Path(str(prefix)+'.bin').exists(), f'Incomplete run: {prefix}'
    cmd = [str(out/'spatial'), str(out/(s['size']+'.model')), s['mode'], str(s['dt_fs']), '5000', '5', str(prefix),
           str(out/f"{s['size']}-seed{s['seed']}.bin"), str(s['perturb']), str(s['field_Oe']*1e-4),
           str(B*s['current']/5e10), '1000', '4000', '.25']
    start = time.monotonic()
    subprocess.run(cmd, check=True, capture_output=True, text=True)
    meta = json.loads(Path(str(prefix)+'.json').read_text())
    assert meta['steps'] == round(5e6/s['dt_fs']) and meta['active_updates'] == round(3e6/s['dt_fs'])
    assert meta['norm_error'] < 1e-12
    row = dict(spec=s, command=cmd, elapsed_seconds=time.monotonic()-start, numerical=meta,
               binary_sha256=digest(out/'spatial'))
    writejson(done, row)
    print(f"Completed {s['name']} ({row['elapsed_seconds']:.1f} s)", flush=True)
    return row


def batch(out, specs, manifest):
    assert json.loads((out/'driven-gates.json').read_text())['pass_']
    writejson(out/(manifest+'-specs.json'), specs)
    with ThreadPoolExecutor(max_workers=3) as pool:
        jobs = [pool.submit(execute, out, s) for s in specs]
        results = [f.result() for f in as_completed(jobs)]
    writejson(out/(manifest+'-runs.json'), results)


PILOTS = ['low_j1_h1_s270', 'low_j1_h-1_s270', 'low_j-1_h-1_s270', 'low_j1_h1_s210',
          'field275', 'field300', 'highp30', 'highm330']


def pilot(out):
    selected = [c for c in cases() if c[0] in PILOTS]
    assert len(selected) == len(PILOTS)
    specs = [spec(size, c, dt) for size in SIZES for c in selected for dt in [5, 2.5]]
    specs += [spec('2x2', c, 5, 'current') for c in selected]
    batch(out, specs, 'pilot')


def pilot_check(out):
    results = []
    for size in SIZES:
        for case in [c for c in cases() if c[0] in PILOTS]:
            names = [spec(size, case, dt)['name'] for dt in [5, 2.5]]
            assert all((out/(n+'-run.json')).exists() for n in names)
            times = [(out/(n+'-snapshot-times.csv')).read_text() for n in names]
            # Floating print rounding may differ; numerical sample grids must agree.
            ta, tb = [[float(x) for x in t.splitlines()[1:]] for t in times]
            assert distance(ta, tb) < 1e-9
            error = distance(*[readbin(out/(n+'.bin')) for n in names])
            results.append(dict(size=size, case=case[0], max_spin_error=error, pass_=error<1e-5))
    writejson(out/'pilot-checks.json', dict(pass_=all(r['pass_'] for r in results), refinements=results))
    print(json.dumps(results, indent=2), flush=True)
    assert all(r['pass_'] for r in results), 'Diagnose/refine numerical failures before the sweep.'


def sweep(out):
    assert json.loads((out/'pilot-checks.json').read_text())['pass_']
    specs = [spec(size, c) for size in SIZES for c in cases()]
    # Independent full-strength control, all 48 cases; no candidate rescaling.
    specs += [spec('2x2', c, 5, 'current') for c in cases()]
    # Perturb selected driven starts on each spatial patch, with refinement.
    specs += [spec(size, c, dt, perturb=.001) for size in SIZES[1:]
              for c in cases() if c[0] in ['low_j1_h1_s210', 'field275', 'field300', 'highm330']
              for dt in [5, 2.5]]
    batch(out, specs, 'sweep')


def observables(path):
    with path.open() as f:
        reader = csv.reader(f)
        header = next(reader)
        columns = {k:array.array('d') for k in header}
        for row in reader:
            for k, value in zip(header, row):
                columns[k].append(float(value))
    return columns


def angle(y, z):
    return math.degrees(math.atan2(y, z)) % 360 if math.hypot(y, z)>1e-10 else None


def angular_error(a, b):
    return abs((a-b+180) % 360-180) if a is not None and b is not None else None


def transition(out):
    # The half-strength control shifts its sampled transition to 400--450 Oe.
    # Verify both bounds and the midpoint; do not retune any physical input.
    assert json.loads((out/'pilot-checks.json').read_text())['pass_']
    results = []
    for field in [400, 425, 450]:
        case = (f'field{field}', 210, field, 5e10, 90)
        names = []
        for dt in [5, 2.5]:
            s = spec('2x2', case, dt)
            execute(out, s);names.append(s['name'])
        error = distance(*[readbin(out/(name+'.bin')) for name in names])
        assert error < 1e-5, (field, error)
        obs = observables(out/(names[0]+'-observables.csv'))
        results.append(dict(field_Oe=field, final_phi_deg=angle(obs['balanced_y'][-1],obs['balanced_z'][-1]),
                            max_spin_refinement_error=error))
    writejson(out/'triangle-transition.json', results)
    print(json.dumps(results, indent=2), flush=True)


def crossing(t, z):
    for i in range(1, len(t)):
        if t[i]>1000 and z[i-1]<0<=z[i]:
            return (t[i-1]-z[i-1]*(t[i]-t[i-1])/(z[i]-z[i-1]))/1000
    return None


def analyze(out):
    assert (out/'sweep-runs.json').exists()
    all_specs = json.loads((out/'sweep-specs.json').read_text())
    rows, refinements, robustness = [], [], []
    # Ref bins store independent full-strength triangle a,b,c spins after time.
    for s in all_specs:
        name = s['name']
        meta = json.loads((out/(name+'-run.json')).read_text())['numerical']
        n, _, _, classes, _, _ = graph(out/(s['size']+'.model'))
        obs = observables(out/(name+'-observables.csv'))
        assert len(obs['time_ps']) == 20001
        assert max(abs(t-i*.25) for i, t in enumerate(obs['time_ps'])) < 1e-8
        for k, values in obs.items():
            if k.startswith('interior_') and n==3:
                assert all(math.isnan(x) for x in values)
            else:
                assert all(math.isfinite(x) for x in values), (name, k)
        snapshots = readbin(out/(name+'.bin'))
        snapshot_times = [float(v) for v in (out/(name+'-snapshot-times.csv')).read_text().splitlines()[1:]]
        assert len(snapshots) == len(snapshot_times)*3*n and all(map(math.isfinite, snapshots))
        assert snapshot_times[0] == 0 and abs(snapshot_times[-1]-5000)<1e-8
        assert distance(snapshots[-3*n:], readbin(out/(name+'-final.bin'))) == 0
        # Independently reconcile dense means with actual full-spin snapshots.
        sample_error = norm_error = 0
        for frame in sorted({0, 1, len(snapshot_times)//2, len(snapshot_times)-2, len(snapshot_times)-1}):
            v = snapshots[frame*3*n:(frame+1)*3*n]
            sample_index = round(snapshot_times[frame]/.25)
            for i in range(n):
                norm_error = max(norm_error, abs(math.sqrt(sum(v[3*i+k]**2 for k in range(3)))-1))
            for k, axis in enumerate('xyz'):
                site = math.fsum(v[3*i+k] for i in range(n))/n
                balanced = sum(math.fsum(v[3*i+k] for i in range(n) if classes[i]==c)/classes.count(c) for c in range(3))/3
                sample_error = max(sample_error, abs(site-obs['site_'+axis][sample_index]),
                                   abs(balanced-obs['balanced_'+axis][sample_index]))
        assert sample_error < 1e-12 and norm_error < 1e-12
        phi = angle(obs['balanced_y'][-1], obs['balanced_z'][-1])
        rawphi = angle(obs['site_y'][-1], obs['site_z'][-1])
        error = angular_error(phi, s['expected'])
        reference = readbin(REF/(s['case']+'_m3_dt5.bin'))
        assert len(reference)==20001*10
        ref_mean = [[math.fsum(reference[i*10+1+k+3*j] for j in range(3))/3 for i in range(20001)] for k in range(3)]
        mean_error = max(distance(obs['balanced_'+axis], ref_mean[k]) for k, axis in enumerate('xyz'))
        angle_errors = [angular_error(angle(obs['balanced_y'][i],obs['balanced_z'][i]), angle(ref_mean[1][i],ref_mean[2][i]))
                        for i in range(200,20001)]
        defined_errors = [e for e in angle_errors if e is not None]
        cross_bal = crossing(obs['time_ps'], obs['balanced_z'])
        cross_raw = crossing(obs['time_ps'], obs['site_z'])
        cross_ref = crossing(obs['time_ps'], ref_mean[2])
        # Retain the saved full-spin states bracketing the first mean reversal.
        crossing_snapshots = None
        if cross_bal is not None:
            before = max(i for i,t in enumerate(snapshot_times) if t<=cross_bal*1000)
            after = min(i for i,t in enumerate(snapshot_times) if t>=cross_bal*1000)
            crossing_snapshots = dict(before_index=before, after_index=after,
                                      before_ps=snapshot_times[before], after_ps=snapshot_times[after])
        row = dict(**s, N=n, populations=[classes.count(c) for c in range(3)], final_balanced_phi=phi,
                   final_site_phi=rawphi, final_angle_error_deg=error, endpoint_matches=error is not None and error<.1,
                   max_mean_component_difference=mean_error, max_defined_angle_difference_deg=max(defined_errors) if defined_errors else None,
                   balanced_crossing_ns=cross_bal, site_crossing_ns=cross_raw, reference_crossing_ns=cross_ref,
                   crossing_difference_ns=cross_bal-cross_ref if cross_bal is not None and cross_ref is not None else None,
                   crossing_snapshots=crossing_snapshots, max_spatial_spread=max(obs['max_sublattice_spread']),
                   final_spatial_spread=obs['max_sublattice_spread'][-1], final_torque_T=meta['max_torque_T'],
                   final_balanced=[obs['balanced_'+axis][-1] for axis in 'xyz'],
                   final_site=[obs['site_'+axis][-1] for axis in 'xyz'],
                   preprojection_error=meta['preprojection_error'], norm_error=meta['norm_error'], snapshot_mean_error=sample_error)
        rows.append(row)
        if s['mode']=='current':
            assert mean_error < 1e-10, (name, mean_error)
        if s['perturb'] and s['dt_fs']==5:
            finer = {**s, 'dt_fs':2.5, 'name':name.replace('-dt5-', '-dt2.5-')}
            err = distance(snapshots, readbin(out/(finer['name']+'.bin')))
            refinements.append(dict(size=s['size'], case=s['case'], max_spin_error=err, pass_=err<1e-5))
            unperturbed = name.replace('-p0.001', '-p0')
            base = readbin(out/(unperturbed+'-final.bin'))
            robustness.append(dict(size=s['size'], case=s['case'], final_spin_difference=distance(base,snapshots[-3*n:])))
        print(f'Analyzed {name}: balanced={phi}, site={rawphi}, endpoint={row["endpoint_matches"]}', flush=True)
    groups = []
    for size, mode in [(s,'triangle') for s in SIZES]+[('2x2','current')]:
        group = [r for r in rows if r['size']==size and r['mode']==mode and r['perturb']==0]
        assert len(group)==48
        low = [r for r in group if r['case'].startswith('low_')]
        high = [r for r in group if r['case'].startswith('high')]
        sweep_rows = sorted([r for r in group if r['case'].startswith('field')], key=lambda r:r['field_Oe'])
        fig4 = next(r for r in low if r['case']=='low_j1_h1_s210')
        groups.append(dict(size=size, mode=mode, N=group[0]['N'], low_field_matches=sum(r['endpoint_matches'] for r in low),
                           high_field_matches=sum(r['endpoint_matches'] for r in high),
                           field_sweep=[{k:r[k] for k in ['field_Oe','final_balanced_phi','final_site_phi','endpoint_matches']} for r in sweep_rows],
                           fig4={k:fig4[k] for k in ['balanced_crossing_ns','site_crossing_ns','crossing_difference_ns','final_balanced','final_site','max_spatial_spread']},
                           max_mean_component_difference=max(r['max_mean_component_difference'] for r in group)))
    writejson(out/'figure-results.json', rows)
    report = dict(physical_cases=48, primary_trajectories=240, perturbation_refinements=refinements,
                  perturbation_recovery=robustness, groups=groups,
                  max_norm_error=max(r['norm_error'] for r in rows), max_preprojection_error=max(r['preprojection_error'] for r in rows),
                  max_snapshot_mean_error=max(r['snapshot_mean_error'] for r in rows))
    writejson(out/'figure-summary.json', report)
    columns = ['size','mode','case','seed','field_Oe','current','dt_fs','perturb','final_balanced_phi','final_site_phi',
               'endpoint_matches','max_mean_component_difference','balanced_crossing_ns','site_crossing_ns',
               'crossing_difference_ns','max_spatial_spread','final_torque_T']
    with (out/'figure-results.csv').open('w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=columns, extrasaction='ignore')
        writer.writeheader();writer.writerows(rows)
    assert all(r['pass_'] for r in refinements), 'Perturbed trajectories require refinement before interpretation.'
    print(json.dumps(report, indent=2), flush=True)


def preserve(out):
    verify_sources(out)
    for name, expected in json.loads((out/'protected.json').read_text()).items():
        p = Path(name)
        assert (digest(p) if p.exists() else None) == expected, name
    writejson(out/'preservation.json', dict(pass_=True, harness={p.name:digest(p) for p in HERE.iterdir() if p.is_file()}))
    print('PASS protected sources/results', flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('output', type=Path)
    parser.add_argument('phase', choices=['prepare', 'check', 'pilot', 'pilot-check', 'sweep', 'transition', 'analyze', 'preserve'])
    args = parser.parse_args()
    out = args.output.resolve()
    if args.phase != 'prepare':
        verify_sources(out)
    {'prepare':prepare, 'check':check, 'pilot':pilot, 'pilot-check':pilot_check,
     'sweep':sweep, 'transition':transition, 'analyze':analyze, 'preserve':preserve}[args.phase](out)
