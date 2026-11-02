"""Figures 2, 4(a), 5 after production GPU/replay gates. Standard library only.

Prepare/static/probes/initial gates: see docs/frozen-field-repair.md.
python3 verification/coupled_figures.py sweep|analyze|replays OUTPUT
"""
import csv
import json
import math
from pathlib import Path
import sys

from three_spin import batch, compare, physical_cases, readbin, run_case


def cases():
    for js in [1, -1]:
        for hs in [1, -1]:
            for seed in [30, 90, 150, 210, 270, 330]:
                expected = {(1, 1): 30, (1, -1): -150, (-1, -1): -30, (-1, 1): 150}[js, hs]
                yield f'low_j{js}_h{hs}_s{seed}', seed, hs*100, js*5e10, expected
    for name, seed, field, current in physical_cases():
        yield name, seed, field, current, 30 if field <= 275 else 90


def sweep(out):
    for name in ['production-gates', 'production-initial', 'dynamics']:
        assert json.loads((out/f'{name}.json').read_text())['pass']
    specs = []
    for label, seed, field, current, expected in cases():
        for mode, steps in [(2, [5, 2.5, 1.25]), (3, [5])]:
            for dt in steps:
                specs.append(dict(name=f'{label}_m{mode}_dt{dt}', case=label, mode=mode,
                                  dt_fs=dt, seed_deg=seed, field_Oe=field, current=current,
                                  duration_ns=5, on_ns=1, width_ns=3, expected=expected))
    results = batch(out, specs)
    (out/'figures-runs.json').write_text(json.dumps(results, indent=2)+'\n')


def crossing(data):
    for n in range(10, len(data), 10):
        z = sum(data[n+3+3*i] for i in range(3))/3
        prev = sum(data[n-10+3+3*i] for i in range(3))/3
        if data[n] > 1e-9 and prev < 0 <= z:
            return (data[n-10]-prev*(data[n]-data[n-10])/(z-prev))*1e9
    raise AssertionError('No zero crossing')


def analyze(out):
    runs = json.loads((out/'figures-runs.json').read_text())
    index = {(r['case'], r['mode'], r['dt_fs']): r for r in runs}
    maxnorm = 0
    for r in runs:
        data = readbin(out/(r['name']+'.bin'))
        assert len(data) == 20001*10
        assert r['updates'] == round(5e6/r['dt_fs'])
        assert r['max_postprojection_norm_error'] < 1e-12
        assert abs((r['final_angle_deg']-r['expected']+180) % 360-180) < .1, r['name']
        for n in range(0, len(data), 10):
            assert all(math.isfinite(v) for v in data[n:n+10])
            assert abs(data[n]-(n//10)*.25e-12) < 1e-22
            for i in range(3):
                maxnorm = max(maxnorm, abs(math.sqrt(sum(data[n+1+3*i+q]**2 for q in range(3)))-1))
    assert maxnorm < 1e-12
    comparisons = []
    for label, *_ in cases():
        a, b, c, ref = [index[label, mode, dt] for mode, dt in [(2, 5), (2, 2.5), (2, 1.25), (3, 5)]]
        e1 = compare(out, a['name'], b['name'])['max_mean_component_difference']
        e2 = compare(out, b['name'], c['name'])['max_mean_component_difference']
        match = compare(out, a['name'], ref['name'])
        # The planned acceptance observable is the mean spin. Also retain the
        # stricter auxiliary per-spin result, including its 275 Oe exception.
        assert max(e1, e2) < 1e-5 and match['max_mean_component_difference'] < 1e-10, (label, e1, e2, match)
        comparisons.append(dict(case=label, coarse_error=e1, fine_error=e2,
                                strict_spin_match_pass=match['max_spin_component_difference'] < 1e-9, **match))
    fig4 = []
    for dt in [5, 2.5, 1.25]:
        r = index['low_j1_h1_s210', 2, dt]
        data = readbin(out/(r['name']+'.bin'))
        def mean_at(k):
            return [sum(data[10*k+1+q+3*i] for i in range(3))/3 for q in range(3)]
        fig4.append(dict(dt_fs=dt, crossing_ns=crossing(data), before=mean_at(3600),
                         at_removal=mean_at(16000), final=mean_at(20000)))
    assert max(r['crossing_ns'] for r in fig4)-min(r['crossing_ns'] for r in fig4) < .001
    result = dict(pass_=True, runs=len(runs), norm_error=maxnorm,
                  max_preprojection_error=max(r['max_preprojection_norm_error'] for r in runs),
                  comparisons=comparisons, fig4=fig4)
    result['pass'] = result.pop('pass_')
    (out/'figures-analysis.json').write_text(json.dumps(result, indent=2)+'\n')
    keys = ['name','case','mode','dt_fs','seed_deg','field_Oe','current','expected','final_angle_deg',
            'updates','max_preprojection_norm_error','max_postprojection_norm_error']
    with (out/'figures-summary.csv').open('w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=keys, extrasaction='ignore')
        writer.writeheader(); writer.writerows(runs)
    print(json.dumps({k:v for k,v in result.items() if k != 'comparisons'}, indent=2))
    print('Maximum halving errors:', max(r['coarse_error'] for r in comparisons), max(r['fine_error'] for r in comparisons))
    print('Maximum reference spin discrepancy:', max(r['max_spin_component_difference'] for r in comparisons))


def replays(out):
    runs = json.loads((out/'figures-runs.json').read_text()); result = []
    labels = [f'low_j{js}_h{hs}_s210' for js in [1,-1] for hs in [1,-1]]
    labels += ['field275', 'field300', 'highp30', 'highm330']
    for label in labels:
        r = next(r for r in runs if r['case'] == label and r['mode'] == 2 and r['dt_fs'] == 5)
        data = readbin(out/(r['name']+'.bin')); start = 4080  # 1.02 ns
        seed = data[start*10+1:start*10+10]; name = f'production_driven_{label}'
        seedfile = name+'.csv'
        with (out/seedfile).open('w', newline='') as f:
            csv.writer(f).writerows([seed[i*3:(i+1)*3] for i in range(3)])
        spec = dict(name=name, mode=2, dt_fs=5, seed_deg=r['seed_deg'], field_Oe=r['field_Oe'],
                    current=r['current'], duration_ns=.005, on_ns=0, width_ns=.005, seed_file=seedfile)
        run_case(out, spec, sample_ps=.005)
        actual = readbin(out/(name+'.bin'))
        error = max(abs(actual[k*50*10+q]-data[(start+k)*10+q]) for k in range(21) for q in range(1,10))
        assert error < 1e-14
        spec['restart_error'] = error; result.append(spec)
    (out/'production-replays.json').write_text(json.dumps(result, indent=2)+'\n')
    print('PASS 8 production replay seeds/segments')


if __name__ == '__main__':
    {'sweep': sweep, 'analyze': analyze, 'replays': replays}[sys.argv[1]](Path(sys.argv[2]).resolve())
