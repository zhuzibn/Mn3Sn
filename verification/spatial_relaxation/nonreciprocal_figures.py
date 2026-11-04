"""Existing receiver-weighted model versus the PRB figure control; standard library.

Usage: python3 verification/spatial_relaxation/nonreciprocal_figures.py OUTPUT PHASE
Run prepare_driven(SOURCE,OUTPUT,'current') between prepare and check.
"""
import argparse
import csv
import json
import math
from pathlib import Path
import shutil
import subprocess
import sys

import driven_figures as d


def archived_check(out):
    """Recheck existing source fixtures when MATLAB is unavailable; never claim a new GPU run."""
    archive = d.ROOT.parent/'simulation_results/Mn3Sn/unequal-anisotropy-switching-20261103-02'
    compiler = json.loads((archive/'compiler.json').read_text())
    assert d.digest(d.HERE/'spatial.cpp')==compiler['source_sha256']
    assert json.loads((archive/'source-hashes.json').read_text())==json.loads((out/'source-hashes.json').read_text())
    fixtures = json.loads((archive/'gpu-gates.json').read_text())
    copied = {}
    for name in ['gpu-gates.json']+[r['name']+'.bin' for r in fixtures]:
        shutil.copy2(archive/name,out/name)
        copied[name] = d.digest(out/name)
    current_replays = []
    for size in d.SIZES:
        for trial in [0,1]:
            original = d.SOURCE/f'{size}-gpu{trial}.bin'
            tag = out/f'{size}-archived-current{trial}'
            command = [str(out/'spatial'),str(out/(size+'.model')),'current','5','.1','.005',str(tag),str(original)]
            subprocess.run(command,check=True,capture_output=True)
            error = d.distance(d.readbin(original),d.readbin(Path(str(tag)+'.bin')))
            assert error<1e-10
            current_replays.append(dict(size=size,trial=trial,max_spin_error=error,command=command,fixture_sha256=d.digest(original)))
    d.writejson(out/'validation-provenance.json',dict(fresh_matlab_gpu_run=False,
                reason='MATLAB launcher unavailable; replaying preserved, source-matched fixtures.',
                driven_archive=str(archive),copied_fixtures=copied,current_undriven_replays=current_replays,
                python=sys.version,platform=sys.platform))
    check(out)


def check(out):
    d.check(out)
    _, _, _, refclasses, refbonds, _ = d.graph(out/'2x2.model')
    results = []
    for size in d.SIZES:
        n, _, _, classes, bonds, _ = d.graph(out/(size+'.model'))
        error = 0
        for i in range(n):
            refi = refclasses.index(classes[i])
            for c in range(3):
                for k in range(2):
                    actual = math.fsum(v[k] for (p,j),v in bonds.items() if p==i and classes[j]==c)
                    expected = math.fsum(v[k] for (p,j),v in refbonds.items() if p==refi and refclasses[j]==c)
                    error = max(error, abs(actual-expected))
        assert error < 1e-14, (size, error)
        pairs = [(i,j) for i,j in bonds if i<j]
        results.append(dict(size=size, N=n, populations=[classes.count(c) for c in range(3)],
                            row_sum_error=error,
                            unequal_exchange_pairs=sum(bonds[i,j][0]!=bonds[j,i][0] for i,j in pairs),
                            nonantisymmetric_dmi_pairs=sum(bonds[i,j][1]!=-bonds[j,i][1] for i,j in pairs)))
    d.writejson(out/'coherent-reduction.json', results)
    print(json.dumps(results, indent=2), flush=True)


def spec(size, case, dt=5, perturb=0):
    return d.spec(size, case, dt, 'current', perturb)


def pilot(out):
    selected = [c for c in d.cases() if c[0] in d.PILOTS]
    d.batch(out, [spec(size,c,dt) for size in d.SIZES for c in selected for dt in [5,2.5]], 'pilot')


def pilot_check(out):
    results = []
    for size in d.SIZES:
        for case in [c for c in d.cases() if c[0] in d.PILOTS]:
            names = [spec(size,case,dt)['name'] for dt in [5,2.5]]
            assert all((out/(name+'-run.json')).exists() for name in names)
            times = [[float(v) for v in (out/(name+'-snapshot-times.csv')).read_text().splitlines()[1:]] for name in names]
            assert d.distance(*times)<1e-9
            error = d.distance(*[d.readbin(out/(name+'.bin')) for name in names])
            results.append(dict(size=size, case=case[0], max_spin_error=error, pass_=error<1e-5))
    d.writejson(out/'pilot-checks.json', dict(pass_=all(r['pass_'] for r in results), refinements=results))
    assert all(r['pass_'] for r in results), 'Refine before interpreting the sweep.'
    print('PASS pilot refinement; maximum spin error', max(r['max_spin_error'] for r in results), flush=True)


def sweep(out):
    assert json.loads((out/'pilot-checks.json').read_text())['pass_']
    specs = [spec(size,c) for size in d.SIZES for c in d.cases()]
    specs += [spec(size,c,dt,.001) for size in d.SIZES[1:] for c in d.cases()
              if c[0] in ['low_j1_h1_s210','field275','field300','highm330'] for dt in [5,2.5]]
    d.batch(out, specs, 'sweep')


def analyze(out):
    assert (out/'sweep-runs.json').exists()
    rows = []
    for s in json.loads((out/'sweep-specs.json').read_text()):
        name = s['name']
        meta = json.loads((out/(name+'-run.json')).read_text())['numerical']
        n, _, _, classes, _, _ = d.graph(out/(s['size']+'.model'))
        obs = d.observables(out/(name+'-observables.csv'))
        assert len(obs['time_ps'])==20001
        assert max(abs(t-i*.25) for i,t in enumerate(obs['time_ps']))<1e-8
        for key,values in obs.items():
            if key.startswith('interior_') and n==3:
                assert all(math.isnan(v) for v in values)
            else:
                assert all(map(math.isfinite,values)), (name,key)
        spins = d.readbin(out/(name+'.bin'))
        times = [float(v) for v in (out/(name+'-snapshot-times.csv')).read_text().splitlines()[1:]]
        assert len(spins)==len(times)*3*n and all(map(math.isfinite,spins))
        assert times[0]==0 and abs(times[-1]-5000)<1e-8
        assert d.distance(spins[-3*n:],d.readbin(out/(name+'-final.bin')))==0
        sample_error = norm_error = 0
        for frame in sorted({0,1,len(times)//2,len(times)-2,len(times)-1}):
            v = spins[frame*3*n:(frame+1)*3*n]
            idx = round(times[frame]/.25)
            for i in range(n):
                norm_error = max(norm_error,abs(math.sqrt(sum(v[3*i+k]**2 for k in range(3)))-1))
            for k,axis in enumerate('xyz'):
                site = math.fsum(v[3*i+k] for i in range(n))/n
                balanced = math.fsum(math.fsum(v[3*i+k] for i in range(n) if classes[i]==c)/classes.count(c) for c in range(3))/3
                sample_error = max(sample_error,abs(site-obs['site_'+axis][idx]),abs(balanced-obs['balanced_'+axis][idx]))
        assert sample_error<1e-12 and norm_error<1e-12
        reference = d.readbin(d.REF/(s['case']+'_m3_dt5.bin'))
        assert len(reference)==20001*10
        means = [[math.fsum(reference[i*10+1+k+3*j] for j in range(3))/3 for i in range(20001)] for k in range(3)]
        mean_error = max(d.distance(obs['balanced_'+axis],means[k]) for k,axis in enumerate('xyz'))
        phi = d.angle(obs['balanced_y'][-1],obs['balanced_z'][-1])
        rawphi = d.angle(obs['site_y'][-1],obs['site_z'][-1])
        angular_error = d.angular_error(phi,s['expected'])
        crossing = d.crossing(obs['time_ps'],obs['balanced_z'])
        refcross = d.crossing(obs['time_ps'],means[2])
        row = dict(**s,N=n,final_balanced_phi=phi,final_site_phi=rawphi,
                   endpoint_matches=angular_error is not None and angular_error<.1,
                   final_angle_error_deg=angular_error,max_mean_component_difference=mean_error,
                   balanced_crossing_ns=crossing,reference_crossing_ns=refcross,
                   site_crossing_ns=d.crossing(obs['time_ps'],obs['site_z']),
                   max_spatial_spread=max(obs['max_sublattice_spread']),final_spatial_spread=obs['max_sublattice_spread'][-1],
                   pre_pulse_balanced=[obs['balanced_'+a][3600] for a in 'xyz'],
                   pre_pulse_site=[obs['site_'+a][3600] for a in 'xyz'],
                   final_balanced=[obs['balanced_'+a][-1] for a in 'xyz'],final_site=[obs['site_'+a][-1] for a in 'xyz'],
                   preprojection_error=meta['preprojection_error'],norm_error=meta['norm_error'],snapshot_mean_error=sample_error)
        if s['size']=='2x2' and not s['perturb']:
            assert mean_error<1e-10, 'Numerical reference control failed.'
        if s['perturb'] and s['dt_fs']==5:
            fine = name.replace('-dt5-','-dt2.5-')
            fine_times = [float(v) for v in (out/(fine+'-snapshot-times.csv')).read_text().splitlines()[1:]]
            assert d.distance(times,fine_times)<1e-9
            row['refinement_error'] = d.distance(spins,d.readbin(out/(fine+'.bin')))
            row['final_recovery_error'] = d.distance(spins[-3*n:],d.readbin(out/(name.replace('-p0.001','-p0')+'-final.bin')))
        rows.append(row)
    groups = []
    for size in d.SIZES:
        group = [r for r in rows if r['size']==size and r['perturb']==0]
        assert len(group)==48
        groups.append(dict(size=size,N=group[0]['N'],
                           low_field_matches=sum(r['endpoint_matches'] for r in group if r['case'].startswith('low_')),
                           high_field_matches=sum(r['endpoint_matches'] for r in group if r['case'].startswith('high')),
                           field_sweep=[{k:r[k] for k in ['field_Oe','final_balanced_phi','final_site_phi','endpoint_matches']}
                                        for r in group if r['case'].startswith('field')],
                           fig4=next(r for r in group if r['case']=='low_j1_h1_s210'),
                           max_mean_component_difference=max(r['max_mean_component_difference'] for r in group),
                           max_spatial_spread=max(r['max_spatial_spread'] for r in group)))
    refinements = [r for r in rows if 'refinement_error' in r]
    report = dict(primary_trajectories=192,groups=groups,perturbations=refinements,
                  max_norm_error=max(r['norm_error'] for r in rows),
                  max_preprojection_error=max(r['preprojection_error'] for r in rows),
                  max_snapshot_mean_error=max(r['snapshot_mean_error'] for r in rows))
    d.writejson(out/'figure-results.json',rows)
    d.writejson(out/'figure-summary.json',report)
    columns = ['size','case','dt_fs','perturb','final_balanced_phi','final_site_phi','endpoint_matches',
               'max_mean_component_difference','balanced_crossing_ns','site_crossing_ns','max_spatial_spread','final_spatial_spread']
    with (out/'figure-results.csv').open('w',newline='') as f:
        w = csv.DictWriter(f,fieldnames=columns,extrasaction='ignore');w.writeheader();w.writerows(rows)
    assert all(r['refinement_error']<1e-5 for r in refinements), 'Perturbed runs require refinement.'
    print(json.dumps(dict(groups=[{k:g[k] for k in ['size','low_field_matches','high_field_matches','max_mean_component_difference','max_spatial_spread']} for g in groups],
                         max_perturbed_refinement=max(r['refinement_error'] for r in refinements),
                         max_final_recovery=max(r['final_recovery_error'] for r in refinements)),indent=2),flush=True)


if __name__=='__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('output',type=Path)
    p.add_argument('phase',choices=['prepare','check','archived-check','pilot','pilot-check','sweep','analyze','preserve'])
    args = p.parse_args();out = args.output.resolve()
    if args.phase!='prepare':
        d.verify_sources(out)
    {'prepare':d.prepare,'check':check,'archived-check':archived_check,'pilot':pilot,'pilot-check':pilot_check,
     'sweep':sweep,'analyze':analyze,'preserve':d.preserve}[args.phase](out)
