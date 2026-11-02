"""Isolated three-spin comparison. Python standard library only; no packages.

Usage: python3 verification/three_spin.py prepare OUTPUT
       python3 verification/three_spin.py preserve OUTPUT
OUTPUT must be a fresh directory under this repository's .sisyphus directory.
Run stages in the order documented in docs/three-spin-comparison.md.
"""
import hashlib
import json
import array
import csv
from concurrent.futures import ThreadPoolExecutor
import math
from pathlib import Path
import shutil
import subprocess
import sys

REPO = Path(__file__).resolve().parents[1]
REFERENCE = REPO.parent / 'Mn3Sn_3atoms'


def digest(path):
    h = hashlib.sha256()
    with path.open('rb') as f:
        for block in iter(lambda: f.read(1024*1024), b''):
            h.update(block)
    return h.hexdigest()


def prepare(out):
    assert out.parent == REPO / '.sisyphus' and not out.exists()
    out.mkdir(parents=True)
    protected = {}
    for label, root in [('native', REPO), ('reference', REFERENCE)]:
        folder = out / label
        folder.mkdir()
        for p in sorted(root.glob('*.m')):
            shutil.copy2(p, folder / p.name)
            protected[str(p)] = digest(p)
        for name in ['finalnew.mat', 'atomtypee.mat', 'final.mat', 'result.txt']:
            p = root / name
            protected[str(p)] = digest(p) if p.exists() else None
    (out / 'protected.json').write_text(json.dumps(protected, indent=2) + '\n')
    # Execute only source initialization excerpts, never either main.m.
    s = (out / 'reference/main.m').read_text().split('\nrk4_4llg();')[0]
    s = s[:s.index('%% calc')]
    s = s.replace('clear ;clc;close all;', '').replace('\ntic\n', '\n')
    (out / 'reference/initialize.m').write_text(s)
    print(out)


def preserve(out):
    manifest = json.loads((out / 'protected.json').read_text())
    for name, expected in manifest.items():
        p = Path(name)
        actual = digest(p) if p.exists() else None
        assert actual == expected, f'Protected source/result changed: {name}'
    (out / 'preservation.json').write_text(json.dumps({'pass': True, 'paths': len(manifest)}, indent=2)+'\n')
    sources={p.name:digest(p) for p in (REPO/'verification').iterdir() if p.is_file()}
    if (out/'compiler.json').exists():
        compiler=json.loads((out/'compiler.json').read_text())
        assert sources['three_spin_replay.cpp']==compiler['source_sha256']
    # Keep the exact harness alongside source snapshots and raw outputs.
    (out/'harness').mkdir(exist_ok=True)
    for name in sources:
        shutil.copy2(REPO/'verification'/name,out/'harness'/name)
    provenance={'harness_sha256':sources,'python':sys.version,
        'native_revision':subprocess.check_output(['git','-C',str(REPO),'rev-parse','HEAD'],text=True).strip(),
        'reference_revision':subprocess.check_output(['git','-C',str(REFERENCE),'rev-parse','HEAD'],text=True).strip(),
        'binary_sha256':{p.name:digest(p) for p in sorted(out.glob('*.bin'))}}
    (out/'run-manifest.json').write_text(json.dumps(provenance,indent=2)+'\n')
    print(f'PASS protected sources/results: {len(manifest)} paths')


def run_case(out, spec, sample_ps=.25):
    mode=spec['mode']
    args=[str(out/'replay'), str(out/(spec['name']+'.bin')),
          str(out/('native.params' if mode in (0, 2) else 'reference.params')), str(mode),
          *[str(spec[k]) for k in ['dt_fs','seed_deg','field_Oe','current','duration_ns','on_ns','width_ns']],
          str(sample_ps), 'physical']
    if 'seed_file' in spec:
        args.append(str(out/spec['seed_file']))
    result=subprocess.run(args, check=True, capture_output=True, text=True)
    metadata={**spec, **json.loads(result.stdout), 'sample_ps':sample_ps, 'command':args}
    (out/(spec['name']+'.json')).write_text(json.dumps(metadata,indent=2)+'\n')
    return metadata


def probes(out):
    static=json.loads((out/'static.json').read_text())
    assert static['pass']
    for name in ['native','reference']:
        p=static[name]
        values=[p['J']/p['mu'],p['D']/p['mu'],p['Kfield']/p['mu'],p['gamma'],p['alpha'],p['b']]
        (out/(name+'.params')).write_text(' '.join(format(v,'.17g') for v in values)+'\n')
    command=['g++','-O3','-std=c++17','-Wall','-Wextra','-ffp-contract=off',
             str(REPO/'verification/three_spin_replay.cpp'),'-o',str(out/'replay')]
    subprocess.run(command,check=True)
    provenance={'command':command,'compiler':subprocess.check_output(['g++','--version'],text=True),
                'source_sha256':digest(REPO/'verification/three_spin_replay.cpp')}
    (out/'compiler.json').write_text(json.dumps(provenance,indent=2)+'\n')
    cases=[]
    for mode in range(4):
        for dt,current,field,seed in [(5,5e10,275,210),(1.25,5e10,1000,30),(5,-5e10,1000,330)]:
            name=f'probe_{mode}_{dt}_{int(current/5e10)}'
            spec=dict(name=name,mode=mode,dt_fs=dt,seed_deg=seed,field_Oe=field,current=current,
                      duration_ns=.005,on_ns=.001,width_ns=.003)
            run_case(out,spec,sample_ps=dt/1000)
            cases.append(spec)
    (out/'probe_cases.json').write_text(json.dumps(cases,indent=2)+'\n')
    print(f'Prepared {len(cases)} short source/replay comparisons')


def readbin(path):
    values=array.array('d')
    with path.open('rb') as f:
        values.frombytes(f.read())
    if sys.byteorder != 'little':
        values.byteswap()
    assert len(values)%10==0
    return values


def compare(out,a,b):
    av=readbin(out/(a+'.bin'));bv=readbin(out/(b+'.bin'))
    assert len(av)==len(bv)
    spin=mean=0.
    for n in range(0,len(av),10):
        assert abs(av[n]-bv[n])<1e-22
        for q in range(3):
            mean=max(mean,abs(sum(av[n+1+q+3*i]-bv[n+1+q+3*i] for i in range(3))/3))
        spin=max(spin,max(abs(av[n+q]-bv[n+q]) for q in range(1,10)))
    return {'max_spin_component_difference':spin,'max_mean_component_difference':mean}


def batch(out,cases):
    assert json.loads((out/'dynamics.json').read_text())['pass']
    with ThreadPoolExecutor(max_workers=4) as pool:
        jobs=[pool.submit(run_case,out,c) for c in cases]
        result=[]
        for job in jobs:
            item=job.result();result.append(item)
            print(f"{item['name']}: angle {item['final_angle_deg']:.8f}; pre-norm {item['max_preprojection_norm_error']:.3g}",flush=True)
    return result


def pilot(out):
    cases=[]
    for mode in range(4):
        for dt in ([5,2.5,1.25,.5] if mode<2 else [5,2.5,1.25]):
            cases.append(dict(name=f'pilot_m{mode}_dt{dt}',mode=mode,dt_fs=dt,seed_deg=210,
                              field_Oe=300,current=5e10,duration_ns=5,on_ns=1,width_ns=3))
    results=batch(out,cases)
    comparisons={}
    for mode in range(4):
        comparisons[f'm{mode}_5_2.5']=compare(out,f'pilot_m{mode}_dt5',f'pilot_m{mode}_dt2.5')
        comparisons[f'm{mode}_2.5_1.25']=compare(out,f'pilot_m{mode}_dt2.5',f'pilot_m{mode}_dt1.25')
    for dt in [5,2.5,1.25]:
        comparisons[f'coupled_reference_dt{dt}']=compare(out,f'pilot_m2_dt{dt}',f'pilot_m3_dt{dt}')
    for key,item in comparisons.items():
        if key.startswith('m2') or key.startswith('m3') or key.startswith('coupled'):
            assert item['max_mean_component_difference']<1e-5,(key,item)
    (out/'pilot.json').write_text(json.dumps({'pass':True,'runs':results,'comparisons':comparisons},indent=2)+'\n')
    print(json.dumps(comparisons,indent=2))


def physical_cases():
    for field in [175,200,225,250,275,300,325,350,375,400,450,500]:
        yield f'field{field}',210,field,5e10
    for sign,current in [('p',5e10),('m',-5e10)]:
        for seed in [30,90,150,210,270,330]:
            yield f'high{sign}{seed}',seed,1000,current


def sweep(out):
    pilot=json.loads((out/'pilot.json').read_text());assert pilot['pass']
    cases=[];previous=[]
    for label,seed,field,current in physical_cases():
        for mode in range(4):
            for dt in ([5,2.5,1.25,.5] if mode<2 else [5,2.5,1.25]):
                if label=='field300':
                    item=next(r for r in pilot['runs'] if r['mode']==mode and r['dt_fs']==dt)
                    previous.append({**item,'case':label})
                else:
                    cases.append(dict(name=f'{label}_m{mode}_dt{dt}',case=label,mode=mode,dt_fs=dt,
                        seed_deg=seed,field_Oe=field,current=current,duration_ns=5,on_ns=1,width_ns=3))
    results=previous+batch(out,cases)
    (out/'sweep.json').write_text(json.dumps(results,indent=2)+'\n')
    print(f'PASS completed {len(results)} full trajectories (includes {len(previous)} pilot runs)')


def analyze(out):
    runs=json.loads((out/'sweep.json').read_text());summaries=[];norm_error=0
    index={(r['case'],r['mode'],r['dt_fs']):r for r in runs}
    for r in runs:
        data=readbin(out/(r['name']+'.bin'));assert len(data)==20001*10
        for n in range(0,len(data),10):
            assert abs(data[n]-(n//10)*.25e-12)<1e-22
            assert all(math.isfinite(v) for v in data[n:n+10])
            for i in range(3):
                norm_error=max(norm_error,abs(math.sqrt(sum(data[n+1+i*3+q]**2 for q in range(3)))-1))
        assert r['max_postprojection_norm_error']<1e-12
        assert r['updates']==round(5e6/r['dt_fs'])
        expected=30 if r['field_Oe']<=275 else 90
        angle_error=abs((r['final_angle_deg']-expected+180)%360-180)
        if r['mode']>=2: assert angle_error<.1,(r['name'],angle_error)
        ref=index[(r['case'],3,1.25)]
        discrepancy=compare(out,r['name'],ref['name'])
        driven=data[16000*10:16001*10]  # State at precisely 4 ns, after last driven update.
        my=sum(driven[2+3*i] for i in range(3))/3;mz=sum(driven[3+3*i] for i in range(3))/3
        summaries.append({k:r[k] for k in ['name','case','mode','dt_fs','field_Oe','current','seed_deg',
                         'updates','final_angle_deg','max_preprojection_norm_error','max_postprojection_norm_error']} |
                         {'angle_error_deg':angle_error,'angle_at_4ns_deg':math.degrees(math.atan2(my,mz)),**discrepancy})
    comparisons=[]
    for label,*_ in physical_cases():
        for mode in range(4):
            for coarse,fine in [(5,2.5),(2.5,1.25)]+([(1.25,.5)] if mode<2 else []):
                a=index[(label,mode,coarse)];b=index[(label,mode,fine)]
                result=compare(out,a['name'],b['name'])
                angle=abs((a['final_angle_deg']-b['final_angle_deg']+180)%360-180)
                if mode>=2:assert result['max_mean_component_difference']<1e-5 and angle<.1
                comparisons.append(dict(case=label,mode=mode,coarse_fs=coarse,fine_fs=fine,
                    final_angle_difference_deg=angle,**result))
    assert norm_error<1e-12
    with (out/'summary.csv').open('w',newline='') as f:
        writer=csv.DictWriter(f,fieldnames=list(summaries[0]));writer.writeheader();writer.writerows(summaries)
    analysis={'pass':True,'full_runs':len(runs),'sampled_norm_error':norm_error,'comparisons':comparisons,
              'max_preprojection_norm_error':max(r['max_preprojection_norm_error'] for r in runs)}
    (out/'analysis.json').write_text(json.dumps(analysis,indent=2)+'\n')
    print(f'PASS {len(runs)} trajectories: physical times, finite spins, norms, coupled destinations/convergence')


def refine(out):
    assert json.loads((out/'pilot.json').read_text())['pass']
    cases=[dict(name=f'refine_m{mode}_dt{dt}',mode=mode,dt_fs=dt,seed_deg=210,
                field_Oe=300,current=5e10,duration_ns=5,on_ns=1,width_ns=3)
           for mode in [0,1] for dt in [.25,.125,.0625]]
    runs=batch(out,cases);comparisons=[]
    for mode in [0,1]:
        previous=f'pilot_m{mode}_dt0.5'
        for dt in [.25,.125,.0625]:
            name=f'refine_m{mode}_dt{dt}'
            comparisons.append(dict(mode=mode,dt_fs=dt,previous=previous,name=name,
                successive=compare(out,previous,name),reference=compare(out,name,'pilot_m3_dt1.25')))
            previous=name
    (out/'refinement.json').write_text(json.dumps({'runs':runs,'comparisons':comparisons},indent=2)+'\n')
    print(json.dumps(comparisons,indent=2))


def replays(out):
    runs=json.loads((out/'sweep.json').read_text());cases=[]
    for case in ['field300','highm330']:
        for mode in range(4):
            r=next(r for r in runs if r['case']==case and r['mode']==mode and r['dt_fs']==5)
            data=readbin(out/(r['name']+'.bin'));start=6000  # Exactly 1.5 ns.
            seed=data[start*10+1:start*10+10];name=f'driven_{case}_m{mode}'
            seedfile=name+'.csv'
            with (out/seedfile).open('w',newline='') as f:
                csv.writer(f).writerows([seed[i*3:(i+1)*3] for i in range(3)])
            spec=dict(name=name,mode=mode,dt_fs=5,seed_deg=r['seed_deg'],field_Oe=r['field_Oe'],
                current=r['current'],duration_ns=.005,on_ns=0,width_ns=.005,seed_file=seedfile)
            run_case(out,spec,sample_ps=.005)
            actual=readbin(out/(name+'.bin'));error=0
            for k in range(21):
                error=max(error,max(abs(actual[k*50*10+q]-data[(start+k)*10+q]) for q in range(1,10)))
            assert error<1e-14,(name,error)
            spec['restart_error']=error;cases.append(spec)
    (out/'driven_cases.json').write_text(json.dumps(cases,indent=2)+'\n')
    print('PASS 8 restarted replay segments match their full trajectories')


if __name__ == '__main__':
    {'prepare': prepare, 'preserve': preserve, 'probes':probes,'pilot':pilot,
     'sweep':sweep,'analyze':analyze,'refine':refine,'replays':replays}[sys.argv[1]](Path(sys.argv[2]).resolve())
