"""Source-checked diagnostic replay; standard-library Python, g++, prepared MATLAB inputs."""
import argparse
import array
from concurrent.futures import ThreadPoolExecutor
import hashlib
import json
import math
from pathlib import Path
import random
import subprocess
import sys

parser = argparse.ArgumentParser()
parser.add_argument('output', type=Path)
parser.add_argument('--check-only', action='store_true')
args = parser.parse_args()
out = args.output.resolve()
here = Path(__file__).resolve().parent
models = sorted(out.glob('*.model'))
assert len(models) == 6, 'First run prepare_spatial into a fresh output directory.'
if not args.check_only:
    assert not list(out.glob('*-5fs-5ns.bin')), 'Use a fresh directory to preserve prior runs.'
subprocess.run(['g++', '-O3', '-std=c++17', '-Wall', '-Wextra', str(here/'spatial.cpp'),
                '-o', str(out/'spatial')], check=True)
manifest = {str(f.relative_to(here.parent.parent)): hashlib.sha256(f.read_bytes()).hexdigest()
            for f in [*here.glob('*.*'), *here.parent.parent.glob('*.m')]}
(out/'diagnostic-hashes.json').write_text(json.dumps(manifest, indent=2))


def read(path):
    a = array.array('d')
    a.frombytes(path.read_bytes())
    return a


def run(command):
    r = subprocess.run(command, check=True, capture_output=True, text=True)
    print(r.stdout.strip(), flush=True)


gates = []
for model in models:
    lines = model.read_text().splitlines()
    h = lines[0].split()
    n = int(h[0])
    J, D, K = map(float, h[3:6])
    e = [list(map(float, line.split()))[:3] for line in lines[1:n+1]]
    bonds = {}
    for line in lines[n+2:]:
        i, j, a, c = map(float, line.split())
        bonds[int(i), int(j)] = (J*.5 if a else 0, math.copysign(D*.5, c) if c else 0)
    for (i, j), (a, c) in bonds.items():
        assert bonds[j, i] == (a, -c)

    def energy(s):
        val = -.5*K*sum(sum(x*y for x, y in zip(s[i], e[i]))**2 for i in range(n))
        for (i, j), (a, c) in bonds.items():
            val += .5*(a*sum(x*y for x, y in zip(s[i], s[j]))
                       - c*(s[i][1]*s[j][2]-s[i][2]*s[j][1]))
        return val

    rng = random.Random(20261102)
    s = [[rng.uniform(-1, 1) for _ in range(3)] for _ in range(n)]
    field = [[K*sum(x*y for x, y in zip(s[i], e[i]))*a for a in e[i]] for i in range(n)]
    for (i, j), (a, c) in bonds.items():
        for k in range(3):
            field[i][k] -= a*s[j][k]
        field[i][1] += c*s[j][2]
        field[i][2] -= c*s[j][1]
    error = 0
    for i in sorted({0, n//3, n//2, n-1}):
        for k in range(3):
            old = s[i][k]
            eps = 1e-5
            s[i][k] = old+eps
            ep = energy(s)
            s[i][k] = old-eps
            em = energy(s)
            s[i][k] = old
            error = max(error, abs((ep-em)/(2*eps)+field[i][k]))
    assert error < 1e-6
    replay = []
    for trial in range(4):
        mode = 'current' if trial < 2 else 'reciprocal'
        tag = f'{model.stem}-replay{trial}'
        gpu = out/f'{model.stem}-gpu{trial}.bin'
        run([str(out/'spatial'), str(model), mode, '5', '.1', '.005', str(out/tag), str(gpu)])
        a, b = read(gpu), read(out/f'{tag}.bin')
        assert len(a) == len(b) == 21*3*n
        err = max(abs(x-y) for x, y in zip(a, b))
        assert err < 1e-10
        replay.append(err)
    gates.append(dict(size=model.stem, gradient_error_T=error, replay_errors=replay))
(out/'verified-gates.json').write_text(json.dumps(gates, indent=2))

if not args.check_only:
    jobs = []
    followups = []
    for model in models:
        for mode in (['current'] if model.stem == '2x2' else ['current', 'reciprocal']):
            tag = f'{model.stem}-{mode}'
            cmd = [str(out/'spatial'), str(model), mode]
            jobs.append(cmd+['5', '5000', '100', str(out/f'{tag}-5fs-5ns')])
            followups.append(cmd+['2.5', '1000', '100', str(out/f'{tag}-2.5fs-1ns')])
            followups.append(cmd+['5', '2000', '100', str(out/f'{tag}-recovery'),
                                 str(out/f'{tag}-5fs-5ns-final.bin'), '.001'])
            followups.append(cmd+['5', '2000', '100', str(out/f'{tag}-perturbed-start'), '-', '.01'])
            if model.stem in ['5x5', '21x22']:
                for dt in ['5', '2.5', '0.5']:
                    followups.append(cmd+[dt, '10', '.1', str(out/f'{tag}-{dt}fs-10ps')])
    (out/'commands.json').write_text(json.dumps(jobs+followups, indent=2))
    for batch in [jobs, followups]:
        with ThreadPoolExecutor(max_workers=3) as pool:
            list(pool.map(run, batch))
subprocess.run([sys.executable, str(here/'summarize.py'), str(out)], check=True)
