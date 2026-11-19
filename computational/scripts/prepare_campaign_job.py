#!/usr/bin/env python3
"""Render the official launcher, then attach the bounded project coordinator."""
import hashlib
import json
import os
import shlex
import subprocess
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT / 'campaign_36h_v2'
REPO = Path('/storage/group/yjz5549/default/RFD4/RFD4-Proteina')
assert (BASE / 'campaign.sqlite').exists()
assert not (BASE / 'submission.json').exists(), 'Refusing duplicate submission.'
env = {k: v for k, v in os.environ.items() if not k.startswith(('SLURM_', 'SBATCH_'))}
env.update(ACCOUNT='yjz5549_cr_default', PARTITION='standard', EXCLUDE_NODES='p-gc-3039',
           PATH=str(ROOT / 'scripts/render_only') + ':' + env['PATH'])
args = ['bash', str(REPO / 'scripts/launch_generate.sh'), '--njobs', '1', '--gpus-per-job', '1',
        '--no-concatenate', '--no-foldseek', '--gres', 'gpu:a100:1', '--cpus-per-task', '2',
        '--mem', '48G', '--time', '36:00:00', '--config-path', str(ROOT / 'configs'),
        '--config-name', 'pilot_generation', f'inputs={ROOT}/inputs/annotated', f'out_dir={BASE}',
        'run_name_=il4ra_1600_36h', 'num_replicates=400', 'compress_outputs=false']
result = subprocess.run(args, env=env, text=True, capture_output=True)
(BASE / 'launcher_render.log').write_text(result.stdout + result.stderr)
assert result.returncode == 97, result.stderr[-2000:]
script = BASE / 'slurm_batch_script.sh'
rendered = script.read_text()
(BASE / 'official_launcher_render.sh').write_text(rendered)
header = [line for line in rendered.splitlines() if (line.startswith('#!') or line.startswith('#SBATCH')) and line != '#SBATCH --requeue']
header += ['#SBATCH --qos=preempt', '#SBATCH --no-requeue']
body = '''
set -euo pipefail
if [[ "$(hostname -s)" == p-gc-3039 ]]; then
    echo 'Refusing the protected node.' >&2
    exit 1
fi
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
ulimit -n "$(ulimit -Hn)" 2>/dev/null || true
# One bounded coordinator serializes generation and independent refolding.
exec bash {environment} python -B {controller}
'''.format(environment=shlex.quote(str(ROOT / 'scripts/environment.sh')),
           controller=shlex.quote(str(ROOT / 'scripts/campaign36.py')))
script.write_text('\n'.join(header) + '\n' + body)
subprocess.run(['bash', '-n', str(script)], check=True)
assert '--exclude=p-gc-3039' in script.read_text()
assert '--array' not in script.read_text()
(BASE / 'slurm_run_outputs').mkdir(exist_ok=True)
print(script.read_text())
