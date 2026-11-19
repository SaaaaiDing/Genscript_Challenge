#!/usr/bin/env python3
"""Submit the reviewed campaign once and retain a content-addressed receipt."""
import ast
import datetime
import hashlib
import json
import os
import shutil
import sqlite3
import subprocess
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT / 'campaign_36h_v2'
script = BASE / 'slurm_batch_script.sh'
receipt_path = BASE / 'submission.json'
assert not receipt_path.exists(), 'Refusing duplicate submission.'
assert '#SBATCH --requeue\n' not in script.read_text()
assert '#SBATCH --qos=preempt' in script.read_text()
assert '#SBATCH --exclude=p-gc-3039' in script.read_text()
ast.parse((ROOT / 'scripts/campaign36.py').read_text())
subprocess.run(['bash', '-n', str(script)], check=True)
env = {k: v for k, v in os.environ.items() if not k.startswith(('SLURM_', 'SBATCH_'))}
result = subprocess.run(['/storage/sys/slurm/bin/sbatch', '--parsable', str(script)], cwd=BASE,
                        env=env, text=True, capture_output=True, check=True)
job = result.stdout.strip().split(';')[0]
assert job.isdigit(), result.stdout
receipt = dict(job_id=job, submitted_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),
               requested_designs=1600, time_limit='36:00:00', qos='preempt', excluded_node='p-gc-3039',
               gpus=1, cpus=2, memory_gib=48, replaces_cancelled_job='57451588',
               script_sha256=hashlib.sha256(script.read_bytes()).hexdigest(),
               controller_sha256=hashlib.sha256((ROOT / 'scripts/campaign36.py').read_bytes()).hexdigest())
receipt_path.write_text(json.dumps(receipt, indent=2) + '\n')
with sqlite3.connect(BASE / 'campaign.sqlite') as conn:
    conn.execute('UPDATE runs SET slurm_job_id=? WHERE phase=?', (job, 'baseline'))
archive = BASE / 'submission_sources'
archive.mkdir(exist_ok=True)
for path in [ROOT / 'scripts/campaign36.py', ROOT / 'scripts/quality_screen.py', ROOT / 'scripts/environment.sh',
             ROOT / 'configs/paired_evaluation.yaml', ROOT / 'inputs/raw/known_il4ra_cdrs.json']:
    shutil.copy2(path, archive / path.name)
print(json.dumps(receipt, indent=2))
