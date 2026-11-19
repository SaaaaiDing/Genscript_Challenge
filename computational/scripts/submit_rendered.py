#!/usr/bin/env python3
"""Submit one already-inspected job and retain the exact job ID and script hash."""
import datetime
import hashlib
import json
import subprocess
import sys
from pathlib import Path
root=Path(__file__).resolve().parents[1]
script=Path(sys.argv[1])
receipt=root/'provenance/submission.json'
if receipt.exists():raise SystemExit('A submission receipt already exists; refusing duplicate submission.')
result=subprocess.run(['/storage/sys/slurm/bin/sbatch','--parsable',str(script)],cwd=script.parent,capture_output=True,text=True,check=True)
job_id=result.stdout.strip().split(';')[0]
assert job_id.isdigit(), result.stdout
data={'job_id':job_id,'submitted_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),
      'script':str(script),'script_sha256':hashlib.sha256(script.read_bytes()).hexdigest(),
      'qos':'preempt','excluded_node':'p-gc-3039','gpus':1,'cpus':2,'memory_gib':24,
      'requested_designs':9,'representative_smoke_designs':1,'exploratory_designs':8,'design_problems':4,'replicates_per_problem':2,
      'status':'submitted_not_yet_confirmed_running'}
receipt.write_text(json.dumps(data,indent=2)+'\n')
print(json.dumps(data,indent=2))
