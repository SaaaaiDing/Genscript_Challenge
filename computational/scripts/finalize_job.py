#!/usr/bin/env python3
"""Inspect and harden only the project-owned rendered job before submission."""
import sys
from pathlib import Path
p=Path(sys.argv[1])
s=p.read_text()
# The upstream launcher strips inherited SLURM_* variables, including its extra-options
# variable. Add the reviewed QoS explicitly to this project-owned rendered script.
assert '#SBATCH --qos' not in s
s=s.replace('#SBATCH --partition standard', '#SBATCH --partition standard\n#SBATCH --qos=preempt')
s=s.replace('#SBATCH --requeue', '#SBATCH --no-requeue')
required=['--qos=preempt','--exclude=p-gc-3039','--gres=gpu:a100:1','--cpus-per-task 2','--mem 24G']
assert all(x in s for x in required), 'Missing resource isolation directive.'
assert '#SBATCH --array' not in s
assert 'af3_daemon.py' not in s and 'af2_daemon.py' not in s, 'Generation must not launch refolding daemons.'
s=s.replace('export OMP_NUM_THREADS=4','export OMP_NUM_THREADS=1').replace('export OPENBLAS_NUM_THREADS=4','export OPENBLAS_NUM_THREADS=1')
lines=s.splitlines()
last=max(i for i,line in enumerate(lines) if line.startswith('#SBATCH'))
guard='''
set -euo pipefail
if [[ "$(hostname -s)" == p-gc-3039 ]]; then
    echo 'Refusing to use the protected node.' >&2
    exit 1
fi
TASK_ROOT=/storage/group/yjz5549/default/GenScript_challenge/computational
test ! -e "$TASK_ROOT/STOP"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export PYTHONPATH="/storage/group/yjz5549/default/RFD4/RFD4-Proteina/src:${PYTHONPATH:-}"
bash "$TASK_ROOT/scripts/environment.sh" bash skills/rfd4-proteina-setup/scripts/preflight.sh
bash "$TASK_ROOT/scripts/environment.sh" rfd4-generate \\
    --config-path "$TASK_ROOT/configs" --config-name pilot_generation \\
    "inputs=$TASK_ROOT/inputs/annotated/il4ra_p1_h3_16.cif" \\
    "out_dir=$TASK_ROOT/runs/representative_smoke" \\
    "hf_ckpt_path=/storage/group/yjz5549/default/RFD4/huggingface/hub/models--nvidia--RFD4-Proteina/snapshots/a06e2509c7033b010fc43614e9fd840303588502" \\
    num_replicates=1 compress_outputs=false
bash "$TASK_ROOT/scripts/environment.sh" python -B "$TASK_ROOT/scripts/validate_smoke.py"
test ! -e "$TASK_ROOT/STOP"
'''
lines.insert(last+1,guard)
p.write_text('\n'.join(lines)+'\n')
print('Job resource directives, serial generation, and GPU preflight verified.')
