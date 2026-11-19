#!/usr/bin/env bash
set -euo pipefail
ROOT=/storage/group/yjz5549/default/GenScript_challenge/computational
REPO=/storage/group/yjz5549/default/RFD4/RFD4-Proteina
RUN="$ROOT/runs/initial_designs"
test -f "$ROOT/provenance/paired_input_validation.json"
test ! -e "$ROOT/provenance/submission.json"
export ACCOUNT=yjz5549_cr_default PARTITION=standard
export SLURM_EXTRA_OPTIONS='--qos=preempt'
export EXCLUDE_NODES=p-gc-3039
export PATH="$ROOT/scripts/render_only:$PATH"
set +e
bash "$ROOT/scripts/environment.sh" bash "$REPO/scripts/launch_generate.sh" \
  --njobs 1 --gpus-per-job 1 --no-concatenate --no-foldseek \
  --gres gpu:a100:1 --cpus-per-task 2 --mem 24G --time 01:00:00 \
  --config-path "$ROOT/configs" --config-name pilot_generation \
  "inputs=$ROOT/inputs/annotated" "out_dir=$RUN" \
  "hf_ckpt_path=/storage/group/yjz5549/default/RFD4/huggingface/hub/models--nvidia--RFD4-Proteina/snapshots/a06e2509c7033b010fc43614e9fd840303588502" \
  num_replicates=2 compress_outputs=false \
  > "$ROOT/logs/render_initial_job.log" 2>&1
render_status=$?
set -e
test "$render_status" -eq 97
test -f "$RUN/slurm_batch_script.sh"
python3 "$ROOT/scripts/finalize_job.py" "$RUN/slurm_batch_script.sh"
bash -n "$RUN/slurm_batch_script.sh"
# Remove the renderer and inherited allocation variables before the real submit.
export PATH="${PATH#"$ROOT/scripts/render_only:"}"
while IFS= read -r variable; do unset "$variable"; done < <(compgen -e | rg '^(SLURM|SBATCH)_')
bash "$ROOT/scripts/environment.sh" python -B "$ROOT/scripts/submit_rendered.py" "$RUN/slurm_batch_script.sh"
