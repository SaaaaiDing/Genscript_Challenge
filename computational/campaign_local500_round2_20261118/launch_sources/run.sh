#!/usr/bin/env bash
set -euo pipefail
BASE=/storage/group/yjz5549/default/GenScript_challenge/computational/campaign_local500_round2_20261118
exec 9>"$BASE/controller.lock"
flock -n 9 || { echo 'Round 2 is already running.' >&2; exit 1; }
ulimit -n "$(ulimit -Hn)" 2>/dev/null || true
exec srun --jobid=57427235 --overlap --exact --nodes=1 --ntasks=1 \
  --cpus-per-task=2 --gres=gpu:a100:1 --mem=48G --time=16:00:00 \
  --job-name=il4ra_second500 \
  bash /storage/group/yjz5549/default/GenScript_challenge/computational/scripts/environment.sh \
  python -B /storage/group/yjz5549/default/GenScript_challenge/computational/scripts/campaign_round2_500.py
