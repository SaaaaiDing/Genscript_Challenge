#!/usr/bin/env bash
set -euo pipefail
TASK_RUN=/storage/group/yjz5549/default/GenScript_challenge/computational
REPO=/storage/group/yjz5549/default/RFD4/RFD4-Proteina
export CLUSTER=local
export PYTHONDONTWRITEBYTECODE=1 PYTHONUNBUFFERED=1
export HF_HOME=/storage/group/yjz5549/default/RFD4/huggingface
export HF_HUB_OFFLINE=1 TRANSFORMERS_OFFLINE=1 WANDB_MODE=disabled
export XDG_CACHE_HOME="$TASK_RUN/cache"
export TORCH_HOME="$TASK_RUN/cache/torch"
export TORCHINDUCTOR_CACHE_DIR="$TASK_RUN/cache/torchinductor"
export TRITON_CACHE_DIR="$TASK_RUN/cache/triton"
export CUDA_CACHE_PATH="$TASK_RUN/cache/cuda"
export MPLCONFIGDIR="$TASK_RUN/cache/matplotlib"
export TMPDIR="$TASK_RUN/tmp"
export APPTAINER_CACHEDIR="$TASK_RUN/cache/apptainer"
export APPTAINER_TMPDIR="$TASK_RUN/tmp/apptainer"
export AF3_JAX_CACHE_DIR="$TASK_RUN/cache/of3_jax"
export AF3_SOCKET_DIR="$TASK_RUN/sockets"
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
export PATH="/storage/group/yjz5549/default/RFD4/bin:$PATH"
mkdir -p "$TMPDIR" "$APPTAINER_TMPDIR" "$AF3_SOCKET_DIR"
cd "$REPO"
source .venv/bin/activate
source scripts/load_local_env.sh
exec "$@"
