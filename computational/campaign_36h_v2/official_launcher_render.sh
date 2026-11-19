#!/bin/bash
#SBATCH --account yjz5549_cr_default
#SBATCH --job-name il4ra_1600_36h
#SBATCH --partition standard
#SBATCH --time 36:00:00
#SBATCH --nodes 1
#SBATCH --gres=gpu:a100:1
#SBATCH --ntasks-per-node 1
#SBATCH --cpus-per-task 2
#SBATCH --mem 48G
#SBATCH --open-mode append
#SBATCH --requeue
#SBATCH --exclude=p-gc-3039
#SBATCH --output ./slurm_run_outputs/slurm.out
#SBATCH --error ./slurm_run_outputs/slurm.err

# Raise the open-file soft limit to the hard cap on hosts with low defaults.
ulimit -n "$(ulimit -Hn)" 2>/dev/null || true
cd /storage/group/yjz5549/default/RFD4/RFD4-Proteina

# Environment setup (all clusters)
export PYTHONPATH=/storage/group/yjz5549/default/RFD4/RFD4-Proteina:$PYTHONPATH
export PYTHONFAULTHANDLER=1
export OMP_NUM_THREADS=4
export OPENBLAS_NUM_THREADS=4

# Job array configuration
ARRAY_TASK_ID=${SLURM_ARRAY_TASK_ID:-0}
GPUS_PER_JOB=1
NJOBS=1

# --exclusive below gives each step its own GPU, but caps it at one CPU unless asked.
CPUS_PER_STEP=$(( ${SLURM_CPUS_ON_NODE:-8} / GPUS_PER_JOB ))
CPUS_PER_STEP=$(( CPUS_PER_STEP >= 2 ? CPUS_PER_STEP : 2 ))  # rule of thumb: never < 2 CPUs/GPU

echo "Array task: $ARRAY_TASK_ID, GPUs per job: $GPUS_PER_JOB, total jobs: $NJOBS, CPUs per step: $CPUS_PER_STEP"

# Cluster identity — read by exec_core.sh (shebang) to source .ipd/.env.${CLUSTER}
export CLUSTER=local
export CONTAINER_TYPE=none

export MASTER_PORT=$((1024 + RANDOM % 64512))
export MASTER_ADDR=$(scontrol show hostnames "$SLURM_JOB_NODELIST" | head -n 1)

# Verify scontrol works (needed for auto-requeue)
scontrol show job $SLURM_JOB_ID > /dev/null 2>&1 && echo '[slurm] scontrol works - auto-requeue enabled' || echo '[slurm] WARNING: scontrol not working - auto-requeue disabled'

# Run generation - each GPU gets its own job_id slice of the dataset
GEN_PIDS=()
for LOCAL_ID in $(seq 0 $((GPUS_PER_JOB - 1))); do
    JOB_ID=$((ARRAY_TASK_ID * GPUS_PER_JOB + LOCAL_ID))
    if [ $JOB_ID -lt $NJOBS ]; then
        echo "Running generation with job_id=$JOB_ID on GPU $LOCAL_ID"
        AF3_DAEMON_RANK=$LOCAL_ID AF2_DAEMON_RANK=$LOCAL_ID srun --exclusive --ntasks=1 --gpus=1 \
            --cpus-per-task=$CPUS_PER_STEP --cpu-bind=none \
            --kill-on-bad-exit /storage/group/yjz5549/default/RFD4/RFD4-Proteina/src/rfproteina/generate.py \
            --config-path /storage/group/yjz5549/default/GenScript_challenge/computational/configs --config-name pilot_generation \
             inputs=/storage/group/yjz5549/default/GenScript_challenge/computational/inputs/annotated out_dir=/storage/group/yjz5549/default/GenScript_challenge/computational/campaign_36h_v2 run_name_=il4ra_1600_36h num_replicates=400 compress_outputs=false \
            out_dir=/storage/group/yjz5549/default/GenScript_challenge/computational/campaign_36h_v2 \
            hardware.ngpus_per_node_=1 \
            hardware.nnodes_=1 \
            job_id=$JOB_ID \
            gen_njobs=$NJOBS &
        GEN_PIDS+=($!)
    fi
done
wait "${GEN_PIDS[@]}"
