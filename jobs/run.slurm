#!/bin/bash
#SBATCH -J juliaserial
#SBATCH --account=arpitdua
#SBATCH --partition=normal_q
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --time=3-00:00:00
#SBATCH --output=logs/%A.%a.out
#SBATCH --error=logs/%A.%a.err

set -euo pipefail
mkdir -p output logs

PARAMS=$(sed -n "$((SLURM_ARRAY_TASK_ID + ${OFFSET:-0}))p" jobs/new_failed_params.txt)

module reset
module load Julia/1.11.3-linux-x86_64

echo "PARAMS: $PARAMS"
julia jobs/run.jl $PARAMS
