#!/bin/bash
#SBATCH -J juliaserial
#SBATCH --account=arpitdua
#SBATCH --partition=normal_q
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --time=0-01:00:00
#SBATCH --output=logs/%A.%a.out
#SBATCH --error=logs/%A.%a.err

set -euo pipefail
mkdir -p output logs

module reset
module load Julia/1.11.3-linux-x86_64

julia jobs/collate.jl output/IsingML_fineq_sweep
