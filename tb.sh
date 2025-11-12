#!/bin/bash
#SBATCH -p stats
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4000
#SBATCH --time=48:00:00
#SBATCH --array 1-20
#SBATCH --mail-user=alicia.gill@warwick.ac.uk
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --output=/storage/stats/maundl/Paper/Output/tb_%a.stdout
#SBATCH --error=/storage/stats/maundl/Paper/Error/tb_%a.stderr

module purge
module load GCC/11.3.0 OpenMPI/4.1.4 R/4.2.1

export TASK_ID=$SLURM_ARRAY_TASK_ID

Rscript /storage/stats/maundl/Paper/tb.R