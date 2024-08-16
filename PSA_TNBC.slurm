#!/bin/bash
#SBATCH --job-name=PSA_simulations
#SBATCH --array=0-1
#SBATCH --output=output_logs/output_%a.log
#SBATCH --error=error_logs/error_%a.log
#SBATCH --mem-per-cpu=2G
#SBATCH --time=00:20:00
#SBATCH --mail-user=joelne@umich.edu
#SBATCH --mail-type=ALL

SAVE_FOLDER=/nfs/turbo/umms-ukarvind/joelne/PSA_simulations/TNBC
NUM_PER_PROCESS=5
START_INDEX=$(($SLURM_ARRAY_TASK_ID))*NUM_PER_PROCESS+1
STOP_INDEX=$(($SLURM_ARRAY_TASK_ID+1))*NUM_PER_PROCESS

echo $START_INDEX
echo $STOP_INDEX

module load matlab/R2023b

matlab -nojvm -nosplash -batch "addpath(genpath('~/repositories/QSPIO-TNBC/')); PSA_iteration_TNBC(${START_INDEX},${STOP_INDEX},${SAVE_FOLDER}); exit"
