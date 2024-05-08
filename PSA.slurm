#!/bin/bash
#SBATCH --job-name=PSA_simulations
#SBATCH --array=0-20
#SBATCH --output=output_%a.log  # Specifies the output file, where %A is the job ID and %a is the task ID
#SBATCH --mem-per-cpu=4G
#SBATCH --time=00:15:00
#SBATCH --mail-user=joelne@umich.edu
#SBATCH --mail-type=ALL

module load python3.10-anaconda
module load matlab/R2023b

source /home/${USER}/.bashrc
conda activate matlab_QSP

python PSA_NSCLC_iteration.py --slurm_task_id $SLURM_ARRAY_TASK_ID --n_param_sets 30 --write_directory /nfs/turbo/umms-ukarvind/joelne/PSA_simulations
