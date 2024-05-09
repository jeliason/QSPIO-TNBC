#!/bin/bash
#SBATCH --job-name=PSA_simulations
#SBATCH --array=0-20
#SBATCH --output=output_logs/output_%a.log
#SBATCH --error=error_logs/error_%a.log
#SBATCH --mem-per-cpu=4G
#SBATCH --time=00:15:00
#SBATCH --mail-user=joelne@umich.edu
#SBATCH --mail-type=ALL

module load python3.10-anaconda
module load matlab/R2023b

source /home/${USER}/.bashrc
conda activate matlab_QSP

python PSA_NSCLC_iteration.py --task_index $SLURM_ARRAY_TASK_ID --n_param_sets 30 --write_directory /nfs/turbo/umms-ukarvind/joelne/PSA_simulations

# Check if there are any errors in the error log
if [ -s error_logs/error_${SLURM_ARRAY_TASK_ID}.log ]; then
    echo "Python script encountered errors. Error log contents:"
    cat error_logs/error_${SLURM_ARRAY_TASK_ID}.log >> error_logs/python_errors_${SLURM_ARRAY_TASK_ID}.log
fi

# Remove the original error log file
rm error_logs/error_${SLURM_ARRAY_TASK_ID}.log
