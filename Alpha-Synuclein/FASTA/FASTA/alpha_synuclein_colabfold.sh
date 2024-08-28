#!/bin/bash
#SBATCH --nodes=1 --ntasks-per-node=1 --partition=gpu # Pod cluster's GPU queue
#SBATCH --gres=gpu:1
#SBATCH --time=05:00:00
#SBATCH --job-name=alpha_synuclein_testing
#SBATCH --mail-user=leifgriem@ucsb.edu # uncomment these two lines and include email if desired
#SBATCH --mail-type=END,FAIL    # Send email at begin and end of job

cd $SLURM_SUBMIT_DIR
export PATH="/sw/alphafold/localcolabfold/colabfold-conda/bin:$PATH"
module load cuda/11.2

srun --gres=gpu:1 colabfold_batch alpha_synuclein_fragments.fasta alpha_synuclein_15aa_consecutive 
