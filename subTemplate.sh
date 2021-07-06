#!/bin/bash
#SBATCH --account=ACCOUNT
#SBATCH --gres=gpu:v100:1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=1-0:00
#SBATCH --mail-user=EMAIL_ADDRESS
#SBATCH --mail-type=END

# LOAD MODULES RELEVANT FOR YOUR INSTALLATION
module load StdEnv/2016.4
module load pgi/19.4
module load cuda/10.0.130

#export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export LD_LIBRARY_PATH=MMB_PATH # put your MMB path here
source ~/tinkermd/bin/activate # activate to your virtual environment

python ./main.py --run_num=0
