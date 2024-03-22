#!/bin/bash 
#SBATCH --nodes=1 
#SBATCH --ntasks-per-node=20
#SBATCH --mem=400GB 
#SBATCH --time=48:00:00 
#SBATCH --partition=open
#SBATCH --output STDOUT.out
#SBATCH --error STDERR.out
#SBATCH --mail-type=end   # send email when job ends
#SBATCH --mail-type=fail  # send email if job fails
#SBATCH --mail-user=dzg5526@psu.edu

module purge
module load matlab
module load abaqus

cd $SLURM_SUBMIT_DIR
matlab -nodisplay -nosplash < ResonatorDesign.m > log.matlab

# Jobs managment
# sbatch shellfile.sh
# squeue --user abc123
# scancel 123456
