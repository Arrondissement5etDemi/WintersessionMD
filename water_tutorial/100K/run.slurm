#!/bin/bash
#SBATCH --job-name=1         # create a short name for your job
#SBATCH --nodes=1                # node count
#SBATCH --ntasks=1               # total number of tasks across all nodes
#SBATCH --cpus-per-task=1        # cpu-cores per task (>1 if multi-threaded tasks)
#SBATCH --mem-per-cpu=4G         # memory per cpu-core (4G is default)
##SBATCH --constraint=gpu80
##SBATCH --gpu-mps
##SBATCH --gres=gpu:4             # number of gpus per node
#SBATCH --time=00:30:00          # total run time limit (HH:MM:SS)
##SBATCH --mail-type=begin        # send email when job begins
##SBATCH --mail-type=end          # send email when job ends
##SBATCH --mail-user=yifanl@princeton.edu

module purge
module load anaconda3/2022.5
conda activate my-lammps-env

lmp -in spce.lmp
