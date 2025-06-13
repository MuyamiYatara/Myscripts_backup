#!/bin/bash -l 


#SBATCH -A hmt03
#SBATCH --time=72:00:00
##SBATCH --exclude=hpcc[154,155,156,114]
#SBATCH --nodelist=node020
#SBATCH --nodes=1
#SBATCH --partition=regular6430
#SBATCH --cpus-per-task=1
#SBATCH --ntasks-per-node=32
#SBATCH --job-name=NaN
#SBATCH --output=./log
#SBATCH --error=./errormsg
export OMP_NUM_THREADS=1 
export MKL_NUM_THREADS=1 
export MV2_ENABLE_AFFINITY=0
echo "The current job ID is $SLURM_JOB_ID"
echo "Running on $SLURM_JOB_NUM_NODES nodes:"
echo $SLURM_JOB_NODELIST
echo "Using $SLURM_NTASKS_PER_NODE tasks per node"
echo "A total of $SLURM_NTASKS tasks is used"

ulimit -s unlimited
ulimit -c unlimited
module purge
module load oneapi22.3
module load mkl/mkl2022.2.1
ulimit -s unlimited
sleep 1d