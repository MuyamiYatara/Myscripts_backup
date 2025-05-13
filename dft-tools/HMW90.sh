#!/bin/bash -l
#
#SBATCH -N 2
#SBATCH -p regular
#SBATCH --ntasks-per-node=52
#SBATCH --job-name=NaN
#SBATCH --output=./log
#SBATCH --error=./err
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
#module load oneapi2022/mkl/2022.1.0
module load pmix/2.2.2
module load parallel_studio/2020.2.254
#module load intel/20.2.254
#module load openmpi3/3.1.4
#module load ohpc
module swap intel gnu8
module load intelmpi/2020.2.254
module list
#module load prun/1.3
#module load vasp6/6.1

export LD_LIBRARY_PATH=/usr/local/lib64:$LD_LIBRARY_PATH
echo $LD_LIBRARY_PATH
#srun --mpi=pmi2 vasp_ncl
#sleep 10000
mpirun /home/users/shenyc/vasp5.4.4/wannier90-1.2/wannier90.x ./wannier90

