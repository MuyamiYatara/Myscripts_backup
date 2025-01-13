#!/bin/bash
#
#SBATCH -N 1
#SBATCH -n 10
#SBATCH -A hmt03  
#SBATCH -p amdgpu40g
#SBATCH --job-name=Fe3GaD03
#SBATCH --output=./log 
#SBATCH --error=./err
#SBATCH --gres=gpu:1
#

#sinfo -n node[005] -o "%C"

ulimit -s unlimited
ulimit -c unlimited
#module purge
#module load oneapi22.3
module load vasp/vasp6.3.0_wannier90-1.2_amd_gpu
module load nvhpc/22.9
#module load mkl/mkl2022.2.1
### 一些提示性输出
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export MV2_ENABLE_AFFINITY=0
echo ”The current job ID is $SLURM_JOB_ID”
echo ”Running on $SLURM_JOB_NUM_NODES nodes:”
echo $SLURM_JOB_NODELIST
echo ”Using $SLURM_NTASKS_PER_NODE tasks per node”
echo ”A total of $SLURM_NTASKS tasks is used”
### 对任务执行的内存不做限制
#ulimit -s unlimited
#ulimit -c unlimited
### 加载任务所需要的库  /home/ycshen/wannier_tools_test/bin/wt.x 
export LD_LIBRARY_PATH=/usr/local/lib64:$LD_LIBRARY_PATH

mpirun -np 1 vasp_ncl




  
