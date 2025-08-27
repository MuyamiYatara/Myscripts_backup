#!/bin/bash
#SBATCH -N 1 
#SBATCH -n 56
###SBATCH --exclusive
###SBATCH --gres=gpu:0
###SBATCH --nodelist=node017
#SBATCH --partition=regular
###SBATCH --ntasks-per-node=50
#SBATCH --job-name=openmx
#SBATCH -A hmt03
#SBATCH --output=./log
#SBATCH --error=./err

#加载执行任务需要的模块
#source /.bashrc
module load openmx3.9
module load cuda11.8
module load oneapi22.3
module load mkl/mkl2022.2.1
module load nvhpc/22.11
###module openmx3.9

export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export MV2_ENABLE_AFFINITY=0
echo "The current job ID is $SLURM_JOB_ID"
echo "Running on $SLURM_JOB_NUM_NODES nodes:"
echo $SLURM_JOB_NODELIST
echo "Using $SLURM_NTASKS_PER_NODE tasks per node"
echo "A total of $SLURM_NTASKS tasks is used"
ulimit -s unlimited

#export PATH=/home/apps/openmx/openmx3.9/work:$PATH
#ncpu=$(($SLURM_NTASKS_PER_NODE*$SLURM_JOB_NUM_NODES))
#echo $ncpu
NP=$SLURM_NTASKS
echo "np is $NP"
# 输出路径下的内容
cd /data/home/apps/openmx/openmx3.9/DFT_DATA19
echo "The current path is `pwd`"
cd $SLURM_SUBMIT_DIR
mpirun -np $NP openmx openmx_in.dat 
mkdir -p HS_forValleyProj
cp -r openmx.scfout HS_forValleyProj/
analysis_example openmx.scfout > HS.out
cp openmx_in.dat MABANDS.dat
cat >> MABANDS.dat << EOF
scf.restart    on

Band.dispersion              on
Band.Nkpath                3
<Band.kpath
15   0.0000 0.0000 0.0000       0.5000 0.0000 0.0000       G  M
15  0.5000 0.0000 0.0000       0.3333333 0.3333333 0.0000   M  K
15   0.3333333 0.3333333 0.0000       0.0000 0.0000 0.0000   K  G
Band.kpath>
EOF
mpirun -np $NP openmx MABANDS.dat 
rm -rf *cube*

echo Job ended at `date`
