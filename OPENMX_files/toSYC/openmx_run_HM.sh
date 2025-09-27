#!/bin/bash -l
#
#SBATCH -N 1
#SBATCH -p regular
#SBATCH --ntasks-per-node=52
#SBATCH --job-name=openmx_test
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

mpirun /home/users/shenyc/openmx/openmx3.9/source/openmx openmx.dat


# Band Structure Calculation
mkdir -p BD
cp -r openmx.scfout BD
# analysis_example openmx.scfout > HS.out
cp openmx.dat BD/MABANDS.dat
cd BD
cat >> MABANDS.dat << EOF
scf.restart    on

Band.dispersion              on
Band.Nkpath                8
<Band.kpath
50   0.0000000000   0.0000000000   0.0000000000      0.5000000000   0.0000000000   0.5000000000     G   X
50   0.5000000000   0.0000000000   0.5000000000      0.5000000000   0.2500000000   0.7500000000     X   W
50   0.5000000000   0.2500000000   0.7500000000      0.3750000000   0.3750000000   0.7500000000     W   K 
50   0.3750000000   0.3750000000   0.7500000000      0.0000000000   0.0000000000   0.0000000000     K   G
50   0.0000000000   0.0000000000   0.0000000000      0.5000000000   0.5000000000   0.5000000000     G   L
50   0.5000000000   0.5000000000   0.5000000000      0.6250000000   0.2500000000   0.6250000000     L   U          
50   0.6250000000   0.2500000000   0.6250000000      0.5000000000   0.2500000000   0.7500000000     U   W 
50   0.5000000000   0.2500000000   0.7500000000      0.5000000000   0.5000000000   0.5000000000     W   L
Band.kpath>
EOF
mpirun /home/users/shenyc/openmx/openmx3.9/source/openmx MABANDS.dat 
rm -rf *cube*

echo Job ended at `date`


