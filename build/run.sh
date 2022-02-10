#!/bin/bash
#SBATCH -p amd_256
#SBATCH -N 1
#SBATCH -n 64
source /public3/soft/module/module.sh
module load gcc/9.3.0-lcc-public3 mpi/intel/17.0.5-cjj-public3 cmake/3.17.0-public3
export PATH=/public3/home/sc72221/rly/libtorch1.5.0/bin$PATH
export LD_LIBRARY_PATH=/public3/home/sc72221/rly/libtorch1.5.0/lib:$LD_LIBRARY_PATH
export INCLUDE=/public3/home/sc72221/rly/libtorch1.5.0/include:$INCLUDE
srun -n 40 /public3/home/sc72221/rly/upload/build/MPF model_1026_3.pt
