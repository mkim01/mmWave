#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --time=5:00:00
#SBATCH --mem=32GB
#SBATCH --job-name=myTest
#SBATCH --mail-type=END
#SBATCH --mail-user=mk7291@nyu.edu
#SBATCH --output=mk_%j.out
  
module purge
module load matlab/2017b

RUNDIR=$SCRATCH/mmWave/Simulations
cd $RUNDIR

matlab -nodisplay -nodesktop -r "run SimulationLOS, exit"
