#!/bin/bash

# use sbatch to lauch script

BASE_DIR=/scratch/VR0052/ESA2014/surf/experiments

#SBATCH --job-name=nodes_evaluated.sh
#SBATCH --account="VR0052"
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=10-12:00:00
#SBATCH --mem=240GB
#SBATCH --mail-user simon.gog@unimelb.edu.au
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --reservation=VR0052_14

module load gcc
module load cmake

$BASE_DIR/nodes_evaluated.sh

## --exclusive
