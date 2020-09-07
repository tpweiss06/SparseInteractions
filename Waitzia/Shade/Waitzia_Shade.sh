#!/bin/bash -l

#SBATCH --account=commbayes
#SBATCH --time=48:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=6
#SBATCH --cpus-per-task=1
#SBATCH --mail-type=ALL
#SBATCH --mail-user=cweissle@uwyo.edu
#SBATCH --job-name=WaitziaShade

# Set the parameter combination to use and generate names of R scripts and log files
Rscript=Waitzia_Shade.R
LogFile=Waitzia_Shade.log

# Change to the relevant working directory
cd /project/commbayes/SparseInteractions

# Load R and MPI
module load gcc/7.3.0 r/3.5.3 swset/2018.05 r-rstan/2.18.2-py27

R < $Rscript > $LogFile --no-save 
