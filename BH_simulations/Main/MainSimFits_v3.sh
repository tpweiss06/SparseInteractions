#!/bin/bash -l

#SBATCH --account=commbayes
#SBATCH --time=72:00:00
#SBATCH --nodes=20
#SBATCH --ntasks-per-node=32
#SBATCH --cpus-per-task=1
#SBATCH --mail-type=ALL
#SBATCH --mail-user=cweissle@uwyo.edu
#SBATCH --job-name=SparseInteractions_Rmpi

# Set the parameter combination to use and generate names of R scripts and log files
Rscript=MainSimFits_v3.R
LogFile=MainSimFits_v3.log

# Change to the relevant working directory
cd /project/commbayes/SparseInteractions/BH_sims/

# Load R and MPI
module load gcc/7.3.0 r/3.5.3 swset/2018.05 r-rstan/2.18.2-py27 openmpi/3.1.0 r-rmpi/0.6-9-r353-py27

mpirun -np 1 R CMD BATCH --no-restore --no-save --quiet  $Rscript $LogFile

