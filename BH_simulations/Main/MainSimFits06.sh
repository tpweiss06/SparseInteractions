#!/bin/bash -l

#SBATCH --account=commbayes
#SBATCH --time=144:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=3
#SBATCH --cpus-per-task=1
#SBATCH --mail-type=ALL
#SBATCH --mail-user=cweissle@uwyo.edu
#SBATCH --job-name=SparseInteractionsTrial
#SBATCH --mem-per-cpu=128000

# Set the parameter combination to use and generate names of R scripts and log files
Rscript=MainSimFits06.R
LogFile=MainSimFits06.log

# Change to the relevant working directory
cd /project/commbayes/SparseInteractions/BH_sims/

# Load R and MPI
module load gcc/7.3.0 r/3.5.3 swset/2018.05  gcc/7.3.0 r-rstan/2.18.2-py27

R < $Rscript > $LogFile --no-save 
