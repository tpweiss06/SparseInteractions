#!/bin/bash -l

#SBATCH --account=commbayes
#SBATCH --time=144:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=3
#SBATCH --cpus-per-task=1
#SBATCH --mail-type=ALL
#SBATCH --mail-user=cweissle@uwyo.edu
#SBATCH --job-name=SparseInteractionsMainSim
#SBATCH --mem-per-cpu=1536M

# Set the parameter combination to use and generate names of R scripts and log files
Rscript=MainSimFits.R
LogFile=MainSimFits.log

# Change to the relevant working directory
cd /project/commbayes/SparseInteractions/BH_sims/

# Load R and MPI
module load gcc/7.3.0 r/3.5.3 swset/2018.05  gcc/7.3.0 r-rstan/2.18.2-py27

mpirun -np 1 R CMD BATCH --no-restore --no-save --quiet  $Rscript $LogFile
