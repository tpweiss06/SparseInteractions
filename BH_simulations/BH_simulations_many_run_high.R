
# Set the number of processors and number of simulations to be run
n_processors <- 32 # number of nodes
n_sims <- 1200


# Set the working directory and load necessary data and libraries
setwd("/project/commbayes/SparseInteractions/BH_sims/")
library(parallel)
library(Rmpi)

# Set parameter values for the simulations
#

#data <- read.csv("MyDataFile.csv")

# Write a function to pass to to the worker nodes. The function should perform
#    a single simulation, and will then be repeated enough times to produce
#    the desired number of simulations.

# CW: not sure if this is necessary if I'm just calling an existing function?
Sim_function <- function(i){
  output <- Run.simulation(a.range = 4.5)
  return(output)
}

# Create the cluster and run the simulations
cl <- makeCluster(n_processors - 1, type = "MPI")

# Export any necessary R objects from the preprocessing stage to the indivdidual
#    processors (e.g. Do the processors need to know the value of any of the parameters
#    set up before?). The names of the objects need to be stored together as character
#    vectors because you are basically just telling R what to look for in the environment.

# ObjectsToImport <- c("R", "alpha", "k", "maxt", "Nt_init")
# clusterExport(cl, ObjectsToImport) # values and vectors not functions

# Run any commands necessary on the processors before running the simulation. This 
#    is where you can set the working directory of the processors, source any
#    necessary files, or load any necessary libraries
clusterEvalQ(cl, setwd("/project/commbayes/SparseInteractions/BH_sims/"))
clusterEvalQ(cl, library(dplyr))
clusterEvalQ(cl, library(tidyr))
clusterEvalQ(cl, library(tibble))
clusterEvalQ(cl, source("BH_simulations_many_functions.R")) # put the functions used on each node

# Run the simulations on the cluster. Here x corresponds to the first argument
#    of the SimFunc function. The cluster will evaluate the SimFunc for each element
#    of the vector x below. I usually just use a vector of 1 to the number of simulations,
#    but you could get creative with this too, if you wanted. The output of the SimFunc
#    function is stored in a list which you can immediately save, or process within this
#    same script.
simulations <- clusterApply(cl, x = 1:n_sims, fun = Sim_function) # run it x amount of times


# Save as an r dataframe
save(simulations, file=('BH_simulations_high_1200.RData'))


