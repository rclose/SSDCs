### Run simulations testing performance of richness estimators for MEE MS ###

source("./functions/richness_estimator_functions.R")

source("./subscripts/11.1 Simulation Experiment 1 - systematically varying richness and evenness.R", echo = T)

source("./subscripts/11.2 Simulation Experiment 2 - standardising estimators by coverage.R", echo = T)

source("./subscripts/11.3 Plot results of Simulation Experiment 1.R", echo = T)

source("./subscripts/11.4 Simulation Experiment 3 - randomly varying evenness and richness.R", echo = T)

source("./subscripts/11.5 Simulations 4 - Testing Good's u in simulations.R", echo = T)

save.image(paste("./saved_workspaces/saved_workspace_richness_simulation.RData", sep = ""))
# load(paste("./saved_workspaces/saved_workspace_richness_simulation.RData", sep = ""))
