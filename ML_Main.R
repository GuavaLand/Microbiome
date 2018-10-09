source('DensityMediated_varyingCoefficientLinear.R')
source('ML_FindSS.R')

#Number of of species
N = 10

##############################
#1. Simulation
#2. Find Steady State Density

simulated_result = growthFunction(N)
dat = simulated_result[[1]]