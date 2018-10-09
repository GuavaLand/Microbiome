library(plyr)

source('ML_GetCommunityParam.R')
source('ML_CommunitySimulation.R')
source('ML_FindSS.R')


#Number of of species
N = 10

##############################
#1. Get simulation parameters: alpha, c0, l, init
#Apply binary mask to init and loop:
  #2. Simulation
  #3. Find Steady State Density

#1. Get simulation parameters
simulationParam <- getCommunityParam(N)
alpha <- simulationParam$alpha
c0 <- simulationParam$c0
l <- simulationParam$l
init <- simulationParam$init

#rows of mask: # of data points
M = 2
mask <- matrix(sample(c(0,1),M*N,replace=TRUE), nrow = M, ncol = N)

#count unique number of rows
#cus the way we generate mask may create repeated rows
nrow(count(mask))

init_mask = t(t(mask) * init)

#2. Do simulation and save result to dat_list
dat_list <- apply(init_mask, 1, function(x){growthFunction(N,alpha,c0,l,x)})

#3. Find steady state
#SS <- matrix()
#SS <- apply(dat[,2:(N+1)],MARGIN = 2,findSS)


