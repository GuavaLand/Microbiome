source('Noise_GetCommunityParam.R')
source('Noise_CommunitySimulation.R')
source('ML_FindSSDensity.R')
source('ML_PredictionAccuracy.R')
source('ML_GetBinaryInitialState.R')


#Number of of species
N = 6

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

#mask as presence/absence of each species
mask <- getBinaryInitialState(N, 2^N)
#mask1 = mask$mask1 #for training
#mask2 = mask$mask2 #for testing


init_mask = t(t(mask) * init)
#init_mask_2 = t(t(mask2) * init)

#2. Do simulation and save result to dat_list
dat_list <- apply(init_mask, 1, function(x){growthFunction(N,alpha,c0,l,x)})

#3. Find steady state
matrixToSS <- function(densityMatrix){
  s <- apply(densityMatrix[,2:(N+1)], MARGIN = 2, findSSDensity)
  return(s)
}

#output list: each member is SS density of every species
SS <- lapply(dat_list, matrixToSS)
#rbind all members in the list to form matrix
SS <- do.call(rbind, SS)
#If SS density is too small, set to 0
SS[which(SS < 0.001)] = 0

########    consider modulate prediction part  ##########

init_mask_df = as.data.frame(init_mask)
SS_df = as.data.frame(SS)

if (any(is.na(SS_df))) {
  #Row index where the row contains NA
  NA_containing_rows = unique(which(is.na(SS_df),arr.ind = TRUE)[,1])
  #Remove NA containing rows
  init_mask_df = init_mask_df[-NA_containing_rows,]
  SS_df = SS_df[-NA_containing_rows,]
}

CommunityDensity = as.matrix(apply(SS_df,1,sum))
colnames(CommunityDensity) = c('CommunityDensity')

plot(CommunityDensity)
