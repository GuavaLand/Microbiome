library(plyr)
library(randomForest)
source('ML_GetCommunityParam.R')
source('ML_CommunitySimulation.R')
source('ML_FindSSDensity.R')
source('ML_PredictionAccuracy.R')
source('ML_GetBinaryInitialState.R')


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

#mask as presence/absence of each species
mask <- getBinaryInitialState(N, 0.8)
mask1 = mask$mask1
mask2 = mask$mask2


init_mask = t(t(mask1) * init)
init_mask_2 = t(t(mask2) * init)

#2. Do simulation and save result to dat_list
dat_list <- apply(init_mask, 1, function(x){growthFunction(N,alpha,c0,l,x)})

#3. Find steady state
matrixToSS <- function(densityMatrix){
  s <- apply(densityMatrix[,2:(N+1)], MARGIN = 2, findSSDensity)
  return(s)
}

#output list: each member is SS density of every species
SS <- lapply(dat_list[c(1:325,327:819)], matrixToSS)
#rbind all members in the list to form matrix
SS <- do.call(rbind, SS)
#If SS density is too small, set to 0
SS[which(SS < 0.001)] = 0


############################################
# input: init_mask
# ML
# output: SS
#
#each row is an entry of data (determined by M)
#each col is a species
#data are density
############################################

########    consider modulate prediction part  ##########

train_x = as.data.frame(init_mask[c(1:325,327:819),])
train_y = as.data.frame(SS)

if (any(is.na(train_y))) {
  #Row index where the row contains NA
  NA_containing_rows = unique(which(is.na(train_y),arr.ind = TRUE)[,1])
  #Remove NA containing rows
  train_x = train_x[-NA_containing_rows,]
  train_y = train_y[-NA_containing_rows,]
}

#Train ONE rf of each column in train_y ~ train_x
rf_species = list()
for (i in 1:N) {
  #for species i, if train_x starts with 0, train_y will be 0 
  #exclude 0
  null_begining = which(train_x[,i] == 0)
  train_x_filtered = train_x[-null_begining,]
  train_y_filtered = train_y[-null_begining,]
  
  data = cbind(train_x_filtered,y = train_y_filtered[,i])
  rf_species[[i]] = randomForest(y ~ ., data = data)
}

###########################################
#1. take test_sample_size samples from mask 2
#2. simulate to get 'real' final state
#3. use rf_species to predict final state
#4. calculate prediction accuracy
#5. save as a csv to wd

test_sample_size = 50

#1. take test sample
actual_sample  = init_mask_2[50:(49 +test_sample_size),]

#2. get real final state
actual_dat_list = apply(actual_sample, 1, function(x){growthFunction(N,alpha,c0,l,x)})

#output list: each member is SS density of every species
actual_SS <- lapply(actual_dat_list, matrixToSS)
#rbind all members in the list to form matrix
actual_SS <- do.call(rbind, actual_SS)
#If SS density is too small, set to 0
actual_SS[which(actual_SS < 0.001)] = 0

#3. use rf_species to predict final state
predicted_SS = matrix(nrow = test_sample_size, ncol = N)

for (row in 1:test_sample_size) {
  data_point = t(as.data.frame(actual_sample[row,]))
  for (modID in 1:length(rf_species)) {
    predicted_SS[row,modID] = predict(rf_species[[modID]],data_point)
  }
}

#4. calculate accuracy
actual_SS[is.na(actual_SS)] <- 0
difference_score = sum(predicted_SS - actual_SS)^2
