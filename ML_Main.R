library(plyr)
library(randomForest)
source('ML_GetCommunityParam.R')
source('ML_CommunitySimulation.R')
source('ML_FindSSDensity.R')


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
M = 1000
mask <- matrix(sample(c(0,1),M*N,replace=TRUE), nrow = M, ncol = N)

#count unique number of rows
#cus the way we generate mask may create repeated rows
nrow(count(mask))

init_mask = t(t(mask) * init)

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


############################################
# input: init_mask
# ML
# output: SS
#
#each row is an entry of data (determined by M)
#each col is a species
#data are density
############################################

train_x = as.data.frame(init_mask)
train_y = as.data.frame(SS)

if (any(is.na(train_y))) {
  #Row index where the row contains NA
  NA_containing_rows = unique(which(is.na(train_y),arr.ind = TRUE)[,1])
  #Remove NA containing rows
  train_x = train_x[-NA_containing_rows,]
  train_y = train_y[-NA_containing_rows,]
}

#Train ONE rf train_x to each column in train_y
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

