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

init_mask_df = as.data.frame(init_mask)
SS_df = as.data.frame(SS)

#Of the input and output, those do not reach SS (contains NA in output) should be removed

if (any(is.na(SS_df))) {
  #Row index where the row contains NA
  NA_containing_rows = unique(which(is.na(SS_df),arr.ind = TRUE)[,1])
  #Remove NA containing rows
  init_mask_df = init_mask_df[-NA_containing_rows,]
  SS_df = SS_df[-NA_containing_rows,]
}

#Decide how many sample to train model
train_size = 500
train_x = init_mask_df[1:train_size,]
train_y = SS_df[1:train_size,]

#Train ONE rf of each column in train_y ~ train_x, store models in rf_species
rf_species = list()
#Track how many samples used to train a model
train_size_deduct_0 = as.data.frame(matrix(nrow = N, ncol = 1))
for (i in 1:N) {
  #for species i, if train_x starts with 0, train_y will be 0 
  #exclude 0
  void_init = which(train_x[,i] == 0)
  train_x_filtered = train_x[-void_init,]
  train_y_filtered = train_y[-void_init,]
  train_size_deduct_0[i,1] = nrow(train_x_filtered)
  
  data = cbind(train_x_filtered,y = train_y_filtered[,i])
  rf_species[[i]] = randomForest(y ~ ., data = data)
}

###########################################
#1. take test_sample_size samples from mask 2
#2. simulate to get 'real' final state
##2.1 find and remove communities that cannot reach SS before use sample for prediction
#3. use rf_species to predict final state
##3.1 before using rf_species, if a species start with 0, the final state is 0.
##that is, only use rf_species when a species has some init state. This is also how we train the model.
#4. calculate prediction accuracy
#5. save as a csv to wd

test_sample_size = 50

#1. take test sample
actual_sample  = init_mask_df[(train_size+1):(train_size + test_sample_size),]

#2. get real final state
actual_dat_list = apply(actual_sample, 1, function(x){growthFunction(N,alpha,c0,l,x)})

#output list: each member is SS density of every species
actual_SS <- lapply(actual_dat_list, matrixToSS)
#rbind all members in the list to form matrix
actual_SS <- do.call(rbind, actual_SS)
#If SS density is too small, set to 0
actual_SS[which(actual_SS < 0.001)] = 0


#2.1 If we see NA in actual_SS, remove the row and remove corresponding row in actual_sample (initial state)
##Because of the assumption that the model is only used to predict community that will reach SS
if (any(is.na(actual_SS))) {
  #Row index where the row contains NA
  NA_containing_rows = unique(which(is.na(actual_SS),arr.ind = TRUE)[,1])
  #Remove NA containing rows
  actual_sample = actual_sample[-NA_containing_rows,]
  actual_SS = actual_SS[-NA_containing_rows,]
}


#3. use rf_species to predict final state
predicted_SS = matrix(nrow = nrow(actual_sample), ncol = N)

for (row in 1:nrow(actual_sample)) {
  data_point = as.data.frame(actual_sample[row,])
  for (modID in 1:length(rf_species)) {
    #3.1 if a species init density is 0, predict its final density to be 0
    if (actual_sample[row,modID] == 0) {
      predicted_SS[row,modID] = 0
    }
    else{
      predicted_SS[row,modID] = predict(rf_species[[modID]],data_point)
    }
  }
}

#4. calculate accuracy
difference_score = colSums((predicted_SS - actual_SS)^2)/nrow(actual_SS)
train_sample_size_difference = cbind(train_size_deduct_0,difference_score)
write.table(train_sample_size_difference, paste(getwd(),'/score',N,'_',train_size,'.csv', sep = ''), sep="\t")

#kk = read.csv('C:\\Source\\Microbiome\\score10.csv',sep = '\t')
