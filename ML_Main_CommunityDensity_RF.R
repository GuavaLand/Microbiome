library(plyr)
library(randomForest)
library(ggplot2)
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

CommunityDensity = as.matrix(apply(SS_df,1,sum))
colnames(CommunityDensity) = c('CommunityDensity')

all_data = init_mask_df
all_data$CommunityDensity = CommunityDensity

model <- function(train_size){ #make model training reusable with a single set of simulation (all_data) 
  #Decide how many sample to train model
  train_size = train_size
  train_data = all_data[1:train_size,]
  
  rf = randomForest(CommunityDensity ~ ., data = train_data)
  
  ###########################################
  #1. take test_sample_size samples from all data
  #2. use rf to predict final state
  #3. calculate prediction accuracy
  #4. save as a csv to cwd
  
  test_sample_size = 50
  
  #1. take test sample
  actual_sample  = all_data[(train_size+1):(train_size + test_sample_size),]
  
  #2. use rf_species to predict final state
  predicted_SS = matrix(nrow = nrow(actual_sample), ncol = 1)
  colnames(predicted_SS) = c('Predicted_Community_Density')
  
  for (row in 1:nrow(actual_sample)) {
    data_point = as.data.frame(actual_sample[row,1:N])
    predicted_SS[row,1] = predict(rf,data_point)
  }
  
  #3. calculate accuracy
  difference_score = colSums((predicted_SS - actual_sample$CommunityDensity)^2)/nrow(actual_sample)
    #Another way of expressing accuracy: mean of prediction error%
  prediction_accuracy_mean = mean(abs(predicted_SS - actual_sample$CommunityDensity)/actual_sample$CommunityDensity)
  prediction_accuracy_sd = sd(abs(predicted_SS - actual_sample$CommunityDensity)/actual_sample$CommunityDensity)
  
  re = list(difference_score=difference_score,
            prediction_accuracy_mean=prediction_accuracy_mean,
            prediction_accuracy_sd=prediction_accuracy_sd,
            actual_sample=actual_sample,
            predicted_SS=predicted_SS)
  
  return(re)
}


#Create matrix to store training sample size and the result difference score and accuracy measure
differenceScoreVsSampleSize = matrix(nrow = 18, ncol = 4)
colnames(differenceScoreVsSampleSize) = c('difference_score',
                                          'prediction_accuracy_mean',
                                          'prediction_accuracy_sd',
                                          'sample_size')

par(mfrow = c(3,3))

row_counter = 0
for (counter in c(c(1:9),seq(10,90,10))) {
  sample_size = counter *10
  returned = model(sample_size)
  #Save all the returned matrices
  difference_score = returned$difference_score
  prediction_accuracy_mean = returned$prediction_accuracy_mean
  prediction_accuracy_sd = returned$prediction_accuracy_sd
  actual_sample = returned$actual_sample
  predicted_SS = returned$predicted_SS
  
  row_counter = row_counter + 1
  differenceScoreVsSampleSize[row_counter,1] = difference_score
  differenceScoreVsSampleSize[row_counter,2] = prediction_accuracy_mean
  differenceScoreVsSampleSize[row_counter,3] = prediction_accuracy_sd
  differenceScoreVsSampleSize[row_counter,4] = sample_size
  
  plot(actual_sample$CommunityDensity,predicted_SS, main = paste(sample_size,'Samples Training:',round(difference_score,4)),
       xlab = 'Actual Community Density',ylab = 'Predicted Community Density')
}

#Plot our defined difference score as a function of training sample size
plot(differenceScoreVsSampleSize[,4],differenceScoreVsSampleSize[,1], main = 'Difference Score of Predicted and Actual \nCommunity Density over Sample Size',
     xlab = 'Training Sample Size', ylab = 'Difference Score of Predicted and Actual')


#Plot our accuracy measure as a function of training sample size (scatter plot with error bar)
differenceScoreVsSampleSize = as.data.frame(differenceScoreVsSampleSize)
p <- ggplot(differenceScoreVsSampleSize, aes(x=sample_size, y=prediction_accuracy_mean)) + 
  geom_line()+
  geom_point()+
  geom_errorbar(aes(ymin=prediction_accuracy_mean-prediction_accuracy_sd, ymax=prediction_accuracy_mean+prediction_accuracy_sd))
#title, labels etc.
p+labs(title="RF Prediction Accuracy as % Difference btw Predicted and Actual", x="Training Sample Size", y = "Mean of Prediction Error%")+
  scale_y_continuous(breaks=seq(-0.5,0.8,0.1))

#write.table(differenceScoreVsSampleSize, paste(getwd(),'/score102.csv', sep = ''), sep="\t")
