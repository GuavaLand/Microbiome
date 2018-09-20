library(randomForest)

#Use Random Forest to make prediction
initialDensityAll = t(t(mask) * init.x)
finalDensityAll = SSMatrix

#Go through each species and train model 
for (species in 1:N) {
  #For species i, find 0 initial condition
  dead <- mask[,i] == 0
  final <- finalDensityAll[-dead,i]
  data <- mask[-dead,]
  #Break down into training set and testing set
  training <- final[1:floor(length(final)*0.8),]
  rf <- randomForest(final ~., data)
  
}

