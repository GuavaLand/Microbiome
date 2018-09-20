library(randomForest)

#Use Random Forest to make prediction
initialDensityAll = t(t(mask) * init.x)
finalDensityAll = SSMatrix

#Go through each species and train model
par(mfrow=c(2,2))
for (species in 1:N) {
  #For species i, find 0 initial condition
  dead <- mask[,species] == 0
  y <- finalDensityAll[!dead,species]
  x <- mask[!dead,]
  #Break down into training set and testing set
  dividerIndex <- floor(length(y)*0.8)
  finalIndex <- length(y)
  train_y <- y[1:dividerIndex]
  train_x <- x[1:dividerIndex,]
  test_y <- y[(dividerIndex + 1):finalIndex]
  test_x <- x[(dividerIndex + 1):finalIndex,]
  
  #Train model for species i
  mod <- randomForest(train_y ~., train_x)
  
  #Make prediction
  predicted_y = predict(mod, test_x)
  
  #Evaluate model by ploting
  plot(predicted_y, test_y, main = species)
  
  
}

