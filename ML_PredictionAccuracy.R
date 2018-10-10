#################################################
#Input: predicted and actual density of species (vector)
#Output: prediction accuracy as a percentage (Euclidean distance percentage)
#################################################
predictionAccuracy <- function(predicted, actual){
  #difference
  difference = predicted - actual
  difference = difference^2
  return(difference)
}