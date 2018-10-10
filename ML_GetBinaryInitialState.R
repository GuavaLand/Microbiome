################################################
#Input: number of species;
####how many percent of binary data to be in mask1 (for training)
#Output: mask1 (*initial density for training) and mask2 (*init density for testing)
###############################################

getBinaryInitialState <- function(N,percent){
  ###########   Approach 1 for generating data points  #######
  #Generate boolean mask of 2^N x N matrix
  repeatBinaryNTimes <- rep(list(c(0,1)),N)
  mask <- expand.grid(repeatBinaryNTimes)
  #Randomize rows
  mask = mask[sample(nrow(mask)),]
  
  xPercentRow = floor(nrow(mask)*percent)
  
  mask1 = mask[1:xPercentRow,]
  mask2 = mask[(xPercentRow+1):nrow(mask),]
  ###########################################################
  
  
  ###########   Approach 2 for generating data points  #######
  ##rows of mask: # of data points
  #M = 1000
  #mask <- matrix(sample(c(0,1),M*N,replace=TRUE), nrow = M, ncol = N)
  
  ##count unique number of rows
  ##cus the way we generate mask may create repeated rows
  #nrow(count(mask))
  ###########################################################
  re = list(mask1 = mask1, mask2 = mask2)
  return(re)
}

