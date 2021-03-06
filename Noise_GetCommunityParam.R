getCommunityParam <- function(N){
  #Define species intrinsic growth rate
  alpha <- runif(N,min = 0,max = 0.5)
  
  #Define the constant in species-species interation coefficient
  c0 <- matrix(runif(N*N, min = -0.1, max = 0.1),nrow = N)
  #Set species self interation to -0.5
  for (i in 1:N) {
    for (j in 1:N) {
      if (i == j) {
        c0[i,j] <-  -0.5
      }
    }
  }
  
  #Define coefficient of linear term in species-species interation coefficient
  ck <- matrix(0,nrow = N,ncol = N)
  #ck <- runif(N*N, min = -0.1, max = 0.02)
  #ck <- matrix(ck, nrow = N)
  l <- list()
  for (i in 1:N) {
    #For i-th matrix, i-th row is 0
    temp <- ck
    temp[i,] <- 0
    
    #In i-th matrix, elements are 0 if k (column) == either j (row) or i (matrix order)
    for (j in 1:N) {
      for (k in 1:N) {
        if (k==j | k == i) {
          temp[j,k] <- 0
        }
      }
    }
    
    ##Control the prevalence of thrid party effects
    #for (element in 1:length(temp)) {
    #  dice <- runif(1)
    #  if (dice > 0.1) { #what percent of to assign 0
    #    temp[element] <- 0
    #  }
    #}
    l[[i]] <- temp
  }
  
  #Define initial abundance between 0.1 and 1, to 1 decimal place
  init <- runif(N, min = 0, max = 0.01)
  
  re = list(alpha = alpha, c0 = c0, l = l, init = init)
  return(re)
  
}