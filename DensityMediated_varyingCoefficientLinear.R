library(deSolve)

#Define GLV with varying coefficient
glv <- function(t, x, params){
  with(as.list(params, c(x)),{
    dx <- alpha*x + x*as.vector(c0%*%x)+x*t(t(x)%*%(do.call(cbind, lapply(l, FUN=function(ma) ma%*%x))))
    list(dx)
  })
}
#glv1 as original GLV
#c0 is beta
glv1 <- function(t, x, params){
  with(as.list(params, c(x)),{
    dx <- alpha*x + x*as.vector(c0%*%x)
    list(dx)
  })
}


#Define integration method
n.integrate <- function(time, init.x, model, params){
  as.data.frame(lsoda(init.x, time, model, params))
}


#Define community size
N <- 30

#Define species intrinsic growth rate
alpha <- runif(N)

#Define the constant in species-species interation coefficient
c0 <- matrix(runif(N*N, min = -1, max = 0),nrow = N)
#Set species self interation to -0.5
for (i in 1:N) {
  for (j in 1:N) {
    if (i == j) {
      c0[i,j] <-  -0.5
    }
  }
}

#Define coefficient of linear term in species-species interation coefficient
ck <- runif(N*N, min = -1, max = 0.2)
ck <- matrix(ck, nrow = N)
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
  
  #Control the prevalence of thrid party effects
  #for (element in 1:length(temp)) {
  #  dice <- runif(1)
  #  if (dice > -1) { #what percent of to assign 0
  #    temp[element] <- 0
  #  }
  #}
  l[[i]] <- temp
}


#Define initial abundance between 0.1 and 1, to 1 decimal place
init.x <- floor(runif(N, min = 1, max = 10))/10

#Solve the ode
dat <- n.integrate(0:100, init.x, glv, list(alpha=alpha, c0=c0, l=l))
dat1 <- n.integrate(0:100, init.x, glv1, list(alpha=alpha, c0=c0))
#Plot

matplot(x=dat$time, y=dat[,-1], typ='b', xlab='time', ylab='Absolute abundance', main=paste('Modified GLV-density', N,'species'))
matplot(x=dat1$time, y=dat1[,-1], typ='b', xlab='time', ylab='Absolute abundance', main=paste('Original GLV-density', N,'species'))



###############################################################################
#Generate 2^N communities where species in each community (max possible N) are present or absent
###############################################################################

#Generate boolean mask of 2^N x N matrix
#repeatBinaryNTimes <- rep(list(c(0,1)),N)
#mask <- expand.grid(repeatBinaryNTimes)
#When N is 50, the mask is too large

#Generate random boolean mask of 100 x N
##############
#Implementation without checking duplicate observations
##############
#mask <- as.data.frame(matrix(sample(c(0,1),N*300, replace = TRUE), ncol = N))

##############
#Implementation checking duplicate observations (slow!)
##############
##First row of mask:
#mask <- as.data.frame(matrix(sample(c(0,1),50, replace = TRUE), ncol = 50))
#for (row in 2:100) {
#  tempNextRow <-  sample(c(0,1),50, replace = TRUE)
#  #Assume this has not appeared before
#  appeared <- FALSE
#  #Run through 1:row-1 to check if appeared. If yes, generate a new one and do again
#  while (TRUE) {
#    for (subrow in 1:(row-1)) {
#      if (all(mask[subrow,] == tempNextRow)) {
#        appeared <- TRUE
#        break
#      }
#    }
#    if (appeared == FALSE) {
#      #Break out of while loop
#      break
#    }
#    tempNextRow <-  sample(c(0,1),50, replace = TRUE)
#  }
#  mask[row,] <- tempNextRow
#}


#Loop through 2^N(or 300) to apply each row in mask to initial abundance, solve ode, and retrieve final abundance
#SSMatrix <- as.data.frame(matrix(nrow = 300, ncol = N))
#colnm <- c(1:N)
#colnames(SSMatrix) <- colnm
#for (i in 1:300) {
#  init <- init.x*mask[i,]
#  init <- as.numeric(init)
#  dat <- n.integrate(0:100, init, glv, list(alpha=alpha, c0=c0, l=l))
#  SSMatrix[i,] <- dat[nrow(dat),2:(N+1)]
#}



