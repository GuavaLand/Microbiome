########################################################
#Assuming that the naive model generates predicted data
#HOI is the observed data
#Run regression on predicted data vs. observed data
########################################################

########################################################
#Naive model
########################################################
#Define Naive model
nm <- function(t, x, params){
  with(as.list(params, c(x)),{
    dx <- x + x*(c0%*%x)
    list(dx)
  })
}

#Define integration method
n.integrate <- function(time, init.x, model, params){
  as.data.frame(lsoda(init.x, time, model, params))
}

naiveModel <- function(N, init.x, c0){
  ##Define intial density
  #init.x <- floor(runif(N)*10)/10
  
  ##Define c0 as 1/N
  #c0 <- matrix(rep(-1/N, N^2), nrow = N)
  ##Define intra-species effect
  #diag(c0) <- -0.5
  
  #Solve ode
  dat <- n.integrate(1:500, init.x, nm, list(c0 = c0))
  
  #plot
  matplot(x = dat$time, y = dat[,-1], type = 'b', xlab = 'time', ylab = 'Absolute abundance', main = paste('Naive Model',N,'Species'))
  
  return(dat)
}

########################################################
#HOI model
########################################################
#Define GLV with varying coefficient
glv <- function(t, x, params){
  with(as.list(params, c(x)),{
    dx <- alpha*x + x*as.vector(c0%*%x)+x*t(t(x)%*%(do.call(cbind, lapply(l, FUN=function(ma) ma%*%x))))
    list(dx)
  })
}

growthFunction <- function(N, init.x, c0){
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
  
  #Solve the ode
  dat <- n.integrate(0:500, init.x, glv, list(alpha=alpha, c0=c0, l=l))
  
  #Plot
  matplot(x=dat$time, y=dat[,-1], typ='b', xlab='time', ylab='Absolute abundance', main=paste('Modified GLV-density', N,'species'))
  
  return(dat)
}



#Find species steady state density
SSDensity <- function(vec){
  for (i in 1:(length(vec) - 10)) {
    #From start, search for 3 points each 5 unit time apart where delta < 10^-3
    thisSlope <- abs(vec[i+1] - vec[i])
    nextSlope <- abs((vec[i+5] - vec[i])/5)
    lastSlope <- abs(vec[i+10] - vec[i+9])
    
    if (thisSlope < 10^-3 & nextSlope < 10^-3 & lastSlope < 10^-3) {
      #Verification: make sure no change after plateau: extrapolate
      #Extrapolate by thisSlope and calculate value at end of integration time
      ExtrEndTimeValue = vec[i] + nextSlope*(length(vec) - i)
      if ((ExtrEndTimeValue < (vec[length(vec)]) + 0.01)& (ExtrEndTimeValue > (vec[length(vec)]) - 0.01)) {
        return(vec[i])
      }
      #There is a plateau in the middle, but in the end density changes again: find next plateau and verify agin
    }
  }
  #Til the end no steady state found
  return(NaN)
}



#######################################################
#Common parameter for Naive model and HOI model
#######################################################
#Define intial density
init.x <- floor(runif(N)*10)/10

#Define c0 as 1/N
c0 <- matrix(rep(-1/N, N^2), nrow = N)
#Define intra-species effect
diag(c0) <- -0.5

#Run both models for community size of 10 to 50
for (n in 1:5) {
  #Community size
  N <- n * 10
  
  predicted_dat <- naiveModel(N, init.x, c0)
  observed_dat <- growthFunction(N, init.x, c0)
  
  #Now find steady state density of each species in predicted_dat and observed_dat
  for (i in 2:ncol(predicted_dat)) {
    vec <- predicted_dat[,i]
    SSDensity[vec]
  }
  
  for (j in 2:ncol(observed_dat)) {
    ovec <- observed_dat[,j]
    SSDensity[ovec]
  }
  
  
}