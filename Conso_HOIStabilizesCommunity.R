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
#N <- 30

growthFunction <- function(N){
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
  dat <- n.integrate(0:500, init.x, glv, list(alpha=alpha, c0=c0, l=l))
  dat1 <- n.integrate(0:500, init.x, glv1, list(alpha=alpha, c0=c0))
  
  #Plot
  matplot(x=dat$time, y=dat[,-1], typ='b', xlab='time', ylab='Absolute abundance', main=paste('Modified GLV-density', N,'species'))
  matplot(x=dat1$time, y=dat1[,-1], typ='b', xlab='time', ylab='Absolute abundance', main=paste('Original GLV-density', N,'species'))
  
  returnList <- list(dat,dat1)
  
  return(returnList)
}

######################################################################################
#Given density data (vector) over time (1:length(vector)), find steady state time
#If there is no steady state, return NaN
######################################################################################

FindSS <- function(vec){
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
        return(i)
      }
      #There is a plateau in the middle, but in the end density changes again: find next plateau and verify agin
      
    }
  }
  #Til the end no steady state found
  return(NaN)
}

###############################################################
#Simulate species growth function given number of species N
#Find the steady state for each species and store in matrix HOI_SSN or GLV_SSN, based on model used for simulation
#Save HOI_SSN and GLV_SSN as csv at working directory
###############################################################
numTrial = 1000

HOI_list <- list()
GLV_list <- list()

for (n in 1:2) {
  #n is a counter. N is the number of species in community
  N <- n*5
  
  HOI_SSN <-  matrix(nrow = numTrial, ncol = N)
  GLV_SSN <-  matrix(nrow = numTrial, ncol = N)
  
  for (trial in 1:numTrial) {
    repeat {
      #densityData is a list of dat (HOI model) and dat1 (GLV)
      densityData <- growthFunction(N)
      #Check if any NaN in dat
      dat <- densityData[[1]]
      datNA <- any(is.na(dat))
      #Check if any NaN in dat1
      dat1 <- densityData[[2]]
      dat1NA <- any(is.na(dat1))
      #If any of datNA and dat1NA is true, growth function should run again
      if(!(datNA|dat1NA)){
        break
      }
    }
    
    #Modified GLV
    for (species in 2:ncol(dat)) {
      HOI_SSN[trial,species-1] <-  FindSS(dat[,species])
    }
    
    #Original GLV
    for (oSpecies in 2:ncol(dat1)) {
      GLV_SSN[trial,oSpecies-1] <-  FindSS(dat1[,oSpecies])
    }
  }
  
  HOI_list[[n]] <- HOI_SSN
  GLV_list[[n]] <- GLV_SSN
}

###########################################################
#Find community steady state time and plot histogram
###########################################################

CommunitySteadyStateDF <- data.frame()

#Loop through HOI_SS10, HOI_SS20...GLV_SS10, GLV_SS20...
for (n in 1:length(HOI_list)) {
  
  HOI = HOI_list[[n]]
  GLV = GLV_list[[n]]
  
  numOfEntries = nrow(HOI)
  
  #CommunitySteadyStateTime to store the time when the whole community has reached steady state
  CommunitySteadyStateTime <- matrix(nrow = numOfEntries, ncol = 3)
  colnames(CommunitySteadyStateTime) <- c('CommunitySize', 'SteadyStateTimeHOI', 'SteadyStateTimeGLV')
  
  for (i in 1:numOfEntries) {
    HOINA <- any(is.na(HOI[i,]))
    GLVNA <- any(is.na(GLV[i,]))
    #If there is NA in a particular row, for wither HOI or GLV, this community combination did not reach SS within integration time. Move to next row
    if (!HOINA & !GLVNA) {
      #If no NA in both HOI and GLV (for that row, both model reached steady state)
      CommunitySteadyStateTime[i,1] <- n * 10
      CommunitySteadyStateTime[i,2] <- max(HOI[i,2:ncol(HOI)])
      CommunitySteadyStateTime[i,3] <- max(GLV[i,2:ncol(GLV)])
    }
    else{
      CommunitySteadyStateTime[i,1] <- n * 10
      CommunitySteadyStateTime[i,2] <- NA
      CommunitySteadyStateTime[i,3] <- NA
    }
  }
  #Finished looping through an HOI_SSN and GLV_SSN, store the CommunitySteadyStateTime vector to a list
  CommunitySteadyStateDF <-  rbind(CommunitySteadyStateDF,as.data.frame(CommunitySteadyStateTime))
}

#Calculate the difference of community SS time between GLV model and and HOI model 
CommunitySteadyStateDF$Difference <- CommunitySteadyStateDF$SteadyStateTimeGLV - CommunitySteadyStateDF$SteadyStateTimeHOI

#plot
par(mfrow = c(3,2))
#plot histogram for Community Steady State time difference for each of the community size 
hist(CommunitySteadyStateDF[CommunitySteadyStateDF$CommunitySize == 10,4], main = 'How Much Longer Does GLV Reaches SS Than HOI 10 Species', xlab = 'Time Difference', ylim = c(0,500), breaks = 18)
hist(CommunitySteadyStateDF[CommunitySteadyStateDF$CommunitySize == 20,4], main = 'How Much Longer Does GLV Reaches SS Than HOI 20 Species', xlab = 'Time Difference', ylim = c(0,500), breaks = 18)
#hist(CommunitySteadyStateDF[CommunitySteadyStateDF$CommunitySize == 30,4], main = 'How Much Longer Does GLV Reaches SS Than HOI 30 Species', xlab = 'Time Difference', ylim = c(0,500), breaks = 18)
#hist(CommunitySteadyStateDF[CommunitySteadyStateDF$CommunitySize == 40,4], main = 'How Much Longer Does GLV Reaches SS Than HOI 40 Species', xlab = 'Time Difference', ylim = c(0,500), breaks = 18)
#hist(CommunitySteadyStateDF[CommunitySteadyStateDF$CommunitySize == 50,4], main = 'How Much Longer Does GLV Reaches SS Than HOI 50 Species', xlab = 'Time Difference', ylim = c(0,500), breaks = 18)
#plot boxplotboxplot(SteadyStateTime~CommunitySize ,data = CommunitySteadyStateDF)Differences
boxplot(Difference~CommunitySize ,data = CommunitySteadyStateDF, main = 'How Much Longer Does GLV Reaches SS Than HOI', xlab = 'Community Size', ylab = 'Time Difference')
