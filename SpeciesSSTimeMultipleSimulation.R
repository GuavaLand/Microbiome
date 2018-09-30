###############################################################
#Simulate species growth function given number of species N
#Find the steady state for each species and store in matrix HOI_SSN or GLV_SSN, based on model used for simulation
#Save HOI_SSN and GLV_SSN as csv at working directory
###############################################################



source('FindSS.R')
#source('ML_DensityMediated_RandomForest.R')
source('DensityMediated_varyingCoefficientLinear.R')

numTrial = 1000

N <-  10
HOI_SS10 <-  matrix(nrow = numTrial, ncol = N)
GLV_SS10 <-  matrix(nrow = numTrial, ncol = N)

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
    HOI_SS10[trial,species-1] <-  FindSS(dat[,species])
  }
  
  #Original GLV
  for (oSpecies in 2:ncol(dat1)) {
    GLV_SS10[trial,oSpecies-1] <-  FindSS(dat1[,oSpecies])
  }
}
write.csv(GLV_SS10, 'GLV_SS10.csv')
write.csv(HOI_SS10, 'HOI_SS10.csv')


N <-  20
HOI_SS20 <-  matrix(nrow = numTrial, ncol = N)
GLV_SS20 <-  matrix(nrow = numTrial, ncol = N)

for (trial in 1:numTrial) {
  repeat {
    #densityData is a list of dat (modified model) and dat1 (GLV)
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
    HOI_SS20[trial,species-1] <-  FindSS(dat[,species])
  }
  
  #Original GLV
  for (oSpecies in 2:ncol(dat1)) {
    GLV_SS20[trial,oSpecies-1] <-  FindSS(dat1[,oSpecies])
  }
}
write.csv(GLV_SS20, 'GLV_SS20.csv')
write.csv(HOI_SS20, 'HOI_SS20.csv')


N <-  30
HOI_SS30 <-  matrix(nrow = numTrial, ncol = N)
GLV_SS30 <-  matrix(nrow = numTrial, ncol = N)

for (trial in 1:numTrial) {
  repeat {
    #densityData is a list of dat (modified model) and dat1 (GLV)
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
    HOI_SS30[trial,species-1] <-  FindSS(dat[,species])
  }
  
  #Original GLV
  for (oSpecies in 2:ncol(dat1)) {
    GLV_SS30[trial,oSpecies-1] <-  FindSS(dat1[,oSpecies])
  }
}
write.csv(GLV_SS30, 'GLV_SS30.csv')
write.csv(HOI_SS30, 'HOI_SS30.csv')

N <-  40
HOI_SS40 <-  matrix(nrow = numTrial, ncol = N)
GLV_SS40 <-  matrix(nrow = numTrial, ncol = N)

for (trial in 1:numTrial) {
  repeat {
    #densityData is a list of dat (modified model) and dat1 (GLV)
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
    HOI_SS40[trial,species-1] <-  FindSS(dat[,species])
  }
  
  #Original GLV
  for (oSpecies in 2:ncol(dat1)) {
    GLV_SS40[trial,oSpecies-1] <-  FindSS(dat1[,oSpecies])
  }
}
write.csv(GLV_SS40, 'GLV_SS40.csv')
write.csv(HOI_SS40, 'HOI_SS40.csv')


N <-  50
HOI_SS50 <-  matrix(nrow = numTrial, ncol = N)
GLV_SS50 <-  matrix(nrow = numTrial, ncol = N)

for (trial in 1:numTrial) {
  repeat {
    #densityData is a list of dat (modified model) and dat1 (GLV)
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
    HOI_SS50[trial,species-1] <-  FindSS(dat[,species])
  }
  
  #Original GLV
  for (oSpecies in 2:ncol(dat1)) {
    GLV_SS50[trial,oSpecies-1] <-  FindSS(dat1[,oSpecies])
  }
}
write.csv(GLV_SS50, 'GLV_SS50.csv')
write.csv(HOI_SS50, 'HOI_SS50.csv')