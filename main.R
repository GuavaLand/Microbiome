
source('FindSS.R')
#source('ML_DensityMediated_RandomForest.R')

N <-  10
numTrial = 1000
HOI_SS10 <-  matrix(nrow = numTrial, ncol = N)
GLV_SS10 <-  matrix(nrow = numTrial, ncol = N)

for (trial in 1:numTrial) {
  source('DensityMediated_varyingCoefficientLinear.R')
  #Modified GLV
  for (species in 2:ncol(dat)) {
    HOI_SS10[trial,species-1] <-  FindSS(dat[,species])
  }
  
  #Original GLV
  for (oSpecies in 2:ncol(dat1)) {
    GLV_SS10[trial,oSpecies-1] <-  FindSS(dat1[,oSpecies])
  }
}



N <-  20
HOI_SS20 <-  matrix(nrow = numTrial, ncol = N)
GLV_SS20 <-  matrix(nrow = numTrial, ncol = N)

for (trial in 1:numTrial) {
  source('DensityMediated_varyingCoefficientLinear.R')
  #Modified GLV
  for (species in 2:ncol(dat)) {
    HOI_SS20[trial,species-1] <-  FindSS(dat[,species])
  }
  
  #Original GLV
  for (oSpecies in 2:ncol(dat1)) {
    GLV_SS20[trial,oSpecies-1] <-  FindSS(dat1[,oSpecies])
  }
}


N <-  30
HOI_SS30 <-  matrix(nrow = numTrial, ncol = N)
GLV_SS30 <-  matrix(nrow = numTrial, ncol = N)

for (trial in 1:numTrial) {
  source('DensityMediated_varyingCoefficientLinear.R')
  #Modified GLV
  for (species in 2:ncol(dat)) {
    HOI_SS30[trial,species-1] <-  FindSS(dat[,species])
  }
  
  #Original GLV
  for (oSpecies in 2:ncol(dat1)) {
    GLV_SS30[trial,oSpecies-1] <-  FindSS(dat1[,oSpecies])
  }
}