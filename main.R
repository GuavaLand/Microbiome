
source('FindSS.R')
#source('ML_DensityMediated_RandomForest.R')

N <-  10
numTrial = 1000
HOI_SS10 <-  matrix(nrow = numTrial, ncol = N)
GLV_SS10 <-  matrix(nrow = numTrial, ncol = N)

for (i in 1:numTrial) {
  source('DensityMediated_varyingCoefficientLinear.R')
  #Modified GLV
  for (j in 2:ncol(dat)) {
    HOI_SS10[1,j-1] <-  FindSS(dat[,j])
  }
  
  #Original GLV
  for (j in 2:ncol(dat1)) {
    GLV_SS10[1,j-1] <-  FindSS(dat[,j])
  }
}


N <-  20
HOI_SS20 <-  matrix(nrow = numTrial, ncol = N)
GLV_SS20 <-  matrix(nrow = numTrial, ncol = N)

for (i in 1:numTrial) {
  source('DensityMediated_varyingCoefficientLinear.R')
  #Modified GLV
  for (j in 2:ncol(dat)) {
    HOI_SS20[1,j-1] <-  FindSS(dat[,j])
  }
  
  #Original GLV
  for (j in 2:ncol(dat1)) {
    GLV_SS20[1,j-1] <-  FindSS(dat[,j])
  }
}