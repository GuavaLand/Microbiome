library(deSolve)

#Creating random matrix of 12 traits of 100 species
traitMatrix <-  sample(x = c(0,1), 12, replace = TRUE)

for (i in c(2:100)) {
  vec <-  sample(x = c(0,1), 12, replace = TRUE)
  traitMatrix <-  rbind(traitMatrix, vec)
}
rownames(traitMatrix) <- paste("species",c(1:100),sep = "")
colnames(traitMatrix) <- c('MoleSecr1','MoleSecr2','MoleSecr3','MoleUsa1','MoleUsa2','MoleUsa3',
                           'AntiSecr1','AntiSecr2','AntiSecr3','AntiProt1','AntiProt2','AntiProt3')

#Generate relationship matrix of 100 species
a = matrix(rep(0,10000),nrow=100,ncol=100)
for (row in 1:100) {
  for (col in 1:100) {
    if (row == col) {
      a[row,col] = -0.5 # auto-interactions are negative (intra-species competition)
    } else {
      speciesX = traitMatrix[row,]
      speciesY = traitMatrix[col,]
      score <- interationScore(speciesX,speciesY)
      a[row,col] <- score
      a[col,row] <- score
    }
  }
}

#Creating random growth rate of 100 species
b <- runif(100)

#Parameters of Lotka-Volterra equations
parms <- cbind(b,a)
parms=cbind(rep(100,100),parms)

tstart=0                   # time (start)
tend=30                    # time (end) of integration
tstep=0.1                  # time step (resolution)

y<-runif(100)                # initial species abundances 

#Simulation for coupled system
times<-seq(tstart, tend, by=tstep)

out<-lsoda(y, times, glvmat, parms)



#Generalized Lotka-Volterra
glvmat<-function(t, y, parms){
  N=parms[1,1]  # species number
  b=parms[,2]   # vector of growth rates
  a=parms[,3:ncol(parms)] # interaction matrix
  dydt <- y*(b+a %*% y)
  list(dydt)
}

#Interation score of speciesX and speciesY based on traitMatrix
interationScore <- function(speciesX, speciesY){
  score = 0
  #if any of molecule secreted by Y is usable by X, score += 0.12
  if (speciesY['MoleSecr1'] == 1 & speciesX['MoleUsa1'] == 1) {
    score = score + 0.12
  }
  if (speciesY['MoleSecr2'] == 1 & speciesX['MoleUsa2'] == 1) {
    score = score + 0.12
  }
  if (speciesY['MoleSecr3'] == 1 & speciesX['MoleUsa3'] == 1) {
    score = score + 0.12
  }
  #vice-versa
  if (speciesX['MoleSecr1'] == 1 & speciesY['MoleUsa1'] == 1) {
    score = score + 0.12
  }
  if (speciesX['MoleSecr2'] == 1 & speciesY['MoleUsa2'] == 1) {
    score = score + 0.12
  }
  if (speciesX['MoleSecr3'] == 1 & speciesY['MoleUsa3'] == 1) {
    score = score + 0.12
  }
  #if any of antibody secreted by Y can protect X, score += 0.15
  if (speciesY['AntiSecr1'] == 1 & speciesX['AntiProt1'] == 1) {
    score = score + 0.15
  }
  if (speciesY['AntiSecr2'] == 1 & speciesX['AntiProt2'] == 1) {
    score = score + 0.15
  }
  if (speciesY['AntiSecr3'] == 1 & speciesX['AntiProt3'] == 1) {
    score = score + 0.15
  }
  #vice-versa
  if (speciesX['AntiSecr1'] == 1 & speciesY['AntiProt1'] == 1) {
    score = score + 0.15
  }
  if (speciesX['AntiSecr2'] == 1 & speciesY['AntiProt2'] == 1) {
    score = score + 0.15
  }
  if (speciesX['AntiSecr3'] == 1 & speciesY['AntiProt3'] == 1) {
    score = score + 0.15
  }
  return(score)
}

