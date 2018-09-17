library(deSolve)

#Generalize Lotka-Voltera
glv <- function(t, x, params){
  with(as.list(params, c(x)),{
    dx <- alpha*x + x*(c0%*%x)
    list(dx)
  })
}

#Define integration method
n.integrate <- function(time, init.x, model, params){
  as.data.frame(lsoda(init.x, time, model, params))
}

#Number of species
N <- 10

#Create matrix of HOI (effect on growth rate)
#N matrices in total, each specifying on species effect on other pairs
HOI <- list()

for (i in 1:N) {
  #Species i will have HOI with N-1 species pairwise combination
  temp <- matrix(runif(N*N, min = -0.9, max = 0.1),nrow = N)
  #HOI effect of i in j-k pair and k-j pair should be the same
  for (ro in 1:N) {
    for (co in ro:N) {
      temp[ro,co] <- temp[co,ro]
    }
  }
  #Name row and col based on species ID
  rownames(temp) <- 1:N
  colnames(temp) <- 1:N
  #i's HOI with i-otherSpecies pair is 0
  temp[i,] <- 0
  temp[,i] <- 0
  #There is no HOI on intra-species interation
  diag(temp) <- 0
  
  HOI[[i]] <- temp
}

#Define species intrinsic growth rate
alpha <- runif(N)

#Define pairwise interaction
c0 <- matrix(runif(N*N, min = -1, max = 0),nrow = N)
#Set intra-species interaction to -0.5
diag(c0) <- -0.5
#Equate i-j interaction to j-i interaction
for (ro in 1:N) {
  for (co in ro:N) {
    c0[ro,co] <- c0[co,ro]
  }
}

#Adjust pairwise interaction based on HOI matrix


#Define initial density
init.x <- floor(runif(N, min=0.1, max=1)*10)/10
