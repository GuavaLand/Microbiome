library(deSolve)

#Modified Lotka-Voltera (trait-mediated)
#In the generalized LV, c0 describes species j's influence on i's growth rate.
#We incorporate a new term into c0 to describe HOI of third species.
#HOI is defined as a list of length N, with each member an N*N matrix. i-th matrix describes species i's effect on row-column pair.
#Therefore to get the whole community's HOI effect on a particular pair, 
#we just need to element-sum the list and find the value at corresponding row and column.
#Before doing element-sum of HOI, we need to consider when a i-th species may be absent.
#When this happens, HOI i-th matrix should mutiply by 0 before doing element-sum of HOI.
#So when defining the function, we mutiply i-th matrix with (-exp(-5*xi)+1) (mapply function).
#(-exp(-5*x)+1) is a special function that = 0 when x = 0 and approaches 1 when x increases (90% of 1 when x reaches 0.46).
#This ensures that HOI effect of i-th species only exists when i exists, and the effect remains unchanged as i's density increases.
glv <- function(t, x, params){
  with(as.list(params, c(x)),{
    dx <- alpha*x + x*((c0+Reduce('+', mapply('*',HOI,(-exp(-5*x)+1),SIMPLIFY = FALSE)))%*%x)
    list(dx)
  })
}

#Original Generalized Lotka-Voltera
glv1 <- function(t, x, params){
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
N <- 50

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

#Define initial density
init.x <- floor(runif(N, min=0.1, max=1)*10)/10

#Solve the ode
dat <- n.integrate(1:100, init.x, glv, list(alpha=alpha, c0=c0, HOI=HOI))
dat1 <- n.integrate(1:100, init.x, glv1, list(alpha=alpha, c0=c0))

#Plot
matplot(x=dat$time, y=dat[,-1], typ='b', xlab='time', ylab='Absolute abundance', main='Modified GLV')
matplot(x=dat1$time, y=dat1[,-1], typ='b', xlab='time', ylab='Absolute abundance', main='Original GLV')
