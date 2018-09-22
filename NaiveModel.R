library(deSolve)

#Define community size
N = 50

#Define Naive model
nm <- function(t, x, params){
  with(as.list(params, c(x)),{
    dx <- x - x*as.vector(c0%*%x)
    list(dx)
  })
}

#Define integration method
n.integrate <- function(time, init.x, model, params){
  as.data.frame(lsoda(init.x, time, model, params))
}

#Define intial density
init.x <- floor(runif(N)*10)/10

#Define c0 as 1/N
c0 <- matrix(rep(1/N, N^2), nrow = N)
#Excluding the species itself
diag(c0) <- 0


