library(deSolve)

#Define community size
#N = 10

#Define Naive model
nm <- function(t, x, params){
  with(as.list(params, c(x)),{
    dx <- x + x*(c0%*%x)
    list(dx)
  })
}

#Define integration method
nm.integrate <- function(time, init.x, model, params){
  as.data.frame(lsoda(init.x, time, model, params))
}

naiveModel <- function(N){
  #Define intial density
  init.x <- floor(runif(N)*10)/10
  
  #Define c0 as 1/N
  c0 <- matrix(rep(-1/N, N^2), nrow = N)
  #Define intra-species effect
  diag(c0) <- -0.5
  
  #Solve ode
  dat <- nm.integrate(1:100, init.x, nm, list(c0 = c0))
  
  #plot
  matplot(x = dat$time, y = dat[,-1], type = 'b', xlab = 'time', ylab = 'Absolute abundance', main = paste('Naive Model',N,'Species'))
  
  return(dat)
}

