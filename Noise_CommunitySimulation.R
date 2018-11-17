library(deSolve)

#Define GLV with varying coefficient
glv <- function(t, x, params){
  with(as.list(params, c(x)),{
    dx <- alpha*x + x*as.vector(c0%*%x)+x*(t(do.call(cbind, lapply(l, FUN=function(ma) ma%*%x)))%*%x)
    list(dx)
  })
}

#Define integration method
n.integrate <- function(time, init, model, params){
  as.data.frame(lsoda(init, time, model, params))
}

growthFunction <- function(N, alpha, c0, l, init){
  
  #Solve the ode
  dat <- n.integrate(0:300, init, glv, list(alpha=alpha, c0=c0, l=l))
  
  #Plot
  matplot(x=dat$time, y=dat[,-1], typ='b', xlab='time', ylab='Absolute abundance', main=paste('Modified GLV-density', N,'species'))
  
  return(dat)
}