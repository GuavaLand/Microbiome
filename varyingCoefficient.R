library(deSolve)

#Define GLV with varying coefficient
glv <- function(t, x, params){
  with(as.list(params, c(x)),{
    dx <- alpha*x + x*(c0%*%x)+x*((ck1%*%x)%*%(ck2%*%x))
    list(dx)
  })
}

#Define integration method
n.integrate <- function(time, init.x, model, params){
  as.data.frame(ode(init.x, time, model, params))
}


#Define community size
N <- 10

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
ck1 <- sample(c(1),N,replace = TRUE)
ck2 <- runif(N, min = -1, max = 0.5)

#Define initial abundance between 0.1 and 1, to 1 decimal place
init.x <- floor(runif(N, min = 1, max = 10))/10

#Solve the ode
dat <- n.integrate(0:10, init.x, glv, list(alpha=alpha, c0=c0, ck1 = ck1, ck2 = ck2))

#Plot
matplot(x=dat$time, y=dat[,-1], typ='b', xlab='time', ylab='Absolute abundance')