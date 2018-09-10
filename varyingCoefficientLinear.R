library(deSolve)

#Define GLV with varying coefficient
glv <- function(t, x, params){
  with(as.list(params, c(x)),{
    dx <- alpha*x + x*(c0%*%x)+x*(x%*%(cbind(l[[1]]%*%x,l[[2]]%*%x,l[[3]]%*%x,
                                             l[[4]]%*%x,l[[5]]%*%x,l[[6]]%*%x,
                                             l[[7]]%*%x,l[[8]]%*%x,l[[9]]%*%x,l[[10]]%*%x)))
    list(dx)
  })
}

#Define integration method
n.integrate <- function(time, init.x, model, params){
  as.data.frame(lsoda(init.x, time, model, params))
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
ck <- runif(N*N, min = -1, max = 0.2)
ck <- matrix(ck, nrow = N)
#Modify ck where k == j
for (i in 1:N) {
  for (j in 1:N) {
    if (i==j) {
      ck[i,j] <- 0
    }
  }
}
#Modify ck where k == i
temp <- ck
temp[,1] <- 0
l <- list()
l[[1]] <- temp
for (i in 2:N) {
  temp1 <- ck
  temp1[,i] <- 0
  l[[i]] <- temp1
}

#Define initial abundance between 0.1 and 1, to 1 decimal place
init.x <- floor(runif(N, min = 1, max = 10))/10

#Solve the ode
dat <- n.integrate(0:10, init.x, glv, list(alpha=alpha, c0=c0, l=l))

#Plot
matplot(x=dat$time, y=dat[,-1], typ='b', xlab='time', ylab='Absolute abundance')



###############################################################################
#Generate 2^N communities where species in each community (max possible N) are present or absent
###############################################################################

#Generate boolean mask of 2^N x N matrix
repeatBinaryNTimes <- rep(list(c(0,1)),N)
mask <- expand.grid(repeatBinaryNTimes)

#Loop through 2^N to apply each row in mask to initial abundance, solve ode, and retrieve final abundance
SSMatrix <- as.data.frame(matrix(nrow = 2^N, ncol = N))
colnm <- c(1:10)
colnames(SSMatrix) <- colnm
for (i in 1:2^N) {
  init <- init.x*mask[i,]
  init <- as.numeric(init)
  dat <- n.integrate(0:10, init, glv, list(alpha=alpha, c0=c0, ck1 = ck1, ck2 = ck2))
  SSMatrix[i,] <- dat[nrow(dat),2:(N+1)]
}


###############################################################################
#Machine learning
###############################################################################
# library(randomForest)
# for(i in 1:10){
#   dead <- which(mask[,i] == 0)
#   final <- SSMatrix[-dead,i]
#   mod <- randomForest(final ~., mask[-dead,])
#   plot(mod$y,mod$predicted,main=i)
# }





