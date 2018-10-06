library(deSolve)
library(ggplot2)

#Define GLV with varying coefficient
glv <- function(t, x, params){
  with(as.list(params, c(x)),{
    dx <- alpha*x + x*as.vector(c0%*%x)+x*t(t(x)%*%(do.call(cbind, lapply(l, FUN=function(ma) ma%*%x))))
    list(dx)
  })
}

#Define integration method
n.integrate <- function(time, init.x, model, params){
  as.data.frame(lsoda(init.x, time, model, params))
}

#Define community size
N <- 3

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
l <- list()
for (i in 1:N) {
 #For i-th matrix, i-th row is 0
 temp <- ck
 temp[i,] <- 0
 
 #In i-th matrix, elements are 0 if k (column) == either j (row) or i (matrix order)
 for (j in 1:N) {
   for (k in 1:N) {
     if (k==j | k == i) {
       temp[j,k] <- 0
     }
   }
 }
 
 #Control the prevalence of thrid party effects
 #for (element in 1:length(temp)) {
 #  dice <- runif(1)
 #  if (dice > -1) { #what percent of to assign 0
 #    temp[element] <- 0
 #  }
 #}
 l[[i]] <- temp
}


#Define initial abundance between 0.1 and 1, to 1 decimal place
init.x <- floor(runif(N, min = 1, max = 10))/10

#Solve the ode
dat <- n.integrate(0:20, init.x, glv, list(alpha=alpha, c0=c0, l=l))

#Plot
matplot(x=dat$time, y=dat[,-1], typ='b', xlab='time', ylab='Absolute abundance', main=paste('Modified GLV-density', N,'species'))

#################################
#Find the components of the rate
#################################
dat_density <- dat[,2:ncol(dat)]
#First term of growth rate: multiply alpha to x
term1 <- apply(dat_density,1,function(x) alpha*x)
term1 <- as.data.frame(t(term1))
colnames(term1) <- paste('term1_species', 1:N)
term1$time <- dat$time


#Second term of growth rate
term2 <- apply(dat_density,1,function(x) x*as.vector(c0%*%x))
term2 <- as.data.frame(t(term2))
colnames(term2) <- paste('term2_species', 1:N)
term2$time <- dat$time

#Third term of growth rate
term3 <- apply(dat_density,1,function(x) x*t(t(x)%*%(do.call(cbind, lapply(l, FUN=function(ma) ma%*%x)))))
term3 <- as.data.frame(t(term3))
colnames(term3) <- paste('term3_species', 1:N)
term3$time <- dat$time

#Total growth rate
rate <- term1 + term2 + term3
colnames(rate) <- paste('rate_species', 1:N)
rate$time <- dat$time

##################################
#plot
##################################
newdf <- data.frame(term1[,1], term2[,1], term3[,1], dat$time)
ggplot(newdf,aes(x = dat.time)) + ylim(0,2)+ geom_area(stat = 'bin',binwidth = 50)
