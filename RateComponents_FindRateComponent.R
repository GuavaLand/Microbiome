#################################
#Find the components of the rate
#Input simulated density over time, species ID as column; alpha, c0, l
#Outpu matrices of broken down growth rate
#################################

findRateComponent <- function(dat, alpha, c0, l){
  #Create local variable to store input
  dat = dat
  alpha = alpha
  c0 = c0
  l = l
  
  N = ncol(dat) - 1
  
  #Plot dat
  matplot(x=dat$time, y=dat[,-1], typ='b', xlab='time', ylab='Absolute abundance', main=paste('Modified GLV-density', N,'species'))
  
  dat_density <- dat[,2:ncol(dat)]
  #First term of growth rate: multiply alpha to x
  term1 <- apply(dat_density,1,function(x) alpha*x)
  term1 <- as.data.frame(t(term1))
  colnames(term1) <- paste('species', 1:N, sep='')
  term1$time <- dat$time
  term1$term <- factor(rep(1,nrow(term1)))
  
  
  #Second term of growth rate
  term2 <- apply(dat_density,1,function(x) x*as.vector(c0%*%x))
  term2 <- as.data.frame(t(term2))
  colnames(term2) <- paste('species', 1:N, sep='')
  term2$time <- dat$time
  term2$term <- factor(rep(2,nrow(term2)))
  
  #Third term of growth rate
  term3 <- apply(dat_density,1,function(x) x*t(t(x)%*%(do.call(cbind, lapply(l, FUN=function(ma) ma%*%x)))))
  term3 <- as.data.frame(t(term3))
  colnames(term3) <- paste('species', 1:N, sep='')
  term3$time <- dat$time
  term3$term <- factor(rep(3,nrow(term3)))
  
  returnList <- list(term1 = term1, term2 = term2, term3 = term3)
  
  return(returnList)
  
}



