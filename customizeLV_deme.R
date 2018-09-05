library(deSolve)

#############################Functions#######################################

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
  #if any of molecule secreted by Y is usable by X, score += 0.09
  #if the molecule secreted by both Y and X, score -= 0.12
  if (speciesY['MoleSecr1'] == 1) {
    if (speciesX['MoleUsa1'] == 1) {
      score = score + 0.09
    }
    if (speciesX['MoleSecr1'] == 1) {
      score = score -0.12
    }
  }
  if (speciesY['MoleSecr2'] == 1) {
    if (speciesX['MoleUsa2'] == 1) {
      score = score + 0.09
    }
    if (speciesX['MoleSecr2'] == 1) {
      score = score -0.12
    }
  }
  if (speciesY['MoleSecr3'] == 1) {
    if (speciesX['MoleUsa3'] == 1) {
      score = score + 0.09
    }
    if (speciesX['MoleSecr3'] == 1) {
      score = score -0.12
    }
  }
  #if any of molecule usable by Y is secreted by X, score += 0.09
  #if the molecule usable by both Y and X, score -= 0.12
  if (speciesY['MoleUsa1'] == 1) {
    if (speciesX['MoleSecr1'] == 1) {
      score = score + 0.09
    }
    if (speciesX['MoleUsa1'] == 1) {
      score = score -0.12
    }
  }
  if (speciesY['MoleUsa2'] == 1) {
    if (speciesX['MoleSecr2'] == 1) {
      score = score + 0.09
    }
    if (speciesX['MoleUsa2'] == 1) {
      score = score -0.12
    }
  }
  if (speciesY['MoleUsa3'] == 1) {
    if (speciesX['MoleSecr3'] == 1) {
      score = score + 0.09
    }
    if (speciesX['MoleUsa3'] == 1) {
      score = score -0.12
    }
  }
  #if any of antibody secreted by Y can protect X, score -= 0.15
  if (speciesY['AntiSecr1'] == 1 & speciesX['AntiProt1'] == 1) {
    score = score - 0.15
  }
  if (speciesY['AntiSecr2'] == 1 & speciesX['AntiProt2'] == 1) {
    score = score - 0.15
  }
  if (speciesY['AntiSecr3'] == 1 & speciesX['AntiProt3'] == 1) {
    score = score - 0.15
  }
  #vice-versa
  if (speciesX['AntiSecr1'] == 1 & speciesY['AntiProt1'] == 1) {
    score = score - 0.15
  }
  if (speciesX['AntiSecr2'] == 1 & speciesY['AntiProt2'] == 1) {
    score = score - 0.15
  }
  if (speciesX['AntiSecr3'] == 1 & speciesY['AntiProt3'] == 1) {
    score = score - 0.15
  }
  return(score)
}

############################################################################

#number of species
N = 10

#Creating random matrix of 12 traits of N species
traitMatrix <-  sample(x = c(0,1), 12, replace = TRUE)

for (i in c(2:N)) {
  vec <-  sample(x = c(0,1), 12, replace = TRUE)
  traitMatrix <-  rbind(traitMatrix, vec)
}
rownames(traitMatrix) <- paste("species",c(1:N),sep = "")
colnames(traitMatrix) <- c('MoleSecr1','MoleSecr2','MoleSecr3','MoleUsa1','MoleUsa2','MoleUsa3',
                           'AntiSecr1','AntiSecr2','AntiSecr3','AntiProt1','AntiProt2','AntiProt3')

#Generate relationship matrix of N species
a = matrix(rep(0,N*N),nrow=N,ncol=N)
for (row in 1:N) {
  for (col in row:N) {
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

#Creating random growth rate of N species
b <- runif(N)

#Parameters of Lotka-Volterra equations
parms <- cbind(b,a)
parms=cbind(rep(N,N),parms)

tstart=0                   # time (start)
tend=30                    # time (end) of integration
tstep=0.1                  # time step (resolution)

init<-runif(N)                # initial species abundances 

#Simulation for coupled system
times<-seq(tstart, tend, by=tstep)

#Generate mask of all boolean combination of N species
repeatBinaryNTimes <- rep(list(c(0,1)),N)
mask <- expand.grid(repeatBinaryNTimes)

#Setting presence/absence of each species based on maskand generate steady state matrix
SSMatrix <- matrix(nrow = 2^N, ncol = N)
for (i in 1:2^N) {
  y <- init*mask[i,]
  y <- as.numeric(y)
  out<-lsoda(y, times, glvmat, parms)
  SSMatrix[i,] <- out[nrow(out),2:(N+1)]
}


####### Figure

col.vec = seq(0,1,1/N)
my.colors = hsv(col.vec)

# plot time series
plot(out[,1],out[,2],col=my.colors[1],lty=1, type="l",ylim=range(out[,-1]))
for(i in 2:N){
  lines(out[,1],out[,(i+1)],col=my.colors[i])
}
legend("topleft", as.character(c(1:N)), lty = rep(1,N), col = my.colors, merge = TRUE, bg = "white", text.col="black")


## Simulation results
library(hier.part)
SSMatrix = combos(10)$binary
Smat <- as.data.frame(SSMatrix)
Outcomes <- SSMatrix
for(i in 1:nrow(Outcomes)){
  y <- SSMatrix[i,]
  out <- lsoda(y, times, glvmat, parms)
  Outcomes[i,] <- out[nrow(out),-1]
}

## Machine Learning Stuff (sloppy)
plot(sort(Outcomes[,1]),ylim=c(0,2.5),cex=1)
for(i in 1:10){
  points(sort(Outcomes[,i]),col=i)
}
final <- Outcomes[,i]

cors <- vector()
par(mfrow=c(3,3))
for(i in 1:10){
  dead <- which(mask[,i] == 0)
  final <- Outcomes[-dead,i]
  mod <- randomForest(final ~., mask[-dead,])
  plot(mod$y,mod$predicted,main=i)
  #lmMod <- lm(final ~., Smat[-dead,])
  
  # tempSmat <- Smat[-dead,-i]
  # svmMod <- svm(final ~., tempSmat,cross = 5)
  # train <- sample(length(final),400)
  # sTrain <- svm(final[train] ~., tempSmat[train,])
  # sTest <- predict(sTrain,tempSmat[-train,])
}


tempSmat <- mask[-dead,]

train <- sample(512,50)

mod <- randomForest(final[train] ~., mask2[train,])

future <- predict(mod,mask2[-train,])

plot(final[-train],future,xlab='real',ylab='predicted')
mtext(cor(final[-train],future))

finals <- rep(final,5)
Z <- rbind(tempSmat,tempSmat,tempSmat,tempSmat,tempSmat)

noise <- rnorm(mean=1,sd=.1,length(finals))
finals <- finals*noise
mood <- randomForest(finals ~., Z)

noiseTrain <- sample(length(final),200)
train <- noiseTrain
for(i in 2:5){
  train <- c(train,train+length(final))
}

mood2 <- randomForest(finals[train] ~., Z[train,])
pred2 <- predict(mood2,Z[-train,])
cor(pred2,finals[-train])
plot(pred2,finals[-train])
abline(0,1)
summary(abs(finals[-train] - pred2)^2)

mood2 <- svm(finals[train] ~., Z[train,])
pred2 <- predict(mood2,Z[-train,])
cor(pred2,finals[-train])
summary(abs(finals[-train] - pred2)^2)
