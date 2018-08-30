#####################################
# Run a simulation with the generalized 
# Lotka-Volterra model
#
# Authors: Didier Gonze & Karoline Faust 
#
# Date: 24/10/2011

library(deSolve)

# Parameters

N = 4                 	     # number of species

b=runif(N)                   # intrinsic growth rate vector
a=matrix(nrow=N,ncol=N)      # initialize species interaction matrix
for (i in 1:N){              
  for(j in 1:i){
    if(i==j){
      a[i,j]=-0.5;         # auto-interactions are negative (intra-species competition)      
      a[j,i]=-0.5
    }else{
      a[i,j] = (runif(1) - 0.5)    
      a[j,i] = a[i,j]     # enforce symmetry
    }
  }
}

print(a)

tstart=0                   # time (start)
tend=30                    # time (end) of integration
tstep=0.1                  # time step (resolution)

y<-runif(N)                # initial species abundances 

# parms as matrix
parms=cbind(b,a)
parms=cbind(rep(N,N),parms)

#### Simulation for the coupled system

times<-seq(tstart, tend, by=tstep)

out<-lsoda(y, times, glvmat, parms)


####### Figure

col.vec = seq(0,1,1/N)
my.colors = hsv(col.vec)

# plot time series
plot(out[,1],out[,2],col=my.colors[1],lty=1, type="l")
for(i in 2:N){
  lines(out[,1],out[,(i+1)],col=my.colors[i])
}
legend("topright", as.character(c(1:N)), lty = rep(1,N), col = my.colors, merge = TRUE, bg = "white", text.col="black")


# ==============================================================
# Equations (generalized Lotka-Volterra)
# ==============================================================

# matrix formulation of the ODE set
# t: current simulation time
# y: vector with current values of state variables (initial conditions)
# parms: parameter values
# 
glvmat<-function(t, y, parms){
  N=parms[1,1]  # species number
  b=parms[,2]   # vector of growth rates
  a=parms[,3:(N+2)] # interaction matrix
  dydt <- y*(b+a %*% y)
  list(dydt)
}
