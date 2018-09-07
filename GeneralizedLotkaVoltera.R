library(deSolve)
#Generalized Lotka-Voltera
#Assuming 3 species, X1, X2 and X3
#Growth rate of each species proportional to species number
#Coefficient a function of all species number
#dx1/dt = (a1 +b1X1 + c1x2 + d1X3)X1
#dx2/dt = (a2 +b2X1 + c2x2 + d2X3)X2
#dx3/dt = (a3 +b3X1 + c3x2 + d3X3)X3

glv <- function()