library(ggplot2)

###############################################################
#This function receives three rate components (separate matrices) of all species (column) over time (row)
#and plot rate component over time for each species (N plots)
###############################################################

plotRateComponent <- function(term1, term2, term3){
  
  ########################################
  #cus term3 dies out quickly, if want to look at only the beginning, run the following section
  term1 <- term1[1:60,] #look at rate from begining to 60 second
  term2 <- term2[1:60,]
  term3 <- term3[1:60,]
  ##########################################
  
  
  #number of species
  N = ncol(term1) - 2
  
  #total rate
  rate = term1[,1:N] + term2[,1:N] + term3[,1:N]
  rate$time = term1$time
  
  #Iterate through all species
  for (n in 1:N) {
    print(ggplot(term1,aes_string(x = 'time',y = paste('species',n,sep='')))+
      geom_area(aes(fill = term),alpha=0.5) + geom_area(data=term2,aes(fill = term),alpha=0.5)+
      geom_area(data = term3,aes(fill = term),alpha = 0.5)+
      geom_line(data = rate)) #superimpose with overall growth rate 
  }
  
  
}