HOI_SS10 = read.csv('HOI_SS10.csv')
HOI_SS20 = read.csv('HOI_SS20.csv')
HOI_SS30 = read.csv('HOI_SS30.csv')
HOI_SS40 = read.csv('HOI_SS40.csv')
HOI_SS50 = read.csv('HOI_SS50.csv')

HOI_list = list(HOI_SS10,HOI_SS20,HOI_SS30,HOI_SS40,HOI_SS50)

CommunitySteadyStateDF <- data.frame()
#The order of HOI matrices stored in the list
HOI_i <- 10

#Loop through HOI_SS10, HOI_SS20...
for (HOI in HOI_list) {
  numOfEntries = nrow(HOI)
  #CommunitySteadyStateTime to store the time when the whole community has reached steady state
  CommunitySteadyStateTime <- matrix(nrow = numOfEntries, ncol = 2)
  colnames(CommunitySteadyStateTime) <- c('CommunitySize', 'SteadyStateTime')
  for (i in 1:numOfEntries) {
    #If there is NA in a particular row, this community combination did not reach SS within integration time. Move to next row
    if (!any(is.na(HOI[i,]))) {
      CommunitySteadyStateTime[i,1] <- HOI_i
      CommunitySteadyStateTime[i,2] <- max(HOI[i,2:ncol(HOI)])
    }
    else{
      CommunitySteadyStateTime[i,1] <- HOI_i
      CommunitySteadyStateTime[i,2] <- NA
    }
  }
  #Finished looping through an HOI matrix, store the CommunitySteadyStateTime vector to a list
  CommunitySteadyStateDF <-  rbind(CommunitySteadyStateDF,as.data.frame(CommunitySteadyStateTime))
  HOI_i <- HOI_i+10
}






