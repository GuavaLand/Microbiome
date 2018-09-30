#Read HOI_SSN and GLV_SSN matrices: 1000 simulation (row) of time when species reached steady state (column). N refers to community size of simulation
HOI_SS10 = read.csv('HOI_SS10.csv')
HOI_SS20 = read.csv('HOI_SS20.csv')
HOI_SS30 = read.csv('HOI_SS30.csv')
HOI_SS40 = read.csv('HOI_SS40.csv')
HOI_SS50 = read.csv('HOI_SS50.csv')
GLV_SS10 = read.csv('GLV_SS10.csv')
GLV_SS20 = read.csv('GLV_SS20.csv')
GLV_SS30 = read.csv('GLV_SS30.csv')
GLV_SS40 = read.csv('GLV_SS40.csv')
GLV_SS50 = read.csv('GLV_SS50.csv')

HOI_list = list(HOI_SS10,HOI_SS20,HOI_SS30,HOI_SS40,HOI_SS50)
GLV_list = list(GLV_SS10,GLV_SS20,GLV_SS30,GLV_SS40,GLV_SS50)

CommunitySteadyStateDF <- data.frame()

#Loop through HOI_SS10, HOI_SS20...GLV_SS10, GLV_SS20...
for (item in 1:length(HOI_list)) {
  
  HOI = HOI_list[[item]]
  GLV = GLV_list[[item]]
  
  numOfEntries = nrow(HOI)
  
  #CommunitySteadyStateTime to store the time when the whole community has reached steady state
  CommunitySteadyStateTime <- matrix(nrow = numOfEntries, ncol = 3)
  colnames(CommunitySteadyStateTime) <- c('CommunitySize', 'SteadyStateTimeHOI', 'SteadyStateTimeGLV')
  
  for (i in 1:numOfEntries) {
    HOINA <- any(is.na(HOI[i,]))
    GLVNA <- any(is.na(GLV[i,]))
    #If there is NA in a particular row, for wither HOI or GLV, this community combination did not reach SS within integration time. Move to next row
    if (!HOINA & !GLVNA) {
      #If no NA in both HOI and GLV (for that row, both model reached steady state)
      CommunitySteadyStateTime[i,1] <- item * 10
      CommunitySteadyStateTime[i,2] <- max(HOI[i,2:ncol(HOI)])
      CommunitySteadyStateTime[i,3] <- max(GLV[i,2:ncol(GLV)])
    }
    else{
      CommunitySteadyStateTime[i,1] <- item * 10
      CommunitySteadyStateTime[i,2] <- NA
      CommunitySteadyStateTime[i,3] <- NA
    }
  }
  #Finished looping through an HOI_SSN and GLV_SSN, store the CommunitySteadyStateTime vector to a list
  CommunitySteadyStateDF <-  rbind(CommunitySteadyStateDF,as.data.frame(CommunitySteadyStateTime))
}


##plot
#par(mfrow = c(3,2))
##plot histogram for Community Steady State time for each of the community size 
#hist(CommunitySteadyStateDF[CommunitySteadyStateDF$CommunitySize == 10,2], main = 'Community Steady State Time (HOI Model) for 10 Species', xlab = 'Community SS Time', ylim = c(0,500))
#hist(CommunitySteadyStateDF[CommunitySteadyStateDF$CommunitySize == 20,2], main = 'Community Steady State Time (HOI Model) for 20 Species', xlab = 'Community SS Time', ylim = c(0,500))
#hist(CommunitySteadyStateDF[CommunitySteadyStateDF$CommunitySize == 30,2], main = 'Community Steady State Time (HOI Model) for 30 Species', xlab = 'Community SS Time', ylim = c(0,500))
#hist(CommunitySteadyStateDF[CommunitySteadyStateDF$CommunitySize == 40,2], main = 'Community Steady State Time (HOI Model) for 40 Species', xlab = 'Community SS Time', ylim = c(0,500))
#hist(CommunitySteadyStateDF[CommunitySteadyStateDF$CommunitySize == 50,2], main = 'Community Steady State Time (HOI Model) for 50 Species', xlab = 'Community SS Time', ylim = c(0,500))
##plot boxplotboxplot(SteadyStateTime~CommunitySize ,data = CommunitySteadyStateDF)
#boxplot(SteadyStateTime~CommunitySize ,data = CommunitySteadyStateDF)
