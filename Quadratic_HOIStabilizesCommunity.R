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

#Calculate the difference of community SS time between GLV model and and HOI model 
CommunitySteadyStateDF$Difference <- CommunitySteadyStateDF$SteadyStateTimeGLV - CommunitySteadyStateDF$SteadyStateTimeHOI

#plot
par(mfrow = c(3,2))
#plot histogram for Community Steady State time difference for each of the community size 
hist(CommunitySteadyStateDF[CommunitySteadyStateDF$CommunitySize == 10,4], main = 'How Much Longer Does GLV Reaches SS Than HOI 10 Species', xlab = 'Time Difference', ylim = c(0,500), breaks = 18)
hist(CommunitySteadyStateDF[CommunitySteadyStateDF$CommunitySize == 20,4], main = 'How Much Longer Does GLV Reaches SS Than HOI 20 Species', xlab = 'Time Difference', ylim = c(0,500), breaks = 18)
hist(CommunitySteadyStateDF[CommunitySteadyStateDF$CommunitySize == 30,4], main = 'How Much Longer Does GLV Reaches SS Than HOI 30 Species', xlab = 'Time Difference', ylim = c(0,500), breaks = 18)
hist(CommunitySteadyStateDF[CommunitySteadyStateDF$CommunitySize == 40,4], main = 'How Much Longer Does GLV Reaches SS Than HOI 40 Species', xlab = 'Time Difference', ylim = c(0,500), breaks = 18)
hist(CommunitySteadyStateDF[CommunitySteadyStateDF$CommunitySize == 50,4], main = 'How Much Longer Does GLV Reaches SS Than HOI 50 Species', xlab = 'Time Difference', ylim = c(0,500), breaks = 18)
#plot boxplotboxplot(SteadyStateTime~CommunitySize ,data = CommunitySteadyStateDF)Differences
boxplot(Difference~CommunitySize ,data = CommunitySteadyStateDF, main = 'How Much Longer Does GLV Reaches SS Than HOI', xlab = 'Community Size', ylab = 'Time Difference')
