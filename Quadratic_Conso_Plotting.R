library(ggplot2)

dat = as.data.frame(read.csv('Quadratic_Conso_Difference_matrix.csv'))
dat = dat[,2:15]

par(mfrow = c(4,2))
for (i in 1:ncol(dat)) {
  hist(dat[,i],main = paste('Time Difference for', i*5, 'Species'), xlab = 'Time Difference')
}
boxplot(dat, use.cols = TRUE, main = 'Time Different for Community of Different Size', xlab = 'Community Size', ylab = 'Frequency')

#plot scatter plot and error bar
colmean = colMeans(dat,na.rm = TRUE)
colsd = apply(dat, MARGIN=2, FUN=sd, na.rm=TRUE)
dat2 = t(rbind(colmean,colsd))
dat2 = as.data.frame(dat2)
#rename X5 to X05 
row.names(dat2)[1] = "X05"
#plotting
p<- ggplot(dat2, aes(x=row.names(dat2), y=colmean)) + 
  geom_line() +
  geom_point()+
  geom_errorbar(aes(ymin=colmean-colsd, ymax=colmean+colsd), width=.2,
                position=position_dodge(0.05))
print(p)
#title, labels etc.
p+labs(title="Time Difference for Community of Different Size", x="Community Size", y = "Time Difference")
  
