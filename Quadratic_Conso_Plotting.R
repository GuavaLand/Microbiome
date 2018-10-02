dat = as.data.frame(read.csv('Quadratic_Conso_Difference_matrix.csv'))
dat = dat[,2:15]

par(mfrow = c(4,2))
for (i in 1:ncol(dat)) {
  hist(dat[,i],main = paste('Time Difference for', i*5, 'Species'), xlab = 'Time Difference')
}
boxplot(dat, use.cols = TRUE, main = 'Time Different for Community of Different Size', xlab = 'Community Size', ylab = 'Frequency')
