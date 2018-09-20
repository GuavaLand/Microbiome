library(randomForest)

for(i in 1:10){
 dead <- which(mask[,i] == 0)
 final <- SSMatrix[-dead,i]
 mod <- randomForest(final ~., mask[-dead,])
 plot(mod$y,mod$predicted,main=i)
}