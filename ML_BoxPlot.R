dat500 = read.csv('C:\\Source\\Microbiome\\score10_500.csv',sep = '\t')
dat600 = read.csv('C:\\Source\\Microbiome\\score10_600.csv',sep = '\t')
dat700 = read.csv('C:\\Source\\Microbiome\\score10_700.csv',sep = '\t')
dat800 = read.csv('C:\\Source\\Microbiome\\score10_800.csv',sep = '\t')
dat900 = read.csv('C:\\Source\\Microbiome\\score10_900.csv',sep = '\t')
dat1000 = read.csv('C:\\Source\\Microbiome\\score10_1000.csv',sep = '\t')

dat = rbind(dat500,dat600,dat700,dat800,dat900,dat1000)

rownames(dat) = 1:60