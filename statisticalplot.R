library(ggplot2)
library(grid)
fileplot = paste(getwd(),"/Dropbox/R/Migrationmodelsim/ModelTest/multiplot.R",sep = "")
source(file = fileplot)


# boxplot of p-value and power of test for different
# p = rep(0,40)
# power = rep(0,40)



p = NA
power = NA
for(j in 20:24){
  file1 = paste(getwd(),"/Desktop/data/Modeltestgroup",j,sep = "")
for(i in 1:100){
  if(j == 20) file = paste(file1,"/Testanalysis",i,".Rdata",sep = "")
  else file = paste(file1,"/analysis",i,".Rdata",sep = "")
  load(file = file)
  # print(paste("pvalue",i," = ",pvalue,"; poweroftest",i," = ",poweroftest,sep = "")) 
  i <- 100*(j-20)+i
  p[i] <- pvalue
  power[i] <- poweroftest
}
  print(length(p))
  print(length(power))
}
group = rep(c(1,0.8,0.5,5,0.3),each = 100)
# group = rep(c(50,0.1,0.2,0.15,0.3,0.35,0.4,0.4,10,1,0.8,0.9,1),each = 10)
p2 = cbind(group, p , power)
p2 = as.data.frame(p2)
 plot_p <- ggplot(p2, aes(x=factor(group), y=p)) +geom_boxplot()+xlab("Migration rate")+ylab("p-value")+#geom_hline(yintercept = 0.05)
 geom_hline(aes(yintercept=0.05)) + geom_text(aes(0,0.05,label = 0.05, vjust = 0))
   # + geom_dotplot(binaxis='y', stackdir='center', dotsize=1,binwidth = )
# p <- ggplot(p2, aes(x=factor(group), y=power)) +geom_boxplot()
plot_p

plot_power <- ggplot(p2, aes(x=factor(group), y=power)) +geom_boxplot()+xlab("Migration rate")+ylab("p-value")#geom_hline(yintercept = 0.05)
 # geom_hline(aes(yintercept=0.05)) + geom_text(aes(0,0.05,label = 0.05, vjust = 0))
# + geom_dotplot(binaxis='y', stackdir='center', dotsize=1,binwidth = )
# p <- ggplot(p2, aes(x=factor(group), y=power)) +geom_boxplot()
plot_power

multiplot(plot_p, plot_power, cols=2)


# boxplot for K in different groups
K = NULL
for(j in c(8,10,9,11,2,7)){
  file1 = paste(getwd(),"/Desktop/data/Modeltestgroup",j,sep = "")
  for(i in 1:10){
    file = paste(file1,"/analysis",i,".Rdata",sep = "")
    load(file = file)
    # print(paste("pvalue",i," = ",pvalue,"; poweroftest",i," = ",poweroftest,sep = "")) 
    endmc = 1000
    K1 = 1E+120 * (opt == 1) + pmin(1E+120,out$K_DD1[(endmc + 2):(2 * endmc + 1)]) * (opt == 2) + pmin(1E+120,out$K_DD2[(endmc + 2):(2 * endmc + 1)]) * (opt == 3)
    K = c(K,K1)
  }
  # print(p)
}
M0 = rep(c( 0.1, 0.15,0.2, 0.3,10,50),each = 10000)

p2 = cbind(M0,K)
p2 = p2[ p2[,2] < 70]
p3 = matrix(p2, ncol = 2)
p3 = p3[-which(is.na(p3[,2])),]
p4 = as.data.frame(p3)
colnames(p4) = c("M0", "K")
p <- ggplot(p4, aes(x=factor(M0), y=K)) +geom_boxplot()+xlab("M0")+ylab("K")
# + geom_dotplot(binaxis='y', stackdir='center', dotsize=1,binwidth = )
# p <- ggplot(p2, aes(x=factor(group), y=power)) +geom_boxplot()
p

