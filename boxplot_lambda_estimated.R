lambda_matrix = NULL
p = list()
for(i in c(c(1:32),c(34:94),c(96:247),c(249:370),c(372:500) )){
file1 = paste("/Volumes/Liang/Research/data/DDItest_lambda500/DDI",i,".Rdata",sep = "")
load(file = file1)
lambda_matrix = rbind(lambda_matrix, fit$par)
}
df_lambda = data.frame(value = lambda_matrix, id = "lambda")
plot_p <- ggplot(df_lambda, aes(x=id, y=value)) +geom_boxplot()+
  geom_hline(aes(yintercept=0.8))

p[[1]] = plot_p
p[[2]] = ggplot(df_lambda, aes(x=value)) +geom_histogram()
multiplot(plotlist=p,cols=1)
