lambda_matrix = NULL
for(i in c(c(1:30),c(32:48),  c(50:100))){
file1 = paste("/Volumes/Liang/Research/data/DDItest_lambda/DDI",i,".Rdata",sep = "")
load(file = file1)
lambda_matrix = rbind(lambda_matrix, fit$par)
}
df_lambda = data.frame(value = lambda_matrix, id = "lambda")
plot_p <- ggplot(df_lambda, aes(x=id, y=value)) +geom_boxplot()+
  geom_hline(aes(yintercept=0.8))

plot_p