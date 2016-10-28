source('~/Dropbox/R/MMgit/event_matrix.R')
# source('~/event_matrix.R')
library(Matrix)
library(ggplot2)
index_loc = event_matrix(4)
file = paste("/Users/mac/Dropbox/R/cluster/4Modeltestgroup1/out10sim.Rdata",sep = "")
# file = paste(getwd(), "/out10sim.Rdata",sep = "")
load(file = file)

Ntable = matrix(0,nrow = nrow(result$Ntable),ncol = 4)
for(i in 1:4){
Ntable[,i] = rowSums(result$Ntable[,which(index_loc[i,] ==1)])
}
total <- length(result$t)-1
df <- setNames(data.frame(1:total, rep(0, total),rep(0, total),rep(0, total),rep(0, total)), c("time","num_a","num_b","num_c","num_d"))
time_sim = result$t[-length(result$t)]
 df$num_a <- Ntable[,1]
 df$num_b <- Ntable[,2]
 df$num_c <- Ntable[,3]
 df$num_d <- Ntable[,4]
 # num = cbind(df$num_a,df$num_b,df$num_c,df$num_d)
 df$time <- time_sim
 df_notime = stack(df,select = -time)
 df_notime$time = rep(df$time,4)
 
 # mydf = setNames(data.frame(t(df[,-1])), t(df[,1]))
 # mydf = t(mydf)
 # total = 100
for(i in 1:total) {
  sub_df <- subset(df_notime, df_notime$time <= df_notime$time[i])
  # num = cbind(sub_df$num_a,sub_df$num_b,sub_df$num_c,sub_df$num_d)
  # simul_plot <- matplot(sub_df$time, num, type = "l",lty = 1)# + 
  
  simul_plot <- qplot(time, values, data = sub_df,group = ind, colour = ind, geom = "path") +
    labs(x = "time", y = "number of lineages", title = "Migration simulation for 4 locations") + # + ylim(c(0,0.4)) +
    geom_hline(yintercept = 25, colour = "red", linetype = "longdash")
  ggsave(plot = simul_plot, filename = paste(sprintf("/Users/mac/Dropbox/R/animation/images/Mig_sim_%02d",i),".png", sep = ""), limitsize = FALSE)
  # ggsave(plot = simul_plot, filename = paste(sprintf(getwd(), "Mig_sim_%02d",i),".png", sep = ""), limitsize = FALSE)
  rm(sub_df)
  # dev.off()
}
dev.off()
# ffmpeg -r 10 -i /Users/mac/Dropbox/R/animation/images/Mig_sim_%02d.png -b:v 50M /Users/mac/Dropbox/R/animation/Mig_4_sim_video_each_loc2.mp4
