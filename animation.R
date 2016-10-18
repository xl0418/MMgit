library(ggplot2)
set.seed(2016)
index_Score <- function(){
  # Picking 2 points randomly on the stick at the same time
  x <- runif(n = 2, min = 0, max = 1) 
  a <- min(x) # first point
  b <- max(x) # second point
  # pieces of the stick with their respective length
  pieces <- c(a, b-a, 1-b)
  cond1 <- sum(pieces[c(1,2)]) > pieces[3] # condition # 1
  cond2 <- sum(pieces[c(1,3)]) > pieces[2] # condition # 2
  cond3 <- sum(pieces[c(3,2)]) > pieces[1] # condition # 3
  combine_conds <- ifelse(cond1 & cond2 & cond3, 1, 0) # if all 3 conditions are satisfied
  return(combine_conds)
}

cnt <- c()
total <- 1000
for(k in 1:total) cnt = c(cnt, index_Score())
df <- setNames(data.frame(1:total, rep(0, total)), c("Incrmt","Probs"))
for (i in 1:total)  df$Probs[i] <- sum(cnt[1:i])/i
for(i in 1:total) {
  sub_df <- subset(df, df$Incrmt <= i)
  simul_plot <- qplot(Incrmt, Probs, data = sub_df, geom = "path") + 
    labs(x = "iterations", y = "Probabilities", title = "Monte Carlo Simulation") + ylim(c(0,0.4)) + 
    geom_hline(yintercept = 0.25, colour = "red", linetype = "longdash")
  ggsave(plot = simul_plot, filename = paste(sprintf("/Users/mac/Dropbox/R/animation/images/brokenstick_%02d",i),".png", sep = ""), limitsize = FALSE)
  rm(sub_df)
  # dev.off()
}
dev.off()
# ffmpeg -r 10 -i /Users/mac/Dropbox/R/animation/images/brokenstick_%02d.png -b:v 20M /Users/mac/Dropbox/R/animation/BrokenStick_video.mp4
