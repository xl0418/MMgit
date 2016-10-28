library(DDD)
file3 = paste("/Volumes/Liang/Research/data/loc_M0_pics/all_tree.pdf",sep = "")
file4 = paste("/Volumes/Liang/Research/data/loc_M0_pics/all_ltt.pdf",sep = "")

pdf(file = file4, width = 8.3, height = 11.7) 

plotlist_tree = list()
plotlist_ltt = list()
m = matrix(1:15,5,3)
layout(m)
for(i in 1:15){
  if(i<5){
  file = paste("/Volumes/Liang/Research/data/loc_M0_pics/out",i,"sim.txt",sep = "")
  }
  else   file = paste("/Volumes/Liang/Research/data/loc_M0_pics/out",i,"sim.Rdata",sep = "")
  load(file = file)
  L = result$L
  phy = L2phylo(L)
  phy$tip.label = L[which(L[,4] == -1),5]
  
  # loc = result$loctable
  # loc_lable = paste(loc[,1],loc[,2],loc[,3])
  # phy$tip.label = loc_lable
  plot(phy)
  # ltt.plot(result$tes,log="y")
  plotlist_tree[[i]] <- recordPlot()
  # plotlist_ltt[[i]] <- recordPlot()
}
# file3 = paste("/Volumes/Liang/Research/data/loc_M0_pics/all.pdf",sep = "")
# file4 = paste("/Volumes/Liang/Research/data/loc_M0_pics/all1.pdf",sep = "")

# pdf(file = file3, width = 8.3, height = 11.7) 
# multiplot(plotlist = plotlist_tree, cols = 3)
# multiplot(plotlist = plotlist_ltt, cols = 3)
dev.off()

# l = as.list(c(1:15))
# ggsave(filename = file4, arrangeGrob(grobs = l))
# grid.arrange(plotlist_tree, ncol = 3, main = "Main title")



par(mfrow = c(5, 3))
par(cex = 0.6)
par(mar = c(0, 0, 0, 0), oma = c(4, 4, 0.5, 0.5))
par(tcl = -0.25)
par(mgp = c(2, 0.6, 0))
Index = matrix(1:15,5,3)
B = c(Index[1,],Index[2,],Index[3,],Index[4,],Index[5,])
for (i in 1:15) {
  j = B[i]
  if(j<5){
    file = paste("/Volumes/Liang/Research/data/loc_M0_pics/out",j,"sim.txt",sep = "")
  }
  else   file = paste("/Volumes/Liang/Research/data/loc_M0_pics/out",j,"sim.Rdata",sep = "")
  load(file = file)
  L = result$L
  phy = L2phylo(L)
  phy$tip.label = L[which(L[,4] == -1),5]
  
  # plot(phy)
  ltt.plot(result$tes)
  # plotlist_tree[[j]] <- recordPlot()
  
  # mtext(letters[i], side = 3, line = -1, adj = 0.1, cex = 0.6, col = "grey40")
  if (i %in% c(13, 14, 15))
  {  group = c(2,3,4)
  text = paste("Group",group[i-12])
    mtext(text = text, side = 1, cex = 0.5, line = 1.5,col = "grey20")
    # axis(1, col = "grey40", col.axis = "grey20", at = seq(5,10,15,20))
    
    }
  # axis(1, col = "grey40", col.axis = "grey20", at = seq(0.6,1.2, 0.2))
  if (i %in% c(1, 4,7,10,13))
  { M0 = rep(c(5,1,0.5,0.3,0),each = 3)
    text = paste("M0 =",M0[i])
    mtext(text = text, side = 2, cex = 0.5, line = 1.5,col = "grey20")
    # axis(2, col = "grey40", col.axis = "grey20", at = seq(0,50,100))
    
  }
 #  axis(2, col = "grey40", col.axis = "grey20", at = seq(0.6,1.2, 0.2))
 box(col = "grey60")
   }
mtext("Multiple Locations", side = 1, outer = TRUE, cex = 0.7, line = 2.5,col = "grey20")
mtext("Migration Rate", side = 2, outer = TRUE, cex = 0.7, line = 2.5,col = "grey20")



for (j in 1:15) {
  if(j<5){
    file = paste("/Volumes/Liang/Research/data/loc_M0_pics/out",j,"sim.txt",sep = "")
  }
  else   file = paste("/Volumes/Liang/Research/data/loc_M0_pics/out",j,"sim.Rdata",sep = "")
  load(file = file)
  print(result$tes)
}