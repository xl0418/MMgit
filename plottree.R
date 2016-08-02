library(DDD)
library(plotly)
library(MASS)
library(plot3D)
plottree <- function(x,num){
if(x == 1){
  for(i in 1:100){
  file = paste(getwd(),"/Dropbox/R/cluster/Modeltestgroup",num,"/out",i,"sim.Rdata",sep = "")
  load(file = file)
  pdffilename = paste(getwd(),"/Dropbox/R/cluster/Modeltestgroup",num,"/tree",i,".pdf",sep = "")
  pdf(pdffilename,paper = "a4r", width = 29.7, height = 21)
  L = result$L
  phy = L2phylo(L)
  phy$tip.label = L[which(L[,4] == -1),5]
  plot(phy)
  ltt.plot(result$tes)
  
  try(dev.off())

}
else (x==2){
  data = NA
  for(i in 1:100){
    file = paste(getwd(),"/Dropbox/R/cluster/Modeltestgroup",num,"/out",i,"sim.Rdata",sep = "")
    load(file = file)
   # pdffilename = paste(getwd(),"/Dropbox/R/cluster/Modeltestgroup",num,"/tree",i,".pdf",sep = "")
   # pdf(pdffilename,paper = "a4r", width = 29.7, height = 21)
    L = result$L
    tes = L2phylo(L,dropextinct = T)
    brts= -unname(sort(branching.times(tes),decreasing = T))
    print(length(brts))
    data0 = cbind(i,brts,c(2:(length(brts)+1)))
    data = rbind(data, data0)
   # try(dev.off())
   
    
  }
  kd <- with(MASS::geyser, MASS::kde2d(duration, waiting, n = 50))
  with(kd, plot_ly(x = data[,1], y = data[,2], z = data[,3], type = "surface"))
  
  
  plot_ly(z = data, type = "surface")
  
  
  x = data[,1]
  y = data[,2]
  z = data[,3]
  scatter3D(x, y, z, clab = c("Number of", "lineages"))
  
}
  
  
}
}