library(DDD)
library(plotly)
library(MASS)
library(plot3D)
library(scatterplot3d)
library(lattice)
library(TeachingDemos)

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
}
else if(x==2){
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
  data <- data[-1,]
  time = data[order(data[,2]),2]
  timeu = time[-1,]
  data_lineage = timeu
  for(i in 1:100){
    file = paste(getwd(),"/Dropbox/R/cluster/Modeltestgroup",num,"/out",i,"sim.Rdata",sep = "")
    load(file = file)
    # pdffilename = paste(getwd(),"/Dropbox/R/cluster/Modeltestgroup",num,"/tree",i,".pdf",sep = "")
    # pdf(pdffilename,paper = "a4r", width = 29.7, height = 21)
    L = result$L
    tes = L2phylo(L,dropextinct = T)
    brts= -unname(sort(branching.times(tes),decreasing = T))
    M1 = match(brts,timeu)
    M1[1] = 1
    M11 = diff(M1)
    M13 = length(timeu)-max(M1)+1
    M12 = c(M11,M13)
    N1 = rep(2:(length(brts)+1),M12)
    data_lineage = cbind(data_lineage,N1)
  }
  x = data_lineage[,1]
  y = c(1:100)
  z = data_lineage[,2:101]
  scatter3D(x, y, z,surf = TRUE, clab = c("Number of", "lineages"))
 
    }
}