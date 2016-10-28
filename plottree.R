library(DDD)
library(MASS)
library(rgl)
library(stringr)
plottree <- function(x,num_loc,num){
  # LTT plot
if(x == 1){
  for(i in 1:100){
  file = paste("/Volumes/Liang/Research/data/",num_loc,"Modeltestgroup",num,"/out",i,"sim.Rdata",sep = "")
  # file = paste("/Users/mac/Dropbox/R/cluster/Modeltestgroupmu/out",i,"sim.Rdata",sep = "")
  load(file = file)
  pdffilename = paste("/Volumes/Liang/Research/data/",num_loc,"Modeltestgroup",num,"/tree",i,".pdf",sep = "")
  # pdffilename = paste("/Users/mac/Dropbox/R/cluster/Modeltestgroupmu/tree",i,".pdf",sep = "")
  pdf(pdffilename,paper = "a4r", width = 29.7, height = 21)
  L = result$L
  phy = L2phylo(L)
  phy$tip.label = L[which(L[,4] == -1),5]
  
  # loc = result$loctable
  # loc_lable = paste(loc[,1],loc[,2],loc[,3])
  # phy$tip.label = loc_lable
  plot(phy)
  ltt.plot(result$tes,log="y")
  try(dev.off())

  }
}
  
  # 3D surface for simulations
else if(x==2){
  data = NA
  for(i in 1:100){
    # file = paste("/Volumes/Liang/Research/data/Modeltestgroup",num,"/out",i,"sim.Rdata",sep = "")
    file = paste("/Users/mac/Dropbox/R/cluster/Modeltestgroupmu/out",i,"sim.Rdata",sep = "")
    load(file = file)
   # pdffilename = paste(getwd(),"/Dropbox/R/cluster/Modeltestgroup",num,"/tree",i,".pdf",sep = "")
   # pdf(pdffilename,paper = "a4r", width = 29.7, height = 21)
    L = result$L
    tes = L2phylo(L,dropextinct = T)
    brts= -unname(sort(branching.times(tes),decreasing = T))
    # print(length(brts))
    data0 = cbind(i,brts,c(2:(length(brts)+1)))
    data = rbind(data, data0)
   # try(dev.off())
    
  }
  data <- data[-1,]
  time = data[order(data[,2]),2]
  timeu = unique(time)
  data_lineage = timeu
  for(i in 1:100){
    # file = paste("/Volumes/Liang/Research/data/Modeltestgroup",num,"/out",i,"sim.Rdata",sep = "")
    file = paste("/Users/mac/Dropbox/R/cluster/Modeltestgroupmu/out",i,"sim.Rdata",sep = "")
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
  zlim <- range(z)
  zlen <- zlim[2] - zlim[1] + 1
  colorlut <- terrain.colors(zlen,alpha=0) # height color lookup table
  col <- colorlut[ z-zlim[1]+1 ] # assign colors to heights for each point
  open3d()
  rgl.surface(x, y, z, color=col, alpha=0.75, back="lines")
  axes3d()
  }
}