library(ggplot2)
library(gridExtra)
library(DAISIE)
library(nLTT)
source('~/Dropbox/R/MMgit/multiplot.R')
p = list()
par(mfrow = c(3, 5))
par(cex = 0.6)
par(mar = c(0, 0, 0, 0), oma = c(4, 4, 0.5, 0.5))
par(tcl = -0.25)
par(mgp = c(2, 0.6, 0))
count1 = 0
count2 = 0
for(j in c(2,3,4)){
  loc = event_matrix(j)
  colsum = colSums(loc)
  count2 = count2 +1
  print(count2)
  if(j == 2) num_group = c(1:5)
  else num_group = c(10:6)
  for(num in num_group){
    L_extant_all =NULL
    for (i in 1:100) {
      if(j ==2 & num <5){
        file = paste("/Volumes/Liang/Research/data/",j,"Modeltestgroup",num,"/out",i,"sim.txt",sep = "")
        
      }
      else file = paste("/Volumes/Liang/Research/data/",j,"Modeltestgroup",num,"/out",i,"sim.Rdata",sep = "")
      
      load(file = file)
      
      L = result$L
      L_extant = L[which(L[,4]==-1),]
      L_extant_all = rbind(L_extant_all,L_extant)
      
    }
    L_extant_all[,1] = 20 - L_extant_all[,1]
    age = c(1:20)
    for(i in age){
    L_extant_all_groups = L_extant_all[which(L_extant_all[,1] <age[i]),]
    x = data_lineage[,1]
    z = data_lineage[,2:101]
    
    data_average_z <- apply(z, 1, median)
    data_q0.025_z <- apply(z, 1 , quantile, 0.025)
    data_q0.25_z <- apply(z, 1, quantile, 0.25)
    data_q0.75_z <- apply(z, 1, quantile, 0.75)
    data_q0.975_z <- apply(z, 1, quantile, 0.975)
    data_lower_z <- apply(z,1,min)
    data_upper_z <- apply(z,1,max)
    
    lineage_stat = cbind(x,data_average_z,data_q0.025_z,data_q0.25_z,data_q0.75_z,data_q0.975_z,data_lower_z,data_upper_z)
   
    }
     colnames(lineage_stat) = c("time", "median","0.025","0.25","0.75","0.975","min","max")
    time = min(lineage_stat[,1])
    suppressWarnings(plot(NULL, NULL, xlim = rev(c(0, time)), 
                          ylim = c(1, max(data_upper_z)), ylab = "No of species ", 
                          bty = "l", xaxs = "i", xlab = "Time before present", 
                          main = "nLTT", log = "y", cex.lab = 1.5, 
                          cex.main = 1.2, cex.axis = 1.2))
    polygon(c(lineage_stat[, "time"], rev(lineage_stat[, "time"])), c(lineage_stat[, "min"], rev(lineage_stat[, "max"])), col = "light grey", border = NA,angle = 45)
    
    polygon(c(lineage_stat[, "time"], rev(lineage_stat[, "time"])), c(lineage_stat[, "0.025"] , rev(lineage_stat[, "0.975"] )), col = "dark grey", border = NA,angle = 45)
    polygon(c(lineage_stat[, "time"], rev(lineage_stat[, "time"])), c(lineage_stat[, "0.25"] , rev(lineage_stat[, "0.75"] )), col = "light blue", border = NA,angle = 45)
    lines(lineage_stat[, "time"], data_average_z , lwd = 2, col = "dodgerblue1")
    
    
    # ggplot(fre_data, aes(number, fill = loc)) + geom_histogram(alpha = 0.5, aes(y = ..density..), position = 'identity')
    count1 = count1+1
  }
  
  # distribution of labels
  # hist(L3,col = "lightblue", border = "pink",main = "",breaks = seq(1,range(L3)[2],by = 0.5))
  
  #distribution of number of locations
  # hist(L4,col = "lightblue", border = "pink",main = "",breaks = seq(1,j,by = 0.5))
  
  #  if (j == 4)
  # {  M0 = c(5,1,0.5,0.3,0)
  # text = paste("M0 =",M0[11-num])
  # mtext(text = text, side = 1, cex = 0.5, line = 1.5,col = "grey20")
  # 
  # }
  # if (count1 %in% c(1,6,11))
  # { group = c(2,3,4)
  # text = paste("Group ",group[count2])
  # mtext(text = text, side = 2, cex = 0.5, line = 1.5,col = "grey20")
  # }
}
multiplot(plotlist=L,cols=3,labs=list("Number", "Frequency"),labpos=list(c(0.5,0.01), c(0.01,0.5)))
# mtext("Label", side = 1, outer = TRUE, cex = 0.7, line = 2.5,col = "grey20")
# mtext("Frequency", side = 2, outer = TRUE, cex = 0.7, line = 2.5,col = "grey20")
