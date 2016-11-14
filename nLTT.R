library(ggplot2)
library(gridExtra)
library(DAISIE)
library(nLTT)
source('~/Dropbox/R/MMgit/multiplot.R')
source('~/Dropbox/R/MMgit/event_matrix.R', echo=TRUE)
p = list()
par(mfrow = c(3, 5))
par(cex = 0.6)
par(mar = c(0, 0, 0, 0), oma = c(4, 4, 0.5, 0.5))
par(tcl = -0.25)
par(mgp = c(2, 0.6, 0))
count1 = 1
count2 = 0
group = rep(c(2,3,4),each = 5)
for(j in c(2,3,4)){
  loc = event_matrix(j)
  colsum = colSums(loc)
  count2 = count2 +1
  print(count2)
  if(j == 2) num_group = c(1:5)
  else num_group = c(10:6)
  for(num in num_group){
    L3 = NULL
    fre = NULL
    L = list()
    data = NA
    
    for (i in 1:100) {
      file = paste("/Volumes/Liang/Research/data/",j,"Modeltestgroup",num,"/out",i,"sim.Rdata",sep = "")
    
      load(file = file)
      
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
    file = paste("/Volumes/Liang/Research/data/",j,"Modeltestgroup",num,"/out",i,"sim.Rdata",sep = "")
      
       load(file = file)
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
    z = data_lineage[,2:101]
    
    data_average_z <- apply(z, 1, median)
    data_q0.025_z <- apply(z, 1 , quantile, 0.025)
    data_q0.25_z <- apply(z, 1, quantile, 0.25)
    data_q0.75_z <- apply(z, 1, quantile, 0.75)
    data_q0.975_z <- apply(z, 1, quantile, 0.975)
    data_lower_z <- apply(z,1,min)
    data_upper_z <- apply(z,1,max)
    
    lineage_stat = cbind(x,data_average_z,data_q0.025_z,data_q0.25_z,data_q0.75_z,data_q0.975_z,data_lower_z,data_upper_z)
    colnames(lineage_stat) = c("time", "median","0.025","0.25","0.75","0.975","min","max")
   time = min(lineage_stat[,1])
   df_lineage = data.frame(lineage_stat)
   df_min_max = data.frame(id = "min_max", value = 1, x = c(df_lineage$time,rev(df_lineage$time)), y = c(df_lineage$min,rev(df_lineage$max)))
   df_0025 = data.frame(id = "0025", value = 2, x = c(df_lineage$time,rev(df_lineage$time)), y = c(df_lineage$X0.025,rev(df_lineage$X0.975)))
   df_025 = data.frame(id = "025", value = 3, x = c(df_lineage$time,rev(df_lineage$time)), y = c(df_lineage$X0.25,rev(df_lineage$X0.75)))
   df_lineage_all = rbind(df_min_max,df_025,df_0025)
   
   if(count1 %in% c(5,10,15)){
     x_lab = paste("group =",group[count1])
   p[[count1]] <- ggplot(df_min_max, aes(x = x, y = y)) +
     theme(axis.title.y=element_blank(),legend.position="none")+
   geom_polygon(data = df_min_max, aes(  group = id),fill = "light gray", alpha = 0.8)+
     geom_polygon(data = df_0025, aes( group = id),fill = "dark gray", alpha = 0.8)+
     geom_polygon(data = df_025, aes( group = id), fill = "gray27", alpha = 0.8)+ylab("")+xlab(x_lab)+
     geom_line(data = df_lineage, aes(time, median), col = "black")
   }
   else
   {
     p[[count1]] <- ggplot(df_min_max, aes(x = x, y = y)) +
       theme(axis.title.x=element_blank(),axis.text.x = element_blank() ,axis.title.y=element_blank(),legend.position="none")+
       geom_polygon(data = df_min_max, aes(  group = id),fill = "light gray", alpha = 0.8)+
       geom_polygon(data = df_0025, aes( group = id),fill = "dark gray", alpha = 0.8)+
       geom_polygon(data = df_025, aes( group = id), fill = "gray27", alpha = 0.8)+
       geom_line(data = df_lineage, aes(time, median), col = "black")
   }
     count1 = count1+1
    }
   
  }
multiplot(plotlist=p,cols=3,labs=list("Time", "Number of lineages"),labpos=list(c(0.5,0.01), c(0.01,0.5)))

