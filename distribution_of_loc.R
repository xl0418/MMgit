library(ggplot2)
library(gridExtra)
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
    L3 = NULL
    fre = NULL
    for (i in 1:100) {
      if(j ==2 & num <5){
        file = paste("/Volumes/Liang/Research/data/",j,"Modeltestgroup",num,"/out",i,"sim.txt",sep = "")
        
      }
      else file = paste("/Volumes/Liang/Research/data/",j,"Modeltestgroup",num,"/out",i,"sim.Rdata",sep = "")
      load(file = file)
      L = result$L
      L1 = which(L[,4] == -1)
      L2 = L[L1, 5]
      fre_num_loc = colsum[L2]
      
      fre_each = c(length(fre_num_loc[fre_num_loc == 1]),length(fre_num_loc[fre_num_loc == 2]),length(fre_num_loc[fre_num_loc == 3]),length(fre_num_loc[fre_num_loc == 4]))
      fre = rbind(fre,fre_each)
      L3 = c(L3,L2)
      fre_1 = data.frame(number = fre[,1])
      fre_1$loc = '1'
      fre_2 = data.frame(number = fre[,2])
      fre_2$loc = '2'
      fre_3 = data.frame(number = fre[,3])
      fre_3$loc = '3'
      fre_4 = data.frame(number = fre[,4])
      fre_4$loc = '4'
      fre_data = rbind(fre_1,fre_2,fre_3,fre_4)
    }
    if(j == 2) {
      loc_name = c('1','2')
      fre_data = fre_data[fre_data$loc %in% loc_name ,]
    }
    if(j == 3) {
      loc_name = c('1','2','3')
      fre_data = fre_data[fre_data$loc %in% loc_name ,]
    }
    if(j == 4) {
      loc_name = c('1','2','3','4')
      fre_data = fre_data[fre_data$loc %in% loc_name ,]
    }
    # ggplot(fre_data, aes(number, fill = loc)) + geom_histogram(alpha = 0.5, aes(y = ..density..), position = 'identity')
    count1 = count1+1
    
    if(count1 %in% c(1:5)){
      p[[count1]] <- ggplot(fre_data, aes(number, fill = loc)) +
        geom_histogram(binwidth=.5, alpha=.5, position="identity")+
      # geom_density(alpha = 0.2)+ 
        theme(axis.title.x=element_blank(),axis.title.y=element_blank(),legend.position="none")

      
    }
    if( count1 == 11){
      p[[count1]] <- ggplot(fre_data, aes(number, fill = loc)) +# geom_density(alpha = 0.2)+
        geom_histogram(binwidth=.5, alpha=.5, position="identity")+
        theme(axis.title.x=element_blank(),axis.title.y=element_blank(),legend.position = c(0.85,0.65))
      
    }
  
    else{
      p[[count1]] <- ggplot(fre_data, aes(number, fill = loc)) + # geom_density(alpha = 0.2)+ 
        geom_histogram(binwidth=.5, alpha=.5, position="identity")+
        theme(axis.title.x=element_blank(),axis.title.y=element_blank(),legend.position="none")
    }
    L4 = colsum[L3]
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
}
multiplot(plotlist=p,cols=3,labs=list("Number", "Frequency"),labpos=list(c(0.5,0.01), c(0.01,0.5)))
# mtext("Label", side = 1, outer = TRUE, cex = 0.7, line = 2.5,col = "grey20")
# mtext("Frequency", side = 2, outer = TRUE, cex = 0.7, line = 2.5,col = "grey20")
