library(ggplot2)
library(gridExtra)
source('~/Dropbox/R/MMgit/multiplot.R')

# dis_loc_plot = function(group, )
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
    L_age_loc_all = NULL
    for (i in 1:100) {
      file = paste("/Volumes/Liang/Research/data/",j,"Modeltestgroup",num,"/out",i,"sim.Rdata",sep = "")
      load(file = file)
      L = result$L
      L1 = which(L[,4] == -1)
      L2 = L[L1, 5]
      fre_num_loc = colsum[L2]
      L_age_loc = cbind(L[L1,1],fre_num_loc)
      L_age_loc_all = rbind(L_age_loc_all,L_age_loc)
    }
      fre_1 = data.frame(number = L_age_loc_all[which(L_age_loc_all[,2]==1),1])
      fre_1$loc = '1'
      if(length(which(L_age_loc_all[,2]==2)) == 0) fre_2 = NULL
     else{ fre_2 = data.frame(number = L_age_loc_all[L_age_loc_all[,2]==2,1])
      fre_2$loc = '2'
     }
      
      if(length(which(L_age_loc_all[,2]==3)) == 0) fre_3 = NULL
      else{ fre_3 = data.frame(number = L_age_loc_all[L_age_loc_all[,2]==3,1])
      fre_3$loc = '3'
      }
      
      if(length(which(L_age_loc_all[,2]==4)) == 0) fre_4 = NULL
      else{ fre_4 = data.frame(number = L_age_loc_all[L_age_loc_all[,2]==4,1])
      fre_4$loc = '4'
      }
     
      fre_data = rbind(fre_1,fre_2,fre_3,fre_4)
    
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
   
  }
}
multiplot(plotlist=p,cols=3,labs=list("Number", "Frequency"),labpos=list(c(0.5,0.01), c(0.01,0.5)))
# mtext("Label", side = 1, outer = TRUE, cex = 0.7, line = 2.5,col = "grey20")
# mtext("Frequency", side = 2, outer = TRUE, cex = 0.7, line = 2.5,col = "grey20")
