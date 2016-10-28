

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
for (i in 1:100) {
  if(j ==2 & num <5){
    file = paste("/Volumes/Liang/Research/data/",j,"Modeltestgroup",num,"/out",i,"sim.txt",sep = "")
    
  }
  else file = paste("/Volumes/Liang/Research/data/",j,"Modeltestgroup",num,"/out",i,"sim.Rdata",sep = "")
  load(file = file)
  L = result$L
  L1 = which(L[,4] == -1)
  L2 = L[L1, 5]
  L3 = c(L3,L2)
 
}
   
    
    L4 = colsum[L3]
    count1 = count1+1
    # distribution of labels
    # hist(L3,col = "lightblue", border = "pink",main = "",breaks = seq(1,range(L3)[2],by = 0.5))
    
    #distribution of number of locations
    hist(L4,col = "lightblue", border = "pink",main = "",breaks = seq(1,j,by = 0.5))
   
     if (j == 4)
    {  M0 = c(5,1,0.5,0.3,0)
    text = paste("M0 =",M0[11-num])
    mtext(text = text, side = 1, cex = 0.5, line = 1.5,col = "grey20")

    }
    if (count1 %in% c(1,6,11))
    { group = c(2,3,4)
    text = paste("Group ",group[count2])
    mtext(text = text, side = 2, cex = 0.5, line = 1.5,col = "grey20")
    }
  }
}
mtext("Label", side = 1, outer = TRUE, cex = 0.7, line = 2.5,col = "grey20")
mtext("Frequency", side = 2, outer = TRUE, cex = 0.7, line = 2.5,col = "grey20")