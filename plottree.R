for(i in 1:100){
  file = paste(getwd(),"/Dropbox/R/cluster/Modeltestgroup28/out",i,"sim.Rdata",sep = "")
  load(file = file)
  pdffilename = paste(getwd(),"/Dropbox/R/cluster/Modeltestgroup28/tree",i,".pdf",sep = "")
  pdf(pdffilename,paper = "a4r", width = 29.7, height = 21)
  L = result$L
  phy = L2phylo(L)
  phy$tip.label = L[which(L[,4] == -1),5]
  plot(phy)
  ltt.plot(result$tes)
  
  try(dev.off())

}
