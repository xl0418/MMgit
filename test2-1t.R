library(DDD)
source("/home/p274981/MM_sim.R")
library(expoRkit)
if (length(commandArgs(trailingOnly = TRUE)) == 0) {
  stop("Please supply an argument for this script, e [1,N]")
}

i <- commandArgs(trailingOnly = TRUE)[1]
print(i)
test2<-function(result,pars,endmc,i){
  brtsMM = branching.times(result[[1]])
  file = paste("analysis",i,".Rdata",sep = "")
  out<- DDD::dd_LR(brts = brtsMM,initparsoptDD = pars,initparsoptCR = c(0.5,0.2),missnumspec = 0,endmc = endmc,outputfilename = file)
}

parsN = c(1,1,0)
pars = c(0.8,0.4,20)
pars[3] = pars[1]*pars[3] /(pars[1]-pars[2])
endmc = 1000


  file1 = paste("/home/p274981/out",i,"sim.Rdata",sep = "")
  load(file = file1)
  test2(result = result, pars = pars,endmc = endmc,i=i)

