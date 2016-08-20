# library(ggplot2)
# library(reshape)
library(DDD)
source("/Users/mac/Dropbox/R/MigrationModelsim/number2binary.R")
MM_sim<-function(parsN,age=100,pars ,M = 1,  lambdaab0=0.2, M0=(M == 1)*0.4,model=1){
if(sum(parsN) ==2 | sum(parsN)==1){
  done = 0
  while(done == 0)
  {
  lambda0 = pars[1]
  mu0 = pars[2]
  K = pars[3]
  Ka= lambda0*K/(lambda0 - mu0)
  Kb = Ka
  
  t <- 0
  # Na : number of species in location A
  # Nb : number of species in location B
  # Nab: number of species in location AB
  # age : time
  # Ki : carrying capacity of location i
  # lambda0 : initial speciation rate
  # lambdaab0 : allopatric speciation parameter 
  # mu0 : initial extinction rate  
  # mui : initial extinction rate in location i
  mua=mu0
  mub=mu0
  # coni : contraction rate in location i
  cona=mua
  conb=mub 
  # M0 : initial migration rate
  # string of events
  #B<- c("spec A","spec B","allo spec","ext A","ext B", "ext AB", "mig AB", "mig BA", "con A", "con B", "sym spec A", "sym spec B")
  B<- c(1:11)
  Na = parsN[1]
  Nb = parsN[2]
  Nab = parsN[3]
  N=Na+Nb+Nab
  Ntable = cbind(Na,Nb,Nab)
  i=0
  # L : Ltable used in L2phylo function of DDD package 
  # L = data structure for lineages,
  # . L[,1] = branching times
  # . L[,2] = index of parent species
  # . L[,3] = index of daughter species
  # . L[,4] = time of extinction
  # j = index running through L
  
  loctable= NULL
  
  L = NULL
  for(j in 1:N){
    loc = which(Ntable[1,]!=0)
    L = rbind(L, c(0, 1-j, (-1)^j*j, -1, loc[1]))
    L = matrix(L, ncol = 5)
    if(is.element(loc,1:2)[1]){
      loc1= matrix(0,1,2)
      loc1[1,loc[1]] = 1 
    }
    else {
      loc1 = matrix(1,1,2)
    }
    loctable = rbind(loctable, loc1)
    linlist = cbind(L[,3], L[,5])
    Ntable[1,loc[1]] = Ntable[1,loc[1]] -1 
    newL=j
  }
  # print(loctable)
  
  Ntable = cbind(Na,Nb,Nab)
  while(t[i+1]< age){
    i<-i+1
   # print(t[i])
    Na = Ntable[,1]
    Nb = Ntable[,2]
    Nab = Ntable[,3]
    # speciation event in A
    if(model==1){
    lambdaNa=max(lambda0*(1-(Na[i]+Nab[i])/Ka),0)
    }
    if(model==0) {lambdaNa=lambda0}
    birth_event_A=lambdaNa*Na[i]
    # speciation event in B
    if(model==1){
    lambdaNb=max(lambda0*(1-(Nb[i]+Nab[i])/Kb),0)
    }
    if(model==0) {lambdaNb=lambda0}
    birth_event_B=lambdaNb*Nb[i]
    # Extinction event in A
    death_event_A=mua*Na[i]
    # Extinction event in B
    death_event_B=mub*Nb[i]
    # Migration 
    if(model==1){
    Mab=max(M0*(1-(Nb[i]+Nab[i])/Kb),0)
    Mba=max(M0*(1-(Na[i]+Nab[i])/Ka),0)
    }
    if(model==0){
      Mab=M0
      Mba=M0
    }
   # print(paste("Mab+Mba =",Mab+Mba))
    migration_event_AB=Mab*Na[i]
    migration_event_BA=Mba*Nb[i]
    # Allopatric speciation event in AB
     if (Mab+Mba == 0) {lambdaab= 0}
     else lambdaab=max(lambdaab0/(Mab+Mba),0)
    birth_event_AB=lambdaab*Nab[i]
    # speciation event in AB to A
    birth_event_ABa=lambdaNa*Nab[i]
    # speciation event in AB to B
    birth_event_ABb=lambdaNb*Nab[i]
    #Contraction of A
    contraction_event_A=cona*Nab[i]
    #Contraction of A
    contraction_event_B=conb*Nab[i]
    # Probabilities for each event
    probs= c(birth_event_A,birth_event_B,birth_event_ABa,birth_event_ABb,death_event_A,death_event_B,contraction_event_A,contraction_event_B,birth_event_AB,migration_event_AB,migration_event_BA)
    
    #Total rate
    TR=sum(probs)
    if(TR==0) break
    else{
      A<-DDD::sample2(B,1, prob = probs)
      t[i+1]=t[i]+rexp(1,rate=TR)
      if(t[i+1]>age) break
      loc2 = matrix(0,1,2)
      # Sympatric speciation 
      if (is.element(A,1:4)){
        b1<-number2binary(A,4)
        b2<-2-b1[4]
        Ntable=rbind(Ntable,Ntable[i,])
        Ntable[i+1,b2] = Ntable[i,b2]+1
        newL = newL + 1;
        list0 = matrix(linlist,ncol = 2)
        b3 <- sign(A-2)+2
       # print(b3)
        list1 = linlist[list0[,2]== b3]
       # print(list1)
        list2 = matrix(list1, ncol = 2)
       # print(list2)
        linlist1 = list2[,1]
       # print(linlist1)
        ranL= DDD::sample2(linlist1,1)
        L = rbind(L,c(t[i+1],ranL,sign(ranL) * newL,-1,b2))
        linlist = rbind(linlist,c(sign(ranL) * newL,b2)) 
        loc2[1,b2] = 1
        loctable = rbind(loctable, loc2)
      }
      
      #Extinction
      else if(is.element(A,5:8)) {
        b1<-number2binary(A-4,4)
        b2<-2-b1[4]
        Ntable=rbind(Ntable,Ntable[i,])
        b3 <- sign(A-6)+2
        list0 = matrix(linlist,ncol = 2)
        list1 = linlist[list0[,2]== b3]
        list2 = matrix(list1, ncol = 2)
        linlist1 = list2[,1]
        ranL= DDD::sample2(linlist1,1)
       if(A==5 | A==6){
         Ntable[i+1,b3] = Ntable[i,b3]-1
        L[abs(ranL),4] = t[i+1]
        if(length(L[L[,4]== -1]) == 0) break
        else {
          v = which(linlist[,1] == ranL)
          linlist = linlist[-v,,drop=FALSE]
          linlist = linlist[order(linlist[,1]),]
          linlist = matrix(linlist,ncol=2)
        }
       }
        else{
          Ntable[i+1,b3] = Ntable[i,b3]-1
          Ntable[i+1,3-b2] = Ntable[i,3-b2]+1
          L[abs(ranL),5] <- 9-A
          v = which(linlist[,1] == ranL)
          linlist[v,2] = 9-A
          linlist = linlist[order(linlist[,1]),]
          linlist = matrix(linlist,ncol=2)
          loctable[abs(ranL),A-6] = 0
        }
      }
      
      
      # Allopatric speciation in AB
      else if(A== 9) {
        Ntable=rbind(Ntable,Ntable[i,])
        Ntable[i+1,1] = Ntable[i,1]+1
        Ntable[i+1,2] = Ntable[i,2]+1
        Ntable[i+1,3] = Ntable[i,3]-1
        newL = newL + 1
        list0 = matrix(linlist,ncol = 2)
        list1 = linlist[list0[,2]==3]
        list2 = matrix(list1, ncol = 2)
        linlist1 = list2[,1]
        ranL= DDD::sample2(linlist1,1)
        L[abs(ranL),5] <- 1
        loctable[abs(ranL), 2] = 0
        v = which(linlist[,1] == ranL)
        linlist[v,2] = 1
        L = rbind(L,c(t[i+1],ranL,sign(ranL) * newL,-1,2))
        linlist = rbind(linlist,c(sign(ranL) * newL,2))  
        loc2[1,2] = 1 
        loctable = rbind(loctable,loc2)
      }
     
      #Migration
      else {
        b1<-number2binary(A-9,4)
        b2<-2-b1[4]
        Ntable=rbind(Ntable,Ntable[i,])
        Ntable[i+1,b2]=Ntable[i,b2]-1
        Ntable[i+1,3]=Ntable[i,3]+1
        linlist = matrix(linlist,ncol = 2)
        list1 = linlist[linlist[,2]==b2]
        list2 = matrix(list1, ncol = 2)
        linlist1 = list2[,1]
        ranL1 = DDD::sample2(linlist1,1)
        L[abs(ranL1),5] <- 3
        v = which(linlist[,1] == ranL1)
        linlist[v,2] = 3
        linlist = linlist[order(linlist[,1]),]
        linlist = matrix(linlist,ncol=2)
        loctable[abs(ranL1),] = c(1,1)
      }
      if(sum(linlist[,1] < 0) == 0 | sum(linlist[,1] > 0) == 0) 
        { #print("break")
        break
      }
      # print(Na[i+1])
      N[i+1]=sum(Ntable[i+1,])
      }
  }
  
  if(sum(linlist[,1] < 0) == 0 | sum(linlist[,1] > 0) == 0)
  {
    done = 0
  } else {
    done = 1
  }
  }  
#   print(Ntable)
#   print(linlist)
  
#   print(head(linlist))
#   print(head(L))
  Ltable = cbind(L, loctable)
  # print(head(Ltable))
  
  if(dim(L)[1]==1)
    print(paste("It goes extinct at the begining!"))
  else{
  L[,1]=age-L[,1]
  notmin1 = which(L[,4] != -1)
  L[notmin1,4] = age - c(L[notmin1,4])
  L[which(L[,4] == age + 1),4] = -1
 # print(L)
  tes = L2phylo(L,dropextinct = T)
  tas = L2phylo(L,dropextinct = F)
  out = list(tes = tes,tas = tas,L = L)
#    file = paste("/Users/mac/Dropbox/R/MigrationModelsim/ModelTest/",parsN[1],parsN[2],parsN[3],"&",pars[1],pars[2],pars[3],"&",age,"sim.txt")
#    save(out,file = file)
   return(out)
  }
  }
  else print(paste("please input 1 or 2 species in total!"))
}
print(paste("Please input Na, Nb, Nab, age with no more than 2 species in total!"))