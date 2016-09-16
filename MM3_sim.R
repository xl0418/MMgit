# library(ggplot2)
# library(reshape)
library(DDD)
source("/Users/mac/Dropbox/R/MigrationModelsim/number2binary.R")
MM_sim<-function(parsN,age=100,pars ,M = 1,  lambda_allo=0.2, M0=(M == 1)*0.4){
if(sum(parsN) ==2 | sum(parsN)==1){
  done = 0
  while(done == 0)
  {
  lambda0 = pars[1]
  mu0 = pars[2]
  K = pars[3]
  Ka= lambda0*K/(lambda0 - mu0)
  Kb = Ka
  Kc = Ka
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
  muc=mu0
  
  lambdaab0 = lambda_allo
  lambdaac0 = lambda_allo
  lambdabc0 = lambda_allo
  # M0 : initial migration rate
  # string of events
  #B<- c("spec A","spec B","allo spec","ext A","ext B", "ext AB", "mig AB", "mig BA", "con A", "con B", "sym spec A", "sym spec B")
  B<- c(1:36)
  Na = parsN[1]
  Nb = parsN[2]
  Nc = parsN[3]
  N=Na+Nb+Nc
  Ntable = cbind(Na,Nb,Nc)
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
  
  Nab = Nac = Nbc =Nabc = 0
  Ntable = cbind(Na,Nb,Nc,Nab,Nac,Nbc,Nabc)
  while(t[i+1]< age){
    i<-i+1
   # print(t[i])
    Na = Ntable[,1]
    Nb = Ntable[,2]
    Nc = Ntable[,3]
    Nab = Ntable[,4]
    Nac = Ntable[,5]
    Nbc = Ntable[,6]
    Nabc = Ntable[,7]
    # speciation event in A
    lambdaNa=max(lambda0*(1-(Na[i]+Nab[i]+Nac[i]+Nabc[i])/Ka),0)
    birth_event_A=lambdaNa*Na[i]
    birth_event_ABa=lambdaNa*Nab[i]
    birth_event_ACa=lambdaNa*Nac[i]
    birth_event_ABCa=lambdaNa*Nabc[i]
    
    # speciation event in B
    lambdaNb=max(lambda0*(1-(Nb[i]+Nab[i]+Nbc[i]+Nabc[i])/Kb),0)
    birth_event_B=lambdaNb*Nb[i]
    birth_event_BCb=lambdaNa*Nbc[i]
    birth_event_BAb=lambdaNa*Nab[i]
    birth_event_ABCb=lambdaNa*Nabc[i]
    
    # speciation event in C
    lambdaNc=max(lambda0*(1-(Nc[i]+Nac[i]+Nbc[i]+Nabc[i])/Kc),0)
    birth_event_C=lambdaNc*Nc[i]
    birth_event_CAc=lambdaNc*Nac[i]
    birth_event_CBc=lambdaNc*Nbc[i]
    birth_event_ABCc=lambdaNc*Nabc[i]
    
    # Extinction event in A
    death_event_A=mua*Na[i]
    contraction_event_ABa=mua*Nab[i]
    contraction_event_ACa=mua*Nac[i]
    contraction_event_ABCa=mua*Nabc[i]
    
    
    # Extinction event in B
    death_event_B=mub*Nb[i]
    contraction_event_BCb=mua*Nbc[i]
    contraction_event_BAb=mua*Nab[i]
    contraction_event_ABCb=mua*Nabc[i]
    
    # Extinction event in C
    death_event_C=muc*Nc[i]
    contraction_event_CAc=muc*Nac[i]
    contraction_event_CBc=muc*Nbc[i]
    contraction_event_ABCc=muc*Nabc[i]
    
    
    # Migration 
    
    Mab=max(M0*(1-(Nb[i]+Nab[i]+Nbc[i]+Nabc[i])/Kb),0)
    migration_event_AB=Mab*Na[i]
    Mac=max(M0*(1-(Nc[i]+Nac[i]+Nbc[i]+Nabc[i])/Kc),0)
    migration_event_AC=Mac*Na[i]
    
    Mba=max(M0*(1-(Na[i]+Nab[i]+Nac[i]+Nabc[i])/Ka),0)
    migration_event_BA=Mba*Nb[i]
    Mbc=max(M0*(1-(Nc[i]+Nbc[i]+Nac[i]+Nabc[i])/Kc),0)
    migration_event_BC=Mbc*Nb[i]
    
    Mca=max(M0*(1-(Na[i]+Nab[i]+Nac[i]+Nabc[i])/Ka),0)
    migration_event_CA=Mca*Nc[i]
    Mcb=max(M0*(1-(Nb[i]+Nbc[i]+Nab[i]+Nabc[i])/Kb),0)
    migration_event_CB=Mcb*Nc[i]
    
    Mabc=max(M0*(1-(Nc[i]+Nbc[i]+Nac[i]+Nabc[i])/Kc),0)
    migration_event_ABC=Mabc*Nab[i]
    Macb=max(M0*(1-(Nb[i]+Nbc[i]+Nab[i]+Nabc[i])/Kb),0)
    migration_event_ACB=Macb*Nac[i]
    Mbca=max(M0*(1-(Na[i]+Nab[i]+Nac[i]+Nabc[i])/Ka),0)
    migration_event_BCA=Mbca*Nbc[i]
    
    # Allopatric speciation event in AB
     if (Mab+Mba == 0) {lambdaab= 0}
     else lambdaab=max(lambdaab0/(Mab+Mba),0)
    birth_event_AB=lambdaab*Nab[i]
    
    if (Mac+Mca == 0) {lambdaac= 0}
    else lambdaac=max(lambdaac0/(Mac+Mca),0)
    birth_event_CA=lambdaac*Nac[i]
    
    if (Mbc+Mcb == 0) {lambdabc= 0}
    else lambdabc=max(lambdabc0/(Mbc+Mcb),0)
    birth_event_BC=lambdabc*Nbc[i]
   
   
    # Probabilities for each event
    # birth_event_A = 1
    # birth_event_B = 2
    # birth_event_C = 3 
    # birth_event_ABa = 4
    # birth_event_BCb = 5
    # birth_event_CAc = 6
    # birth_event_ACa = 7
    # birth_event_BAb = 8
    # birth_event_CBc = 9
    # birth_event_ABCa = 10
    # birth_event_ABCb = 11
    # birth_event_ABCc = 12
    # death_event_A = 13
    # death_event_B = 14
    # death_event_C = 15
    # contraction_event_ABa = 16
    # contraction_event_BCb = 17
    # contraction_event_CAc = 18
    # contraction_event_ACa = 19
    # contraction_event_BAb = 20
    # contraction_event_CBc = 21
    # contraction_event_ABCa = 22
    # contraction_event_ABCb = 23
    # contraction_event_ABCc = 24
    # birth_event_AB = 25
    # birth_event_BC = 26
    # birth_event_CA = 27
    # migration_event_AB = 28
    # migration_event_AC = 29
    # migration_event_BA = 30
    # migration_event_BC = 31
    # migration_event_CA = 32
    # migration_event_CB = 33
    # migration_event_ABC = 34
    # migration_event_ACB = 35
    # migration_event_BCA = 36
    
    
    probs= c(birth_event_A,birth_event_B,birth_event_C,birth_event_ABa,birth_event_BCb,birth_event_CAc,birth_event_ACa,birth_event_BCb,birth_event_CBc,birth_event_ABCa,birth_event_ABCb,birth_event_ABCc,death_event_A,death_event_B,death_event_C,contraction_event_ABa,contraction_event_BCb,contraction_event_CAc,contraction_event_ACa,contraction_event_BAb,contraction_event_CBc,contraction_event_ABCa,contraction_event_ABCb,contraction_event_ABCc,birth_event_AB,birth_event_BC,birth_event_CA,migration_event_AB,migration_event_AC,migration_event_BA,migration_event_BC,migration_event_CA,migration_event_CB,migration_event_ABC,migration_event_ACB,migration_event_BCA)
    
    #Total rate
    TR=sum(probs)
    if(TR==0) break
    else{
      A<-DDD::sample2(B,1, prob = probs)
      t[i+1]=t[i]+rexp(1,rate=TR)
      if(t[i+1]>age) break
      loc2 = matrix(0,1,2)
      # Sympatric speciation 
      if (is.element(A,1:12)){
        B = c(1,2,3,4,4,5,5,6,6,7,7,7)
        b1<-A%%3
        if(b1 == 0) b1 = 3
        Ntable=rbind(Ntable,Ntable[i,])
        Ntable[i+1,b1] = Ntable[i,b1]+1
        newL = newL + 1;
        list0 = matrix(linlist,ncol = 2)
        b3 <- B[A]
       # print(b3)
        list1 = linlist[list0[,2]== b3]
       # print(list1)
        list2 = matrix(list1, ncol = 2)
       # print(list2)
        linlist1 = list2[,1]
       # print(linlist1)
        ranL= DDD::sample2(linlist1,1)
        L = rbind(L,c(t[i+1],ranL,sign(ranL) * newL,-1,b1))
        linlist = rbind(linlist,c(sign(ranL) * newL,b1)) 
        # loc2[1,b1] = 1
        # loctable = rbind(loctable, loc2)
      }
      
      #Extinction
      else if(is.element(A,13:24)) {
        B = c(1,2,3,4,4,5,5,6,6,7,7,7)
        b1<-A%%3
        if(b1 == 0) b1 = 3
        Ntable=rbind(Ntable,Ntable[i,])
        b3 <- B[A-12]
        list0 = matrix(linlist,ncol = 2)
        list1 = linlist[list0[,2]== b3]
        list2 = matrix(list1, ncol = 2)
        linlist1 = list2[,1]
        ranL= DDD::sample2(linlist1,1)
       if(is.element(A,13:15)){
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
          B1 = c(0,0,0,2,1,1,3,3,2,6,5,4)
          b2 <- B1[A-12]
          Ntable[i+1,b3] = Ntable[i,b3]-1
          Ntable[i+1,b2] = Ntable[i,b2]+1
          L[abs(ranL),5] <- 9-A
          v = which(linlist[,1] == ranL)
          linlist[v,2] = b2
          linlist = linlist[order(linlist[,1]),]
          linlist = matrix(linlist,ncol=2)
          # loctable[abs(ranL),A-6] = 0
        }
      }
      
      
      # Allopatric speciation in AB
      else if(is.element(A,25:27)) {
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