# library(ggplot2)
# library(reshape)
library(DDD)
source("/Users/mac/Dropbox/R/MMgit/Nindex.R")
source("/Users/mac/Dropbox/R/MMgit/event_matrix.R")
MM4_sim<-function(n,parsN,age=100,pars ,M = 1,  lambda_allo0=0.2, M0=(M == 1)*10){
if(sum(parsN) ==2 | sum(parsN)==1){
  done = 0
  while(done == 0)
  {
  lambda0 = pars[1]
  mu0 = pars[2]
  K = pars[3]
  K= lambda0*K/(lambda0 - mu0)
  K_loc = rep(K,n)
 
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
  mu = rep(mu0,n)
  
  # M0 : initial migration rate
  # string of events
  #B<- c("spec A","spec B","allo spec","ext A","ext B", "ext AB", "mig AB", "mig BA", "con A", "con B", "sym spec A", "sym spec B")
  
  #number of events
  #sym spec
  n = 4
  num_ss = 0
  for(i in 1:(n-1))
  {
    num_ss = num_ss+choose(n-1,i)
  }
  num_ss = n*(1+num_ss)
  
  #ext 
  num_ext = num_ss
  
  #allo spec
  num_as = choose(n,2)
  
  #mig 
  num_mig = 0
  for(i in 1:(n-1)){
    num_mig = num_mig + choose(n,i)*(n-i)
  }
  
  #number of events
  num_event = num_ss + num_ext + num_as + num_mig
  
  B<- c(1:num_event)
  
  Na = parsN[1]
  Nb = parsN[2]
  Nc = parsN[3]
  Nd = parsN[4]
  Nab=Nac=Nad=Nbc=Nbd= Ncd=Nabc=Nabd=Nacd = Nbcd=Nabcd=0
  
  
  N=sum(parsN)
  Ntable_index = Nindex(n)
  Nlength = sum(Ntable_index)
  Ntable = matrix(0,nrow =1, ncol =Nlength)
  Ntable[1:n] = parsN
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
    loc1= matrix(0,1,4)
    loc1[1,loc[1]] = 1 
    loctable = rbind(loctable, loc1)
    linlist = cbind(L[,3], L[,5])
    Ntable[1,loc[1]] = Ntable[1,loc[1]] -1 
    newL=j
  }
  # print(loctable)
  Ntable = cbind(Na,Nb,Nc,Nd,Nab,Nac,Nad,Nbc,Nbd,Ncd,Nabc,Nabd,Nacd,Nbcd,Nabcd)
  
  spec_num = sum(Nindex(n)[1,])
   N_loc = matrix(0,nrow = n,ncol = spec_num)
   N_loc_col = matrix(0,nrow = n,ncol = spec_num)
     Ndistribution = event_matrix(n)
  
     #index for sym speciation
     for(j in 1:n){
     N_loc_col[j,] = which(Ndistribution[j,] == 1,arr.ind = TRUE)#index of each loc in sym spec table
     }
     B_symspec = c(N_loc_col)
     #index for extinction
     B_ext = NULL
     x = Ndistribution
     y = split(x, rep(1:ncol(x), each = nrow(x)))
     for(j in 1:n){
       x1 = x
       x1[j,] = x1[j,]-1
       y1 = split(x1, rep(1:ncol(x1), each = nrow(x1)))
       z = match(y1,y)
       z = z[!is.na(z)]
       B_ext = rbind(B_ext, z)
     }
     B_ext = c(B_ext)
     
     #index for allo speciation
     lambda_allo = NULL
     N_allo = colSums(Ndistribution)
     allo_index = which(N_allo == 2)
     allo_col = which(Ndistribution[,allo_index] == 1,arr.ind = TRUE)
     B_allospec = matrix(allo_col[,1],2)
     B_allodau1 = c(B_allospec[1,])
     B_allodau2 = c(B_allospec[2,])
     
     #index for migration
     B_mig_target = which(Ndistribution == 0, arr.ind = TRUE)
     B_mig_from = c(B_mig_target[,2])
     B_mig_to = c(B_mig_target[,1])
    
     B_mig_bec = NULL
     for(j in 1:length(B_mig_from)){
       Ndis_aftermig = Ndistribution
       Ndis_aftermig[B_mig_to[j],B_mig_from[j]] = 1
       loc_aftermig = as.vector(Ndis_aftermig[,B_mig_from[j]])
       B_mig_bec = c(B_mig_bec,which(sapply( y ,function(x)all(x==loc_aftermig))))
     }
     B_mig_bec = as.numeric(B_mig_bec)
  i = 0
  while(t[i+1]< age){
    i<-i+1
   # print(t[i])
    Na = Ntable[,1]
    Nb = Ntable[,2]
    Nc = Ntable[,3]
    Nd = Ntable[,4]
    Nab = Ntable[,5]
    Nac = Ntable[,6]
    Nad = Ntable[,7]
    Nbc = Ntable[,8]
    Nbd = Ntable[,9]
    Ncd = Ntable[,10]
    Nabc = Ntable[,11]
    Nabd = Ntable[,12]
    Nacd = Ntable[,13]
    Nbcd = Ntable[,14]
    Nabcd = Ntable[,15]
    
    
    num_1 = 0
    for(m in 1:(n-1))
    {
      num_1 = num_1 +choose(n-1,m)
    }
    num_1 = num_1 +1
    
    probs = rep(0,num_event)
    
    # speciation event & extinction event
    lambda_sym = rep(0,n)
    mu = rep(mu0,n)
    
    sym_spec_event = matrix(0,nrow = n,ncol = spec_num)
    ext_event = matrix(0,nrow = n,ncol = spec_num)
    
    for(j in 1:n){
      N_loc[j,] = Ntable[i,which(Ndistribution[j,]==1)]  #number of each loc in sym spec table
      lambda_sym[j]=max(lambda0*(1-sum(N_loc[j,])/K_loc[j]),0)
      sym_spec_event[j,] = lambda_sym[j]*N_loc[j,] 
      ext_event[j,] = mu[j] *N_loc[j,]
    }
    prob_spec_sym = c(sym_spec_event)
    prob_ext = c(ext_event)
    
  
    # Migration 
    prob_mig = NULL
    Mig_dir = rep(0,n)
    for(j in 1:n){
      Mig_dir[j] = max(M0*((1-sum(N_loc[j,])/K_loc[j])),0)
    }
    for(j in 1:(ncol(Ntable)-1)){
      tar = which(Ndistribution[,j]==0)
      prob_mig_each = Ntable[i,j]*Mig_dir[tar]
      prob_mig_each = matrix(prob_mig_each,ncol = length(prob_mig_each))
      prob_mig = cbind(prob_mig,prob_mig_each)
      }
    
    
    # Allopatric speciation event
    lambda_allo = NULL
    N_allo = colSums(Ndistribution)
    allo_index = which(N_allo == 2)
    for(j in allo_index){
      lambda_allo_each = max(lambda_allo0/sum(Mig_dir[which(Ndistribution[,j] == 1)]), 0 )
      lambda_allo = cbind(lambda_allo_each,lambda_allo)
    }
    prob_spec_allo = lambda_allo
   
    
    
    
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
    
    
    probs= c(prob_spec_sym,prob_ext,prob_spec_allo,prob_mig)
    probs_part1 = length(prob_spec_sym)
    probs_part2 = length(prob_ext)
    probs_part3 = length(prob_spec_allo)
    probs_part4 = length(prob_mig)
    #Total rate
    TR=sum(probs)
    if(TR==0) break
    else{
      # print(linlist)
      # print(Ntable)
      # print(which(probs < 0))
      A<-DDD::sample2(B,1, prob = probs)
      t[i+1]=t[i]+rexp(1,rate=TR)
      if(t[i+1]>age) break
      loc2 = matrix(0,1,n)
     
      # Sympatric speciation 
      if (is.element(A,1:prob_part1)){
        b1<-A%%n
        if(b1 == 0) b1 = n
        Ntable=rbind(Ntable,Ntable[i,])
        Ntable[i+1,b1] = Ntable[i,b1]+1
        newL = newL + 1;
        list0 = matrix(linlist,ncol = 2)
        b3 <- B_symspec[A]
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
        loc2[1,b1] = 1
        loctable = rbind(loctable, loc2)
      }
      
      #Extinction
      else if(is.element(A,(probs_part1+1):(probs_part1+probs_part2))) {
        b1<-A%%n
        if(b1 == 0) b1 = 3
        Ntable=rbind(Ntable,Ntable[i,])
        b3 <- B1[A-length(prob_spec_sym)]
        list0 = matrix(linlist,ncol = 2)
        list1 = linlist[list0[,2]== b3]
        list2 = matrix(list1, ncol = 2)
        linlist1 = list2[,1]
        ranL= DDD::sample2(linlist1,1)
        loctable[abs(ranL),b1] = 0
       if(is.element(A,(length(prob_spec_sym)+1):(length(prob_spec_sym)+n))){
         Ntable[i+1,b3] = max(Ntable[i,b3]-1,0)
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
          b2 <- B_ext[A-length(prob_spec_sym)-n]
          Ntable[i+1,b3] = Ntable[i,b3]-1
          Ntable[i+1,b2] = Ntable[i,b2]+1
          L[abs(ranL),5] <- b2
          v = which(linlist[,1] == ranL)
          linlist[v,2] = b2
          linlist = linlist[order(linlist[,1]),]
          linlist = matrix(linlist,ncol=2)
          # loctable[abs(ranL),A-6] = 0
        }
      }
      
      
      # Allopatric speciation in AB
      else if(is.element(A,((probs_part1+probs_part2)+1):(probs_part1+probs_part2+probs_part3))) {
        A1 = A - (probs_part1+probs_part2)
        Ntable=rbind(Ntable,Ntable[i,])
        Ntable[i+1,B_allodau1[A1]] = Ntable[i,B_allodau1[A1]]+1
        Ntable[i+1,B_allodau2[A1]] = Ntable[i,B_allodau2[A1]]+1
        Ntable[i+1,allo_index[A1]] = Ntable[i,allo_index[A1]]-1
        newL = newL + 1
        list0 = matrix(linlist,ncol = 2)
        list1 = linlist[list0[,2]==allo_index[A1]]
        list2 = matrix(list1, ncol = 2)
        linlist1 = list2[,1]
        ranL= DDD::sample2(linlist1,1)
        # A2 = c(1,1,2)
        # A3 = c(2,3,3)
        L[abs(ranL),5] <- B_allodau1[A1]
        loctable[abs(ranL), B_allodau2[A1]] = 0
        v = which(linlist[,1] == ranL)
        linlist[v,2] = B_allodau1[A1]
        L = rbind(L,c(t[i+1],ranL,sign(ranL) * newL,-1,B_allodau2[A1]))
        linlist = rbind(linlist,c(sign(ranL) * newL,B_allodau2[A1]))  
         loc2[1,B_allodau1[A1]] = 1 
         loctable = rbind(loctable,loc2)
      }
     
      #Migration
      else {
        Mig1 = B_mig_from # c(1,1,2,2,3,3,4,5,6)
        Mig2 =  B_mig_bec  #c(4,5,4,6,5,6,7,7,7)
        Mig3 = B_mig_to # c(2,3,1,3,1,2,3,2,1)
        Am = A - (probs_part1+probs_part2+probs_part3)
        b2 = Mig1[Am]
        b3 = Mig2[Am]
        Ntable=rbind(Ntable,Ntable[i,])
        Ntable[i+1,b2]=Ntable[i,b2]-1
        Ntable[i+1,b3]=Ntable[i,b3]+1
        linlist = matrix(linlist,ncol = 2)
        list1 = linlist[linlist[,2]==b2]
        list2 = matrix(list1, ncol = 2)
        linlist1 = list2[,1]
        ranL1 = DDD::sample2(linlist1,1)
        L[abs(ranL1),5] <- b3
        v = which(linlist[,1] == ranL1)
        linlist[v,2] = b3
        linlist = linlist[order(linlist[,1]),]
        linlist = matrix(linlist,ncol=2)
         loctable[abs(ranL1),Mig3[Am]] = 1
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
  
  # Ltable = cbind(L, loctable)
  
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
  out = list(tes = tes,tas = tas,L = L, loctable = loctable)
#    file = paste("/Users/mac/Dropbox/R/MigrationModelsim/ModelTest/",parsN[1],parsN[2],parsN[3],"&",pars[1],pars[2],pars[3],"&",age,"sim.txt")
#    save(out,file = file)
   return(out)
  }
  }
  else print(paste("please input 1 or 2 species in total!"))
}
print(paste("Please input Na, Nb, Nab, age with no more than 2 species in total!"))