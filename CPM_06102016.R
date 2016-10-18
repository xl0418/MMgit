library(iterators)
# library(doParallel)
# library(foreach)
library(Matrix)

lambda_na <- NA
lambda_nb = NA
Mab = NA
Mba = NA
lambda_ab = NA
VDelta = NA
E1 = NA
E2 = NA
E3 = NA
E4 = NA
E5 = NA
E6 = NA
E7 = NA
E8 = NA
E9 = NA
E10 = NA
E11 = NA
E12 = NA
E13 = NA
E14 = NA
E15 = NA
E16 = NA
mu_na = NA
mu_nb = NA
C_a = NA
C_b = NA
K_a = NA
K_b = NA
fun = NA

CPM <- function (k,N,pars,M_0=0.4) {
  timestart<-Sys.time()
  lambda0 = pars[1]
  mu0 = pars[2]
  mu_na <- mu0
  mu_nb <- mu0
  C_a <- mu0
  C_b <- mu0
  K = pars[3]
  K_a <- lambda0*K/(lambda0 - mu0)
  K_b = K_a
  
  som <- expand.grid(0:k,0:k,0:k)         #compute all the combinations of species in phylogeny
  sol <- expand.grid(0:N,0:N,0:N)          #compute all the combinations of missing species into three groups
  sol <- cbind(sol,sum=apply(sol,1,sum))     #compute the sum of each row 
  sol <- sol[sol$sum<=N,1:3]                  #extract the rows with sum less than or equal to N
  sol <- sol[order(sol[,1],sol[,2],sol[,3]),]       #order the combinations
  som <- cbind(som,sum=apply(som,1,sum))     #compute the sum of each row 
  som <- som[som$sum ==k,1:3]                  #extract the rows with the sum equal to k
  som <- som[order(som[,1],som[,2],som[,3]),]    #order the combinations
  som1 <- as.matrix(som)
  sol1 <- as.matrix(sol)
  mat <- Matrix(0,nrow=nrow(som1)*nrow(sol1),ncol=ncol(som1)+ncol(sol1),sparse=TRUE)
  
  
  som2 = paste(som1[,1],som1[,2],som1[,3])
  sol2 = paste(sol1[,1],sol1[,2],sol1[,3])
  
  mat1 = expand.grid(sol2,som2)
  mat2 = paste(mat1[,1],mat1[,2])
  
  out <- Matrix(0,nrow=nrow(mat),ncol=nrow(mat),sparse=TRUE)             #build the matrix
  dimnames(out) <- list(mat2,mat2)                     #name the rows and columes by the combinations
  
  timemid1<-Sys.time()
  # print(paste("name time", timemid1-timestart))
  
  #define the matrix	
  A <- dimnames(out)[[1]]
  B <- dimnames(out)[[2]]
  
  N0 = 0
  for(i in 1:(N+1)){
    N0 = i*(i+1)+N0
  }
  N1 = N0/2
  
  N2<-choose(k+2,2)
  
  dimM = N1 * N2
  
  rka = rep(som1[,1],each = N1)
  rkb = rep(som1[,2],each = N1)
  rkab = rep(som1[,3],each = N1)
  rna = rep(sol1[,1],N2)
  rnb = rep(sol1[,2],N2)
  rnab = rep(sol1[,3],N2)
  
  E1 <- cbind(rna,rnb-1,rnab,rka-1,rkb,rkab+1) #allo spe
  E2 <- cbind(rna,rnb,rnab,rka-1,rkb,rkab+1)   #contrac in B
  E3 <- cbind(rna-1,rnb,rnab,rka,rkb-1,rkab+1) #allo spe
  E4 <- cbind(rna,rnb,rnab,rka,rkb-1,rkab+1)   #contrac in A
  E5 <- cbind(rna-1,rnb-1,rnab+1,rka,rkb,rkab) #allo spe
  E6 <- cbind(rna-1,rnb,rnab,rka,rkb,rkab)     #spe in A
  E7 <- cbind(rna-1,rnb,rnab+1,rka,rkb,rkab)   #contrac in B
  E8 <- cbind(rna,rnb-1,rnab,rka,rkb,rkab)     #spe in B
  E9 <- cbind(rna,rnb-1,rnab+1,rka,rkb,rkab)   #contrac in A
  E10 <- cbind(rna,rnb,rnab,rka,rkb,rkab)      #nothing
  E11 <- cbind(rna,rnb+1,rnab-1,rka,rkb,rkab)  #mig B to A
  E12 <- cbind(rna,rnb+1,rnab,rka,rkb,rkab)    #ext in B
  E13 <- cbind(rna+1,rnb,rnab-1,rka,rkb,rkab)  #mig A to B
  E14 <- cbind(rna+1,rnb,rnab,rka,rkb,rkab)    #ext in A
  E15 <- cbind(rna,rnb,rnab,rka,rkb+1,rkab-1)  #mig B to A
  E16 <- cbind(rna,rnb,rnab,rka+1,rkb,rkab-1)  #mig A to B
  
  N_a<- cbind(E1[,1]+E1[,3]+E1[,4]+E1[,6],E2[,1]+E2[,3]+E2[,4]+E2[,6],E3[,1]+E3[,3]+E3[,4]+E3[,6],E4[,1]+E4[,3]+E4[,4]+E4[,6],E5[,1]+E5[,3]+E5[,4]+E5[,6],E6[,1]+E6[,3]+E6[,4]+E6[,6],E7[,1]+E7[,3]+E7[,4]+E7[,6],E8[,1]+E8[,3]+E8[,4]+E8[,6],E9[,1]+E9[,3]+E9[,4]+E9[,6],E10[,1]+E10[,3]+E10[,4]+E10[,6],E11[,1]+E11[,3]+E11[,4]+E11[,6],E12[,1]+E12[,3]+E12[,4]+E12[,6],E13[,1]+E13[,3]+E13[,4]+E13[,6],E14[,1]+E14[,3]+E14[,4]+E14[,6],E15[,1]+E15[,3]+E15[,4]+E15[,6],E16[,1]+E16[,3]+E16[,4]+E16[,6])
  N_b<- cbind(E1[,2]+E1[,3]+E1[,5]+E1[,6],E2[,2]+E2[,3]+E2[,5]+E2[,6],E3[,2]+E3[,3]+E3[,5]+E3[,6],E4[,2]+E4[,3]+E4[,5]+E4[,6],E5[,2]+E5[,3]+E5[,5]+E5[,6],E6[,2]+E6[,3]+E6[,5]+E6[,6],E7[,2]+E7[,3]+E7[,5]+E7[,6],E8[,2]+E8[,3]+E8[,5]+E8[,6],E9[,2]+E9[,3]+E9[,5]+E9[,6],E10[,2]+E10[,3]+E10[,5]+E10[,6],E11[,2]+E11[,3]+E11[,5]+E11[,6],E12[,2]+E12[,3]+E12[,5]+E12[,6],E13[,2]+E13[,3]+E13[,5]+E13[,6],E14[,2]+E14[,3]+E14[,5]+E14[,6],E15[,2]+E15[,3]+E15[,5]+E15[,6],E16[,2]+E16[,3]+E16[,5]+E16[,6])
  
  
  lambda_na <- lambda0*(1-N_a/K_a)
  lambda_na[which(lambda_na <0) ] = 0
  lambda_nb <- lambda0*(1-N_b/K_b)
  lambda_nb[which(lambda_nb <0) ] = 0
  Mab <- M_0*(1-N_b/K_b)
  Mba <- M_0*(1-N_a/K_a)
  Mab[which(Mab <0) ] = 0
  Mba[which(Mba <0) ] = 0
  
  lambda_ab <- matrix(0.2,nrow = nrow(Mab), ncol = ncol(Mab))   # 2/(Mab+Mba)
  
  
  
  VDelta <- matrix(0,nrow= dimM,ncol= 16) 
  
  
  VDelta[,10] <- -(mu_na+lambda_na[,10] +Mab[,10] )*(E10[,1] +E10[,4])-(mu_nb+lambda_nb[,10]+Mba[,10])*(E10[,2] +E10[,5])-(C_a+C_b)*(E10[,3] +E10[,6])-((lambda_na[,10]+lambda_nb[,10])*(E10[,3] +E10[,6] )+lambda_ab[,10]*(E10[,3] +2*E10[,6] ))
  VDelta[,14] <-mu_na*E14[,1] 
  VDelta[,6] <-lambda_na[,6]*(E6[,1] +2*E6[,4] +E6[,3] +2*E6[,6] ) 
  VDelta[,13] <-Mab[,13]*E13[,1] 
   VDelta[,12] <-mu_nb*E12[,2]
  VDelta[,8] <-lambda_nb[,8]*(E8[,2] +2*E8[,5] +E8[,3] +2*E8[,6]) 
   VDelta[,11] <-Mba[,11]*E11[,2] 
   VDelta[,9] <-C_a*E9[,3] 
   VDelta[,7] <-C_b*E7[,3] 
   VDelta[,5] <-lambda_ab[,5]*E5[,3]
   VDelta[,16] <- Mab[,16]*E16[,4] 
   VDelta[,15] <- Mba[,15]*E15[,5] 
   VDelta[,2] <- C_b*E2[,6] 
   VDelta[,4] <- C_a*E4[,6] 
   VDelta[,1] <- lambda_ab[,1]*E1[,6]
  VDelta[,3] <- lambda_ab[,3]*E3[,6]
  
  
  E1name = paste(E1[,1],E1[,2],E1[,3],E1[,4],E1[,5],E1[,6],sep = " ")
  E2name = paste(E2[,1],E2[,2],E2[,3],E2[,4],E2[,5],E2[,6], sep = " ")
  E3name = paste(E3[,1],E3[,2],E3[,3],E3[,4],E3[,5],E3[,6], sep = " ")
  E4name = paste(E4[,1],E4[,2],E4[,3],E4[,4],E4[,5],E4[,6], sep = " ")
  E5name = paste(E5[,1],E5[,2],E5[,3],E5[,4],E5[,5],E5[,6], sep = " ")
  E6name = paste(E6[,1],E6[,2],E6[,3],E6[,4],E6[,5],E6[,6], sep = " ")
  E7name = paste(E7[,1],E7[,2],E7[,3],E7[,4],E7[,5],E7[,6], sep = " ")
  E8name = paste(E8[,1],E8[,2],E8[,3],E8[,4],E8[,5],E8[,6], sep = " ")
  E9name = paste(E9[,1],E9[,2],E9[,3],E9[,4],E9[,5],E9[,6], sep = " ")
  E10name = paste(E10[,1],E10[,2],E10[,3],E10[,4],E10[,5],E10[,6], sep = " ")
  E11name = paste(E11[,1],E11[,2],E11[,3],E11[,4],E11[,5],E11[,6], sep = " ")
  E12name = paste(E12[,1],E12[,2],E12[,3],E12[,4],E12[,5],E12[,6], sep = " ")
  E13name = paste(E13[,1],E13[,2],E13[,3],E13[,4],E13[,5],E13[,6], sep = " ")
  E14name = paste(E14[,1],E14[,2],E14[,3],E14[,4],E14[,5],E14[,6], sep = " ")
  E15name = paste(E15[,1],E15[,2],E15[,3],E15[,4],E15[,5],E15[,6], sep = " ")
  E16name = paste(E16[,1],E16[,2],E16[,3],E16[,4],E16[,5],E16[,6], sep = " ")
  
  colname = cbind(E1name, E2name, E3name, E4name, E5name, E6name,E7name, E8name, E9name, E10name, E11name, E12name, E13name, E14name, E15name, E16name )
  colname = t(colname)
  Pcol = which(colname %in% B)
  Pcol1 = ceiling(Pcol/16)
  Pout = cbind(mat2[Pcol1],colname[Pcol])
  P1 = match(Pout[,1],mat2)
  P2 = match(Pout[,2],mat2)
  P3 = cbind(P1,P2)
  VDelta = t(VDelta)
  out[P3] = VDelta[Pcol]
  timemid2<-Sys.time()
  # print(paste("value time", timemid2-timemid1))
  
  # out1<-as.matrix(out)
  return(out)
  # file = paste("/Users/mac/Dropbox/R/likelihoodalgorithm/",k,N,"PM_new.Rdata",sep = "")
  # save(out1, file= file)
  
#     out1<-as.matrix(out)
#     file = paste("/Users/mac/Dropbox/R/likelihoodalgorithm/",k,N,"PMnew2.xls",sep = "")
#     write.csv(out1,file=file,fileEncoding = "",row.names=FALSE,col.names = FALSE)
}


