library(Matrix)
library(iterators)
library(doParallel)
library(foreach)

lambda_na <- NA
lambda_nb = NA
Mab = NA
Mba = NA
lambda_ab = NA
VDelta = NA
E1 = NA
E2 = NA
E3 = NA
K_a = NA
K_b = NA
fun = NA

# 
# V2Delta <- function(i,ka,kb,kab,na,nb,nab,lambda_na,lambda_nb,lambda_ab,S1,S2,S3,DDelta1,DDelta2,DDelta3){
#   if(i ==2)  VDelta<- S1*lambda_na * ka + S3*DDelta1*lambda_na * kab
#   else if(i ==3) VDelta<- S2*lambda_nb * kb + S3*DDelta2*lambda_nb * kab
#   else if(i ==1) VDelta<- S3*DDelta3*lambda_ab * kab
#   else print("warning: i is out of range")
# }

CTM <- function (k,N,pars,M_0=0.4, lambda_ab = 0.2) {
  timemid1<-Sys.time()
  lambda0 = pars[1]
  mu0 = pars[2]
  mu_na = mu0
  mu_nb = mu0
  C_a = mu0
  C_b = mu0
  K = pars[3]
  K_a <- lambda0*K/(lambda0 - mu0)
  K_b = K_a  
  som <- expand.grid(0:k,0:k,0:k)         #compute all the combinations of species in phylogeny
  k1 <- k+1
  son <- expand.grid(0:k1,0:k1,0:k1)
  son <- cbind(son,sum=apply(son,1,sum))     #compute the sum of each row 
  son <- son[son$sum ==k1,1:3]                  #extract the rows with sum less than or equal to N
  son <- son[order(son[,1],son[,2],son[,3]),] 
  sol <- expand.grid(0:N,0:N,0:N)          #compute all the combinations of missing species into three groups
  sol <- cbind(sol,sum=apply(sol,1,sum))     #compute the sum of each row 
  sol <- sol[sol$sum<=N,1:3]                  #extract the rows with sum less than or equal to N
  sol <- sol[order(sol[,1],sol[,2],sol[,3]),]       #order the combinations
  som <- cbind(som,sum=apply(som,1,sum))     #compute the sum of each row 
  som <- som[som$sum ==k,1:3]                  #extract the rows with the sum equal to k
  som <- som[order(som[,1],som[,2],som[,3]),]    #order the combinations
  som1 <- as.matrix(som)
  sol1 <- as.matrix(sol)
  son1 <- as.matrix(son)
  mat <- Matrix(0,nrow=nrow(som1)*nrow(sol1),ncol=ncol(som1)+ncol(sol1),sparse=TRUE)
  mat1 <- Matrix(0,nrow=nrow(son1)*nrow(sol1),ncol=ncol(son1)+ncol(sol1),sparse=TRUE)
  
  som2 = paste(som1[,1],som1[,2],som1[,3])
  son2 = paste(son1[,1],son1[,2],son1[,3])
  sol2 = paste(sol1[,1],sol1[,2],sol1[,3])
  
  matm1 = expand.grid(sol2,som2)
  matm2 = paste(matm1[,1],matm1[,2])    #column name
  matn1 = expand.grid(sol2,son2)
  matn2 = paste(matn1[,1],matn1[,2])    #row name
  
  out <- Matrix(0,nrow=nrow(mat1),ncol=nrow(mat),sparse=TRUE)             #build the matrix
  dimnames(out) <- list(matn2,matm2)                     #name the rows and columes by the combinations
  
  #define the matrix	
  A <- dimnames(out)[[1]]
  B <- dimnames(out)[[2]]
  
  
  
  N0 = 0
  for(i in 1:(N+1)){
    N0 = i*(i+1)+N0
  }
  N1 = N0/2
  
  N2<-choose(k+2,2)
  N3<-choose(k+3,2)
  
  cdimM = N1 * N2
  rdimM = N1 * N3
  
  
  rka = rep(son1[,1],each = N1)
  rkb = rep(son1[,2],each = N1)
  rkab = rep(son1[,3],each = N1)
  rna = rep(sol1[,1],N3)
  rnb = rep(sol1[,2],N3)
  rnab = rep(sol1[,3],N3)
  
  cka = rep(som1[,1],each = N1)
  ckb = rep(som1[,2],each = N1)
  ckab = rep(som1[,3],each = N1)
  cna = rep(sol1[,1],N2)
  cnb = rep(sol1[,2],N2)
  cnab = rep(sol1[,3],N2)
  
  E1 <- cbind(rna,rnb,rnab,rka-1,rkb,rkab) #sym spe in A
  E2 <- cbind(rna,rnb,rnab,rka,rkb-1,rkab)   #sym spe in B
  E3 <- cbind(rna,rnb,rnab,rka-1,rkb-1,rkab+1) #allo spe
  
  
  N_a<- cbind(E1[,1]+E1[,3]+E1[,4]+E1[,6],E2[,1]+E2[,3]+E2[,4]+E2[,6],E3[,1]+E3[,3]+E3[,4]+E3[,6])
  N_b<- cbind(E1[,2]+E1[,3]+E1[,5]+E1[,6],E2[,2]+E2[,3]+E2[,5]+E2[,6],E3[,2]+E3[,3]+E3[,5]+E3[,6])
  
  lambda_na <- lambda0*(1-N_a/K_a)
  lambda_na[which(lambda_na <0) ] = 0
  lambda_nb <- lambda0*(1-N_b/K_b)
  lambda_nb[which(lambda_nb <0) ] = 0
  
#   Mab <<- M_0*(1-N_b/K_b)
#   Mba <<- M_0*(1-N_a/K_a)
  
  VDelta <- matrix(0,nrow= rdimM,ncol= 3) 
  
  S1 <- E1[,4] /(E1[,4]+E1[,5]+E1[,6])
  S2 <- E2[,4] /(E2[,4]+E2[,5]+E2[,6])
  S31 <- E1[,6] /(E1[,4]+E1[,5]+E1[,6])
  S32 <- E2[,6] /(E2[,4]+E2[,5]+E2[,6])
  S33 <- E3[,6] /(E3[,4]+E3[,5]+E3[,6])
  
  lambda <- lambda_na+lambda_nb+lambda_ab
  DDelta1<- lambda_na[,1] /lambda[,1]
  DDelta2<- lambda_nb[,2] /lambda[,2]
  DDelta3<- lambda_ab /lambda[,3]
  
  
     VDelta[,1] <- S1*lambda_na[,1] * E1[,4] + S31 *DDelta1*lambda_na[,1] * E1[,6]
    VDelta[,2] <-   S2*lambda_nb[,2] * E2[,5] + S32*DDelta2*lambda_nb[,2] * E2[,6]
    VDelta[,3] <-   S33*DDelta3*lambda_ab * E3[,6]
    
  
  E1name = paste(E1[,1],E1[,2],E1[,3],E1[,4],E1[,5],E1[,6],sep = " ")
  E2name = paste(E2[,1],E2[,2],E2[,3],E2[,4],E2[,5],E2[,6], sep = " ")
  E3name = paste(E3[,1],E3[,2],E3[,3],E3[,4],E3[,5],E3[,6], sep = " ")
  
  
  colname = cbind(E1name, E2name, E3name)
  colname = t(colname)
  Pcol = which(colname %in% B)
  Pcol1 = ceiling(Pcol/3)
  Pout = cbind(matn2[Pcol1],colname[Pcol])
  P1 = match(Pout[,1],matn2)
  P2 = match(Pout[,2],matm2)
  P3 = cbind(P1,P2)
  VDelta = t(VDelta)
  out[P3] = VDelta[Pcol]
  timemid2<-Sys.time()
  # print(paste("value time", timemid2-timemid1))
  out
#     out1<-as.matrix(out)
#     file = paste("/Users/mac/Dropbox/R/likelihoodalgorithm/",k,N,"TMpara.xls",sep = "")
#     write.csv(out1,file=file,fileEncoding = "",row.names=FALSE,col.names = FALSE)
}
