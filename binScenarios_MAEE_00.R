#########################################
# Simulation program for binary outcomes
# Balanced designs
# GEE with MAEE
#########################################
# n: number of clusters
# m: number of divisions per cluster
# k: number of participants per division
# l: number of evaluations per participant
# beta - marginal means
#        beta_0 - the marginal mean P0 in the control arm
#        beta_1 - the marginal mean P1 in the intervention arm
# alpha - correlations 
#         alpha_0 - correlation between different evaluations from the same participant
#         alpha_1 - correlation between evaluations from different participants but within the same division
#         alpha_2 - correlation between evaluations from different divisions but within the same cluster
#########################################

binSIMULATE<-function(n,m,k,l,beta,alpha){
  
  # Source program code
  source("../../../binGEN.R")
  source("../../../binMAEE.R")
  
  # Store results
  nrep<-1000
  nerror<-0
  results<-matrix(NA,nrep,16)
  colnames(results)<-c("nSIM","pzeroflag","pminflag","poneflag","pmaxflag",
                       "beta","MB","ROB","KC","MD","FG","MBN","AVG","K","p","error")
  results[,14]<-n
  results[,15]<-2
  results[,16]<-0
  
  results0<-matrix(NA,nrep,8)
  results1<-matrix(NA,nrep,8)
  results2<-matrix(NA,nrep,8)
  colnames(results0)<-c("alpha0","ROB","KC","MD","FG","MBN","AVG","error")
  colnames(results1)<-c("alpha1","ROB","KC","MD","FG","MBN","AVG","error")
  colnames(results2)<-c("alpha2","ROB","KC","MD","FG","MBN","AVG","error")
  
  # Create X matrix
  X<-cbind(rep(1,n*m*k*l),c(rep(0,n*m*k*l/2),rep(1,n*m*k*l/2)))
  
  # Create Z matrix
  CREATEZ<-function(n,m,k,l){
    rho<-1
    nu<-2
    lambda<-3
    em1<-(1-rho)*diag(m*k*l)
    em2<-(rho-nu)*kronecker(diag(m*k),matrix(1,l,l))
    em3<-(nu-lambda)*kronecker(diag(m),matrix(1,k*l,k*l))
    em4<-lambda*matrix(1,m*k*l,m*k*l)
    POS<-em1+em2+em3+em4
    zrow<-diag(3)
    z_c<-NULL
    for(jp in 1:(m*k*l-1)){
      for(t in (jp+1):(m*k*l)){
        z_c<-rbind(z_c,zrow[POS[jp,t],])
      }
    }
    Z<-z_c
    for(i in 1:n){
      Z<-rbind(Z,z_c)
    }
    rm(z_c)
    return(Z)
  }
  Z<-CREATEZ(n,m,k,l) # large matrix
  
  # Create id
  id<-rep(1:n,each=m*k*l)
  clsize<-rep(m*k*l,n)
  
  # For Loop
  j<-0
  for(i in 1:nrep){
    # Generate outcome and check; store number of simulated data sets
    passGEN<-1
    nSIM<-0
    while(passGEN==1){
      nSIM<-nSIM+1
      j<-j+1
      set.seed(j)
      y<-try(binGEN(n,m,k,l,beta,alpha),silent=TRUE)
      if(class(y)[1]=="try-error"){passGEN<-1;next}
      passGEN<-0
    }
    y<-c(y)
    results[i,1]<-nSIM
    
    # MAEE program
    out=try(binMAEE(y=y,X=X,id=id,n=clsize,Z=Z,maxiter=25,epsilon=0.001,printrange="NO",shrink="ALPHA",makevone="NO"),silent=TRUE)
    if(class(out)[1]=="try-error" | is.null(out)==1){
      nerror<-nerror+1
      results[i,16]<-1
      next
    }
    
    # Save results
    results[i,(6:12)]<-out$outbeta[2,]
    results0[i,(1:6)]<-out$outalpha[1,]
    results1[i,(1:6)]<-out$outalpha[2,]
    results2[i,(1:6)]<-out$outalpha[3,]
    
    # Loop control
    if(i %% 20 == 0){print(paste("Iteration",i,"NERROR =",nerror))}
  }
  
  results[,13]<-rowMeans(results[,9:10], na.rm=TRUE)
  results0[,7]<-rowMeans(results0[,3:4], na.rm=TRUE)
  results1[,7]<-rowMeans(results1[,3:4], na.rm=TRUE)
  results2[,7]<-rowMeans(results2[,3:4], na.rm=TRUE)
  
  results0[,8]<-results[,16]
  results1[,8]<-results[,16]
  results2[,8]<-results[,16]
  
  # save data sets and clear workspace
  if(beta[2]-beta[1]==0){
    name<-"Size"
  }else{
    name<-"Power"}
  save(results,file=paste0(name,"_beta.RData"))
  save(results0,file=paste0(name,"_alpha0.RData"))
  save(results1,file=paste0(name,"_alpha1.RData"))
  save(results2,file=paste0(name,"_alpha2.RData"))
}

# ---------------------------------------------------------------------------------

# Power

# (1)
setwd('binResults/binRData/1row')
binSIMULATE(14, 2, 3, 5, c(0.2, 0.5), c(0.4, 0.1, 0.03)) #14

# (2)
setwd('../2row')
binSIMULATE(14, 2, 3, 10, c(0.2, 0.5), c(0.4, 0.1, 0.03)) #14

# (3)
setwd('../3row')
binSIMULATE(14, 2, 4, 5, c(0.2, 0.5), c(0.4, 0.1, 0.03)) #14

# (4)
setwd('../4row')
binSIMULATE(12, 3, 3, 5, c(0.2, 0.5), c(0.4, 0.1, 0.03)) #12

# (5)
setwd('../5row')
binSIMULATE(10, 2, 3, 5, c(0.2, 0.5), c(0.15, 0.08, 0.02)) #10

# (6)
setwd('../6row')
binSIMULATE(10, 2, 3, 10, c(0.2, 0.5), c(0.15, 0.08, 0.02)) #10

# (7)
setwd('../7row')
binSIMULATE(10, 2, 4, 5, c(0.2, 0.5), c(0.15, 0.08, 0.02)) #10

# (8)
setwd('../8row')
binSIMULATE(8, 3, 3, 5, c(0.2, 0.5), c(0.15, 0.08, 0.02)) #8

# (9)
setwd('../9row')
binSIMULATE(8, 2, 3, 5, c(0.2, 0.5), c(0.1, 0.02, 0.01)) #8

# (10)
setwd('../10row')
binSIMULATE(8, 3, 3, 5, c(0.2, 0.5), c(0.1, 0.02, 0.01)) #8

# (11)
setwd('../11row')
binSIMULATE(8, 3, 3, 5, c(0.2, 0.5), c(0.05, 0.05, 0.02)) #8

# (12)
setwd('../12row')
binSIMULATE(22, 2, 3, 5, c(0.1, 0.3), c(0.4, 0.1, 0.03)) #22

# (13)
setwd('../13row')
binSIMULATE(20, 2, 3, 10, c(0.1, 0.3), c(0.4, 0.1, 0.03)) #20

# (14)
setwd('../14row')
binSIMULATE(20, 2, 4, 5, c(0.1, 0.3), c(0.4, 0.1, 0.03)) #20

# (15)
setwd('../15row')
binSIMULATE(16, 3, 3, 5, c(0.1, 0.3), c(0.4, 0.1, 0.03)) #16

# (16)
setwd('../16row')
binSIMULATE(16, 2, 3, 5, c(0.1, 0.3), c(0.15, 0.08, 0.02)) #16

# (17)
setwd('../17row')
binSIMULATE(14, 2, 3, 10, c(0.1, 0.3), c(0.15, 0.08, 0.02)) #14

# (18)
setwd('../18row')
binSIMULATE(14, 2, 4, 5, c(0.1, 0.3), c(0.15, 0.08, 0.02)) #14

# (19)
setwd('../19row')
binSIMULATE(12, 3, 3, 5, c(0.1, 0.3), c(0.15, 0.08, 0.02)) #12

# (20)
setwd('../20row')
binSIMULATE(12, 2, 3, 5, c(0.1, 0.3), c(0.1, 0.02, 0.01)) #12

# (21)
setwd('../21row')
binSIMULATE(10, 3, 3, 5, c(0.1, 0.3), c(0.1, 0.02, 0.01)) #10

# (22)
setwd('../22row')
binSIMULATE(10, 3, 3, 5, c(0.1, 0.3), c(0.05, 0.05, 0.02)) #10

# (23)
setwd('../23row')
binSIMULATE(26, 2, 4, 5, c(0.5, 0.7), c(0.4, 0.1, 0.03)) #26

# (24)
setwd('../24row')
binSIMULATE(16, 3, 3, 5, c(0.5, 0.7), c(0.15, 0.08, 0.02)) #16

# (25)
setwd('../25row')
binSIMULATE(12, 2, 4, 5, c(0.5, 0.7), c(0.1, 0.02, 0.01)) #12

# (26)
setwd('../26row')
binSIMULATE(14, 3, 3, 5, c(0.5, 0.7), c(0.05, 0.05, 0.02)) #14

# (27)
setwd('../27row')
binSIMULATE(30, 3, 3, 5, c(0.8, 0.9), c(0.15, 0.08, 0.02)) #30

# (28)
setwd('../28row')
binSIMULATE(22, 2, 4, 5, c(0.8, 0.9), c(0.1, 0.02, 0.01)) #22

# (29)
setwd('../29row')
binSIMULATE(28, 2, 4, 5, c(0.8, 0.9), c(0.05, 0.05, 0.02)) #28

# (30)
setwd('../30row')
binSIMULATE(24, 3, 3, 5, c(0.8, 0.9), c(0.05, 0.05, 0.02)) #24

# Size

# (1)
setwd('../1row')
binSIMULATE(14, 2, 3, 5, c(0.2, 0.2), c(0.4, 0.1, 0.03)) #14

# (2)
setwd('../2row')
binSIMULATE(14, 2, 3, 10, c(0.2, 0.2), c(0.4, 0.1, 0.03)) #14

# (3)
setwd('../3row')
binSIMULATE(14, 2, 4, 5, c(0.2, 0.2), c(0.4, 0.1, 0.03)) #14

# (4)
setwd('../4row')
binSIMULATE(12, 3, 3, 5, c(0.2, 0.2), c(0.4, 0.1, 0.03)) #12

# (5)
setwd('../5row')
binSIMULATE(10, 2, 3, 5, c(0.2, 0.2), c(0.15, 0.08, 0.02)) #10

# (6)
setwd('../6row')
binSIMULATE(10, 2, 3, 10, c(0.2, 0.2), c(0.15, 0.08, 0.02)) #10

# (7)
setwd('../7row')
binSIMULATE(10, 2, 4, 5, c(0.2, 0.2), c(0.15, 0.08, 0.02)) #10

# (8)
setwd('../8row')
binSIMULATE(8, 3, 3, 5, c(0.2, 0.2), c(0.15, 0.08, 0.02)) #8

# (9)
setwd('../9row')
binSIMULATE(8, 2, 3, 5, c(0.2, 0.2), c(0.1, 0.02, 0.01)) #8

# (10)
setwd('../10row')
binSIMULATE(8, 3, 3, 5, c(0.2, 0.2), c(0.1, 0.02, 0.01)) #8

# (11)
setwd('../11row')
binSIMULATE(8, 3, 3, 5, c(0.2, 0.2), c(0.05, 0.05, 0.02)) #8

# (12)
setwd('../12row')
binSIMULATE(22, 2, 3, 5, c(0.1, 0.1), c(0.4, 0.1, 0.03)) #22

# (13)
setwd('../13row')
binSIMULATE(20, 2, 3, 10, c(0.1, 0.1), c(0.4, 0.1, 0.03)) #20

# (14)
setwd('../14row')
binSIMULATE(20, 2, 4, 5, c(0.1, 0.1), c(0.4, 0.1, 0.03)) #20

# (15)
setwd('../15row')
binSIMULATE(16, 3, 3, 5, c(0.1, 0.1), c(0.4, 0.1, 0.03)) #16

# (16)
setwd('../16row')
binSIMULATE(16, 2, 3, 5, c(0.1, 0.1), c(0.15, 0.08, 0.02)) #16

# (17)
setwd('../17row')
binSIMULATE(14, 2, 3, 10, c(0.1, 0.1), c(0.15, 0.08, 0.02)) #14

# (18)
setwd('../18row')
binSIMULATE(14, 2, 4, 5, c(0.1, 0.1), c(0.15, 0.08, 0.02)) #14

# (19)
setwd('../19row')
binSIMULATE(12, 3, 3, 5, c(0.1, 0.1), c(0.15, 0.08, 0.02)) #12

# (20)
setwd('../20row')
binSIMULATE(12, 2, 3, 5, c(0.1, 0.1), c(0.1, 0.02, 0.01)) #12

# (21)
setwd('../21row')
binSIMULATE(10, 3, 3, 5, c(0.1, 0.1), c(0.1, 0.02, 0.01)) #10

# (22)
setwd('../22row')
binSIMULATE(10, 3, 3, 5, c(0.1, 0.1), c(0.05, 0.05, 0.02)) #10

# (23)
setwd('../23row')
binSIMULATE(26, 2, 4, 5, c(0.5, 0.5), c(0.4, 0.1, 0.03)) #26

# (24)
setwd('../24row')
binSIMULATE(16, 3, 3, 5, c(0.5, 0.5), c(0.15, 0.08, 0.02)) #16

# (25)
setwd('../25row')
binSIMULATE(12, 2, 4, 5, c(0.5, 0.5), c(0.1, 0.02, 0.01)) #12

# (26)
setwd('../26row')
binSIMULATE(14, 3, 3, 5, c(0.5, 0.5), c(0.05, 0.05, 0.02)) #14

# (27)
setwd('../27row')
binSIMULATE(30, 3, 3, 5, c(0.8, 0.8), c(0.15, 0.08, 0.02)) #30

# (28)
setwd('../28row')
binSIMULATE(22, 2, 4, 5, c(0.8, 0.8), c(0.1, 0.02, 0.01)) #22

# (29)
setwd('../29row')
binSIMULATE(28, 2, 4, 5, c(0.8, 0.8), c(0.05, 0.05, 0.02)) #28

# (30)
setwd('../30row')
binSIMULATE(24, 3, 3, 5, c(0.8, 0.8), c(0.05, 0.05, 0.02)) #24

setwd('../../..')


