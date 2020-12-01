########################################################################################################
#             Examples of binary simulations under four-level CRTs
#
# Menu:
# 0. Design binary simulations
# 1. GEE analyses using the true correlation structure with MAEE under balanced four-level CRTs
# 2. GEE analyses using an independence working correlation matrix under balanced four-level CRTs
# 3. GEE analyses using the true correlation structure with MAEE under unbalanced four-level CRTs
# 4. GEE analyses using an independence working correlation matrix under unbalanced four-level CRTs
#
# Parameters:
# n: Number of clusters
# m: Number of divisions per cluster
# l: Number of subjects per division
# k: Number of evaluations per subject, under balanced designs;
#    mean number of evaluations per subject, under unbalanced designs
# k_cv: coefficient of variance of k, under unbalanced designs
#       (not required under balanced designs)
# beta - effect size on the link function scale
# alpha - correlations 
#         alpha_0 - correlation between different evaluations from the same subject
#         alpha_1 - correlation between evaluations from different subjects but within the same division
#         alpha_2 - correlation between evaluations from different divisions but within the same cluster
########################################################################################################


##################################################################################
######### 0. Design binary simulations
##################################################################################

# This binSIMULATE differs from the actual simulation program.
# This is only used to calculate the power analytically.
binSIMULATE <- function(m,l,k,beta,alpha,pc=0.5,gamma=0.2){
  b <- log(beta[2]/(1-beta[2]))-log(beta[1]/(1-beta[1]))
  lambda4 <- 1+(k-1)*alpha[1]+k*(l-1)*alpha[2]+k*l*(m-1)*alpha[3]
  sigma2 <- lambda4/(m*l*k)*(1/(pc*beta[1]*(1-beta[1]))+1/((1-pc)*beta[2]*(1-beta[2])))
  n <- (qnorm(0.05/2)+qnorm(gamma))^2/b^2*sigma2
  n <- 2*ceiling(n/2)
  nt <- (qt(0.05/2, df=n-2)+qt(gamma, n-2))^2/b^2*sigma2
  while (n<nt){
    n <- n+2
    nt <- (qt(0.05/2, df=n-2)+qt(gamma, n-2))^2/b^2*sigma2
  }
  t_power <- pt(qt(0.05/2, df=n-2)+abs(b)*sqrt(n/sigma2), df=n-2)
  return(data.frame(p0=beta[1], p1=beta[2], alpha0=alpha[1], alpha1=alpha[2], alpha2=alpha[3], 
                    n=n, m=m, l=l, k=k, t_test=t_power))
}

# --------------------------------------------------------------------------------

# Power
binSIMULATE(2, 3, 5, c(0.2, 0.5), c(0.4, 0.1, 0.03)) # n=14, power=0.817


##################################################################################
######### 1. GEE analyses using the true correlation structure with MAEE
#########    under balanced four-level CRTs
##################################################################################
binSIMULATE<-function(n,m,l,k,beta,alpha){
  
  # Source program code
  source("binGEN.R")
  source("binMAEE.R")
  
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
  X<-cbind(rep(1,n*m*l*k),c(rep(0,n*m*l*k/2),rep(1,n*m*l*k/2)))
  
  # Create Z matrix
  CREATEZ<-function(n,m,l,k){
    rho<-1
    nu<-2
    lambda<-3
    em1<-(1-rho)*diag(m*l*k)
    em2<-(rho-nu)*kronecker(diag(m*l),matrix(1,k,k))
    em3<-(nu-lambda)*kronecker(diag(m),matrix(1,l*k,l*k))
    em4<-lambda*matrix(1,m*l*k,m*l*k)
    POS<-em1+em2+em3+em4
    zrow<-diag(3)
    z_c<-NULL
    for(jp in 1:(m*l*k-1)){
      for(t in (jp+1):(m*l*k)){
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
  Z<-CREATEZ(n,m,l,k) # large matrix
  
  # Create id
  id<-rep(1:n,each=m*l*k)
  clsize<-rep(m*l*k,n)
  
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
      y<-try(binGEN(n,m,l,k,beta,alpha),silent=TRUE)
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

# --------------------------------------------------------------------------------

setwd('binRData/example1')
# Power
binSIMULATE(14, 2, 3, 5, c(0.2, 0.5), c(0.4, 0.1, 0.03))
# Size
binSIMULATE(14, 2, 3, 5, c(0.2, 0.2), c(0.4, 0.1, 0.03))


##################################################################################
######### 2. GEE analyses using an independence working correlation matrix
#########    under balanced four-level CRTs
##################################################################################


binSIMULATE<-function(n,m,l,k,beta,alpha){
  
  # Source program code
  source("binGEN.R")
  source("LogitBinBCV.R")
  
  require(gee)
  
  # Store results
  nrep<-1000
  nerror<-0
  results<-matrix(NA,nrep,16)
  colnames(results)<-c("nSIM","pzeroflag","pminflag","poneflag","pmaxflag",
                       "beta","MB","ROB","KC","MD","FG","MBN","AVG","K","p","error")
  results[,14]<-n
  results[,15]<-2
  results[,16]<-0
  
  # Create X matrix
  X<-cbind(rep(1,n*m*l*k),c(rep(0,n*m*l*k/2),rep(1,n*m*l*k/2)))
  
  # Create id
  id<-rep(1:n,each=m*l*k)
  clsize<-rep(m*l*k,n)
  
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
      y<-try(binGEN(n,m,l,k,beta,alpha),silent=TRUE)
      if(class(y)[1]=="try-error"){passGEN<-1;next}
      passGEN<-0
    }
    y<-c(y)
    results[i,1]<-nSIM
    
    # GEE program
    mydata<-data.frame(cbind(y,X,id))
    outfit=try(gee(y~-1+X,id=id,data=mydata,family=binomial(link="logit"),corstr="independence"),silent=TRUE)
    if(class(outfit)[1]=="try-error" | is.null(outfit)==1){
      nerror<-nerror+1
      results[i,16]<-1
      next
    }
    beta_est=as.numeric(outfit$coefficients)
    out=try(LogitBinBCV(y=y,X=X,beta=beta_est,id=id),silent=TRUE)
    if(class(out)[1]=="try-error" | is.null(out)==1){
      nerror<-nerror+1
      results[i,16]<-1
      next
    }
    
    # Save results
    results[i,(6:12)]<-out$outbeta[2,]
    
    # Loop control
    if(i %% 20 == 0){print(paste("Iteration",i,"NERROR =",nerror))}
  }
  
  results[,13]<-rowMeans(results[,9:10], na.rm=TRUE)
  
  # save data sets and clear workspace
  if(beta[2]-beta[1]==0){
    name<-"Size"
  }else{
    name<-"Power"}
  save(results,file=paste0(name,"_beta_ind.RData"))
}

# ---------------------------------------------------------------------------------

setwd('binRData/example2')
# Power
binSIMULATE(14, 2, 3, 5, c(0.2, 0.5), c(0.4, 0.1, 0.03))
# Size
binSIMULATE(14, 2, 3, 5, c(0.2, 0.2), c(0.4, 0.1, 0.03))


##################################################################################
######### 3. GEE analyses using the true correlation structure with MAEE
#########    under unbalanced four-level CRTs
##################################################################################
binSIMULATE<-function(n,m,l,k,k_cv,beta,alpha){
  
  # Source program code
  source("binGEN_var.R")
  source("binMAEE.R")
  
  # Store results
  nrep<-1000
  nerror<-0
  results<-matrix(NA,nrep,17)
  colnames(results)<-c("nSIM","pzeroflag","pminflag","poneflag","pmaxflag",
                       "beta","MB","ROB","KC","MD","FG","MBN","AVG","K","p","error","clsize.mean")
  results[,14]<-n
  results[,15]<-2
  results[,16]<-0
  
  results0<-matrix(NA,nrep,8)
  results1<-matrix(NA,nrep,8)
  results2<-matrix(NA,nrep,8)
  colnames(results0)<-c("alpha0","ROB","KC","MD","FG","MBN","AVG","error")
  colnames(results1)<-c("alpha1","ROB","KC","MD","FG","MBN","AVG","error")
  colnames(results2)<-c("alpha2","ROB","KC","MD","FG","MBN","AVG","error")
  
  # Create Z matrix
  CREATEZ<-function(n,m,l,kv){
    BEGINEND=function(n){
      last=cumsum(n)
      first=last-n+1
      return(cbind(first,last))
    }
    
    Z<- NULL
    for(i in 1:n){
      cv<-kv[((i-1)*m*l+1):(i*m*l)]
      rho<-1
      nu<-2
      lambda<-3
      locr<-BEGINEND(cv)
      em<-matrix(lambda,nrow=sum(cv),ncol=sum(cv))
      for (j in 1:m){
        em[locr[((j-1)*l+1),1]:locr[(j*l),2],locr[((j-1)*l+1),1]:locr[(j*l),2]]<-nu
      }
      for (jj in 1:length(cv)){
        em[locr[jj,1]:locr[jj,2],locr[jj,1]:locr[jj,2]]<-rho
      }
      diag(em)<-1
      zrow<-diag(3)
      z_c<-NULL
      for(jp in 1:(sum(cv)-1)){
        for(t in (jp+1):(sum(cv))){
          z_c<-rbind(z_c,zrow[em[jp,t],])
        }
      }
      Z<-rbind(Z,z_c)
    }
    rm(z_c)
    return(Z)
  }
  
  # For Loop
  j<-0
  gamma_cv=k_cv
  gamma_a=gamma_cv^(-2)
  gamma_b=1/(k*gamma_cv^2)
  for(i in 1:nrep){
    # Generate cluster size
    set.seed(j)
    clsize<-round(rgamma(n*m*l,shape=gamma_a,rate=gamma_b))
    clsize[clsize<=2]=2
    cs.c<-sum(clsize[1:(length(clsize)/2)])
    cs.t<-sum(clsize[(length(clsize)/2+1):length(clsize)])
    results[i,17]<-mean(clsize)
    
    # Create X matrix
    X<-cbind(rep(1,sum(clsize)),c(rep(0,cs.c),rep(1,cs.t)))
    
    # Create id
    locid<-rep(0,n)
    for (ii in 1:n){
      locid[ii]<-sum(clsize[((ii-1)*m*l+1):(ii*m*l)])
    }
    id<-rep(1:n,locid)
    
    # Create Z matrix
    Z<-CREATEZ(n,m,l,clsize) # large matrix
    
    # Generate outcome and check; store number of simulated data sets
    passGEN<-1
    nSIM<-0
    while(passGEN==1){
      nSIM<-nSIM+1
      j<-j+1
      set.seed(j)
      y<-try(binGEN_var(n,m,l,clsize,beta,alpha),silent=TRUE)
      if(class(y)[1]=="try-error"){passGEN<-1;next}
      passGEN<-0
    }
    results[i,1]<-nSIM
    
    # MAEE program
    out=try(binMAEE(y=y,X=X,id=id,n=locid,Z=Z,maxiter=25,epsilon=0.001,printrange="NO",shrink="ALPHA",makevone="NO"),silent=TRUE)
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
  save(results,file=paste0(name,"_beta_var.RData"))
  save(results0,file=paste0(name,"_alpha0_var.RData"))
  save(results1,file=paste0(name,"_alpha1_var.RData"))
  save(results2,file=paste0(name,"_alpha2_var.RData"))
}

# --------------------------------------------------------------------------------

setwd('binRData/example3')
# Power
binSIMULATE(14, 2, 3, 5, 0.25, c(0.2, 0.5), c(0.4, 0.1, 0.03))
# Size
binSIMULATE(14, 2, 3, 5, 0.25, c(0.2, 0.2), c(0.4, 0.1, 0.03))


##################################################################################
######### 4. GEE analyses using an independence working correlation matrix
#########    under unbalanced four-level CRTs
##################################################################################

binSIMULATE<-function(n,m,l,k,k_cv,beta,alpha){
  
  # Source program code
  source("binGEN_var.R")
  source("LogitBinBCV.R")
  
  require(gee)
  
  # Store results
  nrep<-1000
  nerror<-0
  results<-matrix(NA,nrep,17)
  colnames(results)<-c("nSIM","pzeroflag","pminflag","poneflag","pmaxflag",
                       "beta","MB","ROB","KC","MD","FG","MBN","AVG","K","p","error","clsize.mean")
  results[,14]<-n
  results[,15]<-2
  results[,16]<-0
  
  # For Loop
  j<-0
  gamma_cv=k_cv
  gamma_a=gamma_cv^(-2)
  gamma_b=1/(k*gamma_cv^2)
  for(i in 1:nrep){
    # Generate cluster size
    set.seed(j)
    clsize<-round(rgamma(n*m*l,shape=gamma_a,rate=gamma_b))
    clsize[clsize<=2]=2
    cs.c<-sum(clsize[1:(length(clsize)/2)])
    cs.t<-sum(clsize[(length(clsize)/2+1):length(clsize)])
    results[i,17]<-mean(clsize)
    
    # Create X matrix
    X<-cbind(rep(1,sum(clsize)),c(rep(0,cs.c),rep(1,cs.t)))
    
    # Create id
    locid<-rep(0,n)
    for (ii in 1:n){
      locid[ii]<-sum(clsize[((ii-1)*m*l+1):(ii*m*l)])
    }
    id<-rep(1:n,locid)
    
    # Generate outcome and check; store number of simulated data sets
    passGEN<-1
    nSIM<-0
    while(passGEN==1){
      nSIM<-nSIM+1
      j<-j+1
      set.seed(j)
      y<-try(binGEN_var(n,m,l,clsize,beta,alpha),silent=TRUE)
      if(class(y)[1]=="try-error"){passGEN<-1;next}
      passGEN<-0
    }
    #y<-c(y)
    results[i,1]<-nSIM
    
    # GEE program
    mydata<-data.frame(cbind(y,X,id))
    outfit=try(gee(y~-1+X,id=id,data=mydata,family=binomial(link="logit"),corstr="independence"),silent=TRUE)
    if(class(outfit)[1]=="try-error" | is.null(outfit)==1){
      nerror<-nerror+1
      results[i,16]<-1
      next
    }
    beta_est=as.numeric(outfit$coefficients)
    out=try(LogitBinBCV(y=y,X=X,beta=beta_est,id=id),silent=TRUE)
    if(class(out)[1]=="try-error" | is.null(out)==1){
      nerror<-nerror+1
      results[i,16]<-1
      next
    }
    
    # Save results
    results[i,(6:12)]<-out$outbeta[2,]
    
    # Loop control
    if(i %% 20 == 0){print(paste("Iteration",i,"NERROR =",nerror))}
  }
  
  results[,13]<-rowMeans(results[,9:10], na.rm=TRUE)
  
  # save data sets and clear workspace
  if(beta[2]-beta[1]==0){
    name<-"Size"
  }else{
    name<-"Power"}
  save(results,file=paste0(name,"_beta_ind_var.RData"))
}

# --------------------------------------------------------------------------------

setwd('binRData/example4')
# Power
binSIMULATE(14, 2, 3, 5, 0.25, c(0.2, 0.5), c(0.4, 0.1, 0.03))
# Size
binSIMULATE(14, 2, 3, 5, 0.25, c(0.2, 0.2), c(0.4, 0.1, 0.03))


