##########################################################################
# Simulate correlated binary outcomes in a balanced four-level CRT

# Ref: Qaqish, B. F. (2003). A family of multivariate binary distributions 
# for simulating correlated binary variables. Biometrika 90, 455-463.

# INPUT
# n: Number of clusters
# m: Number of divisions per cluster
# k: Number of participants per division
# l: Number of evaluations per participant
# beta: Vector of means=c(mu_0, mu_1)
#       mu_0=population mean of the control arm
#       mu_1=population mean of the intervention arm
# alpha: Vector of correlations 
#        alpha_0=within-participant correlation
#        alpha_1=within-division correlation
#        alpha_2=inter-division correlation
##########################################################################

binGEN<-function(n,m,k,l,beta,alpha){
  
  ########################################################
  # Create extended nested exchangeable correlation matrix.
  ########################################################
  
  enxch<-function(alpha){
    rho<-alpha[1]
    nu<-alpha[2]
    lambda<-alpha[3]
    em1<-(1-rho)*diag(m*k*l)
    em2<-(rho-nu)*kronecker(diag(m*k),matrix(1,l,l))
    em3<-(nu-lambda)*kronecker(diag(m),matrix(1,k*l,k*l))
    em4<-lambda*matrix(1,m*k*l,m*k*l)
    return(em1+em2+em3+em4)
  }
  
  ########################################################
  # a[1:n, 1:n] is the input covariance matrix of Y[1:n].
  # Returns  b[1:n,1:n] such that b[1:t-1, t] are the 
  # slopes for regression of y[t] on y[1:t-1], for t=2:n.
  # Diagonals and lower half of b[,] are copied from a[,].
  # a[,] is assumed +ve definite symmetric, not checked.
  ########################################################
  
  allreg<-function(a){
    n<-nrow(a)
    b<-a
    for(t in 2:n){
      t1<-t-1
      gt<-a[1:t1,1:t1]
      st<-a[1:t1,t]
      bt<-solve(gt,st)
      b[1:t1,t]<-bt
    }
    return(b)
  }
  
  ########################################################
  # returns variance matrix of binary variables with mean
  # vector u[] and corr matrix r[,].
  ########################################################
  
  cor2var<-function(r,u){
    s<-diag(sqrt(u*(1-u)))
    return(s%*%r%*%s)
  }
  
  ########################################################
  # r[1:n, 1:n] is the corr mtx
  # u[1:n] is the mean of a binary vector
  # checks that pairwise corrs are in-range for the given u[]
  # only upper half of r[,] is checked.
  # return 0 if ok
  # return 1 if out of range
  ########################################################
  
  chkbinc<-function(r,u){
    n<-length(u)
    s<-sqrt(u*(1-u))
    for(i in 1:(n-1)){
      for(j in (i+1):n){
        uij<-u[i]*u[j]+r[i,j]*s[i]*s[j]
        ok<-((uij <= min(u[i], u[j])) & (uij >= max(0, u[i]+u[j]-1)))
        if(!ok) {return(1)}
      }
    }
    return(0)
  }
  
  ########################################################
  # Multivariate Binary Simulation by Linear Regression.
  # Simulate a single vector.
  # Returns a simulated binary random vector y[1:n] with mean 
  # u[1:n] and regression coefs matrix b[1:n,1:n] (obtained 
  # by calling allreg() above).
  # y[] and u[] are column vectors.
  # Returns -1 if the cond. linear family not reproducible
  ########################################################
  
  mbslr1<-function(b,u){
    n<-nrow(b)
    y<-rep(-1,n)
    y[1]<-rbinom(1,1,u[1])
    for(i in 2:n){
      i1<-i-1
      r<-y[1:i1]-u[1:i1]              # residuals
      ci<-u[i]+sum(r*b[1:i1,i])       # cond.mean
      if(ci < 0 | ci > 1){
        stop(paste("mbslr1: ERROR:",ci))
        return(-1)
      }
      y[i]<-rbinom(1,1,ci)
    }
    return(y)
  }
  
  # Simulate correlated binary outcomes
  y<-NULL
  r<-enxch(alpha)
  for(i in 1:2){
    u_c<-beta[i]
    u<-rep(u_c,m*k*l)
    v<-cor2var(r,u)                   # v <- cov matrix
    oor<-chkbinc(r,u)                 # corr out of range?
    if(oor){
      stop("ERROR: Corr out of range for given mean")
    }
    b<-allreg(v)     # prepare coeffs
    for (j in 1:(n/2)){
      y<-cbind(y,mbslr1(b,u))         # simulate data matrix
    }
  }
  
  # Return simulated data matrix
  return(y)
}

