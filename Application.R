##############################################
####### RESHAPE (binary) Application #########
##############################################
library(ggplot2)
library(directlabels)
library(cowplot)

# Function of sample size calculations for four-level designs, when a nominal type I error rate is fixed at
# 5%, with binary outcomes under different link functions
# INPUT:
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
# pc: percentage of clusters assigned to the control arm
# gamma: type II error rate; power = 1-gamma
# link: 1 = logit; 2 = identity;  3 = log
# randomization: randomization level
binSIMULATE <- function(m,k,l,beta,alpha,pc=0.5,gamma=0.2,link=1,randomization=4){
  if (link==1){
    b <- log(beta[2]/(1-beta[2]))-log(beta[1]/(1-beta[1]))
  } else if (link==2){
    b <- beta[2]-beta[1]
  } else if (link==3){
    b <- log(beta[2])-log(beta[1])
  }
  
  lambda1 <- 1-alpha[1]
  lambda2 <- 1+(l-1)*alpha[1]-l*alpha[2]
  lambda3 <- 1+(l-1)*alpha[1]+l*(k-1)*alpha[2]-l*k*alpha[3]
  lambda4 <- 1+(l-1)*alpha[1]+l*(k-1)*alpha[2]+l*k*(m-1)*alpha[3]
  
  if (randomization==4){
    if (link==1){
      sigma2 <- lambda4/(m*k*l)*(1/(pc*beta[1]*(1-beta[1]))+1/((1-pc)*beta[2]*(1-beta[2])))
    } else if (link==2){
      sigma2 <- lambda4/(m*k*l)*((beta[1]*(1-beta[1]))/pc+(beta[2]*(1-beta[2]))/(1-pc))
    } else if (link==3){
      sigma2 <- lambda4/(m*k*l)*((1-beta[1])/(pc*beta[1])+(1-beta[2])/((1-pc)*beta[2]))
    }
  } else if (randomization==3){
    if (link==1){
      sigma2 <- lambda3/(m*k*l)*(1/(pc*beta[1]*(1-beta[1]))+1/((1-pc)*beta[2]*(1-beta[2]))) + 
        (lambda4-lambda3)/(m*k*l)*(1/sqrt(beta[1]*(1-beta[1]))-1/sqrt(beta[2]*(1-beta[2])))^2
    } else if (link==2){
      sigma2 <- lambda3/(m*k*l)*((beta[1]*(1-beta[1]))/pc+(beta[2]*(1-beta[2]))/(1-pc)) +
        (lambda4-lambda3)/(m*k*l)*(sqrt(beta[1]*(1-beta[1]))-sqrt(beta[2]*(1-beta[2])))^2
    } else if (link==3){
      sigma2 <- lambda3/(m*k*l)*((1-beta[1])/(pc*beta[1])+(1-beta[2])/((1-pc)*beta[2])) +
        (lambda4-lambda3)/(m*k*l)*(sqrt((1-beta[1])/beta[1])-sqrt((1-beta[2])/beta[2]))^2
    }
  } else if (randomization==2){
    if (link==1){
      sigma2 <- lambda2/(m*k*l)*(1/(pc*beta[1]*(1-beta[1]))+1/((1-pc)*beta[2]*(1-beta[2]))) + 
        (lambda4-lambda2)/(m*k*l)*(1/sqrt(beta[1]*(1-beta[1]))-1/sqrt(beta[2]*(1-beta[2])))^2
    } else if (link==2){
      sigma2 <- lambda2/(m*k*l)*((beta[1]*(1-beta[1]))/pc+(beta[2]*(1-beta[2]))/(1-pc)) +
        (lambda4-lambda2)/(m*k*l)*(sqrt(beta[1]*(1-beta[1]))-sqrt(beta[2]*(1-beta[2])))^2
    } else if (link==3){
      sigma2 <- lambda2/(m*k*l)*((1-beta[1])/(pc*beta[1])+(1-beta[2])/((1-pc)*beta[2])) +
        (lambda4-lambda2)/(m*k*l)*(sqrt((1-beta[1])/beta[1])-sqrt((1-beta[2])/beta[2]))^2
    }
  } else if (randomization==1){
    if (link==1){
      sigma2 <- lambda1/(m*k*l)*(1/(pc*beta[1]*(1-beta[1]))+1/((1-pc)*beta[2]*(1-beta[2]))) + 
        (lambda4-lambda1)/(m*k*l)*(1/sqrt(beta[1]*(1-beta[1]))-1/sqrt(beta[2]*(1-beta[2])))^2
    } else if (link==2){
      sigma2 <- lambda1/(m*k*l)*((beta[1]*(1-beta[1]))/pc+(beta[2]*(1-beta[2]))/(1-pc)) +
        (lambda4-lambda1)/(m*k*l)*(sqrt(beta[1]*(1-beta[1]))-sqrt(beta[2]*(1-beta[2])))^2
    } else if (link==3){
      sigma2 <- lambda1/(m*k*l)*((1-beta[1])/(pc*beta[1])+(1-beta[2])/((1-pc)*beta[2])) +
        (lambda4-lambda1)/(m*k*l)*(sqrt((1-beta[1])/beta[1])-sqrt((1-beta[2])/beta[2]))^2
    }
  }
  
  n <- (qnorm(0.05/2)+qnorm(gamma))^2/b^2*sigma2
  n <- 2*ceiling(n/2)
  if (n==2){
    nt <- 4
  } else{
    nt <- (qt(0.05/2, df=n-2)+qt(gamma, n-2))^2/b^2*sigma2
  }
  while (n<nt){
    n <- n+2
    nt <- (qt(0.05/2, df=n-2)+qt(gamma, n-2))^2/b^2*sigma2
  }
  t_power <- pt(qt(0.05/2, df=n-2)+abs(b)*sqrt(n/sigma2), df=n-2)
  return(data.frame(p0=beta[1], p1=beta[2], alpha0=alpha[1], alpha1=alpha[2], alpha2=alpha[3], 
                    n=n, m=m, k=k, l=l, t_test=t_power))
}

# Function of sample size calculations using the design effect method, for four-level designs when
# randomization is carried out at the fourth level and a nominal type I error rate is fixed at 5%,
# with binary outcomes under the canonical logit link function
# INPUT:
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
# pc: percentage of clusters assigned to the control arm
# gamma: type II error rate; power = 1-gamma
binDE <- function(m,k,l,beta,alpha,pc=0.5,gamma=0.2){
  b <- log(beta[2]/(1-beta[2]))-log(beta[1]/(1-beta[1]))
  lambda4 <- 1+(l-1)*alpha[1]+l*(k-1)*alpha[2]+l*k*(m-1)*alpha[3]
  sigma2 <- (1/(pc*beta[1]*(1-beta[1]))+1/((1-pc)*beta[2]*(1-beta[2])))
  nmkl <- (qnorm(0.05/2)+qnorm(gamma))^2/b^2*sigma2
  nmkl <- 2*ceiling(nmkl/2)
  df_adjusted <- nmkl/(m*k*l)*lambda4-2
  ntmlk <- (qt(0.05/2, df=df_adjusted)+qt(gamma, df_adjusted))^2/b^2*sigma2
  while (nmkl<ntmlk){
    nmkl <- nmkl+2
    df_adjusted <- nmkl/(m*k*l)*lambda4-2
    nt <- (qt(0.05/2, df=df_adjusted)+qt(gamma, df_adjusted))^2/b^2*sigma2
  }
  t_power_i <- pt(qt(0.05/2, df=df_adjusted)+abs(b)*sqrt(nmkl/sigma2), df=df_adjusted)
  nmlk_de <- nmkl*lambda4
  n_de <- nmlk_de/(m*k*l)
  n_de <- 2*ceiling(n_de/2)
  sigma2_c <- lambda4/(m*k*l)*(1/(pc*beta[1]*(1-beta[1]))+1/((1-pc)*beta[2]*(1-beta[2])))
  t_power <- pt(qt(0.05/2, df=n_de-2)+abs(b)*sqrt(n_de/sigma2_c), df=n_de-2)
  return(data.frame(p0=beta[1], p1=beta[2], alpha0=alpha[1], alpha1=alpha[2], alpha2=alpha[3], 
                    n=n_de, m=m, k=k, l=l, t_test=t_power, lambda4, nmkl, t_ind=t_power_i))
}

# Function of power calculations for four-level designs, with binary outcomes under different link functions
# INPUT:
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
# pc: percentage of clusters assigned to the control arm
# link: 1 = logit; 2 = identity;  3 = log
# randomization: randomization level
binPOWER <- function(n,m,k,l,beta,alpha,pc=0.5,link=1,randomization=4){
  if (link==1){
    b <- log(beta[2]/(1-beta[2]))-log(beta[1]/(1-beta[1]))
  } else if (link==2){
    b <- beta[2]-beta[1]
  } else if (link==3){
    b <- log(beta[2])-log(beta[1])
  }
  
  lambda1 <- 1-alpha[1]
  lambda2 <- 1+(l-1)*alpha[1]-l*alpha[2]
  lambda3 <- 1+(l-1)*alpha[1]+l*(k-1)*alpha[2]-l*k*alpha[3]
  lambda4 <- 1+(l-1)*alpha[1]+l*(k-1)*alpha[2]+l*k*(m-1)*alpha[3]
  
  if (randomization==4){
    if (link==1){
      sigma2 <- lambda4/(m*k*l)*(1/(pc*beta[1]*(1-beta[1]))+1/((1-pc)*beta[2]*(1-beta[2])))
    } else if (link==2){
      sigma2 <- lambda4/(m*k*l)*((beta[1]*(1-beta[1]))/pc+(beta[2]*(1-beta[2]))/(1-pc))
    } else if (link==3){
      sigma2 <- lambda4/(m*k*l)*((1-beta[1])/(pc*beta[1])+(1-beta[2])/((1-pc)*beta[2]))
    }
  } else if (randomization==3){
    if (link==1){
      sigma2 <- lambda3/(m*k*l)*(1/(pc*beta[1]*(1-beta[1]))+1/((1-pc)*beta[2]*(1-beta[2]))) + 
        (lambda4-lambda3)/(m*k*l)*(1/sqrt(beta[1]*(1-beta[1]))-1/sqrt(beta[2]*(1-beta[2])))^2
    } else if (link==2){
      sigma2 <- lambda3/(m*k*l)*((beta[1]*(1-beta[1]))/pc+(beta[2]*(1-beta[2]))/(1-pc)) +
        (lambda4-lambda3)/(m*k*l)*(sqrt(beta[1]*(1-beta[1]))-sqrt(beta[2]*(1-beta[2])))^2
    } else if (link==3){
      sigma2 <- lambda3/(m*k*l)*((1-beta[1])/(pc*beta[1])+(1-beta[2])/((1-pc)*beta[2])) +
        (lambda4-lambda3)/(m*k*l)*(sqrt((1-beta[1])/beta[1])-sqrt((1-beta[2])/beta[2]))^2
    }
  } else if (randomization==2){
    if (link==1){
      sigma2 <- lambda2/(m*k*l)*(1/(pc*beta[1]*(1-beta[1]))+1/((1-pc)*beta[2]*(1-beta[2]))) + 
        (lambda4-lambda2)/(m*k*l)*(1/sqrt(beta[1]*(1-beta[1]))-1/sqrt(beta[2]*(1-beta[2])))^2
    } else if (link==2){
      sigma2 <- lambda2/(m*k*l)*((beta[1]*(1-beta[1]))/pc+(beta[2]*(1-beta[2]))/(1-pc)) +
        (lambda4-lambda2)/(m*k*l)*(sqrt(beta[1]*(1-beta[1]))-sqrt(beta[2]*(1-beta[2])))^2
    } else if (link==3){
      sigma2 <- lambda2/(m*k*l)*((1-beta[1])/(pc*beta[1])+(1-beta[2])/((1-pc)*beta[2])) +
        (lambda4-lambda2)/(m*k*l)*(sqrt((1-beta[1])/beta[1])-sqrt((1-beta[2])/beta[2]))^2
    }
  } else if (randomization==1){
    if (link==1){
      sigma2 <- lambda1/(m*k*l)*(1/(pc*beta[1]*(1-beta[1]))+1/((1-pc)*beta[2]*(1-beta[2]))) + 
        (lambda4-lambda1)/(m*k*l)*(1/sqrt(beta[1]*(1-beta[1]))-1/sqrt(beta[2]*(1-beta[2])))^2
    } else if (link==2){
      sigma2 <- lambda1/(m*k*l)*((beta[1]*(1-beta[1]))/pc+(beta[2]*(1-beta[2]))/(1-pc)) +
        (lambda4-lambda1)/(m*k*l)*(sqrt(beta[1]*(1-beta[1]))-sqrt(beta[2]*(1-beta[2])))^2
    } else if (link==3){
      sigma2 <- lambda1/(m*k*l)*((1-beta[1])/(pc*beta[1])+(1-beta[2])/((1-pc)*beta[2])) +
        (lambda4-lambda1)/(m*k*l)*(sqrt((1-beta[1])/beta[1])-sqrt((1-beta[2])/beta[2]))^2
    }
  }
  
  t_power <- pt(qt(0.05/2, df=n-2)+abs(b)*sqrt(n/sigma2), df=n-2)
  return(data.frame(alpha1=alpha[2], alpha2=alpha[3], t_test=t_power))
}

# Function of power contours for four-level designs, with binary outcomes under different link functions
# INPUT:
# n: number of clusters
# m: number of divisions per cluster
# k: number of participants per division
# l: number of evaluations per participant
# a0: correlation between different evaluations from the same participant
# a1_range: range of correlation between evaluations from different participants but within the same division
# a2_range: range of correlation between evaluations from different divisions but within the same cluster
# beta - marginal means
#        beta_0 - the marginal mean P0 in the control arm
#        beta_1 - the marginal mean P1 in the intervention arm
# pc: percentage of clusters assigned to the control arm
# link: 1 = logit; 2 = identity;  3 = log
# randomization: randomization level
binPLOT <- function(n,m,k,l,a0,a1_range,a2_range,beta,pc=0.5,link=1,randomization=4){
  bin_plot<-NULL
  hali<-NULL
  for (a1 in seq(a1_range[1], a1_range[2], 0.001)){
    for (a2 in seq(a2_range[1], a2_range[2], 0.0005)){
      hali<-rbind(hali, binPOWER(n,m,k,l,beta,c(a0, a1, a2),pc,link,randomization))
    }
  }
  fig<-ggplot()  +
    theme_bw() +
    ggtitle(bquote(alpha[0] == ~.(a0))) +
    xlab(expression(alpha[1])) +
    ylab(expression(alpha[2])) +
    scale_x_continuous(breaks=seq(a1_range[1], a1_range[2], 0.01)) +
    scale_y_continuous(breaks=seq(a2_range[1], a2_range[2], 0.01)) +
    stat_contour(data = hali, aes(x = alpha1, y = alpha2, z = t_test, colour = ..level..), 
                 breaks = round(quantile(hali$t_test, seq(0, 1, 0.1)), 2), size = 1) +
    scale_color_continuous(name = "Power") +
    theme(plot.title = element_text(hjust = 0.5),
          legend.justification=c(1, 0), legend.position=c(1, 0))
  bin_plot<-direct.label(fig, "bottom.pieces")
  return(bin_plot)
}

# ---------------------------------------------------------------------------------

# logit link with randomization at level 4
binSIMULATE(3, 3, 36, c(0.785, 0.88), c(0.05, 0.04, 0.03)) # n = 22, power = 0.8265288
binDE(3, 3, 36, c(0.785, 0.88), c(0.05, 0.04, 0.03)) # n = 22, power = 0.8265288

# logit link with randomization at lower levels
binSIMULATE(3, 3, 36, c(0.785, 0.88), c(0.05, 0.04, 0.03), 0.5, 0.2, 1, 3) # n = 8, power = 0.9177897
binSIMULATE(3, 3, 36, c(0.785, 0.88), c(0.05, 0.04, 0.03), 0.5, 0.2, 1, 2) # n = 6 power = 0.9283434
binSIMULATE(3, 3, 36, c(0.785, 0.88), c(0.05, 0.04, 0.03), 0.5, 0.2, 1, 1) # n = 6 power = 0.966874

# identity link with randomization at different levels
binSIMULATE(3, 3, 36, c(0.785, 0.88), c(0.05, 0.04, 0.03), 0.5, 0.2, 2, 4) # n = 20, power = 0.8009576
binSIMULATE(3, 3, 36, c(0.785, 0.88), c(0.05, 0.04, 0.03), 0.5, 0.2, 2, 3) # n = 8, power = 0.9265738
binSIMULATE(3, 3, 36, c(0.785, 0.88), c(0.05, 0.04, 0.03), 0.5, 0.2, 2, 2) # n = 6, power = 0.9357038
binSIMULATE(3, 3, 36, c(0.785, 0.88), c(0.05, 0.04, 0.03), 0.5, 0.2, 2, 1) # n = 6, power = 0.9704446

# log link with randomization at different levels
binSIMULATE(3, 3, 36, c(0.785, 0.88), c(0.05, 0.04, 0.03), 0.5, 0.2, 3, 4) # n = 22, power = 0.8291002
binSIMULATE(3, 3, 36, c(0.785, 0.88), c(0.05, 0.04, 0.03), 0.5, 0.2, 3, 3) # n = 8, power = 0.9055332
binSIMULATE(3, 3, 36, c(0.785, 0.88), c(0.05, 0.04, 0.03), 0.5, 0.2, 3, 2) # n = 6, power = 0.9063791
binSIMULATE(3, 3, 36, c(0.785, 0.88), c(0.05, 0.04, 0.03), 0.5, 0.2, 3, 1) # n = 6, power = 0.9511328

### power contours
# logit link (Figure 4 in the manuscript)
reshape_1_1 <- binPLOT(22, 3, 3, 36, 0.025, c(0, 0.07), c(0, 0.04), c(0.785, 0.88), 0.5, 1, 4)
reshape_1_2 <- binPLOT(22, 3, 3, 36, 0.05, c(0, 0.07), c(0, 0.04), c(0.785, 0.88), 0.5, 1, 4)
reshape_1_3 <- binPLOT(22, 3, 3, 36, 0.1, c(0, 0.07), c(0, 0.04), c(0.785, 0.88), 0.5, 1, 4)
plot_grid(reshape_1_1, reshape_1_2, reshape_1_3, ncol=3)

# identity link (Web Figure 9 in the supporting information)
reshape_2_1 <- binPLOT(20, 3, 3, 36, 0.025, c(0, 0.07), c(0, 0.04), c(0.785, 0.88), 0.5, 2, 4)
reshape_2_2 <- binPLOT(20, 3, 3, 36, 0.05, c(0, 0.07), c(0, 0.04), c(0.785, 0.88), 0.5, 2, 4)
reshape_2_3 <- binPLOT(20, 3, 3, 36, 0.1, c(0, 0.07), c(0, 0.04), c(0.785, 0.88), 0.5, 2, 4)
plot_grid(reshape_2_1, reshape_2_2, reshape_2_3, ncol=3)

# log link (Web Figure 10 in the supporting information)
reshape_3_1 <- binPLOT(22, 3, 3, 36, 0.025, c(0, 0.07), c(0, 0.04), c(0.785, 0.88), 0.5, 3, 4)
reshape_3_2 <- binPLOT(22, 3, 3, 36, 0.05, c(0, 0.07), c(0, 0.04), c(0.785, 0.88), 0.5, 3, 4)
reshape_3_3 <- binPLOT(22, 3, 3, 36, 0.1, c(0, 0.07), c(0, 0.04), c(0.785, 0.88), 0.5, 3, 4)
plot_grid(reshape_3_1, reshape_3_2, reshape_3_3, ncol=3)



##############################################
####### HALI (continuous) Application ########
##############################################

# Function of sample size calculations for four-level designs, when a nominal type I error rate is fixed at
# 5%, with continuous outcomes under the canonical identity link function
# INPUT:
# m: number of divisions per cluster
# k: number of participants per division
# l: number of evaluations per participant
# beta - marginal means
#        beta_0 - the marginal mean in the control arm
#        beta_1 - the marginal mean in the intervention arm
# var: dispersion parameter
# alpha - correlations 
#         alpha_0 - correlation between different evaluations from the same participant
#         alpha_1 - correlation between evaluations from different participants but within the same division
#         alpha_2 - correlation between evaluations from different divisions but within the same cluster
# pc: percentage of clusters assigned to the control arm
# gamma: type II error rate; power = 1-gamma
# randomization: randomization level
conSIMULATE <- function(m,k,l,beta,var,alpha,pc=0.5,gamma=0.2,randomization=4){
  b <- beta[2]-beta[1]
  
  lambda1 <- 1-alpha[1]
  lambda2 <- 1+(l-1)*alpha[1]-l*alpha[2]
  lambda3 <- 1+(l-1)*alpha[1]+l*(k-1)*alpha[2]-l*k*alpha[3]
  lambda4 <- 1+(l-1)*alpha[1]+l*(k-1)*alpha[2]+l*k*(m-1)*alpha[3]
  
  if (randomization==4){
    sigma2 <- var*lambda4/(pc*(1-pc)*m*k*l)
  } else if (randomization==3){
    sigma2 <- var*lambda3/(pc*(1-pc)*m*k*l)
  } else if (randomization==2){
    sigma2 <- var*lambda2/(pc*(1-pc)*m*k*l)
  } else if (randomization==1){
    sigma2 <- var*lambda1/(pc*(1-pc)*m*k*l)
  }
  
  n <- (qnorm(0.05/2)+qnorm(gamma))^2/b^2*sigma2
  n <- 2*ceiling(n/2)
  if (n==2){
    nt <- 4
  } else{
    nt <- (qt(0.05/2, df=n-2)+qt(gamma, n-2))^2/b^2*sigma2
  }
  while (n<nt){
    n <- n+2
    nt <- (qt(0.05/2, df=n-2)+qt(gamma, n-2))^2/b^2*sigma2
  }
  t_power <- pt(qt(0.05/2, df=n-2)+abs(b)*sqrt(n/sigma2), df=n-2)
  return(data.frame(p0=beta[1], p1=beta[2], alpha0=alpha[1], alpha1=alpha[2], alpha2=alpha[3], 
                    n=n, m=m, k=k, l=l, t_test=t_power))
}

# Function of power calculations for four-level designs, with continuous outcomes under the canonical
# identity link function
# INPUT:
# n: number of clusters
# m: number of divisions per cluster
# k: number of participants per division
# l: number of evaluations per participant
# beta - marginal means
#        beta_0 - the marginal mean P0 in the control arm
#        beta_1 - the marginal mean P1 in the intervention arm
# var: dispersion parameter
# alpha - correlations 
#         alpha_0 - correlation between different evaluations from the same participant
#         alpha_1 - correlation between evaluations from different participants but within the same division
#         alpha_2 - correlation between evaluations from different divisions but within the same cluster
# pc: percentage of clusters assigned to the control arm
# randomization: randomization level
conPOWER <- function(n,m,k,l,beta,var,alpha,pc=0.5,randomization=4){
  b <- beta[2]-beta[1]
  
  lambda1 <- 1-alpha[1]
  lambda2 <- 1+(l-1)*alpha[1]-l*alpha[2]
  lambda3 <- 1+(l-1)*alpha[1]+l*(k-1)*alpha[2]-l*k*alpha[3]
  lambda4 <- 1+(l-1)*alpha[1]+l*(k-1)*alpha[2]+l*k*(m-1)*alpha[3]
  
  if (randomization==4){
    sigma2 <- var*lambda4/(pc*(1-pc)*m*k*l)
  } else if (randomization==3){
    sigma2 <- var*lambda3/(pc*(1-pc)*m*k*l)
  } else if (randomization==2){
    sigma2 <- var*lambda2/(pc*(1-pc)*m*k*l)
  } else if (randomization==1){
    sigma2 <- var*lambda1/(pc*(1-pc)*m*k*l)
  }
  
  t_power <- pt(qt(0.05/2, df=n-2)+abs(b)*sqrt(n/sigma2), df=n-2)
  return(data.frame(alpha1=alpha[2], alpha2=alpha[3], t_test=t_power))
}

# Function of power contours for four-level designs, with with continuous outcomes under the canonical
# identity link function
# INPUT:
# n: number of clusters
# m: number of divisions per cluster
# k: number of participants per division
# l: number of evaluations per participant
# beta - marginal means
#        beta_0 - the marginal mean P0 in the control arm
#        beta_1 - the marginal mean P1 in the intervention arm
# var: dispersion parameter
# a0: correlation between different evaluations from the same participant
# a1_range: range of correlation between evaluations from different participants but within the same division
# a2_range: range of correlation between evaluations from different divisions but within the same cluster
# pc: percentage of clusters assigned to the control arm
# randomization: randomization level
conPLOT <- function(n,m,k,l,beta,var,a0,a1_range,a2_range,pc=0.5,randomization=4){
  con_plot<-NULL
  hali<-NULL
  for (a1 in seq(a1_range[1], a1_range[2], 0.001)){
    for (a2 in seq(a2_range[1], a2_range[2], 0.0005)){
      hali<-rbind(hali, conPOWER(n,m,k,l,beta,var,c(a0, a1, a2),pc,randomization))
    }
  }
  fig<-ggplot()  +
    theme_bw() +
    ggtitle(bquote(alpha[0] == ~.(a0))) +
    xlab(expression(alpha[1])) +
    ylab(expression(alpha[2])) +
    ylim(c(a2_range[1], a2_range[2])) +
    scale_x_continuous(breaks=seq(a1_range[1], a1_range[2], 0.01)) +
    stat_contour(data = hali, aes(x = alpha1, y = alpha2, z = t_test, colour = ..level..), 
                 breaks = round(quantile(hali$t_test, seq(0, 1, 0.1)), 2), size = 1) +
    scale_color_continuous(name = "Power") +
    theme(plot.title = element_text(hjust = 0.5),
          legend.justification=c(1, 0), legend.position=c(1, 0))
  con_plot<-direct.label(fig, "bottom.pieces")
  return(con_plot)
}

# ---------------------------------------------------------------------------------

# Effect size of 0.19 standard deviation (SD) with randomization at level 4
conSIMULATE(4, 25, 2, c(11, 11.19), 1, c(0.445, 0.104, 0.008), 0.5, 0.2, 4) # n = 36, power = 0.8087343

# Effect size of 0.25 standard deviation (SD) with randomization at level 4
conSIMULATE(4, 25, 2, c(11, 11.25), 1, c(0.445, 0.104, 0.008), 0.5, 0.2, 4) # n = 22, power = 0.8143021

# Effect size of 0.19 standard deviation (SD) with randomization at lower levels
conSIMULATE(4, 25, 2, c(11, 11.19), 1, c(0.445, 0.104, 0.008), 0.5, 0.2, 3) # n = 30, power = 0.8240136
conSIMULATE(4, 25, 2, c(11, 11.19), 1, c(0.445, 0.104, 0.008), 0.5, 0.2, 2) # n = 8, power = 0.8151828

# Effect size of 0.25 standard deviation (SD) with randomization at lower levels
conSIMULATE(4, 25, 2, c(11, 11.25), 1, c(0.445, 0.104, 0.008), 0.5, 0.2, 3) # n = 18, power = 0.8175463
conSIMULATE(4, 25, 2, c(11, 11.25), 1, c(0.445, 0.104, 0.008), 0.5, 0.2, 2) # n = 6, power = 0.8366901

### power contours
# Effect size of 0.19 standard deviation (SD) and n = 36 (Figure 5 in the manuscript)
hali_1_1 <- conPLOT(36, 4, 25, 2, c(11, 11.19), 1, 0.4, c(0.05, 0.15), c(0, 0.02))
hali_1_2 <- conPLOT(36, 4, 25, 2, c(11, 11.19), 1, 0.445, c(0.05, 0.15), c(0, 0.02))
hali_1_3 <- conPLOT(36, 4, 25, 2, c(11, 11.19), 1, 0.5, c(0.05, 0.15), c(0, 0.02))
plot_grid(hali_1_1, hali_1_2, hali_1_3, ncol=3)

# Effect size of 0.25 standard deviation (SD) and n = 22 (Web Figure 11 in the supporting information)
hali_2_1 <- conPLOT(22, 4, 25, 2, c(11, 11.25), 1, 0.4, c(0.05, 0.15), c(0, 0.02))
hali_2_2 <- conPLOT(22, 4, 25, 2, c(11, 11.25), 1, 0.445, c(0.05, 0.15), c(0, 0.02))
hali_2_3 <- conPLOT(22, 4, 25, 2, c(11, 11.25), 1, 0.5, c(0.05, 0.15), c(0, 0.02))
plot_grid(hali_2_1, hali_2_2, hali_2_3, ncol=3)

# Effect size of 0.25 standard deviation (SD) and n = 26 (Web Figure 12 in the supporting information)
hali_3_1 <- conPLOT(26, 4, 25, 2, c(11, 11.25), 1, 0.4, c(0.05, 0.15), c(0, 0.02))
hali_3_2 <- conPLOT(26, 4, 25, 2, c(11, 11.25), 1, 0.445, c(0.05, 0.15), c(0, 0.02))
hali_3_3 <- conPLOT(26, 4, 25, 2, c(11, 11.25), 1, 0.5, c(0.05, 0.15), c(0, 0.02))
plot_grid(hali_3_1, hali_3_2, hali_3_3, ncol=3)


