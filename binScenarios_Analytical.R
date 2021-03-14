############################################################################################
############################## Design binary simulations ###################################
############################################################################################
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
# gamma: type II error rate; power = 1-gamma
############################################################################################
library(openxlsx)

# Function of sample size calculations for four-level CRTs, when randomization is carried out
# at the fourth level and a nominal type I error rate is fixed at 5%, with binary outcomes 
# under the canonical logit link function
# (This binSIMULATE differs from the actual simulation program.
# This is only used to calculate the power analytically)
binSIMULATE <- function(m,k,l,beta,alpha,pc=0.5,gamma=0.2){
  b <- log(beta[2]/(1-beta[2]))-log(beta[1]/(1-beta[1]))
  lambda4 <- 1+(l-1)*alpha[1]+l*(k-1)*alpha[2]+l*k*(m-1)*alpha[3]
  sigma2 <- lambda4/(m*k*l)*(1/(pc*beta[1]*(1-beta[1]))+1/((1-pc)*beta[2]*(1-beta[2])))
  n <- (qnorm(0.05/2)+qnorm(gamma))^2/b^2*sigma2
  n <- 2*ceiling(n/2)
  nt <- (qt(0.05/2, df=n-2)+qt(gamma, n-2))^2/b^2*sigma2
  while (n<nt){
    n <- n+2
    nt <- (qt(0.05/2, df=n-2)+qt(gamma, n-2))^2/b^2*sigma2
  }
  t_power <- pt(qt(0.05/2, df=n-2)+abs(b)*sqrt(n/sigma2), df=n-2)
  return(data.frame(p0=beta[1], p1=beta[2], alpha0=alpha[1], alpha1=alpha[2], alpha2=alpha[3], 
                    n=n, m=m, k=k, l=l, t_test=t_power))
}

# ---------------------------------------------------------------------------------

# Analytical power
tPower <- NULL

tPower <- rbind(tPower, binSIMULATE(2, 3, 5, c(0.2, 0.5), c(0.4, 0.1, 0.03))) # n=14
tPower <- rbind(tPower, binSIMULATE(2, 3, 10, c(0.2, 0.5), c(0.4, 0.1, 0.03))) # n=14
tPower <- rbind(tPower, binSIMULATE(2, 4, 5, c(0.2, 0.5), c(0.4, 0.1, 0.03))) # n=14
tPower <- rbind(tPower, binSIMULATE(3, 3, 5, c(0.2, 0.5), c(0.4, 0.1, 0.03))) # n=12
tPower <- rbind(tPower, binSIMULATE(2, 3, 5, c(0.2, 0.5), c(0.15, 0.08, 0.02))) # n=10
tPower <- rbind(tPower, binSIMULATE(2, 3, 10, c(0.2, 0.5), c(0.15, 0.08, 0.02))) # n=10
tPower <- rbind(tPower, binSIMULATE(2, 4, 5, c(0.2, 0.5), c(0.15, 0.08, 0.02))) # n=10
tPower <- rbind(tPower, binSIMULATE(3, 3, 5, c(0.2, 0.5), c(0.15, 0.08, 0.02))) # n=8
tPower <- rbind(tPower, binSIMULATE(2, 3, 5, c(0.2, 0.5), c(0.1, 0.02, 0.01))) # n=8
tPower <- rbind(tPower, binSIMULATE(3, 3, 5, c(0.2, 0.5), c(0.1, 0.02, 0.01))) # n=8
tPower <- rbind(tPower, binSIMULATE(3, 3, 5, c(0.2, 0.5), c(0.05, 0.05, 0.02))) # n=8

tPower <- rbind(tPower, binSIMULATE(2, 3, 5, c(0.1, 0.3), c(0.4, 0.1, 0.03))) # n=22
tPower <- rbind(tPower, binSIMULATE(2, 3, 10, c(0.1, 0.3), c(0.4, 0.1, 0.03))) # n=20
tPower <- rbind(tPower, binSIMULATE(2, 4, 5, c(0.1, 0.3), c(0.4, 0.1, 0.03))) # n=20
tPower <- rbind(tPower, binSIMULATE(3, 3, 5, c(0.1, 0.3), c(0.4, 0.1, 0.03))) # n=16
tPower <- rbind(tPower, binSIMULATE(2, 3, 5, c(0.1, 0.3), c(0.15, 0.08, 0.02))) # n=16
tPower <- rbind(tPower, binSIMULATE(2, 3, 10, c(0.1, 0.3), c(0.15, 0.08, 0.02))) # n=14
tPower <- rbind(tPower, binSIMULATE(2, 4, 5, c(0.1, 0.3), c(0.15, 0.08, 0.02))) # n=14
tPower <- rbind(tPower, binSIMULATE(3, 3, 5, c(0.1, 0.3), c(0.15, 0.08, 0.02))) # n=12
tPower <- rbind(tPower, binSIMULATE(2, 3, 5, c(0.1, 0.3), c(0.1, 0.02, 0.01))) # n=12
tPower <- rbind(tPower, binSIMULATE(3, 3, 5, c(0.1, 0.3), c(0.1, 0.02, 0.01))) # n=10
tPower <- rbind(tPower, binSIMULATE(3, 3, 5, c(0.1, 0.3), c(0.05, 0.05, 0.02))) # n=10

tPower <- rbind(tPower, binSIMULATE(2, 4, 5, c(0.5, 0.7), c(0.4, 0.1, 0.03))) # n=26
tPower <- rbind(tPower, binSIMULATE(3, 3, 5, c(0.5, 0.7), c(0.15, 0.08, 0.02))) # n=16
tPower <- rbind(tPower, binSIMULATE(2, 4, 5, c(0.5, 0.7), c(0.1, 0.02, 0.01))) # n=12
tPower <- rbind(tPower, binSIMULATE(3, 3, 5, c(0.5, 0.7), c(0.05, 0.05, 0.02))) # n=14

tPower <- rbind(tPower, binSIMULATE(3, 3, 5, c(0.8, 0.9), c(0.15, 0.08, 0.02))) # n=30
tPower <- rbind(tPower, binSIMULATE(2, 4, 5, c(0.8, 0.9), c(0.1, 0.02, 0.01))) # n=22
tPower <- rbind(tPower, binSIMULATE(2, 4, 5, c(0.8, 0.9), c(0.05, 0.05, 0.02))) # n=28
tPower <- rbind(tPower, binSIMULATE(3, 3, 5, c(0.8, 0.9), c(0.05, 0.05, 0.02))) # n=24

write.xlsx(tPower, file = "binResults/pred_power_bin.xlsx", row.names = FALSE)


