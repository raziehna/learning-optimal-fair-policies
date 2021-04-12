
rm(list = ls(all=TRUE))
cat("\f")

library(mvtnorm)
library(nloptr) # nloptr()
library(rpart)
library(rpart.plot) 

set.seed(0)

setwd("~/Desktop/Algorithmic Fairness/FairPolicy/github/COMPAS/")

source("compute_effects_compas.R")
source("optimization_nloptr_compas.R")
source("compas_Qlearning.R")
source("compas_Qlearning_unfair.R") 

# ++++++++++++++++++++++
# Generate data
# ++++++++++++++++++++++

dat <- read.csv("compas_preprocessed.csv", sep = ",")

n = nrow(dat)

# (X) +++++++++++++++++++++++
# X1: age, X2: gender
# Male (X2 = 1), Female (X2 = 0)
# Young (age < 35) (X1 = 0), Old (X1 = 1) 
X1 = dat$age 
X2 = dat$gender
# X1[X1 ==0] = -1
# X2[X2 == 0] = -1


# (S) +++++++++++++++++++++++ 
# S: race
# Black (S = 1), White (S = 0)
S = dat$race
# S[S == 0] = -1


# (M) +++++++++++++++++++++++ 
# M: criminal history (M | X, S) 
# M = 1 if priors_count is greater than 1 
M = dat$priors 
# M[M == 0] = -1

fmla_M = as.formula(M ~ X1 + X2 + S + SX1 + SX2)
inter_M_S = data.frame(X1, X2) # in M model what interacts with S? 


# (A) +++++++++++++++++++++++
# A: compas score (A | X, S, M)
# A = 0 if decile_score is less than 8
A = dat$compas
# A[A == 0] = -1   

fmla_A = as.formula(A ~ X1 + X2 + S + M + SX1 + SX2 + MS + MX1 + MX2)
inter_A_S = data.frame(X1, X2, M) 
inter_A_M = data.frame(X1, X2, S) 

# (Y) +++++++++++++++++++++++ 
# Y: utility (Y | S, X, M, A)

# interactions
SX1 = S*X1 
SX2 = S*X2
MX1 = M*X1 
MX2 = M*X2
AS = A*S 
AM = A*M 
AX1 = A*X1 
AX2 = A*X2
MS = M*S

## Utility
Y_recid = dat$recid
sink("compas_experiment8.txt")

obsA.vec <- fairA.vec <- unfairA.vec <- c()
obsY.vec <- fairY.vec <- unfairY.vec <- c()
obsBR.vec <- fairBR.vec <- unfairBR.vec <- c()
obsWR.vec <- fairWR.vec <- unfairWR.vec <- c()
theta.vec <- c(4.00,3.50,3.40,3.30,3.20,3.10,3.00,2.80,2.70,2.60,2.58,2.56,2.54,2.52,2.5,2.48,2.46,2.44,2.42,2.4,2.38,2.36,2.34,2.32,2.30,2.20,2.10,2.0,1.8,1.6,1.4,1.2,1.0)

for(i in theta.vec){
  
C1 = rep(-1.0,n) #-abs( rnorm(n,1.5,1) )
theta = i
cat("theta : ", theta, "\n")
C2 = rep(-1.0*theta,n) #-abs( rnorm(n,1.5*theta,1) )
C3 = rep(1.0,n) #abs( rnorm(n,1.5,1) )

Y = A*(C1) + (1-A)*(Y_recid*( C2 ) + (1-Y_recid)*( C3 ))

fmla_Y = as.formula(Y ~ X1 + X2 + S + M + A + SX1 + SX2 + AS + AM + MS + MX1 + MX2 + AX1 + AX2)
inter_Y_S = data.frame(A, M, X1, X2)  
inter_Y_A = data.frame(M, S, X1, X2)
inter_Y_M = data.frame(A, S, X1, X2) 

# data 
dat = data.frame(X1, X2, S, M, A, SX1, SX2, MS, MX1, MX2, AX1, AX2, AS, AM, MS, Y)

# ++++++++++++++
# regular coeff 
# ++++++++++++++
beta_S = mean(dat$S)
beta_M = glm(fmla_M, family=binomial, dat)$coefficients
beta_A = glm(fmla_A, family=binomial, dat)$coefficients
beta_Y = lm(fmla_Y, dat)$coefficients

cat("eff_SY (true) : ", (eff_SY = pse_SY_SM(dat, beta_S, beta_M, inter_M_S)) ,"\n")

cat("eff_SA (true) : ", (eff_SA = pse_SA_SM(dat, beta_S, beta_M, inter_M_S)) , "\n")


# ++++++++++++++
# unfair coeff 
# ++++++++++++++
var = 1
tau_cont_u = Inf
tau_cont_l = -Inf
tau_bin_u = Inf 
tau_bin_l = -Inf
funcSY = pse_SY_SM
funcSA = pse_SA_SM
beta_start = NULL

#result_uf = optimize_nloptr(dat, fmla_Y, fmla_M, funcSY, funcSA, inter_M_S, inter_Y_S, 
#                            var, tau_cont_u, tau_cont_l, tau_bin_u, tau_bin_l, beta_start) 

#result_uf = optimize_nloptr(dat, fmla_M, funcSY, funcSA, inter_M_S, 
#                            tau_cont_u, tau_cont_l, tau_bin_u, tau_bin_l, beta_start)
#beta_M_uf = result_uf$beta_M
#beta_S_uf = result_uf$beta_S

beta_M_uf = beta_M
beta_S_uf = beta_S
beta_Y_uf = beta_Y

#cat("eff_SY (uf) : ", (eff_SY_uf = pse_SY_SM(dat, beta_S_uf, beta_M_uf, inter_M_S)) ,"\n")

#cat("eff_SA (uf) : ", (eff_SA_uf = pse_SA_SM(dat, beta_S_uf, beta_M_uf, inter_M_S)) ,"\n")


# ++++++++++++++
# fair coeff
# ++++++++++++++
var = 1
tau_cont_u = 0.05
tau_cont_l = -0.05 
tau_bin_u = 1.05
tau_bin_l = 0.95
funcSY = pse_SY_SM 
funcSA = pse_SA_SM
beta_start = NULL

result_f = optimize_nloptr(dat, fmla_M, funcSY, funcSA, inter_M_S, 
                           tau_cont_u, tau_cont_l, tau_bin_u, tau_bin_l, beta_start)
beta_M_f = result_f$beta_M
beta_S_f = result_f$beta_S

cat("eff_SY (f) : ", (eff_SY_f = pse_SY_SM(dat, beta_S_f, beta_M_f, inter_M_S)) ,"\n")

cat("eff_SA (f) : ", (eff_SA_f = pse_SA_SM(dat, beta_S_f, beta_M_f, inter_M_S)) , "\n")


# +++++++++++++++++++++++
# Q-learning 
# +++++++++++++++++++++++

# (Unfair Policy)
Aopt_uf_q = Qlearning_unfair(dat, fmla_Y, beta_Y_uf, inter_Y_A)
Yopt_uf_q = Y_Qlearning_unfair(dat, fmla_Y, beta_Y_uf, inter_Y_A, Aopt_uf_q$Aopt)

# (Fair Policy)
Aopt_f_q = Qlearning(dat, fmla_M, fmla_A, fmla_Y, beta_S_f, beta_M_f, beta_A, beta_Y, inter_A_M, inter_Y_M)
#(dat, fmla_M, fmla_Y, beta_M_f, beta_Y, inter_Y_A, inter_Y_M)
Yopt_f_q = Y_Qlearning_unfair(dat, fmla_Y, beta_Y_uf, inter_Y_A, Aopt_f_q$Aopt)

# Comparisons
cat(" E(Y) given: \n",
    "   optimal unfair decisions:", mean(Yopt_uf_q), "\n",
    "   optimal fair decisions:", mean(Yopt_f_q), "\n",
    "   original decisions:", mean(dat$Y), "\n")
cat(" Tables (0,1) : \n",
"A: ", table(dat$A), "\n",
"Aopt fair: ", table(Aopt_f_q), "\n",
"Aopt unfair: ", table(Aopt_uf_q), "\n"
)

obsA <- table(dat$A)
fairA <- table(Aopt_f_q)
unfairA <- table(Aopt_uf_q)

### added 5/08
obsBR <- table(dat$A,dat$S)[2,2] / obsA[2]
obsWR <- table(dat$A,dat$S)[2,1] / obsA[2]
fairBR <- ifelse(all(Aopt_f_q[[1]]==0),0,table(Aopt_f_q[[1]],dat$S)[2,2] / fairA[2])
fairWR <- ifelse(all(Aopt_f_q[[1]]==0),0,table(Aopt_f_q[[1]],dat$S)[2,1] / fairA[2])
unfairBR <- ifelse(all(Aopt_uf_q[[1]]==0),0,table(Aopt_uf_q[[1]],dat$S)[2,2] / unfairA[2])
unfairWR <- ifelse(all(Aopt_uf_q[[1]]==0),0,table(Aopt_uf_q[[1]],dat$S)[2,1] / unfairA[2])
cat(" Breakdown by race : \n",
    "obsBR: ", obsBR, "\n",
    "fairBR: ", fairBR, "\n",
    "unfairBR: ", unfairBR, "\n",
    "fairWR: ", fairWR, "\n",
    "unfairWR: ", unfairWR, "\n"
)
obsBR.vec <- c(obsBR,obsBR.vec)
fairBR.vec <- c(fairBR,fairBR.vec)
unfairBR.vec <- c(unfairBR,unfairBR.vec)
obsWR.vec <- c(obsWR,obsWR.vec)
fairWR.vec <- c(fairWR,fairWR.vec)
unfairWR.vec <- c(unfairWR,unfairWR.vec)
###

obsA.vec <- c(obsA[2],obsA.vec)
fairA.vec <- c(fairA[2],fairA.vec)
unfairA.vec <- c(unfairA[2],unfairA.vec)

obsY.vec <- c(mean(dat$Y),obsY.vec)
fairY.vec <- c(mean(Yopt_f_q),fairY.vec)
unfairY.vec <- c(mean(Yopt_uf_q),unfairY.vec)

} ## end for loop
sink()

#theta_vec <- c(1,2,4,6,8,10,12,14)
#ratios_f <- c(0/5278, 0/5267, 0/5278, 3199/5278, 3199/5278, 4857/5278, 4857/5278, 4857/5278) ## |A==1|/|A==0| incarceration rate
#ratios_uf <- c(0/5278, 0/5267, 2724/5278, 4066/5278, 4218/5278, 4577/5278, 5033/5278, 5110/5278)
#ratios_orig <- c(1524/5278, 1524/5278, 1524/5278, 1524/5278, 1524/5278, 1524/5278, 1524/5278, 1524/5278)

fairA.vec[which(is.na(fairA.vec))] <- 0
unfairA.vec[which(is.na(unfairA.vec))] <- 0
ratios_f <- fairA.vec / rep(n,length(theta.vec))
ratios_uf <- unfairA.vec / rep(n,length(theta.vec))
ratios_orig <- obsA.vec / rep(n,length(theta.vec))

## overall incarceration rate

plot(sort(theta.vec),ratios_f,type = "l", lwd = 2, lty = "dashed", xaxt = "n", xlab = "Utility Parameter", ylab = "Incarceration Rate", cex.lab=1.4, cex.axis=1.4, cex.main=1.4, cex.sub=1.4)
#axis(1, at=1:length(theta.vec), labels=sort(theta.vec))
axis(1, at=1:4, labels=c(1,2,3,4))

lines(sort(theta.vec),ratios_uf, type = "l", lwd = 2)
lines(sort(theta.vec),ratios_orig, lty = "dotted", lwd = 2)
legend(locator(1), legend = c("Optimal unfair policy","Optimal Fair policy","Observed policy"), 
       lty = c(1, 2, 3), lwd = c(2,2,2), cex = 1.0, pt.cex = 2.0, y.intersp = 0.1, x.intersp = 0.25, bty="n", inset=c(.57,-.2))
#legend("bottomleft", legend = c("Optimal unfair policy","Optimal Fair policy","Observed policy"), 
#       lty = c(1, 2, 3), lwd = c(2,2,2), cex = 1.0, pt.cex = 2.0, y.intersp = 0.1, x.intersp = 0.25, bty="n", inset=c(.57,-.2))
#cex = 1.4 y.int = 0.3

plot(sort(theta.vec),fairY.vec,type="l",lwd=2,lty="dashed",xaxt = "n", xaxt = "n", xlab = "Utility Parameter", ylab = "Expected Utility", cex.lab=1.4, cex.axis=1.4, cex.main=1.4, cex.sub=1.4)
axis(1, at=1:4, labels=c(1,2,3,4))
lines(sort(theta.vec),unfairY.vec,type="l",lwd=2)
lines(sort(theta.vec),obsY.vec,type="l",lwd=2,lty="dotted")
legend("topleft", legend = c("Optimal unfair policy","Optimal Fair policy","Observed policy"), 
       lty = c(1, 2, 3), lwd = c(2,2,2), cex = 1.0, pt.cex = 2.0, y.intersp = 0.1, x.intersp = 0.25, bty="n", inset=c(.57,-.2))

plot(sort(theta.vec),unfairY.vec-fairY.vec,type="l",lwd=2,xaxt = "n")
axis(1, at=1:4, labels=c(1,2,3,4))


### new plots...

#main = "Group-level Incarceration Rates as a Function of Utility Parameter in the COMPAS Data"
plot(sort(theta.vec),unfairBR.vec,type="l",lwd=2,lty="dashed",col="blue",xaxt = "n", xaxt = "n",xlab = "Utility Parameter", ylab = "Incarceration Rate",cex.lab=1.4, cex.axis=1.4, cex.main=1.4, cex.sub=1.4)
axis(1, at=1:4, labels=c(1,2,3,4))
lines(sort(theta.vec),unfairWR.vec,type="l",lwd=2,lty="dashed",col="red")

lines(sort(theta.vec),fairBR.vec,type="l",lwd=2,col="blue")
lines(sort(theta.vec),fairWR.vec,type="l",lwd=2,col="red")

lines(sort(theta.vec),obsBR.vec,type="l",lwd=2,lty="dotted",col="blue")
lines(sort(theta.vec),1-obsBR.vec,type="l",lwd=2,lty="dotted",col="red")

legend(locator(1), legend = c("Optimal unrestricted policy","Optimal fair policy","Observed policy", "African-American","Caucasian"), 
       lty = c(2, 1, 3,NA,NA), lwd = c(2,2,2,NA,NA), pch=c(NA,NA,NA,22,22), pt.bg = c(NA,NA,NA,"blue","red"),cex = 1.0, pt.cex = 2.0, y.intersp = 0.1, x.intersp = 0.25, bty="n")#, inset=c(.57,-.7))

