

# Notes:
# Make sure you name interactions in alphabetical order: e.g. AS not SA 

# Graph: X S M A W B Y 
# baseline: X, S 
# mediator: M 
# decision 1: A 
# intermediate: W 
# decision 2: B 
# final outcome: Y  

rm(list = ls(all=TRUE)) 
cat("\f") 

n = 5000 # sample size
iter = 100 # bootstrap iterations

# +++++++++++++++++++++++++++++++++
# Load libraries 
# +++++++++++++++++++++++++++++++++
library(mvtnorm) 
library(nloptr)
library(gtools)
library(rpart) 
library(rpart.plot)

set.seed(0)

setwd("~/Desktop/Algorithmic Fairness/FairPolicy/github/1stage_sims/")

source("compute_effects.R")
source("optimization_nloptr.R")
source("policy_Qlearning.R")
source("policy_Qlearning_unfair.R")
source("policy_Vsearch.R")
source("policy_Vsearch_unfair.R")
source("policy_Gestimation.R")
source("policy_Gestimation_unfair.R")
source("bootstrap.R")


# +++++++++++++++++++++++++++++++++
# Generate data
# +++++++++++++++++++++++++++++++++

# (X) +++++++++++++
# generate X
X = rmvnorm(n, mean = rep(0, 3), sigma = diag(3))

X1 = abs(X[, 1])
X2 = X[, 2]
X3 = X[, 3]

# (S) +++++++++++++
# generate S 
S = rbinom(n, 1, 0.3)

# (M) +++++++++++++
# genereta M | X, S
SX1 = S*X1
SX2 = S*X2
SX3 = S*X3
pM = 1/(1 + exp(-1 + X1 + X2 + X3 + S + 3*SX1 + SX2 + SX3)) 
M = rbinom(n, 1, pM)

fmla_M = as.formula(M ~ X1 + X2 + X3 + S + SX1 + SX2 + SX3)
inter_M_S = data.frame(X1, X2, X3) # in M model what interacts with S? 

# (A) +++++++++++++
# generate A | X, S, M 
MS = M*S
MX1 = M*X1 
MX2 = M*X2
MX3 = M*X3
pA = 1/(1 + exp(1 - X1 + X2 + S + M - SX1 + SX2 + MS - 3*MX1 + 0.5*MX2 ))
A = rbinom(n, 1, pA) 

fmla_A = as.formula(A ~ X1 + X2 + S + M + SX1 + SX2 + MS + MX1 + MX2)
inter_A_S = data.frame(X1, X2, M) # in A model what interacts with S? 
inter_A_M = data.frame(X1, X2, S) # in A model what interacts with M?

# (Y) +++++++++++++
# generate Y | X, S, M, A
AX1 = A*X1
AX2 = A*X2
AX3 = A*X3
AS = A*S 
AM = A*M
epsY = rnorm(n, 0, 1)
Y = -2 + X1 + X2 + S + M + A -3*SX2 + MS + AS + AM + AX2 + AX3 + epsY 

fmla_Y = as.formula(Y ~ X1 + X2 + S + M + A + SX2 + MS + AS + AM + AX2 + AX3)
inter_Y_S = data.frame(X2, M, A)
inter_Y_M = data.frame(A, S)
inter_Y_A = data.frame(S, M, X2, X3)

# data +++++++++++++
dat = data.frame(X1, X2, X3, S, M, A, SX1, SX2, SX3, MS, MX1, MX2, MX3, AS, AM, AX1, AX2, AX3, Y)

Y_vec = c()
for (i in 1:iter){
 cat(i, "\n")
 idx = sample.int(n, n, replace = TRUE)
 Y_vec = c(Y_vec, mean(dat[idx, ]$Y))
}


# +++++++++++++++++++++++++++++++++
# regular coeff 
# +++++++++++++++++++++++++++++++++
beta_S = mean(dat$S)
beta_M = glm(fmla_M, family=binomial, dat)$coefficients
beta_A = glm(fmla_A, family=binomial, dat)$coefficients
beta_Y = lm(fmla_Y, dat)$coefficients

# (eff_SY = pse_SY_SMY(dat, beta_Y, beta_M, inter_M_S, inter_Y_S)) 
(eff_SY = pse_SY_SM(dat, beta_S, beta_M, inter_M_S)) 
(eff_SA = pse_SA_SM(dat, beta_S, beta_M, inter_M_S)) 


# +++++++++++++++++++++++++++++++++
# unfair coeff
# +++++++++++++++++++++++++++++++++
beta_M_uf = beta_M 
beta_Y_uf = beta_Y


# +++++++++++++++++++++++++++++++++
# fair coeff
# +++++++++++++++++++++++++++++++++
tau_cont_u = 0.1
tau_cont_l = -0.1 
tau_bin_u = 1.05
tau_bin_l = 0.95
funcSY = pse_SY_SM
funcSA = pse_SA_SM
beta_start = NULL

result_f = optimize_nloptr(dat, fmla_M, funcSY, funcSA, inter_M_S, 
                           tau_cont_u, tau_cont_l, tau_bin_u, tau_bin_l, beta_start)
beta_S_f = result_f$beta_S 
beta_M_f = result_f$beta_M

(eff_SY_f = pse_SY_SM(dat, beta_S_f, beta_M_f, inter_M_S))
(eff_SA_f = pse_SA_SM(dat, beta_S_f, beta_M_f, inter_M_S))

Yopt_vec_uf_q = c()
Yopt_vec_f_q = c()
Yopt_vec_uf_v = c()
Yopt_vec_f_v = c()
Yopt_vec_uf_g = c()
Yopt_vec_f_g = c()

for (i in 1:iter){
 cat(i, "\n")
 idx = sample.int(n, n, replace = TRUE)
 
# +++++++++++++++++++++++++++++++++
# Q-learning   
# +++++++++++++++++++++++++++++++++

# (Unfair Policy)
policy_uf_q = Qlearning_unfair(dat[idx, ], fmla_Y, beta_Y_uf, inter_Y_A[idx, ])
Aopt_uf_q = policy_uf_q$Aopt
Yopt_uf_q = Y_Qlearning_unfair(dat[idx, ], fmla_Y, beta_Y_uf, inter_Y_A[idx, ], Aopt_uf_q)
Yopt_vec_uf_q = c(Yopt_vec_uf_q, mean(Yopt_uf_q))

# boot_uf_q = bootstrap_policy(iter, method="Q-learning", policy="Unfair")
# Aopt_mat_uf_q = boot_uf_q$Aopt
# Bopt_mat_uf_q = boot_uf_q$Bopt
# Yhat_vec_uf_q = boot_uf_q$Yhat_vec
# Yopt_vec_uf_q = boot_uf_q$Yopt_vec

# (Fair Policy)
policy_f_q = Qlearning(dat[idx, ], fmla_M, fmla_A, fmla_Y, beta_S_f, beta_M_f, beta_A, beta_Y, inter_A_M, inter_Y_M)
Aopt_f_q = policy_f_q$Aopt # optimal_fair_qlearning
Yopt_f_q = Y_Qlearning_unfair(dat[idx, ], fmla_Y, beta_Y_uf, inter_Y_A[idx, ], Aopt_f_q)
Yopt_vec_f_q = c(Yopt_vec_f_q, mean(Yopt_f_q))

# boot_f_q = bootstrap_policy(iter, method="Q-learning", policy="Fair")
# Aopt_mat_f_q = boot_f_q$Aopt
# Bopt_mat_f_q = boot_f_q$Bopt
# Yhat_vec_f_q = boot_f_q$Yhat_vec
# Yopt_vec_f_q = boot_f_q$Yopt_vec

# +++++++++++++++++++++++++++++++++
# Value-search
# +++++++++++++++++++++++++++++++++

# Define policies 
policy <- function(alpha, dat){
 K = length(alpha)
 fw_A = list() 
 fw_B = list()
 count = 0
 for (i in 1:K){
  for (j in 1:K){
   for (k in 1:K){
    for (m in 1:K){
     count = count + 1
     # A
     pA = 1/(1 + exp(1 + alpha[m]*X1 + alpha[k]*X2 + alpha[k]*S + alpha[k]*M + alpha[m]*S*X1 + 
                     alpha[k]*S*X2 + alpha[k]*M*S + alpha[j]*M*X1 + alpha[i]*M*X2)) 
     A = rbinom(n, 1, pA) 
     fw_A[[count]] = A
    }
   }
  }
 }
 return(list(fw_A = fw_A))
}

alpha = seq(-3, 3, by=1) 
fw = policy(alpha, dat[idx, ])

# (Unfair Policy)
result_uf_v = Vsearch_unfair(dat[idx, ], fmla_A, fmla_Y, beta_A, beta_Y, fw)
Aopt_uf_v = result_uf_v$Aopt
Policy_uf_v = result_uf_v$Value_max
Yopt_uf_v = Y_Vsearch_unfair(dat[idx, ], fmla_A, fmla_Y, beta_A, beta_Y, Aopt_uf_v)
Yopt_vec_uf_v = c(Yopt_vec_uf_v, mean(Yopt_uf_v))

# boot_uf_v = bootstrap_policy(iter, method="V-search", policy="Unfair")
# Aopt_mat_uf_v = boot_uf_v$Aopt
# Bopt_mat_uf_v = boot_uf_v$Bopt
# Yhat_vec_uf_v = boot_uf_v$Yhat_vec
# Yopt_vec_uf_v = boot_uf_v$Yopt_vec


# (Fair Policy)
result_f_v = Vsearch(dat[idx, ], fmla_M, fmla_A, fmla_Y, beta_S_f, beta_M_f, beta_A, beta_Y,
                      inter_M_S[idx, ], inter_A_M[idx, ], inter_A_S[idx, ], fw)
Aopt_f_v = result_f_v$Aopt
Policy_f_v = result_f_v$Value_max
Yopt_f_v = Y_Vsearch_unfair(dat[idx, ], fmla_A, fmla_Y, beta_A, beta_Y, Aopt_f_v)
Yopt_vec_f_v = c(Yopt_vec_f_v, mean(Yopt_f_v))

# boot_f_v = bootstrap_policy(iter, method="V-search", policy="Fair")
# Aopt_mat_f_v = boot_f_v$Aopt
# Bopt_mat_f_v = boot_f_v$Bopt
# Yhat_vec_f_v = boot_f_v$Yhat_vec
# Yopt_vec_f_v = boot_f_v$Yopt_vec

# +++++++++++++++++++++++++++++++++
# G-estimation
# +++++++++++++++++++++++++++++++++

k = 100 # number of draws from p*(m | . )
W = as.matrix(cbind( dat[idx, 1:5] ))  # remember to change it in Y_Gestimation

# (Unfair Policy) 
result_uf_g = Gestimation_unfair(dat[idx, ], fmla_A, fmla_Y, beta_A, beta_Y, inter_Y_A[idx, ], W)
Aopt_uf_g = result_uf_g$Aopt
psi1_uf = result_uf_g$psi_1
psi2_uf = result_uf_g$psi_2
Yopt_uf_g = Y_Gestimation_unfair(dat[idx, ], fmla_Y, beta_Y, inter_Y_A[idx, ], Aopt_uf_g, psi1_uf, psi2_uf, W)   
Yopt_vec_uf_g = c(Yopt_vec_uf_g, mean(Yopt_uf_g))

# boot_g_uf = bootstrap_policy(iter_boots, n, dat, "G-estimation", "Unfair")
# Yopt_uf_g_samples = boot_g_uf$Yopt

# (Fair Policy) 
# result_f_g = Gestimation(dat, fmla_A, fmla_M, fmla_Y, beta_S, beta_A, beta_M, beta_Y,
#             inter_Y_A, inter_Y_M, inter_Y_S, inter_A_M, inter_A_S, k)
# Aopt_f_g = result_f_g$Aopt
# psi1_f = result_f_g$psi_1 
# psi2_f = result_f_g$psi_2 
# Yopt_f_g = Y_Gestimation_unfair(dat, fmla_Y, beta_Y, inter_Y_A, Aopt_f_g, psi1_uf, psi2_uf, W) 

Aopt_f_g = Gestimation(dat[idx, ], fmla_M, beta_S_f, beta_M_f, inter_M_S[idx, ], psi1_uf, psi2_uf, W)
Yopt_f_g = Y_Gestimation_unfair(dat[idx, ], fmla_Y, beta_Y, inter_Y_A[idx, ], Aopt_f_g, psi1_uf, psi2_uf, W)
Yopt_vec_f_g = c(Yopt_vec_f_g, mean(Yopt_f_g))
}

# +++++++++++++++++++++++++++++++++
# Save resutls 
# +++++++++++++++++++++++++++++++++

# write.csv(Aopt_mat_f_q, "Aopt_mat_f_q.csv", row.names = F)
# write.csv(Aopt_mat_uf_q, "Aopt_mat_uf_q.csv", row.names = F)
# write.csv(Aopt_mat_f_v, "Aopt_mat_f_v.csv", row.names = F)
# write.csv(Aopt_mat_uf_v, "Aopt_mat_uf_v.csv", row.names = F)
# 
# write.csv(Bopt_mat_f_q, "Bopt_mat_f_q.csv", row.names = F)
# write.csv(Bopt_mat_uf_q, "Bopt_mat_uf_q.csv", row.names = F)
# write.csv(Bopt_mat_f_v, "Bopt_mat_f_v.csv", row.names = F)
# write.csv(Bopt_mat_uf_v, "Bopt_mat_uf_v.csv", row.names = F)
# 
write.csv(Yopt_vec_f_q, "Yopt_vec_f_q.csv", row.names = F)
write.csv(Yopt_vec_uf_q, "Yopt_vec_uf_q.csv", row.names = F)
write.csv(Yopt_vec_f_v, "Yopt_vec_f_v.csv", row.names = F)
write.csv(Yopt_vec_uf_v, "Yopt_vec_uf_v.csv", row.names = F)
write.csv(Yopt_vec_f_g, "Yopt_vec_f_g.csv", row.names = F)
write.csv(Yopt_vec_uf_g, "Yopt_vec_uf_g.csv", row.names = F)


# +++++++++++++++++++++++++++++++++ 
# Comparisons 
# +++++++++++++++++++++++++++++++++

cat(" Q-learning - E(Y) given : \n",
    "   optimal unfair decisions:", mean(Yopt_uf_q), "\n", 
    "   optimal fair decisions:", mean(Yopt_f_q), "\n", 
    "   original decisions:", mean(Y), "\n \n", 
    # 
    " Value-search - E(Y) given: \n", 
    "   optimal unfair decisions:", mean(Yopt_uf_v), "\n",
    "   optimal fair decisions:", mean(Yopt_f_v), "\n", 
    "   original decisions:", mean(Y), "\n \n", 
    # 
    " G-estimation - E(Y) given: \n", 
    "   optimal unfair decisions:", mean(Yopt_uf_g), "\n",
    "   optimal fair decisions:", mean(Yopt_f_g), "\n", 
    "   original decisions:", mean(Y), "\n \n")  


# # Maybe remove outliers  
# cat(" Q-learning - E(Y) given : \n",
#     "   optimal unfair decisions:", Yopt_vec_uf_q[1], "+-", 1.96*sd(Yopt_vec_uf_q[-1])/sqrt(iter), "\n",
#     "   optimal fair decisions:", Yopt_vec_f_q[1], "+-", 1.96*sd(Yopt_vec_f_q[-1])/sqrt(iter), "\n",
#     "   original decisions:", Yhat_vec_uf_q[1], "+-", 1.96*sd(Yhat_vec_uf_q[-1])/sqrt(iter), "\n \n",
#     # 
#     " Value-search - E(Y) given: \n",
#     "   optimal unfair decisions:", Yopt_vec_uf_v[1], "+-", 1.96*sd(Yopt_vec_uf_v[-1])/sqrt(iter), "\n",
#     "   optimal fair decisions:", Yopt_vec_f_v[1], "+-", 1.96*sd(Yopt_vec_f_v[-1])/sqrt(iter), "\n",
#     "   original decisions:", Yhat_vec_uf_v[1], "+-", 1.96*sd(Yhat_vec_uf_v[-1])/sqrt(iter), "\n") 

idx = 100

Yopt_vec_f_g[idx] 
1.96*sd(Yopt_vec_f_g[-idx])/sqrt(iter)

Yopt_vec_uf_g[idx] 
1.96*sd(Yopt_vec_uf_g[-idx])/sqrt(iter)

Yopt_vec_f_q[idx] 
1.96*sd(Yopt_vec_f_q[-idx])/sqrt(iter)

Yopt_vec_uf_q[idx] 
1.96*sd(Yopt_vec_uf_q[-idx])/sqrt(iter)

Yopt_vec_f_v[idx] 
1.96*sd(Yopt_vec_f_v[-idx])/sqrt(iter)

Yopt_vec_uf_v[idx] 
1.96*sd(Yopt_vec_uf_v[-idx])/sqrt(iter)

