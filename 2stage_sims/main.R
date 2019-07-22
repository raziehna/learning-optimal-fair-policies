

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

setwd("~/Desktop/Algorithmic Fairness/FairPolicy/github/2stage_sims/")

source("compute_effects.R")
source("optimization_nloptr.R")
source("policy_Qlearning.R")
source("policy_Qlearning_unfair.R")


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
S = rbinom(n, 1, 0.5)

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


# (W) +++++++++++++
# generate W | X, S, M, A
AX1 = A*X1
AX2 = A*X2
AX3 = A*X3
AS = A*S 
AM = A*M
pW = 1/(1 + exp(-2 + X1 + X2 + S + M + A + SX2 + MS + AS + AM)) 
W = rbinom(n, 1, pW)

fmla_W = as.formula(W ~ X1 + X2 + S + M + A + SX2 + MS + AS + AM)
inter_W_S = data.frame(X2, M, A)
inter_W_M = data.frame(A, S)
inter_W_A = data.frame(S, M)


# (B) +++++++++++++
# generate B | X, S, M, A, W
SW = W*S
AW = W*A
WX1 = W*X1
WX2 = W*X2
WX3 = W*X3
MW = M*W

pB = 1/(1 + exp(1 - X1 + X2 + S + M + A + W - SX1 + SX2 + MS - 3*MX1 + 0.5*MX2 - AS - AX1 - AX2))
B = rbinom(n, 1, pB)   

fmla_B = as.formula(B ~ X1 + X2 + S + M + A + W + SX1 + SX2 + MS + MX1 + MX2 + AS + AX1 + AX2)
inter_B_S = data.frame(X1, X2, M, A)
inter_B_M = data.frame(S, X1, X2)
inter_B_A = data.frame(X1, X2, S)


# Y +++++++++++++
# generate Y | X, S, M, A, W, B
BX1 = B*X1 
BX2 = B*X2
BX3 = B*X3
BS = B*S 
BM = B*M 
AB = B*A 
BW = B*W
epsY = rnorm(n, 0, 1)
Y = 2.5 + X1 + X2 + S + M + A + W + B + SX1 + SX2 + MS + 
  AS + AM - 2*AW + MW + SW - 3*BX1 + 2*BX2 - BM + WX1 + epsY 

fmla_Y = as.formula(Y ~ X1 + X2 + S + M + A + W + B + SX1 + SX2 + MS + 
                      AS + AM + AW + MW + SW + BX1 + BX2 + BM + WX1)
inter_Y_S = data.frame(X1, X2, M, A, W)   
inter_Y_M = data.frame(S, A, W, B)
inter_Y_A = data.frame(S, M, W) # in Y model what interacts with A? 
inter_Y_B = data.frame(X1, X2, M)
inter_Y_W = data.frame(A, M, S, X1)


# data +++++++++++++
dat = data.frame(X1, X2, X3, S, M, A, W, B, SX1, SX2, SX3, MS, MX1, MX2, MX3, AS, AM, AX1, 
                 AX2, AX3, AW, MW, SW, WX1, WX2, WX3, BS, AB, BM, BW, BX1, BX2, BX3, Y)


# +++++++++++++++++++++++++++++++++
# regular coeff 
# +++++++++++++++++++++++++++++++++
beta_S = mean(dat$S)
beta_M = glm(fmla_M, family=binomial, dat)$coefficients
beta_A = glm(fmla_A, family=binomial, dat)$coefficients
beta_W = lm(fmla_W, dat)$coefficients
beta_B = glm(fmla_B, family=binomial, dat)$coefficients
beta_Y = lm(fmla_Y, dat)$coefficients

# (eff_SY = pse_SY_SMY(dat, beta_Y, beta_M, inter_M_S, inter_Y_S)) 
(eff_SY = pse_SY_SM(dat, beta_S, beta_M, inter_M_S)) 
(eff_SB = pse_SB_SM(dat, beta_S, beta_M, inter_M_S))
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
funcSB = pse_SB_SM
funcSA = pse_SA_SM
beta_start = NULL

result_f = optimize_nloptr(dat, fmla_M, funcSY, funcSB, funcSA, inter_M_S, tau_cont_u, tau_cont_l, tau_bin_u, tau_bin_l, beta_start)

beta_S_f = result_f$beta_S 
beta_M_f = result_f$beta_M

(eff_SY_f = pse_SY_SM(dat, beta_S_f, beta_M_f, inter_M_S))
(eff_SB_f = pse_SB_SM(dat, beta_S_f, beta_M_f, inter_M_S))
(eff_SA_f = pse_SA_SM(dat, beta_S_f, beta_M_f, inter_M_S))


# +++++++++++++++++++++++++++++++++
# Q-learning
# +++++++++++++++++++++++++++++++++

# (Unfair Policy)
policy_uf_q = Qlearning_unfair(dat, fmla_Y, fmla_W, beta_Y_uf, beta_W,
                               inter_Y_B, inter_Y_W, inter_Y_A, inter_W_A)
Aopt_uf_q = policy_uf_q$Aopt
Bopt_uf_q = policy_uf_q$Bopt
Yopt_uf_q = Y_Qlearning_unfair(dat, fmla_Y, beta_Y_uf, inter_Y_B, inter_Y_A, Aopt_uf_q, Bopt_uf_q)

# (Fair Policy)
policy_f_q = Qlearning(dat, fmla_M, fmla_A, fmla_W, fmla_B, fmla_Y, beta_S_f, beta_M_f, beta_A, beta_W, beta_Y,
                       inter_A_M, inter_W_M, inter_W_A, inter_Y_B, inter_Y_M)
Aopt_f_q = policy_f_q$Aopt # optimal_fair_qlearning
Bopt_f_q = policy_f_q$Bopt
Yopt_f_q = Y_Qlearning_unfair(dat, fmla_Y, beta_Y_uf, inter_Y_B, inter_Y_A, Aopt_f_q, Bopt_f_q)


# +++++++++++++++++++++++++++++++++ 
# Comparisons 
# +++++++++++++++++++++++++++++++++

cat(" Q-learning - E(Y) given : \n",
    "   optimal unfair decisions:", mean(Yopt_uf_q), "\n",
    "   optimal fair decisions:", mean(Yopt_f_q), "\n \n")


