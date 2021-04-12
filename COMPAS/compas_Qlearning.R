
# +++++++++++++++++++++++++++++++++++++++
# Q-learning
# +++++++++++++++++++++++++++++++++++++++

Qlearning <- function(dat, fmla_M, fmla_A, fmla_Y, beta_S, beta_M, beta_A, beta_Y, inter_A_M, inter_Y_M){
  
  Xm = as.data.frame(model.matrix(fmla_M, data=model.frame(dat)))
  Xa = as.data.frame(model.matrix(fmla_A, data=model.frame(dat)))
  Xy = as.data.frame(model.matrix(fmla_Y, data=model.frame(dat)))
  A = dat$A
  
  # S 
  pS = beta_S
  
  # M, S = 1 
  Xm_S1 = edit_dat_S(Xm, inter_M_S, s = 1) 
  pM1_S1 = 1/(1+exp(-as.matrix(Xm_S1)%*%beta_M))
  # pM1_S1[pM1_S1 < 0.005] = 0.005
  # pM1_S1[pM1_S1 > 0.995] = 0.995 
  
  # M, S = 0 
  Xm_S0 = edit_dat_S(Xm, inter_M_S, s = 0)
  pM1_S0 = 1/(1+exp(-as.matrix(Xm_S0)%*%beta_M))
  # pM1_S0[pM1_S0 < 0.005] = 0.005 
  # pM1_S0[pM1_S0 > 0.995] = 0.995 

  
  # (K = 2) ++++++++++++++++++++
  
  # A, M = 1, S = 1
  Xa_M1S1 = edit_dat_MS(Xa, inter_A_M, inter_A_S, m = 1, s = 1)
  pA1_M1S1 = 1/(1+exp(-as.matrix(Xa_M1S1)%*%beta_A))
  # pA1_M1S1[pA1_M1S1 < 0.005] = 0.005
  # pA1_M1S1[pA1_M1S1 > 0.995] = 0.995 
  
  # A, M = 0, S = 1
  Xa_M0S1 = edit_dat_MS(Xa, inter_A_M, inter_A_S, m = 0, s = 1)
  pA1_M0S1 = 1/(1+exp(-as.matrix(Xa_M0S1)%*%beta_A))
  # pA1_M0S1[pA1_M0S1 < 0.005] = 0.005
  # pA1_M0S1[pA1_M0S1 > 0.995] = 0.995
  
  # A, M = 1, S = 0
  Xa_M1S0 = edit_dat_MS(Xa, inter_A_M, inter_A_S, m = 1, s = 0)
  pA1_M1S0 = 1/(1+exp(-as.matrix(Xa_M1S0)%*%beta_A))
  # pA1_M1S0[pA1_M1S0 < 0.005] = 0.005
  # pA1_M1S0[pA1_M1S0 > 0.995] = 0.995 
  
  # A, M = 0, S = 0
  Xa_M0S0 = edit_dat_MS(Xa, inter_A_M, inter_A_S, m = 0, s = 0)
  pA1_M0S0 = 1/(1+exp(-as.matrix(Xa_M0S0)%*%beta_A))
  # pA1_M0S0[pA1_M0S0 < 0.005] = 0.005
  # pA1_M0S0[pA1_M0S0 > 0.995] = 0.995
  
  # p(ai | m = 1, s = 1, ...)
  pA_M1S1 = pA1_M1S1 
  pA_M1S1[A == 0] = 1 - pA1_M1S1[A == 0]  
  # p(ai | m = 0, s = 1, ...)
  pA_M0S1 = pA1_M0S1 
  pA_M0S1[A == 0] = 1 - pA1_M0S1[A == 0]
  # p(ai | m = 1, s = 0, ...)
  pA_M1S0 = pA1_M1S0
  pA_M1S0[A == 0] = 1 - pA1_M1S0[A == 0]  
  # p(ai | m = 0, s = 0, ...)
  pA_M0S0 = pA1_M0S0
  pA_M0S0[A == 0] = 1 - pA1_M0S0[A == 0]
  
  # \sum_{m, s} p(S, M, A | X) 
  denom = pA_M1S1 * pM1_S1 * pS + 
          pA_M0S1 * ( 1- pM1_S1) * pS +  
          pA_M1S0 * pM1_S0 * (1 - pS) + 
          pA_M0S0 * ( 1- pM1_S0) * (1 - pS)  
  
  # E[Y | ., A = 1] 
  Xy_A1M1S1 = edit_dat_AMS(Xy, inter_Y_A, inter_Y_M, inter_Y_S, a = 1, m = 1, s = 1) 
  Xy_A1M0S1 = edit_dat_AMS(Xy, inter_Y_A, inter_Y_M, inter_Y_S, a = 1, m = 0, s = 1) 
  Xy_A1M1S0 = edit_dat_AMS(Xy, inter_Y_A, inter_Y_M, inter_Y_S, a = 1, m = 1, s = 0) 
  Xy_A1M0S0 = edit_dat_AMS(Xy, inter_Y_A, inter_Y_M, inter_Y_S, a = 1, m = 0, s = 0)
  
  Y_A1M1S1 = as.matrix(Xy_A1M1S1)%*%beta_Y
  Y_A1M0S1 = as.matrix(Xy_A1M0S1)%*%beta_Y
  Y_A1M1S0 = as.matrix(Xy_A1M1S0)%*%beta_Y
  Y_A1M0S0 = as.matrix(Xy_A1M0S0)%*%beta_Y
  Y_A1 = (1/denom)*{ Y_A1M1S1 * pA_M1S1 * pM1_S1 * pS + 
                     Y_A1M0S1 * pA_M0S1 * ( 1 - pM1_S1) * pS + 
                     Y_A1M1S0 * pA_M1S0 * pM1_S0 * (1 - pS) + 
                     Y_A1M0S0 * pA_M0S0 * ( 1 - pM1_S0) * (1 - pS) }

  # E[Y | ., A = 0]
  Xy_A0M1S1 = edit_dat_AMS(Xy, inter_Y_A, inter_Y_M, inter_Y_S, a = 0, m = 1, s = 1) 
  Xy_A0M0S1 = edit_dat_AMS(Xy, inter_Y_A, inter_Y_M, inter_Y_S, a = 0, m = 0, s = 1) 
  Xy_A0M1S0 = edit_dat_AMS(Xy, inter_Y_A, inter_Y_M, inter_Y_S, a = 0, m = 1, s = 0) 
  Xy_A0M0S0 = edit_dat_AMS(Xy, inter_Y_A, inter_Y_M, inter_Y_S, a = 0, m = 0, s = 0)
  
  Y_A0M1S1 = as.matrix(Xy_A0M1S1)%*%beta_Y
  Y_A0M0S1 = as.matrix(Xy_A0M0S1)%*%beta_Y
  Y_A0M1S0 = as.matrix(Xy_A0M1S0)%*%beta_Y
  Y_A0M0S0 = as.matrix(Xy_A0M0S0)%*%beta_Y
  Y_A0 = (1/denom)*{ Y_A0M1S1 * pA_M1S1 * pM1_S1 * pS + 
                     Y_A0M0S1 * pA_M0S1 * ( 1 - pM1_S1) * pS + 
                     Y_A0M1S0 * pA_M1S0 * pM1_S0 * (1 - pS) + 
                     Y_A0M0S0 * pA_M0S0 * ( 1 - pM1_S0) * (1 - pS) }

  # Optimal B
  Aopt = (Y_A1 > Y_A0) + 0 
  table(Aopt) 
  
  return(list(Aopt=Aopt)) 
}


