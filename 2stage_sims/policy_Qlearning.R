
# +++++++++++++++++++++++++++++++++++++++
# Q-learning
# +++++++++++++++++++++++++++++++++++++++

Qlearning <- function(dat, fmla_M, fmla_A, fmla_W, fmla_B, fmla_Y, beta_S, beta_M, beta_A, beta_W, beta_Y, 
                      inter_A_M, inter_W_M, inter_W_A, inter_Y_B, inter_Y_M){
  
  Xm = as.data.frame(model.matrix(fmla_M, data=model.frame(dat)))
  Xa = as.data.frame(model.matrix(fmla_A, data=model.frame(dat)))
  Xw = as.data.frame(model.matrix(fmla_W, data=model.frame(dat)))
  Xb = as.data.frame(model.matrix(fmla_B, data=model.frame(dat)))
  Xy = as.data.frame(model.matrix(fmla_Y, data=model.frame(dat)))
  A = dat$A
  W = dat$W
  
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
  
  
  # W, M = 1, S = 1
  Xw_M1S1 = edit_dat_MS(Xw, inter_W_M, inter_W_S, m = 1, s = 1)
  pW1_M1S1 = 1/(1+exp(-as.matrix(Xw_M1S1)%*%beta_W))
  # pW1_M1S1[pW1_M1S1 < 0.005] = 0.005
  # pW1_M1S1[pW1_M1S1 > 0.995] = 0.995 
  
  # W, M = 0, S = 1
  Xw_M0S1 = edit_dat_MS(Xw, inter_W_M, inter_W_S, m = 0, s = 1)
  pW1_M0S1 = 1/(1+exp(-as.matrix(Xw_M0S1)%*%beta_W))
  # pW1_M0S1[pW1_M0S1 < 0.005] = 0.005
  # pW1_M0S1[pW1_M0S1 > 0.995] = 0.995 
  
  # W, M = 1, S = 0
  Xw_M1S0 = edit_dat_MS(Xw, inter_W_M, inter_W_S, m = 1, s = 0)
  pW1_M1S0 = 1/(1+exp(-as.matrix(Xw_M1S0)%*%beta_W))
  # pW1_M1S0[pW1_M1S0 < 0.005] = 0.005
  # pW1_M1S0[pW1_M1S0 > 0.995] = 0.995 
  
  # W, M = 0, S = 0
  Xw_M0S0 = edit_dat_MS(Xw, inter_W_M, inter_W_S, m = 0, s = 0)
  pW1_M0S0 = 1/(1+exp(-as.matrix(Xw_M0S0)%*%beta_W))
  # pW1_M0S0[pW1_M0S0 < 0.005] = 0.005
  # pW1_M0S0[pW1_M0S0 > 0.995] = 0.995 
  
  # p(wi | m = 1, s = 1, ...)
  pW_M1S1 = pW1_M1S1 
  pW_M1S1[W == 0] = 1 - pW1_M1S1[W == 0]
  # p(wi | m = 0, s = 1, ...)
  pW_M0S1 = pW1_M0S1 
  pW_M0S1[W == 0] = 1 - pW1_M0S1[W == 0]
  # p(wi | m = 1, s = 0, ...)
  pW_M1S0 = pW1_M1S0
  pW_M1S0[W == 0] = 1 - pW1_M1S0[W == 0]
  # p(wi | m = 0, s = 0, ...)
  pW_M0S0 = pW1_M0S0 
  pW_M0S0[W == 0] = 1 - pW1_M0S0[W == 0]
  
  
  # B 
  # B, M = 1, S = 1
  Xb_M1S1 = edit_dat_MS(Xb, inter_B_M, inter_B_S, m = 1, s = 1)
  pB1_M1S1 = 1/(1+exp(-as.matrix(Xb_M1S1)%*%beta_B))
  # pB1_M1S1[pB1_M1S1 < 0.005] = 0.005
  # pB1_M1S1[pB1_M1S1 > 0.995] = 0.995 
  
  # B, M = 0, S = 1
  Xb_M0S1 = edit_dat_MS(Xb, inter_B_M, inter_B_S, m = 0, s = 1)
  pB1_M0S1 = 1/(1+exp(-as.matrix(Xb_M0S1)%*%beta_B))
  # pB1_M0S1[pB1_M0S1 < 0.005] = 0.005
  # pB1_M0S1[pB1_M0S1 > 0.995] = 0.995 
  
  # B, M = 1, S = 0
  Xb_M1S0 = edit_dat_MS(Xb, inter_B_M, inter_B_S, m = 1, s = 0)
  pB1_M1S0 = 1/(1+exp(-as.matrix(Xb_M1S0)%*%beta_B))
  # pB1_M1S0[pB1_M1S0 < 0.005] = 0.005
  # pB1_M1S0[pB1_M1S0 > 0.995] = 0.995 
  
  # B, M = 0, S = 0
  Xb_M0S0 = edit_dat_MS(Xb, inter_B_M, inter_B_S, m = 0, s = 0)
  pB1_M0S0 = 1/(1+exp(-as.matrix(Xb_M0S0)%*%beta_B))
  # pB1_M0S0[pB1_M0S0 < 0.005] = 0.005
  # pB1_M0S0[pB1_M0S0 > 0.995] = 0.995 
  
  # p(bi | m = 1, s = 1, ...)
  pB_M1S1 = pB1_M1S1 
  pB_M1S1[B == 0] = 1 - pB1_M1S1[B == 0]
  # p(bi | m = 0, s = 1, ...)
  pB_M0S1 = pB1_M0S1 
  pB_M0S1[B == 0] = 1 - pB1_M0S1[B == 0]
  # p(bi | m = 1, s = 0, ...)
  pB_M1S0 = pB1_M1S0
  pB_M1S0[B == 0] = 1 - pB1_M1S0[B == 0]
  # p(bi | m = 0, s = 0, ...)
  pB_M0S0 = pB1_M0S0 
  pB_M0S0[B == 0] = 1 - pB1_M0S0[B == 0]
  
  # \sum_{m, s} p(S, M, A, W, B | X) 
  denom = pB_M1S1 * pW_M1S1 * pA_M1S1 * pM1_S1 * pS + 
          pB_M0S1 * pW_M0S1 * pA_M0S1 * ( 1- pM1_S1) * pS +  
          pB_M1S0 * pW_M1S0 * pA_M1S0 * pM1_S0 * (1 - pS) + 
          pB_M0S0 * pW_M0S0 * pA_M0S0 * ( 1- pM1_S0) * (1 - pS)  
  
  # E[Y | ., B = 1] 
  Xy_B1M1S1 = edit_dat_BMS(Xy, inter_Y_B, inter_Y_M, inter_Y_S, b = 1, m = 1, s = 1) 
  Xy_B1M0S1 = edit_dat_BMS(Xy, inter_Y_B, inter_Y_M, inter_Y_S, b = 1, m = 0, s = 1) 
  Xy_B1M1S0 = edit_dat_BMS(Xy, inter_Y_B, inter_Y_M, inter_Y_S, b = 1, m = 1, s = 0) 
  Xy_B1M0S0 = edit_dat_BMS(Xy, inter_Y_B, inter_Y_M, inter_Y_S, b = 1, m = 0, s = 0)
  
  Y_B1M1S1 = as.matrix(Xy_B1M1S1)%*%beta_Y
  Y_B1M0S1 = as.matrix(Xy_B1M0S1)%*%beta_Y
  Y_B1M1S0 = as.matrix(Xy_B1M1S0)%*%beta_Y
  Y_B1M0S0 = as.matrix(Xy_B1M0S0)%*%beta_Y
  Y_B1 = (1/denom)*{ Y_B1M1S1 * pB_M1S1 * pW_M1S1 * pA_M1S1 * pM1_S1 * pS + 
                     Y_B1M0S1 * pB_M0S1 * pW_M0S1 * pA_M0S1 * ( 1- pM1_S1) * pS + 
                     Y_B1M1S0 * pB_M1S0 * pW_M1S0 * pA_M1S0 * pM1_S0 * (1 - pS) + 
                     Y_B1M0S0 * pB_M0S0 * pW_M0S0 * pA_M0S0 * ( 1- pM1_S0) * (1 - pS) }

  # E[Y | ., B = 0]
  Xy_B0M1S1 = edit_dat_BMS(Xy, inter_Y_B, inter_Y_M, inter_Y_S, b = 0, m = 1, s = 1) 
  Xy_B0M0S1 = edit_dat_BMS(Xy, inter_Y_B, inter_Y_M, inter_Y_S, b = 0, m = 0, s = 1) 
  Xy_B0M1S0 = edit_dat_BMS(Xy, inter_Y_B, inter_Y_M, inter_Y_S, b = 0, m = 1, s = 0) 
  Xy_B0M0S0 = edit_dat_BMS(Xy, inter_Y_B, inter_Y_M, inter_Y_S, b = 0, m = 0, s = 0)
  
  Y_B0M1S1 = as.matrix(Xy_B0M1S1)%*%beta_Y
  Y_B0M0S1 = as.matrix(Xy_B0M0S1)%*%beta_Y
  Y_B0M1S0 = as.matrix(Xy_B0M1S0)%*%beta_Y
  Y_B0M0S0 = as.matrix(Xy_B0M0S0)%*%beta_Y
  Y_B0 = (1/denom)*{ Y_B0M1S1 * pB_M1S1 * pW_M1S1 * pA_M1S1 * pM1_S1 * pS + 
                     Y_B0M0S1 * pB_M0S1 * pW_M0S1 * pA_M0S1 * ( 1 - pM1_S1) * pS + 
                     Y_B0M1S0 * pB_M1S0 * pW_M1S0 * pA_M1S0 * pM1_S0 * (1 - pS) + 
                     Y_B0M0S0 * pB_M0S0 * pW_M0S0 * pA_M0S0 * ( 1 - pM1_S0) * (1 - pS) }

  # Optimal B
  Bopt = (Y_B1 > Y_B0) + 0 
  table(Bopt) 
  
  # (K = 1) ++++++++++++++++++++
 
  # Y | a = 1
  Xy_BoptS1M1W1A1 = edit_dat_BSMWA(Xy, inter_Y_B, inter_Y_S, inter_Y_M, inter_Y_W, inter_Y_A, b = Bopt, s = 1, m = 1, w = 1, a = 1) 
  Xy_BoptS1M1W0A1 = edit_dat_BSMWA(Xy, inter_Y_B, inter_Y_S, inter_Y_M, inter_Y_W, inter_Y_A, b = Bopt, s = 1, m = 1, w = 0, a = 1)
  Xy_BoptS1M0W1A1 = edit_dat_BSMWA(Xy, inter_Y_B, inter_Y_S, inter_Y_M, inter_Y_W, inter_Y_A, b = Bopt, s = 1, m = 0, w = 1, a = 1)
  Xy_BoptS1M0W0A1 = edit_dat_BSMWA(Xy, inter_Y_B, inter_Y_S, inter_Y_M, inter_Y_W, inter_Y_A, b = Bopt, s = 1, m = 0, w = 0, a = 1)
  Xy_BoptS0M1W1A1 = edit_dat_BSMWA(Xy, inter_Y_B, inter_Y_S, inter_Y_M, inter_Y_W, inter_Y_A, b = Bopt, s = 0, m = 1, w = 1, a = 1) 
  Xy_BoptS0M1W0A1 = edit_dat_BSMWA(Xy, inter_Y_B, inter_Y_S, inter_Y_M, inter_Y_W, inter_Y_A, b = Bopt, s = 0, m = 1, w = 0, a = 1)
  Xy_BoptS0M0W1A1 = edit_dat_BSMWA(Xy, inter_Y_B, inter_Y_S, inter_Y_M, inter_Y_W, inter_Y_A, b = Bopt, s = 0, m = 0, w = 1, a = 1)
  Xy_BoptS0M0W0A1 = edit_dat_BSMWA(Xy, inter_Y_B, inter_Y_S, inter_Y_M, inter_Y_W, inter_Y_A, b = Bopt, s = 0, m = 0, w = 0, a = 1)
  
  # W | a = 1
  Xw_A1M1S1 = edit_dat_AMS(Xw, inter_W_A, inter_W_M, inter_W_S, a = 1, m = 1, s = 1)
  pW1_A1M1S1 = as.matrix(Xw_A1M1S1)%*%beta_W
  
  Xw_A1M1S0 = edit_dat_AMS(Xw, inter_W_A, inter_W_M, inter_W_S, a = 1, m = 1, s = 0)
  pW1_A1M1S0 = as.matrix(Xw_A1M1S0)%*%beta_W
  
  Xw_A1M0S1 = edit_dat_AMS(Xw, inter_W_A, inter_W_M, inter_W_S, a = 1, m = 0, s = 1)
  pW1_A1M0S1 = as.matrix(Xw_A1M0S1)%*%beta_W
  
  Xw_A1M0S0 = edit_dat_AMS(Xw, inter_W_A, inter_W_M, inter_W_S, a = 1, m = 0, s = 0)
  pW1_A1M0S0 = as.matrix(Xw_A1M0S0)%*%beta_W

  Y_A1 = (as.matrix(Xy_BoptS1M1W1A1)%*%beta_Y) * pW1_A1M1S1 * pM1_S1 * pS + 
         (as.matrix(Xy_BoptS1M1W0A1)%*%beta_Y) * (1 - pW1_A1M1S1) * pM1_S1 * pS + 
         (as.matrix(Xy_BoptS1M0W1A1)%*%beta_Y) * pW1_A1M0S1 * (1 - pM1_S1) * pS + 
         (as.matrix(Xy_BoptS1M0W0A1)%*%beta_Y) * (1 - pW1_A1M0S1) * (1 - pM1_S1) * pS +  
         (as.matrix(Xy_BoptS0M1W1A1)%*%beta_Y) * pW1_A1M1S1 * pM1_S1 * (1 - pS) + 
         (as.matrix(Xy_BoptS0M1W0A1)%*%beta_Y) * (1 - pW1_A1M1S1) * pM1_S1 * (1 - pS) + 
         (as.matrix(Xy_BoptS0M0W1A1)%*%beta_Y) * pW1_A1M0S1 * (1 - pM1_S1) * (1 - pS) + 
         (as.matrix(Xy_BoptS0M0W0A1)%*%beta_Y) * (1 - pW1_A1M0S1) * ( 1 - pM1_S1) * ( 1 - pS) 
  
  
  # Y | a = 0
  Xy_BoptS1M1W1A0 = edit_dat_BSMWA(Xy, inter_Y_B, inter_Y_S, inter_Y_M, inter_Y_W, inter_Y_A, b = Bopt, s = 1, m = 1, w = 1, a = 0) 
  Xy_BoptS1M1W0A0 = edit_dat_BSMWA(Xy, inter_Y_B, inter_Y_S, inter_Y_M, inter_Y_W, inter_Y_A, b = Bopt, s = 1, m = 1, w = 0, a = 0)
  Xy_BoptS1M0W1A0 = edit_dat_BSMWA(Xy, inter_Y_B, inter_Y_S, inter_Y_M, inter_Y_W, inter_Y_A, b = Bopt, s = 1, m = 0, w = 1, a = 0)
  Xy_BoptS1M0W0A0 = edit_dat_BSMWA(Xy, inter_Y_B, inter_Y_S, inter_Y_M, inter_Y_W, inter_Y_A, b = Bopt, s = 1, m = 0, w = 0, a = 0)
  Xy_BoptS0M1W1A0 = edit_dat_BSMWA(Xy, inter_Y_B, inter_Y_S, inter_Y_M, inter_Y_W, inter_Y_A, b = Bopt, s = 0, m = 1, w = 1, a = 0) 
  Xy_BoptS0M1W0A0 = edit_dat_BSMWA(Xy, inter_Y_B, inter_Y_S, inter_Y_M, inter_Y_W, inter_Y_A, b = Bopt, s = 0, m = 1, w = 0, a = 0)
  Xy_BoptS0M0W1A0 = edit_dat_BSMWA(Xy, inter_Y_B, inter_Y_S, inter_Y_M, inter_Y_W, inter_Y_A, b = Bopt, s = 0, m = 0, w = 1, a = 0)
  Xy_BoptS0M0W0A0 = edit_dat_BSMWA(Xy, inter_Y_B, inter_Y_S, inter_Y_M, inter_Y_W, inter_Y_A, b = Bopt, s = 0, m = 0, w = 0, a = 0)
  
  # W | a = 0
  Xw_A0M1S1 = edit_dat_AMS(Xw, inter_W_A, inter_W_M, inter_W_S, a = 0, m = 1, s = 1)
  pW1_A0M1S1 = as.matrix(Xw_A0M1S1)%*%beta_W
  
  Xw_A0M1S0 = edit_dat_AMS(Xw, inter_W_A, inter_W_M, inter_W_S, a = 0, m = 1, s = 0)
  pW1_A0M1S0 = as.matrix(Xw_A0M1S0)%*%beta_W
  
  Xw_A0M0S1 = edit_dat_AMS(Xw, inter_W_A, inter_W_M, inter_W_S, a = 0, m = 0, s = 1)
  pW1_A0M0S1 = as.matrix(Xw_A0M0S1)%*%beta_W
  
  Xw_A0M0S0 = edit_dat_AMS(Xw, inter_W_A, inter_W_M, inter_W_S, a = 0, m = 0, s = 0)
  pW1_A0M0S0 = as.matrix(Xw_A0M0S0)%*%beta_W
  
  Y_A0 = (as.matrix(Xy_BoptS1M1W1A0)%*%beta_Y) * pW1_A0M1S1 * pM1_S1 * pS + 
   (as.matrix(Xy_BoptS1M1W0A1)%*%beta_Y) * (1 - pW1_A0M1S1) * pM1_S1 * pS + 
   (as.matrix(Xy_BoptS1M0W1A1)%*%beta_Y) * pW1_A0M0S1 * (1 - pM1_S1) * pS + 
   (as.matrix(Xy_BoptS1M0W0A1)%*%beta_Y) * (1 - pW1_A0M0S1) * (1 - pM1_S1) * pS +  
   (as.matrix(Xy_BoptS0M1W1A1)%*%beta_Y) * pW1_A0M1S1 * pM1_S1 * (1 - pS) + 
   (as.matrix(Xy_BoptS0M1W0A1)%*%beta_Y) * (1 - pW1_A0M1S1) * pM1_S1 * (1 - pS) + 
   (as.matrix(Xy_BoptS0M0W1A1)%*%beta_Y) * pW1_A0M0S1 * (1 - pM1_S1) * (1 - pS) + 
   (as.matrix(Xy_BoptS0M0W0A1)%*%beta_Y) * (1 - pW1_A0M0S1) * ( 1 - pM1_S1) * ( 1 - pS)
  
  # Optimal A
  Aopt = (Y_A1 > Y_A0) + 0
  table(Aopt)
  
  # # Ystar 
  # Xy_BoptM1 = edit_dY_BM(Xy, inter_Y_B, inter_Y_M, b = Bopt, m = 1) 
  # Xy_BoptM0 = edit_dY_BM(Xy, inter_Y_B, inter_Y_M, b = Bopt, m = 0) 
  # Y_BoptM1 = as.matrix(Xy_BoptM1)%*%beta_Y
  # Y_BoptM0 = as.matrix(Xy_BoptM0)%*%beta_Y
  # Ystar = (1/denom)*{ Y_BoptM1 * pW_M1 * pA_M1 * pM + Y_BoptM0 * pW_M0 * pA_M0 * (1-pM) }
  # 
  # dat$Ystar = Ystar 
  # beta_Ystar = lm(fmla_Ystar, dat)$coefficients
  # Xystar = as.data.frame(model.matrix(fmla_Ystar, data=model.frame(dat)))
  # 
  # # E[Ystar | ., A = 1] 
  # # A = 1, M = 1
  # Xystar_A1M1 = edit_dY_AM(Xystar, inter_Ystar_A, inter_Ystar_M, a = 1, m = 1)
  # Ystar_A1M1 = (as.matrix(Xystar_A1M1)%*%beta_Ystar)*pM 
  # # A = 1, M = 0
  # Xystar_A1M0 = edit_dY_AM(Xystar, inter_Ystar_A, inter_Ystar_M, a = 1, m = 0)
  # Ystar_A1M0 = (as.matrix(Xystar_A1M0)%*%beta_Ystar)*(1 - pM)
  # Ystar_A1 = Ystar_A1M1 + Ystar_A1M0
  # 
  # # E[Ystar | ., A = 0] 
  # # A = 0, M = 1
  # Xystar_A0M1 = edit_dY_AM(Xystar, inter_Ystar_A, inter_Ystar_M, a = 0, m = 1)
  # Ystar_A0M1 = (as.matrix(Xystar_A0M1)%*%beta_Ystar)*pM
  # # A = 0, M = 0
  # Xystar_A0M0 = edit_dY_AM(Xystar, inter_Ystar_A, inter_Ystar_M, a = 0, m = 0)
  # Ystar_A0M0 = (as.matrix(Xystar_A0M0)%*%beta_Ystar)*(1 - pM)
  # Ystar_A0 = Ystar_A0M1 + Ystar_A0M0
  # 
  # # Optimal A
  # Aopt = (Ystar_A1 > Ystar_A0) + 0
  # table(Aopt)
  
  return(list(Bopt=Bopt, Aopt=Aopt)) 
}


