

# +++++++++++++++++++++++++++++++++++++++
# Q-learning
# +++++++++++++++++++++++++++++++++++++++
Qlearning_unfair <- function(dat, fmla_Y, fmla_W, beta_Y, beta_W, 
                             inter_Y_B, inter_Y_W, inter_Y_A, inter_W_A){
 
  Xy = as.data.frame(model.matrix(fmla_Y, data=model.frame(dat)))
  
  
  # (K = 2) ++++++++++++++++++++

  # E[Y | ., B = 1] 
  Xy_B1 = edit_dat_B(Xy, inter_Y_B, b = 1)
  Y_B1 = as.matrix(Xy_B1)%*%beta_Y
  
  # E[Y | ., B = 0] 
  Xy_B0 = edit_dat_B(Xy, inter_Y_B, b = 0)
  Y_B0 = as.matrix(Xy_B0)%*%beta_Y
  
  # Optimal B
  Bopt = (Y_B1 > Y_B0) + 0
  table(Bopt)
  
  # (K = 1) ++++++++++++++++++++
  
  Xw = as.data.frame(model.matrix(fmla_W, data=model.frame(dat)))
  Xw_A1 = edit_dat_A(Xw, inter_W_A, a = 1) 
  pW_A1 = 1/(1+exp(-as.matrix(Xw_A1)%*%beta_W))
  
  Xw_A0 = edit_dat_A(Xw, inter_W_A, a = 0) 
  pW_A0 = 1/(1+exp(-as.matrix(Xw_A0)%*%beta_W))

  
  # Y | a = 1
  Xy_BoptW1A1 = edit_dat_BWA(Xy, inter_Y_B, inter_Y_W, inter_Y_A, b = Bopt, w = 1, a = 1) 
  Xy_BoptW0A1 = edit_dat_BWA(Xy, inter_Y_B, inter_Y_W, inter_Y_A, b = Bopt, w = 0, a = 1)
  
  Y_A1 = (as.matrix(Xy_BoptW1A1)%*%beta_Y)*pW_A1 + 
    (as.matrix(Xy_BoptW0A1)%*%beta_Y)*(1 - pW_A1) 
  
  # Y | a = 0
  Xy_BoptW1A0 = edit_dat_BWA(Xy, inter_Y_B, inter_Y_W, inter_Y_A, b = Bopt, w = 1, a = 0) 
  Xy_BoptW0A0 = edit_dat_BWA(Xy, inter_Y_B, inter_Y_W, inter_Y_A, b = Bopt, w = 0, a = 0)
  
  Y_A0 = (as.matrix(Xy_BoptW1A0)%*%beta_Y)*pW_A0 + 
  (as.matrix(Xy_BoptW0A0)%*%beta_Y)*(1 - pW_A0) 
  
  # Optimal A
  Aopt = (Y_A1 > Y_A0) + 0
  table(Aopt)
  
  # 
  # 
  # # Ystar 
  # Xy_Bopt = edit_dY_B(Xy, inter_Y_B, b = Bopt)
  # Ystar = as.matrix(Xy_Bopt)%*%beta_Y 
  # dat$Ystar = Ystar
  # beta_Ystar = lm(fmla_Ystar, dat)$coefficients
  # Xystar = as.data.frame(model.matrix(fmla_Ystar, data=model.frame(dat)))
  # 
  # # E[Ystar | ., A = 1] 
  # Xystar_A1 = edit_dY_A(Xystar, inter_Ystar_A, a = 1) 
  # Ystar_A1 = as.matrix(Xystar_A1)%*%beta_Ystar
  # 
  # # E[Ystar | ., A = 0] 
  # Xystar_A0 = edit_dY_A(Xystar, inter_Ystar_A, a = 0) 
  # Ystar_A0 = as.matrix(Xystar_A0)%*%beta_Ystar
  # 
  # # Optimal A
  # Aopt = (Ystar_A1 > Ystar_A0) + 0
  # table(Aopt)
  
  return(list(Bopt=Bopt, Aopt=Aopt))
}


Y_Qlearning_unfair <- function(dat, fmla_Y, beta_Y, inter_Y_B, inter_Y_A, Aopt, Bopt){ 
  
  # compute E[Y | everything + optimal decisions]
  Xy = as.data.frame(model.matrix(fmla_Y, data=model.frame(dat)))

  # E[Y | ., A = Aopt, B = Bopt] 
  Xy_Bopt = edit_dat_B(Xy, inter_Y_B, b = Bopt)
  Xy_AoptBopt = edit_dat_A(Xy_Bopt, inter_Y_A, a = Aopt)
  
  Y_AoptBopt = as.matrix(Xy_AoptBopt)%*%beta_Y
  
  return(Y_AoptBopt)
}

