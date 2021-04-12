

# +++++++++++++++++++++++++++++++++++++++
# G-estimation 
# +++++++++++++++++++++++++++++++++++++++
Gestimation_unfair <- function(dat, fmla_A, fmla_Y, beta_A, beta_Y, inter_Y_A, W){
  
  Xa = as.data.frame(model.matrix(fmla_A, data=model.frame(dat)))
  Xy = as.data.frame(model.matrix(fmla_Y, data=model.frame(dat)))
  A = dat$A
  Y = dat$Y

  # p(A = 1 | x, s, m)
  pA1 = 1/(1+exp(-as.matrix(Xa)%*%beta_A))
  pA1[pA1 < 0.005] = 0.005
  pA1[pA1 > 0.995] = 0.995
  
  # Calculate E[Y | x, s, m] 
  # A = 1
  Xy_A1 = edit_dat_A(Xy, inter_Y_A, a = 1)
  Y_A1 = (as.matrix(Xy_A1)%*%beta_Y)*pA1
  # A = 0
  Xy_A0 = edit_dat_A(Xy, inter_Y_A, a = 0)
  Y_A0 = (as.matrix(Xy_A0)%*%beta_Y)*(1 - pA1) 
  
  Y_hat = Y_A0 + Y_A1 
  
  # W should be a function of (x, s, m)
  # W = as.matrix(cbind( sqrt(abs(dat[, 1:5])) ))  # remember to change it in Y_Gestimation
  
  C = as.matrix((A - pA1)*(Y - Y_hat))
  B1 = -as.matrix(A - pA1)^2
  B2 = W*as.vector(B1) 
  
  # (B*hw)*psi = -C*hw 
  hw1 = as.matrix(rep(1, nrow(dat)))
  hw2 = W
  hw = cbind(hw1, hw2)
  C = -t(hw)%*%C
  B1 = t(hw)%*%B1 
  B2 = t(hw)%*%B2
  B = cbind(B1, B2) 
  
  # solve: B*psi = C 
  # gamma = psi_1*A + psi_2*A*W(x, s, m)
  (psi = solve(B, C))
  psi_1 = psi[1]
  psi_2 = psi[-1]
  
  # find optimal treatment 
  Aopt = as.numeric(psi_1 + W%*%psi_2 > 0 ) 
  table(Aopt)
 
  return(list(Aopt=Aopt, 
              psi_1=psi_1, 
              psi_2=psi_2))
}


Y_Gestimation_unfair <- function(dat, fmla_Y, beta_Y, inter_Y_A, Aopt, psi_1, psi_2, W){
 
  Xy = as.data.frame(model.matrix(fmla_Y, data=model.frame(dat))) 
  
  # E[Y | x, s, m, a = 0]
  Xy_A0 = edit_dat_A(Xy, inter_Y_A, a = 0)
  Y_A0 = as.matrix(Xy_A0)%*%beta_Y
 
  # W = as.matrix(cbind( sqrt(abs(dat[, 1:5])) ))
  Y_Aopt = as.matrix(Aopt)*psi_1 + (W%*%psi_2)*as.matrix(Aopt) + Y_A0
  
  return(Y_Aopt)
}

