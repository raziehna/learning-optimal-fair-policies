

# +++++++++++++++++++++++++++++++++++++++
# Vsearch
# +++++++++++++++++++++++++++++++++++++++
Vsearch <- function(dat, fmla_M, fmla_A, fmla_Y, beta_S, beta_M, beta_A, beta_Y, 
                    inter_M_S, inter_A_M, inter_A_S, fw){
  
  fw_A = fw$fw_A
  
  Xm = as.data.frame(model.matrix(fmla_M, data=model.frame(dat)))
  Xa = as.data.frame(model.matrix(fmla_A, data=model.frame(dat)))
  Xy = as.data.frame(model.matrix(fmla_Y, data=model.frame(dat)))
  A = dat$A
  Y = dat$Y

  # S 
  pS = beta_S
 
  # M, S = 1 
  Xm_S1 = edit_dat_S(Xm, inter_M_S, s = 1) 
  pM1_S1 = 1/(1+exp(-as.matrix(Xm_S1)%*%beta_M))
  pM1_S1[pM1_S1 < 0.005] = 0.005 
  pM1_S1[pM1_S1 > 0.995] = 0.995 
  
  # M, S = 0 
  Xm_S0 = edit_dat_S(Xm, inter_M_S, s = 0)
  pM1_S0 = 1/(1+exp(-as.matrix(Xm_S0)%*%beta_M))
  pM1_S0[pM1_S0 < 0.005] = 0.005 
  pM1_S0[pM1_S0 > 0.995] = 0.995  
  
  # E[Y(d)] where expectaion is with respect to p*(x, s, m) 
  Value_d = c()
  for (i in 1:length(fw_A)){
    dA = fw_A[[i]]
   
    # I(A = f(w_A)) 
    I_dA = (A == dA) + 0
    
    # p(dA = 1 | m, s, . ) 
    Xa_M1S1 = edit_dat_MS(Xa, inter_A_M, inter_A_S, m = 1, s = 1)
    PdA1_M1S1 = 1/(1+exp(-as.matrix(Xa_M1S1)%*%beta_A))
    Xa_M1S0 = edit_dat_MS(Xa, inter_A_M, inter_A_S, m = 1, s = 0)
    PdA1_M1S0 = 1/(1+exp(-as.matrix(Xa_M1S0)%*%beta_A))
    Xa_M0S1 = edit_dat_MS(Xa, inter_A_M, inter_A_S, m = 0, s = 1)
    PdA1_M0S1 = 1/(1+exp(-as.matrix(Xa_M0S1)%*%beta_A))
    Xa_M0S0 = edit_dat_MS(Xa, inter_A_M, inter_A_S, m = 0, s = 0)
    PdA1_M0S0 = 1/(1+exp(-as.matrix(Xa_M0S0)%*%beta_A))
    
    # p(dA | m, s, . ) 
    PdA_M1S1 = PdA1_M1S1
    PdA_M1S1[dA == 0] = 1 - PdA1_M1S1[dA == 0] 
    PdA_M1S1[PdA_M1S1 < 0.005] = 0.005
    PdA_M1S1[PdA_M1S1 > 0.995] = 0.995

    PdA_M1S0 = PdA1_M1S0
    PdA_M1S0[dA == 0] = 1 - PdA1_M1S0[dA == 0] 
    PdA_M1S0[PdA_M1S0 < 0.005] = 0.005
    PdA_M1S0[PdA_M1S0 > 0.995] = 0.995
    
    PdA_M0S1 = PdA1_M0S1
    PdA_M0S1[dA == 0] = 1 - PdA1_M0S1[dA == 0] 
    PdA_M0S1[PdA_M0S1 < 0.005] = 0.005
    PdA_M0S1[PdA_M0S1 > 0.995] = 0.995
    
    PdA_M0S0 = PdA1_M0S0
    PdA_M0S0[dA == 0] = 1 - PdA1_M0S0[dA == 0] 
    PdA_M0S0[PdA_M0S0 < 0.005] = 0.005
    PdA_M0S0[PdA_M0S0 > 0.995] = 0.995
    
    
    E_Yd_M1S1 =  mean((I_dA*Y/PdA_M1S1)*pM1_S1*pS)
    E_Yd_M1S0 =  mean((I_dA*Y/PdA_M1S0)*pM1_S0*(1 - pS))
    E_Yd_M0S1 =  mean((I_dA*Y/PdA_M0S1)*(1- pM1_S1)*pS)
    E_Yd_M0S0 =  mean((I_dA*Y/PdA_M0S0)*(1- pM1_S0)*(1 - pS))
   
    # expected value of Y(d) with respect to p(x, s)
    E_Yd = E_Yd_M1S1 + E_Yd_M1S0 + E_Yd_M0S1 + E_Yd_M0S0
    
    # Value(d)
    Value_d = c(Value_d, E_Yd)
  }
  
  (Value_max = which.max(Value_d))
  Aopt = fw_A[[Value_max]]
  table(Aopt)
  
  return(list(Aopt=Aopt, Value_max=Value_max))
}







