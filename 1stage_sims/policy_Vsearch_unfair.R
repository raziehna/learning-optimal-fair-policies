
# +++++++++++++++++++++++++++++++++++++++
# Vsearch
# +++++++++++++++++++++++++++++++++++++++
Vsearch_unfair <- function(dat, fmla_A, fmla_Y, beta_A, beta_Y, fw){
  
  fw_A = fw$fw_A
 
  Xa = as.data.frame(model.matrix(fmla_A, data=model.frame(dat)))
  Xy = as.data.frame(model.matrix(fmla_Y, data=model.frame(dat)))
  A = dat$A
  Y = dat$Y
  
  Value_d = c()
  for (i in 1:length(fw_A)){
    dA = fw_A[[i]]
    
    # I(A = f(w_A)) 
    I_dA = (A == dA) + 0
    
    # p(dA = 1 | . ) 
    PdA1 = 1/(1+exp(-as.matrix(Xa)%*%beta_A))
    # p(dA | . ) 
    PdA = PdA1
    PdA[dA == 0] = 1 - PdA1[dA == 0] 
    PdA[PdA < 0.05] = 0.05
    PdA[PdA > 0.995] = 0.995
    
    # Value(d)
    s = sum(I_dA/PdA)
    Value_d = c(Value_d, (1/s)*sum((I_dA*Y)/PdA))
  }
  
  (Value_max = which.max(Value_d))
  Aopt = fw_A[[Value_max]]
  table(Aopt)
  
  return(list(Aopt=Aopt, Value_max=Value_max))
}


Y_Vsearch_unfair <- function(dat, fmla_A, fmla_Y, beta_A, beta_Y, Aopt){
 
  Xa = as.data.frame(model.matrix(fmla_A, data=model.frame(dat)))
  Xy = as.data.frame(model.matrix(fmla_Y, data=model.frame(dat)))
  A = dat$A
  Y = dat$Y
  
  dA = Aopt
  
  # I(A = dA)
  I_dA = (A == dA) + 0
  
  # p(dA = 1 | . )   
  PdA1 = 1/(1+exp(-as.matrix(Xa)%*%beta_A))
  # p(dA | . ) 
  PdA = PdA1
  PdA[dA == 0] = 1 - PdA1[dA == 0] 
  PdA[PdA < 0.005] = 0.005
  PdA[PdA > 0.995] = 0.995
  
  # Value(d)
  # s = sum((I_dA*I_dB)/(PdA*PdB))
  Y_opt = (I_dA*Y)/PdA
  
  return(Y_opt)
}





