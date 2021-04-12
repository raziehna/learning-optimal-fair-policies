

# +++++++++++++++++++++++++++++++++++++++
# Q-learning
# +++++++++++++++++++++++++++++++++++++++
Qlearning_unfair <- function(dat, fmla_Y, beta_Y, inter_Y_A){
 
  Xy = as.data.frame(model.matrix(fmla_Y, data=model.frame(dat)))
 
  # E[Y | ., A = 1] 
  Xy_A1 = edit_dat_A(Xy, inter_Y_A, a = 1)
  Y_A1 = as.matrix(Xy_A1)%*%beta_Y
  
  # E[Y | ., A = 0] 
  Xy_A0 = edit_dat_A(Xy, inter_Y_A, a = 0)
  Y_A0 = as.matrix(Xy_A0)%*%beta_Y
  
  # Optimal A
  Aopt = (Y_A1 > Y_A0) + 0
  table(Aopt)
 
  return(list(Aopt=Aopt))
}


Y_Qlearning_unfair <- function(dat, fmla_Y, beta_Y, inter_Y_A, Aopt){ 
  
  # compute E[Y | everything + optimal decisions]
  Xy = as.data.frame(model.matrix(fmla_Y, data=model.frame(dat)))

  # E[Y | ., A = Aopt] 
  Xy_Aopt = edit_dat_A(Xy, inter_Y_A, a = Aopt)
  
  Y_Aopt = as.matrix(Xy_Aopt)%*%beta_Y
  
  return(Y_Aopt)
}

