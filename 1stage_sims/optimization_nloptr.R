

# -------------------------------------- 
# Optimize a constrained logistic model
# --------------------------------------
optimize_nloptr <- function(dat, fmla_M, funcSY, funcSA, inter_M_S, 
                            tau_cont_u, tau_cont_l, tau_bin_u, tau_bin_l, beta_start){
  
  Xm = as.matrix(model.matrix(fmla_M, data=model.frame(dat)))
  S = as.matrix(dat$S)
  M = as.matrix(dat$M)
  
  if (length(beta_start) == 0){
    beta_start = rep(0.1, 1 + ncol(Xm))
    names(beta_start) = c("s", colnames(Xm))
  }
  
  # Define the negative log likelihood function
  eval_f <- function(beta){
    beta_S = beta[1]
    beta_M = beta[2:length(beta)]
    names(beta_M) = colnames(Xm)

    f = - sum(S == 1)*log(beta_S) - sum(S == 0)*log(1 - beta_S) + M*log(1+exp(-Xm%*%beta_M)) + (1-M)*log(1+exp(Xm%*%beta_M)) 
    f = sum(f)/nrow(dat)
    return(f)
  }
  
  # Define the inequlity constraint 
  eval_g_ineq <- function(beta){
    beta_S = beta[1]
    beta_M = beta[2:length(beta)]
    names(beta_M) = colnames(Xm)
    
    # eff_SY - tau_u < 0, tau_l - eff_SY < 0
    eff_SY = funcSY(dat, beta_S, beta_M, inter_M_S)
    eval_SY =  c(eff_SY - tau_cont_u, tau_cont_l - eff_SY)
    
    # eff_SA - tau_u < 0, tau_l - eff_SA < 0
    eff_SA = funcSA(dat, beta_S, beta_M, inter_M_S)
    eval_SA =  c(eff_SA - tau_bin_u, tau_bin_l - eff_SA)
    
    # combine all four constraints 
    eval_g = c(eval_SY, eval_SA)
    return(eval_g)
  }

  # Solve the optimization problem
  mle = nloptr(x0=beta_start, 
               eval_f=eval_f, 
               eval_g_ineq=eval_g_ineq,
               opts = list("algorithm"="NLOPT_LN_COBYLA","xtol_rel"=1.0e-8, "maxeval"=5000)
               )
  
  # Returnt the parameters
  beta = mle$solution
  beta_S = beta[1]
  beta_M = beta[2:length(beta)]
  names(beta_M) = colnames(Xm)
  
  return(list(beta_S=beta_S, beta_M=beta_M))
}
