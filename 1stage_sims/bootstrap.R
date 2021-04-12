

# BOOTSTRAP
bootstrap_policy <- function(iter, method, policy){
 # set.seed(0)
 
 n = nrow(dat)

 Aopt_mat = matrix(0, nrow = n, ncol = (iter+1))
 Bopt_mat = matrix(0, nrow = n, ncol = (iter+1))
 Yhat_vec = c()
 Yopt_vec = c()
 
 if (method == "Q-learning"){
  if (policy == "Unfair"){
   # (Unfair Policy)
   policy_uf_q = Qlearning_unfair(dat, fmla_Y, fmla_W, beta_Y_uf, beta_W, 
                                  inter_Y_B, inter_Y_W, inter_Y_A, inter_W_A)
   Aopt = policy_uf_q$Aopt
   Bopt = policy_uf_q$Bopt
   Yopt = Y_Qlearning_unfair(dat, fmla_Y, beta_Y_uf, inter_Y_B, inter_Y_A, Aopt, Bopt)

  }else{ 
   # (Fair Policy)
   policy_f_q = Qlearning(dat, fmla_M, fmla_A, fmla_W, fmla_B, fmla_Y, beta_S_f, beta_M_f, beta_A, 
                          beta_W, beta_Y, inter_A_M, inter_W_M, inter_W_A, inter_Y_B, inter_Y_M)
   Aopt = policy_f_q$Aopt 
   Bopt = policy_f_q$Bopt
   Yopt = Y_Qlearning_unfair(dat, fmla_Y, beta_Y_uf, inter_Y_B, inter_Y_A, Aopt, Bopt)
  }
 }

 if (method == "V-search"){
  fw = policy(alpha, dat)
  if (policy == "Unfair"){
   # (Unfair Policy)
   result_uf_v = Vsearch_unfair(dat, fmla_A, fmla_B, fmla_Y, beta_A, beta_B, beta_Y, fw)
   Aopt = result_uf_v$Aopt
   Bopt = result_uf_v$Bopt
   Yopt = Y_Vsearch_unfair(dat, fmla_A, fmla_B, fmla_Y, beta_A, beta_B, beta_Y, Aopt, Bopt)
  }else{
   # (Fair Policy)
   result_f_v = Vsearch(dat, fmla_M, fmla_A, fmla_B, fmla_Y, beta_S_f, beta_M_f, beta_A, beta_B, 
                        beta_Y, inter_M_S, inter_A_M, inter_A_S, inter_B_M, inter_B_S, fw)
   Aopt = result_f_v$Aopt
   Bopt = result_f_v$Bopt
   Yopt = Y_Vsearch_unfair(dat, fmla_A, fmla_B, fmla_Y, beta_A, beta_B, beta_Y, Aopt, Bopt)
  }
 }

 Aopt_mat[, 1] = Aopt 
 Bopt_mat[, 1] = Bopt
 Yhat_vec = c(Yhat_vec, mean(dat$Y))
 Yopt_vec = c(Yopt_vec, mean(Yopt))
 
 cat("\n \n")
 for (i in 1:iter){
  cat("Bootstrap # ", i, "\n")
  idx = sample.int(n, n, replace = TRUE)
  
  if (method == "Q-learning"){
   if (policy == "Unfair"){
    # (Unfair Policy)
    policy_uf_q = Qlearning_unfair(dat[idx, ], fmla_Y, fmla_W, beta_Y_uf, beta_W, 
                                   inter_Y_B, inter_Y_W, inter_Y_A, inter_W_A)
    Aopt = policy_uf_q$Aopt
    Bopt = policy_uf_q$Bopt
    Yopt = Y_Qlearning_unfair(dat[idx, ], fmla_Y, beta_Y_uf, inter_Y_B, inter_Y_A, Aopt, Bopt)
   }else{ 
    # (Fair Policy)
    policy_f_q = Qlearning(dat[idx, ], fmla_M, fmla_A, fmla_W, fmla_B, fmla_Y, beta_S_f, beta_M_f, beta_A, 
                           beta_W, beta_Y, inter_A_M, inter_W_M, inter_W_A, inter_Y_B, inter_Y_M)
    Aopt = policy_f_q$Aopt 
    Bopt = policy_f_q$Bopt
    Yopt = Y_Qlearning_unfair(dat[idx, ], fmla_Y, beta_Y_uf, inter_Y_B, inter_Y_A, Aopt, Bopt)
   }
  }
  if (method == "V-search"){
   fw = policy(alpha, dat[idx, ])  
   if (policy == "Unfair"){
    # (Unfair Policy)
    result_uf_v = Vsearch_unfair(dat[idx, ], fmla_A, fmla_B, fmla_Y, beta_A, beta_B, beta_Y, fw)
    Aopt = result_uf_v$Aopt
    Bopt = result_uf_v$Bopt
    Yopt = Y_Vsearch_unfair(dat[idx, ], fmla_A, fmla_B, fmla_Y, beta_A, beta_B, beta_Y, Aopt, Bopt)
   }else{
    # (Fair Policy)
    result_f_v = Vsearch(dat[idx, ], fmla_M, fmla_A, fmla_B, fmla_Y, beta_S_f, beta_M_f, beta_A, beta_B, 
                         beta_Y, inter_M_S, inter_A_M, inter_A_S, inter_B_M, inter_B_S, fw)
    Aopt = result_f_v$Aopt
    Bopt = result_f_v$Bopt
    Yopt = Y_Vsearch_unfair(dat[idx, ], fmla_A, fmla_B, fmla_Y, beta_A, beta_B, beta_Y, Aopt, Bopt)
   }
  }
  
  Aopt_mat[, (i+1)] = Aopt
  Bopt_mat[, (i+1)] = Bopt
  Yhat_vec = c(Yhat_vec, mean(dat[idx, ]$Y))
  Yopt_vec = c(Yopt_vec, mean(Yopt))
 }
 
 cat("\n \n")
 return(list(Aopt = Aopt_mat, 
             Bopt = Bopt_mat, 
             Yhat_vec = Yhat_vec, 
             Yopt_vec = Yopt_vec))
}






