

# +++++++++++++++++++++++++++++++++++++++
# G-estimation
# +++++++++++++++++++++++++++++++++++++++

Gestimation <- function(dat, fmla_M, beta_S, beta_M, inter_M_S, psi_1, psi_2, W){

 W = as.data.frame(cbind( sqrt(abs(dat[, 1:5])) ))

 # m = 1, s = 1
 W$M = 1
 W$S = 1
 W_11 = as.matrix(W)

 # W_MS
 # m = 1, s = 0
 W$M = 1
 W$S = 0
 W_10 = as.matrix(W)

 # m = 0, s = 1
 W$M = 0
 W$S = 1
 W_01 = as.matrix(W)

 # m = 0, s = 0
 W$M = 0
 W$S = 0
 W_00 = as.matrix(W)


 # Work on S model
 # S
 pS = beta_S

 # Woek on M model
 Xm = as.data.frame(model.matrix(fmla_M, data=model.frame(dat)))
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

 A = rep(1, n)
 gamma_star_A1 = as.matrix(A)*psi_1 +
                {(W_11%*%psi_2)*as.matrix(A)*pM1_S1*pS +
                 (W_10%*%psi_2)*as.matrix(A)*pM1_S0*(1 - pS) +
                  (W_01%*%psi_2)*as.matrix(A)*( 1 - pM1_S1)*pS +
                   (W_00%*%psi_2)*as.matrix(A)*( 1 - pM1_S1)*( 1- pS)}

 A = rep(0, n)
 gamma_star_A0 = as.matrix(A)*psi_1 +
 {(W_11%*%psi_2)*as.matrix(A)*pM1_S1*pS +
   (W_10%*%psi_2)*as.matrix(A)*pM1_S0*(1 - pS) +
   (W_01%*%psi_2)*as.matrix(A)*( 1 - pM1_S1)*pS +
   (W_00%*%psi_2)*as.matrix(A)*( 1 - pM1_S1)*( 1- pS)}


 # Optimal A
 Aopt = (gamma_star_A1 > gamma_star_A0) + 0
 table(Aopt)

 return(Aopt)
}




# # +++++++++++++++++++++++++++++++++++++++
# # G-estimation 
# # +++++++++++++++++++++++++++++++++++++++
# Gestimation <- function(dat, fmla_A, fmla_M, fmla_Y, beta_S, beta_A, beta_M, beta_Y,
#                         inter_Y_A, inter_Y_M, inter_Y_S, inter_A_M, inter_A_S, k,){
# 
#   n = nrow(dat)
# 
#   ks = 1
#   k = 1
# 
#   # Work on S model
#   pS = beta_S
#   s_star = rbinom(ks*n, 1, pS)
# 
#   # create 2 replica of dat and row-bind them + interaction data frames
#   inter_A_S_new = data.frame(inter_A_S[rep(seq_len(nrow(inter_A_S)), each = ks), ])
#   inter_Y_S_new = data.frame(inter_Y_S[rep(seq_len(nrow(inter_Y_S)), each = ks), ])
#   colnames(inter_A_S_new) = colnames(inter_A_S)
#   colnames(inter_Y_S_new) = colnames(inter_Y_S)
# 
#   # update interactions with s*
#   dat_S = dat[rep(seq_len(nrow(dat)), each = ks), ]
#   dat_S = edit_dat_S(dat_S, inter_A_S_new, s = s_star)
#   dat_S = edit_dat_S(dat_S, inter_Y_S_new, s = s_star)
# 
# 
#   # Work on M model
#   Xm = as.data.frame(model.matrix(fmla_M, data=model.frame(dat)))
#   pM = 1/(1+exp(-as.matrix(Xm)%*%beta_M))
#   pM[pM < 0.005] = 0.005
#   pM[pM > 0.995] = 0.995
#   m_star = c()
#   for (i in 1:n) m_star = c(m_star, rbinom(k*ks, 1, pM[i]))
# 
#   # create 100 replica of dat and row-bind them + interaction data frames
#   inter_A_M_new = data.frame(inter_A_M[rep(seq_len(nrow(inter_A_M)), each = k*ks), ])
#   inter_Y_M_new = data.frame(inter_Y_M[rep(seq_len(nrow(inter_Y_M)), each = k*ks), ])
#   inter_Y_A_new = data.frame(inter_Y_A[rep(seq_len(nrow(inter_Y_A)), each = k*ks), ])
#   colnames(inter_A_M_new) = colnames(inter_A_M)
#   colnames(inter_Y_M_new) = colnames(inter_Y_M)
#   colnames(inter_Y_A_new) = colnames(inter_Y_A)
# 
#   # update interactions with m*
#   dat_M = dat_S[rep(seq_len(nrow(dat_S)), each = k), ]
#   dat_M = edit_dat_M(dat_M, inter_A_M_new, m = m_star)
#   dat_new = edit_dat_M(dat_M, inter_Y_M_new, m = m_star)
# 
#   # ++++++++++++++++++++++
# 
#   Xa = as.data.frame(model.matrix(fmla_A, data=model.frame(dat_new)))
#   Xy = as.data.frame(model.matrix(fmla_Y, data=model.frame(dat_new)))
#   A = dat_new$A
#   Y = dat_new$Y
# 
#   # p(A = 1 | x, s*, m*)
#   pA1 = 1/(1+exp(-as.matrix(Xa)%*%beta_A))
#   pA1[pA1 < 0.005] = 0.005
#   pA1[pA1 > 0.995] = 0.995
# 
# 
#   # Calculate E[Y | x, s*, m*]
#   # A = 1
#   Xy_A1 = edit_dat_A(Xy, inter_Y_A_new, a = 1)
#   Y_A1 = (as.matrix(Xy_A1)%*%beta_Y)*pA1
#   # A = 0
#   Xy_A0 = edit_dat_A(Xy, inter_Y_A_new, a = 0)
#   Y_A0 = (as.matrix(Xy_A0)%*%beta_Y)*(1 - pA1)
# 
#   Y_hat = Y_A0 + Y_A1
# 
#   # W should be a function of (x, s, m)
#   W = as.matrix(cbind( sqrt(abs(dat_new[, 1:5])) ))  # remember to change it in Y_Gestimation
#   C = as.matrix((A - pA1)*(Y - Y_hat))
#   B1 = -as.matrix(A - pA1)^2
#   B2 = W*as.vector(B1)
# 
#   # (B*hw)*psi = -C*hw
#   hw1 = as.matrix(rep(1, nrow(dat_new)))
#   hw2 = W
#   hw = cbind(hw1, hw2)
#   C = -t(hw)%*%C
#   B1 = t(hw)%*%B1
#   B2 = t(hw)%*%B2
#   B = cbind(B1, B2)
# 
#   # solve: B*psi = C
#   # gamma = psi_1*A + psi_2*A*W(x, s, m)
#   (psi = solve(B, C))
#   psi_1 = psi[1]
#   psi_2 = psi[-1]
# 
# 
#   # Find optimal treatment
#   W_new = as.matrix(cbind( sqrt(abs(dat_new[, 1:5])) ))
#   Xm_new = as.data.frame(model.matrix(fmla_M, data=model.frame(dat_new)))
#   pM_new = 1/(1+exp(-as.matrix(Xm_new)%*%beta_M))
#   pM_new[pM_new < 0.005] = 0.005
#   pM_new[pM_new > 0.995] = 0.995
# 
# 
#   # \gamma(., m = 1)*(pM)
#   W_new[, 5] = 1
#   gamma_M1 = (psi_1 + W_new%*%psi_2)*pM_new
# 
#   # \gamma(., m = 0)*(1 - pM)
#   W_new[, 5] = 0
#   gamma_M0 = (psi_1 + W_new%*%psi_2)*(1 - pM_new)
# 
#   # sum out m from gamma
#   gamma = gamma_M1 + gamma_M0
# 
#   Aopt = (gamma > 0 ) + 0
#   table(Aopt)
# 
#   # find expected outcome under Aopt
#   # Xy = as.data.frame(model.matrix(fmla_Y, data=model.frame(dat_new)))
#   # inter_Y_A_new = data.frame(inter_Y_A[rep(seq_len(nrow(inter_Y_A)), each = k*2), ])
#   # colnames(inter_Y_A_new) = colnames(inter_Y_A)
# 
#   # E[Y | x, s, m, a = 0]
#   Xy_new = as.data.frame(model.matrix(fmla_Y, data=model.frame(dat_new)))
#   Xy_A0 = edit_dat_A(Xy_new, inter_Y_A, a = 0)
#   Y_A0 = as.matrix(Xy_A0)%*%beta_Y
# 
#   W_new = as.matrix(cbind( sqrt(abs(dat_new[, 1:5])) ))
#   Y_Aopt = as.matrix(Aopt)*psi_1 + (W_new%*%psi_2)*as.matrix(Aopt) + Y_A0
# 
#   mean(Y_Aopt)
# 
#   return(list(Aopt=Aopt,
#               psi_1=psi_1,
#               psi_2=psi_2))
# }



