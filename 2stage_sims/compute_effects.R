
# ------------------------------------------------------- 
# PSE of S on Y using S and M models
# ------------------------------------------------------- 

# g*h for all g_i in g and h_i in h
add_interactions <- function(g, h){
  gh_inter = h
  for (i in 1:ncol(g)){
    for (j in 1:ncol(h)){
      gh_inter[paste0(colnames(g)[i], colnames(h)[j])] = g[, i]*h[, j]
    }
  }
  names = colnames(gh_inter)[-c(1:ncol(h))]
  gh_inter = data.frame(gh_inter[, -c(1:ncol(h))])
  colnames(gh_inter) = names
  return(gh_inter)
}


# edit_dM
# edit_dY_S
edit_dat_S <- function(dat, inter_set, s){
  dat$S = s 
  g = data.frame(S = dat$S) 
  h = inter_set
  if (length(h)!=0){
   inter_S = add_interactions(g, h)
   for (i in 1:length(inter_S)){
    col_name = colnames(inter_S)[i]
    if (nchar(col_name) == 2){
     colnames(inter_S)[i] = paste(sort(unlist(strsplit(col_name, ""))), collapse = "")
    }
   }
   for (i in 1:ncol(inter_S)){
    col_id = which(colnames(dat) %in% colnames(inter_S[i]))
    dat[, col_id] = inter_S[i]
   }
  }
  return(dat)
}


# edit_dA_M
# edit_dW_M
edit_dat_M <- function(dat, inter_set, m){
  dat$M = m
  g = data.frame(M = dat$M) 
  h = inter_set
  if (length(h)!=0){
    inter_M = add_interactions(g, h)
    for (i in 1:length(inter_M)){
      col_name = colnames(inter_M)[i]
      if (nchar(col_name) == 2){
        colnames(inter_M)[i] = paste(sort(unlist(strsplit(col_name, ""))), collapse = "")
      }
    }
    for (i in 1:ncol(inter_M)){
     col_id = which(colnames(dat) %in% colnames(inter_M[i]))
     dat[, col_id] = inter_M[i]
    }
  }
  return(dat)
}


# edit_dY_A
edit_dat_A <- function(dat, inter_set, a){ 
 dat$A = a
 g = data.frame(A = dat$A) 
 h = inter_set
 if (length(h)!=0){
  inter_A = add_interactions(g, h)
  for (i in 1:length(inter_A)){
   col_name = colnames(inter_A)[i]
   if (nchar(col_name) == 2){
    colnames(inter_A)[i] = paste(sort(unlist(strsplit(col_name, ""))), collapse = "")
   }
  } 
  for (i in 1:ncol(inter_A)){
   col_id = which(colnames(dat) %in% colnames(inter_A[i]))
   dat[, col_id] = inter_A[i]
  }
 }
 return(dat)
}


# edit_dY_B
edit_dat_B <- function(dat, inter_set, b){ 
 dat$B = b
 g = data.frame(B = dat$B) 
 h = inter_set
 if (length(h)!=0){
  inter_B = add_interactions(g, h)
  for (i in 1:length(inter_B)){
   col_name = colnames(inter_B)[i]
   if (nchar(col_name) == 2){
    colnames(inter_B)[i] = paste(sort(unlist(strsplit(col_name, ""))), collapse = "")
   }
  }
  for (i in 1:ncol(inter_B)){
   col_id = which(colnames(dat) %in% colnames(inter_B[i]))
   dat[, col_id] = inter_B[i]
  }
 }
 return(dat)
}


# edit_dA_MS
# edit_dW_MS
# edit_dB_MS
edit_dat_MS <- function(dA, inter_A_M, inter_A_S, m, s){
 dA$M = m
 g2 = data.frame(M = dA$M) 
 h2 = inter_A_M
 if (length(h2)!=0){
  inter_M = add_interactions(g2, h2)
  for (i in 1:length(inter_M)){
   col_name = colnames(inter_M)[i]
   if (nchar(col_name) == 2){
    colnames(inter_M)[i] = paste(sort(unlist(strsplit(col_name, ""))), collapse = "")
   }
  }
  for (i in 1:ncol(inter_M)){
   col_id = which(colnames(dA) %in% colnames(inter_M[i]))
   dA[, col_id] = inter_M[i]
  }
 }
 
 dA$S = s
 g1 = data.frame(S = dA$S) 
 h1 = inter_A_S
 if (length(h1)!=0){
  inter_S = add_interactions(g1, h1)
  for (i in 1:length(inter_S)){
   col_name = colnames(inter_S)[i]
   if (nchar(col_name) == 2){
    colnames(inter_S)[i] = paste(sort(unlist(strsplit(col_name, ""))), collapse = "")
   }
  }
  for (i in 1:ncol(inter_S)){
   col_id = which(colnames(dA) %in% colnames(inter_S[i]))
   dA[, col_id] = inter_S[i]
  }
 }
 
 if ("MS" %in% colnames(dA)){
  dA$MS = m*s
 }
return(dA)
}


# edit_dY_AM
edit_dat_AM <- function(dY, inter_Y_A, inter_Y_M, a, m){ 
  dY$A = a
  g1 = data.frame(A = dY$A) 
  h1 = inter_Y_A
  if (length(h1)!=0){
    inter_A = add_interactions(g1, h1)
    for (i in 1:length(inter_A)){
      col_name = colnames(inter_A)[i]
      if (nchar(col_name) == 2){
        colnames(inter_A)[i] = paste(sort(unlist(strsplit(col_name, ""))), collapse = "")
      }
    }
    for (i in 1:ncol(inter_A)){
     col_id = which(colnames(dY) %in% colnames(inter_A[i]))
     dY[, col_id] = inter_A[i]
    }
  }
  
  dY$M = m
  g2 = data.frame(M = dY$M) 
  h2 = inter_Y_M
  if (length(h2)!=0){
    inter_M = add_interactions(g2, h2)
    for (i in 1:length(inter_M)){
      col_name = colnames(inter_M)[i]
      if (nchar(col_name) == 2){
        colnames(inter_M)[i] = paste(sort(unlist(strsplit(col_name, ""))), collapse = "")
      }
    }
    for (i in 1:ncol(inter_M)){
     col_id = which(colnames(dY) %in% colnames(inter_M[i]))
     dY[, col_id] = inter_M[i]
    }
  }
  if ("AM" %in% colnames(dY)){
   dY$AM = a*m
  }
  return(dY)
}


# edit_dY_BM
edit_dat_BM <- function(dY, inter_Y_B, inter_Y_M, b, m){ 
  dY$B = b
  g1 = data.frame(B = dY$B) 
  h1 = inter_Y_B
  if (length(h1)!=0){
    inter_B = add_interactions(g1, h1)
    for (i in 1:length(inter_B)){
      col_name = colnames(inter_B)[i]
      if (nchar(col_name) == 2){
        colnames(inter_B)[i] = paste(sort(unlist(strsplit(col_name, ""))), collapse = "")
      }
    }
    for (i in 1:ncol(inter_B)){
     col_id = which(colnames(dY) %in% colnames(inter_B[i]))
     dY[, col_id] = inter_B[i]
    }
  }
  
  dY$M = m
  g2 = data.frame(M = dY$M) 
  h2 = inter_Y_M
  if (length(h2)!=0){
    inter_M = add_interactions(g2, h2)
    for (i in 1:length(inter_M)){
      col_name = colnames(inter_M)[i]
      if (nchar(col_name) == 2){
        colnames(inter_M)[i] = paste(sort(unlist(strsplit(col_name, ""))), collapse = "")
      }
    }
    for (i in 1:ncol(inter_M)){
     col_id = which(colnames(dY) %in% colnames(inter_M[i]))
     dY[, col_id] = inter_M[i]
    }
  }
  
  if ("BM" %in% colnames(dY)){
   dY$BM = b*m
  }
  return(dY)
}


# edit_dW_AMS
edit_dat_AMS <- function(dW, inter_W_A, inter_W_M, inter_W_S, a, m, s){
 dW$A = a
 g3 = data.frame(A = dW$A) 
 h3 = inter_W_A
 if (length(h3)!=0){
  inter_A = add_interactions(g3, h3)
  for (i in 1:length(inter_A)){
   col_name = colnames(inter_A)[i]
   if (nchar(col_name) == 2){
    colnames(inter_A)[i] = paste(sort(unlist(strsplit(col_name, ""))), collapse = "")
   }
  }
  for (i in 1:ncol(inter_A)){
   col_id = which(colnames(dW) %in% colnames(inter_A[i]))
   dW[, col_id] = inter_A[i]
  }
 }
 
 dW$M = m
 g2 = data.frame(M = dW$M) 
 h2 = inter_W_M
 if (length(h2)!=0){
  inter_M = add_interactions(g2, h2)
  for (i in 1:length(inter_M)){
   col_name = colnames(inter_M)[i]
   if (nchar(col_name) == 2){
    colnames(inter_M)[i] = paste(sort(unlist(strsplit(col_name, ""))), collapse = "")
   }
  }
  for (i in 1:ncol(inter_M)){
   col_id = which(colnames(dW) %in% colnames(inter_M[i]))
   dW[, col_id] = inter_M[i]
  }
 }
 
 dW$S = s
 g1 = data.frame(S = dW$S) 
 h1 = inter_W_S
 if (length(h1)!=0){
  inter_S = add_interactions(g1, h1)
  for (i in 1:length(inter_S)){
   col_name = colnames(inter_S)[i]
   if (nchar(col_name) == 2){
    colnames(inter_S)[i] = paste(sort(unlist(strsplit(col_name, ""))), collapse = "")
   }
  }
  for (i in 1:ncol(inter_S)){
   col_id = which(colnames(dW) %in% colnames(inter_S[i]))
   dW[, col_id] = inter_S[i]
  }
 }
 if ("AM" %in% colnames(dW)) dW$AM = a*m
 if ("AS" %in% colnames(dW)) dW$AS = a*s
 if ("MS" %in% colnames(dW)) dW$MS = m*s
 
 return(dW)
}


# edit_dY_BMS
edit_dat_BMS <- function(dY, inter_Y_B, inter_Y_M, inter_Y_S, b, m, s){ 
 dY$B = b
 g1 = data.frame(B = dY$B) 
 h1 = inter_Y_B
 if (length(h1)!=0){
  inter_B = add_interactions(g1, h1)
  for (i in 1:length(inter_B)){
   col_name = colnames(inter_B)[i]
   if (nchar(col_name) == 2){
    colnames(inter_B)[i] = paste(sort(unlist(strsplit(col_name, ""))), collapse = "")
   }
  }
  for (i in 1:ncol(inter_B)){
   col_id = which(colnames(dY) %in% colnames(inter_B[i]))
   dY[, col_id] = inter_B[i]
  }
 }
 
 dY$M = m
 g2 = data.frame(M = dY$M) 
 h2 = inter_Y_M
 if (length(h2)!=0){
  inter_M = add_interactions(g2, h2)
  for (i in 1:length(inter_M)){
   col_name = colnames(inter_M)[i]
   if (nchar(col_name) == 2){
    colnames(inter_M)[i] = paste(sort(unlist(strsplit(col_name, ""))), collapse = "")
   }
  }
  for (i in 1:ncol(inter_M)){
   col_id = which(colnames(dY) %in% colnames(inter_M[i]))
   dY[, col_id] = inter_M[i]
  }
 }
 
 dY$S = s
 g3 = data.frame(S = dY$S) 
 h3 = inter_Y_S
 if (length(h3)!=0){
  inter_S = add_interactions(g3, h3)
  for (i in 1:length(inter_S)){
   col_name = colnames(inter_S)[i]
   if (nchar(col_name) == 2){
    colnames(inter_S)[i] = paste(sort(unlist(strsplit(col_name, ""))), collapse = "")
   }
  }
  for (i in 1:ncol(inter_S)){
   col_id = which(colnames(dY) %in% colnames(inter_S[i]))
   dY[, col_id] = inter_S[i]
  }
 }
 
 if ("BM" %in% colnames(dY)) dY$AM = b*m
 if ("BS" %in% colnames(dY)) dY$AS = b*s
 if ("MS" %in% colnames(dY)) dY$MS = m*s
 
 return(dY)
}


# edit_dY_BWA
edit_dat_BWA <- function(dY, inter_Y_B, inter_Y_W, inter_Y_A, b, w, a){ 
 dY$A = a
 g4 = data.frame(A = dY$A) 
 h4 = inter_Y_A
 if (length(h4)!=0){
  inter_A = add_interactions(g4, h4)
  for (i in 1:length(inter_A)){
   col_name = colnames(inter_A)[i]
   if (nchar(col_name) == 2){
    colnames(inter_A)[i] = paste(sort(unlist(strsplit(col_name, ""))), collapse = "")
   }
  }
  for (i in 1:ncol(inter_A)){
   col_id = which(colnames(dY) %in% colnames(inter_A[i]))
   dY[, col_id] = inter_A[i]
  }
 }
 
  dY$B = b
  g1 = data.frame(B = dY$B) 
  h1 = inter_Y_B
  if (length(h1)!=0){
    inter_B = add_interactions(g1, h1)
    for (i in 1:length(inter_B)){
      col_name = colnames(inter_B)[i]
      if (nchar(col_name) == 2){
        colnames(inter_B)[i] = paste(sort(unlist(strsplit(col_name, ""))), collapse = "")
      }
    }
    for (i in 1:ncol(inter_B)){
     col_id = which(colnames(dY) %in% colnames(inter_B[i]))
     dY[, col_id] = inter_B[i]
    }
  }
  
  dY$W = w
  g3 = data.frame(W = dY$W) 
  h3 = inter_Y_W
  if (length(h3)!=0){
    inter_W = add_interactions(g3, h3)
    for (i in 1:length(inter_W)){
      col_name = colnames(inter_W)[i]
      if (nchar(col_name) == 2){
        colnames(inter_W)[i] = paste(sort(unlist(strsplit(col_name, ""))), collapse = "")
      }
    }
    for (i in 1:ncol(inter_W)){
     col_id = which(colnames(dY) %in% colnames(inter_W[i]))
     dY[, col_id] = inter_W[i]
    }
  }
  
  if ("AB" %in% colnames(dY)) dY$AB = a*b
  if ("AW" %in% colnames(dY)) dY$AW = a*w
  if ("BW" %in% colnames(dY)) dY$BW = b*w
  
  return(dY)
}



# edit_dY_BSMWA
edit_dat_BSMWA <- function(dY, inter_Y_B, inter_Y_S, inter_Y_M, inter_Y_W, inter_Y_A, b, s, m, w, a){ 
 dY$A = a
 g4 = data.frame(A = dY$A) 
 h4 = inter_Y_A
 if (length(h4)!=0){
  inter_A = add_interactions(g4, h4)
  for (i in 1:length(inter_A)){
   col_name = colnames(inter_A)[i]
   if (nchar(col_name) == 2){
    colnames(inter_A)[i] = paste(sort(unlist(strsplit(col_name, ""))), collapse = "")
   }
  }
  for (i in 1:ncol(inter_A)){
   col_id = which(colnames(dY) %in% colnames(inter_A[i]))
   dY[, col_id] = inter_A[i]
  }
 }
 
 dY$B = b
 g1 = data.frame(B = dY$B) 
 h1 = inter_Y_B
 if (length(h1)!=0){
  inter_B = add_interactions(g1, h1)
  for (i in 1:length(inter_B)){
   col_name = colnames(inter_B)[i]
   if (nchar(col_name) == 2){
    colnames(inter_B)[i] = paste(sort(unlist(strsplit(col_name, ""))), collapse = "")
   }
  }
  for (i in 1:ncol(inter_B)){
   col_id = which(colnames(dY) %in% colnames(inter_B[i]))
   dY[, col_id] = inter_B[i]
  }
 }
 
 dY$M = m
 g2 = data.frame(M = dY$M) 
 h2 = inter_Y_M
 if (length(h2)!=0){
  inter_M = add_interactions(g2, h2)
  for (i in 1:length(inter_M)){
   col_name = colnames(inter_M)[i]
   if (nchar(col_name) == 2){
    colnames(inter_M)[i] = paste(sort(unlist(strsplit(col_name, ""))), collapse = "")
   }
  }
  for (i in 1:ncol(inter_M)){
   col_id = which(colnames(dY) %in% colnames(inter_M[i]))
   dY[, col_id] = inter_M[i]
  }
 }
 
 dY$S = s
 g5 = data.frame(S = dY$S) 
 h5 = inter_Y_S
 if (length(h5)!=0){
  inter_S = add_interactions(g5, h5)
  for (i in 1:length(inter_S)){
   col_name = colnames(inter_S)[i]
   if (nchar(col_name) == 2){
    colnames(inter_S)[i] = paste(sort(unlist(strsplit(col_name, ""))), collapse = "")
   }
  }
  for (i in 1:ncol(inter_S)){
   col_id = which(colnames(dY) %in% colnames(inter_S[i]))
   dY[, col_id] = inter_S[i]
  }
 }
 
 dY$W = w
 g3 = data.frame(W = dY$W) 
 h3 = inter_Y_W
 if (length(h3)!=0){
  inter_W = add_interactions(g3, h3)
  for (i in 1:length(inter_W)){
   col_name = colnames(inter_W)[i]
   if (nchar(col_name) == 2){
    colnames(inter_W)[i] = paste(sort(unlist(strsplit(col_name, ""))), collapse = "")
   }
  }
  for (i in 1:ncol(inter_W)){
   col_id = which(colnames(dY) %in% colnames(inter_W[i]))
   dY[, col_id] = inter_W[i]
  }
 }
 
 edit_dat_BSMWA
 
 if ("AB" %in% colnames(dY)) dY$AB = a*b
 if ("AM" %in% colnames(dY)) dY$AM = a*m
 if ("AS" %in% colnames(dY)) dY$AS = a*s
 if ("AW" %in% colnames(dY)) dY$AW = a*w
 if ("BM" %in% colnames(dY)) dY$BM = b*m
 if ("BS" %in% colnames(dY)) dY$BS = b*s
 if ("BW" %in% colnames(dY)) dY$BW = b*w
 if ("MS" %in% colnames(dY)) dY$MS = m*s
 if ("MW" %in% colnames(dY)) dY$MW = m*w
 if ("SW" %in% colnames(dY)) dY$SW = s*w
 
 return(dY)
}


# ------------------------------------------------------- 
# PSE of S on Y using S and M and Y models
# -------------------------------------------------------
# if ("SX1" %in% attributes(beta_M)$names){
#   dM$SX1 = dat$X1
# }
pse_SY_SMY <- function(dat, beta_Y, beta_M, inter_M_S, inter_Y_S){
  n = nrow(dat)
  
  # S
  p_S1 = mean(dat$S)
  p_S0 = 1 - p_S1
  
  # Identity vecotrs 
  I_S1 = dat$S
  I_S0 = 1 - I_S1
  
  # Work on Y model
  # Y_hat: means that use E[Y | .] instead of Y
  idx_Y = match(attributes(beta_Y)$names[-1], colnames(dat))
  dY = as.data.frame(cbind(intercept = 1, dat[,idx_Y]))
  # S = 1 
  dY_S1 = edit_dY_S(dY, inter_Y_S, s = 1)
  Yhat_S1 = as.matrix(dY_S1)%*%beta_Y
  # S = 0 
  dY_S0 = edit_dY_S(dY, inter_Y_S, s = 0)
  Yhat_S0 = as.matrix(dY_S0)%*%beta_Y
  
  # Work on M model 
  idx_M = match(attributes(beta_M)$names[-1], colnames(dat))
  dM = cbind(intercept = 1, dat[,idx_M])
  # S = 1
  dM_S1 = edit_dM(dM, inter_M_S, s = 1) 
  p_M1_S1 = (1+exp(-as.matrix(dM_S1)%*%beta_M))^(-1)
  
  p_M_S1 = p_M1_S1 
  p_M_S1[M == 0] = 1 - p_M1_S1[M == 0] 
  
  # S = 0
  dM_S0 = edit_dM(dM, inter_M_S, s = 0)
  p_M1_S0 = (1+exp(-as.matrix(dM_S0)%*%beta_M))^(-1)
  
  p_M_S0 = p_M1_S0
  p_M_S0[M == 0] = 1 - p_M1_S0[M == 0] 
  
  # eff = mean( ( (I_S1*Y)/p_S1 )*(p_m_S0/p_m_S1) - I_S0*Y/p_S0 )
  
  w = sum(((I_S1)/p_S1)*(p_M_S0/p_M_S1))
  eff = sum(((I_S1*Yhat_S1)/p_S1)*(p_M_S0/p_M_S1))/w - mean(I_S0*Yhat_S0/p_S0)
  
  return(eff)
}



# ------------------------------------------------------- 
# PSE of S on Y using S and M models
# -------------------------------------------------------
pse_SY_SM <- function(dat, beta_S, beta_M, inter_M_S){
 n = nrow(dat)
 
 # S
 p_S1 = beta_S
 p_S0 = 1 - p_S1
 
 # Identity vecotrs 
 I_S1 = dat$S
 I_S0 = 1 - I_S1
 
 # Work on M model 
 idx_M = match(attributes(beta_M)$names[-1], colnames(dat))
 dM = cbind(intercept = 1, dat[,idx_M])
 # S = 1
 dM_S1 = edit_dat_S(dM, inter_M_S, s = 1) 
 p_M1_S1 = (1+exp(-as.matrix(dM_S1)%*%beta_M))^(-1)
 
 p_M_S1 = p_M1_S1 
 p_M_S1[M == 0] = 1 - p_M1_S1[M == 0] 
 
 # S = 0
 dM_S0 = edit_dat_S(dM, inter_M_S, s = 0)
 p_M1_S0 = (1+exp(-as.matrix(dM_S0)%*%beta_M))^(-1)
 
 p_M_S0 = p_M1_S0
 p_M_S0[M == 0] = 1 - p_M1_S0[M == 0] 
 
 # eff = mean( ( (I_S1*Y)/p_S1 )*(p_m_S0/p_m_S1) - I_S0*Y/p_S0 )
 
 w = sum(((I_S1)/p_S1)*(p_M_S0/p_M_S1))
 eff = sum(((I_S1*Y)/p_S1)*(p_M_S0/p_M_S1))/w - mean(I_S0*Y/p_S0)
 
 return(eff)
}




# ------------------------------------------------------- 
# PSE of S on B using S and M models
# -------------------------------------------------------
pse_SB_SM <- function(dat, beta_S, beta_M, inter_M_S){
  n = nrow(dat)
  
  # S
  # p_S1 = mean(dat$S)
  p_S1 = beta_S
  p_S0 = 1 - p_S1
  
  # Identity vecotrs 
  I_S1 = dat$S
  I_S0 = 1 - I_S1
  
  # B vector 
  B = dat$B
  
  # Work on M model 
  idx_M = match(attributes(beta_M)$names[-1], colnames(dat))
  dM = cbind(intercept = 1, dat[,idx_M])
  # S = 1
  dM_S1 = edit_dat_S(dM, inter_M_S, s = 1)
  p_M1_S1 = (1+exp(-as.matrix(dM_S1)%*%beta_M))^(-1)
  
  p_M_S1 = p_M1_S1 
  p_M_S1[M == 0] = 1 - p_M1_S1[M == 0] 
  
  # S = 0
  dM_S0 = edit_dat_S(dM, inter_M_S, s = 0)
  p_M1_S0 = (1+exp(-as.matrix(dM_S0)%*%beta_M))^(-1) 
  
  p_M_S0 = p_M1_S0
  p_M_S0[M == 0] = 1 - p_M1_S0[M == 0] 
  
  # eff = mean( ( (I_S1*Y)/p_S1 )*(p_m_S0/p_m_S1) - I_S0*Y/p_S0 )
  
  w = sum(((I_S1)/p_S1)*(p_M_S0/p_M_S1))
  pB_S1 = sum(((I_S1*B)/p_S1)*(p_M_S0/p_M_S1))/w
  pB_S0 = mean(I_S0*B/p_S0)
  
  eff = (pB_S1 / (1 - pB_S1)) / (pB_S0/(1 - pB_S0))
  
  return(eff)
}


# ------------------------------------------------------- 
# PSE of S on A using S and M models
# -------------------------------------------------------
pse_SA_SM <- function(dat, beta_S, beta_M, inter_M_S){
  n = nrow(dat)
  
  # S
  # p_S1 = mean(dat$S)
  p_S1 = beta_S
  p_S0 = 1 - p_S1
  
  # Identity vecotrs 
  I_S1 = dat$S
  I_S0 = 1 - I_S1
  
  # A vector 
  A = dat$A
  
  # Work on M model 
  idx_M = match(attributes(beta_M)$names[-1], colnames(dat))
  dM = cbind(intercept = 1, dat[,idx_M])
  # S = 1
  dM_S1 = edit_dat_S(dM, inter_M_S, s = 1)
  p_M1_S1 = (1+exp(-as.matrix(dM_S1)%*%beta_M))^(-1)
  
  p_M_S1 = p_M1_S1 
  p_M_S1[M == 0] = 1 - p_M1_S1[M == 0]
  
  # S = 0
  dM_S0 = edit_dat_S(dM, inter_M_S, s = 0)
  p_M1_S0 = (1+exp(-as.matrix(dM_S0)%*%beta_M))^(-1)
  
  p_M_S0 = p_M1_S0
  p_M_S0[M == 0] = 1 - p_M1_S0[M == 0] 
  
  # eff = mean( ( (I_S1*Y)/p_S1 )*(p_m_S0/p_m_S1) - I_S0*Y/p_S0 )
  
  w = sum(((I_S1)/p_S1)*(p_M_S0/p_M_S1))
  pA_S1 = sum(((I_S1*A)/p_S1)*(p_M_S0/p_M_S1))/w 
  pA_S0 = mean(I_S0*A/p_S0)
  
  eff = (pA_S1 / (1 - pA_S1)) / (pA_S0/(1 - pA_S0))
  
  return(eff)
}

