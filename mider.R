
#X <- read.csv(file, header = FALSE)
mider <- function(X, taumax, ert_crit){
  ntotal <- ncol(X)
  n <- nrow(X)
  taumax <- 10
  fraction  = 0.1*(log10(n)-1)
  q <- 1
  threshold <- 1
  ert_crit <- 3
  
  MI <- H2 <- array(0, dim = c(ntotal, ntotal, taumax + 1))
  
  for (tau in 0:taumax){
    for (i in 1:ntotal){
      for (j in 1:ntotal){
        pb <- 5
        #print("here0")
        calc <- estimateH2(X[1:(n-tau), i], X[(tau+1):n, j], pb, q)
        MI[i, j, tau+1] <- calc$mutinfo; fracn <- calc$fracn; H2[i, j, tau+1] <- calc$H2
        #if (tau == 0){print(calc$H2)}
        while (fracn < fraction){
          pb = pb+1
          calc <- estimateH2(X[1:(n-tau), i], X[(tau+1):n, j], pb, q)
          MI[i, j, tau+1] <- calc$mutinfo; fracn <- calc$fracn; H2[i, j, tau+1] <- calc$H2
        }
      }
    }
  }
  
  #H2[1,,]
  
  #### Entropy of each variable (H1)
  
  H1 <- diag(H2[, , 1])
  
  
  #### EMC distance
  EMC_dist_tau = exp(-MI)
  dist <- sapply(1:ntotal, function(i) sapply(1:ntotal, function(j) min(EMC_dist_tau[j, i, ])))
  taumin <- sapply(1:ntotal, function(i) sapply(1:ntotal, function(j) which.min(EMC_dist_tau[j, i, ])))
  
  
  ### conditional entropies
  
  cond_entr2 <- cond_entr2_A <- cond_entr2_B <- entropy_error <- rel_entr_reduc_2 <- matrix(0, ntotal, ntotal)
  
  for (m in 1:ntotal){
    for (n in 1:ntotal){
      cond_entr2_A[m,n] <- H1[m] - MI[m,n,1]
      cond_entr2_B[m,n] <- H2[n, m, 1] - H1[n]
      cond_entr2[m, n] <- cond_entr2_A[m,n]
      entropy_error[m, n] <- abs(cond_entr2_A[m, n] - cond_entr2_B[m, n])
      rel_entr_reduc_2[m, n] <- MI[m, n, 1]/H1[m]
    }
  }
  
  #### ERT 1st round
  
  if (ert_crit >= 1) {
    spec_ind_2 <- apply(cond_entr2 + 1000*diag(ntotal), 2, which.min)
    ERT_array_2 <- cbind(1:ntotal, spec_ind_2)
  }
  
  
  ### ERT 2nd round
  
  if (ert_crit >= 2) {
    
    MI3 <- H3 <- array(0, dim = c(ntotal, ntotal, ntotal))
    
    for (i in 1:ntotal){
      for (j in 1:ntotal){
        for (k in 1:ntotal){
          pb <- 5
          calc <- estimateH3(X[, i], X[, j], X[, k],  pb, q)
          H3[i, j, k] <- calc$H3; fracn <- calc$fracn
          #if (tau == 0){print(calc$H2)}
          while (fracn < fraction){
            pb = pb+1
            calc <- estimateH3(X[, i], X[, j], X[, k],  pb, q)
            H3[i, j, k] <- calc$H3; fracn <- calc$fracn
          }
        }
      }
    }
    
    for (i in 1:ntotal){
      for (j in 1:ntotal){
        for (k in 1:ntotal){
          MI3[i, j, k] <- -H1[i]-H1[j]+2*H2[i,j, 1]+H2[i,k, 1]+H2[j,k, 1]-2*H3[i,j,k]
        }
      }
    }
    
    cond_entr3 <- rel_entr_reduc3 <- array(0, dim = c(ntotal, ntotal, ntotal))
    for (m in 1:ntotal){
      for (n in 1:ntotal){
        for (p in 1:ntotal){
          cond_entr3[m, n, p] <- H3[n, p, m] - H2[n,p,1]
          rel_entr_reduc3[m, n, p] <- (cond_entr2[m, n] - cond_entr3[m, n, p] )/ H1[m]
        }
      }
    }
    
    diag_ones_3 <- array(0, dim = c(ntotal, ntotal, ntotal))
    for(i in 1:ntotal){
      diag_ones_3[i, i, ] <- diag_ones_3[i, , i] <- diag_ones_3[, i, i] <- 1
    }
    
    cond_entr3m <- cond_entr3 + 1000*diag_ones_3
    min_diff_3 <- spec_ind_3 <- matrix(0, ntotal, 1)
    
    for (i in 1:ntotal){
      min_diff_3[i] <- min(cond_entr3m[ERT_array_2[i, 1], ERT_array_2[i, 2], ])
      spec_ind_3[i] <- which.min(cond_entr3m[ERT_array_2[i, 1], ERT_array_2[i, 2], ])
    }
    
    ERT_array_3 <- cbind(1:ntotal, spec_ind_2, spec_ind_3)
  }
  
  
  ### ERT 3rd round
  
  if (ert_crit >= 3) {
    
    H4 <- array(0, dim = c(ntotal, ntotal, ntotal, ntotal))
    
    for (i in 1:ntotal){
      for (j in 1:ntotal){
        for (k in 1:ntotal){
          for (l in 1:ntotal){
            pb <- 5
            calc <- estimateH4(X[, i], X[, j], X[, k], X[, l], pb, q)
            H4[i, j, k, l] <- calc$H4; fracn <- calc$fracn
            #if (tau == 0){print(calc$H2)}
            while (fracn < fraction){
              pb = pb+1
              calc <- estimateH4(X[, i], X[, j], X[, k], X[, l], pb, q)
              H4[i, j, k, l] <- calc$H4; fracn <- calc$fracn
            }
          }
        }
      }
    }
    
    cond_entr4 <- rel_entr_reduc4 <- array(0, dim = c(ntotal, ntotal, ntotal, ntotal))
    for (m in 1:ntotal){
      for (n in 1:ntotal){
        for (p in 1:ntotal){
          for (r in 1:ntotal){
            cond_entr4[m, n, p, r] <- H4[n, p, r, m] - H3[n, p ,r]
            rel_entr_reduc4[m, n, p, r] <- (cond_entr3[m, n, p] - cond_entr4[m, n, p, r] )/ H1[m]
          }
        }
      }
    }
    
    diag_ones_4 <- array(0, dim = c(ntotal, ntotal, ntotal, ntotal))
    for(i in 1:ntotal){
      diag_ones_4[i, i, , ] <- diag_ones_4[i, , i, ] <- diag_ones_4[i, , , i] <- 1
      diag_ones_4[ , i, i, ] <- diag_ones_4[, i, , i] <- diag_ones_4[ , , i, i] <- 1
    }
    
    cond_entr4m <- cond_entr4 + 1000*diag_ones_4
    min_diff_4 <- spec_ind_4<- matrix(0, ntotal, 1)
    
    for (i in 1:ntotal){
      min_diff_4[i] <- min(cond_entr4m[ERT_array_3[i, 1], ERT_array_3[i, 2], ERT_array_3[i, 3], ])
      spec_ind_4[i] <- which.min(cond_entr4m[ERT_array_3[i, 1], ERT_array_3[i, 2], ERT_array_3[i, 3], ])
    }
    
    ERT_array_4 <- cbind(1:ntotal, spec_ind_2, spec_ind_3, spec_ind_4)
  }
  
  
  ###### transfer entropy
  n <- nrow(X)
  T <- T3 <- T4 <- matrix(0, ntotal, ntotal)
  T1 <- T2 <- matrix(0, ntotal, 1)
  
  for(i in 1:ntotal){
    for (j in 1:ntotal){
      pb <- 5
      calc <- estimateH2(X[(1+taumin[i,j]):n, j], X[1:(n-taumin[i, j]), j], pb, q)
      T1[j] <- calc$H2; fracn <- calc$fracn
      while (fracn < fraction){
        pb = pb+1
        calc <- estimateH2(X[(1+taumin[i,j]):n, j], X[1:(n-taumin[i, j]), j], pb, q)
        T1[j] <- calc$H2; fracn <- calc$fracn
      }
      
      pb <- 5
      calc <- estimateH2(X[1:(n-taumin[i, j]), j], X[1:(n-taumin[i, j]), j], pb, q)
      T2[j] <- calc$H2; fracn <- calc$fracn
      while (fracn < fraction){
        pb = pb+1
        calc <- estimateH2(X[1:(n-taumin[i, j]), j], X[1:(n-taumin[i, j]), j], pb, q)
        T2[j] <- calc$H2; fracn <- calc$fracn
      }
      
      pb <- 5
      calc <- estimateH3(X[(1+taumin[i,j]):n, j], X[1:(n-taumin[i, j]), j], X[1:(n-taumin[i, j]), i], pb, q)
      T3[i, j] <- calc$H3; fracn <- calc$fracn
      while (fracn < fraction){
        pb = pb+1
        calc <- estimateH3(X[(1+taumin[i,j]):n, j], X[1:(n-taumin[i, j]), j], X[1:(n-taumin[i, j]), i], pb, q)
        T3[i, j] <- calc$H3; fracn <- calc$fracn
      }
      
      pb <- 5
      calc <- estimateH2(X[1:(n-taumin[i, j]), j], X[1:(n-taumin[i, j]), i], pb, q)
      T4[i, j] <- calc$H2; fracn <- calc$fracn
      while (fracn < fraction){
        pb = pb+1
        calc <- estimateH2(X[1:(n-taumin[i, j]), j], X[1:(n-taumin[i, j]), i], pb, q)
        T4[i, j] <- calc$H2; fracn <- calc$fracn
      }
    }
  }
  
  for (i in 1:ntotal){
    for (j in 1:ntotal){
      T[i, j] <- T1[j] - T2[j] - T3[j,i] + T4[j,i]
    }
  }
  
  ###### connections
  
  con_array <- matrix(0, ntotal, ntotal)
  
  if (ert_crit >= 2){
    for (i in 1:ntotal){
      con_array[i, ERT_array_2[i,2]] = rel_entr_reduc_2[i, ERT_array_2[i,2]]
    }
  }
  
  if (ert_crit >= 3){
    for (i in 1:ntotal){
      con_array[i, ERT_array_3[i,3]] = rel_entr_reduc3[i, ERT_array_2[i,2], 
                                                       ERT_array_3[i,3]]
    }
  }
  
  if (ert_crit >= 4){
    for (i in 1:ntotal){
      con_array[i, ERT_array_4[i,4]] = rel_entr_reduc4[i, ERT_array_2[i,2], 
                                                       ERT_array_3[i,3],
                                                       ERT_array_4[i,4]]
    }
  }
  
  
  if (threshold==1){
    maxMI <- max(con_array)
    threshold <- ifelse(maxMI < 0.3, 0, ifelse(maxMI > 0.7, 0.2, 0.5*(maxMI-0.3)))
  }
  
  #print(threshold); #print(con_array)
  #threshold <- 0
  con_array[con_array < threshold] <- 0
  
  return(list("con_array" = con_array,
              "T" = T,
              "dist" = dist))
}

