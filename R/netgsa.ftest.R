netgsa.ftest <-
function(s2g, s2e, D, DtD, DDInv, n_vec, B, beta_hat){
  ncond = length(D)
  npath = dim(B)[1] 
  p = dim(B)[2]
  Ip = diag(rep(1, p))
  
  teststat2 = matrix(0, npath, 1)
  df2 = matrix(0, npath, 1)
  
  ####C = inverse of Psi ' W^{-1} Psi is calculated in the main code
  Cmat = lapply(1:ncond, function(k) (s2g*Ip + s2e * DDInv[[k]])/n_vec[k])
  Cmat = as.matrix(bdiag(Cmat))
  
  for (rr in 1:npath){  	   
    ##obtain the contrast matrix L for a given pathway
    L = get.contrast(D, B[rr,]) 
    
    ##get the test statistic
    LCL = L %*% Cmat %*% t(L)
    LCL_inv = solve(LCL)
    q = rankMatrix(L)		
    teststat2[rr] = (t(unlist(beta_hat)) %*% t(L) %*% LCL_inv %*% L %*% unlist(beta_hat))/q
    
    ##find projection P and diagonal D such that LCL' = P'DP
    tmp <- eigen(LCL)
    D_diag <- diag(tmp$values)		
    
    ## Calculate the first-order derivatives of C wrt to parameters s2g and s2e	
    g1Mat = as.matrix(bdiag(lapply(1:ncond, function(ix) Ip/n_vec[ix])))
    g2Mat = as.matrix(bdiag(lapply(1:ncond, function(ix) DDInv[[ix]]/n_vec[ix])))
    
    ##calculate the empirical covariance matrix Kmat
    Sigma = lapply(1:ncond, function(ix) s2e * Ip + s2g * DtD[[ix]] )
    SigmaInv = lapply(1:ncond, function(ix) chol2inv(chol(Sigma[[ix]])) )
    SigmaInvD = lapply(1:ncond, function(ix) SigmaInv[[ix]] %*% DtD[[ix]] )
    SinvSinv = lapply(1:ncond, function(ix) matTr(SigmaInv[[ix]] %*% SigmaInv[[ix]]))
    SinvDSinvD = lapply(1:ncond, function(ix) matTr(SigmaInvD[[ix]] %*% SigmaInvD[[ix]]))
    SinvSinvD = lapply(1:ncond, function(ix) matTr(SigmaInv[[ix]] %*% SigmaInvD[[ix]]))
    EH11 = (1/2) * Reduce("+", SinvDSinvD)
    EH12 = (1/2) * Reduce("+", SinvSinvD)
    EH22 = (1/2) * Reduce("+", SinvSinv)
    
    Kmat = matrix(c(EH11, EH12, EH12, EH22), 2, 2, byrow = TRUE)
    KmatInv = solve(Kmat)
    
    Em <- 0				
    for(m in 1:q){
      lm <- matrix(L[m,], nrow=1)		## lm is the mth row of L
      gm1 <- lm %*% g1Mat %*% t(lm)
      gm2 <- lm %*% g2Mat %*% t(lm)
      gm <- c(gm1, gm2)  ## gm is the gradient of lm C lm' wrt theta
      vm <- (2*D_diag[m,m]^2) / (t(gm) %*% KmatInv %*% gm )
      if(vm > 2){
        Em <- Em + (vm/(vm-2))
      }
    }
    
    ##The first degree of freedom is q. We need to calculate the 2nd degree of freedom df2. 
    if(Em > q){
      df2[rr] = (2*Em) / (Em - q)
    }
    
    ## NOTE: matlab only accepts positive integers as df's for F-dist.
    ## Therefore, I had set the denom df to 1 if it is zero, otherwise round it to integer values
    df2[rr] <- (df2[rr] >= 1) * ceiling(df2[rr]) + (df2[rr] < 1) * 1
  }
  
  pvals = 1 - pf(abs(teststat2), q, df2) + pf(-abs(teststat2), q, df2)  
  
  return(list(teststat = teststat2, df = df2, p.value = pvals))
}
