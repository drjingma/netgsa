netgsa.ftest <-
function(s2g, s2e, D, DDt, DtDInv, n_vec, B, beta){
  ncond = length(n_vec)
  npath = nrow(B) 
  teststat = matrix(0, npath, 1)
  df = matrix(0, npath, 1)
  p = length(beta[[1]])
  
  ##calculate the empirical covariance matrix Kmat
  ## Sigma is W in the notes. 
  Sigma = lapply(DDt, lapply, function(ix) s2e * diag(1, nrow(ix)) + s2g * ix)
  SigmaInv = lapply(Sigma, lapply, function(ix) chol2inv(chol(ix)))
  SigmaInvD = lapply(1:ncond, function(j) mapply(function(a,b) crossprod(a,b), SigmaInv[[j]], DDt[[j]], SIMPLIFY = FALSE))
  SinvSinv = lapply(SigmaInv, sapply, function(ix) matTr(ix, ix))
  SinvDSinvD = lapply(SigmaInvD, sapply, function(ix) matTr(ix, ix))
  SinvSinvD = lapply(1:ncond, function(j) mapply(function(a,b) matTr(a,b), SigmaInv[[j]], SigmaInvD[[j]]))
  EH11 = 0.5*sum(sapply(SinvDSinvD,sum)) 
  EH12 = 0.5*sum(sapply(SinvSinvD,sum))
  EH22 = 0.5*sum(sapply(SinvSinv,sum))
  
  Kmat = matrix(c(EH11, EH12, EH12, EH22), 2, 2, byrow = TRUE)
  KmatInv = solve(Kmat)
  
  for (rr in 1:npath){  	   
    ##obtain the contrast matrix L for a given pathway
    LN_list = lapply(1:ncond, function(j) crossprod(B[rr,], as.matrix(bdiag(D[[j]]))) * B[rr,])
    
    ##calculate ll' and l(Lambda'Lambda)^{-1}l' for each row of L
    llt <- sapply(1:ncond,function(j) tcrossprod(LN_list[[j]])/n_vec[j]) 
    lDtDlt <- sapply(1:ncond ,function(j) LN_list[[j]]%*%as.matrix(bdiag(DtDInv[[j]]))%*%t(LN_list[[j]])/n_vec[j])
    Lbeta_full <- sapply(1:ncond, function(j) crossprod(t(LN_list[[j]]), beta[[j]]))
    Lbeta <- Lbeta_full[-1] - Lbeta_full[-ncond]
    
    ##construct LCL' matrix
    LCL <- matrix(0, ncond-1, ncond-1)
    L_mat <- matrix(0, ncond-1, ncond*p) 
    for (j in 1:(ncond-1)){
      L_mat[j, ((j-1)*p + 1): ((j+1)*p)] = c(-LN_list[[j]], LN_list[[j+1]]) 
      
      LCL[j,j] <- (s2e*lDtDlt[j] + s2g*llt[j]) + (s2e*lDtDlt[j+1] + s2g*llt[j+1])
      if (j<ncond-1){
        LCL[j,j+1] <- - (s2e*lDtDlt[j+1] + s2g*llt[j+1])
        LCL[j+1,j] <- LCL[j,j+1]
      }
      if (j>1){
        LCL[j-1,j] <- - (s2e*lDtDlt[j] + s2g*llt[j])
        LCL[j,j-1] <- LCL[j-1,j]
      } 
    }
    q = rankMatrix(L_mat)		
    
    ##get the test statistic
    teststat[rr] = (t(Lbeta) %*% solve(LCL) %*% Lbeta)/q
    
    ##find projection P and diagonal D such that LCL' = P'DP
    D_diag <- eigen(LCL)$values
    
    gm <- lapply(1:q, function(j) c(llt[j]+llt[j+1], lDtDlt[j] + lDtDlt[j+1]))
    vm <- sapply(1:q, function(j) (2*D_diag[j]^2) / (t(gm[[j]]) %*% KmatInv %*% gm[[j]] ))
    Em <- sum(sapply(vm, function(a) ifelse(a>2, a/(a-2), 0)))
    
    ##The first degree of freedom is q. The second degree of freedom is df. 
    if(Em > q){
      df[rr] = (2*Em) / (Em - q)
    }
    
    ## NOTE: matlab only accepts positive integers as df's for F-dist.
    ## Therefore, I had set the denom df to 1 if it is zero, otherwise round it to integer values
    df[rr] <- (df[rr] >= 1) * ceiling(df[rr]) + (df[rr] < 1) * 1
  }
  
  pvals = 1 - pf(abs(teststat), q, df) + pf(-abs(teststat), q, df)  
  
  return(list(teststat = teststat, df = df, p.value = pvals))
}
