netgsa.onesample <-
function(s2g, s2e, D, DDt, DDInv, n_vec, B, beta) {
  # Here D, DDt and DDInv are lists of block diagonal matrices. 
  npath = nrow(B)
  p = ncol(B)
  n1 = n_vec
  ncond = 1
  
  ##Initializing the vector for degrees of freedom for the test statistics. 
  teststat = matrix(0, npath, 1)
  num.tstat = matrix(0, npath, 1)
  df = matrix(0, npath, 1)
  
  ##----------------- 
  ##CALCULATING DEGREES OF FREEDOM & TEST STATS 
  ##----------------- 
  #matrices needed in calculatoin of degrees of freedom
  ## Sigma is W in the notes. 
  Sigma = lapply(DDt, function(ix) s2e * diag(1, nrow(ix)) + s2g * ix)
  SigmaInv = lapply(Sigma, function(ix) chol2inv(chol(ix)))
  SigmaInvD = mapply(function(a,b) crossprod(a,b), SigmaInv, DDt, SIMPLIFY = FALSE)
  SinvSinv = sapply(SigmaInv, function(ix) matTr(ix, ix))
  SinvDSinvD = sapply(SigmaInvD, function(ix) matTr(ix, ix))
  SinvSinvD = mapply(function(a,b) matTr(a,b), SigmaInv, SigmaInvD)
  EH11 = 0.5*sum(SinvDSinvD)
  EH12 = 0.5*sum(SinvSinvD)
  EH22 = 0.5*sum(SinvSinv)
  Kmat = matrix(c(EH11, EH12, EH12, EH22), 2, 2, byrow = TRUE)
  
  ## Kmat will be the expected information matrix. 
  ## Here we need its inverse in calculating the degrees of freedom. 
  KmatInv = solve(Kmat)
  
  ##Building the "contrast" matrix L, see Result in the paper 
  LN = crossprod(t(B), as.matrix(bdiag(D))) * B
  
  ##----------------- 
  for (rr in 1:npath) {
    Lrow1 = t(as.matrix(LN[rr, ])) #single row 

    g1 = (1/n1) * tcrossprod(Lrow1)
    g2 = (1/n1) * Lrow1 %*% as.matrix(bdiag(DDInv)) %*% t(Lrow1)
    g = matrix(c(g1, g2), 2, 1)
    LC11Lprime = s2g * g1 + s2e * g2
    
    #test statistic 
    num.tstat[rr] = crossprod(t(Lrow1), beta)
    teststat[rr] = num.tstat[rr]/sqrt(LC11Lprime)
    
    #calculating df based on the Satterthwaite approximation method  
    #using the formula nu=(2*LCL'^2)/g'Kg with K being the empirical covariance matrix.
    #NOTE: If df2<2, it is set to 2 
    
    df[rr] = 2 * (LC11Lprime)^2/(crossprod(g, KmatInv) %*% g)
    if (df[rr] < 2) 
      df[rr] = 2
    
  }
  pvals = 1 - pt(abs(teststat), df) + pt(-abs(teststat), df)
  
  return(list(teststat = teststat, df = df, p.value = pvals))
}
