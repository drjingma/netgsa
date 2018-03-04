netgsa.onesample <-
function(s2g, s2e, D, DtD, DDInv, n_vec, B, beta) {
  ncond = length(D)
  npath = dim(B)[1]
  p = dim(B)[2]
  Ip = matrix(0, p, p)
  diag(Ip) = 1
  
  n1 = n_vec[1]
  teststat = matrix(0, npath, 1)
  num.tstat = matrix(0, npath, 1)
  
  ##Initializing the vector for degrees of freedom for the test statistics. 
  df = matrix(0, npath, 1)
  
  ##Building the "contrast" matrix L, see Result in the paper 
  LN = (B %*% D[[1]]) * B
  
  ##----------------- 
  ##CALCULATING DEGREES OF FREEDOM & TEST STATS 
  ##----------------- 
  #matrices needed in calculatoin of degrees of freedom
  Sigma = lapply(1:ncond, function(ix) s2e * Ip + s2g * DtD[[ix]])
  SigmaInv = lapply(1:ncond, function(ix) chol2inv(chol(Sigma[[ix]])))
  SigmaInvD = lapply(1:ncond, function(ix) SigmaInv[[ix]] %*% DtD[[ix]])
  SinvSinv = lapply(1:ncond, function(ix) matTr(SigmaInv[[ix]] %*% SigmaInv[[ix]]))
  SinvDSinvD = lapply(1:ncond, function(ix) matTr(SigmaInvD[[ix]] %*% SigmaInvD[[ix]]))
  SinvSinvD = lapply(1:ncond, function(ix) matTr(SigmaInv[[ix]] %*% SigmaInvD[[ix]]))
  EH11 = (1/2) * Reduce("+",SinvDSinvD)
  EH12 = (1/2) * Reduce("+",SinvSinvD)
  EH22 = (1/2) * Reduce("+",SinvSinv)
  
  #In this version of the code, K matrix is calculated directly! 
  ## Kmat will be the expected information matrix. Here we need its inverse  
  ## in calculating the degrees of freedom. 
  Kmat = matrix(c(EH11, EH12, EH12, EH22), 2, 2, byrow = TRUE)
  KmatInv = solve(Kmat)
  
  #These matrices are needed in the calculation of test statistics 
  mctildi = s2e * DDInv[[1]] + s2g * Ip

  ##----------------- 
  for (rr in 1:npath) {
    Lrow = t(as.matrix(LN[rr, ])) #single row of L3 
    Lrow1 = t(as.matrix(Lrow[, 1:p]))

    LC11Lprime = (1/n1) * Lrow1 %*% mctildi %*% t(Lrow1)
    
    g1 = (1/n1) * (Lrow1 %*% t(Lrow1))
    g2 = (1/n1) * Lrow1 %*% DDInv[[1]] %*% t(Lrow1)
    g = matrix(c(g1, g2), 2, 1)
    
    #test statistic 
    num.tstat[rr] = Lrow1 %*% beta[[1]]
    teststat[rr] = num.tstat[rr]/sqrt(LC11Lprime)
    
    #calculating df based on the Satterthwaite approximation method  
    #using the formula nu=(2*LCL'^2)/g'Kg with K being the empirical covariance matrix.
    #NOTE: If df2<2, it is set to 2 
    
    df[rr] = 2 * (LC11Lprime)^2/(t(g) %*% KmatInv %*% g)
    if (df[rr] < 2) 
      df[rr] = 2
    
  }
  pvals = 1 - pt(abs(teststat), df) + pt(-abs(teststat), df)
  
  return(list(teststat = teststat, df = df, p.value = pvals))
}
