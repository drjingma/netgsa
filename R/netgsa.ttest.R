netgsa.ttest <-
function(s2g, s2e, D, DDt, DtDInv, n_vec, B, beta) {
  ncond = length(n_vec)
  p = ncol(B)
  n1 = n_vec[1]
  n2 = n_vec[2]

  ##Initializing test statistics, its numerator and degrees of freedom. 
  npath = nrow(B)
  teststat = matrix(0, npath, 1)
  num.tstat = matrix(0, npath, 1)
  df = matrix(0, npath, 1)
  
  ## There matrices are needed in calculating DOF
  ## Sigma is W in the notes. 
  Sigma = lapply(DDt, lapply, function(ix) s2e * diag(1, nrow(ix)) + s2g * ix)
  SigmaInv = lapply(Sigma, lapply, function(ix) chol2inv(cholCpp(ix)))
  SigmaInvD = lapply(1:ncond, function(j) mapply(function(a,b) crossprodCpp(a,b), SigmaInv[[j]], DDt[[j]], SIMPLIFY = FALSE))
  SinvSinv = lapply(SigmaInv, sapply, function(ix) matTr(ix, ix))
  SinvDSinvD = lapply(SigmaInvD, sapply, function(ix) matTr(ix, ix))
  SinvSinvD = lapply(1:ncond, function(j) mapply(function(a,b) matTr(a,b), SigmaInv[[j]], SigmaInvD[[j]]))
  EH11 = 0.5*sum(sapply(SinvDSinvD,sum)) #This accommodates unequal number of blocks across conditions.
  EH12 = 0.5*sum(sapply(SinvSinvD,sum))
  EH22 = 0.5*sum(sapply(SinvSinv,sum))
  Kmat = matrix(c(EH11, EH12, EH12, EH22), 2, 2, byrow = TRUE)
  ## The Kmat will be the expected information matrix. 
  ## Here we need its inverse in calculating the degrees of freedom. 
  KmatInv = solveCpp(Kmat)
  
  ##Building the "contrast" vector; see Result in the paper 
  LN = lapply(1:ncond, function(j) crossprodCpp(t(B), as.matrix(bdiag(D[[j]]))) * B)
  ##----------------- 
  ##CALCULATING DEGREES OF FREEDOM & TEST STATS 
  ##----------------- 
  
  for (rr in 1:npath) {
    llt <- sapply(1:ncond,function(j) crossprodCpp(LN[[j]][rr,])/n_vec[j]) 
    lDtDlt <- sapply(1:ncond ,function(j) t(LN[[j]][rr,])%*Cpp%as.matrix(bdiag(DtDInv[[j]]))%*Cpp%LN[[j]][rr,]/n_vec[j])
    Lbeta_full <- sapply(1:ncond, function(j) crossprodCpp(LN[[j]][rr,], beta[[j]]))
    
    g = matrix(c(sum(llt), sum(lDtDlt)), 2, 1)
    LC11Lprime = s2g * g[1] + s2e * g[2]

    #test statistic 
    num.tstat[rr] =  Lbeta_full[2] - Lbeta_full[1]
    teststat[rr] = num.tstat[rr]/sqrt(LC11Lprime)
    
    #calculating df based on the Satterthwaite approximation method  
    #using the formula nu=(2*LCL'^2)/g'Kg with K being the empirical covariance matrix.
    #NOTE: If df2<2, it is set to 2 
    df[rr] = 2 * LC11Lprime^2/(t(g) %*Cpp% KmatInv %*Cpp% g)
    if (df[rr] < 2) 
      df[rr] = 2
    
  }
  pvals = 1 - pt(abs(teststat), df) + pt(-abs(teststat), df)
  
  return(list(teststat = teststat, df = df, p.value = pvals))
}
