he_method = function(n, D, DtD, resid, control=NULL){
  #resid / D / DtD are all length K lists, resid_m is pxn matrix
  if (is.null(control)) {
    sampling = T
    lklMethod = "REHE"
    p_sample = 0.75
    ratio = 0.2
  } else {
    sampling = ifelse(is.null(control$sampling), TRUE, control$sampling)
    lklMethod = control$lklMethod
    p_sample = ifelse(is.null(control$p_sample), 0.75, control$p_sample)  
    ratio = ifelse(is.null(control$ratio), 0.2, control$ratio)
  }
  
  resid_m = do.call(cbind, resid)
  p = dim(resid_m)[1]
  ncond = length(D)

  rm(resid_m)
  
  if (!sampling)
  {
    ##------------------
    ## full-data estimation
    ##---------------------------------------
    XTX = matrix(0, 2, 2)
    XTX[1,1] = n*p
    XTX[1, 2]<-XTX[2, 1]<- sum(sapply(1:ncond, function(k)dim(resid[[k]])[2]*sum(diag(DtD[[k]]))))
    XTX[2,2] = sum(sapply(1:ncond, function(k) dim(resid[[k]])[2]*sum((DtD[[k]][lower.tri(DtD[[k]], diag = T)])^2)))
    
    s0 = matrix(unlist(sapply(1:p, function(i){c(1, rep(0, p-i))})), ncol=1)
    cond_k = function(k){
      nk=dim(resid[[k]])[2]
      
      si = matrix(DtD[[k]][lower.tri(DtD[[k]], diag = T)], ncol=1)
      resid_value = function(i){
        resid_sq = resid[[k]][,i, drop=F] %*% t(resid[[k]][,i, drop=F])
        return(resid_sq[lower.tri(resid_sq, diag = T)])
      }
      y_sum = matrix(rowSums(sapply(1:nk, resid_value)), ncol=1)
      return(c(t(s0)%*%y_sum, t(si)%*%y_sum))
    }
    XTY = rowSums(sapply(1:ncond, cond_k))
    
    if (lklMethod == 'HE'){
      res = solve(XTX) %*% XTY
      s2e = ifelse(res[1]<0, 0, res[1])
      s2g = ifelse(res[-1]<0, 0, res[-1])
    } else if (lklMethod == 'REHE'){
      res_qp = solve.QP(Dmat = XTX, dvec=XTY, Amat = diag(1,ncond), bvec = rep(0, ncond), meq=0, factorized=FALSE)
      s2e = res_qp$solution[1]
      s2g = res_qp$solution[-1]
    }
    
  }
  
  if (sampling){
    ## uniform sample 20% of observations; then do weighted sampling by norm for p genes
    
    # condition K:
    cond_k_sample = function(k){
      nk = dim(resid[[k]])[2]
      nk_index = sample(1:nk, size = round(ratio*nk, 0), replace = F)
      nk_s = length(nk_index)
      resid_k_s = resid[[k]][, nk_index]
      return(resid_k_s)
    }
    
    resid_s = lapply(1:ncond, cond_k_sample)
    temp = do.call(cbind, resid_s)
    n_s = dim(temp)[2]
    
    weight_normp = sqrt(colSums((t(temp)-colMeans(temp))^2))  
    index_p = sample(1:p, size = round(p*(p_sample), 0), replace=F, prob = weight_normp)
    resid_s = lapply(1:ncond, function(k) resid_s[[k]][index_p,])
    p_s = length(index_p)
    rm(temp)
    
    resid_value_ps = function(k, i){
      resid_sq = resid_s[[k]][,i, drop=F] %*% t(resid_s[[k]][,i, drop=F])
      return(resid_sq[lower.tri(resid_sq, diag = T)])                 
    }
    
    
    #D1D1 need modified     
    #D_s = lapply(1:ncond, function(k) D[[k]][index_p, index_p])
    DtD_s = lapply(1:ncond, function(k)  DtD[[k]][index_p, index_p])
    
    ###
    XTX = matrix(0, 2, 2)
    XTX[1,1] = n_s*p_s
    XTX[1, 2]<-XTX[2, 1]<- sum(sapply(1:ncond, function(k) dim(resid_s[[k]])[2]*sum(diag(DtD_s[[k]]))))
    XTX[2,2] = sum(sapply(1:ncond, function(k) dim(resid_s[[k]])[2]*sum((DtD_s[[k]][lower.tri(DtD_s[[k]], diag = T)])^2)))
    
    s0 = matrix(unlist(sapply(1:p_s, function(i){c(1, rep(0, p_s-i))})), ncol=1)
    cond_k = function(k){
      nk = dim(resid_s[[k]])[2]
      si = matrix(DtD_s[[k]][lower.tri(DtD_s[[k]], diag = T)], ncol=1)
      resid_value = function(i){
        resid_sq = resid_s[[k]][,i, drop=F] %*% t(resid_s[[k]][,i, drop=F])
        return(resid_sq[lower.tri(resid_sq, diag = T)])
      }
      y_sum = matrix(rowSums(sapply(1:nk, resid_value)), ncol=1)
      return(c(t(s0)%*%y_sum, t(si)%*%y_sum))
    }
    XTY = rowSums(sapply(1:ncond, cond_k))
    
    if (lklMethod == 'HE'){
      res = solve(XTX) %*% XTY
      s2e = ifelse(res[1]<0, 0, res[1])
      s2g = ifelse(res[-1]<0, 0, res[-1])
    } else if (lklMethod == 'REHE'){
      res_qp = solve.QP(Dmat = XTX, dvec=XTY, Amat = diag(1,ncond), bvec = rep(0, ncond), meq=0, factorized=FALSE)
      s2e = res_qp$solution[1]
      s2g = res_qp$solution[-1]
    }
  }
  
  return(list(s2e = s2e, s2g = s2g)) 
}
