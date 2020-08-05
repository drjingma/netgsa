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
  
  ncond = length(D)
  p = nrow(resid[[1]])
  
  n_k = sapply(resid, ncol)


  if (!sampling)
  {
    ##------------------
    ## full-data estimation
    ##-------------------------------------
    I_p = diag(1, p)
    mat_index = lower.tri(I_p, diag = T)
    XTX = matrix(0, 2, 2)
    XTX[1,1] = n*p
    XTX[1, 2]<-XTX[2, 1]<- sum(sapply(1:ncond, function(k)n_k[k]*sum(diag(DtD[[k]]))))
    XTX[2, 2] = sum(sapply(1:ncond, function(k) n_k[k]*sum(DtD[[k]][mat_index]^2)))
    
    
    
    cond_k = function(k){
      nk=n_k[k]
      
      tmp = sapply(1:nk, function(i){
        resid_sq = resid[[k]][,i, drop=F] %*Cpp% t(resid[[k]][,i, drop=F])
        XY_1 = sum(diag(resid_sq))
        XY_2 = sum((resid_sq*DtD[[k]])[mat_index])
        return(c(XY_1, XY_2))
      })

      return(rowSums(tmp))
    }
    XTY = matrix(rowSums(sapply(1:ncond, cond_k)), c(2,1))
    
    if (lklMethod == 'HE'){
      res = solveCpp(XTX) %*Cpp% XTY
      s2e = ifelse(res[1]<0, 0, res[1])
      s2g = ifelse(res[-1]<0, 0, res[-1])
    } else if (lklMethod == 'REHE'){
      res_qp = solve.QP(Dmat = XTX, dvec=XTY, Amat = diag(1,2), bvec = rep(0, 2), meq=0, factorized=FALSE)
      res_qp$solution[res_qp$solution<0]<-0
      s2e = res_qp$solution[1]
      s2g = res_qp$solution[-1]
    }
    
  }
  
  # for subsampling, by default is repeated subsampling of 50 replications
  if (sampling){
    ## uniform sample observations and genes, but seperate sampleof gene for each subject
    #if(!is.null(control$sample_seeds)){set.seed(control$sample_seed[1])} ##MH Added Sample seed - REMOVED useless b/c cannot set seed for "p"

    # condition K:
    cond_k_sample = function(k){
      nk = n_k[k]
      nk_index = sample(1:nk, size = round(ratio*nk, 0), replace = T)
      nk_s = length(nk_index)
      resid_k_s = resid[[k]][, nk_index]
      return(resid_k_s)
    }
    
    one_subsample =function(rep_subsample){
      resid_s = lapply(1:ncond, cond_k_sample)
      temp = do.call(cbind, resid_s)
      n_s_k = sapply( resid_s, ncol)
      n_s = dim(temp)[2]
      #weight_normp = sqrt(colSums((t(temp)-colMeans(temp))^2))  
      weight_normp = rep(1/p, p)
      
      ## for repeated subsampling, cannot set a same seed to different replications
      #if(!is.null(control$sample_seeds)){set.seed(control$sample_seeds[2])} ##MH Added Sample seed
      mat_index = lower.tri(diag(0, round(p*(p_sample), 0)), diag=T)
      cond_k_stat = lapply(1:ncond, function(k) 
        sapply(1:n_s_k[k], function(i){
          index_p = sample(1:p, size = round(p*(p_sample), 0), replace=T, prob = weight_normp)
          resid_sq = resid_s[[k]][index_p,i, drop=F] %*Cpp% t(resid_s[[k]][index_p,i, drop=F])
          DTD_tmp = DtD[[k]][index_p, index_p]
          XY_1 = sum(diag(resid_sq))
          XY_2 = sum((resid_sq*DTD_tmp)[mat_index])
          XX12 = sum(diag(DTD_tmp))
          XX22 = sum((DTD_tmp^2)[mat_index])
          return(c(XY_1, XY_2, XX12, XX22))
        })
      )
      
      tmp = rowSums(do.call(cbind, cond_k_stat))
      
      XTX = matrix(0, 2, 2)
      XTX[1,1] <- round(p*(p_sample), 0)*n_s
      
      XTX[1,2] <-XTX[2,1] <-tmp[3]
      
      XTX[2,2] <- tmp[4]
      
      XTY = tmp[1:2]
      
      
      if (lklMethod == 'HE'){
        res = solveCpp(XTX) %*Cpp% XTY
        s2e = ifelse(res[1]<0, 0, res[1])
        s2g = ifelse(res[-1]<0, 0, res[-1])
      } else if (lklMethod == 'REHE'){
        res_qp = solve.QP(Dmat = XTX, dvec=XTY, Amat = diag(1,2), bvec = rep(0, 2), meq=0, factorized=FALSE)
        res_qp$solution[res_qp$solution<0]<-0
        s2e = res_qp$solution[1]
        s2g = res_qp$solution[-1]
      }
      
      return(c(s2e = s2e, s2g = s2g)) 
      
    }
    
    subsample_res = sapply(1:50, one_subsample)
    subsample_res_summary = apply(subsample_res, FUN=mean, 1)
    s2e = subsample_res_summary[1]
    s2g = subsample_res_summary[2]
    }
    
    
    

  
  return(list(s2e = s2e, s2g = s2g, test='NULL')) 
}
