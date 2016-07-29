approxVarEst <-
function(se0, sg0, D, r, n_vec, control=NULL){
    #D and r are both lists. 
    if (is.null(control)) {
      tol = 0.01
      s2profile = "se"
      lklMethod = "REML"
    } else {
      tol = control$tol
      s2profile = control$s2profile
      lklMethod = control$lklMethod
    }
    
    matTr <- function(z) sum(diag(z))
    
    p = nrow(D[[1]])
    Ip = diag(rep(1, p))
    ncond = length(D)
    
    n = sum(n_vec)
    N = n * p
    
    ##residual matrices
    Rm = vector("list", ncond)
    for (loop_cond in 1:ncond){
      Rm[[loop_cond]] = matrix(0, p, p)
      tmp = r[[loop_cond]]
      for (i in 1:n_vec[loop_cond]){
        Rm[[loop_cond]] = Rm[[loop_cond]] +  tmp[, i] %o% tmp[, i]
      }
    }
    DtD = lapply(1:ncond, function(k) D[[k]] %*% t(D[[k]]))
    
    gap = 1
    sg = sg0
    se = se0
    cnt = 0
    ## Whether it's ML or REML
    lklConst = ifelse(lklMethod == "REML", N - ncond*p, N)
    
    while (gap > tol) {
      sg0 = sg
      se0 = se
      
      # print(cnt)
      cnt = cnt + 1
      
      SS = lapply(1:ncond, function(k) solve(sg * DtD[[k]] + se*Ip) %*% Rm[[k]] )
      tmp0 = lapply(SS, matTr)
      tmp0 = Reduce("+", tmp0)
      
      if (s2profile == "se") {     
        se =  tmp0 / lklConst
        tmp2 = vector("list", ncond)
        for (loop_cond in 1:ncond){
          tmp1 = cov(t(r[[loop_cond]]))
          tmp = min(diag(tmp1))
          tmp = min(se, tmp)
          tmp1 = tmp1 - diag(rep(tmp, p))
          tmp2[[loop_cond]] = solve(D[[loop_cond]]) %*% tmp1  
        }
        
        tmp3 = lapply(tmp2, function(A) mean(diag(A)))      
        sg = Reduce("+", tmp3)/ncond
        tau = sg/se     
      } else {     
        sg = tmp0/lklConst
        tmp2 = vector("list", ncond)
        for (loop_cond in 1:ncond){
          tmp1 = cov(t(r[[loop_cond]]))
          tmp = min(diag(tmp1))
          tmp = min(sg, tmp)
          tmp1 = tmp1 - diag(rep(tmp, p))
          tmp2[[loop_cond]] = solve(D[[loop_cond]]) %*% tmp1  
        }
        
        tmp3 = lapply(tmp2, function(A) mean(diag(A)))  
        se = Reduce("+", tmp3)/ncond
        tau = se/sg
      }
      
      gap = abs(sg - sg0) + abs(se - se0)
    }
    
    return(list(s2e = se, s2g = sg, tau = tau, finalgap = gap, niter = cnt))
  }
