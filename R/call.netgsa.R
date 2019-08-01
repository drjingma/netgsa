call.netgsa <-
  function(
    D,    	        #Influence matrices in a list (each list is a list of block digonal matrices)
    x, 			      	#the p x n data matrix
    y, 			      	#vector of class indicators of length n
    B, 		    	  	#indicator matix for pathways (npath x p)
    varEstCntrl     #parameters to pass for variance estimation
  ){
    p = nrow(x)     #No. of genes
    n = length(y)   #No. of samples in total
    npath = nrow(B) #No. of gene sets to consider

    y = as.integer(as.factor(y))
    ncond = length(unique(y))
    n_vec = as.numeric(table(y))

    xx = vector("list", ncond)  
    for (i in 1:ncond){
      xx[[i]] = x[, (y==i)]
    }
    
    DInv = lapply(D, lapply, solve)
    DDt = lapply(D, lapply, tcrossprod) # D %*% t(D)
    DtD = lapply(D, lapply, crossprod) # t(D) %*% D
    DtDInv = lapply(DtD, lapply, solve)
    
    #Initialzing pvalues of each test
    pvals = matrix(0, npath, 1)
    
    ##------------------
    ##ESTIMATION OF BETA
    ##------------------
    beta = vector("list", ncond)
    for (i in 1:ncond){
      beta[[i]] = (1/n_vec[i]) * as.matrix(bdiag(DInv[[i]])) %*% rowSums(xx[[i]])
    }
    
    ##-----------------
    ##ESTIMATION OF VAR COMPONENTS
    ##-----------------
    #First calculate residuals: 
    resid = vector("list", ncond)
    for (i in 1:ncond){
      resid[[i]] = xx[[i]] - replicate(n_vec[i], rowMeans(xx[[i]]))
    }
    temp = do.call(cbind, resid)
    
    #initial estimates
    sg0 = (mean(apply(temp, 2, sd)))^2
    se0 = (sd(diag(temp)))^2
    
    ##Estimation of variance components can be based on Restricted Haseman-Elston (REHE), or 
    ## simple Haseman-Elston (HE), and HE/REHE can use full data method or sampling method;
    ## or profile likeelihood with Newton's method using either ML or REML
    if (p<5000){
      D_list <- lapply(D, function(m) as.matrix(bdiag(m)))
      if (varEstCntrl$lklMethod %in% c('REML', 'ML')){
        approxEst = approxVarEst(se0, sg0, D_list, resid, n_vec, control = varEstCntrl) 
        se0 = approxEst$s2e
        sg0 = approxEst$s2g
        if (se0 < 1e-08 || sg0 < 1e-08){
          sg0 = (mean(apply(temp, 2, sd)))^2
          se0 = (sd(diag(temp)))^2
        }
        S = profile.newton.se(sg0/se0, D_list, resid, control = varEstCntrl)
        s2e = S$s2e
        s2g = S$s2g  	
      }
      if (varEstCntrl$lklMethod %in% c('HE', 'REHE')){
        DDt_list <- lapply(DDt, function(m) as.matrix(bdiag(m)))
        S = he_method(n, D_list, DDt_list, resid, control = varEstCntrl )
        s2e = S$s2e
        s2g = S$s2g
      } 
    }  else {
      s2e = se0
      s2g = sg0  	
    } 
    
    
    if (ncond==1){
      ##Use the one-sample T test
      res = netgsa.onesample(s2g, s2e, D[[1]], DDt[[1]], DtDInv[[1]], n_vec[1], B, beta[[1]])
      output = list(beta = beta, teststat = res$teststat, df = res$df, p.value = res$p.value, s2.epsilon = s2e, s2.gamma = s2g)      
    } else if (ncond==2){
      ##Use the two-sample T test
      res = netgsa.ttest(s2g, s2e, D, DDt, DtDInv, n_vec, B, beta)
      output = list(beta = beta, teststat = res$teststat, df = res$df, p.value = res$p.value, s2.epsilon = s2e, s2.gamma = s2g)      
    } else {
      ##Use the F test
      res = netgsa.ftest(s2g, s2e, D, DDt, DtDInv, n_vec, B, beta)
      output = list(beta = beta, teststat = res$teststat, df = res$df, p.value = res$p.value, s2.epsilon = s2e, s2.gamma = s2g) 
    } 
    
    return(output)
  }
