call.netgsa <-
function(
    D,    	        #Influence matrices in a list
    x, 			      	#the p x n data matrix
    y, 			      	#vector of class indicators of length n
    B, 		    	  	#indicator matix for pathways (npath x p)
    varEstCntrl     #parameters to pass for variance estimation
  ){
    p = nrow(x) #No. of genes
    n = length(y) #No. of samples in total
    npath = nrow(B) #No. of gene sets to consider
    
    y = as.integer(as.factor(y))
    ncond = length(unique(y))
    n_vec = as.numeric(table(y))
    
    ##The identity matrix 
    Ip = matrix(0, p, p)
    diag(Ip) = 1
    
    xx = vector("list", ncond)  
    for (i in 1:ncond){
      xx[[i]] = x[, (y==i)]
    }
    
    DInv = lapply(D, solve)
    DtD = lapply(1:ncond, function(k) D[[k]] %*% t(D[[k]]))
    tDD = lapply(1:ncond, function(k) t(D[[k]]) %*% D[[k]]) 
    DDInv = lapply(tDD, solve)
    
    #Initialzing pvalues of each test
    pvals = matrix(0, npath, 1)
    
    ##------------------
    ##ESTIMATION OF BETA
    ##------------------
    beta = vector("list", ncond)
    for (i in 1:ncond){
      beta[[i]] = (1/n_vec[i]) * DInv[[i]] %*% apply(xx[[i]], 1, sum)
    }
    
    ##-----------------
    ##ESTIMATION OF VAR COMPONENTS
    ##-----------------
    #First calculate residuals: 
    resid = vector("list", ncond)
    for (i in 1:ncond){
      resid[[i]] = xx[[i]] - replicate(n_vec[i], c(D[[i]] %*% beta[[i]]))
    }
    temp = do.call(cbind, resid)
    
    #initial estimates
    sg0 = (mean(apply(temp, 2, sd)))^2
    se0 = (sd(diag(temp)))^2
    
    ##Estimation of variance components is based on profile likelihood with Newton's method.
    ##It can be done using either ML or REML.
    if (p<5000){
      approxEst = approxVarEst(se0, sg0, D, resid, n_vec, control = varEstCntrl) 
      se0 = approxEst$s2e
      sg0 = approxEst$s2g
      S = profile.newton.se(sg0/se0, D, resid, control = varEstCntrl)
      s2e = S$s2e
      s2g = S$s2g  	
    } else {
      s2e = se0
      s2g = sg0  	
    } 
    
    if (ncond==2){
      ##Use the T test
      res = netgsa.ttest(s2g, s2e, D, DtD, DDInv,n_vec, B, beta)
      output = list(beta = beta, teststat = res$teststat, df = res$df, p.value = res$p.value, s2.epsilon = s2e, s2.gamma = s2g)      
    } else {
      ##Use the F test
      res = netgsa.ftest(s2g, s2e, D, DtD, DDInv, n_vec, B, beta)
      output = list(beta = beta, teststat = res$teststat, df = res$df, p.value = res$p.value, s2.epsilon = s2e, s2.gamma = s2g) 
    } 
    
    return(output)
  }
