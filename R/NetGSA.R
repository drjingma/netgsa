NetGSA <-
function(
    A,  #Adjacency matrix in a list	    	
    x, 	#The p x n data matrix		      	
    y,    #vector of class indicators of length n
    B,  	#indicator matrix for pathways (npath x p)	    	  	
    lklMethod = c("REML","ML"),
    directed = FALSE,          
    eta = 1e-1,           
    lim4kappa = 500       
  ){
    this.call <- match.call()
    lklMethod <- match.arg(lklMethod)
    
    p = dim(x)[1] #No. of genes
    n = length(y) #No. of samples in total
    
    if (dim(x)[2] != n) {
      stop("The dimensions of the data matrix and class vector don't match!")
    }
    
    if (dim(B)[2] != p) {
      stop("The dimensions of the data matrix and indicator matrix don't match!")
    }
    
    if (min(apply(B,1,sum))==0){
      stop("Empty pathways detected!")
    }
    
    if (length(unique(y)) < 2) {
      stop("There should be at least 2 unique classes in the class indicator!")
    }
    
    if (min(sapply(lapply(A, abs), sum))==0) {
      stop("No network interactions were found!")
    }
    
    ##-----------------
    ##Determine whether the network is DAG
    ##Assume A1 and A2 are of the same type (directed or undirected)
    isNetDAG = FALSE
    if (directed) {    
      gA <- graph.adjacency(A[[1]], mode="directed")
      isNetDAG = is.dag(gA)
    }  
    
    if (p > 5000) {
      warning("netGSA may be slow for datasets with large number of genes.")
    }
    
    ##-----------------
    ##setting up control parameters for the var estimation procedures
    varEstCntrl = list(lklMethod = lklMethod,                    
                       s2profile = "se",   
                       lb = 0.5,           
                       ub = 100,           
                       tol = 0.01)         
    
    ##-----------------
    ##Find influence matrices based on adjacency matrices A1 and A2
    ##Check if the influence matrices are well conditioned. Otherwise update eta.
    if (directed){
      D = lapply(A, function(m) adj2inf(AA=m, isDAG = isNetDAG, eta = eta))
      
      tmp = min(sapply(D, kappa)) 
      while ((tmp> lim4kappa) && !isNetDAG) {
        eta = eta * 2
        warning(paste("Influence matrix is ill-conditioned, using eta =", eta))
        D = lapply(A, function(m) adj2inf(AA=m, isDAG = isNetDAG, eta = eta))
        tmp = min(sapply(D, kappa)) 
      }
      
      DD = lapply(D, function(m) m %*% t(m))
      
      tmp = min(sapply(DD, kappa))    
      while ((tmp > lim4kappa) && !isNetDAG) {
        eta = eta * 2
        warning(paste("Influence matrix is ill-conditioned, using eta =", eta))  
        D = lapply(A, function(m) adj2inf(AA=m, isDAG = isNetDAG, eta = eta))
        DD = lapply(D, function(m) m %*% t(m))
        tmp = min(sapply(DD, kappa))
      }
      
    } else {
      #Undirected gaussian graphical model
      #Need to normalize the matrix differently
      Ip = diag( rep(1,p) )
      D = lapply(A, function(m) t(chol(solve(Ip - m))) ) 
    }
    
    output = call.netgsa(D, x, y, B, varEstCntrl)
    
    return(output)
  }
