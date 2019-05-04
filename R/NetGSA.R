NetGSA <-
  function(
    A,     	
    x, 	
    group,    
    pathways, 	
    lklMethod=c("REML","ML", "HE", "REHE"),
    sampling = FALSE,
    sample_n = NULL,
    sample_p = NULL, 
    minsize=5,
    eta=0.1,           
    lim4kappa=500
  ){
    this.call <- match.call()
    lklMethod <- match.arg(lklMethod)
    
    p <- dim(x)[1] #No. of genes
    n <- length(group) #No. of samples in total
    
    if (is.null(rownames(x))){
      stop('Data matrix must have rownames!')
    }
    
    if (is.null(rownames(A[[1]])) || is.null(rownames(A[[2]]))){
      stop('Adjacency matrix must have rownames!')
    }
    
    if (!identical(rownames(x),colnames(pathways))){
      stop('Genes in the data matrix and the pathway indicator matrix must be in the same order!')
    }
    
    if (p > 3000) {
      warning("netGSA may be slow for datasets with large number of genes.")
    }
    
    if (dim(x)[2] != n) {
      stop("The dimensions of the data matrix and class vector don't match!")
    }
    
    if (n<10){
      warning("The sample size is too small! Use NetGSA at your discretion!")
    }
    
    #If any of these are out of order, reorder
    outoforder <- vapply(A, function(Ai, data) { 
                          if(setequal(rownames(Ai), rownames(data)))   stop("Adjacency matrices and data do not contain same list of genes")
                          else if(all(rownames(Ai) == rownames(data))) return(FALSE)
                          else                                         return(TRUE)
                        }, FUN.VALUE = logical(1), data = x)
    if(any(outoforder)){
      order <- rownames(A[[1]])
      A <- lapply(A, function(Ai) Ai[order, order])
      x <- x[order,]
    }
    
    ##-----------------
    ##setting up control parameters for the var estimation procedures
    varEstCntrl = list(lklMethod = lklMethod,                    
                       s2profile = "se",
                       sampling = sampling,
                       ratio = sample_n,
                       p_sample = sample_p,
                       lb = 0.5,           
                       ub = 100,           
                       tol = 0.01)         
    
    A_c <- A
    if (min(sapply(lapply(A_c, abs), sum))==0) {
      warning("No network interactions were found! Check your networks!")
    }
    
    ##-----------------
    ##Determine whether the network is DAG
    ##Assume A_c[[1]] and A_c[[2]] are of the same type (directed or undirected)
    isNetDAG <- FALSE
    gA <- igraph::graph_from_adjacency_matrix((abs(A_c[[1]])>1e-06), mode="directed")
    isNetDAG <- igraph::is_dag(gA)
    
    p_c <- p
    x_c <- x[match(rownames(A_c[[1]]),rownames(x)),]
    
    ##Find influence matrices based on adjacency matrices in A_c
    ##Check if the influence matrices are well conditioned. Otherwise update eta.
    if (isNetDAG){
      D <- lapply(A_c, function(m) adj2inf(AA=m, isDAG = isNetDAG, eta = eta))
      tmp <- min(sapply(D, kappa)) 
      while ((tmp> lim4kappa) && !isNetDAG) {
        eta <- eta * 2
        warning(paste("Influence matrix is ill-conditioned, using eta =", eta))
        D <- lapply(A_c, function(m) adj2inf(AA=m, isDAG = isNetDAG, eta = eta))
        tmp <- min(sapply(D, kappa)) 
      }
      
      DD <- lapply(D, function(m) m %*% t(m))
      tmp <- min(sapply(DD, kappa))    
      while ((tmp > lim4kappa) && !isNetDAG) {
        eta <- eta * 2
        warning(paste("Influence matrix is ill-conditioned, using eta =", eta))  
        D <- lapply(A_c, function(m) adj2inf(AA=m, isDAG = isNetDAG, eta = eta))
        DD <- lapply(D, function(m) m %*% t(m))
        tmp <- min(sapply(DD, kappa))
      }
      
    } else {
      #Undirected gaussian graphical model
      Ip <- diag( rep(1,p_c) )
      D <- lapply(A_c, function(m) t(chol(pseudoinverse(Ip - m)))) 
    }
    output <- call.netgsa(D, x_c, group, pathways, varEstCntrl)
    
    ## Update the format of the output to be consistent with other methods. 
    out <- data.frame('pathway'= rownames(pathways),
                      'pSize' = rowSums(pathways),
                      'pval' = output$p.value,
                      'pFdr' = p.adjust(output$p.value,"BH"))
    
    out <- out[order(out$pFdr),]
    rownames(out) <- NULL
    
    return(list(results=out,beta=output$beta,s2.epsilon=output$s2.epsilon,s2.gamma=output$s2.gamma))
  }
