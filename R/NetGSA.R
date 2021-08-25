NetGSA <-
  function(
    A,     # a list of adj matrices
    x, 	
    group,    
    pathways, 	
    lklMethod="REHE",
    sampling = FALSE,
    sample_n = NULL,
    sample_p = NULL, 
    minsize=5,
    eta=0.1,           
    lim4kappa=500
  ){
    this.call <- match.call()
    lklMethod <- match.arg(lklMethod, c("REML","ML", "HE", "REHE"))
    
    edgelist <- A[["edgelist"]]
    A[["edgelist"]] <- NULL
    
    p <- dim(x)[1] #No. of genes
    n <- length(group) #No. of samples in total
    
    # Assume the adj matrix is always block diagonal; It has one block in the worst case. 
    # Also add dimnames onto matrix. This automatically works with clustering or no clustering
    A_mat <- lapply(A, function(a) {a_full <- as.matrix(bdiag(a))
                                    dimnames(a_full) <- list(do.call(c, lapply(a, rownames)), do.call(c, lapply(a, colnames)))
                                    return(a_full)
                                    })
    
    #If any of these are out of order, reorder
    outoforder <- vapply(A_mat, function(Ai, data) { 
                            if(!setequal(rownames(Ai), rownames(data)))  stop("Adjacency matrices and data do not contain same list of genes")
                            else if(all(rownames(Ai) == rownames(data))) return(FALSE)
                            else                                         return(TRUE)
                          }, FUN.VALUE = logical(1), data = x)
    if(any(outoforder)){
      order <- rownames(A_mat[[1]])
      A_mat <- lapply(A_mat, function(Ai) Ai[order, order])
      x <- x[order,]
      pathways <- pathways[,order, drop = FALSE]
    }
    
    
    if (max(sapply(lapply(A_mat,abs),sum))==0) {
      warning("No network interactions were found! Check your networks!")
    }
   
    if (is.null(rownames(x))){
      stop('Data matrix must have rownames!')
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
    ##-----------------
    ##setting up control parameters for the var estimation procedures
    varEstCntrl = list(lklMethod = lklMethod,                    
                       s2profile = "se",
                       sampling = sampling,
                       ratio = sample_n,
                       p_sample = sample_p,
                       lb = 0.5,           
                       ub = 100,           
                       tol = 0.01,
                       maxIter = 100)         
    
    ##-----------------
    ##Determine whether the network is DAG
    isNetDAG <- FALSE
    gA <- igraph::graph_from_adjacency_matrix((abs(A_mat[[1]])>1e-06), mode="directed")
    isNetDAG <- igraph::is_dag(gA)
    
    ##Find influence matrices based on adjacency matrices in A
    ##Check if the influence matrices are well conditioned. Otherwise update eta.
    if (isNetDAG){
      D <- lapply(A, lapply, function(m) adj2inf(AA=m, isDAG = isNetDAG, eta = eta))
      tmp <- min(sapply(D,function(m) kappa(as.matrix(bdiag(m)))))
      while ((tmp > lim4kappa) && !isNetDAG) {
        eta <- eta * 2
        warning(paste("Influence matrix is ill-conditioned, using eta =", eta))
        D <- lapply(A, lapply, function(m) adj2inf(AA=m, isDAG = isNetDAG, eta = eta))
        tmp <- min(sapply(D,function(m) kappa(as.matrix(bdiag(m)))))
      }
      
      DD <- lapply(D, lapply, tcrossprod)
      tmp <- min(sapply(DD,function(m) kappa(as.matrix(bdiag(m)))))
      while ((tmp > lim4kappa) && !isNetDAG) {
        eta <- eta * 2
        warning(paste("Influence matrix is ill-conditioned, using eta =", eta))  
        D <- lapply(A, lapply, function(m) adj2inf(AA=m, isDAG = isNetDAG, eta = eta))
        DD <- lapply(D, lapply, tcrossprod)
        tmp <- min(sapply(DD,function(m) kappa(as.matrix(bdiag(m)))))
      }
      
    } else {
      #Undirected gaussian graphical model
      D <- lapply(A, lapply, function(b) t(cholCpp(pseudoinverse(diag(1,nrow(b)) - b)))) 
    }
    
    output <- call.netgsa(D, x, group, pathways, varEstCntrl)
    
    ## Update the format of the output to be consistent with other methods. 
    out <- data.frame('pathway'= rownames(pathways),
                      'pSize' = rowSums(pathways),
                      'pval' = output$p.value,
                      'pFdr' = p.adjust(output$p.value,"BH"),
                      'teststat' = output$teststat)
    
    out <- out[order(out$pFdr),]
    rownames(out) <- NULL
    #Test individual genes for plotting
    if(length(unique(group)) > 1){
      gene_tests <- setNames(genefilter::rowFtests(x, factor(group)), c("teststat", "pval"))
    } else{
      gene_tests <- setNames(genefilter::rowttests(x)[,c("statistic", "p.value")], c("teststat", "pval"))
    }
    #Copy b/c setkeyv then orders the rownames of x, not just for gene_tests. This is some quite odd data.table behavior
    gene_tests[, c("gene", "pFdr")] <- list(copy(rownames(x)), p.adjust(gene_tests[["pval"]], "BH"))
    setDT(gene_tests); setkeyv(gene_tests, "gene")
    
    #Graph object for plotting
    graph_obj <- list(edgelist = edgelist, pathways = reshapePathways(pathways), gene.tests = gene_tests)
    netgsa_obj <- list(results=out,beta=output$beta,s2.epsilon=output$s2.epsilon,s2.gamma=output$s2.gamma, graph = graph_obj)
    class(netgsa_obj) <- "NetGSA"
    
    return(netgsa_obj)
  }



# Helper functions -----------------------------------------------------------------

reshapePathways <- function(pathways){
  . <- pathway <- gene <- NULL #Added to avoid data.table note in R CMD check
  res <- data.table::setDT(reshape2::melt(pathways, varnames = c("pathway", "gene")))
  res[, c("pathway", "gene") := .(as.character(pathway), as.character(gene))]
  return(res[res$value == 1,][,c("pathway", "gene")])
}

`%*Cpp%` <- function(x, y){
  matMult(x,y)
}
