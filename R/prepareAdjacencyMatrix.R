prepareAdjacencyMatrix <-
  function(x,             
           group,         
           pathways,      
           import_from_kegg=FALSE,
           file_e=c(NA, file_e), 
           file_ne=c(NULL, file_ne, NA),
           estimate_network=FALSE,
           lambda_c=1,
           eta=0.5,
           minsize=5,
           fileEncoding=""
  ) {
    this.call <- match.call()

    if (!is.na(file_e)){
      if (is.character(file_e)) {
        file_e <- if (nzchar(fileEncoding))
          file(file_e, "rt", encoding = fileEncoding)
        else file(file_e, "rt")
        on.exit(close(file_e))
      } 
    }
    if (!is.null(file_ne) && !is.na(file_ne)){
      if (is.character(file_ne)) {
        file_ne <- if (nzchar(fileEncoding))
          file(file_ne, "rt", encoding = fileEncoding)
        else file(file_ne, "rt")
        on.exit(close(file_ne))
      }
    }
    
    ## first match variables in each pathway with variables in the data
    pathways <- lapply(pathways, function(a) intersect(a,rownames(x)))
    if (min(sapply(pathways, length))<minsize){
      pathways <- pathways[-which(sapply(pathways, length)<=minsize)]      
    }
    npath <- length(pathways)
    ncond <- length(unique(group))
    n <- as.numeric(table(group))
    p <- nrow(x)
    
    ## make sure the pathways are organized as a matrix
    listB <- lapply(1:npath, function(a) as.numeric(!is.na(match(rownames(x),pathways[[a]]))))
    B <- do.call(rbind,listB)
    colnames(B) <- rownames(x)
    rownames(B) <- names(pathways)
    
    if (is.na(file_e)){
      if (!import_from_kegg){
        warning('No edges provided!')
        Adj <- NULL
        g <- igraph::make_empty_graph(p, directed = FALSE)
      } 
      if (import_from_kegg){
        paths <- graphite::pathways('hsapiens','kegg')
        nets <- lapply(1:npath, function(a) if(names(pathways)[a] %in% names(paths)){
          graphite::pathwayGraph(paths[[which(names(paths)==names(pathways)[a])]])})
        
        gg <- lapply(1:npath, function(a){
          if (is.null(nets[[a]])){
            NULL
          } else {
            g <- graph::subGraph(intersect(pathways[[a]],nodes(nets[[a]])), nets[[a]])
            g <- igraph::igraph.from.graphNEL(g)
            igraph::get.edgelist(g)
          }
        }
        )
        el <- do.call(rbind, gg)
        g <- igraph::graph_from_edgelist(el, directed = FALSE)
        if (length(V(g)$name) < nrow(x)){
          g <- igraph::add_vertices(g, length(setdiff(rownames(x), V(g)$name)), name=setdiff(rownames(x), V(g)$name))
        }
        Adj <- as.matrix(igraph::get.adjacency(as.undirected(g, mode="collapse")))
        
        if (max(Adj)>1){
          Adj <- 1 *(Adj!=0)
        }
      }
    }
    
    if (!is.na(file_e)){
      ## user uploaded
      dat <- read.csv(file_e) 
      is.directed <- (as.character(dat$direction[1])=='directed')
      Adj <- matrix(0, p, p)
      
      if (length(intersect(rownames(x),c(dat$src,dat$dest)))==0){
        Adj <- matrix(0, p, p)
      } else {
        g <- igraph::graph_from_data_frame(dat, directed=is.directed)
        
        if (length(V(g)$name)>nrow(x)){
          ## subgraph
          g <- igraph::induced_subgraph(g, which(V(g)$name %in% rownames(x) == 1))
        }
        
        if (length(V(g)$name) < nrow(x)){
          g <- igraph::add_vertices(g, length(setdiff(rownames(x), V(g)$name)), name=setdiff(rownames(x), V(g)$name))
        }
        
        if (is_dag(g)){
          reOrder <- igraph::topo_sort(g,"in");
          Adj <- 1*Adj[reOrder, reOrder]
        } else {
          # treated as undirected
          Adj <- as.matrix(igraph::get.adjacency(g, type="both"))
          Adj <- Adj + t(Adj)
          if (max(Adj)>1){
            Adj <- 1*(Adj>0)
          }
        }
      }
    }
    
    if (is.null(file_ne)){
      Zero_Adj <- NULL
      if (is.null(Adj)){
        warning('No network information provided!')
      }
    } else if (is.na(file_ne)){
      if (is.null(Adj)){
        warning('No network information provided!')
        Zero_Adj <- NULL
      } else {
        Zero_Adj <- matrix(1, p, p) - Adj - diag(1, p)
      }
    } else {
      dat_ne <- read.csv(file_ne) 
      is.directed <- (as.character(dat_ne$direction[1])=='directed')
      g_ne <- igraph::graph_from_data_frame(dat_ne, directed=is.directed)
      if (length(V(g_ne)$name)>nrow(x)){
        ## subgraph
        g_ne <- igraph::induced_subgraph(g_ne, which(V(g_ne)$name %in% rownames(x) == 1))
      }
      
      if (length(V(g_ne)$name) < nrow(x)){
        g_ne <- igraph::add_vertices(g_ne, length(setdiff(rownames(x), V(g_ne)$name)), name=setdiff(rownames(x), V(g_ne)$name))
      } 
      
      Zero_Adj <- as.matrix(igraph::get.adjacency(g_ne, type="both"))
      if (max(Zero_Adj)>1){
        Zero_Adj <- 1*(Zero_Adj>0)
      }          
    }
    
    if (is.null(Adj)){
      x_c <- x
    } else {
      x_c <- x[match(rownames(Adj),rownames(x)),]
      B <- B[,match(rownames(Adj),rownames(x))]
      if (!is.null(Zero_Adj) && max(Zero_Adj)==1){
        Zero_Adj <- Zero_Adj[match(rownames(Adj),rownames(Zero_Adj)),match(rownames(Adj),rownames(Zero_Adj))]
      }
    }
    
    ## estimate_network
    current_data <- vector("list", ncond)
    Amat <- vector("list", ncond)
    
    if (estimate_network){
      for (k in 1:ncond){
        current_data[[k]] <- x_c[,(group==k)]
        ## Estimate the partial correlation matrices
        ## Full estimation might be slow if p is large.
        ## Use cluster graphical lasso. If input network information is inconsistent with cluster structure, discard them.
        if (is_dag(g)){
          fit <- netEst.dir(current_data[[k]], one=Adj, zero=Zero_Adj,lambda=lambda_c*sqrt(log(p)/n[k]))
        } else {
          fit <- netEst.undir(current_data[[k]], one=Adj, zero=Zero_Adj,lambda=lambda_c*sqrt(log(p)/n[k]), rho=0.1*sqrt(log(p)/n[k]), eta=eta)
        }
        Amat[[k]] <- fit$Adj
      } 
    }
    
    return(list(Amat=Amat, Adj=Adj, Zero_Adj=Zero_Adj, B=B))
  }
