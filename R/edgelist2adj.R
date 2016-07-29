edgelist2adj <-
function(file, vertex.names, mode=c("directed", "undirected")) {
    this.call <- match.call()
    mode <- match.arg(mode)
    is.directed <- ifelse(mode=="directed", TRUE, FALSE)
    
    dat = read.table(file, header=TRUE) 
    el = as.matrix(dat) 
    # el[,1] = as.character(el[,1])
    # el[,2] = as.character(el[,2])
    # g = graph_from_edgelist(el, directed=is.directed)
    g = graph.data.frame(el, directed=is.directed)
    
    extra.vertices <- setdiff(vertex.names, V(g)$name)
    if (length(extra.vertices)>0){
      g <- add.vertices(g, nv = length(extra.vertices), name = extra.vertices)
    }
    Adj = as.matrix(get.adjacency(g, type="both"))
    if (is.directed){
     Adj = t(Adj) 
     if (max(Adj[upper.tri(Adj)])>0){
       stop('The adjacency matrix is not lower triangular! Check ordering of the variables!')
     }
    }
    
    reOrder <- match(vertex.names,rownames(Adj))
    Adj <- Adj[reOrder, reOrder]
    
    return(Adj)
  }
