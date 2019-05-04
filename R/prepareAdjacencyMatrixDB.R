prepareAdjacencyMatrixDB <-
  function(X, 
           genes,
           group,         
           databases = NULL,
           file_e=NULL,  
           file_ne=NULL,
           lambda_c=1,
           eta=0.5
  ) {

    user_info       <- checkUserEdges(file_e, file_ne)
    if(any(! c(user_info[["user_edges"]][["src"]], user_info[["user_edges"]][["dest"]], 
               user_info[["user_non_edges"]][["src"]], user_info[["user_non_edges"]][["dest"]]) %in% genes)){
      stop("Some genes in user supplied edges or non-edges were not specified in gene list")
    }

    # What if we don't find any edges? Should probably throw warning
    if(!is.null(databases)){
      edgeL_info      <- obtainEdgeList(genes, databases)
      edgeL           <- edgeL_info[["edgelist"]]
      edgeL_freq      <- if(nrow(edgeL) == 0) {warning("No edges found in databases") 
                                               NULL}
                         else                  edgeL[, .(frequency = .N), by = .(base_gene_src, base_id_src, base_gene_dest, base_id_dest)]
    } else{
      edgeL_info  <- NULL
      edgeL_freq  <- NULL
    }

    
    edgeL_freq_usr  <- addUserEdges(edgeL_freq, user_info[["user_edges"]])
    
    net_info        <- convertEdgeListToZeroOne(edgeL_freq_usr, user_info[["user_non_edges"]], genes, edgeL_info[["genes_not_in_dbs"]])
    
    if(net_info[["directed"]]) {reorder <- rownames(net_info[["ones"]]); X <- X[reorder, ]}
 
    n <- table(group)
    p <- nrow(X)

    network <- lapply(unique(group), estimateNetwork, X, group, net_info, n, p, lambda_c, eta)
    network <- setNames(network, unique(group))
    network_reshaped  <- Reduce(function(x,y) Map(function(x2, y2) {if(identical(class(x2), "list")) c(x2,       list(y2))
                                                                    else                             c(list(x2), list(y2))}, x, y), network)
    network_reshaped2 <- lapply(network_reshaped, setNames, unique(group))

    return(network_reshaped2)
    
  }

# Check user edges and non-edges. 
  checkUserEdges <- function(edgelist, non_edges){
    
    # Checking edgelist
    if(!is.null(edgelist)){
      
      # Read in as data.table
      edgelist  <- fread(edgelist)
  
      # If not 5 or 4 columns specified, throw error
      if(ncol(edgelist) == 5) {
        setnames(edgelist, names(edgelist)[1:5], c("base_gene_src", "base_id_src", "base_gene_dest", "base_id_dest", "frequency"))
      }else if(ncol(edgelist) == 4){
        setnames(edgelist, names(edgelist)[1:4], c("base_gene_src", "base_id_src", "base_gene_dest", "base_id_dest"))
      }else{
        stop(paste0(ncol(edgelist), " column(s) specified in user edgelist. Edgelist should have either 5 or 4 columns"))
      }
        edgelist[, c("base_gene_src", "base_id_src", "base_gene_dest", "base_id_dest") := lapply(.SD, as.character), .SDcols = c("base_gene_src", "base_id_src", "base_gene_dest", "base_id_dest")]
      
      # Make sure frequency is either not existant or numeric value > 0 (I suppose can be 0.5, 1, 1.2, etc)
      if("frequency" %in% names(edgelist)){
        if(any(is.na(edgelist[["frequency"]])))       stop("NA's detected in \"frequency\" column of user edgelist. If specified, \"frequency\" should be a numeric value > 0")
        if(any(!is.numeric(edgelist[["frequency"]]))) stop("Non-numeric values detected in \"frequency\" column of user edgelist. If specified, \"frequency\" should be a numeric value > 0")
        if(any(edgelist[["frequency"]] <=0))          stop("Numeric values <= 0 detected in \"frequency\" column of user edgelist. Please specify numeric values > 0")
      } else{
        edgelist[, frequency := 1]
      }
      
    }
    
    # Checking non_edges
    if(!is.null(non_edges)){
      
      # Read in as data.table
      non_edges <- as.data.table(fread(non_edges))
      
      # If not 4 columns specified, throw error
      if(ncol(non_edges) == 4){
        # Rename first 4 columns
        setnames(non_edges, names(non_edges)[1:4], c("base_gene_src", "base_id_src", "base_gene_dest", "base_id_dest"))
      } else{
        stop(paste0(ncol(non_edges), " column(s) specified in user non-edgelist. Non-edgelist should have exactly 4 columns"))
      }
        non_edges[, c("base_gene_src", "base_id_src", "base_gene_dest", "base_id_dest") := lapply(.SD, as.character), .SDcols = c("base_gene_src", "base_id_src", "base_gene_dest", "base_id_dest")]
    }
    
    # Error if user specifies same edge and non-edge
    if(!is.null(edgelist) & !is.null(non_edges)){
      error_edges <- edgelist[non_edges, on = .(base_gene_src, base_id_src, base_gene_dest, base_id_dest), nomatch = 0L]
      
      if(nrow(error_edges) != 0) {error_src  <- paste0(error_edges[["base_id_src"]], ":", error_edges[["base_gene_src"]])
                                  error_dest <- paste0(error_edges[["base_id_dest"]], ":", error_edges[["base_gene_dest"]])
                                  stop(paste0("Edges identified in both user edgelist and user non-edgelist. Please specify each edge in either the edgelist or the non-edgelist: \n", 
                                               paste0(error_src, " -> ", error_dest, collapse = ", ")))}
    }
    
    # Return list of edited edgelist and non_edges
    return(list(user_edges     = edgelist,
                user_non_edges = non_edges))
  }



# Add user edges. Error if they specify same edge more than once. Also, warning and use user edge if an edge is both user specified and in our database
  addUserEdges <- function(non_user_edges, user_edges){
    
    if(is.null(user_edges) & !is.null(non_user_edges))      return(non_user_edges)
    else if(is.null(user_edges) & is.null(non_user_edges))  return(NULL)
    else if(!is.null(user_edges) & is.null(non_user_edges)) return(user_edges)
    else{
    # Otherwise we have user edges & non_user_edges so do all the checks
    
    # Throw error if user specifies same edge more than once
    specified_n_edges <- user_edges[, .(.N), by = .(base_gene_src, base_id_src, base_gene_dest, base_id_dest)]
    if(any(specified_n_edges[["N"]] > 1)) { error_idxs        <- which(specified_n_edges[["N"]] > 1)
                                            error_src_edges   <- paste0(specified_n_edges[["base_id_src"]][error_idxs], ":", specified_n_edges[["base_gene_src"]][error_idxs]) 
                                            error_dest_edges  <- paste0(specified_n_edges[["base_id_dest"]][error_idxs], ":", specified_n_edges[["base_gene_dest"]][error_idxs]) 
                                            stop(paste0("The following edges were specified multiple times in the user edgelist. Please only enter each edge once :\n",
                                                         paste0(error_src_edges, " -> ", error_dest_edges, collapse = "\n")))}
    
    
    all_edges       <- rbindlist(list("no" = non_user_edges, "yes" = user_edges), use.names = TRUE, fill = TRUE, idcol = "user_specified")
    
    # Warning if user specified edge & in non_user_edges. Use user edge instead
    idxs            <- all_edges[, {if(.N > 1) { warning(paste0("Edge ", base_id_src, ":", base_gene_src, " -> ", base_id_dest, ":",base_gene_dest, " identified in database and in user input. Using user defined frequency instead."))
                                                 .I[.SD == "yes"]
                                                } else{.I}}, by = .(base_gene_src, base_id_src, base_gene_dest, base_id_dest), .SDcols = "user_specified"][["V1"]]
    
    
    return(all_edges[idxs][, .(base_gene_src, base_id_src, base_gene_dest, base_id_dest, frequency)])
    }
    
    
  }



# Convert edgelist to zero / one matrices. Edgelist should include user edges. Also output frequency of edges for later development work. Determines for you whether or not it is directed.
# If genes are not in the databases, we don't make them non-edges. We are essentially saying "we have no info" about them.
  convertEdgeListToZeroOne <- function(edgelist, non_edges, genes, genes_not_in_dbs){
    #This is the only place we need to paste names together. I like having them as separate variables
    genes            <- paste0(names(genes), ":", genes)
    genes_not_in_dbs <- if(!is.null(genes_not_in_dbs)) paste0(names(genes_not_in_dbs), ":", genes_not_in_dbs)
    edgelist         <- if(!is.null(edgelist))         edgelist[, c("base_gene_src", "base_gene_dest")  := .(paste0(base_id_src, ":",  base_gene_src),
                                                                       paste0(base_id_dest, ":", base_gene_dest))][, c("base_id_src", "base_id_dest") := NULL]
    non_edges        <- if(!is.null(non_edges))        non_edges[, c("base_gene_src", "base_gene_dest") := .(paste0(base_id_src, ":",  base_gene_src),
                                                                       paste0(base_id_dest, ":", base_gene_dest))][, c("base_id_src", "base_id_dest") := NULL]
    
    # Start by removing any user specified non_edges if we see them
    if(!is.null(non_edges)){
      #Warn if conflict between non_edges and edgelist. Say we'll use user non-edges
      warn_edges <- edgelist[non_edges, on = .(base_gene_src, base_gene_dest), nomatch = 0L]
      
      if(nrow(warn_edges) != 0) {
        warning("Edges identified in both database edgelist and user non-edgelist. User non-edgelist will be used instead : \n", 
                                        paste0(paste0(warn_edges[["base_gene_src"]], " -> ", warn_edges[["base_gene_dest"]]), collapse = ", "))
        #If we have overlapping edges, remove from edgelist
        edgelist <- edgelist[!non_edges, on = .(base_gene_src, base_gene_dest)]
        
      }
    }
    
    # If we have edges (including user specified), make the zero / one matrices based on those edges
    if(!is.null(edgelist)){
      ones_freq                                          <- matrix(0, nrow = length(genes), ncol = length(genes), dimnames = list(genes, genes))
      ones_freq[cbind(edgelist[["base_gene_src"]], 
                      edgelist[["base_gene_dest"]])]     <- edgelist[["frequency"]]
      diag(ones_freq)                                    <- 0
    
      # Check if edges are directed
      directed                                           <- igraph::is_dag(igraph::graph_from_edgelist(as.matrix(edgelist[, c("base_gene_src", "base_gene_dest"), with = FALSE])))
      
      # If directed, then no undirected edges so can simply return the ones_weighted, ones, and zeros. If undirected need to edit a little
      # Make any directed edges undirected and give same weight
      if(!directed){
        #Faster than x + t(x)
        additional_dir_edges <- edgelist[, is_undirected := edgelist[edgelist, .(is_undirected = max(!is.na(x.base_gene_src))), on = .(base_gene_src = base_gene_dest, base_gene_dest = base_gene_src), by = .EACHI][["is_undirected"]]][is_undirected == 0]
        
        #Removing non_edges so we don't accidently overwrite these. We remove both directions from the user specified non_edges
        if(!is.null(non_edges)) additional_dir_edges <- additional_dir_edges[!non_edges, on = .("base_gene_src" = base_gene_dest, "base_gene_dest" = base_gene_src)]
        
        ones_freq[cbind(additional_dir_edges[["base_gene_dest"]], additional_dir_edges[["base_gene_src"]])] <- additional_dir_edges[["frequency"]]
        
        #Also assume user non_edge undirected if graph is undirected. This is b/c ones MUST be symmetric
        if(!is.null(non_edges)) ones_freq[cbind(non_edges[["base_gene_dest"]], non_edges[["base_gene_src"]])] <- 0
        
        if(!all(ones_freq == t(ones_freq))) stop("Undirected graph, but ones matrix is not symmetric")
    
      } else{
        #If directed, calculate correct order, error check, and then reOrder ones_freq and zeros
        preReOrder    <- igraph::topo_sort(igraph::graph_from_edgelist(as.matrix(edgelist[, c("base_gene_src", "base_gene_dest"), with = FALSE])), mode = "in")
        reOrder       <- c(names(preReOrder), setdiff(rownames(ones_freq), names(preReOrder)))
        to_test       <- ones_freq[reOrder, reOrder]
        stopifnot(all(rownames(to_test) %in% rownames(ones_freq)) & all(colnames(to_test) %in% colnames(ones_freq)))
        stopifnot(all(to_test[upper.tri(to_test)] == 0))
        
        ones_freq     <- to_test
        
        #Test if everything is lower.tri. If it is then it is Wald ordering. If not, error
        if(any(ones_freq[upper.tri(ones_freq)] != 0)) stop("Directed graph detected, but incorrect Wald ordering. Please arrange gene list and data matrix to be lower triangular")
      }
      
      ones <- 1*(ones_freq> 0)
      
      #Lastly, make 0s matrix. We should have done the ones correctly so zeros is easy
      zeros                                              <- 1*(ones_freq==0)
      diag(zeros)                                        <- 0
      zeros[, genes_not_in_dbs]                          <- 0
      zeros[genes_not_in_dbs, ]                          <- 0
      #If we overwrote some of the user non-edges put them back in
      zeros[cbind(non_edges[["base_gene_src"]],
                  non_edges[["base_gene_dest"]])]        <- 1
      
      if(any( (ones + zeros) > 1 )) stop("Overlapping information identified")
      
      
    } 
    # Otherwise, we dont have edges so set ones to be all 0 (assumes we searched through DB). Directed FALSE. Zeros is all ones
    else{
      ones_freq   <- ones <- matrix(0, nrow = length(genes), ncol = length(genes), dimnames = list(genes, genes))
      zeros       <- matrix(1, nrow = length(genes), ncol = length(genes), dimnames = list(genes, genes))
      diag(zeros) <- 0
      directed <- FALSE
    }
    
    return(list(ones_freq = ones_freq, ones = ones, zeros = zeros, directed = directed))
  
  }

#Simple wrapper to estimateNetwork for given grp value in group
estimateNetwork <- function(grp, X, group, net_info, n, p, lambda_c, eta) {
                    if(is.null(rownames(X))) dimnames(net_info$zeros) <- dimnames(net_info$ones) <- dimnames(net_info$ones_freq) <- NULL
                       
                    if(net_info$directed){
                      return(netEst.dir(X[, group == grp], zero = net_info$zeros, one = net_info$ones, lambda=lambda_c*sqrt(log(p)/n[grp])))
                    }else{
                      return(netEst.undir(X[, group == grp], zero = net_info$zeros, one = net_info$ones, lambda=lambda_c*sqrt(log(p)/n[grp]), rho=0.1*sqrt(log(p)/n[grp]), eta=eta))
                    }
}

