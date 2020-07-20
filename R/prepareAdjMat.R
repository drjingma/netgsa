prepareAdjMat <-
  function(x, 
           group,         
           databases=NULL,
           cluster=TRUE,
           file_e=NULL,  
           file_ne=NULL,
           lambda_c=1,
           penalize_diag = TRUE,
           eta=0.5
  ) {
    lambda_c = sort(lambda_c, decreasing = TRUE)
    genes <- setNames(gsub(".*:(.*)", "\\1", rownames(x)), gsub("(.*):.*", "\\1", rownames(x)))
                      
    user_info       <- checkUserEdges(file_e, file_ne)

    if(any(! c(user_info[["user_edges"]][["src"]], user_info[["user_edges"]][["dest"]], 
               user_info[["user_non_edges"]][["src"]], user_info[["user_non_edges"]][["dest"]]) %in% genes)){
      stop("Some genes in user supplied edges or non-edges were not specified in gene list")
    }
    
    if (!all(class(databases) %in% c("obtainedEdgeList", "character", "NULL"))){
      stop(paste0("Don't know how to handle databases object of type: ", class(databases)))
    }

    # What if we don't find any edges? Should probably throw warning
    if(!is.null(databases)){
      
      edgeL_info <- if(identical(class(databases), "character")){
                        obtainEdgeList(rownames(x), databases)
                    } else if(identical(class(databases), "obtainedEdgeList")){
                        copy(databases) #To avoid overwriting if data.table
                    }
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
    
    if(net_info[["directed"]]) {reorder <- rownames(x);
                                net_info[c("ones","ones_freq", "zeros")] <- lapply(net_info[c("ones","ones_freq", "zeros")], function(x) x[reorder, reorder])
    }
 
    n <- table(group)
    p <- nrow(x)

    #If user doesn't specify clustering choose automatically based on p
    if(is.null(cluster)){
      if(p > 2500) cluster <- TRUE
      else         cluster <- FALSE
    }
  
    #Can only cluster if we have information
    #Always cluster by C.C.
    if(!all(net_info$ones == 0)){
      net_clusters <- obtainClusters(net_info$ones, rownames(x), cluster)
      stopifnot(all(names(net_clusters) == rownames(x))) #clusters returned must be in same order as data matrix x
    } else{
      net_clusters <- NULL
    }
    

    # Each row has same lambda, its just rescaled for each cluster based on # of genes in that cluster and # of obs for group
    network <- lapply(unique(group), function(g) {estimateNetwork(grp = g, X = x, group = group, net_info = net_info, 
                                                                  n = n, p = p, lambda_c = lambda_c, eta = eta, 
                                                                  net_clusters = net_clusters, penalize_diag = penalize_diag)})
    network <- rbindlist(network)

    if(!net_info$directed){
      empcov <- cov(t(x)) #Could either do this quick calc or store one for each group/lambda combo   
      network[, c("bic", "df") := {calcBIC(invcov[[1]], empcov, n[group])}, by = .(group, lambda)]
      network <- network[,.SD[which.min(bic)], by = group]
    }
    
    # For directed, should only return 1 row per condition. We don't test multiple lambdas
    stopifnot(nrow(network) == length(unique(group)))
    ugrp <- unique(group)
    Adj_mats    <- setNames(lapply(ugrp, function(grp) network[group == grp][["Adj"]][[1]]), ugrp)
    invcov_mats <- setNames(lapply(ugrp, function(grp) network[group == grp][["invcov"]][[1]]), ugrp)
    lambda_l    <- setNames(lapply(ugrp, function(grp) network[group == grp][["lambda_l"]][[1]]), ugrp)

    result <- list(Adj = c(Adj_mats, list("edgelist" = edgeL_freq_usr)), invcov = invcov_mats, lambda = lambda_l)
    return(result)
    
}


# Helper functions -----------------------------------------------------------------

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
  idxs            <- all_edges[, {if(.N > 1) { #warning(paste0("Edge ", base_id_src, ":", base_gene_src, " -> ", base_id_dest, ":",base_gene_dest, " identified in database and in user input. Using user defined frequency instead.\n"))
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
      #If we have overlapping edges, remove from edgelist. Edgelist should not include any user edges that are also in the user non-edges (would have errored out)
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
      #Faster than x + t(x). Getting directed edges
      additional_dir_edges <- edgelist[, is_undirected := edgelist[edgelist, .(is_undirected = max(!is.na(x.base_gene_src))), on = .(base_gene_src = base_gene_dest, base_gene_dest = base_gene_src), by = .EACHI][["is_undirected"]]][is_undirected == 0]
      
      #Removing non_edges so we don't accidently overwrite these. We remove both directions from the user specified non_edges
      if(!is.null(non_edges)) additional_dir_edges <- additional_dir_edges[!non_edges, on = .("base_gene_src" = base_gene_dest, "base_gene_dest" = base_gene_src)]
      
      ones_freq[cbind(additional_dir_edges[["base_gene_dest"]], additional_dir_edges[["base_gene_src"]])] <- additional_dir_edges[["frequency"]]

      #Also assume user non_edge undirected if graph is undirected. This is b/c ones MUST be symmetric
      if(!is.null(non_edges)) ones_freq[cbind(non_edges[["base_gene_dest"]], non_edges[["base_gene_src"]])] <- 0
      
      if(!all(1*(ones_freq>0) == t(1*(ones_freq>0)))) stop("Undirected graph, but ones matrix is not symmetric")
      ## Some of the edges reported (e.g. in breast cancer example) point to themselves. That is we see for example
      ## and edge between ENTREZID:781 <-> ENTREZID:781, same for ENTREZID:4338 and ENTREZID:4854.
      ## This is why the nrow(edgelist) + nrow(additional_dir_edges) do not equal sum(ones_freq>0)
  
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
    zeros       <- matrix(0, nrow = length(genes), ncol = length(genes), dimnames = list(genes, genes))
    #diag(zeros) <- 0
    directed <- FALSE
  }
  
  return(list(ones_freq = ones_freq, ones = ones, zeros = zeros, directed = directed))

}

combineSingletons <- function(c_vec){
  members_size_1 <- as.numeric(names(Filter(function(x) x == 1, table(c_vec))))
  if(length(members_size_1) != 0) c_vec[c_vec %in% members_size_1] <- c_vec[1]
  c_vec_refactor <- setNames(as.numeric(factor(c_vec)), names(c_vec))
  return(c_vec_refactor)
}

#Function to get the clusters for a 0-1 Adjacency matrix
obtainClusters <- function(A, order, cluster){

  #Converting to undirected for purposes of clustering
  A_graph <- igraph::as.undirected(igraph::graph_from_adjacency_matrix(A))
  ccs <- igraph::components(A_graph)$membership
  if (!cluster){ return(combineSingletons(ccs))} #Always cluster by cc's for speed up. If actually want to cluster, loop over ccs
  
  ccs_tbl <- table(ccs)
  ccs_to_cluster <- as.numeric(names(ccs_tbl[which(ccs_tbl > 1000)]))

  cluster_methods_list  <- list(igraph::cluster_walktrap, igraph::cluster_leading_eigen, igraph::cluster_fast_greedy, igraph::cluster_label_prop, igraph::cluster_infomap, igraph::cluster_louvain)
  
  #Within big connected components, perform clustering
  clustered_ccs <- lapply(ccs_to_cluster, function(c){
    subg <- igraph::induced_subgraph(A_graph, names(ccs[ccs == c]))
    
    cluster_method_info   <- rbindlist(lapply(cluster_methods_list, analyzeClusterMethod, graph = subg), use.names = TRUE, fill = TRUE)

    #Choose method with lowest edge-loss among methods with max cluster size < 1000. If tied choose one with lowest max cluster. If still tied, choose random one
    #   If all > 1000 choose with smallest max cluster size
    if(all(cluster_method_info$max_cluster_size > 1000)){
      best_membership <- cluster_method_info[sample(which(max_cluster_size == min(max_cluster_size)), 1), membership][[1]] #Random if ties
    } else{
      best_membership <- cluster_method_info[max_cluster_size < 1000][which(perc_lost_edges == min(perc_lost_edges))][sample(which(max_cluster_size == min(max_cluster_size)),1), membership][[1]] #Random if ties
    }

    #Make names unique
    return(setNames(paste0(c, "_", best_membership), names(best_membership)))
  })
 
  #Leftover ccs get put into their own clusters
  leftover_ccs <- ccs[!ccs %in% ccs_to_cluster]

  final_clusters <- combineSingletons(c(leftover_ccs, do.call(c, clustered_ccs)))
  #Must be in same order as data rows
  return(final_clusters[order])
}

#Simple function to return important info from each cluster method
analyzeClusterMethod <- function(cluster_method, graph){
  x <- cluster_method(graph)
  
  data.table(algorithm                 = igraph::algorithm(x),
             max_cluster_size          = max(tabulate(igraph::membership(x))),
             modularity                = igraph::modularity(x),
             n_between_clust_edges     = sum(igraph::crossing(x, graph)),
             n_total_edges             = igraph::ecount(graph),
             perc_lost_edges           = sum(igraph::crossing(x, graph)) / igraph::ecount(graph),
             membership                = list(igraph::membership(x)))
}


#Loops through clusters, subsets data and zero/one matrices and estimates the network.
#For clusters of size 1 it returns a matrix(0,1,1)
netEstClusts <- function(grp, X, group, net_info, n, lambda_c, eta, net_clusters, penalize_diag){
  net_info_clusts <- lapply(unique(net_clusters), function(clust){
                                #Subsetting to clusters
                                idxs       <- net_clusters == clust
                                #Make sure this is same returned type as netEst.undir for multiple lambdas.
                                if(sum(idxs) == 1) {return(list(Adj = rep(list(matrix(0, sum(idxs), sum(idxs), dimnames = list(names(net_clusters[idxs]), names(net_clusters[idxs])))), length(lambda_c)), 
                                                                invcov = rep(list(matrix(0,sum(idxs), sum(idxs), dimnames = list(names(net_clusters[idxs]), names(net_clusters[idxs])))), length(lambda_c)), 
                                                                lambda = rep(list(NA), length(lambda_c))))}

                                X_new      <- X[idxs, group == grp]
                                zero_new   <- net_info$zeros[idxs, idxs]
                                one_new    <- net_info$ones[idxs, idxs]
                                
                                ### Return directed
                                if(net_info$directed){
                                  return(netEst.dir(X_new, zero = zero_new, one = one_new))
                                  
                                ### Return undirected
                                }else{
                                  n = n[grp]; p = sum(idxs)
                                  lambda_try = lambda_c*sqrt(log(p)/n)
                                  #tuning_params <- getUndirTuningParams(X_new, zero = zero_new, one = one_new, lambda_c = lambda_c, n = n[grp], p = sum(idxs), eta = eta)
                                  return(netEst.undir(X_new, zero = zero_new, one = one_new, lambda=lambda_try, rho=(0.1 * sqrt(log(p)/n)), penalize_diag = penalize_diag, eta=eta))
                                }
                              })
  
  #For now return as list
  # list of length lambda. So will give length(lambda) rows in our data.table
  if(!net_info$directed){
    A_list <- lapply(1:length(lambda_c), function(i) lapply(net_info_clusts, function(x) x[["Adj"]][[i]])) 
    invcov <- lapply(1:length(lambda_c), function(i) lapply(net_info_clusts, function(x) x[["invcov"]][[i]]))
    lambda <- lapply(1:length(lambda_c), function(i) lapply(net_info_clusts, function(x) x[["lambda"]][[i]]))
  } else{
    A_list <- list(lapply(net_info_clusts, function(x) x[["Adj"]]))
    invcov <- list(lapply(net_info_clusts, function(x) x[["infmat"]])) #If directed, we call invcov even though its infmat for consistency w/ undirected
    lambda <- list(lapply(net_info_clusts, function(x) x[["lambda"]]))
    lambda_c <- NA #Directed uses built in lambdas, we never use lambda_c
  }

                   
  # Data.table can store columns of lists
  res = data.table(group = grp, lambda = lambda_c, Adj = A_list, invcov = invcov, lambda_l = lambda)
  stopifnot(nrow(res) == length(lambda_c))
  
  return(res)
}
  
#Simple wrapper to estimateNetwork for given grp value in group
estimateNetwork <- function(grp, X, group, net_info, n, p, lambda_c, eta, net_clusters = NULL, penalize_diag) {
  if(is.null(rownames(X))) dimnames(net_info$zeros) <- dimnames(net_info$ones) <- dimnames(net_info$ones_freq) <- NULL
  
  #No clustering
  if(is.null(net_clusters)){
    if(net_info$directed){
      net_estd = netEst.dir(X[, group == grp], zero = net_info$zeros, one = net_info$ones)
      res = data.table(group = grp, lambda = NA, Adj = list(list(net_estd$Adj)), invcov = list(list(net_estd$infmat)), lambda_l = list(list(net_estd$lambda)))
      return(res)
    }else{
      n2 = n[grp]
      lambda_try = lambda_c*sqrt(log(p)/n2)
      net_estd   = netEst.undir(X[, group == grp], zero = net_info$zeros, one = net_info$ones, lambda=lambda_try, rho=(0.1 * sqrt(log(p)/n2)), penalize_diag = penalize_diag, eta=eta)
      A_list     = lapply(1:length(lambda_c), function(i) list(net_estd[["Adj"]][[i]]))
      invcov     = lapply(1:length(lambda_c), function(i) list(net_estd[["invcov"]][[i]]))
      lambda     = lapply(1:length(lambda_c), function(i) list(net_estd[["lambda"]][[i]]))
      res = data.table(group = grp, lambda = lambda_c, Adj = A_list, invcov = invcov, lambda_l = lambda)
      return(res)
    }
    #Clustering
  } else{ 
    obj <- netEstClusts(grp, X, group, net_info, n, lambda_c, eta, net_clusters, penalize_diag = penalize_diag)
    return(obj)
  }
  
}

# Calculate BIC
calcBIC <- function(invcov_list, empcov, n, eps = 1e-08){
  siginv <- as.matrix(bdiag(invcov_list))
  stopifnot(class(invcov_list) == "list")
  rownames(siginv) <- colnames(siginv) <- do.call(c, lapply(invcov_list, rownames))
  siginv <- siginv[rownames(empcov), colnames(empcov)]
  
  no.edge = sum(abs(siginv) > eps) - ncol(siginv)
  #Faster matrix trace using entrywise product
  bic = sum(t(empcov) * siginv) - determinant(siginv, logarithm = T)$modulus + log(n) * no.edge/(2 * n)
  return(list(bic = bic, df = no.edge))
}
