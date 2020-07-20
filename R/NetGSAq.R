NetGSAq <-
  function(
    x, 	
    group,    
    pathways, 	
    lambda_c = 1,
    file_e = NULL,
    file_ne = NULL,
    lklMethod="REHE",
    cluster = TRUE,
    sampling = TRUE,
    sample_n = NULL,
    sample_p = NULL, 
    minsize=5,
    eta=0.1,           
    lim4kappa=500
  ){
    dbs <- unique(as.character(graphite::pathwayDatabases()$database))
    adj_mats <- prepareAdjMat(x, group, databases = dbs, cluster = cluster, file_e = file_e, file_ne = file_ne, lambda_c = lambda_c, eta = eta)
    net_est  <- NetGSA(adj_mats$Adj, x, group, pathways, lklMethod = lklMethod, sampling = sampling, sample_n = sample_n,
                       sample_p = sample_p, minsize = minsize, eta = eta, lim4kappa = lim4kappa)
    return(net_est)
}
