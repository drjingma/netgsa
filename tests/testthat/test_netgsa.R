
##---example-----
test_that("data dimension",{
  # paths <- pathways('hsapiens','kegg')
  cat('Loading data...\n')
  system.file("doc/breastcancer2012.rda",package="netgsa")
  expect_equal(nrow(x), 2598)
})

###------------------BIC--------------------------
# cat('Run the small example ... \n')
# ## Consider genes from the "ErbB signaling pathway" and "Jak-STAT signaling pathway"

# ncond <- length(unique(group))
# Amat <- vector("list",ncond)
# for (k in 1:ncond){
#   data_c <- sx[,(group==k)]
#   fitBIC <- bic.netEst.undir(data_c,one=out$Adj,zero=out$Zero_Adj,
#                              lambda=seq(1,10)*sqrt(log(p)/ncol(data_c)),eta=0.1)
#   fit <- netEst.undir(data_c,one=out$Adj,zero=out$Zero_Adj,
#                       lambda=which.min(fitBIC$BIC)*sqrt(log(p)/ncol(data_c)),eta=0.1)
#   Amat[[k]] <- fit$Adj
# }
# test <- NetGSA(Amat, sx, group, pathways = out$B, lklMethod = 'REHE')


###------------------DAG--------------------------
## NetGSA can also work with directed acyclic graphs (DAGs).
# e.g. the "Adrenergic signaling in cardiomyocytes" pathway from KEGG is a DAG.

# load KEGG network
# paths <- pathways('hsapiens','kegg')
# 
# # get pathway topology using pathwayGraph() from the graphite package
# pg <- pathwayGraph(paths[[which(names(paths)=="Adrenergic signaling in cardiomyocytes")]])
# 
# # extract genes in both the pathway and the data
# genenames <- intersect(nodes(pg), rownames(x))
# p <- length(genenames)
# gx <- x[match(genenames,rownames(x)),]
# 
# # extract the subgraph
# g <- subGraph(genenames, pg)
# g <- igraph.from.graphNEL(g)
# print(is_dag(g))
# 
# # reorder the variables and get the adjacency matrix
# reOrder <- topo_sort(g,"in")
# Adj <- as.matrix(get.adjacency(g))
# Adj <- Adj[reOrder,reOrder]
# B <- matrix(rep(1,p),nrow=1)
# rownames(B) <- "Adrenergic signaling in cardiomyocytes"
# 
# Amat <- vector("list", ncond)
# for (k in 1:2){
#   data_c <- gx[,which(group==k)]
#   Amat[[k]] <- netEst.dir(data_c, one = Adj)$Adj
# }
# test <- NetGSA(Amat, gx, group, pathways = B, lklMethod = 'REHE')
