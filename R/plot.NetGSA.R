plot.NetGSA <- function(x, graph_layout = NULL, rescale_node = c(2,10), rescale_label = c(0.5,0.6)){
  edges_pathways_list <- makePathwayEdges(x$graph$edgelist, x$graph$pathways)
  edges_pathways      <- edges_pathways_list[["edges_pathways"]]
  edges_all           <- edges_pathways_list[["edges_all"]]
  
  cytoscape_open      <- tryCatch({httr::GET("http://localhost:1234/v1/version")$status_code == 200}, error = function(e){return(FALSE)})
  
  fdrCutoffCols <- list("0_005" = list(cutoffs = c(0, 0.05), cols = c("#FF0000", "#F74A4A")), "005_01" = list(cutoffs = c(0.05001, 1), cols = c("#F78F8F", "#F2B8B8")), 
                        "01_02" = list(cutoffs = c(0.1001, 0.2), cols = c("#F7D7D7", "#FCE8E8")), "02_1" = list(cutoffs = c(0.2001, 1), cols = c( "#F7F5F5", "#F7F5F5")))
  
  if(cytoscape_open){
    cyto <- plot_cytoscape.NetGSA(edges_pathways, edges_all, x$graph$pathways, x$results, x$graph$gene.tests, fdrCutoffCols = fdrCutoffCols, graph_layout = graph_layout)
    legend_title <- "Cytoscape plot legend"
  } else{
    warning("For better visualization results, please install and open Cytoscape")
    cyto <- NULL
    legend_title <- "Igraph plot legend"
  }
  plotCytoLegend(x$results, fdrCutoffCols, legend_title, igraph = ! cytoscape_open)
  plot_igraph.NetGSA(edges_pathways, x$results, x$graph$pathways, fdrCutoffCols = fdrCutoffCols, nodes_layout = cyto[["node_locations"]], cytoscape_open = cytoscape_open, graph_layout = graph_layout, rescale_node = rescale_node, rescale_label = rescale_label)
  #Reset plot layout
  par(mfrow=c(1,1))
}


plot_cytoscape.NetGSA <- function(edges_pathways, edges_all, pathway_gene_map, pathway_results, gene_results, fdrCutoffCols, graph_layout = NULL, title = "Pathway Network"){
  title <-  gsub(" ", "\\ ", title, fixed = TRUE)
  #Using Cytoscape
  network_ids         <- createNestedNetwork(edges_pathways = edges_pathways, edges_all = edges_all, pathway_vertices = pathway_gene_map, main = title)
  
  #Also adding edge weights incase we want that
  RCy3::setCurrentNetwork(network_ids$networks[1])
  edge_df             <- setnames(copy(edges_pathways)[, id := paste0(src_pathway, " (pp) ", dest_pathway)], "weight_sum", "weight")
  RCy3::loadTableData(data = pathway_results, data.key.column = "pathway", table = "node", table.key.column = "name", network = network_ids$networks[1])
  RCy3::loadTableData(data = edge_df, data.key.column = "id", table = "edge", table.key.column = "name", network = network_ids$networks[1])
  if(is.null(graph_layout)){
    graph_layout = 'force-directed defaultSpringCoefficient=0.00000004 defaultSpringLength=100'
  }
  RCy3::layoutNetwork(graph_layout) ##Ooohhh this one is good!
  node_locations      <- getCytoscapeXYCoords(network_ids$networks[1])
  
  #Adding mappings
  teststat_lim <- max(abs(min(pathway_results$teststat)), abs(max(pathway_results$teststat)))
  #Mappings are constant for a given style. So these are really editing the "default" visual style mapping. This is done for ALL networks that use "default". Not just the network specified
  RCy3::copyVisualStyle("default", "pathway_style")
  RCy3::copyVisualStyle("default", "gene_style") #For use later
  
  #Setting up "pathway_style"
  RCy3::setNodeColorMapping("teststat", c(-teststat_lim, 0, teststat_lim), c("#FFA500", "#FFFFFF", "#0000FF"), mapping.type = "continuous", style.name = "pathway_style")
  RCy3::setNodeBorderColorMapping("pFdr", unname(do.call(c, lapply(fdrCutoffCols, "[[", "cutoffs"))), unname(do.call(c, lapply(fdrCutoffCols, "[[", "cols"))), mapping.type = "continuous", style.name = "pathway_style")
  RCy3::setNodeWidthMapping("pSize", c(min(pathway_results$pSize), max(pathway_results$pSize)), c(20, 120), style.name = "pathway_style")
  RCy3::setNodeHeightMapping("pSize", c(min(pathway_results$pSize), max(pathway_results$pSize)), c(20, 120), style.name = "pathway_style")
  RCy3::setNodeTooltipMapping("pathway", style.name = "pathway_style")
  
  #These are for specific nodes. Only for network specified
  RCy3::setNodeBorderWidthBypass(pathway_results$pathway, 10, network = network_ids$networks[1])
  RCy3::setNodeShapeBypass(pathway_results$pathway, "ellipse", network = network_ids$networks[1]) #I was getting weird error because some of the nodes werent added!
  RCy3::setVisualStyle("pathway_style")
  
  #Significant Pathways
  sigNodes <- pathway_results$pathway[pathway_results$pFdr <= 0.05]
  if(length(sigNodes) != 0){
    sigNet <- RCy3::createSubnetwork(nodes=sigNodes, nodes.by.col = "pathway", network = network_ids$networks[1], subnetwork.name="Significant Pathways")
  }
  
  return(list(node_locations = node_locations))
}


plot_igraph.NetGSA <- function(edges_pathways, pathway_results, pathway_gene_map, fdrCutoffCols, nodes_layout = NULL, cytoscape_open = FALSE, graph_layout = NULL, rescale_node = c(2,10), rescale_label = c(0.5,0.6)){
  if(is.null(nodes_layout)) nodes_list <- pathway_results[,"pathway", drop = FALSE]
  else                      nodes_list <- nodes_layout
  #Using igraph if Cytoscape not available
  ig                      <- igraph::graph_from_data_frame(edges_pathways, directed = FALSE, vertices = nodes_list)
  pathway_results_ordered <- pathway_results[match(igraph::V(ig)$name, pathway_results$pathway), ] #Make sure this is in correct order
  
  ## Color vertex outline based on FDR
  #Define color ranges of interest
  fdr_ramps               <- list("0_0.05" = colorRamp(fdrCutoffCols[["0_005"]][["cols"]]), "0.05_0.1" = colorRamp(fdrCutoffCols[["005_01"]][["cols"]]), "0.1_0.2" = colorRamp(fdrCutoffCols[["01_02"]][["cols"]]), "0.2_1" = colorRamp(fdrCutoffCols[["02_1"]][["cols"]])) #Coloring p1
  fdr_ramps_vals          <- vapply(pathway_results_ordered$pFdr, function(x){
    if( x < 0.05) return(fdr_ramps[["0_0.05"]](x))
    else if (x < 0.1) return(fdr_ramps[["0.05_0.1"]](x))
    else if (x < 0.2) return(fdr_ramps[["0.1_0.2"]](x))
    else return(fdr_ramps[["0.2_1"]](x))
  }, FUN.VALUE = numeric(3))
  fdr_rgb_vals            <- rgb(t(fdr_ramps_vals)/255) #Coloring p2
  
  ## Color vertex based on teststat
  teststat_lim            <- max(abs(min(pathway_results_ordered$teststat)), abs(max(pathway_results_ordered$teststat)))
  teststat_cols           <- c("orange", "white", "blue") #Orange is small vals (-), blue is large vals (+)
  teststat_ramp           <- colorRamp(teststat_cols)
  std_teststat            <- (pathway_results_ordered$teststat - (-teststat_lim)) / (teststat_lim - (-teststat_lim))
  teststat_rgb_cols       <- rgb(teststat_ramp(std_teststat) /255)
  
  #Vertex size - scale to between 5 and 40
  pSize_0_1_rescale       <- (pathway_results_ordered$pSize - (min(pathway_results_ordered$pSize))) / (max(pathway_results_ordered$pSize) - min(pathway_results_ordered$pSize))
  vertex_sizes            <- pSize_0_1_rescale*(rescale_node[2] - rescale_node[1]) + rescale_node[1]
  
  label_sizes             <- pSize_0_1_rescale*(rescale_label[2] - rescale_label[1]) + rescale_label[1]
  
  if (cytoscape_open)              l <- NULL #Change to layout of cytoscape
  else if (!is.null(graph_layout)) l <- graph_layout(ig)
  else                             l <- igraph::layout_on_sphere(ig)
  
  #Cytoscape open, use coords and plot ig as well as legend on side
  if(cytoscape_open){
    layout_mat <- matrix(c(1,1,2,3), nrow = 2)
    layout(layout_mat, widths = c(8,2))
    par(mai=c(0, 0, 0.5, 0))
    plot(ig, vertex.color = teststat_rgb_cols, vertex.frame.color = fdr_rgb_vals, vertex.size = vertex_sizes, vertex.label.cex = 0.5, layout = l) 
    
    teststat_image <- as.raster(matrix(rev(colorRampPalette(teststat_cols)(20)), ncol=1))
    plot(c(0,0.5),c(0,1),type = 'n', axes = F,xlab = '', ylab = '')
    text(x=0.225, y = c(0,0.5,1), labels = round(c(-teststat_lim, 0, teststat_lim), digits = 1))
    rasterImage(teststat_image, 0, 0, 0.15,1)
    mtext("Node Color\nTest statistic", side = 3, at = c(0.075))
    
    fdr_image <- as.raster(matrix(rev(c(colorRampPalette(fdrCutoffCols[["0_005"]][["cols"]])(20), colorRampPalette(fdrCutoffCols[["005_01"]][["cols"]])(20), colorRampPalette(fdrCutoffCols[["01_02"]][["cols"]])(40), colorRampPalette(fdrCutoffCols[["02_1"]][["cols"]])(320))), ncol=1))
    plot(c(0,0.5),c(0,1),type = 'n', axes = F,xlab = '', ylab = '')
    text(x=0.225, y = c(0,0.05, 0.1, 0.2, 1), labels = c(0,0.05, 0.1, 0.2, 1))
    rasterImage(fdr_image, 0, 0, 0.15,1)
    mtext("Node Border\nFDR adjusted p-value", side = 3, at = c(0.075))
    
  } else{
    #plotCytoLegend(pathway_results, title = "Igraph color legend", igraph = TRUE) #Plot legend in base R because looks better
    igraph::rglplot(ig, vertex.color = teststat_rgb_cols, vertex.size = vertex_sizes, vertex.label.cex = label_sizes, edge.color = "grey98", label.dist = 0, layout = l) #Label size not supported in RGL. No other options to make labels more visible. Label.degree doesnt put them on axis like I thought
  }
}

plotCytoLegend <- function(pathway_results, fdrCutoffCols, title = "Cytoscape plot legend", igraph = FALSE){
  fdr_cols                <- c("red", "white") #Red is small vals (sig) white is large vals (insig)
  teststat_cols           <- c("orange", "white", "blue") #Orange is small vals (-), blue is large vals (+)
  teststat_lim            <- max(abs(min(pathway_results$teststat)), abs(max(pathway_results$teststat)))
  
  
  par(mar=c(2.5,2.5,1,1))
  if (igraph){
    layout_mat <- matrix(c(1,2), nrow = 2)
    layout(layout_mat, heights = c(2,7))
  } else{
    layout_mat <- matrix(c(1,2,3), nrow = 3)
    layout(layout_mat, heights = c(1,3,3))
  }
  
  
  plot.new()
  text(0.5,0.5,title,cex=2,font=2)
  
  teststat_image <- as.raster(matrix(colorRampPalette(teststat_cols)(20), nrow=1))
  plot(c(0,1),c(0,0.5),type = 'n', axes = F,xlab = '', ylab = '')
  text(y=0.175, x = c(0,0.5,1), labels = round(c(-teststat_lim, 0, teststat_lim), digits = 2))
  rasterImage(teststat_image, 0, 0.25, 1,0.5)
  mtext("Node Color\nTest statistic", side = 3, at = c(0.5))
  
  if (!igraph){
    fdr_image <- as.raster(matrix(c(colorRampPalette(fdrCutoffCols[["0_005"]][["cols"]])(20), colorRampPalette(fdrCutoffCols[["005_01"]][["cols"]])(20), colorRampPalette(fdrCutoffCols[["01_02"]][["cols"]])(40), colorRampPalette(fdrCutoffCols[["02_1"]][["cols"]])(320)), nrow=1))
    plot(c(0,1),c(0,0.5),type = 'n', axes = F,xlab = '', ylab = '')
    text(y=0.175, x = c(0,0.05, 0.1, 0.2, 1), labels = c(0,0.05, 0.1, 0.2, 1))
    rasterImage(fdr_image, 0, 0.25, 1,0.5)
    mtext("Node Border\nFDR adjusted p-value", side = 3, at = c(0.5))
  }
}

#Returns cytoscape X & Y coordinates for given network
getCytoscapeXYCoords <- function(network){
  RCy3::setCurrentNetwork(network)
  nodes_all <- RCy3::getAllNodes()
  #For some reason igraph flips the y. So we need to do -y from cytoscape to get same plot in igraph
  return(data.table(nodes = nodes_all, x = RCy3::getNodeProperty(nodes_all, "NODE_X_LOCATION"), y = -RCy3::getNodeProperty(nodes_all, "NODE_Y_LOCATION")))
}




formatPathways <- function(x, pways, graph_layout = NULL){
  gene_results <- x$graph$gene.tests
  cytoscape_open      <- tryCatch({httr::GET("http://localhost:1234/v1/version")$status_code == 200}, error = function(e){return(FALSE)})
  if(!cytoscape_open) stop("formatPathways is only compatible with Cytoscape plots")
  #Layout of each pathway - Slow & prone to crashing
  for (pway in pways){
      if (! pway %in% RCy3::getNetworkList()) {
        stop(paste0("Network does not exist: ", pway))
      }
      #RCy3 has an error in code so doing manually
      cmd <- paste0("network get attribute network=\"", pway, "\" namespace=\"default\" columnList=\"SUID\"")
      net <- RCy3::commandsPOST(cmd)[[1]]
      RCy3::setCurrentNetwork(net)
      curr_nodes <- unlist(RCy3::getAllNodes())
      #Set and edit our new visual style
      RCy3::loadTableData(data = gene_results[J(curr_nodes), ], data.key.column = "gene", table = "node", table.key.column = "name", network = net)
      #Setting up "gene_style"
      RCy3::setNodeColorMapping("pFdr", c(0,1), c("#FF0000", "#FFFFFF"), mapping.type = "continuous", style.name = "gene_style")
      RCy3::setNodeTooltipMapping("gene", style.name = "gene_style")
      RCy3::setNodeShapeDefault("ellipse", style.name = "gene_style") #I was getting weird error because some of the nodes werent added!
      #setNodeBorderWidthDefault(10, style.name = "gene_style")
      RCy3::setNodeWidthDefault(75, style.name = "gene_style")
      RCy3::setNodeHeightDefault(75, style.name = "gene_style")
      RCy3::setVisualStyle("gene_style", network = net)
      if(is.null(graph_layout)){
        layout_str = 'force-directed defaultSpringCoefficient=0.00000004 defaultSpringLength=100'
      }
      RCy3::layoutNetwork(layout_str)
  }
}



zoomPathway <- function(x, pway, graph_layout = NULL){
  if (!pway %in% x$graph$pathways$pathway) stop(paste0("Pathway: \"", pway, "\" not found in list of pathways"))
  edges_pathways_list <- makePathwayEdges(x$graph$edgelist, x$graph$pathways)
  edges_all           <- edges_pathways_list[["edges_all"]]
  pathway_gene_edges  <- edges_all[src_pathway ==  pway & dest_pathway == pway]
  pathway_graph       <- igraph::graph_from_data_frame(pathway_gene_edges[, .(base_gene_src, base_gene_dest)], directed = FALSE, vertices = unique(x$graph$pathways[pathway == pway,][["gene"]]))
  
  if(!is.null(graph_layout)) l <- graph_layout(pathway_graph)
  else                       l <- igraph::layout_with_graphopt(pathway_graph, spring.length = 300, spring.constant = 0.00000004)
  plot(pathway_graph, main = pway, layout = l)
}



# Helper functions -----------------------------------------------------------------
#Return just the edges between pathways & the edges with pathways merged on (for the nested network)
#Delete self edges
makePathwayEdges <- function(gene_edges, pathway_gene_map){
  res <- gene_edges[pathway_gene_map, .(frequency, base_gene_src, base_gene_dest, src_pathway = i.pathway), on = .(base_gene_src = gene), allow.cartesian = TRUE, nomatch = 0L][
    pathway_gene_map, .(frequency, base_gene_src, base_gene_dest, src_pathway, dest_pathway = i.pathway), on = .(base_gene_dest = gene), allow.cartesian = TRUE, nomatch = 0L]
  return(list(edges_all = res, edges_pathways = res[src_pathway != dest_pathway, .(weight_sum = sum(frequency)), by = .(src_pathway, dest_pathway)]))
}

createNestedNetwork <- function(edges_pathways, edges_all, pathway_vertices, main){
  
  #Copy so doesn't update outside of function
  edges_pathways   <- copy(edges_pathways)[, c("src_pathway", "dest_pathway") := lapply(.SD, function(x) gsub(" ", "\\ ", trimws(x), fixed = TRUE)), .SDcols = c("src_pathway", "dest_pathway")]
  edges_all        <- copy(edges_all)[, c("base_gene_src", "base_gene_dest", "src_pathway", "dest_pathway") := lapply(.SD, function(x) gsub(" ", "\\ ", trimws(x), fixed = TRUE)), .SDcols = c("base_gene_src", "base_gene_dest", "src_pathway", "dest_pathway")]
  pathway_vertices <- copy(pathway_vertices)[, c("pathway", "gene") := lapply(.SD, function(x) gsub(" ", "\\ ", trimws(x), fixed = TRUE))]
  
  pathways_edges_nnf     <- paste0(main, " ", edges_pathways$src_pathway, " pp ", edges_pathways$dest_pathway)
  no_edge_vertices       <- setdiff(unique(pathway_vertices$pathway), unique(c(edges_pathways$src_pathway, edges_pathways$dest_pathway)))
  no_edge_vertices_nnf   <- if (length(no_edge_vertices) == 0) NULL
                            else paste0(main, " ", no_edge_vertices)
  
  
  #Ones we want to link
  pathways_gene_vertices <- paste0(pathway_vertices$pathway, " ", pathway_vertices$gene)
  gene_edges_nnf_in      <- edges_all[src_pathway==dest_pathway]
  gene_edges_nnf_in2     <- paste0(gene_edges_nnf_in$src_pathway, " ", gene_edges_nnf_in$base_gene_src, " pp ", gene_edges_nnf_in$base_gene_dest)
  
  #Making NNF File
  all_nnf  <- c(pathways_edges_nnf, no_edge_vertices_nnf, pathways_gene_vertices, gene_edges_nnf_in2)
  temp_loc <- paste0(tempfile(),".nnf")
  f        <- base::file(temp_loc) 
  writeLines(all_nnf, f)
  close(f)
  
  nested_nets <- RCy3::importNetworkFromFile(temp_loc) ##First SUID returned if for first network (in this case one with pathways)
  unlink(temp_loc) #Delete file
  
  return(nested_nets)
}

