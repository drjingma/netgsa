preparePathways <-
  function(db=c("kegg", "MSigDB"), 
           type=c("H","C1","C2","C3","C4","C5","C6","C7"),
           genename= c("EntrezID", "symbol")) {
    this.call <- match.call()
    db <- match.arg(db)
    m.type <- match.arg(type)
    g.id <- match.arg(genename)
    
    if (db=="kegg" && is.null(g.id)){g.id="EntrezID"}
    
    if (db=="MSigDB" && is.null(m.type)){m.type="C2";g.id=="symbol"} 

    if (db=="kegg"){
      paths <- graphite::pathways('hsapiens','kegg')
      if (g.id=="symbol"){
        paths<- lapply(paths, function(p) graphite::convertIdentifiers(p,g.id))
      }
      genesets <- lapply(paths, nodes)
    } else {
      ## these gene sets are defined based on gene symbols;
      path.id <- grep(m.type,names(MSigDB))
      genesets <- MSigDB[[path.id]] 
    }
    
    return(genesets)
  }
