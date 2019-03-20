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
        paths<- lapply(paths, function(p) convertIdentifiers(p,g.id))
      }
      genesets <- lapply(paths, nodes)
      genesets <- lapply(genesets, function(a) gsub("ENTREZID:","",a))
    } else {
      m_df = msigdbr::msigdbr(species = "Homo sapiens", category = m.type)
      m_t2g = m_df %>% dplyr::select(.data$gs_name, .data$entrez_gene) %>% as.data.frame()
      if (g.id=="symbol"){
        m_t2g = m_df %>% dplyr::select(.data$gs_name, .data$gene_symbol) %>% as.data.frame()
      }
      gs_df <- m_t2g %>% dplyr::group_by(.data$gs_name) %>% 
        dplyr::summarise(entrez_gene = paste0(.data$entrez_gene, collapse=","))
      genesets <- lapply(gs_df$entrez_gene, function(s) unlist(strsplit(s,",")))
      names(genesets) <- gs_df$gs_name
    }
    
    return(genesets)
  }
