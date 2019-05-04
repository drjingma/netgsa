obtainEdgeList <- function(genes, databases){


  metabolites_ids   <- as.data.table(graphite:::metabolites())
  genes_with_id     <- reshapeIDs(genes, metabolites_ids)
  
  # Creating database here so we only need to process once. Also doing after genes_with_id since this takes the longest. Want to fail as early as possible
  full_edgelist     <- stackDatabases(databases)
  
  # Obtaining unique ids used by all chosen databases. We will convert to these IDs.
  # E.g. if databases only use "ENTREZID", "CHEBI", and "UNIPROT", this return a list of list(org.Hs.eg.db = c("ENTREZID", "UNIPROT"), metabolites = "CHEBI")
  databases_ids     <- findDatabaseIDs(full_edgelist, metabolites_ids)
  
  # Convert proteins
  ids_to_convert_using_orgHsegdb <- genes_with_id[conversion_type == "org.Hs.eg.db"]
  spl_ids                        <- split(ids_to_convert_using_orgHsegdb, ids_to_convert_using_orgHsegdb[["id_type"]])
  converted_proteins             <- rbindlist(lapply(spl_ids, convert_protein_groups, databases_ids[["org.Hs.eg.db"]]), use.names = TRUE)
  
  # Convert metabolites
  metabs_ids            <- genes_with_id[conversion_type == "graphite_metabolites"]
  converted_metabolites <- if(nrow(metabs_ids) == 0) NULL
                           else                      rbindlist(Map(convert_single_metabolite, metabs_ids[["id_type"]], metabs_ids[["gene"]], list(databases_ids[["metabolites"]]), list(metabolites_ids)), use.names = TRUE)
  
  # All finalized_ids should be character
  finalized_ids <- rbindlist(list(converted_proteins, converted_metabolites), use.names = TRUE, fill = TRUE)

  #Merge twice separately
  src_subs             <- setnames(full_edgelist[finalized_ids, on = .(src = converted_gene,  src_type  = converted_id), nomatch = 0L], c("base_gene", "base_id"), c("base_gene_src",  "base_id_src"))
  full_subs            <- setnames(     src_subs[finalized_ids, on = .(dest = converted_gene, dest_type = converted_id), nomatch = 0L], c("base_gene", "base_id"), c("base_gene_dest", "base_id_dest"))
  
  # Make all directed and stack on
  full_edgelist_subs_dir  <- rbindlist(list(full_subs, full_subs[direction == "undirected", .(src = dest,                     src_type = dest_type,       dest = src,                     dest_type = src_type, 
                                                                                              base_gene_src = base_gene_dest, base_id_src = base_id_dest, base_gene_dest = base_gene_src, base_id_dest = base_id_src,
                                                                                              database = database)]), use.names = TRUE, fill = TRUE)
  # Take unique
  # UNIT TEST IDEA: make sure when grouping by ID & database that there are no obs with > 1 occurrence. E.g. an edge shouldnt be in a database > 1 time.
  full_edgelist_subs_dir_uniq  <- unique(full_edgelist_subs_dir[, c("database", "base_gene_src", "base_id_src", "base_gene_dest", "base_id_dest"), with = FALSE])
  
  ## Lastly, we need to know which genes never appear in the full stack edgelist. If they dont appear then we likely dont have enough info
  ## to tell whether or not they should be 0 or 1.
  genes_not_in_dbs <- findGenesNotInDb(full_edgelist, finalized_ids)
  
  return(list(edgelist = full_edgelist_subs_dir_uniq, genes_not_in_dbs = genes_not_in_dbs))
}






# Wrapper to stack databases of interest 
  stackDatabases <- function(databases){
    
    not_in_graphite <- ! databases %in% unique(as.character(graphite::pathwayDatabases()[["database"]]))
    if(any(not_in_graphite)) stop(paste0("Database(s): ", paste0(databases[not_in_graphite], collapse = ", "), 
                                                          " are not found within graphite. Please choose databases from graphite::pathwayDatabases()$databases"))
    
    db_stacked <- rbindlist(lapply(databases, function(db){
                                              rbindlist(lapply(pathways("hsapiens", db), graphite::edges, which = "mixed", stringsAsFactors = FALSE), use.names = TRUE, idcol = "Pathway")}) 
                            %>% setNames(., databases), use.names = TRUE, idcol = "database")
    return(unique(db_stacked[,.(database, src, src_type, dest, dest_type, direction)]))
  }


# Reshape IDs into format that's easier to merge onto database
  reshapeIDs <- function(genes, metabolites){
    
    # Make sure all genes have names
    if(is.null(names(genes)) | any(names(genes) == "")) stop("Please select identfier type of all genes. You can do this by setting the name of each element to the identifier type of that gene")
    
    # Make sure their gene is gene type we know about
    unconv_ids <- !names(genes) %in% c(names(metabolites), AnnotationDbi::keytypes(org.Hs.eg.db))
    
    if(any(unconv_ids)) stop(paste0("Don't know how to convert ID type(s): ", paste(names(genes)[unconv_ids], collapse = ", ")))
  
    
    # All ID types are valid, drop duped IDs and make sure actual genes are valid for specific ID type
    genes_cl <- checkIDs(genes, metabolites)
    conversion_type = ifelse(names(genes_cl) %in% names(metabolites),                      "graphite_metabolites",
                      ifelse(names(genes_cl) %in% AnnotationDbi::keytypes(org.Hs.eg.db),   "org.Hs.eg.db", NA))
    
    data.table(id_type = names(genes_cl), gene = genes_cl, conversion_type = conversion_type)
    
  }

# Checking to make sure IDs are valid. That is, make sure we can find the IDs in our databases
  checkIDs <- function(genes, metabolites){
    
    # Error if cannot find IDs in database
    spl_genes <- split(genes, names(genes))
  
    not_found_ids <- lapply(spl_genes, function(gene_subs, metabolites){
                                              id_type       <- names(gene_subs)[1]
                                              valid_id_list <- if(      id_type %in% names(metabolites)    )                     metabolites[[id_type]]
                                                               else if( id_type %in% AnnotationDbi::keytypes(org.Hs.eg.db))      AnnotationDbi::keys(org.Hs.eg.db, id_type)
                                              
                                              not_found     <- gene_subs[! gene_subs %in% valid_id_list]
                                              if( length(not_found) == 0 )  return(NULL)
                                              else                          return(paste0(id_type, ": ", paste(not_found, collapse = ", ")))
                        }, metabolites = metabolites)
    
    error_ids <- Filter(function(x) !is.null(x), not_found_ids)
    
    if(length(error_ids) != 0){
      prefix <- "The following IDs were not found in the keylist and thus are not able to be converted: \n  "
      stop(paste0(prefix, do.call(paste0, list(error_ids, collapse = "\n  "))))
    }
    
    # Warn if duplicated ID and ID type
    duped_id <- duplicated(genes) & duplicated(names(genes))
    if(any(duped_id)) {stop(paste0("The following duplicate IDs were detected. Please remove and rerun: ", paste0(paste0(names(genes[duped_id]), " = ", genes[duped_id]), collapse = ", ")))}
    
    
    return(genes[!duped_id])
  }


# Finding set of ids used across chosen databases
  findDatabaseIDs <- function(databases, metabolites){
    uniq_ids                   <- databases[, unique(c(src_type, dest_type))]
    
    # IDs are converted either with metabolites from graphite OR org.Hs.eg.db
    metabolite_or_orgHsegdb_id <- ifelse(uniq_ids %in% AnnotationDbi::keytypes(org.Hs.eg.db), "org.Hs.eg.db",
                                  ifelse(uniq_ids %in% names(metabolites),                    "metabolites", NA))
    return(split(uniq_ids, metabolite_or_orgHsegdb_id))
  }


# Function to convert one set of proteins with a common base ID to all the ids we want to convert to
  convert_protein_groups <- function(id_table, convert_to){
  
    # Should only ever be of length one b/c we split by id_type
    keytype                                         <- unique(id_table[["id_type"]])
    
    # If keytype == convert_to simply return original data reshaped / renamed
    if(all(keytype == convert_to))                     return(data.table(base_gene = id_table[["gene"]], converted_gene = id_table[["gene"]], base_id = id_table[["id_type"]], converted_id = id_table[["id_type"]]))
    
    # We select c(keytype, convert_to) so that we can select base_id twice if its in list we want to convert to. This is possible b/c data.table allows duplicate column names
    converted                                       <- setDT(AnnotationDbi::select(org.Hs.eg.db, keys = id_table[["gene"]], columns = convert_to, keytype = keytype))[, c(keytype, convert_to), with = FALSE]
    converted                                       <- unique(data.table::melt(converted, id.vars = keytype, variable.name = "converted_id", value.name = "converted_gene", variable.factor = FALSE, value.factor = FALSE, na.rm = TRUE))
    converted[, c("base_id")]                       <- keytype
    names(converted)[names(converted) == keytype]   <- "base_gene"
    
    unconvertables                                  <- which(is.na(converted[["converted_gene"]]))
    if(length(unconvertables) != 0)                    warning(paste0("\nCould not convert the following genes :\n",
                                                                      paste0(paste0(keytype, " = ", converted[unconvertables, ][["base_gene"]], " to ", converted[unconvertables, ][["converted_id"]]), collapse = "\n")))
    
    return(converted)
  }


# Function to convert a single metabolite ID to all the ids we want to convert to 
  convert_single_metabolite <- function(base_id, base_gene, convert_to, metabolites){
  
    if(all(base_id == convert_to))                                           return(data.table(base_gene = base_gene, converted_id = base_id, converted_gene = base_gene, base_id = base_id))
    converted_ids_metabs                                                  <- unique(data.table::melt(metabolites[which(metabolites[[base_id]] == base_gene), c(base_id, convert_to), with = FALSE], id.vars = base_id, na.rm = TRUE, 
                                                                                                     variable.name = "converted_id", value.name = "converted_gene", variable.factor = FALSE, value.factor = FALSE))
  
    # Throwing warning about unconvertable IDs
    unconv_ids                                                            <- if(grepl("KEGG", base_id)) setdiff(convert_to[!grepl("KEGG", convert_to)], converted_ids_metabs[["converted_id"]])
                                                                             else                       setdiff(convert_to, converted_ids_metabs[["converted_id"]])
    if(length(unconv_ids != 0)){                                             warning(paste0("\nCould not convert gene :\n", base_id, " = ", base_gene, " to ", paste0(unconv_ids, collapse = ", ")))}
    
    if(nrow(converted_ids_metabs) == 0){                                     return(NULL)}
    
    names(converted_ids_metabs)[names(converted_ids_metabs) == base_id]   <- "base_gene"
    converted_ids_metabs[["base_id"]]                                     <- base_id
  
    return(converted_ids_metabs)
  }

#Given edgelist from stackDatabases and idlist from convertID functions, get the genes not in the database ever
  findGenesNotInDb <- function(edgelist, idlist){
    
    src_m  <- edgelist[idlist, on = .(src  = converted_gene, src_type  = converted_id), .(in_db = max(!is.na(dest)),  base_gene, base_id), nomatch = NA, by = .EACHI]
    dest_m <- edgelist[idlist, on = .(dest = converted_gene, dest_type = converted_id), .(in_db = max(!is.na(src)), base_gene, base_id), nomatch = NA, by = .EACHI]
    
    src_m2  <- src_m[, .(any_in_db_src = max(in_db)), by= .(base_gene, base_id)]
    dest_m2 <- dest_m[, .(any_in_db_dest = max(in_db)), by= .(base_gene, base_id)]
    
    not_in_dbs <- src_m2[dest_m2, on = c("base_gene", "base_id"), .(any_in_db = pmax(any_in_db_src, any_in_db_dest)), by = .EACHI][any_in_db == 0]
    
    if(nrow(not_in_dbs) > 0) return(setNames(not_in_dbs[["base_gene"]], not_in_dbs[["base_id"]]))
    else                     return(NULL)
    
  }



