sc_drugde = function(sce){
  metadata <- colData(sce)
  metadata$sm_name <- relevel(metadata$sm_name, ref = "Dimethyl Sulfoxide")
  
  sce <- computeLibraryFactors(sce)
  sce <- logNormCounts(sce)
  
  df <- t(sce@assays@data$logcounts) %>%
    as.data.frame() %>%
    mutate(cell_type_orig = metadata$cell_type_orig,
           sm_name = metadata$sm_name,
           donor_id = metadata$donor_id,
           plate = metadata$plate_name,
           well = metadata$well) %>%
    dplyr::group_by(cell_type_orig, plate, well, sm_name, donor_id) %>%
    summarise_each(funs(mean))
  
  cells <- unique(as.character(df$cell_type_orig))
  res <- lapply(cells, function(cc){
    . <- subset(df, cell_type_orig == cc)
    block <- .$donor_id
    design <- model.matrix(~sm_name, data = .)
    eset <- t(.[, -c(1:5)])
    fit <- lmFit(eset, design)
    fit <- eBayes(fit, trend = TRUE)
    top_all <- lapply(colnames(design)[2:144], function(i){
      top <- topTable(fit, coef = i, adjust.method = "BH", n = nrow(fit)) %>%
        mutate(gene = rownames(.))
      top$sm_name <- gsub("sm_name", "", i)
      top$cc <- cc
      top
    }) %>%
      do.call(rbind, .)
    list(top_all=top_all, normE = eset)
  }) %>%
    zip_nPure()
  
  toptables <- res$top_all %>% do.call(rbind, .)
  normEs <- res$normE
  names(normEs) <- cells
  list(toptables=toptables, normEs=normEs)
}