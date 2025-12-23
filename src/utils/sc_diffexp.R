sc_diffexp = function(counts, metadata, cell_ann_col = "cluster"){
  # create singleCellExperiment object
  sce <- SingleCellExperiment(
    assays = list(counts = counts[, rownames(metadata)]),
    colData = metadata
  )
  
  sce <- computeLibraryFactors(sce)
  sce <- logNormCounts(sce)
  
  (sce <- prepSCE(sce,
                  kid = cell_ann_col, # subpopulation assignments
                  gid = "condition",  # group IDs (ctrl/stim)
                  sid = "sample_id", # sample IDs (ctrl/stim.1234)
                  drop = FALSE))  # drop all other colData columns
  
  # store cluster + sample IDs as well as number of clusters and samples
  nk <- length(kids <- levels(sce$cluster_id))
  ns <- length(sids <- levels(sce$sample_id))
  names(kids) <- kids; names(sids) <- sids
  
  # nb. of cells per cluster-sample
  df <- t(table(sce$cluster_id, sce$sample_id))
  df
  
  # aggregation of single-cell to pseudo-bulk data
  pb <- aggregateData(sce,
                      assay = "logcounts", fun = "mean",
                      by = c("cluster_id", "sample_id"))
  
  # construct design & contrast matrix
  ei <- metadata(sce)$experiment_info
  ei$subj <- sapply(strsplit(as.character(ei$sample_id), "_"), function(i) i[[1]])
  ei$group_id <- relevel(droplevels(ei$group_id), ref = "Bln")
  mm <- model.matrix(~ 0 + group_id, data = ei)
  # dimnames(mm) <- list(ei$sample_id, colnames(mm))
  dimnames(mm) <- list(ei$sample_id, levels(ei$group_id))
  
  contrast <- makeContrasts(Ag-Bln, levels = mm)
  
  # run DS analysis
  res <- pbDS(pb, design = mm, contrast = contrast, method="limma-trend",
              min_cells = 3)
  
  # access results table
  tbl <- lapply(res$table, function(i){
    do.call(rbind, i)
  }) %>%
    do.call(rbind, .)
  
  tbl2 <- tbl %>%
    filter(gene %in% latephase_exacer_biomarkers) %>%
    group_by(contrast, cluster_id) %>%
    mutate(p_adj.loc = p.adjust(p_val, "BH")) %>%
    arrange(p_val) %>%
    mutate(n = 1:n())
  # tbl2$contrast <- gsub("group_id", "", tbl2$contrast)
  list(toptables=tbl2)
}