test_biomarkers = function(eset,pheno,biomarkers){
  group <- pheno$phenotype
  sex <- pheno$sex
  design <- model.matrix(~group+sex+0)
  colnames(design) <- c("control", "mod", "sev", "sex")
  lmfit <- lmFit(eset[intersect(rownames(eset), biomarkers), ], design)
  cont <- makeContrasts(sev-control, mod-control, sev-mod, levels = design)
  lmfit.cont <- contrasts.fit(lmfit, cont)
  lmfit.cont.ebayes <- eBayes(lmfit.cont)
  top <- lapply(colnames(cont), function(contrast){
    topTable(lmfit.cont.ebayes, coef = contrast,
             adjust.method = "BH", n= nrow(lmfit.cont.ebayes)) %>% 
      as.data.frame() %>% 
      mutate(contrast = contrast) %>% 
      mutate(feature = rownames(.)) %>% 
      mutate(n = 1:n())
  }) %>% 
    do.call(rbind, .)
}
