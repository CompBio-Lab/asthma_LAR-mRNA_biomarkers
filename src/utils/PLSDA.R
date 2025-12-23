suppressPackageStartupMessages({
  library(dplyr)
  library(mixOmics)
})

## -----------------------------
## small helpers
## -----------------------------
zip_nPure <- function(.x, .fields = NULL, .simplify = FALSE) {
  if (length(.x) == 0) return(list())
  if (is.null(.fields)) {
    .fields <- if (is.null(names(.x[[1]]))) seq_along(.x[[1]])
    else stats::setNames(names(.x[[1]]), names(.x[[1]]))
  } else if (is.character(.fields) && is.null(names(.fields))) {
    names(.fields) <- .fields
  }
  out <- lapply(.fields, function(i) lapply(.x, .subset2, i))
  if (.simplify) out <- lapply(out, simplify_if_possible)
  out
}

plsda_block <- function(x, group, biomarkers, keep_mask, y_levels, gse, sampletype, comparison,
                        folds = 5, nrepeat = 20, ncomp = 2) {
  
  genes <- intersect(rownames(x), biomarkers)
  X <- t(x[genes, keep_mask, drop = FALSE])
  Y <- factor(group[keep_mask], levels = y_levels)
  
  result <- mixOmics::plsda(X = X, Y = Y, ncomp = ncomp)
  BPPARAM <- BiocParallel::MulticoreParam(workers = folds)
  cv <- mixOmics::perf(result, validation = "Mfold", folds = folds,
                       nrepeat = nrepeat, auc = TRUE, BPPARAM = BPPARAM)
  
  df <- data.frame(
    value = cv$auc[[paste0("comp", ncomp)]],
    stat  = names(cv$auc[[paste0("comp", ncomp)]]),
    gse = gse,
    sampletype = sampletype,
    comparison = comparison,
    model = "PLSDA"
  )
  
  feat <- data.frame(importance = mixOmics::vip(result)[, ncomp]) %>% 
    mutate(feature = rownames(.),
           gse = gse,
           sampletype = sampletype,
           comparison = comparison,
           model = "PLSDA")
  
  list(df = df, feat = feat)
}

multilevel_plsda_block <- function(eset, demo, biomarkers, exclude_level,
                                   y_levels, gse, sampletype, comparison,
                                   folds = 5, nrepeat = 20, ncomp = 2) {
  
  genes <- intersect(rownames(eset), biomarkers)
  keep <- demo$phenotype != exclude_level
  
  X <- t(eset[genes, keep, drop = FALSE])
  Y <- factor(demo$phenotype[keep], levels = y_levels)
  design <- data.frame(sample = demo$donor[keep])
  
  Xw <- mixOmics::withinVariation(X = X, design = design)
  
  result <- mixOmics::plsda(X = Xw, Y = Y, ncomp = ncomp)
  BPPARAM <- BiocParallel::MulticoreParam(workers = folds)
  cv <- mixOmics::perf(result, validation = "Mfold", folds = folds, 
                       nrepeat = nrepeat, auc = TRUE, BPPARAM = BPPARAM)
  
  df <- data.frame(
    value = cv$auc[[paste0("comp", ncomp)]],
    stat  = names(cv$auc[[paste0("comp", ncomp)]]),
    gse = gse,
    sampletype = sampletype,
    comparison = comparison,
    model = "PLSDA"
  )
  
  feat <- data.frame(importance = mixOmics::vip(result)[, ncomp]) %>% 
    mutate(feature = rownames(.),
                     gse = gse,
                     sampletype = sampletype,
                     comparison = comparison,
                     model = "PLSDA")
  
  list(df = df, feat = feat)
}

## -----------------------------
## streamlined multilevel_plsda
## -----------------------------
multilevel_plsda <- function(eset, demo, biomarkers,
                             folds = 5, nrepeat = 20, ncomp = 2,
                             gse = "GSE19301", sampletype = "blood") {
  
  mat <- table(demo$donor, demo$phenotype)
  keepdonors <- names(which(rowSums(mat == 0) == 0))
  demo <- demo[demo$donor %in% keepdonors, , drop = FALSE]
  eset <- eset[, rownames(demo), drop = FALSE]
  
  blocks <- list(
    list(exclude = "FOLLOW UP", levels = c("QUIET", "EXACERBATION"),
         comparison = "Exacer vs. Quiet"),
    list(exclude = "QUIET", levels = c("FOLLOW UP", "EXACERBATION"),
         comparison = "Exacer vs. followup")
  )
  
  out <- lapply(blocks, function(b) {
    multilevel_plsda_block(
      eset = eset, demo = demo, biomarkers = biomarkers,
      exclude_level = b$exclude, y_levels = b$levels,
      gse = gse, sampletype = sampletype, comparison = b$comparison,
      folds = folds, nrepeat = nrepeat, ncomp = ncomp
    )
  })
  
  list(
    feat = bind_rows(lapply(out, `[[`, "feat")),
    df   = bind_rows(lapply(out, `[[`, "df"))
  )
}

## -----------------------------
## streamlined lar_plsda
## -----------------------------
lar_plsda <- function(eset_list, pheno_list, gse_list, biomarkers,
                      folds = 5, nrepeat = 20, ncomp = 2) {
  
  res <- mapply(function(x, y, gse) {
    x <- x[, rownames(y), drop = FALSE]
    group <- y$phenotype
    sampletype <- unique(y$sampletype)
    
    blocks <- list(
      list(keep = group != "severeAsthma",
           levels = c("healthyControls", "moderateAsthma"),
           comparison = "healthyControls_vs_moderateAsthma"),
      list(keep = group != "moderateAsthma",
           levels = c("healthyControls", "severeAsthma"),
           comparison = "healthyControls_vs_severeAsthma"),
      list(keep = group != "healthyControls",
           levels = c("moderateAsthma", "severeAsthma"),
           comparison = "moderateAsthma_vs_severeAsthma")
    )
    
    out <- lapply(blocks, function(b) {
      plsda_block(
        x = x, group = group, biomarkers = biomarkers,
        keep_mask = b$keep, y_levels = b$levels,
        gse = gse, sampletype = sampletype, comparison = b$comparison,
        folds = folds, nrepeat = nrepeat, ncomp = ncomp
      )
    })
    
    list(
      feat = bind_rows(lapply(out, `[[`, "feat")),
      df   = bind_rows(lapply(out, `[[`, "df"))
    )
  }, x = eset_list[gse_list], y = pheno_list[gse_list], gse = gse_list, SIMPLIFY = FALSE) %>%
    zip_nPure()
  
  res
}
