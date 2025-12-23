## ============================================================
## Random Forest (ranger)
## with parallelized nrepeat (future.apply)
## ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(purrr)
  library(pROC)
  library(ranger)
  library(future)
  library(future.apply)
})

## ------------------------------------------------------------
## Utility: zip_nPure
## ------------------------------------------------------------
zip_nPure <- function(.x, .fields = NULL, .simplify = FALSE) {
  if (length(.x) == 0) return(list())
  if (is.null(.fields)) {
    if (is.null(names(.x[[1]]))) {
      .fields <- seq_along(.x[[1]])
    } else {
      .fields <- stats::setNames(names(.x[[1]]), names(.x[[1]]))
    }
  } else {
    if (is.character(.fields) && is.null(names(.fields))) names(.fields) <- .fields
  }
  out <- lapply(.fields, function(i) lapply(.x, .subset2, i))
  if (.simplify) out <- lapply(out, simplify_if_possible)
  out
}


# rf_cv <- function(X, Y, gse, sampletype, comparison,
#                   folds = 5, nrepeat = 20, seed = 1,
#                   importance = "impurity", num.trees = 1000) {
#   
#   # keep complete cases
#   keep <- complete.cases(X) & !is.na(Y)
#   X <- X[keep, , drop = FALSE]
#   Y <- droplevels(Y[keep])
#   
#   if (nrow(X) < 4 || ncol(X) < 2 || nlevels(Y) != 2) {
#     return(list(
#       df = data.frame(
#         value = NA_real_,
#         stat = c("AUC", "Accuracy"),
#         gse = gse,
#         sampletype = sampletype,
#         comparison = comparison
#       ),
#       feat = data.frame()
#     ))
#   }
#   
#   lvls <- levels(Y)
#   y_bin <- as.integer(Y == lvls[2])  # positive class = 2nd level
#   n <- nrow(X)
#   
#   repeat_out <- future.apply::future_lapply(seq_len(nrepeat), function(r) {
#     set.seed(seed + r)
#     
#     # stratified folds
#     idx1 <- which(Y == lvls[1])
#     idx2 <- which(Y == lvls[2])
#     f1 <- sample(rep_len(seq_len(folds), length(idx1)))
#     f2 <- sample(rep_len(seq_len(folds), length(idx2)))
#     fold_id <- integer(n)
#     fold_id[idx1] <- f1
#     fold_id[idx2] <- f2
#     
#     preds <- rep(NA_real_, n)
#     
#     for (k in seq_len(folds)) {
#       tr <- fold_id != k
#       te <- fold_id == k
#       
#       dat_tr <- data.frame(Y = Y[tr], X[tr, , drop = FALSE])
#       dat_te <- data.frame(Y = Y[te], X[te, , drop = FALSE])
#       
#       fit <- ranger::ranger(
#         Y ~ ., data = dat_tr,
#         probability = TRUE,
#         importance = importance,
#         num.trees = num.trees
#       )
#       
#       preds[te] <- predict(fit, data = dat_te)$predictions[, lvls[2]]
#     }
#     
#     auc <- as.numeric(pROC::auc(pROC::roc(y_bin, preds, quiet = TRUE)))
#     acc <- mean((preds >= 0.5) == (y_bin == 1), na.rm = TRUE)
#     
#     stats_df <- data.frame(
#       value = c(auc, acc),
#       stat  = c("AUC", "Accuracy"),
#       rep_id = r
#     )
#     
#     # importance on full data (per repeat; averaged later)
#     fit_full <- ranger::ranger(
#       Y ~ ., data = data.frame(Y = Y, X),
#       probability = TRUE,
#       importance = importance,
#       num.trees = num.trees
#     )
#     imp_vec <- ranger::importance(fit_full)
#     imp_df <- data.frame(
#       feature = names(imp_vec),
#       importance = as.numeric(imp_vec),
#       rep_id = r
#     )
#     
#     list(stats = stats_df, imp = imp_df)
#   }, future.seed = TRUE)
#   
#   all_stats <- lapply(repeat_out, `[[`, "stats")
#   all_imp   <- lapply(repeat_out, `[[`, "imp")
#   
#   df <- dplyr::bind_rows(all_stats) %>%
#     dplyr::group_by(stat) %>%
#     dplyr::summarise(AUC.mean = mean(value, na.rm = TRUE),
#                      AUC.sd = sd(value, na.rm = TRUE), 
#                      .groups = "drop") %>%
#     dplyr::mutate(
#       gse = gse,
#       sampletype = sampletype,
#       comparison = comparison
#     ) %>% 
#     dplyr::filter(stat == "AUC") %>% 
#     dplyr::select(-stat) %>% 
#     gather(stat, value, AUC.mean:AUC.sd) %>% 
#     dplyr::mutate(model = "RF")
#   
#   feat <- dplyr::bind_rows(all_imp) %>%
#     dplyr::group_by(feature) %>%
#     dplyr::summarise(value = mean(importance, na.rm = TRUE), .groups = "drop") %>%
#     dplyr::mutate(
#       gse = gse,
#       sampletype = sampletype,
#       comparison = comparison,
#       model = "RF"
#     ) %>%
#     dplyr::arrange(dplyr::desc(value))
#   
#   list(df = df, feat = feat)
# }

## ------------------------------------------------------------
## Main: lar_rf
## ------------------------------------------------------------
rf_cv <- function(X, Y, gse, sampletype, comparison,
                                    folds = 5, nrepeat = 20, seed = 1,
                                    importance = "impurity", num.trees = 1000) {
  # clean input
  keep <- complete.cases(X) & !is.na(Y)
  X <- X[keep, , drop = FALSE]
  Y <- droplevels(Y[keep])
  
  if (nrow(X) < 4 || ncol(X) < 2 || nlevels(Y) != 2) {
    return(list(
      df = data.frame(
        value = NA_real_,
        stat = c("AUC", "Accuracy"),
        gse = gse,
        sampletype = sampletype,
        comparison = comparison
      ),
      feat = data.frame()
    ))
  }
  
  lvls <- levels(Y)
  y_bin <- as.integer(Y == lvls[2])  # positive class = 2nd level
  n <- nrow(X)
  
  ## ---- parallelize repeats ----
  repeat_out <- future.apply::future_lapply(seq_len(nrepeat), function(r) {
    set.seed(seed + r)
    
    # stratified fold assignment
    idx1 <- which(Y == lvls[1])
    idx2 <- which(Y == lvls[2])
    f1 <- sample(rep_len(seq_len(folds), length(idx1)))
    f2 <- sample(rep_len(seq_len(folds), length(idx2)))
    fold_id <- integer(n)
    fold_id[idx1] <- f1
    fold_id[idx2] <- f2
    
    preds <- rep(NA_real_, n)
    
    for (k in seq_len(folds)) {
      tr <- fold_id != k
      te <- fold_id == k
      
      dat_tr <- data.frame(Y = Y[tr], X[tr, , drop = FALSE])
      dat_te <- data.frame(Y = Y[te], X[te, , drop = FALSE])
      
      fit <- ranger::ranger(
        Y ~ ., data = dat_tr,
        probability = TRUE,
        importance = importance,
        num.trees = num.trees
      )
      
      preds[te] <- predict(fit, data = dat_te)$predictions[, lvls[2]]
    }
    
    auc <- as.numeric(pROC::auc(pROC::roc(y_bin, preds, quiet = TRUE)))
    acc <- mean((preds >= 0.5) == (y_bin == 1), na.rm = TRUE)
    
    stats_df <- data.frame(
      value = c(auc, acc),
      stat  = c("AUC", "Accuracy"),
      rep_id = r
    )
    
    # variable importance from full-data fit (per repeat; averaged later)
    fit_full <- ranger::ranger(
      Y ~ ., data = data.frame(Y = Y, X),
      probability = TRUE,
      importance = importance,
      num.trees = num.trees
    )
    
    imp_vec <- ranger::importance(fit_full)
    imp_df <- data.frame(
      feature = names(imp_vec),
      importance = as.numeric(imp_vec),
      rep_id = r
    )
    
    list(stats = stats_df, imp = imp_df)
  }, future.seed = TRUE)
  ## -----------------------------
  
  all_stats <- lapply(repeat_out, `[[`, "stats")
  all_imp   <- lapply(repeat_out, `[[`, "imp")
  
  df <- dplyr::bind_rows(all_stats) %>%
    dplyr::group_by(stat) %>%
    dplyr::summarise(AUC.mean = mean(value, na.rm = TRUE),
                     AUC.sd = sd(value, na.rm = TRUE), 
                     .groups = "drop") %>%
    dplyr::mutate(
      gse = gse,
      sampletype = sampletype,
      comparison = comparison
    ) %>% 
    dplyr::filter(stat == "AUC") %>% 
    dplyr::select(-stat) %>% 
    gather(stat, value, AUC.mean:AUC.sd) %>% 
    dplyr::mutate(model = "RF")
  
  feat <- dplyr::bind_rows(all_imp) %>%
    dplyr::group_by(feature) %>%
    dplyr::summarise(importance = mean(importance, na.rm = TRUE), .groups = "drop") %>%
    dplyr::mutate(
      gse = gse,
      sampletype = sampletype,
      comparison = comparison,
      model = "RF"
    ) %>%
    dplyr::arrange(dplyr::desc(importance))
  
  list(df = df, feat = feat)
}
lar_rf <- function(eset_list, pheno_list, gse_list, biomarkers,
                   folds = 5,
                   nrepeat = 20,
                   seed = 1,
                   importance = "impurity",
                   num.trees = 1000,
                   workers = NULL) {
  
  ## Parallel plan (set workers=NULL to let R decide a default; or pass a number)
  ## If you already set a plan elsewhere, you can comment this out.
  if (!is.null(workers)) {
    future::plan(future::multisession, workers = workers)
  } else if (!inherits(future::plan(), "multisession") && !inherits(future::plan(), "multicore")) {
    future::plan(future::multisession)
  }
  
  res <- mapply(function(x, y, gse) {
    # align samples
    x <- x[, rownames(y)]
    group <- y$phenotype
    sampletype <- unique(y$sampletype)
    
    genes <- intersect(rownames(x), biomarkers)
    
    ## 1) healthyControls vs moderateAsthma (exclude severe)
    X1 <- t(x[genes, group != "severeAsthma", drop = FALSE])
    Y1 <- droplevels(group[group != "severeAsthma"])
    out1 <- rf_cv(X1, Y1, gse, sampletype, "healthyControls_vs_moderateAsthma")
    
    ## 2) healthyControls vs severeAsthma (exclude moderate)
    X2 <- t(x[genes, group != "moderateAsthma", drop = FALSE])
    Y2 <- droplevels(group[group != "moderateAsthma"])
    out2 <- rf_cv(X2, Y2, gse, sampletype, "healthyControls_vs_severeAsthma")
    
    ## 3) moderateAsthma vs severeAsthma (exclude healthy)
    X3 <- t(x[genes, group != "healthyControls", drop = FALSE])
    Y3 <- droplevels(group[group != "healthyControls"])
    out3 <- rf_cv(X3, Y3, gse, sampletype, "moderateAsthma_vs_severeAsthma")
    
    list(
      feat = dplyr::bind_rows(out1$feat, out2$feat, out3$feat),
      df   = dplyr::bind_rows(out1$df,   out2$df,   out3$df)
    )
    
  }, x = eset_list[gse_list],
  y = pheno_list[gse_list],
  gse = gse_list,
  SIMPLIFY = FALSE) %>%
    zip_nPure()
  
  return(res)
}

## ============================================================
## Example usage:
## future::plan(multisession, workers = 6) # optional; or pass workers=6 to lar_rf()
## res <- lar_rf(eset_list, pheno_list, gse_list, biomarkers, nrepeat=20, folds=5, workers=6)
## res$df   # list of per-GSE df tables
## res$feat # list of per-GSE feature-importance tables
## ============================================================

# future::plan(multisession, workers = 6)
# all_sev <- lapply(names(latephase_biomarkers)[1:2], function(i){
#   res <- lar_rf(eset_list = esetList, 
#                 pheno_list = phenoDataList, 
#                 gse_list = select_gse, 
#                 biomarkers = latephase_biomarkers[[i]],
#                 folds = 2, nrepeat = 2, seed = 1,
#                 importance = "impurity")
#   list(feat = do.call(rbind, res$feat) %>% mutate(panel = i),
#        df = do.call(rbind, res$df) %>% mutate(panel = i))
# }) %>%
#   zip_nPure()

## ------------------------------------------------------------
## multilevel_rf: replaces multilevel_plsda
## ------------------------------------------------------------
multilevel_rf <- function(eset, demo, biomarkers,
                          folds = 5, nrepeat = 20, seed = 1,
                          importance = "impurity", num.trees = 1000,
                          workers = NULL,
                          gse = "GSE19301",
                          sampletype = "blood") {
  
  suppressPackageStartupMessages(library(mixOmics))
  
  ## parallel plan
  if (!is.null(workers)) {
    future::plan(future::multisession, workers = workers)
  } else if (!inherits(future::plan(), "multisession") &&
             !inherits(future::plan(), "multicore")) {
    future::plan(future::multisession)
  }
  
  ## keep only donors with all phenotypes
  mat <- table(demo$donor, demo$phenotype)
  keepdonors <- names(which(rowSums(mat == 0) == 0))
  demo <- demo[demo$donor %in% keepdonors, , drop = FALSE]
  eset <- eset[, rownames(demo), drop = FALSE]
  
  genes <- intersect(rownames(eset), biomarkers)
  
  ## =========================================================
  ## quiet vs exacerbation (exclude FOLLOW UP)
  ## =========================================================
  X <- t(eset[genes, demo$phenotype != "FOLLOW UP", drop = FALSE])
  Y <- demo$phenotype
  Y <- factor(Y[Y != "FOLLOW UP"], levels = c("QUIET", "EXACERBATION"))
  
  design <- data.frame(sample = demo$donor[demo$phenotype != "FOLLOW UP"])
  
  ## EXACT same multilevel transform as PLS-DA
  Xw <- mixOmics::withinVariation(X = X, design = design)
  
  out1 <- rf_cv(
    X = Xw,
    Y = Y,
    gse = gse,
    sampletype = sampletype,
    comparison = "Exacer vs. Quiet",
    folds = folds,
    nrepeat = nrepeat,
    seed = seed,
    importance = importance,
    num.trees = num.trees
  )
  
  ## =========================================================
  ## exacerbation vs follow-up (exclude QUIET)
  ## =========================================================
  X <- t(eset[genes, demo$phenotype != "QUIET", drop = FALSE])
  Y <- demo$phenotype
  Y <- factor(Y[Y != "QUIET"], levels = c("FOLLOW UP", "EXACERBATION"))
  
  design <- data.frame(sample = demo$donor[demo$phenotype != "QUIET"])
  
  Xw <- mixOmics::withinVariation(X = X, design = design)
  
  out2 <- rf_cv(
    X = Xw,
    Y = Y,
    gse = gse,
    sampletype = sampletype,
    comparison = "Exacer vs. followup",
    folds = folds,
    nrepeat = nrepeat,
    seed = seed + 10000,
    importance = importance,
    num.trees = num.trees
  )
  
  list(
    feat = dplyr::bind_rows(out1$feat, out2$feat),
    df   = dplyr::bind_rows(out1$df,   out2$df)
  )
}

# res <- multilevel_rf(eset = eset, demo = demo, biomarkers = biomarkers,
#   folds = 2, nrepeat = 2, seed = 1,
#   importance = "impurity", num.trees = 1000,
#   workers = NULL,
#   gse = "GSE19301",
#   sampletype = "blood")
# 
# latephase_biomarkers <- readRDS(here::here("results", "latephase_biomarkers.rds"))
# latephase_biomarkers$all <- latephase_exacer_biomarkers
# exacer_plsda <- lapply(names(latephase_biomarkers), function(i){
#   biomarkers <- latephase_biomarkers[[i]]
#   ## all subjects
#   demo <- subset(phenoDataList$GSE19301, severity == "Moderate Persistent")
#   eset <- esetList$GSE19301[intersect(rownames(esetList$GSE19301), biomarkers), rownames(demo)]
#   exacer_all_mod <- multilevel_rf(eset = eset, demo = demo,
#                                      biomarkers = biomarkers)
#   exacer_all_mod$df$comparison <- exacer_all_mod$feat$comparison <- paste(exacer_all_mod$df$comp, "(Moderate Persistent)")
#   exacer_all_mod$df$sex <- exacer_all_mod$feat$sex <- "M+F"
#   demo <- subset(phenoDataList$GSE19301, severity == "Severe Persistent")
#   eset <- esetList$GSE19301[intersect(rownames(esetList$GSE19301), biomarkers), rownames(demo)]
#   exacer_all_sev <- multilevel_plsda(eset = eset, demo = demo,
#                                      biomarkers = biomarkers)
#   exacer_all_sev$df$comparison <- exacer_all_sev$feat$comparison <- paste(exacer_all_sev$df$comp, "(Severe Persistent)")
#   exacer_all_sev$df$sex <- exacer_all_sev$feat$sex <- "M+F"
#   
#   ## females
#   demo <- subset(phenoDataList$GSE19301, severity == "Moderate Persistent" & sex ==  "F")
#   eset <- esetList$GSE19301[intersect(rownames(esetList$GSE19301), biomarkers), rownames(demo)]
#   exacer_f_mod <- multilevel_plsda(eset = eset, demo = demo,
#                                    biomarkers = biomarkers)
#   exacer_f_mod$df$comparison <- exacer_f_mod$feat$comparison <- paste(exacer_f_mod$df$comp, "(Moderate Persistent)")
#   exacer_f_mod$df$sex <- exacer_f_mod$feat$sex <- "F"
#   demo <- subset(phenoDataList$GSE19301, severity == "Severe Persistent" & sex ==  "F")
#   eset <- esetList$GSE19301[intersect(rownames(esetList$GSE19301), biomarkers), rownames(demo)]
#   exacer_f_sev <- multilevel_plsda(eset = eset, demo = demo,
#                                    biomarkers = biomarkers)
#   exacer_f_sev$df$comparison <- exacer_f_sev$feat$comparison <- paste(exacer_f_sev$df$comp, "(Severe Persistent)")
#   exacer_f_sev$df$sex <- exacer_f_sev$feat$sex <- "F"
#   
#   ## males
#   demo <- subset(phenoDataList$GSE19301, severity == "Moderate Persistent" & sex ==  "M")
#   eset <- esetList$GSE19301[intersect(rownames(esetList$GSE19301), biomarkers), rownames(demo)]
#   exacer_m_mod <- multilevel_plsda(eset = eset, demo = demo,
#                                    biomarkers = biomarkers)
#   exacer_m_mod$df$comparison <- exacer_m_mod$feat$comparison <- paste(exacer_m_mod$df$comp, "(Moderate Persistent)")
#   exacer_m_mod$df$sex <- exacer_m_mod$feat$sex <- "M"
#   demo <- subset(phenoDataList$GSE19301, severity == "Severe Persistent" & sex ==  "M")
#   eset <- esetList$GSE19301[intersect(rownames(esetList$GSE19301), biomarkers), rownames(demo)]
#   exacer_m_sev <- multilevel_plsda(eset = eset, demo = demo,
#                                    biomarkers = biomarkers)
#   exacer_m_sev$df$comparison <- exacer_m_sev$feat$comparison <- paste(exacer_m_sev$df$comp, "(Severe Persistent)")
#   exacer_m_sev$df$sex <- exacer_m_sev$feat$sex <- "M"
#   
#   exacer_df <- rbind(exacer_all_mod$df, exacer_all_sev$df,
#                      exacer_f_mod$df, exacer_f_sev$df,
#                      exacer_m_mod$df, exacer_m_sev$df)
#   exacer_feat <- rbind(exacer_all_mod$feat, exacer_all_sev$feat,
#                        exacer_f_mod$feat, exacer_f_sev$feat,
#                        exacer_m_mod$feat, exacer_m_sev$feat)
#   exacer_df$panel <- exacer_feat$panel <- i
#   list(df = exacer_df, feat = exacer_feat)
# }) %>%
#   zip_nPure()
# saveRDS(exacer_plsda, here("results", "exacer_plsda.rds"))