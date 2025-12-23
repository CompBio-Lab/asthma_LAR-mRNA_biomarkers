#' Donor-adjusted (LME) pairwise phenotype p-values within each sex + plot
#'
#' @param gene Character. Gene symbol present as a column in geoPanelGenes.
#' @param severity Character. Severity level present in phenoDataList$GSE19301$severity.
#' @param geoPanelGenes Matrix/data.frame with samples in rows; must include gene, phenotype, sex, donor columns.
#' @param phenoDataList List containing phenoDataList$GSE19301 with a `severity` vector aligned to geoPanelGenes rows.
#' @param gse Character. Name of GSE within phenoDataList (default "GSE19301").
#' @param phenotype_levels Character vector specifying order of phenotype factor.
#' @param adjust_method P-value adjustment method for emmeans::pairs (default "BH").
#' @param facet_scales Facet scales ("fixed" or "free_y"). Default "fixed".
#' @param bracket_step Fraction of y-range used to stack p-value brackets.
#' @return A list with: plot (ggplot), data (df), stats (stat_tbl)
plot_gene_sex_lme <- function(
    gene,
    severity,
    geoPanelGenes,
    phenoDataList,
    gse = "GSE19301",
    phenotype_levels = c("QUIET", "EXACERBATION", "FOLLOW UP"),
    adjust_method = "BH",
    facet_scales = "fixed",
    bracket_step = 0.08
){
  stopifnot(is.character(gene), length(gene) == 1)
  stopifnot(is.character(severity), length(severity) == 1)
  stopifnot(gse %in% names(phenoDataList))
  
  # Required columns
  req_cols <- c("phenotype", "sex", "donor")
  missing_cols <- setdiff(req_cols, colnames(geoPanelGenes))
  if (length(missing_cols) > 0) {
    stop("geoPanelGenes is missing required columns: ", paste(missing_cols, collapse = ", "))
  }
  if (!(gene %in% colnames(geoPanelGenes))) {
    stop("Gene '", gene, "' not found in geoPanelGenes columns.")
  }
  if (!("severity" %in% colnames(phenoDataList[[gse]]))) {
    stop("phenoDataList[['", gse, "']] must contain a 'severity' column.")
  }
  
  # Subset by severity (assumes row alignment between phenoDataList[[gse]] and geoPanelGenes)
  idx <- phenoDataList[[gse]]$severity == severity
  if (sum(idx, na.rm = TRUE) == 0) {
    stop("No samples found for severity = '", severity, "'.")
  }
  
  # Build df
  df <- data.frame(
    expr      = geoPanelGenes[idx, gene],
    phenotype = geoPanelGenes[idx, "phenotype"],
    sex       = geoPanelGenes[idx, "sex"],
    donor     = geoPanelGenes[idx, "donor"],
    stringsAsFactors = FALSE
  ) %>%
    dplyr::mutate(
      phenotype = factor(phenotype, levels = phenotype_levels),
      sex = factor(sex),
      donor = factor(donor)
    ) %>%
    dplyr::filter(!is.na(expr), !is.na(phenotype), !is.na(sex), !is.na(donor))
  
  if (nrow(df) == 0) stop("After filtering NAs, no rows remain to plot/model.")
  
  # Helper: phenotype pairwise tests within one sex (random intercept donor)
  get_sex_pairs <- function(sex_level) {
    d <- df %>% dplyr::filter(sex == sex_level)
    
    # Need >= 2 phenotype levels with data
    if (nrow(d) < 3) return(NULL)
    if (dplyr::n_distinct(d$phenotype) < 2) return(NULL)
    
    fit <- lmerTest::lmer(expr ~ phenotype + (1 | donor), data = d, REML = FALSE)
    
    emm <- emmeans::emmeans(fit, ~ phenotype)
    
    pw <- emmeans::contrast(
      emm,
      method = "pairwise",
      adjust = adjust_method
    ) %>% as.data.frame()
    
    out <- pw %>%
      dplyr::transmute(
        sex = sex_level,
        group1 = gsub(" -.*$", "", contrast),
        group2 = gsub("^.*- ", "", contrast),
        p = p.value
      )
    
    # y positions for brackets
    y0 <- max(d$expr, na.rm = TRUE)
    rng <- diff(range(d$expr, na.rm = TRUE))
    step <- bracket_step * rng
    if (!is.finite(step) || step == 0) step <- 0.1
    
    out$y.position <- y0 + seq_len(nrow(out)) * step
    
    out$p.signif <- as.character(
      symnum(
        out$p,
        cutpoints = c(0, 1e-4, 1e-3, 1e-2, 0.05, 1),
        symbols   = c("****", "***", "**", "*", "ns")
      )
    )
    out
  }
  
  sex_levels <- levels(df$sex)
  stat_tbl <- dplyr::bind_rows(lapply(sex_levels, get_sex_pairs)) %>%
    dplyr::filter(!is.na(p))
  
  # Plot
  p <- ggplot2::ggplot(df, ggplot2::aes(x = phenotype, y = expr)) +
    ggplot2::geom_boxplot(
      ggplot2::aes(fill = sex),
      alpha = 0.18,
      outlier.shape = NA
    ) +
    ggplot2::geom_line(
      ggplot2::aes(group = donor, color = sex),
      alpha = 0.35
    ) +
    ggplot2::geom_point(
      ggplot2::aes(color = sex),
      size = 2.2, alpha = 0.85,
      position = ggplot2::position_jitter(width = 0.10, height = 0)
    ) +
    ggplot2::facet_wrap(~ sex, scales = facet_scales) +
    (if (nrow(stat_tbl) > 0)
      ggpubr::stat_pvalue_manual(
        stat_tbl,
        label = "p.signif",
        xmin = "group1", xmax = "group2",
        y.position = "y.position",
        tip.length = 0.01
      ) else ggplot2::geom_blank()) +
    ggplot2::labs(
      x = "Phenotype",
      y = gene,
      title = paste0(gene, " (", severity, ")")
    ) +
    ggplot2::theme_classic(base_size = 12) +
    ggplot2::theme(
      legend.position = "none",
      axis.text.x = ggplot2::element_text(angle = 30, hjust = 1)
    )
  
  list(plot = p, data = df, stats = stat_tbl)
}
