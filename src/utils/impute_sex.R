#' Impute biological sex using Y-chromosome gene expression
#'
#' This function infers biological sex from gene-expression data by identifying
#' Y-chromosome–encoded genes, performing dimensionality reduction, and clustering
#' samples into two groups corresponding to male and female expression profiles.
#'
#' Sex is inferred by:
#' \enumerate{
#'   \item Mapping input gene symbols to Ensembl and identifying Y-chromosome genes
#'   \item Performing PCA on scaled Y-chromosome expression
#'   \item Clustering samples into two groups using k-means
#'   \item Assigning the cluster with higher average Y-gene expression as male
#' }
#'
#' If known biological sex labels are provided, classification accuracy is also returned.
#'
#' @param eset A numeric gene-expression matrix or data frame with genes as rows
#'   and samples as columns. Row names must correspond to HGNC gene symbols.
#'
#' @param true_sex Optional character vector of true biological sex labels
#'   (values must be \code{"F"} or \code{"M"}). Must be named or ordered to match
#'   the columns of \code{eset}.
#'
#' @return A list with two elements:
#' \describe{
#'   \item{cluster}{A character vector of imputed sex labels (\code{"F"} or \code{"M"})
#'   for each sample. Names correspond to sample IDs.}
#'   \item{acc}{Classification accuracy if \code{true_sex} is provided; otherwise \code{NA}.}
#' }
#'
#' @details
#' This method assumes that Y-chromosome genes are reliably detected in the dataset
#' and that samples contain sufficient expression signal for sex discrimination.
#' Performance may be reduced in datasets with:
#' @examples
#' \dontrun{
#' sex_pred <- impute_sex(eset)
#'
#' sex_pred <- impute_sex(
#'   eset = expr_matrix,
#'   true_sex = metadata$sex
#' )
#'
#' sex_pred$acc
#' }
#'
#' @export
impute_sex <- function(eset, true_sex = NULL,
                       mart = NULL,
                       use_pcs = 2,
                       centers = 2,
                       nstart = 25,
                       min_y_genes = 3,
                       verbose = FALSE) {
  # ---- checks ----
  if (is.null(rownames(eset)) || is.null(colnames(eset))) {
    stop("`eset` must have rownames (HGNC symbols) and colnames (sample IDs).")
  }
  eset <- as.matrix(eset)
  storage.mode(eset) <- "double"
  
  # Handle true_sex alignment (either named vector or in same order as samples)
  if (!is.null(true_sex)) {
    if (!all(true_sex %in% c("F", "M"))) {
      stop("`true_sex` must contain only 'F' and 'M'.")
    }
    if (!is.null(names(true_sex))) {
      true_sex <- true_sex[colnames(eset)]
    } else {
      if (length(true_sex) != ncol(eset)) {
        stop("`true_sex` must be length ncol(eset) or a named vector matching colnames(eset).")
      }
    }
  }
  
  # ---- biomaRt lookup (allow reusing mart object) ----
  if (is.null(mart)) {
    mart <- biomaRt::useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  }
  
  genes_ychr <- biomaRt::getBM(
    attributes = c("chromosome_name", "hgnc_symbol"),
    filters    = "hgnc_symbol",
    values     = rownames(eset),
    mart       = mart
  )
  
  genes_ychr <- genes_ychr |>
    dplyr::filter(.data$chromosome_name == "Y") |>
    dplyr::distinct(.data$hgnc_symbol)
  
  ygenes <- intersect(genes_ychr$hgnc_symbol, rownames(eset))
  
  if (length(ygenes) < min_y_genes) {
    stop(sprintf(
      "Found %d Y-chromosome gene(s) in `eset` (need at least %d).",
      length(ygenes), min_y_genes
    ))
  }
  if (verbose) message("Using ", length(ygenes), " Y-chromosome gene(s) for sex imputation.")
  
  # ---- build sample x gene matrix ----
  sex_eset <- t(eset[ygenes, , drop = FALSE])  # samples x genes
  
  # Remove near-zero variance genes to avoid PCA issues
  vars <- apply(sex_eset, 2, stats::var, na.rm = TRUE)
  keep <- is.finite(vars) & vars > 0
  sex_eset <- sex_eset[, keep, drop = FALSE]
  if (ncol(sex_eset) < min_y_genes) {
    stop("After removing zero-variance Y genes, too few genes remain for PCA.")
  }
  
  # Impute missing values (simple median per gene) to keep PCA stable
  if (anyNA(sex_eset)) {
    for (j in seq_len(ncol(sex_eset))) {
      if (anyNA(sex_eset[, j])) {
        sex_eset[is.na(sex_eset[, j]), j] <- stats::median(sex_eset[, j], na.rm = TRUE)
      }
    }
  }
  
  # ---- PCA + kmeans ----
  pca <- stats::prcomp(sex_eset, scale. = TRUE, center = TRUE)
  
  k <- min(use_pcs, ncol(pca$x))
  x <- pca$x[, seq_len(k), drop = FALSE]
  
  km <- stats::kmeans(x, centers = centers, nstart = nstart)
  cl <- km$cluster
  names(cl) <- rownames(sex_eset)
  
  # ---- determine which cluster is male (higher average scaled Y expression) ----
  z <- scale(sex_eset)
  cluster_means <- tapply(rowMeans(z), cl, mean, na.rm = TRUE)
  male_cluster <- as.integer(names(cluster_means)[which.max(cluster_means)])
  
  pred <- ifelse(cl == male_cluster, "M", "F")
  names(pred) <- names(cl)
  
  # ---- accuracy (optional) ----
  acc <- NA_real_
  cm <- NULL
  if (!is.null(true_sex)) {
    truth <- factor(true_sex, levels = c("F", "M"))
    pred_f <- factor(pred[colnames(eset)], levels = c("F", "M"))
    
    cm <- table(truth, pred_f)
    acc <- sum(diag(cm)) / sum(cm)
  }
  
  list(
    cluster = pred,
    acc = acc,
    confusion = cm,
    y_genes = colnames(sex_eset),
    pca = pca,
    kmeans = km
  )
}

