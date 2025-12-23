#' Plot pathway-level feature importance across comparisons and sample types
#'
#' This function visualizes ranked feature importance scores aggregated at the
#' pathway level. For each comparison and sample type, the top-ranked features
#' are selected, aggregated across comparisons, filtered by a cumulative rank
#' threshold, and displayed as faceted bar plots.
#'
#' The function is designed for use with feature-importance outputs from
#' machine-learning models such as Random Forests or PLS-DA, where importance
#' scores have already been summarized into an annotation table.
#'
#' @param annot A data frame containing feature importance annotations with
#'   at least the following columns:
#'   \describe{
#'     \item{feature}{Gene or feature name}
#'     \item{Pathway}{Pathway or functional annotation}
#'     \item{avg_rank}{Average (or relative) importance rank}
#'     \item{comparison}{Model comparison label (e.g., case vs control)}
#'     \item{sampletype}{Sample type (e.g., PBMC, BALF)}
#'   }
#'
#' @param top_n Integer specifying the number of top-ranked features to retain
#'   per comparison and sample type (default: 30).
#'
#' @param sum_threshold Numeric threshold applied to the summed importance ranks
#'   across comparisons. Features with a cumulative rank less than or equal to
#'   this value are excluded (default: 2).
#'
#' @param include_unassigned Logical indicating whether to include features
#'   annotated as \code{"Unassigned"} pathways. If \code{TRUE}, only unassigned
#'   pathways are shown; if \code{FALSE}, unassigned pathways are excluded.
#' @export
plot_varimp_pathways <- function(
    annot,
    top_n = 30,
    sum_threshold = 2,
    include_unassigned = TRUE
) {
  
  ## ---- pathway filter ----
  pathway_filter <- if (include_unassigned) {
    rlang::expr(Pathway == "Unassigned")
  } else {
    rlang::expr(Pathway != "Unassigned")
  }
  
  ## ---- select features to keep ----
  keep_feat <- annot %>% 
    filter(!!pathway_filter) %>%
    group_by(comparison, sampletype) %>% 
    arrange(desc(avg_rank)) %>% 
    dplyr::slice(1:top_n) %>% 
    group_by(Pathway, feature) %>% 
    summarise(s = sum(avg_rank), .groups = "drop") %>% 
    filter(s > sum_threshold) %>% 
    mutate(pathway_feature = paste(Pathway, feature, sep = "_"))
  
  ## ---- plotting data ----
  p <- annot %>% 
    mutate(pathway_feature = paste(Pathway, feature, sep = "_")) %>%
    filter(pathway_feature %in% keep_feat$pathway_feature) %>%
    mutate(feature = paste0("*", feature, "*")) %>% 
    ggplot(aes(
      x = avg_rank,
      y = tidytext::reorder_within(feature, avg_rank, Pathway),
      fill = comparison
    )) +
    geom_bar(stat = "identity") +
    facet_grid(Pathway ~ sampletype, scales = "free") +
    tidytext::scale_y_reordered() +
    theme_bw() +
    theme(
      axis.title.y = element_blank(),
      axis.text.y = ggtext::element_markdown(),
      legend.position = "right",
      strip.text.y.right = element_text(angle = 0)
    ) +
    xlab("Average Rank")
  
  return(p)
}