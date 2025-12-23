convert_enrichr_to_links <- function(enrichr_result, fdr_cutoff){
  result <- lapply(names(enrichr_result), function(db){
    res <- subset(enrichr_result[[db]], Adjusted.P.value < fdr_cutoff)
    res_genes <- lapply(res$Genes, function(i) strsplit(i, ";")[[1]])
    names(res_genes) <- res$Term
    if(db == "Rummagene_transcription_factors"){
      df <- data.frame(geneset = rep(res$tf, sapply(res_genes, length)),
                       genes = unlist(res_genes))
      
    } else {
      df <- data.frame(geneset = rep(names(res_genes), sapply(res_genes, length)),
                       genes = unlist(res_genes))
    }
    df$db <- db
    df %>% 
      group_by(geneset, genes, db) %>% 
      dplyr::slice(1)
  }) %>% 
    do.call(rbind, .)
}