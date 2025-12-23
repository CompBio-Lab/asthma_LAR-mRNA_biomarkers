#' Calculate desriptive statistics for each group and compare with lm() or lme()
#'
#' takes in data and determines the summary statistics (Mean & SD) and also compares the levels of the groups (binary) variable
#' @param demo is an (nxp) dataset
#' @param groups specifies the column name of a binary variable in demo
#' @param variables is a vector of column names to be compared between the 2 groups
#' @param paired boolean (T/F) - default = FALSE (repeated measures or not)
#' @param pairing is a column name that containing the pairing information
#' @export
descriptiveStat = function(demo, groups, variables, paired = FALSE, pairing = NULL){
  library(dplyr)
  library(tidyr)
  library(broom)
  
  if(all(paired)){
    X <- demo[, c(variables, groups, pairing), drop = FALSE]
    colnames(X) <- c(variables, "Group", "Pairing")
    lvls <- levels(X$Group)
    meanSD <- X %>% gather(Variable, Value, -c(Group, Pairing)) %>% dplyr::group_by(Variable,
                                                                                    Group) %>% dplyr::summarise(MEAN = mean(Value, na.rm = TRUE),
                                                                                                                SD = sd(Value, na.rm = TRUE))
    
    pval0 <- X %>% gather(Variable, Value, -c(Group, Pairing)) %>% dplyr::group_by(Variable) %>%
      nest() %>% dplyr::mutate(model = purrr::map(data,
                                                  ~tryCatch(lme(Value ~ Group, random = ~ 1 | Pairing, data = .), error = function(e) NA)
      ))
    pval <- do.call(rbind, lapply(pval0$model, function(i){
      tryCatch(summary(i)$tTable[2,], error = function(e) NA)
    })) %>%
      data.frame %>% mutate(Variable = variables, term = paste("Group", lvls[2]),
                            BH.FDR = p.adjust(p.value, "BH"))
  } else {
    X <- demo[, c(variables, groups), drop = FALSE]
    colnames(X) <- c(variables, "Group")
    lvls <- levels(X$Group)
    meanSD <- X %>% 
      gather(Variable, Value, -Group) %>% 
      dplyr::group_by(Variable, Group) %>%
      dplyr::summarise(MEAN = mean(Value, na.rm = TRUE),
                       SD = sd(Value, na.rm = TRUE))
    
    pval <- X %>% 
      tidyr::gather(Variable, Value, -Group) %>% 
      dplyr::group_by(Variable) %>%
      tidyr::nest() %>% 
      dplyr::mutate(model = purrr::map(data, ~ .x %>% 
                                         summarise(model = list(broom::tidy(lm(Value ~ Group)))))) %>%
      dplyr::select(-data) %>% 
      tidyr::unnest(model) %>%
      tidyr::unnest(model) %>%
      group_by(Variable) %>% dplyr::slice(2)
    pval$BH.FDR <- p.adjust(pval$p.value, "BH")
  }
  
  return(list(meanSD = meanSD, pval = pval))
}