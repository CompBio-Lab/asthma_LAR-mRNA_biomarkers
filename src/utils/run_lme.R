run_lme = function(data, fdr){
  ## quiet vs. exacerbations
  geoEset_QuietvsExacer <- filter(ungroup(data), phenotype != "FOLLOW UP")
  group <- factor(geoEset_QuietvsExacer$phenotype, levels = c("QUIET", "EXACERBATION"))
  subject <- factor(geoEset_QuietvsExacer$donor)
  sex <- factor(geoEset_QuietvsExacer$sex)
  genes_pval_QuietvsExacer <- list()
  for(i in 4 : ncol(data)){
    gene <- as.vector(as.matrix(geoEset_QuietvsExacer[, i]))
    df <- data.frame(gene, group, subject, sex)
    groupedDat <- groupedData(gene ~ group | subject, data = df)
    fit <- lme(gene ~ group, data = groupedDat, random = ~ 1 | subject, na.action = na.omit)
    genes_pval_QuietvsExacer[[i-3]] <- coef(summary(fit))["groupEXACERBATION",]
  }
  names(genes_pval_QuietvsExacer) <- colnames(data)[-c(1:3)]
  genes_adjpval_QuietvsExacer <- as.data.frame(do.call(rbind, genes_pval_QuietvsExacer))
  genes_adjpval_QuietvsExacer$adj.P.Val <- p.adjust(genes_adjpval_QuietvsExacer$`p-value`, "BH")
  genes_adjpval_QuietvsExacer <- genes_adjpval_QuietvsExacer[order(genes_adjpval_QuietvsExacer$`p-value`), ]
  genes_adjpval_QuietvsExacer$feature <- rownames(genes_adjpval_QuietvsExacer)
  genes_adjpval_QuietvsExacer$GSE <- "GSE19301"
  genes_adjpval_QuietvsExacer$contrast <- "Exacer - Quiet"
  genes_adjpval_QuietvsExacer$tissue <- "blood-PBMCs"
  genes_adjpval_QuietvsExacer$fc <- ifelse(genes_adjpval_QuietvsExacer$Value > 0, "UP", "DOWN")
  
  # exacerbation vs. followup
  geoEset_ExacerVsFollowUp <- filter(ungroup(data), phenotype != "QUIET")
  group <- factor(geoEset_ExacerVsFollowUp$phenotype, levels = c("FOLLOW UP","EXACERBATION"))
  subject <- factor(geoEset_ExacerVsFollowUp$donor)
  sex <- factor(geoEset_ExacerVsFollowUp$sex)
  genes_pval_ExacerVsFollowUp <- list()
  for(i in 4 : ncol(data)){
    gene <- as.vector(as.matrix(geoEset_ExacerVsFollowUp[, i]))
    df <- data.frame(gene, group, subject, sex)
    groupedDat <- groupedData(gene ~ group | subject, data = df)
    fit <- lme(gene ~ group, data = groupedDat, random = ~ 1 | subject, na.action = na.omit)
    genes_pval_ExacerVsFollowUp[[i-2]] <- coef(summary(fit))["groupEXACERBATION",]
  }
  names(genes_pval_ExacerVsFollowUp) <- colnames(data)[-c(1:3)]
  genes_adjpval_ExacerVsFollowUp <- as.data.frame(do.call(rbind, genes_pval_ExacerVsFollowUp))
  genes_adjpval_ExacerVsFollowUp$adj.P.Val <- p.adjust(genes_adjpval_ExacerVsFollowUp$`p-value`, "BH")
  genes_adjpval_ExacerVsFollowUp <- genes_adjpval_ExacerVsFollowUp[order(genes_adjpval_ExacerVsFollowUp$`p-value`), ]
  genes_genSym_ExacerVsFollowUp <- rownames(genes_adjpval_ExacerVsFollowUp)[genes_adjpval_ExacerVsFollowUp$adj.P.Val < fdr]
  length(genes_genSym_ExacerVsFollowUp)
  genes_adjpval_ExacerVsFollowUp$feature <- rownames(genes_adjpval_ExacerVsFollowUp)
  genes_adjpval_ExacerVsFollowUp$GSE <- "GSE19301"
  genes_adjpval_ExacerVsFollowUp$contrast <- "Exacer - Followup"
  genes_adjpval_ExacerVsFollowUp$tissue <- "blood-PBMCs"
  genes_adjpval_ExacerVsFollowUp$fc <- ifelse(genes_adjpval_QuietvsExacer$Value > 0, "UP", "DOWN")
  
  rbind(genes_adjpval_QuietvsExacer, genes_adjpval_ExacerVsFollowUp)
}