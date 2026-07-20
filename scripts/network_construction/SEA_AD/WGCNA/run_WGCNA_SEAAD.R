# this R script is used for generating co expression matrx for 27 patients using wgcna

rm(list = ls())
load("D:\\Research_project\\For Test\\H20.33.001\\coexpression-generation\\data_uni\\data_for_uni_27.RData")

############################################################
# 0. Load packages
############################################################

library(Giotto)
#install.packages("BiocManager")
#BiocManager::install("impute")
BiocManager::install("preprocessCore")
library(WGCNA)

options(stringsAsFactors = FALSE)
allowWGCNAThreads()



############################################################
# Fisher Z
############################################################

FisherZ <- function(cor_mat, eps = 1e-6) {
  cor_mat[cor_mat >= 1] <- 1 - eps
  cor_mat[cor_mat <= -1] <- -1 + eps
  
  z_mat <- 0.5 * log((1 + cor_mat) / (1 - cor_mat))
  return(z_mat)
}

############################################################
# Extract normalized expression matrix from Giotto object
############################################################

extract_giotto_normalized_matrix <- function(gobject) {
  
  if (exists("get_expression_values", mode = "function")) {
    norm_expr <- get_expression_values(
      gobject,
      values = "normalized",
      spat_unit = "cell",
      feat_type = "rna",
      output = "matrix"
    )
  } else if (exists("getExpression", mode = "function")) {
    norm_expr <- getExpression(
      gobject,
      values = "normalized",
      spat_unit = "cell",
      feat_type = "rna",
      output = "matrix"
    )
  } else {
    stop("Cannot find Giotto expression extraction function.")
  }
  
  norm_expr <- as.matrix(norm_expr)
  return(norm_expr)
}

############################################################
# One patient pipeline
# Input: cell × gene raw count matrix
############################################################

make_giotto_wgcna_fisherz_network <- function(count_mat,
                                              genes_use = NULL,
                                              scalefactor = 10000) {
  
  # count_mat is cell × gene
  count_mat <- as.matrix(count_mat)
  
  # Optional: restrict genes
  if (!is.null(genes_use)) {
    genes_use <- intersect(genes_use, colnames(count_mat))
    count_mat <- count_mat[, genes_use, drop = FALSE]
  }
  
  # Remove empty cells
  count_mat <- count_mat[rowSums(count_mat) > 0, , drop = FALSE]
  
  # Remove all-zero genes
  count_mat <- count_mat[, colSums(count_mat) > 0, drop = FALSE]
  
  if (nrow(count_mat) < 3) {
    stop("After filtering, fewer than 3 cells remain.")
  }
  
  if (ncol(count_mat) < 2) {
    stop("After filtering, fewer than 2 genes remain.")
  }
  
  ##########################################################
  # Giotto expects gene × cell expression matrix
  # Therefore transpose before creating Giotto object
  ##########################################################
  
  gobject <- createGiottoObject(
    raw_exprs = t(count_mat)
  )
  
  gobject <- normalizeGiotto(
    gobject,
    scalefactor = scalefactor,
    log_norm = TRUE,
    verbose = TRUE
  )
  
  ##########################################################
  # Extract normalized expression
  # Usually returned as gene × cell
  ##########################################################
  
  norm_expr <- extract_giotto_normalized_matrix(gobject)
  
  ##########################################################
  # Convert normalized expression to cell × gene for WGCNA
  ##########################################################
  
  if (all(colnames(count_mat) %in% rownames(norm_expr))) {
    datExpr <- t(norm_expr[colnames(count_mat), , drop = FALSE])
  } else if (all(colnames(count_mat) %in% colnames(norm_expr))) {
    datExpr <- norm_expr[, colnames(count_mat), drop = FALSE]
  } else {
    stop("Cannot match gene names between count_mat and Giotto normalized matrix.")
  }
  
  ##########################################################
  # WGCNA co-expression calculation
  ##########################################################
  
  raw_corr <- adjacency(
    datExpr = datExpr,
    type = "signed",
    power = 1
  )
  
  raw_corr <- 2 * raw_corr - 1
  raw_corr[is.na(raw_corr)] <- 0
  diag(raw_corr) <- 0
  
  transformed_corr <- FisherZ(raw_corr)
  diag(transformed_corr) <- 0
  
  return(list(
    normalized_expression_cell_by_gene = datExpr,
    raw_correlation = raw_corr,
    fisherz_correlation = transformed_corr
  ))
}


genes_use <- colnames(express_matrix_temp[[1]])

result_list <- lapply(express_matrix_temp, function(count_mat) {
  make_giotto_wgcna_fisherz_network(
    count_mat = count_mat,
    genes_use = genes_use,
    scalefactor = 10000
  )
})

wgcna_raw_ls <- lapply(result_list, function(x) x$raw_correlation)

wgcna_transformed_ls <- lapply(result_list, function(x) x$fisherz_correlation)

saveRDS(wgcna_transformed_ls,file = "baseline_co_expression.RDS")
