mtx <- readRDS("/home/yangiwen/n/workspace/bsnmani/data/bc/Basel_patient_metal_count.RDS")
meta <- read.csv("/home/yangiwen/n/workspace/bsnmani/data/bc/Basel_CleanCore_Age_OS.csv")
protein <- read.csv("/home/yangiwen/n/workspace/bsnmani/data/bc/ProteinMarker_GeneMapping.csv")

# Subset the 253 patients
mtx <- mtx[meta$core]

# Create concat seurat object count matrix
obj_cnt <- NULL
for (patient in names(mtx)) {
  count <- as.data.frame(mtx[[patient]]$count)
  rownames(count) <- count$metal_tag
  count <- count[, -1]
  count <- count[protein$Metal,]
  obj_cnt <- if (is.null(obj_cnt)) count else cbind(obj_cnt, count)
}

# Mark case and control
median_survival_time <- median(meta$OSmonth)
meta[meta$OSmonth >= median_survival_time, "group"] <- "control"
meta[meta$OSmonth < median_survival_time, "group"] <- "case"

# Create seurat object meta
obj_meta <- data.frame(row.names = colnames(obj_cnt))
obj_meta$id <- rownames(obj_meta)
obj_meta$core <- sub("(_[^_]*)$", "", rownames(obj_meta))
obj_meta <- merge(obj_meta, meta, on = "core")
rownames(obj_meta) <- obj_meta$id

# limma
model_matrix <- model.matrix(~0 + group, data = obj_meta)
fit <- limma::lmFit(obj_cnt, model_matrix)
contrasts <- limma::makeContrasts(groupcase - groupcontrol, levels = colnames(coef(fit)))
fit <- limma::contrasts.fit(fit, contrasts = contrasts)
fit <- limma::eBayes(fit)
de_out <- limma::topTable(fit, sort.by = "P", number = nrow(fit))
degs <- data.frame(
  row.names = rownames(de_out),
  score = de_out$t,
  pval = de_out$P.Value,
  adj_pval = de_out$adj.P.Val
)

# Replace metal tag with gene name
rownames(degs) <- protein[match(rownames(de_out), protein$Metal), "Gene.Symbol"]

# EnrichR
degs <- degs[order(-abs(degs$score)),]
subnetwork_genes <- list(
  c("GATA3", "MKI67", "ESR1", "MS4A1", "PGR", "EZH2", "CA9", "SNAI2", "H3F3A", "TWIST1", "MYC", "KRT5", "KRT14", "EGFR", "CASP3", "CDH1", "ERBB2", "KRT8", "KRT19", "KRT18", "KRT7"),
  c("CDH1", "ERBB2", "KRT8", "KRT19", "KRT18", "KRT7", "FN1", "ACTA2", "VIM"),
  c("PTPRC", "CD3E", "CD44", "CD68", "VWF", "FN1", "ACTA2", "VIM"),
  c("TWIST1", "MYC", "KRT5", "KRT14", "EGFR", "CASP3", "CDH1", "ERBB2", "KRT8", "KRT19", "KRT18", "KRT7", "CD44", "CD68", "VWF", "FN1", "ACTA2", "VIM")
)
library(enrichR)
enrichR::setEnrichrSite("Enrichr")
for (i in seq_along(subnetwork_genes)) {
  pos_genes <- intersect(subnetwork_genes[[i]], rownames(degs[degs$score > 0,]))
  neg_genes <- intersect(subnetwork_genes[[i]], rownames(degs[degs$score < 0,]))
  pos_gse <- enrichR::enrichr(pos_genes, c(
    "KEGG_2021_Human"
  ))
  pos_gse <- do.call(rbind, lapply(names(pos_gse), function(n) cbind(pos_gse[[n]], db = n)))
  pos_gse <- pos_gse[pos_gse$P.value < 0.05,]
  # pos_gse <- pos_gse[grep("signaling pathway|Apoptosis|Cell cycle", pos_gse$Term, ignore.case = T),]
  pos_gse$category <- "activated"
  neg_gse <- enrichR::enrichr(neg_genes, c(
    "KEGG_2021_Human"
  ))
  neg_gse <- do.call(rbind, lapply(names(neg_gse), function(n) cbind(neg_gse[[n]], db = n)))
  neg_gse <- neg_gse[neg_gse$P.value < 0.05,]
  # neg_gse <- neg_gse[grep("signaling pathway|Apoptosis|Cell cycle", neg_gse$Term, ignore.case = T),]
  neg_gse$category <- "suppressed"
  gse <- rbind(pos_gse, neg_gse)
  gse <- na.omit(gse)
  write.csv(gse, file.path("/home/yangiwen/n/workspace/bsnmani/output", "gsea", paste0("pathway_", i, ".csv")))
}

for (i in seq_along(subnetwork_genes)) {
  gse <- read.csv(file.path("/home/yangiwen/n/workspace/bsnmani/output", "gsea", paste0("pathway_", i, "_filtered.csv")), row.names = 1)
  gse$hit_count <- as.integer(sub("/.*", "", gse$Overlap))
  gse$gene_ratio <- gse$hit_count / nrow(protein)
  gse$category <- factor(gse$category, levels = c("Suppressed", "Activated"))
  png(file.path("/home/yangiwen/n/workspace/bsnmani/output", "gsea", paste0("subnetwork_", i, ".png")), width = 1792, height = 1024, res = 120)
  print(ggplot2::ggplot(gse, ggplot2::aes(x = gene_ratio, y = Term, size = hit_count, color = P.value)) +
          ggplot2::geom_point(alpha = 0.8) +
          ggplot2::scale_color_gradient(low = "red", high = "blue", name = "p value") +
          ggplot2::scale_size(range = c(3, 10), name = "Gene Count", breaks = sort(unique(gse$hit_count))) +
          ggplot2::facet_wrap(~category, nrow = 1, scales = "free_x") +
          ggplot2::theme_bw() +
          ggplot2::theme(
            strip.text = ggplot2::element_text(face = "bold", size = 12),
            axis.text.y = ggplot2::element_text(size = 18),
            axis.title.y = ggplot2::element_blank(),
            plot.title = ggplot2::element_text(hjust = 0.5)
          ) +
          ggplot2::labs(
            title = "GSEA Enrichment: KEGG Terms by Activation Status",
            x = "Gene Ratio"
          ) +
          ggplot2::guides(
            size = ggplot2::guide_legend(order = 2),
            color = ggplot2::guide_colorbar(order = 1)
          )
  )
  dev.off()
}

