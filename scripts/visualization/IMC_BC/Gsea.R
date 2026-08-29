library(dplyr)
library(tidyr)
library(ggplot2)
library(enrichR)

# --- 1. Inputs: only the gene mapping (for the symbol guard) and the lists ---
protein <- read.csv("ProteinMarker_GeneMapping.csv")

subnetwork_genes <- list(
  c("KRT14", "SNAI2", "KRT5", "EGFR", "TWIST1", "MYC", "CASP3"),
  c("KRT18", "CDH1", "ERBB2", "KRT19", "GATA3", "KRT8", "PGR", "ESR1", "H3F3A", "CA9", "EZH2")
)
# Guard: no symbol can silently fail to match the panel (the CDH11 bug class)
stopifnot(length(setdiff(unlist(subnetwork_genes), protein$Gene.Symbol)) == 0)

enrichR::setEnrichrSite("Enrichr")
dbs <- c("KEGG_2021_Human", "GO_Biological_Process_2025")
save_dir <- "/home/koe/BSNMani_application-main/Github_BSNMani/IMC/Results/GSEA"
if (!dir.exists(save_dir)) dir.create(save_dir, recursive = TRUE, showWarnings = FALSE)

# --- 2. One enrichment query per subnetwork on the FULL gene list ------------
for (i in seq_along(subnetwork_genes)) {
  cat(sprintf("\nProcessing Subnetwork %d (%d genes)...\n",
              i, length(subnetwork_genes[[i]])))

  enr <- enrichR::enrichr(subnetwork_genes[[i]], dbs)

  for (db in dbs) {
    tab <- enr[[db]]
    if (is.null(tab) || nrow(tab) == 0) next
    tab <- tab[tab$Adjusted.P.value < 0.05, ]
    if (nrow(tab) == 0) next
    tab$count <- as.numeric(sub("/.*", "", tab$Overlap))

    db_short <- ifelse(db == "KEGG_2021_Human", "KEGG", "GO")
    out_file <- file.path(save_dir,
                sprintf("Annotation_%s_Subnetwork%d_Pathways.csv", db_short, i))
    write.csv(tab, out_file, row.names = FALSE)
    cat(sprintf("  [Saved] %s: %d significant terms -> %s\n",
                db_short, nrow(tab), basename(out_file)))
  }
}

# Define the base directory for your results
base_dir <- ### set directory
library(dplyr)
library(tidyr)
library(ggplot2)

for (i in 1:2) {
  # 1. Define file paths for both KEGG and GO
  kegg_file <- file.path(base_dir, paste0("Fixed_Filtered_KEGG_Subnetwork", i, "_Pathways.csv"))
  go_file   <- file.path(base_dir, paste0("Fixed_Filtered_GO_Subnetwork",  i, "_Pathways.csv"))

  df_list <- list()

  # 2. Read and label GO data (if it exists)
  if (file.exists(go_file)) {
    gse_go <- read.csv(go_file)
    gse_go$type <- "GO"
    df_list[[length(df_list) + 1]] <- gse_go
  } else {
    warning(paste("GO file not found for subnetwork", i))
  }

  # 3. Read and label KEGG data (if it exists)
  if (file.exists(kegg_file)) {
    gse_kegg <- read.csv(kegg_file)
    gse_kegg$type <- "KEGG"
    df_list[[length(df_list) + 1]] <- gse_kegg
  } else {
    warning(paste("KEGG file not found for subnetwork", i))
  }

  if (length(df_list) == 0) next

  # 4. Combine KEGG and GO datasets
  gse_all <- bind_rows(df_list)

  # 5. Parse the Overlap column to calculate Gene Ratio
  gse_all <- gse_all %>%
    separate(Overlap, into = c("hit_count", "Total"), sep = "/", convert = TRUE, remove = FALSE) %>%
    mutate(gene_ratio = hit_count / Total)

  # 6. Set factors for correct ordering in the plot
  gse_all$category <- factor(gse_all$category, levels = c("activated", "suppressed"))
  gse_all$type     <- factor(gse_all$type,     levels = c("GO", "KEGG"))

  # 6b. GHOST PANELS: guarantee every (type x category) facet cell exists,
  #     even when a direction has no pathways, so both subnetwork figures
  #     share the identical two-column layout. geom_blank() draws nothing
  #     but forces the panel and its x-range into existence.
  needed  <- expand.grid(type     = levels(gse_all$type),
                         category = levels(gse_all$category),
                         stringsAsFactors = FALSE)
  present <- distinct(gse_all, type, category) %>%
    mutate(type = as.character(type), category = as.character(category))
  missing <- anti_join(needed, present, by = c("type", "category"))

  # x-range for the empty column: reuse this figure's own data range so the
  # empty panel's axis looks native. To match Subnetwork 2's suppressed axis
  # EXACTLY instead, replace the next line with:  xr <- c(0.01, 0.08)
  xr <- range(gse_all$gene_ratio, na.rm = TRUE)

  blank_df <- do.call(rbind, lapply(seq_len(nrow(missing)), function(k) {
    ty <- missing$type[k]
    data.frame(type       = factor(ty, levels = levels(gse_all$type)),
               category   = factor(missing$category[k], levels = levels(gse_all$category)),
               Term       = gse_all$Term[gse_all$type == ty][1],  # existing term: adds no y label
               gene_ratio = xr)
  }))
  if (is.null(blank_df)) blank_df <- gse_all[0, c("type", "category", "Term", "gene_ratio")]

  # 7. Generate the plot
  final_plot <- ggplot(gse_all, aes(x = gene_ratio, y = Term)) +
    geom_blank(data = blank_df, aes(x = gene_ratio, y = Term)) +
    geom_point(aes(size = hit_count, color = Adjusted.P.value), alpha = 0.8) +
    # ^ previously color = P.value with legend "p value" while the manuscript
    #   colorbar said p.adjust — now the figure and its legend agree.
    scale_color_gradient(low = "red", high = "blue", name = "p.adjust") +
    scale_size(range = c(3, 10), name = "Gene Count",
               breaks = sort(unique(gse_all$hit_count))) +
    facet_grid(type ~ category, scales = "free") +
    theme_bw() +
    theme(
      strip.text   = element_text(face = "bold", size = 12),
      axis.text.y  = element_text(size = 12),
      axis.text.x  = element_text(size = 10),
      axis.title.y = element_blank(),
      plot.title   = element_text(hjust = 0.5, size = 14, face = "bold"),
      panel.spacing = unit(1.5, "lines")
    ) +
    labs(
      title = "Enrichment: GO and KEGG Terms by Activation Status",
      # ^ the method is over-representation analysis, not GSEA — title corrected.
      x = "Gene Ratio"
    ) +
    guides(
      size  = guide_legend(order = 2),
      color = guide_colorbar(order = 1)
    )

  # 8. Save the combined plot
  out_file <- file.path(base_dir, paste0("Bubble_Plot_Subnetwork_", i, "_GO_KEGG_Combined.png"))
  ggsave(out_file, plot = final_plot, width = 12, height = 8, dpi = 300)

  cat(sprintf("Successfully saved combined plot for Subnetwork %d\n", i))
}

