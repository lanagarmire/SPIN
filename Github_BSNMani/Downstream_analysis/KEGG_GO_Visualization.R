
#############################
## Downstream bioinformatics analysis
#############################

####
# packages preparation
####
rm(list = ls())
library(enrichR)
library(ggplot2)
library(dplyr)
library(tidyr)

####
# enrich analysis Pipeline 
####

## input genes name vector          ****
genes_up = c("CD44","FN1","CASP3")
genes_down = c("CD45","CD68","ACTA2","VIM")

dbs = c("KEGG_2021_Human","GO_Biological_Process_2025")
enrich_up = enrichr(genes_up, dbs)
enrich_down = enrichr(genes_down, dbs)

##  KEGG and res
enrich_kegg_up = enrich_up$KEGG_2021_Human
enrich_kegg_down = enrich_down$KEGG_2021_Human
enrich_go_up = enrich_up$GO_Biological_Process_2025
enrich_go_down = enrich_down$GO_Biological_Process_2025

enrich_kegg_up_filter = enrich_kegg_up[enrich_kegg_up$Adjusted.P.value < 0.05, ]
enrich_kegg_down_filter = enrich_kegg_down[enrich_kegg_down$Adjusted.P.value < 0.05, ]

enrich_go_up_filter = enrich_go_up[enrich_go_up$Adjusted.P.value < 0.05, ]
enrich_go_down_filter = enrich_go_down[enrich_go_down$Adjusted.P.value < 0.05, ]

enrich_kegg_up_filter$count = as.numeric(sub("/.*", "", enrich_kegg_up_filter$Overlap))
enrich_kegg_down_filter$count = as.numeric(sub("/.*", "", enrich_kegg_down_filter$Overlap))
enrich_go_up_filter$count = as.numeric(sub("/.*", "", enrich_go_up_filter$Overlap))
enrich_go_down_filter$count = as.numeric(sub("/.*", "", enrich_go_down_filter$Overlap))

## final df for heatmap
## manually select best terms          ****
go_up_final = enrich_go_up_filter[c(4,3,9),]
go_down_final = enrich_go_down_filter[c(16,2),]
kegg_up_final = enrich_kegg_up_filter[c(1,9,13,7),]
kegg_down_final = enrich_kegg_down_filter[c(1),]



###############
## heatmap pipeline
###############


go_up_final$type = "GO"
go_down_final$type = "GO"
kegg_up_final$type = "KEGG"
kegg_down_final$type = "KEGG"

go_up_final$status = "activated"
go_down_final$status = "suppressed"
kegg_up_final$status = "activated"
kegg_down_final$status = "suppressed"

all_df = bind_rows(go_up_final, go_down_final, kegg_up_final, kegg_down_final)

all_df = all_df %>%
  rename(Term = Term, GeneRatio = Overlap, padj = Adjusted.P.value)

all_df <- all_df %>%
  separate(GeneRatio, into = c("Count", "Total"), sep = "/", convert = TRUE) %>%
  mutate(GeneRatio_numeric = Count / Total)

setwd("D:\\Research_project\\For Test\\helphaowen\\helphaowen_later\\Downstream_res\\new_results")
saveRDS(all_df, "subnetwork2_res.RDS")

rm(list = ls())
setwd("D:\\Research_project\\For Test\\helphaowen\\helphaowen_later\\Downstream_res\\new_results")
all_df = readRDS("subnetwork2_res.RDS")   ****

final_plot = ggplot(all_df, aes(x = GeneRatio_numeric, y = Term)) +
  geom_point(aes(size = count, color = padj)) +
  scale_color_gradient(low = "red", high = "blue", name = "p.adjust") +
  scale_size_continuous(name = "Gene Count") +
  facet_grid(type ~ status, scales = "free") + 
  theme_bw() +
  theme(
    strip.text = element_text(size = 12),
    axis.text.y = element_text(size = 10),
    axis.text.x = element_text(size = 10),
    panel.spacing = unit(1.5, "lines"),
    plot.title = element_text(size = 14, face = "bold")
  ) +
  labs(title = "GSEA Enrichment: GO and KEGG Terms by Activation Status",
       x = "Gene Ratio", y = NULL)
setwd("D:\\Research_project\\For Test\\helphaowen\\helphaowen_later\\Downstream_res\\new_results")

## manually edit file name           ****
ggsave("subnetwork2.png", plot = final_plot, width = 10, height = 6, dpi = 600)




##########################
## draw tripartite plot
##########################
rm(list = ls())
library(dplyr)
setwd("D:\\Research_project\\For Test\\helphaowen\\helphaowen_later\\Downstream_res\\new_results")

rds_files = list.files(pattern = "\\.RDS$")
rds_list = lapply(rds_files, readRDS)
names(rds_list) = tools::file_path_sans_ext(rds_files)

for (i in c(1,2))    #**** i in c(1,2,3,,,, number of subnetworks)
{
  rds_list[[i]]$'subnetwork' = paste0("subnetwork",i)
}
df_all = rbind(rds_list[[1]],rds_list[[2]])
df_all_final = df_all[,c("Term","status","Genes","subnetwork","padj")]
library(writexl)
write_xlsx(df_all_final, "df_all.xlsx")
df_unique = df_all_final %>% distinct(Term, .keep_all = TRUE)

library(writexl)
write_xlsx(df_unique, "df_unique.xlsx")

df_unique = df_unique[order(df_unique$padj), ][1:30, ]
write_xlsx(df_unique, "df_tripartite.xlsx")

library(dplyr)
library(tidyr)

df_long_gene = df_unique %>%
  separate_rows(Genes, sep = ";")

edge_gene_term = df_long_gene %>%
  transmute(
    source = Genes,
    target = Term,
    interaction = "Gene-Term",
    color = "gray",
    width = 1
  )
edge_term_subnet = df_unique %>%
  transmute(
    source = Term,
    target = subnetwork,
    interaction = "Term-Subnet",
    color = ifelse(status == "activated", "red", "blue"),
    width = -log10(padj)
  )
edge_table = bind_rows(edge_gene_term, edge_term_subnet)
write.csv(edge_table, "edges.csv", row.names = FALSE, quote = FALSE)
# Gene
gene_nodes = df_long_gene %>%
  distinct(Genes) %>%
  transmute(id = Genes, type = "Gene")

# Term
term_nodes = df_unique %>%
  distinct(Term) %>%
  transmute(id = Term, type = "Term")

# Subnetworks
subnet_nodes = df_unique %>%
  distinct(subnetwork) %>%
  transmute(id = subnetwork, type = "Subnetwork")

node_table = bind_rows(gene_nodes, term_nodes, subnet_nodes)
write.csv(node_table, "nodes.csv", row.names = FALSE, quote = FALSE)


