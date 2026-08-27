#########################
## draw cell type proportion plot
#########################
rm(list = ls())
setwd("D:\\Research_project\\For Test\\Cell_type_specific_BSNMani\\data")
load("data_for_uni_27.RData")
colnames(cell_loc_temp[["H21.33.006"]])[3] = "Subclass"

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(scales)
})

## ---- Merge Subclass data from all patients ----
all_cells <- dplyr::bind_rows(
  lapply(cell_loc_temp, function(df) {
    df <- as.data.frame(df)  
    df["Subclass"]            
  })
)

## ---- calculate proportion ----
prop_df <- all_cells %>%
  filter(!is.na(Subclass) & Subclass != "") %>%
  count(Subclass, name = "n") %>%
  mutate(Proportion = n / sum(n)) %>%
  arrange(desc(Proportion))

print(head(prop_df))

## ---- draw proportion plot ----
p <- ggplot(prop_df, aes(x = reorder(Subclass, -Proportion), y = Proportion)) +
  geom_col(fill = "#74c0fc") +
  labs(title = "Cell Type Proportion Across All Patients",
       x = NULL, y = "Proportion") +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  theme_minimal(base_size = 18) +
  theme(
    panel.grid       = element_blank(),
    axis.text.x      = element_text(angle = 45, hjust = 1, size = 16),
    axis.text.y      = element_text(size = 16),
    axis.title.y     = element_text(size = 18),
    panel.background = element_rect(fill = "white", colour = NA),  
    plot.background  = element_rect(fill = "white", colour = NA)   
  )
p

setwd("D:\\Research_project\\For Test\\Cell_type_specific_BSNMani")
ggsave("celltype_proportion_all_patients_uni.png", p,
       width = 12, height = 7, dpi = 600, bg = "white")
