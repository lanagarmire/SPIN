############################################################
# Generate hdWGCNA TOM co-expression matrices for 27 patients
# Oligodendrocyte-only data
############################################################

rm(list = ls())

############################################################
# 0. Load data
############################################################

setwd("D:/Research_project/For Test/H20.33.001/Oli/data")

cell_loc_temp <- readRDS("cell_loc_uni_Oli.RDS")
express_matrix_temp <- readRDS("express_uni_Oli.RDS")

############################################################
# 1. Load packages
############################################################

library(Seurat)
library(hdWGCNA)
library(WGCNA)

set.seed(12345)
enableWGCNAThreads(nThreads = 8)

############################################################
# 2. Define paths
############################################################

tom_dir <- "D:/Research_project/For Test/H20.33.001/Oli/data/TOM"

dir.create(tom_dir, recursive = TRUE, showWarnings = FALSE)

gene_names <- readRDS(
  "D:/Research_project/For Test/H20.33.001/coexpression-generation/hdWGCNA/common_genes.RDS"
)

genes_use <- colnames(express_matrix_temp[[1]])

# Check consistency
stopifnot(setequal(gene_names, genes_use))

# Use the common gene order as final 140-gene order
all_genes <- gene_names

############################################################
# 3. Patient IDs
############################################################

patient_ids <- c(
  "H20.33.035","H21.33.038","H20.33.004","H21.33.023","H21.33.001",
  "H21.33.025","H20.33.001","H21.33.032","H21.33.031","H21.33.019",
  "H21.33.016","H20.33.040","H20.33.012","H20.33.025","H21.33.015",
  "H21.33.005","H21.33.040","H21.33.028","H21.33.012","H20.33.015",
  "H20.33.044","H21.33.022","H21.33.011","H21.33.013","H21.33.014",
  "H21.33.021","H21.33.006"
)

############################################################
# 4. Function: expand TOM to 140 * 140
############################################################

expand_tom_to_140 <- function(tom_mat, valid_genes, all_genes) {
  
  full_mat <- matrix(
    0,
    nrow = length(all_genes),
    ncol = length(all_genes),
    dimnames = list(all_genes, all_genes)
  )
  
  common_genes <- intersect(valid_genes, all_genes)
  
  tom_mat <- tom_mat[common_genes, common_genes, drop = FALSE]
  
  full_mat[common_genes, common_genes] <- tom_mat
  
  diag(full_mat) <- 0
  
  return(full_mat)
}

############################################################
# 5. Function: run hdWGCNA for one patient
############################################################

run_hdwgcna_one_patient <- function(count_mat,
                                    loc_df,
                                    patient_id,
                                    genes_use,
                                    wgcna_name = "hdwgcna_spatial",
                                    min_cells = 10,
                                    k = 5,
                                    max_shared = 2,
                                    soft_power = 6) {
  
  count_mat <- as.matrix(count_mat)
  loc_df <- as.data.frame(loc_df)
  
  stopifnot(nrow(count_mat) == nrow(loc_df))
  
  cell_ids <- paste0(patient_id, "_cell_", seq_len(nrow(count_mat)))
  rownames(count_mat) <- cell_ids
  rownames(loc_df) <- cell_ids
  
  seu <- CreateSeuratObject(
    counts = t(count_mat),
    assay = "RNA",
    project = patient_id
  )
  
  seu@meta.data$patient_id <- patient_id
  seu@meta.data$cell_group <- "all_cell"
  
  seu@meta.data$row <- as.numeric(loc_df$sdimy)
  seu@meta.data$col <- as.numeric(loc_df$sdimx)
  seu@meta.data$imagerow <- as.numeric(loc_df$sdimy)
  seu@meta.data$imagecol <- as.numeric(loc_df$sdimx)
  
  seu <- NormalizeData(seu)
  seu <- ScaleData(seu, features = genes_use)
  seu <- RunPCA(seu, features = genes_use)
  
  seu <- SetupForWGCNA(
    seu,
    wgcna_name = wgcna_name,
    features = genes_use
  )
  
  seu <- MetacellsByGroups(
    seurat_obj = seu,
    group.by = "cell_group",
    reduction = "pca",
    k = k,
    max_shared = max_shared,
    ident.group = "cell_group",
    assay = "RNA",
    min_cells = min_cells
  )
  
  seu <- NormalizeMetacells(seu)
  
  seu <- SetDatExpr(
    seu,
    assay = "RNA",
    wgcna_name = wgcna_name
  )
  
  # fixed soft power, avoid Inf soft_power problem
  seu <- ConstructNetwork(
    seu,
    soft_power = soft_power,
    tom_name = patient_id,
    overwrite_tom = TRUE,
    networkType = "signed",
    TOMType = "signed",
    wgcna_name = wgcna_name
  )
  
  return(seu)
}

############################################################
# 6. Run hdWGCNA for all patients
############################################################

failed_patients <- list()

for (pid in patient_ids) {
  
  message("Running hdWGCNA for patient: ", pid)
  
  tryCatch({
    
    seu <- run_hdwgcna_one_patient(
      count_mat = express_matrix_temp[[pid]],
      loc_df = cell_loc_temp[[pid]],
      patient_id = pid,
      genes_use = genes_use,
      wgcna_name = "hdwgcna_spatial",
      min_cells = 10,
      k = 5,
      max_shared = 2,
      soft_power = 6
    )
    
    valid_genes <- colnames(seu@misc$hdwgcna_spatial$datExpr)
    
    saveRDS(
      valid_genes,
      file = file.path(tom_dir, paste0(pid, "_valid_genes.rds"))
    )
    
    rm(seu)
    gc()
    
  }, error = function(e) {
    
    message("Failed patient: ", pid)
    message("Error: ", e$message)
    failed_patients[[pid]] <<- e$message
    
  })
}

saveRDS(
  failed_patients,
  file = "D:/Research_project/For Test/H20.33.001/Oli/data/hdWGCNA_failed_patients_Oli.rds"
)

############################################################
# 7. Read TOM files and build 140 * 140 matrix list
############################################################

tom_files <- list.files(
  path = tom_dir,
  pattern = "_TOM\\.rda$",
  full.names = TRUE
)

hdwgcna_tom_list <- list()

for (tom_file in tom_files) {
  
  pid <- basename(tom_file)
  pid <- sub("_TOM\\.rda$", "", pid)
  
  message("Loading TOM matrix for patient: ", pid)
  
  load(tom_file)  # loads consTomDS
  
  tom_mat <- as.matrix(consTomDS)
  
  valid_gene_file <- file.path(tom_dir, paste0(pid, "_valid_genes.rds"))
  
  if (file.exists(valid_gene_file)) {
    valid_genes <- readRDS(valid_gene_file)
  } else {
    valid_genes <- all_genes[seq_len(nrow(tom_mat))]
    warning("No valid_genes file found for ", pid, 
            ". Using first n genes as fallback. Please verify.")
  }
  
  if (nrow(tom_mat) == length(valid_genes)) {
    rownames(tom_mat) <- valid_genes
    colnames(tom_mat) <- valid_genes
  } else {
    stop(
      "Dimension mismatch for ", pid,
      ": TOM dim = ", nrow(tom_mat),
      ", valid genes = ", length(valid_genes)
    )
  }
  
  tom_mat_140 <- expand_tom_to_140(
    tom_mat = tom_mat,
    valid_genes = valid_genes,
    all_genes = all_genes
  )
  
  hdwgcna_tom_list[[pid]] <- tom_mat_140
  
  rm(consTomDS, tom_mat, tom_mat_140)
  gc()
}

############################################################
# 8. Reorder and save
############################################################

available_patients <- intersect(patient_ids, names(hdwgcna_tom_list))
missing_patients <- setdiff(patient_ids, names(hdwgcna_tom_list))

message("Available patients: ", length(available_patients))
message("Missing patients: ", paste(missing_patients, collapse = ", "))

hdwgcna_tom_list <- hdwgcna_tom_list[available_patients]

# Check dimensions
print(sapply(hdwgcna_tom_list, dim))

saveRDS(
  hdwgcna_tom_list,
  file = "D:/Research_project/For Test/H20.33.001/coexpression-generation/hdWGCNA/hdWGCNA_co_express_list_uni_Oli.RDS"
)


test = as.data.frame(hdwgcna_tom_list[["H21.33.006"]])
