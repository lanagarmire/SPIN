# generate hdWGCNA co-expression matrix for 27 uni patients

setwd("D:\\Research_project\\For Test\\H20.33.001\\coexpression-generation\\data_uni")
load("data_for_uni_27.RData")

###########################################
# essential function
library(Seurat)

# hdWGCNA
BiocManager::install("GeneOverlap")
BiocManager::install("UCell")
remotes::install_github("smorabit/hdWGCNA")
library(hdWGCNA)
library(WGCNA)

set.seed(12345)
enableWGCNAThreads(nThreads = 8)

run_hdwgcna_spatial_one_patient <- function(count_mat,
                                            loc_df,
                                            patient_id,
                                            genes_use,
                                            wgcna_name = "hdwgcna_spatial",
                                            min_cells = 50,
                                            k = 25,
                                            max_shared = 10) {
  
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
  seu@meta.data$Subclass <- loc_df$Subclass
  
  seu@meta.data$row <- as.numeric(loc_df$sdimy)
  seu@meta.data$col <- as.numeric(loc_df$sdimx)
  seu@meta.data$imagerow <- as.numeric(loc_df$sdimy)
  seu@meta.data$imagecol <- as.numeric(loc_df$sdimx)
  
  seu <- NormalizeData(seu)
  
  seu <- ScaleData(
    seu,
    features = genes_use
  )
  
  seu <- RunPCA(
    seu,
    features = genes_use
  )
  
  seu <- SetupForWGCNA(
    seu,
    wgcna_name = wgcna_name,
    features = genes_use
  )
  
  seu <- MetacellsByGroups(
    seurat_obj = seu,
    group.by = "Subclass",
    reduction = "pca",
    k = k,
    max_shared = max_shared,
    ident.group = "Subclass",
    assay = "RNA",
    min_cells = min_cells
  )
  
  seu <- NormalizeMetacells(seu)
  
  seu <- SetDatExpr(
    seu,
    assay = "RNA",
    wgcna_name = wgcna_name
  )
  
  seu <- TestSoftPowers(
    seu,
    networkType = "signed",
    wgcna_name = wgcna_name
  )
  
  seu <- ConstructNetwork(
    seu,
    tom_name = patient_id,
    overwrite_tom = TRUE,
    networkType = "signed",
    TOMType = "signed",
    wgcna_name = wgcna_name
  )
  
  seu <- ModuleEigengenes(
    seu,
    wgcna_name = wgcna_name
  )
  
  seu <- ModuleConnectivity(
    seu,
    wgcna_name = wgcna_name
  )
  
  return(seu)
}

#####################################
# all patients
#####################################

############################################################
# 1. Define paths and basic objects
############################################################

tom_dir <- "D:/Research_project/For Test/H20.33.001/coexpression-generation/data_uni/TOM"

gene_names <- readRDS(
  "D:/Research_project/For Test/H20.33.001/coexpression-generation/hdWGCNA/common_genes.RDS"
)

patient_ids <- c(
  "H20.33.035","H21.33.038","H20.33.004","H21.33.023","H21.33.001",
  "H21.33.025","H20.33.001","H21.33.032","H21.33.031","H21.33.019",
  "H21.33.016","H20.33.040","H20.33.012","H20.33.025","H21.33.015",
  "H21.33.005","H21.33.040","H21.33.028","H21.33.012","H20.33.015",
  "H20.33.044","H21.33.022","H21.33.011","H21.33.013","H21.33.014",
  "H21.33.021","H21.33.006"
)

#patient_ids = c("H21.33.006")
#colnames(cell_loc_temp[["H21.33.006"]])[3] = "Subclass"

genes_use <- colnames(express_matrix_temp[[1]])

dir.create(tom_dir, recursive = TRUE, showWarnings = FALSE)

failed_patients <- list()

############################################################
# 2. Run hdWGCNA for all patients
#    Each patient will generate one *_TOM.rda file
############################################################

for (pid in patient_ids) {
  
  message("Running hdWGCNA for patient: ", pid)
  
  tom_file_expected <- file.path(tom_dir, paste0(pid, "_TOM.rda"))

  if (file.exists(tom_file_expected)) {
    message("TOM already exists. Skipping: ", pid)
    next
  }
  
  tryCatch({
    
    this_min_cells <- ifelse(pid == "H20.33.001", 10, 50)
    
    seu <- run_hdwgcna_spatial_one_patient(
      count_mat = express_matrix_temp[[pid]],
      loc_df = cell_loc_temp[[pid]],
      patient_id = pid,
      genes_use = genes_use,
      wgcna_name = "hdwgcna_spatial",
      min_cells = this_min_cells,
      k = 25,
      max_shared = 10
    )
    
    rm(seu)
    gc()
    
  }, error = function(e) {
    
    message("Failed patient: ", pid)
    message("Error: ", e$message)
    failed_patients[[pid]] <<- e$message
    
  })
}

############################################################
# 3. Read all existing TOM files and build list
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
  
  rownames(tom_mat) <- gene_names
  colnames(tom_mat) <- gene_names
  diag(tom_mat) <- 0
  
  hdwgcna_tom_list[[pid]] <- tom_mat
  
  rm(consTomDS, tom_mat)
  gc()
}

############################################################
# 4. Reorder and save
############################################################

available_patients <- intersect(patient_ids, names(hdwgcna_tom_list))
missing_patients <- setdiff(patient_ids, names(hdwgcna_tom_list))

message("Available patients: ", length(available_patients))
message("Missing patients: ", paste(missing_patients, collapse = ", "))

hdwgcna_tom_list <- hdwgcna_tom_list[available_patients]

saveRDS(
  hdwgcna_tom_list,
  file = "D:/Research_project/For Test/H20.33.001/coexpression-generation/hdWGCNA/hdWGCNA_co_express_list_uni.RDS"
)
