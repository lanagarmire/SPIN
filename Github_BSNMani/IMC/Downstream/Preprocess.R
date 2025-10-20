source("/home/yangiwen/n/workspace/common/Utils.R")

# Yijun's preprocessed protein data
data <- readRDS("/home/yangiwen/n/workspace/bsnmani/data/extern/bc_imc/preprocessed/Basel_patient_metal_count.RDS")

# Yijun's preprocessed clinical data
clinical <- readRDS("/home/yangiwen/n/workspace/bsnmani/data/extern/bc_imc/preprocessed/clinical_data.RDS")

# Spatial coordinates from Jackson's raw data
spatial <- read.csv("/home/yangiwen/n/workspace/bsnmani/data/bc/Basel_SC_locations.csv")

# Cell type cluster from Jackson's raw data
cluster <- read.csv("/home/yangiwen/n/workspace/bsnmani/data/bc/Cluster_labels/Basel_metaclusters.csv")

# Pick a random patient
patient <- "BaselTMA_SP41_257_X3Y1"

# Subset patient protein data
patient_data <- data[[patient]]
patient_data <- as.data.frame(t(patient_data))
# Make first row colnames
colnames(patient_data) <- as.character(unlist(patient_data[1, ]))
patient_data <- patient_data[-1, ]

# Subset spatial coords
patient_spatial <- spatial[spatial$core == patient, ]
rownames(patient_spatial) <- patient_spatial$id

# Store in Seurat object
obj <- Seurat::CreateSeuratObject(t(patient_data), meta.data = patient_spatial)
saveRDS(obj, paste0("../output/", patient, ".rds"))

# Store in separate files
saveRDS(patient_data, "../output/patient_data.rds")
write.csv(patient_spatial, "../output/patient_spatial.csv")

# Draw violin plot for all proteins
png("../output/vln.png", width = 1024, height = 4096, res = 120)
print(Seurat::VlnPlot(obj, features = rownames(obj), pt.size = 0))
dev.off()
