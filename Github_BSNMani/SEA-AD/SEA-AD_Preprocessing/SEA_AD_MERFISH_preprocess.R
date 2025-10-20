##############
## GL environment
## module load Bioinformatics
## module load Rgiotto
##############
rm(list=ls())
library(tidyr)
library(dplyr)
library(readxl)
library(data.table)
library(WGCNA)
library(DescTools)
library(olsrr)
library(gplots)
library(Giotto,lib.loc = "/sw/pkgs/med/Rgiotto/1.1.0")

seaad_dir = "/nfs/turbo/umms-lgarmire/liyijun/BSNMani_ST_Application_Project/SEA_AD"## your data root directory
setwd(seaad_dir)
patient_IDs_merFISH = list.files("raw_data/middle-temporal-gyrus") ##MERFISH patient ID
print(getwd())

############################
### preprocess clinical data
############################
meta = read_excel("raw_data/sea-ad_cohort_donor_metadata_072524.xlsx")
meta_merFISH = meta %>% filter(`Donor ID` %in% patient_IDs_merFISH)

##### select predictors
clinical_concensus_variables_ID = which(sapply(colnames(meta_merFISH),FUN=function(x){grepl(x,pattern="Consensus")}))
clinical_concensus_variables = colnames(meta_merFISH)[clinical_concensus_variables_ID]
meta_merFISH$Clinical_Concensus = 
  apply(meta_merFISH[,clinical_concensus_variables_ID],1,FUN=function(x){diagnosis = clinical_concensus_variables[which(x=="Checked")]; diagnosis = paste(diagnosis,collapse = "; "); if(!is.na(x[length(x)])){diagnosis = paste(diagnosis,x[length(x)],sep="; ")}; return(diagnosis)})
meta_merFISH = meta_merFISH %>% mutate(Age = `Age at Death` - mean(`Age at Death`),
                                       Sex = as.numeric(Sex=="Female"),
                                       Dementia = as.numeric(`Consensus Clinical Dx (choice=Control)`=="Unchecked"),
                                       Dementia_v2 = as.numeric(`Cognitive Status`=="Dementia"),
                                       AD_change = as.numeric(`Overall AD neuropathological Change`=="Not AD"),
                                       Braak = as.numeric(factor(Braak,levels=c("Braak 0","Braak II","Braak III","Braak IV","Braak V","Braak VI")))-1,
                                       CERAD = as.numeric(`CERAD score`!="Absent"),
                                       CAA = as.numeric(`Overall CAA Score`!="Not identified"),
                                       Atherosclerosis = as.numeric(factor(Atherosclerosis,levels=c("None","Mild","Moderate","Severe")))-1,
                                       Arteriolosclerosis = as.numeric(factor(Arteriolosclerosis,levels=c("Mild","Moderate","Severe"))),
                                       brain_weight = as.numeric(`Fresh Brain Weight`),
                                       yrs_education = `Years of education`,
                                       brain_ph = `Brain pH`,
                                       n_microinfarcts_screening = `Total microinfarcts in screening sections`,
                                       MMSE = `Last MMSE Score`,
                                       LATE = as.numeric(factor(LATE,levels=c("Not Identified","LATE Stage 1","LATE Stage 2","LATE Stage 3")))-1)
meta_merFISH$Thal = unname(sapply(meta_merFISH$Thal, FUN=function(x){as.numeric(strsplit(x,split=" ")[[1]][2])}))

rownames(meta_merFISH) = meta_merFISH$`Donor ID`
meta_merFISH = meta_merFISH[patient_IDs_merFISH,]
rownames(meta_merFISH) = NULL
print(getwd())
saveRDS(meta_merFISH,file="preprocessed/merFISH_clinical.RDS")

##### check association of clinical covariates with MMSE; feel free to explore other variables as clinical outcomes as well
mod = lm(MMSE~yrs_education+Age+Dementia+Dementia_v2+Sex+brain_weight+AD_change+Thal+Braak+CERAD+CAA+Atherosclerosis+Arteriolosclerosis+brain_ph+n_microinfarcts_screening+LATE,data=meta_merFISH)
var_selection_res = ols_step_all_possible(mod)
saveRDS(var_selection_res,file="preprocessed/merFISH_clinical_predictor_var_selection.RDS")

## check the stepwise selection performance, limit clinical predictors to 3 or less, since sample size is small, prioritize adjr, rmse, fpe, aic
var_selection_res$result %>% filter(n<4) %>% arrange(desc(adjr)) %>% head()
var_selection_res$result %>% filter(n<4) %>% arrange(rmse) %>% head()
var_selection_res$result %>% filter(n<4) %>% arrange(aic) %>% head()
var_selection_res$result %>% filter(n<4) %>% arrange(fpe) %>% head() ## predictive
##  Sex, Atherosclerosis, Arteriolosclerosis

var_selection_res$result %>% filter(n<3) %>% arrange(desc(adjr)) %>% head()
var_selection_res$result %>% filter(n<3) %>% arrange(rmse) %>% head()
var_selection_res$result %>% filter(n<3) %>% arrange(aic) %>% head()
var_selection_res$result %>% filter(n<3) %>% arrange(fpe) %>% head() ## predictive

############################
### preprocess gene expression data
############################
instr = createGiottoInstructions(python_path="/sw/pkgs/arc/python3.9-anaconda/2021.11/bin")
nCells_vec = c()
GO_obj_ls = list()
genes_ls = list()
raw_exprs_ls = list()
spatial_locs_ls = list()
sample_ID_vec = c()
#### my current code is only using one specimen from each donor (patient), feel free to explore all the specimens.
for(i in 1:length(patient_IDs_merFISH)){
  ## select patient sample
  patient_ID = patient_IDs_merFISH[i]
  patient_samples = list.files(fs::path("raw_data/middle-temporal-gyrus",patient_ID))
  sample_ID = patient_samples[1]
  sample_ID_vec = c(sample_ID_vec,sample_ID)
  
  ## process raw expression
  df = fread(fs::path("raw_data/middle-temporal-gyrus",patient_ID,sample_ID,"cellpose-detected_transcripts",ext="csv"))
  
  ## Modified by TONG 4/5/2025. The first column has no meaning
  df = df %>% select(-1)
  ## Column 1 is not count number
  #colnames(df)[1] = "counts"
  df = df %>% filter(cell_id != -1)
  df$cell_id = as.character(df$cell_id)
  
  ## aggregate raw expression  NEED TO USE SUM(COUNTS) INSTEAD OF MEAN(COUNTS) BECAUSE OF LIMITATION OF SpaceX Model！！
  ## change from counts = mean(counts) to counts=n() By TONG 4/5/2025
  aggregated_df = df %>% group_by(gene,cell_id) %>% summarise(counts=n(),barcode_id=mode(barcode_id),global_x=mean(global_x),global_y=mean(global_y),global_z=mean(global_z),x=mean(x),y=mean(y),fov=mode(fov))
  raw_exprs_df = aggregated_df %>% select(counts,gene,cell_id) %>% spread(key=cell_id,value=counts)
  raw_exprs = as.matrix(raw_exprs_df[,-1])
  rownames(raw_exprs) = raw_exprs_df$gene
  raw_exprs[which(is.na(raw_exprs))] = 0
  
  spatial_locs = aggregated_df %>% ungroup() %>% group_by(cell_id) %>% summarise(sdimx=mean(x),sdimy=mean(y))
  colnames(spatial_locs)[which(colnames(spatial_locs)=="cell_id")] = "cell_ID"

  ## build Giotto object
  GO_obj = createGiottoObject(raw_exprs=raw_exprs,spatial_locs=spatial_locs,instructions=instr)
  
  ## gene filtering
  gene_IDs = rownames(raw_exprs)
  mito_gs = gene_IDs[grep(pattern="(^MT-)|(^mt-)",gene_IDs)]
  print(paste("removing",length(mito_gs),"mitochondrial genes",sep=" "))
  blank_gs = gene_IDs[grep(pattern="(^Blank)|(blank)",gene_IDs)]
  print(paste("removing",length(blank_gs),"blank genes",sep=" "))
  negprb_gs = gene_IDs[grep(gene_IDs,pattern="(^NegPrb)")]
  print(paste("removing",length(negprb_gs),"negative probes",sep=" "))
  gene_ID_keep = setdiff(gene_IDs, Reduce(union,list(mito_gs,blank_gs,negprb_gs)))
  print(gene_ID_keep)
  GO_obj_filtered = subsetGiotto(GO_obj, gene_ids = gene_ID_keep)
  print(str(GO_obj_filtered@raw_exprs))
  
  ## cell filtering
  genes_expressed_per_cell = colSums(as.matrix(GO_obj_filtered@raw_exprs)>0)
  cell_ID_keep = names(genes_expressed_per_cell[which(genes_expressed_per_cell > 20)])
  GO_obj_filtered = subsetGiotto(GO_obj_filtered, cell_ids = cell_ID_keep)

  ## normalize and scale
  GO_obj_filtered = normalizeGiotto(gobject=GO_obj_filtered,scalefactor=10000)
  
  raw_exprs_ls[[i]] = GO_obj_filtered@raw_exprs
  GO_obj_ls[[i]] = GO_obj_filtered@norm_scaled_expr
  genes_ls[[i]] = rownames(GO_obj_filtered@raw_exprs)
  spatial_locs_ls[[i]] = GO_obj_filtered@spatial_locs
  ## By Tong 20/4/2025
  names(spatial_locs_ls)[i] = patient_ID
  #names(GO_obj_ls)[i] = patient_ID
  
  nCells_vec = c(nCells_vec,length(unique(df$cell_id)))
}
saveRDS(raw_exprs_ls, file = "preprocessed/raw_exprs_merFISH.RDS")
saveRDS(GO_obj_ls, file = "preprocessed/norm_scaled_exprs_merFISH.RDS")
saveRDS(genes_ls, file = "preprocessed/patient_genes_merFISH.RDS")
saveRDS(spatial_locs_ls, file = "preprocessed/spatial_locs_merFISH.RDS")
cat("spatial_locs_merFISH file saved successfully")

meta_2 = data.frame(patient=patient_IDs_merFISH, sample=sample_ID_vec, nCells = nCells_vec, platform="merFISH")
saveRDS(meta_2, file="preprocessed/patient_sample_ID.RDS")

############################
### build WGCNA gene co-expression networks
############################
wgcna_raw_ls = wgcna_transformed_ls = list()
for(i in 1:length(patient_IDs_merFISH)){
  ## construct raw network
  raw_corr = adjacency(datExpr = t(GO_obj_ls[[i]][genes_ls[[1]],]), type="signed", power=1)
  raw_corr = 2*raw_corr-1
  diag(raw_corr) = 0
  ## normalize raw network
  transformed_corr = FisherZ(raw_corr)
  
  wgcna_raw_ls[[i]] = raw_corr
  wgcna_transformed_ls[[i]] = transformed_corr
}
saveRDS(wgcna_raw_ls, file = "preprocessed/WGCNA_raw_merFISH.RDS")
saveRDS(wgcna_transformed_ls, file = "preprocessed/WGCNA_transformed_merFISH.RDS")

## visualize one of the networks (without context of how genes are clustered)
col.heat = colorRampPalette(c("blue", "white", "red"))(256)
plot_heatmap = function(x,color,mask=TRUE){
  mat = x%*%t(x)
  if(mask){
    diag(mat)=NA
  }
  heatmap.2(mat, trace="none", scale="none", na.color="grey", Rowv=NULL, Colv=NULL, col=color, density.info="none", dendrogram="none")
}

png(fs::path(paste(patient_IDs_merFISH[1],meta_2$sample[1],"FisherZ",sep="_"),ext="png"))
heatmap.2(wgcna_transformed_ls[[1]], trace="none", scale="none", na.color="grey", Rowv=NULL, Colv=NULL, col=col.heat, density.info="none", dendrogram="none")
dev.off()
