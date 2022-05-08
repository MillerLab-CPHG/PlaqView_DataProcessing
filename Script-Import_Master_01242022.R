#### Preprocessing Master ####
## this contains most of the scenario for data import to seurat for plaqview 
## uncomment the import type and modify from the base template!


## General steps
# 1. read data
# 2. update to latest seurat format
# 3. find author provided labels and save as Author_Provided
# 4. export as UNPROCESSED.rds

#### Libraries ####
suppressPackageStartupMessages({
  library(tidyverse) # CRAN
  library(Seurat) # CRAN
  library(patchwork) # CRAN
  library(readr) # CRAN
  library(SingleR) # BIOCONDUCTOR
  library(tidyverse) # CRAN
  library(monocle3) # SPECIFIC INSTALLATION ON WEBSITE
  library(SeuratData) # satijalab/seurat-data
  library(magrittr)# CRAN
  library(ggrepel)# CRAN
  # library(dyno) # devtools::install_github("dynverse/dyno")
  library(SeuratDisk) # remotes::install_github("mojaveazure/seurat-disk")
  library(celldex) # BiocManager::install("celldex")
  library(data.table) # CRAN
  library(matrixStats)# CRAN
  library(Matrix)# CRAN
  # library(bayNorm) # for transposition of sparase matrix
  library(future)
})

# #### From .txt ####
# plaqviewobj <- read.delim("GSE155512_RAW/GSM4705589_RPE004_matrix.txt",  row.names=1)
# plaqviewobj <- CreateSeuratObject(
#   gset1,
#   project = "plaqview_data",
#   assay = "RNA",
#   min.cells = 0,
#   min.features = 0
# )
# 
# #### From multiple .txt ###$
# gset1 <- read.delim("GSE155512_RAW/GSM4705589_RPE004_matrix.txt",  row.names=1)
# gset2 <- read.delim("GSE155512_RAW/GSM4705590_RPE005_matrix.txt",  row.names=1)
# gset3 <- read.delim("GSE155512_RAW/GSM4705591_RPE006_matrix.txt",  row.names=1)
# 
# gset1 <- CreateSeuratObject(
#   gset1,
#   project = "pan_2020",
#   assay = "RNA",
#   min.cells = 0,
#   min.features = 0
# )
# 
# gset2 <- CreateSeuratObject(
#   gset2,
#   project = "pan_2020",
#   assay = "RNA",
#   min.cells = 0,
#   min.features = 0
# )
# 
# gset3 <- CreateSeuratObject(
#   gset3,
#   project = "pan_2020",
#   assay = "RNA",
#   min.cells = 0,
#   min.features = 0
# )
# 
# plaqview_list <- list(gset1, gset2, gset3)
# plaqview_list <- lapply(X = plaqview_list, FUN = function(x) {
#   x <- NormalizeData(x)
#   x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
# })
# 
# # find integration anchors
# features <- SelectIntegrationFeatures(object.list = plaqview_list)
# anchors <- FindIntegrationAnchors(object.list = plaqview_list, anchor.features = features)
# 
# # integrate the samples  
# plaqviewobj <- IntegrateData(anchorset = anchors)
# DefaultAssay(plaqviewobj) <- "integrated" # set defaul to integrated
# 
# #### From .txt with Metadata #### 
# gset1 <- read.delim("PDGFrB _human_counts_2021_09_16.txt",  row.names=1)
# met <- read.delim("metadata_PDGFrB _human_counts_2021_09_16.txt",  row.names=1)
# 
# plaqviewobj <- CreateSeuratObject(
#   gset1,
#   project = "PlaqView_data",
#   assay = "RNA",
#   metadata = met,
#   names.field = 2,
#   min.cells = 0,
#   min.features = 0
# )

# #### From .h5ad ####
# Convert(source = "global_raw.h5ad", dest = "global_raw.h5seurat")
# plaqviewobj <- LoadH5Seurat(file = "global_raw.h5seurat", assays = "RNA")

# #### From .rds ####
# plaqviewobj <- readRDS(file = "pan_from_mingyao_mouse.rds")
# plaqviewobj <- UpdateSeuratObject(plaqviewobj)

# #### From .robj ####
# load("droplet_Heart_and_Aorta_seurat_tiss.Robj")
# plaqviewobj <- UpdateSeuratObject(tiss) # tiss is the name of the robj loaded, may need to change

# #### From .loom ####
#### From 10x cellranger ####
temp1 <- Read10X(data.dir = "filtered_feature_bc_matrix_healthy/")
temp2 <- Read10X(data.dir = "filtered_feature_bc_matrix_stenosis/")

temp1 <- CreateSeuratObject(counts = temp1, project = "PlaqView_healthyXU", min.cells = 3, min.features = 200)
temp2 <- CreateSeuratObject(counts = temp2, project = "PlaqView_stenosisXU", min.cells = 3, min.features = 200)


# #### read and convert ####
# convertedobj <- Connect(filename = "aneursymal-aortic-tissue-human-blood-vessel-10XV3.loom", mode = "r")
# convertedobj <- as.Seurat(convertedobj)


#### OPTIONAL: Calculate nFeatures ####
# some data are missing this portion, esp. .h5ad
plaqviewobj$nCount = colSums(x = plaqviewobj, slot = "counts")  # nCount_RNA
plaqviewobj$nFeature = colSums(x = GetAssayData(object = plaqviewobj, slot = "counts") > 0)  # nFeatureRNA
plaqviewobj@meta.data$nCount_RNA = colSums(x = plaqviewobj, slot = "counts")  # nCount_RNA
plaqviewobj@meta.data$nFeature_RNA = colSums(x = GetAssayData(object = plaqviewobj, slot = "counts") > 0)  # nFeatureRNA

#### RENAME AUTHOR PROVIDED 'Author_Provided" ####
plaqviewobj$Author_Provided <- plaqviewobj$cell_type

# plaqviewobj$Author_Provided <- "Authors did not provide labels"

#### OUTPUT #### 
saveRDS(plaqviewobj, file = "UNPROCESSED.RDS")



