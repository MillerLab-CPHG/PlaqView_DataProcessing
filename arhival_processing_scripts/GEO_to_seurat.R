#### Library and data loading ----
library(Seurat)
library(tidyverse)
library(readr)

#### READ .TXT FROM GEO ####
gset1 <- read.delim("GSE155512_RAW/GSM4705589_RPE004_matrix.txt",  row.names=1)
gset2 <- read.delim("GSE155512_RAW/GSM4705590_RPE005_matrix.txt",  row.names=1)
gset3 <- read.delim("GSE155512_RAW/GSM4705591_RPE006_matrix.txt",  row.names=1)

gset1 <- CreateSeuratObject(
  gset1,
  project = "pan_2020",
  assay = "RNA",
  min.cells = 0,
  min.features = 0
)

gset2 <- CreateSeuratObject(
  gset2,
  project = "pan_2020",
  assay = "RNA",
  min.cells = 0,
  min.features = 0
)

gset3 <- CreateSeuratObject(
  gset3,
  project = "pan_2020",
  assay = "RNA",
  min.cells = 0,
  min.features = 0
)

plaqview_list <- list(gset1, gset2, gset3)
plaqview_list <- lapply(X = plaqview_list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

# find integration anchors
features <- SelectIntegrationFeatures(object.list = plaqview_list)
anchors <- FindIntegrationAnchors(object.list = plaqview_list, anchor.features = features)

# integrate the samples  
plaqviewobj <- IntegrateData(anchorset = anchors)
DefaultAssay(plaqviewobj) <- "integrated" # set defaul to integrated


#### SAVE .RDS ####

saveRDS(plaqviewobj, file = "X_UNPROCESSED.rds")
plaqviewobj <- readRDS(file = "X_UNPROCESSED.rds")

