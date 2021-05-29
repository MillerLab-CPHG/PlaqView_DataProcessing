#### Library and data loading ----
library(Seurat)
library(tidyverse)
library(readr)

#### READ .TXT FROM GEO ####
gset1 <- read.delim("GSE131776_mouse_scRNAseq_wirka_et_al_GEO.txt",  row.names=1)

gset1 <- CreateSeuratObject(
  gset1,
  project = "wirka_mouse_2019",
  assay = "RNA",
  min.cells = 0,
  min.features = 0
)


#### SAVE .RDS ####

saveRDS(gset1, file = "UNPROCESSED.rds")
gset1 <- readRDS(file = "UNPROCESSED.rds")

