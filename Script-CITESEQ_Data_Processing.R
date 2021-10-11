#### Library and data loading ----
library(Seurat)
library(patchwork)
library(readr)
library(scCATCH)
library(SingleR)
library(tidyverse)
library(monocle3)
library(SeuratData)
library(magrittr)
library(ggrepel)
library(dyno)
library(readxl)
library(SeuratDisk)

original_color_list <-
  {c("rosybrown2",
     "cadetblue1",
     "lemonchiffon3",
     "darkseagreen",
     "skyblue3",
     "thistle3",
     "cadetblue3",
     "darkseagreen1",
     "palevioletred3",
     "palevioletred1",
     "darkseagreen2",
     "rosybrown3",
     "thistle2",
     "lightsteelblue3",
     "salmon1",
     "palevioletred4",
     "lemonchiffon4",
     "cadetblue2"
  )}

color_function <- colorRampPalette(original_color_list)
manual_color_list <- color_function(40) # change this if clusters >40



#### Read Matrices ####

# Load in the RNA UMI matrix

# Note that this dataset also contains ~5% of mouse cells, which we can use as negative
# controls for the protein measurements. For this reason, the gene expression matrix has
# HUMAN_ or MOUSE_ appended to the beginning of each gene.
CITESEQ.RNA <- as.sparse(read.csv(file = "../DataProcessing/data/Fernandez_2019/source_files/citeseq_pbmc_gex_macrophages.csv", sep = ",",
                                  header = TRUE))


# To make life a bit easier going forward, we're going to discard all but the top 100 most
# highly expressed mouse genes, and remove the 'HUMAN_' from the CITE-seq prefix
# cbmc.rna <- CollapseSpeciesExpressionMatrix(cbmc.rna)

# Load in the ADT UMI matrix
CITESEQ.ADT <- as.sparse(read_csv("data/Fernandez_2019/source_files/citeseq_pbmc_adt_macrophages.csv", 
                                  col_types = cols(...1 = col_character())))



# Note that since measurements were made in the same cells, the two matrices have identical
# column names
all.equal(colnames(CITESEQ.RNA), colnames(CITESEQ.ADT))
# here most are the same just an error in the first 
# overriding

# creates a Seurat object based on the scRNA-seq data
CITESEQ.OBJ <- CreateSeuratObject(assay = "RNA", counts = CITESEQ.RNA)


# create a new assay to store ADT information
adt_assay <- CreateAssayObject(counts = CITESEQ.ADT)

# add this assay to the previously created Seurat object
CITESEQ.OBJ[["ADT"]] <- adt_assay

# Validate that the object now contains multiple assays
CITESEQ.OBJ@assays$RNA

# Extract a list of features measured in the ADT assay
rownames(CITESEQ.OBJ[["ADT"]])




