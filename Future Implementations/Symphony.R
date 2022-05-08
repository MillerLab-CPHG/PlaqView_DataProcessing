# install.packages("symphony")

library(symphony)

suppressPackageStartupMessages({
  library(symphony)
  library(singlecellmethods) #devtools::install_github("immunogenomics/singlecellmethods")
  library(tidyverse)
  library(data.table)
  library(matrixStats)
  library(Matrix)
  library(plyr)
  library(dplyr)
  
  # Plotting
  library(ggplot2)
  library(ggthemes)
  library(ggrastr)
  library(RColorBrewer)
  library(patchwork)
})

fig.size <- function (height, width) {
  options(repr.plot.height = height, repr.plot.width = width)
}

plaqviewobj <- readRDS(file = "data/Alencar_2020/Alencar_2020.rds")
ref = readRDS(file = "references/Symphony_ref_data/fibroblast_atlas.rds")

fig.size(5, 5)
p = plotReference(ref,
                  as.density = TRUE,      # plot density or individual cells
                  bins = 14,              # if density, nbins parameter for stat_density_2d
                  bandwidth = 1,        # if density, bandwidth parameter for stat_density_2d
                  title = "Symphony Reference: 10x PBMCs",    # Plot title
                  color.by = 'cell_type', # metadata column name for cell type labels
                  celltype.colors = pbmc_colors, # custom color palette
                  show.legend = TRUE,     # Show cell type legend
                  show.labels = TRUE,     # Show cell type labels
                  show.centroids = FALSE) # Plot soft cluster centroid locations)
p


query = symphony::mapQuery(
  plaqviewobj[['RNA']]@counts,
  plaqviewobj@meta.data,
  ref, 
  vars = 'Author_Provided',  # use column names from your meta_data
  do_normalize = TRUE)

colnames(query$umap) = c('UMAP1', 'UMAP2')
umap_labels = cbind(query$meta_data, query$umap)

