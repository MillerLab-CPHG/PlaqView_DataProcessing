# devtools::install_github("atakanekiz/CIPR-Package", build_vignettes = T)
library(Seurat)



library(CIPR)
plaqviewobj <- NormalizeData(plaqviewobj)

allmarkers <- FindAllMarkers(plaqviewobj)
avgexp <- AverageExpression(plaqviewobj)

# Plot summarizing top scoring references per cluster (logFC comparison)
CIPR(input_dat = allmarkers,
     comp_method = "logfc_dot_product", 
     reference = "immgen", 
     plot_ind = F,
     plot_top = T)

# Plot summarizing top scoring references per cluster (all-genes correlation)
CIPR(input_dat = allmarkers,
     comp_method = "logfc_dot_product", 
     reference = "immgen", 
     plot_ind = F,
     plot_top = T)


# Plots for individual clusters
CIPR(input_dat = allmarkers,
     comp_method = "logfc_dot_product", 
     reference = "immgen", 
     plot_ind = T,
     plot_top = F)

# Limiting the analysis to certain reference subsets
CIPR(input_dat = allmarkers,
     comp_method = "logfc_dot_product", 
     reference = "immgen", 
     plot_ind = F,
     plot_top = T, 
     select_ref_subsets = c("T cell", "B cell", "NK cell"))



