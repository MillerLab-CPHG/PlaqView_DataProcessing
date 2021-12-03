library(SCORPIUS)
library(Seurat)



plaqviewobj <- NormalizeData(plaqviewobj)

expression <- Matrix::t(Matrix(plaqviewobj@assays$RNA@data, sparse = T))
group_name <- plaqviewobj@meta.data$
  
space <- reduce_dimensionality(expression, dist = "spearman", ndim = 15)
draw_trajectory_plot(space, progression_group = group_name, contour = TRUE)
