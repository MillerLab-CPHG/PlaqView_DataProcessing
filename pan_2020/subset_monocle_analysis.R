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


plaqviewobj <- readRDS(file = "Pan_2020.rds")
plaqviewobj.cds <- readRDS(file = "Pan_2020_CDS.rds")

#### MONOCLE3 TRAJECTORY INFERENCE ----
# in previous versions we tried the seurat wrapper it just didnt work
# below we manually wrap the data ourselves

# convert to monocle cds object 
# Extract data, phenotype data, and feature data from the SeuratObject
expressiondata <- plaqviewobj@assays[["RNA"]]@data

cellmd <- plaqviewobj@meta.data

genemd <- data.frame(gene_short_name = row.names(expressiondata), 
                     row.names = row.names(expressiondata))

# Construct monocle cds
plaqviewobj.cds <- new_cell_data_set(expression_data = expressiondata,
                                     cell_metadata = cellmd,
                                     gene_metadata = genemd)
plaqviewobj.cds <- preprocess_cds(plaqviewobj.cds, num_dim = 30) # we used 30 in earlier seurat scripts

# 
# run clustering again (didnt transfer from seurat)
plaqviewobj.cds <- reduce_dimension(plaqviewobj.cds, reduction_method = "UMAP")
plaqviewobj.cds <- cluster_cells(plaqviewobj.cds, reduction_method = "UMAP")


#### TRANSFER SEURAT EMBEDDINGS #####
# Note that these may be calculated on the Integrated object, not the counts
#   and thus will involve fewer genes
temp.cds <- ProjectDim(plaqviewobj, reduction = "pca") # this will be removed
reducedDim(plaqviewobj.cds, type = "PCA") <- temp.cds@reductions$pca@cell.embeddings
plaqviewobj.cds@preprocess_aux$prop_var_expl <- temp.cds@reductions$pca@stdev
plot_pc_variance_explained(plaqviewobj.cds)

# Transfer Seurat UMAP embeddings
plaqviewobj.cds@int_colData@listData$reducedDims$UMAP <- temp.cds@reductions$umap@cell.embeddings

## transfer singleR labels to moncle3 object
colData(plaqviewobj.cds)$assigned_cell_type <- plaqviewobj@meta.data[["SingleR.calls"]] # call this by opening the object

#### MONOCLE3 CONT. ----
# now learn the PATH (trajectory)
plaqviewobj.cds <- learn_graph(plaqviewobj.cds)

# this calls up a shiny app, choose the ROOT NODE
plaqviewobj.cds <- order_cells(plaqviewobj.cds, reduction_method = "UMAP")

# finally, you can visualize the learned path
pdf("monocle3_RNAvelocity_seuratpartition.pdf", width=6, height=6)
plot_cells(plaqviewobj.cds,
           color_cells_by = "assigned_cell_type",
           label_groups_by_cluster=F,
           show_trajectory_graph = T,
           trajectory_graph_segment_size = 1,
           label_leaves=F, # this gives a little node label (outcome)
           label_roots = T,
           label_branch_points = F,
           graph_label_size = 1, # size of # in circle
           group_label_size = 3,
           cell_size = 1,
           alpha = 0.7,
           scale_to_range = T) +
  scale_color_manual(values = manual_color_list) # sync color scheme

dev.off()

mon3 <- plot_cells(plaqviewobj.cds,
                   color_cells_by = "assigned_cell_type",
                   label_groups_by_cluster=F,
                   show_trajectory_graph = T,
                   trajectory_graph_segment_size = 1,
                   label_leaves=F, # this gives a little node label (outcome)
                   label_roots = T,
                   label_branch_points = F,
                   graph_label_size = 1, # size of # in circle
                   group_label_size = 3,
                   cell_size = 1,
                   alpha = 0.7,
                   scale_to_range = T) +
  scale_color_manual(values = manual_color_list) # sync color scheme

saveRDS(mon3, file = "dyno/monocle3.rds")
# now you can show pseudotime
pdf("monocle3_pseudotime_seuratpartition.pdf", width=7, height=6)
plot_cells(plaqviewobj.cds,
           color_cells_by = "pseudotime",
           show_trajectory_graph = F,
           trajectory_graph_segment_size = 1,
           label_leaves=F, # this gives a little node label (outcome)
           label_roots = T,
           label_branch_points = F,
           graph_label_size = 1, # size of # in circle
           group_label_size = 3,
           cell_size = 1,
           alpha = 0.7,
           scale_to_range = T) 
dev.off()


# #### NOT RUN Subset Trajectory & analysis of SMC----
# plaqviewobj.cds_subset <- choose_cells(plaqviewobj.cds) # calls up shiny app
# 
# plot_cells(plaqviewobj.cds_subset,
#            color_cells_by = "pseudotime",
#            show_trajectory_graph = T,
#            trajectory_graph_segment_size = 1,
#            label_leaves=F, # this gives a little node label (outcome)
#            label_roots = T,
#            label_branch_points = F,
#            graph_label_size = 1, # size of # in circle
#            group_label_size = 3,
#            cell_size = 1,
#            alpha = 0.7,
#            scale_to_range = T) 
# 
# #### NOT RUN MORAN's I Test of Autocorrelation ####
# now we can extrapolate genes that are differentially expressed in this region
# Moran’s I is a measure of multi-directional and multi-dimensional spatial autocorrelation. 
# the statistic tells you whether cells at nearby positions on a 
# trajectory will have similar (or dissimilar) +
# expression levels for the gene being tested.
## first lets do the whole dataset
# a special gene module score heatmap (for the whole dataset)
# pr_graph_test_res <- graph_test(plaqviewobj.cds, neighbor_graph="principal_graph", cores=2)
write.csv(pr_graph_test_res, file = "moransI_all_clusters.csv")
pr_deg_ids <- row.names(subset(pr_graph_test_res, q_value < 0.00000001)) # you can adjust the p-value here
head(pr_deg_ids)
gene_module_df <- find_gene_modules(plaqviewobj.cds[pr_deg_ids,], resolution=1e-3)
cell_group_df <- tibble::tibble(cell=row.names(colData(plaqviewobj.cds)), 
                                cell_group=colData(plaqviewobj.cds)$assigned_cell_type)
agg_mat <- aggregate_gene_expression(plaqviewobj.cds, gene_module_df, cell_group_df)
row.names(agg_mat) <- stringr::str_c("Module ", row.names(agg_mat))
pheatmap::pheatmap(agg_mat,
                   scale="column", clustering_method="ward.D2")

# which then can be visualized like so;
# this can show you the different gene modules that can are responsible for changes over pseudotime
plot_cells(plaqviewobj.cds,
           genes=gene_module_df %>% filter(module %in% c(2,3,7)), # specify the module you want to examine
           label_cell_groups=T,
           show_trajectory_graph=F)

subset(gene_module_df, module == 2)

## now lets do the subsets
# pr_graph_test_res.sub <- graph_test(plaqviewobj.cds_subset, neighbor_graph="principal_graph", cores=2)
pr_deg_ids.sub <- row.names(subset(pr_graph_test_res.sub, q_value < 0.00000001))
write.csv(pr_graph_test_res.sub, file = "moransI_subset_cluster.csv")
head(pr_deg_ids.sub)

# collect the trajectory-variable genes into modules
gene_module_df.sub <- find_gene_modules(plaqviewobj.cds_subset[pr_deg_ids.sub,], resolution=1e-3)
# visualize these genes
# here I am just pulling out genes that have high moran's i and might be helpful in the paper
# SELECTED FOR PUBLICATIONS
pdf("monocle3_genesoverpseudotime_seuratpartition_extended.pdf", width=7, height=6)
plot_cells(plaqviewobj.cds_subset, 
           genes=c("MYH11", 'IGFBP2',"PPP1R14A","CNN1", "TNFRSF11B",
                   "C7", "C3",
                   "SERPINF1",  "FBLN1", 
                   "CXCL12", "MMP2", 
                   "FN1"), # this is faceting by the genes that are DE
           show_trajectory_graph=FALSE, 
           label_cell_groups=F, cell_size = 1)

dev.off()

# recluster at higher definition
plaqviewobj.cds_subset = cluster_cells(plaqviewobj.cds_subset, resolution=1e-2)

pdf("monocle3_RNAvelocitySUBSET_seuratpartition.pdf", width=6, height=6)
plot_cells(plaqviewobj.cds_subset, 
           color_cells_by="cluster",
           label_groups_by_cluster=F,
           show_trajectory_graph = T,
           trajectory_graph_segment_size = 1,
           label_leaves=F, # this gives a little node label (outcome)
           label_roots = F,
           label_branch_points = F,
           graph_label_size = 1, # size of # in circle
           group_label_size = 4,
           cell_size = 1,
           alpha = 0.5,
           scale_to_range = T)
dev.off()


#### rename manual label####

plaqviewobj[["manually_annotated_labels"]] <- "ORIGINAL AUTHORS HAVE NOT SUPPLIED ANNOTATIONS, PLEASE USE SINGLER and OTHER LABELING METHODS"


#### DYNO TRAJECTORY INFERENCES ####
object_counts <- Matrix::t(as(as.matrix(plaqviewobj@assays$RNA@counts), 'sparseMatrix'))
object_expression <- Matrix::t(as(as.matrix(plaqviewobj@assays$RNA@data), 'sparseMatrix'))
object_cellinfo <- plaqviewobj@meta.data[["SingleR.labels"]]

plaqviewobj.dyno <- wrap_expression(
  counts = object_counts,
  expression = object_expression)


#### slingshot: construct the model ####
# make sure to call up docker images

model <- infer_trajectory(plaqviewobj.dyno, "slingshot", verbose = T)

#### slingshot: project the model ###
# add dim reduction
model <- model %>% 
  add_dimred(dimred = as.matrix(plaqviewobj@reductions$umap@cell.embeddings),
             expression_source = plaqviewobj.dyno$expression)

pdf("dyno/dyno_slingshot_full.pdf", width=7, height=6)
slingshot <- plot_dimred(
  model, 
  expression_source = plaqviewobj.dyno$expression,
  grouping = object_cellinfo # basically stanford@meta.data[["SingleR.labels"]]
)

saveRDS(slingshot, file = "dyno/slingshot.rds")
slingshot
dev.off()

#### slingshot: show a gene expression
plot_dimred(
  model, 
  expression_source = plaqviewobj.dyno$expression, 
  feature_oi = "FN1"
)

#### scorpius: construct the model ####
# make sure to call up docker images

model <- infer_trajectory(plaqviewobj.dyno, "scorpius")

#### scorpius: project the model ###
# add dim reduction
model <- model %>% 
  add_dimred(dimred = as.matrix(plaqviewobj@reductions$umap@cell.embeddings),
             expression_source = plaqviewobj.dyno$expression)

pdf("dyno/dyno_scorpius_full.pdf", width=7, height=6)
scorpius <- plot_dimred(
  model, 
  expression_source = plaqviewobj.dyno$expression,
  grouping = object_cellinfo # basically stanford@meta.data[["SingleR.labels"]]
)

saveRDS(scorpius, file = "dyno/scorpius.rds")
scorpius
dev.off()

#### PAGA: construct the model ####
model <- infer_trajectory(plaqviewobj.dyno, "projected_paga", verbose = T)

#### PAGA: project the model ###
# add dim reduction
model <- model %>% 
  add_dimred(dimred = as.matrix(plaqviewobj@reductions$umap@cell.embeddings),
             expression_source = plaqviewobj.dyno$expression)

paga <- plot_dimred(
  model, 
  expression_source = plaqviewobj.dyno$expression, 
  grouping = object_cellinfo # basically stanford@meta.data[["SingleR.labels"]]
)
paga
saveRDS(paga, file = "dyno/paga.rds")


#### REDUCE SIZE & OUTPUT ####
plaqviewobj <- DietSeurat(plaqviewobj, counts = T, data = T, dimreducs = c('umap'))

saveRDS(plaqviewobj, file = "Pan_2020.rds")
saveRDS(plaqviewobj.cds, file = "Pan_2020_CDS.rds")


plaqviewobj <- readRDS(file = "Pan_2020.rds")
plaqviewobj.cds <- readRDS(file = "Pan_2020_CDS.rds")

#### DIFF EX GENE LIST ####
Idents(object = plaqviewobj) <- "SingleR.calls"
difflist <- Seurat::FindAllMarkers(plaqviewobj)
write_csv(difflist, file = "differential/diff_by_singleR.csv")

Idents(object = plaqviewobj) <- "manually_annotated_labels"
singlerdifflist <- Seurat::FindAllMarkers(plaqviewobj)
write_csv(difflist, file = "differential/diff_by_author.csv")

Idents(object = plaqviewobj) <- "seurat_clusters"
singlerdifflist <- Seurat::FindAllMarkers(plaqviewobj)
write_csv(difflist, file = "differential/diff_by_seurat.csv")

