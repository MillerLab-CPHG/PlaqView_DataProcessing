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


#### unzip from GEO ####
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


#### SEURAT: QC and reduction---- 
#### INTEGRATED DATASET PROCESSING ####
# Run the standard workflow for visualization and clustering
plaqviewobj <- ScaleData(plaqviewobj, verbose = FALSE)
plaqviewobj <- RunPCA(plaqviewobj, npcs = 30, verbose = FALSE)
plaqviewobj <- RunUMAP(plaqviewobj, reduction = "pca", dims = 1:30)
plaqviewobj <- FindNeighbors(plaqviewobj, reduction = "pca", dims = 1:30)
plaqviewobj <- FindClusters(plaqviewobj, resolution = 0.5)

#### SINGLE-R LABELING (post Seurat automated cell cluster annotation alternative) ----
# BiocManager::install("SingleR")
# here we are using Human Primary Cell Atlas design for blood
# https://bioconductor.org/packages/3.12/data/experiment/vignettes/celldex/inst/doc/userguide.html#2_General-purpose_references
hpca.se <- celldex::HumanPrimaryCellAtlasData() # build the reference
hpca.se

# now run the prediction using the reference
# singleR requires that it be in a 'singlecellexperiment' format
# they are workout agnostic

for_singleR_input <- GetAssayData(plaqviewobj)
pred.plaqviewobj <- SingleR(test = for_singleR_input, 
                         ref = hpca.se, 
                         label = hpca.se$label.main) # reference cell types
pred.plaqviewobj
# summarize distribution
table(pred.plaqviewobj$labels)

# to show annotation confidence map
plotScoreHeatmap(pred.plaqviewobj)

# to show # that are pruned due to low score
summary(is.na(pred.plaqviewobj$pruned.labels))

### to place the singleR predictions into Seurat as a sep unit ###
# seurat.obj[["SingleR.labels"]] <- singler.results$labels
plaqviewobj[["SingleR.labels"]] <- pred.plaqviewobj$labels # this nest under metadata

# Copy over the labels and pruned.labels (Note: any other column of the results could be used as well)
plaqviewobj$SingleR.pruned.calls <- pred.plaqviewobj$pruned.labels
plaqviewobj$SingleR.calls <- pred.plaqviewobj$labels

#### RECODE SINGLE-R NAMES ----
plaqviewobj@meta.data[["SingleR.calls"]] <- recode(plaqviewobj@meta.data[["SingleR.calls"]], Smooth_muscle_cells = "SMC")
plaqviewobj@meta.data[["SingleR.calls"]] <- recode(plaqviewobj@meta.data[["SingleR.calls"]], Endothelial_cells = "EC")
plaqviewobj@meta.data[["SingleR.calls"]] <- recode(plaqviewobj@meta.data[["SingleR.calls"]], NK_cell = "NK")
plaqviewobj@meta.data[["SingleR.calls"]] <- recode(plaqviewobj@meta.data[["SingleR.calls"]], Chondrocytes = "CH")
plaqviewobj@meta.data[["SingleR.calls"]] <- recode(plaqviewobj@meta.data[["SingleR.calls"]], Fibroblasts = "FB")
plaqviewobj@meta.data[["SingleR.calls"]] <- recode(plaqviewobj@meta.data[["SingleR.calls"]], Monocyte = "Mono")
plaqviewobj@meta.data[["SingleR.calls"]] <- recode(plaqviewobj@meta.data[["SingleR.calls"]], B_cell = "B_Cells")
plaqviewobj@meta.data[["SingleR.calls"]] <- recode(plaqviewobj@meta.data[["SingleR.calls"]], Macrophage = "Mø")
plaqviewobj@meta.data[["SingleR.calls"]] <- recode(plaqviewobj@meta.data[["SingleR.calls"]], Tissue_stem_cells = "SC")
plaqviewobj@meta.data[["SingleR.calls"]] <- recode(plaqviewobj@meta.data[["SingleR.calls"]], T_cells = "T_Cells")
plaqviewobj@meta.data[["SingleR.calls"]] <- recode(plaqviewobj@meta.data[["SingleR.calls"]], 'Pre-B_cell_CD34-' = "PreB_CD34-")
plaqviewobj@meta.data[["SingleR.calls"]] <- recode(plaqviewobj@meta.data[["SingleR.calls"]], 'Pro-B_cell_CD34+' = "ProB_CD34+")


table(plaqviewobj@meta.data[["SingleR.calls"]])


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

#### STACKED POPULATION PLOT ####
# Idents(plaqviewobj) <- plaqviewobj@meta.data[["SingleR.calls"]]
# pop1 <- as.data.frame(prop.table(table(Idents(plaqviewobj))))
# pop1$Method <- "SingleR"
# 
# Idents(plaqviewobj) <- plaqviewobj@meta.data[["manually_annotated_labels"]]
# pop2 <- as.data.frame(prop.table(table(Idents(plaqviewobj))))
# pop2$Method <- "ManualCluster"
# 
# plot1 <- ggplot(pop1, aes(y = Freq, x =Method,)) + 
#   geom_bar(position="stack", stat="identity",
#            fill = manual_color_list,
#            width = 0.1) + 
#   xlab("") +
#   ylab("Cell Porportions") +
#   geom_label_repel(aes(label = Var1), colour = "black",
#             position = position_stack(vjust = 0.5),
#             max.overlaps = 15 ,
#             force = 10,
#             force_pull = 3,
#             max.iter = 999999,
#             max.time = 1,
#             min.segment.length = 0.01,
#             xlim = c(1, 3)) +
#   theme(axis.text.y = element_blank(),
#         panel.grid.major = element_blank(), 
#         panel.grid.minor = element_blank(),
#         panel.background = element_blank()) 
#   
# ggsave(plot1, file = "verticle_cellpop_plot_singleR.pdf",
#        width = 5, height = 7)
# 
# 
# plot2 <- ggplot(pop2, aes(y = Freq, x =Method,)) + 
#   geom_bar(position="stack", stat="identity",
#            fill = manual_color_list[1:14],
#            width = 0.1) + 
#   xlab("") +
#   ylab("Cell Porportions") +
#   geom_label_repel(aes(label = Var1), colour = "black",
#                    position = position_stack(vjust = 0.5),
#                    max.overlaps = 15 ,
#                    force = 10,
#                    force_pull = 3,
#                    max.iter = 999999,
#                    max.time = 1,
#                    min.segment.length = 0.001,
#                    xlim = c(1, 3)) +
#   theme(axis.text.y = element_blank(),
#         panel.grid.major = element_blank(), 
#         panel.grid.minor = element_blank(),
#         panel.background = element_blank()) 
# 
# 
# ggsave(plot2, file = "verticle_cellpop_plot_manual.pdf",
#        width = 5, height = 7)
# 
#### scCATCH ####
Idents(object = plaqviewobj) <- "seurat_clusters"

clu_markers <- findmarkergenes(
  plaqviewobj,
  species = "Human",
  cluster = 'All',
  match_CellMatch = FALSE, # set T for large dataset
  cancer = NULL,
  tissue = NULL,
  cell_min_pct = 0.25,
  logfc = 0.25,
  pvalue = 0.05
)


## blood vessell ## 
clu_ann_BV <- scCATCH(clu_markers$clu_markers,
                      species = "Human",
                      cancer = NULL,
                      tissue = "Blood vessel")

bv_annotations <- clu_ann_BV$cell_type
names(bv_annotations) <- levels(plaqviewobj)
bv_annotations <- replace_na(bv_annotations, "Unknown")
plaqviewobj[["scCATCH_BV"]] <- bv_annotations[match(plaqviewobj@meta.data$seurat_clusters, names(bv_annotations))]

## heart ##
clu_ann_HT <- scCATCH(clu_markers$clu_markers,
                      species = "Human",
                      cancer = NULL,
                      tissue = "Heart")
# write.csv(clu_ann, file = "scCATCH_vs_singleR_heart.csv")
bv_annotations <- clu_ann_HT$cell_type
names(bv_annotations) <- levels(plaqviewobj)
bv_annotations <- replace_na(bv_annotations, "Unknown")
plaqviewobj[["scCATCH_Heart"]] <- bv_annotations[match(plaqviewobj@meta.data$seurat_clusters, names(bv_annotations))]


## blood ###
clu_ann_Blood <- scCATCH(clu_markers$clu_markers,
                         species = "Human",
                         cancer = NULL,
                         tissue = "Blood")
# write.csv(clu_ann, file = "scCATCH_vs_singleR_blood.csv")

bv_annotations <- clu_ann_Blood$cell_type
names(bv_annotations) <- levels(plaqviewobj)
bv_annotations <- replace_na(bv_annotations, "Unknown")
plaqviewobj[["scCATCH_Blood"]] <- bv_annotations[match(plaqviewobj@meta.data$seurat_clusters, names(bv_annotations))]

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

#### DIFF EX GENE LIST ####
Idents(object = plaqviewobj) <- "SingleR.calls"
difflist <- Seurat::FindAllMarkers(plaqviewobj)
write_csv(difflist, file = "differential/diff_by_singleR.csv")

Idents(object = plaqviewobj) <- "manually_annotated_labels"
difflist <- Seurat::FindAllMarkers(plaqviewobj)
write_csv(difflist, file = "differential/diff_by_author.csv")

Idents(object = plaqviewobj) <- "seurat_clusters"
difflist <- Seurat::FindAllMarkers(plaqviewobj)
write_csv(difflist, file = "differential/diff_by_seurat.csv")

