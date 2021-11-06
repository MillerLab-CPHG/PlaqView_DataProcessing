#### IMPORTANT: SET TO RESPECTIVE DIRECTORY in FUNCTIONS ####
#### Library and data loading ----
library(Seurat)
library(patchwork)
library(readr)
# library(scCATCH)
library(SingleR)
library(tidyverse)
library(monocle3)
library(SeuratData)
library(magrittr)
library(ggrepel)
library(dyno) # devtools::install_github("dynverse/dyno")
library(readxl)
library(SeuratDisk)
library(celldex) # BiocManager::install("celldex")

library(symphony)
library(tidyverse)
library(data.table)
library(matrixStats)
library(Matrix)
library(bayNorm) # for transposition of sparase matrix

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


#### Function for Human Data ####
human_process <- function(datasetID){
  #### STEP1: READ DATASET DIRECTORY ####
  # you must change this if your source is different
  # get this to the dataprocessing - data folder in the first piece
  path.to.destination <- file.path(paste("~/wm5wt@virginia.edu - Google Drive/My Drive/UVA/Grad School/Projects/PlaqView/DataProcessing/data/",
                                         datasetID, "/source_files", sep=""))
  
  setwd(path.to.destination) 
  
  plaqviewobj <- readRDS(file = "UNPROCESSED.rds")
  plaqviewobj <- UpdateSeuratObject(plaqviewobj)
  
  #### STEP1B: READ REFERENCE ####
  humanatlasref <- LoadH5Seurat(file = "~/wm5wt@virginia.edu - Google Drive/My Drive/UVA/Grad School/Projects/PlaqView/DataProcessing/references/Tabula_sapiens_reference/TS_Vasculature.h5seurat", assays = "RNA")
  
  #### STEP2: SEURAT PROCESS ####
  # Run the standard workflow for visualization and clustering
  plaqviewobj <- NormalizeData(plaqviewobj, normalization.method = "LogNormalize", scale.factor = 10000)
  plaqviewobj <- FindVariableFeatures(plaqviewobj, verbose = T, nfeatures = 2000)
  plaqviewobj <- ScaleData(plaqviewobj, verbose = T)
  
  plaqviewobj <- RunPCA(plaqviewobj, npcs = 30, verbose = FALSE)
  plaqviewobj <- RunUMAP(plaqviewobj, reduction = "pca", dims = 1:20)
  plaqviewobj <- FindNeighbors(plaqviewobj, reduction = "pca", dims = 1:20)
  plaqviewobj <- FindClusters(plaqviewobj, resolution = 0.5)
  
  #### STEP3: SINGLER ----
  # BiocManager::install("SingleR")
  # here we are using Human Primary Cell Atlas design for blood
  # https://bioconductor.org/packages/3.12/data/experiment/vignettes/celldex/inst/doc/userguide.html#2_General-purpose_references
  hpca.se <- HumanPrimaryCellAtlasData()
  

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

  #### STEP3A: RECODE SINGLE-R LABELS ----
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


  ## DEPRECATED## #### STEP3B: scCATCH ####
  # Idents(object = plaqviewobj) <- "seurat_clusters"
  # 
  # clu_markers <- findmarkergenes(
  #   plaqviewobj,
  #   species = "Human",
  #   cluster = 'All',
  #   match_CellMatch = FALSE, # set T for large dataset
  #   cancer = NULL,
  #   tissue = NULL,
  #   cell_min_pct = 0.25,
  #   logfc = 0.25,
  #   pvalue = 0.05
  # )
  # 
  # 
  # ## blood vessell ##
  # clu_ann_BV <- scCATCH(clu_markers$clu_markers,
  #                       species = "Human",
  #                       cancer = NULL,
  #                       tissue = "Blood vessel")
  # 
  # bv_annotations <- clu_ann_BV$cell_type
  # names(bv_annotations) <- levels(plaqviewobj)
  # bv_annotations <- replace_na(bv_annotations, "Unknown")
  # plaqviewobj[["scCATCH_BV"]] <- bv_annotations[match(plaqviewobj@meta.data$seurat_clusters, names(bv_annotations))]
  # 
  # ## heart ##
  # clu_ann_HT <- scCATCH(clu_markers$clu_markers,
  #                       species = "Human",
  #                       cancer = NULL,
  #                       tissue = "Heart")
  # # write.csv(clu_ann, file = "scCATCH_vs_singleR_heart.csv")
  # bv_annotations <- clu_ann_HT$cell_type
  # names(bv_annotations) <- levels(plaqviewobj)
  # bv_annotations <- replace_na(bv_annotations, "Unknown")
  # plaqviewobj[["scCATCH_Heart"]] <- bv_annotations[match(plaqviewobj@meta.data$seurat_clusters, names(bv_annotations))]
  # 
  # 
  # ## blood ###
  # clu_ann_Blood <- scCATCH(clu_markers$clu_markers,
  #                          species = "Human",
  #                          cancer = NULL,
  #                          tissue = "Blood")
  # # write.csv(clu_ann, file = "scCATCH_vs_singleR_blood.csv")
  # 
  # bv_annotations <- clu_ann_Blood$cell_type
  # names(bv_annotations) <- levels(plaqviewobj)
  # bv_annotations <- replace_na(bv_annotations, "Unknown")
  # plaqviewobj[["scCATCH_Blood"]] <- bv_annotations[match(plaqviewobj@meta.data$seurat_clusters, names(bv_annotations))]
  # 
  # #### STEP3B: SYMPHONY ####
  # ref_pbmcs = readRDS('references/Symphony_ref_data/fibroblast_atlas.rds')
  # 
  # query = symphony::mapQuery(plaqviewobj@assays$RNA, plaqviewobj@meta.data, ref_pbmcs, 
  #                  vars = plaqviewobj@meta.data, 
  #                  do_normalize = TRUE)
  # 
  #### STEP3C: SEURAT/TABULA SAPIENS LABELING ####

  #### preprocess ref seurat 
  Idents(humanatlasref) <-  humanatlasref@meta.data[["Annotation"]]
  
  humanatlasref <- NormalizeData(humanatlasref, verbose = T)
  humanatlasref <- FindVariableFeatures(humanatlasref, selection.method = "vst", verbose = T)
  
  DefaultAssay(plaqviewobj) <- 'RNA'
  DefaultAssay(humanatlasref) <- 'RNA'
  
  
  anchors <- FindTransferAnchors(reference = humanatlasref, query = plaqviewobj, 
                                 dims = 1:30)
  predictions <- TransferData(anchorset = anchors, refdata = humanatlasref$Annotation, 
                              dims = 1:30)
  plaqviewobj <- AddMetaData(plaqviewobj, metadata = predictions)
  
  #### rename transferred column metadata 
  plaqviewobj@meta.data[["Seurat_with_Tabula_Ref"]] <- plaqviewobj@meta.data[["predicted.id"]]
  
  # capitalize the lettering
  plaqviewobj@meta.data[["Seurat_with_Tabula_Ref"]] <-str_to_title(plaqviewobj@meta.data[["Seurat_with_Tabula_Ref"]], locale = "en")
  
  # set to active idents
  Idents(plaqviewobj) <- plaqviewobj@meta.data[["Seurat_with_Tabula_Ref"]]
  
  #### STEP3D: RECODE SEURAT/TABULA SAPIENS ####
  plaqviewobj@meta.data[["Seurat_with_Tabula_Ref"]] <- recode(plaqviewobj@meta.data[["Seurat_with_Tabula_Ref"]], 
                                                                   'Smooth Muscle Cell' = "SMCs")
  plaqviewobj@meta.data[["Seurat_with_Tabula_Ref"]] <- recode(plaqviewobj@meta.data[["Seurat_with_Tabula_Ref"]], 
                                                                   'Pancreatic Acinar Cell' = "Panc Acinar Cell")
  plaqviewobj@meta.data[["Seurat_with_Tabula_Ref"]] <- recode(plaqviewobj@meta.data[["Seurat_with_Tabula_Ref"]], 
                                                                   'Fibroblast' = "FB")
  plaqviewobj@meta.data[["Seurat_with_Tabula_Ref"]] <- recode(plaqviewobj@meta.data[["Seurat_with_Tabula_Ref"]], 
                                                                   'Endothelial Cell' = "EC")
  plaqviewobj@meta.data[["Seurat_with_Tabula_Ref"]] <- recode(plaqviewobj@meta.data[["Seurat_with_Tabula_Ref"]], 
                                                                   'Macrophage' = "Mø")
  plaqviewobj@meta.data[["Seurat_with_Tabula_Ref"]] <- recode(plaqviewobj@meta.data[["Seurat_with_Tabula_Ref"]], 
                                                                   'Natural Killer Cell' = "NK")
  Idents(plaqviewobj) <- plaqviewobj@meta.data[["Seurat_with_Tabula_Ref"]]
  
  #### plot the cells 
  DimPlot(
    plaqviewobj,
    reduction = "umap",
    label = TRUE,
    label.size = 5,
    repel = T,
    # repel labels
    pt.size = 1,
    cols = manual_color_list) + # group.by is important, use this to call metadata separation
    theme(legend.position="bottom", 
          legend.box = "vertical") +
    ggtitle("UMAP by Cell Type") +
    theme(plot.title = element_text(hjust =  0.5)) +
    guides(color = guide_legend(nrow = 5))
  
  #### STEP3E: RESERVED SPACE #### 
  #### STEP4: MONOCLE3 TRAJECTORY INFERENCE ----
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
  
  
  #### STEP4A: TRANSFER SEURAT EMBEDDINGS ###
  # Note that these may be calculated on the Integrated object, not the counts
  #   and thus will involve fewer genes
  temp.cds <- ProjectDim(plaqviewobj, reduction = "pca") # this will be removed
  reducedDim(plaqviewobj.cds, type = "PCA") <- temp.cds@reductions$pca@cell.embeddings
  plaqviewobj.cds@preprocess_aux$prop_var_expl <- temp.cds@reductions$pca@stdev
  plot_pc_variance_explained(plaqviewobj.cds)
  
  # Transfer Seurat UMAP embeddings
  plaqviewobj.cds@int_colData@listData$reducedDims$UMAP <- temp.cds@reductions$umap@cell.embeddings
  
  ## transfer singleR labels to moncle3 object
  colData(plaqviewobj.cds)$assigned_cell_type <- 
    plaqviewobj@meta.data[["Seurat_with_Tabula_Ref"]] # call this by opening the object
  
  #### MONOCLE3 CONT. ---
  # now learn the PATH (trajectory)
  plaqviewobj.cds <- learn_graph(plaqviewobj.cds)
  
  # a helper function to identify the root principal points:
  get_earliest_principal_node <- function(cds, assigned_cell_type= "SMCs"){ # change celltype if desired
    cell_ids <- which(colData(cds)[, "assigned_cell_type"] == assigned_cell_type)
    
    closest_vertex <-
      cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
    closest_vertex <- as.matrix(closest_vertex[colnames(cds), ])
    root_pr_nodes <-
      igraph::V(principal_graph(cds)[["UMAP"]])$name[as.numeric(names
                                                                (which.max(table(closest_vertex[cell_ids,]))))]
    
    root_pr_nodes
  }
  
  plaqviewobj.cds <- order_cells(plaqviewobj.cds, root_pr_nodes=get_earliest_principal_node(plaqviewobj.cds),
                                 reduction_method = "UMAP")
  
  # 
  # # this calls up a shiny app, choose the ROOT NODE
  # plaqviewobj.cds <- order_cells(plaqviewobj.cds, reduction_method = "UMAP")
  
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
  
  saveRDS(mon3, file = "../monocle3.rds")
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
  
  
  #### STEP5A: DYNO TRAJECTORY INFERENCES ####
  object_counts <- Matrix::t(Matrix(plaqviewobj@assays$RNA@counts, sparse = T))
  object_expression <- Matrix::t(Matrix(plaqviewobj@assays$RNA@data, sparse = T))
  object_cellinfo <- plaqviewobj@meta.data[["Seurat_with_Tabula_Ref"]]
  
  plaqviewobj.dyno <- wrap_expression(
    counts = object_counts,
    expression = object_expression)
  
  
  # may need to run this to clear out memory
  # rm(list= ls()[!(ls() %in% c('plaqviewobj.dyno','plaqviewobj','object_cellinfo'))])
  # gc()
  
  
  
  # may need to run this to clear out memory
  # rm(list= ls()[!(ls() %in% c('plaqviewobj.dyno','plaqviewobj','object_cellinfo'))])
  # gc()
  
  #### STEP5B: SlingShot ####
  # make sure to call up docker images

  model <- infer_trajectory(plaqviewobj.dyno, "slingshot", verbose = T,
                            cluster_method = 'clara') # THIS IS ADDED TO REDUCE MEM REQ

  # add dim reduction
  model <- model %>%
    add_dimred(dimred = as.matrix(plaqviewobj@reductions$umap@cell.embeddings),
               expression_source = plaqviewobj.dyno$expression)

  pdf("../dyno_slingshot_full.pdf", width=7, height=6)
  slingshot <- plot_dimred(
    model,
    expression_source = plaqviewobj.dyno$expression,
    grouping = object_cellinfo # basically stanford@meta.data[["SingleR.calls"]]
  )

  saveRDS(slingshot, file = "../slingshot.rds")
  slingshot
  dev.off()

  #### STEP5C: scorpius ####
  # make sure to call up docker images

  model <- infer_trajectory(plaqviewobj.dyno, "scorpius", verbose = T,
                            cluster_method = 'clara') # THIS IS ADDED TO REDUCE MEM REQ

  # add dim reduction
  model <- model %>%
    add_dimred(dimred = as.matrix(plaqviewobj@reductions$umap@cell.embeddings),
               expression_source = plaqviewobj.dyno$expression)

  pdf("../dyno_scorpius_full.pdf", width=7, height=6)
  scorpius <- plot_dimred(
    model,
    expression_source = plaqviewobj.dyno$expression,
    grouping = object_cellinfo # basically plaqviewobj@meta.data[["SingleR.calls"]]
  )

  saveRDS(scorpius, file = "../scorpius.rds")
  scorpius
  dev.off()

  #### STEP5D: PAGA ####
  model <- infer_trajectory(plaqviewobj.dyno, "projected_paga", verbose = T,
                            cluster_method = 'clara') # THIS IS ADDED TO REDUCE MEM REQ

  #### PAGA: project the model ###
  # add dim reduction
  model <- model %>%
    add_dimred(dimred = as.matrix(plaqviewobj@reductions$umap@cell.embeddings),
               expression_source = plaqviewobj.dyno$expression)
  
  pdf("../dyno_paga_full.pdf", width=7, height=6)
  paga <- plot_dimred(
    model,
    expression_source = plaqviewobj.dyno$expression,
    grouping = object_cellinfo # basically stanford@meta.data[["SingleR.calls"]]
  )
  paga
  saveRDS(paga, file = "../paga.rds")
  dev.off()
  
  #### Clean-Up Metadata ####
  # show all metadata columns
  names(plaqviewobj@meta.data)
  
  # rename number clusters  
  plaqviewobj@meta.data$Seurat_Clusters <- plaqviewobj@meta.data$seurat_clusters

  # choose which ones to keep for display
  plaqviewobj@meta.data <- 
    plaqviewobj@meta.data[, which(colnames(plaqviewobj@meta.data)
                                  %in% c(
                                    "Seurat_Clusters",
                                    "Author_Provided",
                                    "SingleR.calls",
                                    "Seurat_with_Tabula_Ref"  
                                  ))]
  #### STEP6: REDUCE SIZE & OUTPUT ####
  plaqviewobj <- DietSeurat(plaqviewobj, counts = T, data = T, dimreducs = c('umap'))
  
  final.file.name <- file.path(paste("../", datasetID, ".rds", sep="")) # ../ moves up one level in file
  
  saveRDS(plaqviewobj, file = final.file.name)
  
  plaqviewobj <- readRDS(file = final.file.name)
  
  #### STEP7: DIFF EX GENE LIST ####
  Idents(object = plaqviewobj) <- "SingleR.calls"
  difflist <- Seurat::FindAllMarkers(plaqviewobj)
  write_csv(difflist, file = "../diff_by_singleR.csv")
  
  Idents(object = plaqviewobj) <- "Author_Provided"
  difflist <- Seurat::FindAllMarkers(plaqviewobj)
  write_csv(difflist, file = "../diff_by_author.csv")
  
  Idents(object = plaqviewobj) <- "Seurat_Clusters"
  difflist <- Seurat::FindAllMarkers(plaqviewobj)
  write_csv(difflist, file = "../diff_by_seurat.csv")
  
  Idents(object = plaqviewobj) <- "Seurat_with_Tabula_Ref"
  difflist <- Seurat::FindAllMarkers(plaqviewobj)
  write_csv(difflist, file = "../diff_by_Seurat_with_Tabula_Ref.csv")
  

  
  #### STEP8: PRINT CELL COUNT, RAM gc ####
  
  n <- summary(plaqviewobj$SingleR.calls)[1]
  tab <- summary(as.factor(plaqviewobj$SingleR.calls))
  
  write.csv(tab, file = paste(n, "_cell_count.csv"))
  
  gc()

  
}
#### Function for Mouse Data ####
mouse_process <- function(datasetID){
  path.to.unprocessed <- file.path(paste("data/", datasetID, "/source_files/",
                                         "UNPROCESSED.rds", sep=""))
  path.to.destination <- file.path(paste("data/", datasetID, "/source_files/", sep=""))
  setwd(path.to.destination) 
  plaqviewobj <- readRDS(file = path.to.unprocessed)
  plaqviewobj <- UpdateSeuratObject(plaqviewobj)
  
  #### STEP2: SEURAT PROCESS ####
  # Run the standard workflow for visualization and clustering
  plaqviewobj <- ScaleData(plaqviewobj, verbose = T)
  plaqviewobj <- FindVariableFeatures(plaqviewobj, verbose = T)
  plaqviewobj <- RunPCA(plaqviewobj, npcs = 30, verbose = T)
  plaqviewobj <- RunUMAP(plaqviewobj, reduction = "pca", dims = 1:30)
  plaqviewobj <- FindNeighbors(plaqviewobj, reduction = "pca", dims = 1:30)
  plaqviewobj <- FindClusters(plaqviewobj, resolution = 0.5)
  
  # #### STEP3: SINGLER ----
  # # BiocManager::install("SingleR")
  # # https://bioconductor.org/packages/3.12/data/experiment/vignettes/celldex/inst/doc/userguide.html#2_General-purpose_references
  # hpca.se <- celldex::MouseRNAseqData()
  # hpca.se
  # 
  # # now run the prediction using the reference
  # # singleR requires that it be in a 'singlecellexperiment' format
  # # they are workout agnostic
  # 
  # for_singleR_input <- GetAssayData(plaqviewobj)
  # pred.plaqviewobj <- SingleR(test = for_singleR_input, 
  #                             ref = hpca.se, 
  #                             label = hpca.se$label.main) # reference cell types
  # pred.plaqviewobj
  # # summarize distribution
  # table(pred.plaqviewobj$labels)
  # 
  # # to show annotation confidence map
  # plotScoreHeatmap(pred.plaqviewobj)
  # 
  # # to show # that are pruned due to low score
  # summary(is.na(pred.plaqviewobj$pruned.labels))
  # 
  # ### to place the singleR predictions into Seurat as a sep unit ###
  # # seurat.obj[["SingleR.labels"]] <- singler.results$labels
  # plaqviewobj[["SingleR.labels"]] <- pred.plaqviewobj$labels # this nest under metadata
  # 
  # # Copy over the labels and pruned.labels (Note: any other column of the results could be used as well)
  # plaqviewobj$SingleR.pruned.calls <- pred.plaqviewobj$pruned.labels
  # plaqviewobj$SingleR.calls <- pred.plaqviewobj$labels
  # 
  # #### STEP3A: RECODE SINGLE-R LABELS ----
  # plaqviewobj@meta.data[["SingleR.calls"]] <- recode(plaqviewobj@meta.data[["SingleR.calls"]], Smooth_muscle_cells = "SMC")
  # plaqviewobj@meta.data[["SingleR.calls"]] <- recode(plaqviewobj@meta.data[["SingleR.calls"]], Endothelial_cells = "EC")
  # plaqviewobj@meta.data[["SingleR.calls"]] <- recode(plaqviewobj@meta.data[["SingleR.calls"]], NK_cell = "NK")
  # plaqviewobj@meta.data[["SingleR.calls"]] <- recode(plaqviewobj@meta.data[["SingleR.calls"]], Chondrocytes = "CH")
  # plaqviewobj@meta.data[["SingleR.calls"]] <- recode(plaqviewobj@meta.data[["SingleR.calls"]], Fibroblasts = "FB")
  # plaqviewobj@meta.data[["SingleR.calls"]] <- recode(plaqviewobj@meta.data[["SingleR.calls"]], Monocyte = "Mono")
  # plaqviewobj@meta.data[["SingleR.calls"]] <- recode(plaqviewobj@meta.data[["SingleR.calls"]], B_cell = "B_Cells")
  # plaqviewobj@meta.data[["SingleR.calls"]] <- recode(plaqviewobj@meta.data[["SingleR.calls"]], Macrophage = "Mø")
  # plaqviewobj@meta.data[["SingleR.calls"]] <- recode(plaqviewobj@meta.data[["SingleR.calls"]], Tissue_stem_cells = "SC")
  # plaqviewobj@meta.data[["SingleR.calls"]] <- recode(plaqviewobj@meta.data[["SingleR.calls"]], T_cells = "T_Cells")
  # plaqviewobj@meta.data[["SingleR.calls"]] <- recode(plaqviewobj@meta.data[["SingleR.calls"]], 'Pre-B_cell_CD34-' = "PreB_CD34-")
  # plaqviewobj@meta.data[["SingleR.calls"]] <- recode(plaqviewobj@meta.data[["SingleR.calls"]], 'Pro-B_cell_CD34+' = "ProB_CD34+")
  # 
  # 
  # table(plaqviewobj@meta.data[["SingleR.calls"]])
  # 
  
  # #### STEP3B: scCATCH ####
  # Idents(object = plaqviewobj) <- "seurat_clusters"
  # 
  # clu_markers <- findmarkergenes(
  #   plaqviewobj,
  #   species = 'Mouse',
  #   cluster = 'All',
  #   match_CellMatch = TRUE, # set T for large dataset
  #   cancer = NULL,
  #   tissue = 'Aorta',
  #   cell_min_pct = 0.25,
  #   logfc = 0.25,
  #   pvalue = 0.05
  # )
  # 
  # 
  # ## blood vessell ## 
  # clu_ann_BV <- scCATCH(clu_markers$clu_markers,
  #                       species = "Mouse",
  #                       cancer = NULL,
  #                       tissue = "Blood vessel")
  # 
  # bv_annotations <- clu_ann_BV$cell_type
  # names(bv_annotations) <- levels(plaqviewobj)
  # bv_annotations <- replace_na(bv_annotations, "Unknown")
  # plaqviewobj[["scCATCH_BV"]] <- bv_annotations[match(plaqviewobj@meta.data$seurat_clusters, names(bv_annotations))]
  # 
  # ## heart ##
  # clu_ann_HT <- scCATCH(clu_markers$clu_markers,
  #                       species = "Mouse",
  #                       cancer = NULL,
  #                       tissue = "Heart")
  # # write.csv(clu_ann, file = "scCATCH_vs_singleR_heart.csv")
  # bv_annotations <- clu_ann_HT$cell_type
  # names(bv_annotations) <- levels(plaqviewobj)
  # bv_annotations <- replace_na(bv_annotations, "Unknown")
  # plaqviewobj[["scCATCH_Heart"]] <- bv_annotations[match(plaqviewobj@meta.data$seurat_clusters, names(bv_annotations))]
  # 
  # 
  # ## blood ###
  # clu_ann_Blood <- scCATCH(clu_markers$clu_markers,
  #                          species = "Mouse",
  #                          cancer = NULL,
  #                          tissue = "Blood")
  # # write.csv(clu_ann, file = "scCATCH_vs_singleR_blood.csv")
  # 
  # bv_annotations <- clu_ann_Blood$cell_type
  # names(bv_annotations) <- levels(plaqviewobj)
  # bv_annotations <- replace_na(bv_annotations, "Unknown")
  # plaqviewobj[["scCATCH_Blood"]] <- bv_annotations[match(plaqviewobj@meta.data$seurat_clusters, names(bv_annotations))]
  # 
  # #### STEP3C: rename manual label####
  # 
  # plaqviewobj[["manually_annotated_labels"]] <- plaqviewobj$orig.ident
  
  #### STEP3D: SEURAT/TABULA MURINE LABELING ####
  #### load tabulus reference
  # using only the aorta and heart data
  mouseatlas <- read_rds(file = "~/Google Drive/UVA/Grad School/Projects/PlaqView/DataProcessing/Tabula_muris_reference/TM_aorta_heart.rds")
  mouseatlas <- UpdateSeuratObject(mouseatlas)
  
  #### preprocess ref seurat 
  Idents(mouseatlas) <-  mouseatlas@meta.data[["cell_ontology_class"]]
  DefaultAssay(mouseatlas) <- "RNA" 
  DefaultAssay(plaqviewobj) <- "RNA"
  
  mouseatlas <- NormalizeData(mouseatlas, verbose = T)
  mouseatlas <- FindVariableFeatures(mouseatlas, selection.method = "vst", verbose = T)
  ElbowPlot(plaqviewobj, ndims = 20, reduction = "pca")
  ElbowPlot(mouseatlas, ndims = 20, reduction = "pca")
  
  anchors <- FindTransferAnchors(reference = mouseatlas, query = plaqviewobj, 
                                 dims = 1:10,
                                 approx = FALSE)
  
  predictions <- TransferData(anchorset = anchors, refdata = mouseatlas@meta.data[["cell_ontology_class"]], 
                              dims = 1:10)
  plaqviewobj <- AddMetaData(plaqviewobj, metadata = predictions)
  
  #### rename transferred column metadata 
  plaqviewobj@meta.data[["predicted.id_Tabula"]] <- plaqviewobj@meta.data[["predicted.id"]]
  
  # capitalize the lettering
  plaqviewobj@meta.data[["predicted.id_Tabula"]] <-str_to_title(plaqviewobj@meta.data[["predicted.id_Tabula"]], locale = "en")
  
  # set to active idents
  Idents(plaqviewobj) <- plaqviewobj@meta.data[["predicted.id_Tabula"]]
  
  #### STEP3E: recode tabula muris labels ####
  plaqviewobj@meta.data[["predicted.id_Tabula"]] <- recode(plaqviewobj@meta.data[["predicted.id_Tabula"]], 
                                                           'Smooth Muscle Cell' = "SMCs")
  plaqviewobj@meta.data[["predicted.id_Tabula"]] <- recode(plaqviewobj@meta.data[["predicted.id_Tabula"]], 
                                                           'Pancreatic Acinar Cell' = "Panc Acinar Cell")
  plaqviewobj@meta.data[["predicted.id_Tabula"]] <- recode(plaqviewobj@meta.data[["predicted.id_Tabula"]], 
                                                           'Fibroblast' = "FB")
  plaqviewobj@meta.data[["predicted.id_Tabula"]] <- recode(plaqviewobj@meta.data[["predicted.id_Tabula"]], 
                                                           'Endothelial Cell' = "EC")
  plaqviewobj@meta.data[["predicted.id_Tabula"]] <- recode(plaqviewobj@meta.data[["predicted.id_Tabula"]], 
                                                           'Macrophage' = "Mø")
  plaqviewobj@meta.data[["predicted.id_Tabula"]] <- recode(plaqviewobj@meta.data[["predicted.id_Tabula"]], 
                                                           'Natural Killer Cell' = "NK")
  Idents(plaqviewobj) <- plaqviewobj@meta.data[["predicted.id_Tabula"]]
  
  #### plot the cells 
  DimPlot(
    plaqviewobj,
    reduction = "umap",
    label = TRUE,
    label.size = 5,
    repel = T,
    # repel labels
    pt.size = 1,
    cols = manual_color_list) + # group.by is important, use this to call metadata separation
    theme(legend.position="bottom", 
          legend.box = "vertical") +
    ggtitle("UMAP by Cell Type") +
    theme(plot.title = element_text(hjust =  0.5)) +
    guides(color = guide_legend(nrow = 5))
  
  #### STEP4: MONOCLE3 TRAJECTORY INFERENCE ----
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
  
  
  # run clustering again (didnt transfer from seurat)
  plaqviewobj.cds <- reduce_dimension(plaqviewobj.cds, reduction_method = "UMAP")
  plaqviewobj.cds <- cluster_cells(plaqviewobj.cds, reduction_method = "UMAP")
  
  
  #### STEP4A: TRANSFER SEURAT EMBEDDINGS ###
  # Note that these may be calculated on the Integrated object, not the counts
  #   and thus will involve fewer genes
  temp.cds <- ProjectDim(plaqviewobj, reduction = "pca") # this will be removed
  reducedDim(plaqviewobj.cds, type = "PCA") <- temp.cds@reductions$pca@cell.embeddings
  plaqviewobj.cds@preprocess_aux$prop_var_expl <- temp.cds@reductions$pca@stdev
  plot_pc_variance_explained(plaqviewobj.cds)
  
  # Transfer Seurat UMAP embeddings
  plaqviewobj.cds@int_colData@listData$reducedDims$UMAP <- temp.cds@reductions$umap@cell.embeddings
  
  ## transfer singleR labels to moncle3 object
  colData(plaqviewobj.cds)$assigned_cell_type <- 
    plaqviewobj$predicted.id_Tabula # call this by opening the object
  
  #### MONOCLE3 CONT. ---
  # now learn the PATH (trajectory)
  plaqviewobj.cds <- learn_graph(plaqviewobj.cds)
  
  # a helper function to identify the root principal points:
  get_earliest_principal_node <- function(cds, assigned_cell_type= "Endocardial Cell"){ # change celltype if desired
    cell_ids <- which(colData(cds)[, "assigned_cell_type"] == assigned_cell_type)
    
    closest_vertex <-
      cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
    closest_vertex <- as.matrix(closest_vertex[colnames(cds), ])
    root_pr_nodes <-
      igraph::V(principal_graph(cds)[["UMAP"]])$name[as.numeric(names
                                                                (which.max(table(closest_vertex[cell_ids,]))))]
    
    root_pr_nodes
  }
  
  plaqviewobj.cds <- order_cells(plaqviewobj.cds, root_pr_nodes=get_earliest_principal_node(plaqviewobj.cds),
                                 reduction_method = "UMAP")
  
  # 
  # # this calls up a shiny app, choose the ROOT NODE
  # plaqviewobj.cds <- order_cells(plaqviewobj.cds, reduction_method = "UMAP")
  
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
  
  saveRDS(mon3, file = "../monocle3.rds")
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
  
  
  # #### STEP4B: (NOT RUN) Subset Trajectory & analysis of SMC----
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
  # # #### STEP4C: (NOT RUN) MORAN's I Test of Autocorrelation ####
  # # now we can extrapolate genes that are differentially expressed in this region
  # # Moran’s I is a measure of multi-directional and multi-dimensional spatial autocorrelation.
  # # the statistic tells you whether cells at nearby positions on a
  # # trajectory will have similar (or dissimilar) +
  # # expression levels for the gene being tested.
  # ## first lets do the whole dataset
  # # a special gene module score heatmap (for the whole dataset)
  # # pr_graph_test_res <- graph_test(plaqviewobj.cds, neighbor_graph="principal_graph", cores=2)
  # write.csv(pr_graph_test_res, file = "../moransI_all_clusters.csv")
  # pr_deg_ids <- row.names(subset(pr_graph_test_res, q_value < 0.00000001)) # you can adjust the p-value here
  # head(pr_deg_ids)
  # gene_module_df <- find_gene_modules(plaqviewobj.cds[pr_deg_ids,], resolution=1e-3)
  # cell_group_df <- tibble::tibble(cell=row.names(colData(plaqviewobj.cds)),
  #                                 cell_group=colData(plaqviewobj.cds)$assigned_cell_type)
  # agg_mat <- aggregate_gene_expression(plaqviewobj.cds, gene_module_df, cell_group_df)
  # row.names(agg_mat) <- stringr::str_c("Module ", row.names(agg_mat))
  # pheatmap::pheatmap(agg_mat,
  #                    scale="column", clustering_method="ward.D2")
  # 
  # # which then can be visualized like so;
  # # this can show you the different gene modules that can are responsible for changes over pseudotime
  # plot_cells(plaqviewobj.cds,
  #            genes=gene_module_df %>% filter(module %in% c(2,3,7)), # specify the module you want to examine
  #            label_cell_groups=T,
  #            show_trajectory_graph=F)
  # 
  # subset(gene_module_df, module == 2)
  # 
  # ## now lets do the subsets
  # # pr_graph_test_res.sub <- graph_test(plaqviewobj.cds_subset, neighbor_graph="principal_graph", cores=2)
  # pr_deg_ids.sub <- row.names(subset(pr_graph_test_res.sub, q_value < 0.00000001))
  # write.csv(pr_graph_test_res.sub, file = "moransI_subset_cluster.csv")
  # head(pr_deg_ids.sub)
  # 
  # # collect the trajectory-variable genes into modules
  # gene_module_df.sub <- find_gene_modules(plaqviewobj.cds_subset[pr_deg_ids.sub,], resolution=1e-3)
  # # visualize these genes
  # # here I am just pulling out genes that have high moran's i and might be helpful in the paper
  # # SELECTED FOR PUBLICATIONS
  # pdf("monocle3_genesoverpseudotime_seuratpartition_extended.pdf", width=7, height=6)
  # plot_cells(plaqviewobj.cds_subset,
  #            genes=c("MYH11", 'IGFBP2',"PPP1R14A","CNN1", "TNFRSF11B",
  #                    "C7", "C3",
  #                    "SERPINF1",  "FBLN1",
  #                    "CXCL12", "MMP2",
  #                    "FN1"), # this is faceting by the genes that are DE
  #            show_trajectory_graph=FALSE,
  #            label_cell_groups=F, cell_size = 1)
  # 
  # dev.off()
  # 
  # # recluster at higher definition
  # plaqviewobj.cds_subset = cluster_cells(plaqviewobj.cds_subset, resolution=1e-2)
  # 
  # pdf("monocle3_RNAvelocitySUBSET_seuratpartition.pdf", width=6, height=6)
  # plot_cells(plaqviewobj.cds_subset,
  #            color_cells_by="cluster",
  #            label_groups_by_cluster=F,
  #            show_trajectory_graph = T,
  #            trajectory_graph_segment_size = 1,
  #            label_leaves=F, # this gives a little node label (outcome)
  #            label_roots = F,
  #            label_branch_points = F,
  #            graph_label_size = 1, # size of # in circle
  #            group_label_size = 4,
  #            cell_size = 1,
  #            alpha = 0.5,
  #            scale_to_range = T)
  # dev.off()
  # 
  #### STEP5A: DYNO TRAJECTORY INFERENCES ####
  object_counts <- Matrix::t(as(as.matrix(plaqviewobj@assays$RNA@counts), 'sparseMatrix'))
  object_expression <- Matrix::t(as(as.matrix(plaqviewobj@assays$RNA@data), 'sparseMatrix'))
  object_cellinfo <- plaqviewobj@meta.data[["SingleR.labels"]]
  
  plaqviewobj.dyno <- wrap_expression(
    counts = object_counts,
    expression = object_expression)
  
  
  # may need to run this to clear out memory
  # rm(list= ls()[!(ls() %in% c('plaqviewobj.dyno','plaqviewobj','object_cellinfo'))])
  # gc()
  
  #### STEP5B: SlingShot ####
  # make sure to call up docker images
  
  model <- infer_trajectory(plaqviewobj.dyno, "slingshot", verbose = T,
                            cluster_method = 'clara') # THIS IS ADDED TO REDUCE MEM REQ
  
  # add dim reduction
  model <- model %>%
    add_dimred(dimred = as.matrix(plaqviewobj@reductions$umap@cell.embeddings),
               expression_source = plaqviewobj.dyno$expression)
  
  pdf("dyno_slingshot_full.pdf", width=7, height=6)
  slingshot <- plot_dimred(
    model,
    expression_source = plaqviewobj.dyno$expression,
    grouping = object_cellinfo # basically stanford@meta.data[["SingleR.labels"]]
  )
  
  saveRDS(slingshot, file = "../slingshot.rds")
  slingshot
  dev.off()
  
  #### STEP5C: scorpius ####
  # make sure to call up docker images
  
  model <- infer_trajectory(plaqviewobj.dyno, "scorpius", verbose = T,
                            cluster_method = 'clara') # THIS IS ADDED TO REDUCE MEM REQ
  
  # add dim reduction
  model <- model %>%
    add_dimred(dimred = as.matrix(plaqviewobj@reductions$umap@cell.embeddings),
               expression_source = plaqviewobj.dyno$expression)
  
  pdf("dyno_scorpius_full.pdf", width=7, height=6)
  scorpius <- plot_dimred(
    model,
    expression_source = plaqviewobj.dyno$expression,
    grouping = object_cellinfo # basically plaqviewobj@meta.data[["SingleR.labels"]]
  )
  
  saveRDS(scorpius, file = "../scorpius.rds")
  scorpius
  dev.off()
  
  #### STEP5D: PAGA ####
  model <- infer_trajectory(plaqviewobj.dyno, "projected_paga", verbose = T,
                            cluster_method = 'clara') # THIS IS ADDED TO REDUCE MEM REQ
  
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
  saveRDS(paga, file = "../paga.rds")
  
  
  #### Clean-Up Metadata ####
  
  # show all metadata columns
  names(plaqviewobj@meta.data)
  
  # last chance to rename any 
  plaqviewobj@meta.data$Seurat_Clusters <- plaqviewobj@meta.data$seurat_clusters

  # choose which ones to keep for display
  plaqviewobj@meta.data <- 
    plaqviewobj@meta.data[, which(colnames(plaqviewobj@meta.data)
                                  %in% c(
                                    "Seurat_Clusters",
                                    "scCATCH_Blood",
                                    "scCATCH_BV",
                                    "scCATCH_Heart",
                                    "Author_Provided",
                                    "SingleR.calls",
                                    "Seurat_with_Tabula_Ref"  
                                  ))]
  
  #### STEP6: REDUCE SIZE & OUTPUT ####
  plaqviewobj <- DietSeurat(plaqviewobj, counts = T, data = T, dimreducs = c('umap'))
  
  saveRDS(plaqviewobj, file = "../PROCESSED_NEEDTORENAME.rds")
  # saveRDS(plaqviewobj.cds, file = "../PROCESSED_CDS_NO_FOR_PLAQVIEW_NEEDTORENAME.rds")
  
  plaqviewobj <- readRDS(file = "../PROCESSED_NEEDTORENAME.rds")
  
  #### STEP7: DIFF EX GENE LIST ####
  Idents(object = plaqviewobj) <- "SingleR.calls"
  difflist <- Seurat::FindAllMarkers(plaqviewobj)
  write_csv(difflist, file = "../diff_by_singleR.csv")
  
  Idents(object = plaqviewobj) <- "Author_Provided"
  difflist <- Seurat::FindAllMarkers(plaqviewobj)
  write_csv(difflist, file = "../diff_by_author.csv")
  
  Idents(object = plaqviewobj) <- "Seurat_Clusters"
  difflist <- Seurat::FindAllMarkers(plaqviewobj)
  write_csv(difflist, file = "../diff_by_seurat.csv")
  
  Idents(object = plaqviewobj) <- "Seurat_with_Tabula_Ref"
  difflist <- Seurat::FindAllMarkers(plaqviewobj)
  write_csv(difflist, file = "../diff_by_predicted.id_Tabula.csv")
  
  
  ####
  
  
  
  
  
}
#### Read Current Available Data Sets ####
# this needs to point to the set spreadsheet
all.data <- 
  read_excel("data/Available_datasets.xlsx")

#### Pull Out Human and Mouse Data IDs ####
humanIDs <- all.data$DataID[all.data$Species == "Human"]
mouseIDs <- all.data$DataID[all.data$Species == "Mouse"]

#### Apply Functions for Human Sets ####
# lapply(humanIDs, human_process)

human_process(datasetID = "Li_2020")
human_process(datasetID = "Wirka_2019")

#### Apply Functions for Mouse Sets ####
# tapply: Apply a function to subsets of a vector X and defined the subset by vector Y.
