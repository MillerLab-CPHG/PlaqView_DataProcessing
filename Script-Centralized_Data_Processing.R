#### IMPORTANT: SET TO RESPECTIVE DIRECTORY in FUNCTIONS ####

# point this FROM THE ROOT "./" to DataProcessing's 'data' folder 
root.to.data <- "~/Documents/My Drive/PlaqView_Master/DataProcessing/data/"
root.to.TS.ref <- "~/Documents/My Drive/PlaqView_Master/DataProcessing/references/Tabula_sapiens_reference/TS_Vasculature.h5seurat"
root.to.TS.mouse.ref <- "~/Documents/My Drive/PlaqView_Master/DataProcessing/references/Tabula_muris_reference/updated.TS.muris.RDS"
root.to.mastermetadata <- "~/Documents/My Drive/PlaqView_Master/DataProcessing/Summary-Master_Metadata.csv"

# # this version is for weis macbookpro
# root.to.data <- "~/Documents/My Drive/PlaqView_Master/DataProcessing/data/"
# root.to.TS.ref <- "~/Documents/My Drive/PlaqView_Master/DataProcessing/references/Tabula_sapiens_reference/TS_Vasculature.h5seurat"
# root.to.TS.mouse.ref <- "~/Documents/My Drive/PlaqView_Master/DataProcessing/references/Tabula_muris_reference/updated.TS.muris.RDS"

# this version if outside of docker
# root.to.data <- "~/My Drive (wm5wt@virginia.edu)/PlaqView_Master/DataProcessing/data/"
# root.to.TS.ref <- "~/My Drive (wm5wt@virginia.edu)/PlaqView_Master/DataProcessing/references/Tabula_sapiens_reference/TS_Vasculature.h5seurat"
# root.to.TS.mouse.ref <- "~/My Drive (wm5wt@virginia.edu)/PlaqView_Master/DataProcessing/references/Tabula_muris_reference/updated.TS.muris.RDS"

#### Library and Color Schemes ----
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
original_color_list <- {c("rosybrown2",
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


#### Future: Setting Parallel Computing ####
plan("multicore", workers = 3)
options(future.globals.maxSize= 8000 * 1024^2)

#### Function: Process Data ####
plaqview_data_process <- function(datasetID, species.ref = "Human", mitopercentage = 5){
  #### STEP 1: READ DATASET DIRECTORY ####
  # you must change this if your source is different
  # get this to the dataprocessing - data folder in the first piece
  path.to.destination <- file.path(paste(root.to.data,
                                         datasetID, "/source_files", sep=""))
  
  setwd(path.to.destination) 
  print(path.to.destination)
  plaqviewobj <- readRDS(file = "UNPROCESSED.rds")
  plaqviewobj <- UpdateSeuratObject(plaqviewobj)
  
  #### STEP 1B: READ REFERENCES ####
  ## SingleR References ##
  # here we are using Human Primary Cell Atlas design for blood
  # https://bioconductor.org/packages/3.12/data/experiment/vignettes/celldex/inst/doc/userguide.html#2_General-purpose_references
  if(species.ref == "Human"){
    hpca.se <- HumanPrimaryCellAtlasData()
  } else {
    hpca.se <- celldex::MouseRNAseqData()
  }
  
  ## Tabula Sapien References ##
  
  if(species.ref == "Human"){
    TSref <- LoadH5Seurat(file = file.path(root.to.TS.ref), assays = "RNA")
    Idents(TSref) <-  TSref@meta.data[["Annotation"]]
    
  } else {
    TSref <- readRDS(file = file.path(root.to.TS.mouse.ref))
    Idents(TSref) <-  TSref@meta.data[["cell_ontology_class"]]
    TSref@meta.data[["Annotation"]] <-  TSref@meta.data[["cell_ontology_class"]]
    
      }

  #### STEP 2: SEURAT PROCESS ####
  # Run standard cleanup (remove low feature/too many feature/too many mt) (not always needed)
  plaqviewobj[["percent.mt"]] <- PercentageFeatureSet(plaqviewobj, pattern = "^MT-", assay = "RNA")
  
  # this just tells us the distribution of counts
  VlnPlot(plaqviewobj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
  
  # here we will take off the top 1 percentile nfeatures to rm outlier cells
  plaqviewobj <- subset(plaqviewobj, subset = nFeature_RNA > 200 & nFeature_RNA < quantile(plaqviewobj$nFeature_RNA, .99) & percent.mt < mitopercentage)
  
  # Run the standard workflow for visualization and clustering
  plaqviewobj <- FindVariableFeatures(plaqviewobj, verbose = T, nfeatures = 2000, assay = "RNA")
  
  plaqviewobj <- NormalizeData(plaqviewobj, assay = "RNA")
  
  plaqviewobj <- ScaleData(plaqviewobj, verbose = T)
  
  plaqviewobj <- RunPCA(plaqviewobj, npcs = 30, verbose = FALSE)
  plaqviewobj <- RunUMAP(plaqviewobj, reduction = "pca", dims = 1:20)
  plaqviewobj <- FindNeighbors(plaqviewobj, reduction = "pca", dims = 1:20)
  plaqviewobj <- FindClusters(plaqviewobj, resolution = 0.5)
  
  #### STEP 3: SINGLER ----
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

  #### STEP 3A: RECODE SINGLE-R LABELS ----
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


  # #### STEP 3B: SYMPHONY ####
  # ref_pbmcs = readRDS('references/Symphony_ref_data/fibroblast_atlas.rds')
  # 
  # query = symphony::mapQuery(plaqviewobj@assays$RNA, plaqviewobj@meta.data, ref_pbmcs, 
  #                  vars = plaqviewobj@meta.data, 
  #                  do_normalize = TRUE)
  # 
  #### STEP 3C: SEURAT/TABULA SAPIENS LABELING ####
  #### preprocess references  
  TSref <- NormalizeData(TSref, verbose = T)
  TSref <- FindVariableFeatures(TSref, selection.method = "vst", verbose = T)
  
  DefaultAssay(plaqviewobj) <- 'RNA'
  DefaultAssay(TSref) <- 'RNA'
  
  
  anchors <- FindTransferAnchors(reference = TSref, query = plaqviewobj, 
                                 dims = 1:30)
  
  predictions <- TransferData(anchorset = anchors, refdata = TSref$Annotation, 
                              dims = 1:30)
  
  plaqviewobj <- AddMetaData(plaqviewobj, metadata = predictions)
  
  #### rename transferred column metadata 
  plaqviewobj@meta.data[["Seurat_with_Tabula_Ref"]] <- plaqviewobj@meta.data[["predicted.id"]]
  
  # capitalize the lettering
  plaqviewobj@meta.data[["Seurat_with_Tabula_Ref"]] <-str_to_title(plaqviewobj@meta.data[["Seurat_with_Tabula_Ref"]], locale = "en")
  
  # set to active idents
  Idents(plaqviewobj) <- plaqviewobj@meta.data[["Seurat_with_Tabula_Ref"]]
  
  #### STEP 3D: RECODE SEURAT/TABULA ####
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
  
  #### STEP 4: MONOCLE3 TRAJECTORY INFERENCE ----
  
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
  

  
  ## this is species dependent. for human, we will use seurat/ts
  ## for mouse we will use singleR
  if(species.ref == "Human"){
    
    ## transfer seurat labels to moncle3 object
    colData(plaqviewobj.cds)$assigned_cell_type <- 
      plaqviewobj@meta.data[["Seurat_with_Tabula_Ref"]] # call this by opening the object
    
    most.common.cell <- names(sort(table(plaqviewobj@meta.data[["Seurat_with_Tabula_Ref"]]),decreasing = T)[1])
    
    #### MONOCLE3 CONT. ---
    # now learn the PATH (trajectory)
    plaqviewobj.cds <- learn_graph(plaqviewobj.cds)
    
    get_earliest_principal_node <- function(cds, assigned_cell_type = most.common.cell){ # most common human celltype
      cell_ids <- which(colData(cds)[, "assigned_cell_type"] == assigned_cell_type)
      
      closest_vertex <-
        cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
      closest_vertex <- as.matrix(closest_vertex[colnames(cds), ])
      root_pr_nodes <-
        igraph::V(principal_graph(cds)[["UMAP"]])$name[as.numeric(names
                                                                  (which.max(table(closest_vertex[cell_ids,]))))]
      
      root_pr_nodes
    }
  }else{
    ## transfer seurat labels to moncle3 object
    colData(plaqviewobj.cds)$assigned_cell_type <- 
      plaqviewobj@meta.data[["SingleR.calls"]] # call this by opening the object
    
    most.common.cell <- names(sort(table(plaqviewobj@meta.data[["SingleR.calls"]]),decreasing = T)[1])
    
    #### MONOCLE3 CONT. ---
    # now learn the PATH (trajectory)
    plaqviewobj.cds <- learn_graph(plaqviewobj.cds)
    
    
    get_earliest_principal_node <- function(cds, assigned_cell_type = most.common.cell){ # most common mouse celltype
      cell_ids <- which(colData(cds)[, "assigned_cell_type"] == assigned_cell_type)
      
      closest_vertex <-
        cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
      closest_vertex <- as.matrix(closest_vertex[colnames(cds), ])
      root_pr_nodes <-
        igraph::V(principal_graph(cds)[["UMAP"]])$name[as.numeric(names
                                                                  (which.max(table(closest_vertex[cell_ids,]))))]
      
      root_pr_nodes
    }
  }
  

  plaqviewobj.cds <- order_cells(plaqviewobj.cds, 
                                 root_pr_nodes=get_earliest_principal_node(plaqviewobj.cds),
                                 reduction_method = "UMAP")
  
  #### STEP 5: Clean-Up Metadata ####
  # show all metadata columns
  names(plaqviewobj@meta.data)
  
  # rename number clusters  
  plaqviewobj@meta.data$Seurat_Clusters <- plaqviewobj@meta.data$seurat_clusters
  plaqviewobj.all.metadata.headings <- names(plaqviewobj@meta.data)
  write.csv(plaqviewobj.all.metadata.headings, file = "metadataheadings_can_delete.csv")
  # # choose which ones to keep for display
  # plaqviewobj@meta.data <- 
  #   plaqviewobj@meta.data[, which(colnames(plaqviewobj@meta.data)
  #                                 %in% c(
  #                                   "Seurat_Clusters",
  #                                   "Author_Provided",
  #                                   "SingleR.calls",
  #                                   "Seurat_with_Tabula_Ref"  
  #                                 ))]
  #### STEP 6: REDUCE SIZE & SAVE RDS ####
  # plaqviewobj <- DietSeurat(plaqviewobj, counts = T, data = T, dimreducs = c('umap'))
  
  final.file.name <- file.path(paste("../", datasetID, ".rds", sep="")) # ../ moves up one level in file
  final.file.name.cds <- file.path(paste("../", datasetID, "_cds.rds", sep="")) # ../ moves up one level in file
  
  saveRDS(plaqviewobj, file = final.file.name)
  saveRDS(plaqviewobj.cds, file = final.file.name.cds)
  
  
  # plaqviewobj <- readRDS(file = final.file.name)
  
  #### STEP 7: DIFF EX GENE LIST ####
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

  #### STEP 8: PRINT CELL COUNT, RAM gc ####
  
  tab <- as.data.frame(summary(as.factor(plaqviewobj$Seurat_Clusters)))
  
  write_csv(tab, file = paste(datasetID, "_cell_count.csv"))
  
  gc()

  
  #### STEP 9: FINAL CHECK for ANNOTATIONS ####
  
}


#### Function: Check Metadata ####
plaqview_check.metadata <- function(datasetID, species.ref = "Human", is.deployed = TRUE){
  #### STEP 1: READ DATASET DIRECTORY ####
  # you must change this if your source is different
  # get this to the dataprocessing - data folder in the first piece
  if(is.deployed == TRUE){
    processedfile <- file.path(paste(root.to.data, "Archival-Deployed/Public/",
                                     datasetID, "/", datasetID, ".rds", sep=""))
  }else{
    processedfile <- file.path(paste(root.to.data,
                                     datasetID, "/", datasetID, ".rds", sep=""))
  }
  
  print(processedfile)
  
  # read the file
  plaqviewobj <- readRDS(file = processedfile)
  
  # find, sort and append metadata
  plaqviewobj.all.metadata.headings <- str_sort(names(plaqviewobj@meta.data))
  

  # check to see if contains minimium required metadata
  n <- str_count(plaqviewobj.all.metadata.headings, "Author_Provided")
  
  
  
  } # close function

tryCatch(plaqview_check.metadata(datasetID = "Pan_2020"))
tryCatch(plaqview_check.metadata(datasetID = "Wirka_2019"))

#### Process Human Datasets ####
# plaqview_data_process(datasetID = "Alsaigh_2020")
# plaqview_data_process(datasetID = "Li_2020", mitopercentage = 100, upperfeaturelimit = 6500) # this one intentionally seq mitochondrial reads
# plaqview_data_process(datasetID = "Wirka_2019") # rerun 1-22
# plaqview_data_process(datasetID = "Alencar_2020") # rerun 1-22
# tryCatch(plaqview_data_process(datasetID = "Pan_2020")) # rerun 1-26
# tryCatch(plaqview_data_process(datasetID = "Slender_2021")) # rerun 1-26
# plaqview_data_process(datasetID = "Litvinukova_2020")
# plaqview_data_process(datasetID = "Tucker_2020")
# plaqview_data_process(datasetID = "Zernecke_2020", root.cell.type.human = "Monocytes")
# plaqview_data_process(datasetID = "Litvinukova_2020_adipocyte")
# plaqview_data_process(datasetID = "Litvinukova_2020_atrial")
# plaqview_data_process(datasetID = "Litvinukova_2020_fibroblast")
# plaqview_data_process(datasetID = "Litvinukova_2020_immune")
# plaqview_data_process(datasetID = "Litvinukova_2020_neuronal")
# plaqview_data_process(datasetID = "Litvinukova_2020_skeletal")
# plaqview_data_process(datasetID = "Litvinukova_2020_vascular")
# plaqview_data_process(datasetID = "Litvinukova_2020_ventricular")
# plaqview_data_process(datasetID = "Tabula_sapiens_2021_heart")
# plaqview_data_process(datasetID = "Tabula_sapiens_2021_vasculature")
# plaqview_data_process(datasetID = "Delorey_2021")

plaqview_data_process(datasetID = "Xu_2020")


#### Process Mouse Datasets ####
# plaqview_data_process(datasetID = "Alencar_2020_dual", species.ref = "Mouse") # rerun 1-21
# plaqview_data_process(datasetID = "Alencar_2020_KLF4", species.ref = "Mouse") # rerun 1-21
# plaqview_data_process(datasetID = "Pan_2020_mouse", species.ref = "Mouse") # rerun 1-21
# plaqview_data_process(datasetID = "Wirka_2019_mouse", species.ref = "Mouse") # rerun 1-21
# plaqview_data_process(datasetID = "Zernecke_2020_mouse", species.ref = "Mouse")
# plaqview_data_process(datasetID = "Tabula_muris_2019", species.ref = "Mouse")

# plaqview_data_process(datasetID = "vanKuijk_2022_integrated", species.ref = "Mouse")
# plaqview_data_process(datasetID = "vanKuijk_2022_healthy", species.ref = "Mouse")
# plaqview_data_process(datasetID = "Andueza_2020", species.ref = "Mouse")
# plaqview_data_process(datasetID = "Gu_2019", species.ref = "Mouse")
# plaqview_data_process(datasetID = "Dobnikar_2018", species.ref = "Mouse")

