library(Seurat)
library(dplyr)

# Load the ECKO dataset
ECKO.data <- Read10X(data.dir = "~/Desktop/scRNAseq_Data/KOPOS_filtered_gene_bc_matrices/mm10witheyfp/")

# Examine the memory savings between regular and sparse matrices
dense.size <- object.size(x = as.matrix(x = ECKO.data))
dense.size

# Initialize the Seurat object with the raw (non-normalized data).  Keep all
# genes expressed in >= 3 cells (~0.1% of the data). Keep all cells with at
# least 200 detected genes
ECKO <- CreateSeuratObject(raw.data = ECKO.data, min.cells = 3, min.genes = 200, 
                           project = "10X_ECKO")

# The number of genes and UMIs (nGene and nUMI) are automatically calculated
# for every object by Seurat.  For non-UMI data, nUMI represents the sum of
# the non-normalized values within a cell We calculate the percentage of
# mitochondrial genes here and store it in percent.mito using AddMetaData.
# We use object@raw.data since this represents non-transformed and
# non-log-normalized counts The % of UMI mapping to MT-genes is a common
# scRNA-seq QC metric.
mito.genes <- grep(pattern = "^MT-", x = rownames(x = ECKO@data), value = TRUE)
percent.mito <- Matrix::colSums(ECKO@raw.data[mito.genes, ])/Matrix::colSums(ECKO@raw.data)

# AddMetaData adds columns to object@meta.data, and is a great place to
# stash QC stats
ECKO <- AddMetaData(object = ECKO, metadata = percent.mito, col.name = "percent.mito")
VlnPlot(object = ECKO, features.plot = c("nGene", "nUMI", "percent.mito"), nCol = 3)

# GenePlot is typically used to visualize gene-gene relationships, but can
# be used for anything calculated by the object, i.e. columns in
# object@meta.data, PC scores etc.  Since there is a rare subset of cells
# with an outlier level of high mitochondrial percentage and also low UMI
# content, we filter these as well
par(mfrow = c(1, 2))
GenePlot(object = ECKO, gene1 = "nUMI", gene2 = "percent.mito")
GenePlot(object = ECKO, gene1 = "nUMI", gene2 = "nGene")

# We filter out cells that have unique gene counts over 2,500 or less than
# 200 Note that low.thresholds and high.thresholds are used to define a
# 'gate'.  -Inf and Inf should be used if you don't want a lower or upper
# threshold.
ECKO <- FilterCells(object = ECKO, subset.names = c("nGene", "percent.mito"), 
                    low.thresholds = c(200, -Inf), high.thresholds = c(5000, 0.05))

# Normalizing the data
ECKO <- NormalizeData(object = ECKO, normalization.method = "LogNormalize", 
                      scale.factor = 10000)

# Detection of variable genes across the single cells
ECKO <- FindVariableGenes(object = ECKO, mean.function = ExpMean, dispersion.function = LogVMR, 
                          x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5)

# Scaling the data and removing unwanted sources of variation
ECKO <- ScaleData(object = ECKO, vars.to.regress = c("nUMI", "percent.mito"))

# Perform linear dimensional reduction
ECKO <- RunPCA(object = ECKO, pc.genes = ECKO@var.genes, do.print = TRUE, pcs.print = 1:5, 
               genes.print = 5)

# Examine and visualize PCA results a few different ways
PrintPCA(object = ECKO, pcs.print = 1:5, genes.print = 5, use.full = FALSE)

VizPCA(object = ECKO, pcs.use = 1:5)

PCAPlot(object = ECKO, dim.1 = 1, dim.2 = 2)

# ProjectPCA scores each gene in the dataset (including genes not included
# in the PCA) based on their correlation with the calculated components.
# Though we don't use this further here, it can be used to identify markers
# that are strongly correlated with cellular heterogeneity, but may not have
# passed through variable gene selection.  The results of the projected PCA
# can be explored by setting use.full=T in the functions above
ECKO <- ProjectPCA(object = ECKO, do.print = FALSE)

PCHeatmap(object = ECKO, pc.use = 1, cells.use = 500, do.balanced = TRUE, label.columns = FALSE)

PCHeatmap(object = ECKO, pc.use = 1:12, cells.use = 500, do.balanced = TRUE, 
          label.columns = FALSE, use.full = FALSE)

# Determine statistically significant principal components
# NOTE: This process can take a long time for big datasets, comment out for
# expediency.  More approximate techniques such as those implemented in
# PCElbowPlot() can be used to reduce computation time
ECKO <- JackStraw(object = ECKO, num.replicate = 100, display.progress = FALSE)

JackStrawPlot(object = ECKO, PCs = 1:18)
# PC11 fails to reach significance and therefore we will use the first 10 PCs going forward

PCElbowPlot(object = ECKO)

# save.SNN = T saves the SNN so that the clustering algorithm can be rerun
# using the same graph but with a different resolution value (see docs for
# full details)
ECKO <- FindClusters(object = ECKO, reduction.type = "pca", dims.use = 1:10, 
                     resolution = 1, print.output = 0, save.SNN = TRUE)

PrintFindClustersParams(object = ECKO)

# Run Non-linear dimensional reduction (tSNE)
ECKO <- RunTSNE(object = ECKO, dims.use = 1:10, do.fast = TRUE)

# note that you can set do.label=T to help label individual clusters
TSNEPlot(object = ECKO)

saveRDS(ECKO, file = "~/Projects/datasets/pbmc3k/pbmc_tutorial.rds")

# find all markers of cluster 1
cluster1.markers <- FindMarkers(object = ECKO, ident.1 = 1, min.pct = 0.25)
print(x = head(x = cluster1.markers, n = 5))

# How well does a marker distinguish the cluster?
cluster1.markers <- FindMarkers(object = ECKO, ident.1 = 0, thresh.use = 0.25, 
                                test.use = "roc", only.pos = TRUE)

# find markers for every cluster compared to all remaining cells, report
# only the positive ones
ECKO.markers <- FindAllMarkers(object = ECKO, only.pos = TRUE, min.pct = 0.25, 
                               thresh.use = 0.25)
ECKO.markers %>% group_by(cluster) %>% top_n(2, avg_logFC)


VlnPlot(object = ECKO, features.plot = c("Lrg1", "Pi16", "Igfbp3", "Serpine1", "Igfbp7","Col3a1","Sfrp1","Ephx2","Comp","Col2a1","Ankrd1", "Sfrp2"))

FeaturePlot(object = ECKO, features.plot = c("Acta2", "Cdh5", "Lgals3", "Il1r1"), cols.use = c("grey", "blue"), 
            reduction.use = "tsne")

top10 <- ECKO.markers %>% group_by(cluster) %>% top_n(10, avg_logFC)
# setting slim.col.label to TRUE will print just the cluster IDS instead of
# every cell name
DoHeatmap(object = ECKO, genes.use = top10$gene, slim.col.label = TRUE, remove.key = TRUE)

current.cluster.ids <- c(0, 1, 2, 3, 4, 5)
new.cluster.ids <- c("CD4 T cells", "CD14+ Monocytes", "B cells", "CD8 T cells", 
                     "FCGR3A+ Monocytes", "NK cells", "Dendritic cells", "Megakaryocytes")
pbmc@ident <- plyr::mapvalues(x = pbmc@ident, from = current.cluster.ids, to = new.cluster.ids)
TSNEPlot(object = pbmc, do.label = TRUE, pt.size = 0.5)
