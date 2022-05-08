#### PreReq Codes ####
# below line is commented for shinyapp.io deployment temp
### set this once in terminal before deploying to shinyapps.io ###
# options(repos = BiocManager::repositories())

# enrichR functions
# handcurate db names 
dbs <- c("KEGG_2019_Human",
         "WikiPathways_2019_Human",
         "GO_Biological_Process_2018",
         "ChEA_2016",
         "GWAS_Catalog_2019",
         "ENCODE_and_ChEA_Consensus_TFs_from_ChIP-X",
         "Gene_Perturbations_from_GEO_down",
         "Gene_Perturbations_from_GEO_up")
enrichRdb <- sort(dbs)

# color definitions
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

#### LIBRARIES #### 
# library(BiocManager)
library(shiny)
library(shinythemes)
library(Seurat)
library(shinybusy) #install.packages("shinybusy")
library(enrichR) # install.packages("enrichR")
# library(imager)
library(waiter)
library(DT)
library(readxl)
library(shinyWidgets)
library(shinyjs)
# library(RColorBrewer)
library(rDGIdb) # BiocManager::install("rDGIdb")
library(tidyverse)
# library(rsconnect)
# library(CellChat)

#### WIRKA- PRDM16 ####
data <- readRDS(file = "data/Wirka_2019/Wirka_2019.rds")

feat <- c("PRDM16")

RidgePlot(data,
          cols = manual_color_list,
          group.by = "Author_Provided",
          features =  feat,
          log = T,
          y.max = 1,
          same.y.lims = T) + # a trick to sep long string input
  #theme(legend.position="bottom", legend.box = "vertical") + # group.by is important, use this to call metadata separation
  theme(plot.title = element_text(hjust =  0.5)) +
  guides(color = guide_legend(nrow = 3))

VlnPlot(data,
        cols = manual_color_list,
        group.by = "Author_Provided",
        features =  "PRDM16",
        log = F,
        y.max = 1,
        same.y.lims = T) + # a trick to sep long string input
  #theme(legend.position="bottom", legend.box = "vertical") + # group.by is important, use this to call metadata separation
  theme(plot.title = element_text(hjust =  0.5)) +
  guides(color = guide_legend(nrow = 3))

pdf(file = "PRDM16-wirka.pdf", width = 5, height = 5)
FeaturePlot(data,
          features =  feat, pt.size = 1, order = T
          ) +
  #theme(legend.position="bottom", legend.box = "vertical") + # group.by is important, use this to call metadata separation
  theme(plot.title = element_text(hjust =  0.5))
dev.off()


#### WIRKA- TRIP4 ####
data <- readRDS(file = "data/Archival-Deployed/Public/Wirka_2019/Wirka_2019.rds")

feat <- c("TRIP4")

RidgePlot(data,
          cols = manual_color_list,
          group.by = "Author_Provided",
          features =  feat,
          log = T,
          y.max = 1,
          same.y.lims = T) + # a trick to sep long string input
  #theme(legend.position="bottom", legend.box = "vertical") + # group.by is important, use this to call metadata separation
  theme(plot.title = element_text(hjust =  0.5)) +
  guides(color = guide_legend(nrow = 3))

VlnPlot(data,
        cols = manual_color_list,
        group.by = "Author_Provided",
        features =  "TRIP4",
        log = F,
        y.max = 1,
        same.y.lims = T) + # a trick to sep long string input
  #theme(legend.position="bottom", legend.box = "vertical") + # group.by is important, use this to call metadata separation
  theme(plot.title = element_text(hjust =  0.5)) +
  guides(color = guide_legend(nrow = 3))

pdf(file = "TRIP4-wirka.pdf", width = 5, height = 5)
FeaturePlot(data,
            features =  feat, pt.size = 1, order = T
) +
  #theme(legend.position="bottom", legend.box = "vertical") + # group.by is important, use this to call metadata separation
  theme(plot.title = element_text(hjust =  0.5))
dev.off()



#### WIRKA-TBX2 ####


feat <- c("TBX2")

pdf(file = "TBX2-wirka.pdf", width = 5, height = 5)
FeaturePlot(data,
            features =  feat, pt.size = 1, order = T
) +
  #theme(legend.position="bottom", legend.box = "vertical") + # group.by is important, use this to call metadata separation
  theme(plot.title = element_text(hjust =  0.5))
dev.off()

#### WIRKA-UMAP  ####

pdf(file = "UMAP-wirka.pdf", width = 6, height = 5)

Seurat::DimPlot(data, 
                    cols = manual_color_list,
                    group.by = "Author_Provided",
) +
  #theme(legend.position="bottom", legend.box = "vertical") + # group.by is important, use this to call metadata separation
  theme(plot.title = element_text(hjust =  0.5))
dev.off()




#### PAN- PRDM16 ####
data <- readRDS(file = "data/Pan_2020_mouse/Pan_2020_mouse.rds")

feat <- c("Prdm16")

pdf(file = "PRDM16-Pan.pdf", width = 5, height = 5)
FeaturePlot(data,
            features =  feat, pt.size = 1, order = T) +
  #theme(legend.position="bottom", legend.box = "vertical") + # group.by is important, use this to call metadata separation
  theme(plot.title = element_text(hjust =  0.5))
dev.off()


#### Pan-TBX2 ####
feat <- c("Tbx2")

pdf(file = "TBX2-Pan.pdf", width = 5, height = 5)
FeaturePlot(data,
            features =  feat, pt.size = 1, order = T
) +
  #theme(legend.position="bottom", legend.box = "vertical") + # group.by is important, use this to call metadata separation
  theme(plot.title = element_text(hjust =  0.5))
dev.off()

#### Pan-UMAP  ####

pdf(file = "UMAP-Pan.pdf", width = 6, height = 5)

Seurat::DimPlot(data, 
                cols = manual_color_list,
                group.by = "Author_Provided",
) +
  #theme(legend.position="bottom", legend.box = "vertical") + # group.by is important, use this to call metadata separation
  theme(plot.title = element_text(hjust =  0.5))
dev.off()





#### WIRKA- long list look up ####
data <- readRDS(file = "data/Deployed/Public/Wirka_2019/Wirka_2019.rds")
CACgenes_Wei <- read_excel("~/Documents/Desktop/CACgenes_Wei.xlsx", 
                           col_names = FALSE)
feat <- CACgenes_Wei$...1


pdf(file = "wirka-lookup-dec172021.pdf", width = 10, height = 8)
DotPlot(data,
        # cols = manual_color_list,
        group.by = "Seurat_with_Tabula_Ref",
        features =  feat) + # a trick to sep long string input
  #theme(legend.position="bottom", legend.box = "vertical") + # group.by is important, use this to call metadata separation
  theme(plot.title = element_text(hjust =  0.5)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +  
  guides(color = guide_legend(nrow = 3))


dev.off()



#### PAN- long list look up ####
data <- readRDS(file = "data/Deployed/Public/Pan_2020/Pan_2020.rds")
CACgenes_Wei <- read_excel("~/Documents/Desktop/CACgenes_Wei.xlsx", 
                           col_names = FALSE)
feat <- CACgenes_Wei$...1

pdf(file = "pan-lookup-dec172021.pdf", width = 10, height = 8)
DotPlot(data,
        # cols = manual_color_list,
        group.by = "Seurat_with_Tabula_Ref",
        features =  feat) + # a trick to sep long string input
  #theme(legend.position="bottom", legend.box = "vertical") + # group.by is important, use this to call metadata separation
  theme(plot.title = element_text(hjust =  0.5)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +  
  guides(color = guide_legend(nrow = 3))

dev.off()


#### LOTTE- long list look up ####
data <- readRDS(file = "data/Deployed/NoSharing/Slender_2021/Slender_2021.rds")
CACgenes_Wei <- read_excel("~/Documents/Desktop/CACgenes_Wei.xlsx", 
                           col_names = FALSE)
feat <- CACgenes_Wei$...1

pdf(file = "slender-lookup-dec172021.pdf", width = 10, height = 8)
DotPlot(data,
        # cols = manual_color_list,
        group.by = "Seurat_with_Tabula_Ref",
        features =  feat) + # a trick to sep long string input
  #theme(legend.position="bottom", legend.box = "vertical") + # group.by is important, use this to call metadata separation
  theme(plot.title = element_text(hjust =  0.5)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +  
  guides(color = guide_legend(nrow = 3))


dev.off()



#### ALSAIGH- long list look up ####
data <- readRDS(file = "data/Archiva-Deployed/Public/Alsaigh_2020/Alsaigh_2020.rds")
CACgenes_Wei <- readxl::read_excel("Special Lookup Results for Clint/cacgenes.xlsx", 
                           col_names = FALSE)
feat <- CACgenes_Wei$...1

pdf(file = "alsaigh-lookup-dec01142022.pdf", width = 10, height = 8)
DotPlot(data,
        # cols = manual_color_list,
        group.by = "Seurat_with_Tabula_Ref",
        features =  feat) + # a trick to sep long string input
  #theme(legend.position="bottom", legend.box = "vertical") + # group.by is important, use this to call metadata separation
  theme(plot.title = element_text(hjust =  0.5)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +  
  guides(color = guide_legend(nrow = 3))


dev.off()


