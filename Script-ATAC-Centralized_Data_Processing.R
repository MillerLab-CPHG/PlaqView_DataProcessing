#### LIBRARIES #### 
# library(BiocManager)
library(shiny)
library(shinythemes)
library(Seurat)
library(shinybusy) #install.packages("shinybusy")
library(enrichR) # install.packages("enrichR")
library(waiter)
library(DT)
library(readxl)
library(shinyWidgets)
library(shinyjs)
# library(RColorBrewer)
library(rDGIdb) # BiocManager::install("rDGIdb")
library(tidyverse)
# library(rsconnect)
library(monocle3)
library(ggpubr)
library(gtools)
library(CIPR)
library(ArchR)
library(Signac)
library(parallel)

# library(reactlog)
# library(future)
#
# # tell shiny to log all reactivity
# reactlog_enable()
# 
# # tell shiny to try to paralle compute
# future::plan("multisession")
 


# time code
ptm <- proc.time()
addArchRThreads(threads = 1) # doesnt seem like archr supports multithreading here

#### READ FROM ARCHR ####
proj <- loadArchRProject(path = "data/Turner_2022_ATAC/")

#### RENAME AUTHOR PROVIDED ####
proj$Author_Provided <- proj$Clusters2

saveArchRProject(proj, outputDirectory = "Turner_2022_ATAC_2", overwrite = F)
