#### library ####
library(Seurat)

#### optional processing
#### use to reduce seurat obj sizes!

#### stanford file reduction ####
original <- read_rds(file =  "Wirka_with_TS_05032021.rds")
slim <- DietSeurat(original, counts = T, data = T, dimreducs = c('umap'))

saveRDS(slim, file = "Wirka_2019_slim.rds")
