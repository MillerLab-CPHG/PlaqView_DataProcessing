#### installation ####

# devtools::install_github("omicsCore/scTyper")
# http://htmlpreview.github.io/?https://github.com/omicsCore/scTyper/blob/master/vignettes/Sample_analysis.html

#### test ####
library(scTyper)

plaqviewobj <- readRDS(file = "data/Alencar_2020/Alencar_2020.rds")

### this is designed for cancer data anylsis so we wont go thru this