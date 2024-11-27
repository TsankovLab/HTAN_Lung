data.path <- "./data/"
figures.dir <- "../../figures/"
system(paste0("mkdir -p ", figures.dir))

library(Seurat)
library(ggplot2)
library(ggpubr)

load(paste0(data.path, "spatial_list.Rda"))

spatial_list <- spatial_list[c(1:4, 7:16)]

obj <- spatial_list[[4]]

pdf(paste0(figures.dir, 'FIG_7D.pdf'), useDingbats = F, width = 10, height = 10)
SpatialFeaturePlot(obj, features=c('TAM.SPP1', 'CAF.COL', 'Myofibroblast', 'hypoxia.hallmark1', 'emt.hallmark1', 'mean_nUMI_factorsfact_14'))
dev.off()