# contains feature plots for all genes for all 14 samples - Fig 3D, 3I, 4I, and 5G

data.path <- "./data/"
figures.dir <- "../../figures/spatial_feature_plots/"
system(paste0("mkdir -p ", figures.dir))

library(Seurat)
library(ggplot2)
library(ggpubr)

load(paste0(data.path, "spatial_list.Rda"))

spatial_list <- spatial_list[c(1:4, 7:16)]

genes <- c('NRP1', 'SEMA3A', 'TGFB2', 'TGFBR2', 'SPP1', 'CD44', 'PVR', 'TIGIT')

for (i in 1:length(spatial_list)) {
    pdf(paste0(figures.dir, names(spatial_list)[i], '.pdf'), useDingbats = F, width = 20, height = 20)
    print(SpatialFeaturePlot(spatial_list[[i]], features=genes))
    dev.off()
}