data.path <- "./data/"
figures.dir <- "../../figures/"
system(paste0("mkdir -p ", figures.dir))

source(paste0(data.path, "/R_utils/plotutils.R"))
source(paste0(data.path, "/R_utils/seuratutils.R"))
source(paste0(data.path, "/R_utils/seuratutilsV3.R"))
source(paste0(data.path, "/R_utils/color.R"))

library(ggplot2)
library(ggpubr)
library(ComplexHeatmap)
library(gdata)
library(data.table)

load(paste0(data.path, "spatial_list.Rda"))

spatial_list <- spatial_list[c(1:4, 7:16)]

patient.corr.list <- list()
for (i in 1:length(spatial_list)) { # filtered - no BW16, etc.
  obj <- spatial_list[[i]]@meta.data
  corr <- cor(obj[,c(14:20, 23:24, 78:90)], obj[,c(14:20, 23:24, 78:90)]) # 2v4 mes subtype, high lvl
  patient.corr.list[[i]] <- corr
}

avg.mtx <- Reduce("+", patient.corr.list) / length(patient.corr.list)

# avg corr all slides
range <- 0.3

pdf(file = paste0(figures.dir, "EX_FIG_7F.pdf"), useDingbats = F, width = 6, height = 6)
pheatmap(avg.mtx, breaks = seq(-range, range, length.out = 100), cluster_cols = F, cluster_rows = F)
dev.off()