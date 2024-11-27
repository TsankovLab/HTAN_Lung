data.path <- "./data/"
figures.dir <- "../../figures/"
system(paste0("mkdir -p ", figures.dir))

library(Seurat)
library(ggplot2)
library(ggpubr)
library(reshape)

load(paste0(data.path, "spatial_list.Rda"))

spatial_list <- spatial_list[c(1:4, 7:16)]

obj <- spatial_list[['BW14.A']]@meta.data[c('Mesenchymal', 'Myeloid', 'Cancer', 'mean_nUMI_factorsfact_14')]
obj <- obj[order(obj$mean_nUMI_factorsfact_14), ]

mdata <- melt(obj, id= "mean_nUMI_factorsfact_14")

pdf(paste0(figures.dir, 'EX_FIG_15C.pdf'), useDingbats = F, height = 4, width = 6)
print(ggplot(mdata, aes(mean_nUMI_factorsfact_14, value, col = variable, group = variable)) +
geom_smooth(se=F) + theme_bw())
dev.off()