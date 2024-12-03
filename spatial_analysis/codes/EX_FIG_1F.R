data.path <- "./data/"
figures.dir <- "../../figures/"
system(paste0("mkdir -p ", figures.dir))

library(ggplot2)
library(ggpubr)
library(gplots)
library(data.table)

source(paste0(data.path, "/R_utils/plotutils.R"))
source(paste0(data.path, "/R_utils/seuratutils.R"))
source(paste0(data.path, "/R_utils/seuratutilsV3.R"))
source(paste0(data.path, "/R_utils/color.R"))

load(paste0(data.path, "spatial_list.Rda"))

prop.list <- list()
for (i in 1:length(spatial_list)) {
  meta <- spatial_list[[i]]@meta.data[,14:24] # high lvl
  means <- as.data.frame(t(as.data.frame(colMeans(meta))))
  prop.list[[i]] <- means
}

prop.df <- rbindlist(prop.list)
prop.df <- as.data.frame(prop.df)
rownames(prop.df) <- names(spatial_list)
prop.df <- as.data.frame(t(prop.df))
Dout <- as.matrix(prop.df)
my_color_palette <- material.heat(nrow(Dout))

pdf(paste0(figures.dir, 'EX_FIG_1F.pdf'),9,5, useDingbats=FALSE)
par(mar=c(5, 4.1, 4.1, 2.1))
barplot2(Dout, main="Cellular composition", legend = rownames(Dout),col = my_color_palette, xlim=c(0, ncol(Dout) + 15), las=2,
        legend.text=TRUE,args.legend=list(x=ncol(Dout)+15,y=max(colSums(Dout)),bty = "n"))
dev.off()
