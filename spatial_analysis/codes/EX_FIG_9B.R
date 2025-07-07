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
  # change this for every cell type I'm looking at
  corr <- cor(obj[,c(14:15, 21:24, 59:77)], obj[,c(14:15, 21:24, 59:77)]) # 2v4 bmast, tnk subtype, high lvl
  patient.corr.list[[i]] <- corr
}

avg.mtx <- Reduce("+", patient.corr.list) / length(patient.corr.list)

order <- c("Cancer", "Nonmalig", "CD8.IFN", "T.Exhausted", "CD8.M.E", "CD8.Naive", "B.Follicular", "Plasma.IGHA", "Plasma.IGHG",
"TH17", "TFH","CD4.M.E", "CD4.Naive","Treg", "Myeloid","NK", "T.Delta.Gamma", "T.Cycling", "Plasmacytoid.DC", "Mesenchymal",
"Pericyte", "Endothelial", "Mast", "B.MZ", "B.Cycling")

real_order <- c("Cancer", "Nonmalig", "CD8.IFN", "T.Exhausted", "CD8.GZMK", "CD8.TRM", "B.Follicular", "Plasma.IGHA", "Plasma.IGHG",
"T.Stress", "TFH","CD4.TRM", "CD4.Naive","Treg", "Myeloid","NK.CD56.dim", "NK.CD56.bright", "T.Cycling", "Plasmacytoid.DC", "Mesenchymal",
"Pericyte", "Endothelial", "Mast", "B.MZ", "B.Cycling")

plt.avg.mtx <- avg.mtx
plt.avg.mtx <- plt.avg.mtx[order,order]

colnames(plt.avg.mtx) <- real_order
rownames(plt.avg.mtx) <- real_order

# avg corr all slides
range <- 0.3

pdf(file = paste0(figures.dir, "EX_FIG_9B_mean.pdf"), useDingbats = F, width = 6, height = 6)
pheatmap(plt.avg.mtx, breaks = seq(-range, range, length.out = 100), cluster_cols = F, cluster_rows = F)
dev.off()

patient.corr.list2 <- list()
for (i in 1:length(spatial_list)) { # filtered - no BW16, etc.
  obj <- spatial_list[[i]]@meta.data
  corr <- cor(obj[,c(14:15, 21:24, 59:77)], obj[,c(14:15, 21:24, 59:77)]) # 2v4 bmast, tnk subtype, high lvl
  patient.corr.list2[[i]] <- as.data.frame(t(as.data.frame(unmatrix(corr))))
}

rbind.df <- t(rbindlist(patient.corr.list2))
pvals.list <- apply(rbind.df, 1, function(x) {
  pval <- wilcox.test(x[1:4], x[5:14])$p.value
  pval <- -log10(pval)
})

avg.mtx <- t(avg.mtx) # not for x vs factors
avg.mtx.melt <- reshape2::melt(avg.mtx)
avg.mtx.melt$pvals <- pvals.list

meandiff.list <- apply(rbind.df, 1, function(x) {
  meandiff <- mean(x[1:4], na.rm = T) - mean(x[5:14], na.rm = T)
})

avg.mtx.melt$meandiff <- meandiff.list
avg.mtx.melt$signif <- ""
avg.mtx.melt$signif[avg.mtx.melt$pvals >= 1.3] <- "*"

avg.mtx.melt2 <- avg.mtx.melt
avg.mtx.melt2 <- reshape2::dcast(avg.mtx.melt2, Var1 ~ Var2, value.var = "meandiff") # why does it change...
rownames(avg.mtx.melt2) <- avg.mtx.melt2$Var1
avg.mtx.melt2$Var1 <- NULL

avg.mtx.melt3 <- avg.mtx.melt
avg.mtx.melt3 <- reshape2::dcast(avg.mtx.melt3, Var1 ~ Var2, value.var = "signif")
rownames(avg.mtx.melt3) <- avg.mtx.melt3$Var1
avg.mtx.melt3$Var1 <- NULL

avg.mtx.melt2 <- t(avg.mtx.melt2)
avg.mtx.melt3 <- t(avg.mtx.melt3)

range <- 0.3

avg.mtx.melt2 <- avg.mtx.melt2[order,order]
avg.mtx.melt3 <- avg.mtx.melt3[order,order]

colnames(avg.mtx.melt2) <- real_order
rownames(avg.mtx.melt2) <- real_order
colnames(avg.mtx.melt3) <- real_order
rownames(avg.mtx.melt3) <- real_order

pdf(file = paste0(figures.dir, "EX_FIG_9B_diff.pdf"), useDingbats = F, width = 6, height = 6)
pheatmap(avg.mtx.melt2, display_numbers = avg.mtx.melt3, border_color="gray", breaks = seq(-range, range, length.out = 100), cluster_cols = F, cluster_rows = F)
dev.off()
