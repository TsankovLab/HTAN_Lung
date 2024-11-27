data.path <- "./data/"
figures.dir <- "../../figures/"
system(paste0("mkdir -p ", figures.dir))

library(Seurat)
library(ggplot2)
library(ggpubr)

spatial_signif <- read.csv(file = paste0(data.path, "spatial_list_u50_factor14_highlvl_sigs.corr.df.csv"), check.names = F)

colnames(spatial_signif)[1] <- "cellpair"

subset <- spatial_signif[spatial_signif$cellpair %in% "hypoxia.hallmark1:mean_nUMI_factorsfact_14",]
df.t <- t(subset[,c(2:15)])
df.t <- as.data.frame(df.t)
df.t$p53_status <- c(rep("p53_mut", 4), rep("WT", 10))
df.t$p53_status  <- factor(df.t$p53_status , levels = c("WT", "p53_mut"))
colnames(df.t) <- c("Hallmark hypoxia-NMF7", "p53_status")

pdf(paste0(figures.dir, "FIG_7E_hallmarkhypoxia_nmf.pdf"), width = 3, height = 3, useDingbats=FALSE)
ggboxplot(df.t, x = "p53_status", y = "Hallmark hypoxia-NMF7", color = "p53_status",
                       palette = c("#0091CA", "#D8423D"), add = "jitter",
                       title = "Hallmark hypoxia-NMF7 Spatial Corr") +
stat_compare_means(label.x = 0.9, label.y = .9 * max(as.numeric(df.t[,"Hallmark hypoxia-NMF7"]))) +
theme(plot.title = element_text(size = 12, face = "bold")) + theme(legend.position = "none")
dev.off()

subset <- spatial_signif[spatial_signif$cellpair %in% "emt.hallmark1:mean_nUMI_factorsfact_14",]
df.t <- t(subset[,c(2:15)])
df.t <- as.data.frame(df.t)
df.t$p53_status <- c(rep("p53_mut", 4), rep("WT", 10))
df.t$p53_status  <- factor(df.t$p53_status , levels = c("WT", "p53_mut"))
colnames(df.t) <- c("Hallmark EMT-NMF7", "p53_status")

pdf(paste0(figures.dir, "FIG_7E_hallmarkemt_nmf.pdf"), width = 3, height = 3, useDingbats=FALSE)
ggboxplot(df.t, x = "p53_status", y = "Hallmark EMT-NMF7", color = "p53_status",
                       palette = c("#0091CA", "#D8423D"), add = "jitter",
                       title = "Hallmark EMT-NMF7 Spatial Corr") +
stat_compare_means(label.x = 0.9, label.y = .9 * max(as.numeric(df.t[,"Hallmark EMT-NMF7"]))) +
theme(plot.title = element_text(size = 12, face = "bold")) + theme(legend.position = "none")
dev.off()

subset <- spatial_signif[spatial_signif$cellpair %in% "Endothelial:mean_nUMI_factorsfact_14",]
df.t <- t(subset[,c(2:15)])
df.t <- as.data.frame(df.t)
df.t$p53_status <- c(rep("p53_mut", 4), rep("WT", 10))
df.t$p53_status  <- factor(df.t$p53_status , levels = c("WT", "p53_mut"))
colnames(df.t) <- c("Endothelial-NMF7", "p53_status")

pdf(paste0(figures.dir, "FIG_7E_endo_nmf.pdf"), width = 3, height = 3, useDingbats=FALSE)
ggboxplot(df.t, x = "p53_status", y = "Endothelial-NMF7", color = "p53_status",
                       palette = c("#0091CA", "#D8423D"), add = "jitter",
                       title = "Endothelial-NMF7 Spatial Corr") +
stat_compare_means(label.x = 0.9, label.y = .9 * max(as.numeric(df.t[,"Endothelial-NMF7"]))) +
theme(plot.title = element_text(size = 12, face = "bold")) + theme(legend.position = "none")
dev.off()