data.path <- "./data/"
navin.dir <- paste0(data.path, "navin_annots_final/")
figures.dir <- "../../figures/"
system(paste0("mkdir -p ", figures.dir))

library(Seurat)
library(ggplot2)
library(ggpubr)

load(paste0(data.path, "spatial_list.Rda"))

obj <- spatial_list[[1]]

pdf(paste0(figures.dir, 'FIG_1G_spatial.pdf'), useDingbats = F, width = 10, height = 10)
SpatialFeaturePlot(obj, features=c("Cancer", "Mesenchymal", "Lymphoid", "Endothelial"), min.cutoff=0, max.cutoff=1)
dev.off()

annots_BW11C<-read.csv(paste0(navin.dir, "BW11.C_new.csv"), header = T, fill = T)
rownames(annots_BW11C) <- annots_BW11C[,1]
annots_BW11C[,1] <- NULL
colnames(annots_BW11C) <- c("annots")
obj$navin_annots <- annots_BW11C$annots
obj$navin_annots_2 <- obj$navin_annots
obj$navin_annots_2[obj$navin_annots_2 == 'Adenocarcinoma_lepidic'] = 'Malignant'
obj$navin_annots_2[obj$navin_annots_2 == 'Smooth muscle' | obj$navin_annots_2 == 'Endothelial cell' | obj$navin_annots_2 == 'Fibroblasts'] = 'Stromal'
obj$navin_annots_2[obj$navin_annots_2 == 'Immune cells_lymphocytes'] = 'Immune'

pdf(paste0(figures.dir, 'FIG_1G_spatial_patho.pdf'), useDingbats = F, width = 10, height = 10)
SpatialDimPlot(obj, 'navin_annots_2', cols=c('white', 'blue', 'red', 'green'))
dev.off()

df_list <- list()

for(sample in c('1174.C', '1174.D', '1179.B', 'BW09.A', 'BW09.D', 'BW11.C', 'BW11.D', 'BW14.A', 'BW14.B', 'BW19.C', 'BW19.D', 'BW23.A', 'BW23.B')){

    print(sample)
    obj <- spatial_list[[sample]]
    annots <- read.csv(paste0(navin.dir, sample, "_new.csv"), header = T, fill = T)
    rownames(annots) <- annots[,1]
    annots[,1] <- NULL
    colnames(annots) <- c("annots")
    obj$navin_annots <- annots$annots
    obj$navin_annots_2 <- obj$navin_annots
    obj$navin_annots_2[obj$navin_annots_2 == 'Adenocarcinoma_lepidic'] = 'Malignant'
    obj$navin_annots_2[obj$navin_annots_2 == 'Smooth muscle' | obj$navin_annots_2 == 'Endothelial cell' | obj$navin_annots_2 == 'Fibroblasts'] = 'Stromal'
    obj$navin_annots_2[obj$navin_annots_2 == 'Immune cells_lymphocytes'] = 'Immune'

    df <- data.frame(
        navin_annots_2 = obj@meta.data$navin_annots_2,
        Stromal = obj@meta.data$Pericyte + obj@meta.data$Mesenchymal + obj@meta.data$Endothelial,
        Immune = obj@meta.data$Myeloid + obj@meta.data$Lymphoid + obj@meta.data$NK + obj@meta.data$Plasma + obj@meta.data$Bcell + obj@meta.data$Mast,
        Malignant = obj@meta.data$Cancer
    )
    df <- df[df$navin_annots_2 %in% c("Malignant", "Stromal", "Immune"),]

    df_list[[sample]] <- df
}

merged_df <- do.call(rbind, df_list)

pdf(paste0(figures.dir, "FIG_1G_malig.pdf"), width = 6, height = 6, useDingbats=FALSE)
ggboxplot(merged_df, x = "navin_annots_2", y = "Malignant", color = "navin_annots_2",
                       title = "Malignant") +
stat_compare_means(label.x = 0.9, label.y = .9 * max(as.numeric(merged_df[,"Malignant"]))) +
theme(plot.title = element_text(size = 12, face = "bold")) + theme(legend.position = "none")
dev.off()

pdf(paste0(figures.dir, "FIG_1G_stromal.pdf"), width = 6, height = 6, useDingbats=FALSE)
ggboxplot(merged_df, x = "navin_annots_2", y = "Stromal", color = "navin_annots_2",
                       title = "Stromal") +
stat_compare_means(label.x = 0.9, label.y = .9 * max(as.numeric(merged_df[,"Stromal"]))) +
theme(plot.title = element_text(size = 12, face = "bold")) + theme(legend.position = "none")
dev.off()

pdf(paste0(figures.dir, "FIG_1G_immune.pdf"), width = 6, height = 6, useDingbats=FALSE)
ggboxplot(merged_df, x = "navin_annots_2", y = "Immune", color = "navin_annots_2",
                       title = "Immune") +
stat_compare_means(label.x = 0.9, label.y = .9 * max(as.numeric(merged_df[,"Immune"]))) +
theme(plot.title = element_text(size = 12, face = "bold")) + theme(legend.position = "none")
dev.off()
