library(ggplot2)
library(Seurat)

path <- "/path/to/data/"
setwd(path)

sc_obj <- get(load("htan_lung_scrna.Rda")) # load Seurat scRNA-seq object with all cells

DefaultAssay(sc_obj) <- "RNA"
Idents(sc_obj) <- "cell_class"
Idents(sc_obj) <- factor(Idents(all), levels = rev(c("Epithelial", "Myeloid", "Tcell", "NK", "Endothelial", "Lymphatic",
"CAF", "SmoothMuscle", "Pericyte", "Bcell", "Plasma", "Mast")))

# Dotplot of cell class markers
cell_class_markers <- c("EPCAM", "KRT19", "KRT8", "LYZ", "TYROBP", "FCER1G", "TRAC", "CD3D",
"IL32", "GNLY", "NKG7", "GZMB", "VWF", "CLDN5", "PECAM1", "CCL21", "FABP4", "TFF3",
"LUM", "COL3A1", "COL1A1", "ACTA2", "TAGLN", "MYH11", "HIGD1B", "COX4I2", "PTN",
"CD79A", "MS4A1", "CD79B", "IGLC2", "IGKC", "IGHA1", "TPSB2", "TPSAB1", "CPA3")

# Figure 1C
pdf("Dotplot_cell_classes.pdf", width = 12, height = 4)
  DotPlot(object = sc_obj, features = genes, scale = T) +
  theme(axis.text.x = element_text(angle = 45, hjust=1)) +
  geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.5)
dev.off()

# UMAP of cell classes
DefaultAssay(sc_obj) <- "integrated"

my_color_palette <-  c("#283593", "#FFB605", "#3D4EB2", "#268AE8", "#08B2DB", "#32B37C", "#70BA4C", "#ACCF41",
  "#E1E239", "#FFDD29", "#FF8D05", "#FE5523", "#F44336")

# Figure 1D
DimPlot(sc_obj, cols = my_color_palette)
ggsave("Umap_cell_classes.png", width = 8, height = 7)
DimPlot(sc_obj, group.by = "p53_status")
ggsave("Umap_p53_status.png", width = 8, height = 7)

# Cell class proportions in TP53mut vs. TP53WT
Idents(sc_obj) <- "orig.ident"
sc_obj_unsorted <- subset(sc_obj, idents = c("14MA","14ML","17M","19_Norm","MGH_1170","MGH1172a","MGH1172b",
  "MGH1173a","MGH1173b","MGH1174a","MGH1174b","MGH1175a","MGH1175b","MGH1176","MGH1179a","MGH1181",
  "MGH1182","MGH1183","NSC001BLT_FT","NSC004bBN_FT","NSC005bLT_FT","NSC006bLT_FT_V2","NSC009bLT_FT_T1",
  "NSC011bLT_FT_T1","NSC014bLT_FT1","NSC016bLT_FT_T1","NSC019bLTcd45n","NSC023bLT")) # subset of channels that were not sorted by CD45 or EPCAM expression before sequencing

TP53_wt_tumor_ids <- c("MGH14", "MGH1174", "MGH1176", "MGH1183", "BWH01", "BWH04", "BWH05", "BWH09", "BWH19", "BWH23")
TP53_mut_tumor_ids <- c("MGH1172", "MGH1173", "MGH1175", "MGH1182", "BWH16", "BWH11", "BWH14", "BWH06")

cell_class_composition <- as.data.frame(table(sc_obj_unsorted$tumor_id, sc_obj_unsorted$cell_class))
cell_class_composition <- cell_class_composition[c(TP53_wt_tumor_ids, TP53_mut_tumor_ids),]
cell_class_composition$TP53_status <- c(rep("TP53WT", 10), rep("TP53mut", 8)
cell_class_composition$TP53_status <- factor(cell_class_composition$TP53_status, levels = c("TP53WT", "TP53mut"))

cell_class_composition  <- prop.table(cell_class_composition, 1)

for (cell_class in colnames(cell_class_composition)){
  mean.plot <- ggboxplot(cell_class_composition, x = "TP53_status", y = cell_class, color = "TP53_status",
               add = "jitter", title = cell_class) +
               stat_compare_means(label.x = 0.9, label.y = .9 * max(as.numeric(cell_class_composition[,cell_class]))) +
               theme(plot.title = element_text(size = 12, face = "bold")) + NoLegend()
  mean.plots[[cell_class]] <- mean.plot
}

# Figure 1F, left
pdf("Cell_class_composition_TP53mut_vs_TP53WT.pdf"), width = 17, height = 5)
cowplot::plot_grid(plotlist = mean.plots, ncol=8)
dev.off()

# Malignant cells
sc_obj_malig <- get(load("htan_lung_scrna_malig.Rda")) # load Seurat scRNA-seq object with only malignant cells

# Figure 2A
DefaultAssay(sc_obj_malig) <- "integrated"
DimPlot(sc_obj_malig, group.by = "subset")
ggsave("Dimplot_malig_subsets.png", width = 9, height = 7)

# Endothelial cells
sc_obj_endo <- get(load(file = paste0(path, "htan_lung_scrna_endo.Rda"))) # load Seurat scRNA-seq object with only endothelial cells

# Figure 3A
DefaultAssay(sc_obj_endo) <- "integrated"
DimPlot(sc_obj_endo, group.by = "subset")
ggsave("Dimplot_endo_subsets.png", width = 9, height = 7)

# Endothelial subset proportions in TP53mut vs. TP53WT
endo_subset_composition <- as.data.frame(table(sc_obj_endo$tumor_id, sc_obj_endo$subset))
endo_subset_composition <- endo_subset_composition[c(TP53_wt_tumor_ids, TP53_mut_tumor_ids),]
endo_subset_composition$TP53_status <- c(rep("TP53WT", 10), rep("TP53mut", 8)
endo_subset_composition$TP53_status <- factor(endo_subset_composition$TP53_status, levels = c("TP53WT", "TP53mut"))

endo_subset_composition  <- prop.table(endo_subset_composition, 1)

for (endo_subset in colnames(endo_subset_composition)){
  mean.plot <- ggboxplot(endo_subset_composition, x = "TP53_status", y = endo_subset, color = "TP53_status",
               add = "jitter", title = endo_subset) +
               stat_compare_means(label.x = 0.9, label.y = .9 * max(as.numeric(endo_subset_composition[,endo_subset]))) +
               theme(plot.title = element_text(size = 12, face = "bold")) + NoLegend()
  mean.plots[[endo_subset]] <- mean.plot
}

# Figure 3B
pdf("Endo_subset_composition_TP53mut_vs_TP53WT.pdf", width = 17, height = 5)
cowplot::plot_grid(plotlist = mean.plots, ncol=8)
dev.off()

# Mesenchymal cells
sc_obj_mes <- get(load("htan_lung_scrna_mes.Rda")) # load Seurat scRNA-seq object with only mesenchymal cells

# Figure 3E
DefaultAssay(sc_obj_mes) <- "integrated"
DimPlot(sc_obj_mes, group.by = "subset")
ggsave("Dimplot_mes_subsets.png"), width = 9, height = 7)

# Mesenchymal subset proportions in TP53mut vs. TP53WT
mes_subset_composition <- as.data.frame(table(sc_obj_mes$tumor_id, sc_obj_mes$subset))
mes_subset_composition <- mes_subset_composition[c(TP53_wt_tumor_ids, TP53_mut_tumor_ids),]
mes_subset_composition$TP53_status <- c(rep("TP53WT", 10), rep("TP53mut", 8)
mes_subset_composition$TP53_status <- factor(mes_subset_composition$TP53_status, levels = c("TP53WT", "TP53mut"))

mes_subset_composition  <- prop.table(mes_subset_composition, 1)

for (mes_subset in colnames(mes_subset_composition)){
  mean.plot <- ggboxplot(mes_subset_composition, x = "TP53_status", y = mes_subset, color = "TP53_status",
               add = "jitter", title = mes_subset) +
               stat_compare_means(label.x = 0.9, label.y = .9 * max(as.numeric(mes_subset_composition[,mes_subset]))) +
               theme(plot.title = element_text(size = 12, face = "bold")) + NoLegend()
  mean.plots[[mes_subset]] <- mean.plot
}

# Figure 3H
pdf("Mes_subset_composition_TP53mut_vs_TP53WT.pdf", width = 17, height = 5)
cowplot::plot_grid(plotlist = mean.plots, ncol=8)
dev.off()

# Myeloid cells
sc_obj_myl <- get(load("htan_lung_scrna_myl.Rda")) # load Seurat scRNA-seq object with only myeloid cells

# Figure 4A
DefaultAssay(sc_obj_myl) <- "integrated"
DimPlot(sc_obj_myl, group.by = "subset")
ggsave("Dimplot_myl_subsets.png"), width = 9, height = 7)

# Myeloid subset proportions in TP53mut vs. TP53WT
myl_subset_composition <- as.data.frame(table(sc_obj_myl$tumor_id, sc_obj_myl$subset))
myl_subset_composition <- myl_subset_composition[c(TP53_wt_tumor_ids, TP53_mut_tumor_ids),]
myl_subset_composition$TP53_status <- c(rep("TP53WT", 10), rep("TP53mut", 8)
myl_subset_composition$TP53_status <- factor(myl_subset_composition$TP53_status, levels = c("TP53WT", "TP53mut"))

myl_subset_composition  <- prop.table(myl_subset_composition, 1)

for (myl_subset in colnames(myl_subset_composition)){
  mean.plot <- ggboxplot(myl_subset_composition, x = "TP53_status", y = myl_subset, color = "TP53_status",
               add = "jitter", title = myl_subset) +
               stat_compare_means(label.x = 0.9, label.y = .9 * max(as.numeric(myl_subset_composition[,myl_subset]))) +
               theme(plot.title = element_text(size = 12, face = "bold")) + NoLegend()
  mean.plots[[myl_subset]] <- mean.plot
}

# Figure 4C
pdf("Myl_subset_composition_TP53mut_vs_TP53WT.pdf", width = 17, height = 5)
cowplot::plot_grid(plotlist = mean.plots, ncol=8)
dev.off()

# T/NK cells
sc_obj_t_nk <- get(load("htan_lung_scrna_t_nk.Rda")) # load Seurat scRNA-seq object with only t and nk cells

# Figure 5A
DefaultAssay(sc_obj_t_nk) <- "integrated"
DimPlot(sc_obj_t_nk, group.by = "subset")
ggsave("Dimplot_t_nk_subsets.png"), width = 9, height = 7)

# T and NK subset proportions in TP53mut vs. TP53WT
t_nk_subset_composition <- as.data.frame(table(sc_obj_t_nk$tumor_id, sc_obj_t_nk$subset))
t_nk_subset_composition <- t_nk_subset_composition[c(TP53_wt_tumor_ids, TP53_mut_tumor_ids),]
t_nk_subset_composition$TP53_status <- c(rep("TP53WT", 10), rep("TP53mut", 8)
t_nk_subset_composition$TP53_status <- factor(t_nk_subset_composition$TP53_status, levels = c("TP53WT", "TP53mut"))

t_nk_subset_composition  <- prop.table(t_nk_subset_composition, 1)

for (t_nk_subset in colnames(t_nk_subset_composition)){
  mean.plot <- ggboxplot(t_nk_subset_composition, x = "TP53_status", y = t_nk_subset, color = "TP53_status",
               add = "jitter", title = t_nk_subset) +
               stat_compare_means(label.x = 0.9, label.y = .9 * max(as.numeric(t_nk_subset_composition[,t_nk_subset]))) +
               theme(plot.title = element_text(size = 12, face = "bold")) + NoLegend()
  mean.plots[[t_nk_subset]] <- mean.plot
}

# Figure 5C
pdf("T_NK_subset_composition_TP53mut_vs_TP53WT.pdf", width = 17, height = 5)
cowplot::plot_grid(plotlist = mean.plots, ncol=8)
dev.off()

# Spatial analysis
path <- "/path/to/obj/"
st_obj <- get(load(file = paste0(path, "htan_lung_spatial.Rda"))) # load Seurat object with spatial transcriptomic data
