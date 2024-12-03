data.path <- "./data/"
navin.dir <- paste0(data.path, "navin_annots_final/")
figures.dir <- "../../figures/"
system(paste0("mkdir -p ", figures.dir))

library(Seurat)
library(ggplot2)
library(ggpubr)
library(gplots)
library(ComplexHeatmap)

source(paste0(data.path, "/R_utils/plotutils.R"))
source(paste0(data.path, "/R_utils/seuratutils.R"))
source(paste0(data.path, "/R_utils/seuratutilsV3.R"))
source(paste0(data.path, "/R_utils/color.R"))

load(paste0(data.path, "spatial_list.Rda"))

annot_1174C <- read.csv(paste0(navin.dir, "1174.C_new.csv"), header = T, fill = T)
annot_1174D <- read.csv(paste0(navin.dir, "1174.D_new.csv"), header = T, fill = T)
annot_1179B <- read.csv(paste0(navin.dir, "1179.B_new.csv"), header = T, fill = T)
annot_BW09A <- read.csv(paste0(navin.dir, "BW09.A_new.csv"), header = T, fill = T)
annot_BW11C <- read.csv(paste0(navin.dir, "BW11.C_new.csv"), header = T, fill = T)
annot_BW11D <- read.csv(paste0(navin.dir, "BW11.D_new.csv"), header = T, fill = T)
annot_BW14A <- read.csv(paste0(navin.dir, "BW14.A_new.csv"), header = T, fill = T)
annot_BW23B <- read.csv(paste0(navin.dir, "BW23.B_new.csv"), header = T, fill = T)

annot_list <- list(annot_1174C, annot_1174D, annot_1179B, annot_BW09A, annot_BW11C, annot_BW11D, annot_BW14A, annot_BW23B)

for (i in 1:length(annot_list)) {
  colnames(annot_list[[i]]) <- c("Barcode", "Navin_annot")
}

spatial_list_subset <- spatial_list[c("1174.C", "1174.D", "1179.B", "BW09.A", "BW11.C", "BW11.D", "BW14.A", "BW23.B")]

for (i in 1:length(annot_list)) {
  print(all(annot_list[[i]]$Barcode == rownames(spatial_list_subset[[i]]@meta.data)))
  spatial_list_subset[[i]]$Navin_annot <- annot_list[[i]]$Navin_annot
}

total_spots <- list()

for (i in 1:length(spatial_list_subset)) {
    meta <- spatial_list_subset[[i]]@meta.data
    obj <- meta[,c(14:20, 22:24, 78, ncol(meta))]
    caf_indices <- grep("^CAF", names(meta))
    obj$Fibroblast <- rowSums(meta[, caf_indices], na.rm = TRUE)
    total_spots[[i]] <- obj
}

library(data.table)
df <- rbindlist(total_spots)
df <- df[!(df$Navin_annot %in% ""),]

df <- df %>%
  mutate(
    Navin_annot2 = case_when(
      Navin_annot %in% c("Adenocarcinoma_lepidic", 
                         "Mucinous adenocarcinoma_lepidic", 
                         "Mucinous adenocarcinoma_micropapillary", 
                         "Adenocarcinoma_acinar/papillary", 
                         "Adenocarcinoma_micropapillary", 
                         "Adenocarcinoma_acinar", 
                         "Adenocarcinoma_solid") ~ "Malignant",
      Navin_annot %in% c("Fibroblasts", 
                         "Fibrosis", 
                         "Smooth muscle", 
                         "Vascular smooth muscle") ~ "Mesenchymal",
      Navin_annot %in% c("Immune cells", 
                         "Immune cells_lymphocytes", 
                         "Immune_lymphocytes", 
                         "Immune_TLS_BALT", 
                         "Immune cells_granulocytes", 
                         "Immune cells_BALT/TLS") ~ "Immune",
      Navin_annot %in% c("Non-malignant cells", 
                         "Uninvolved lung parenchyma", 
                         "Non-neoplastic lung") ~ "Nonmalignant epithelial",
      Navin_annot %in% c("Vasculature", 
                         "Endothelial cell", 
                         "Endothelium") ~ "Endothelial",
      TRUE ~ NA_character_  # Assign NA to any categories not covered
    )
  )

df$Navin_annot <- NULL

df.list <- split(df, by = "Navin_annot2")

for (i in 1:length(df.list)) {
  means <- colMeans(df.list[[i]][,1:(ncol(df.list[[i]]) - 1)])
  means <- as.data.frame(means)
  colnames(means) <- names(df.list)[i]
  df.list[[i]] <- as.data.frame(t(means))
}

dim(df.list[[1]])
head(df.list[[1]])

table <- rbindlist(df.list)
table <- as.data.frame(table)
rownames(table) <- names(df.list)
table <- as.matrix(table)
table <- prop.table(table,1)
rowSums(table)
table <- t(table)
Dout <- table

row_order = c("Mast", "Plasma", "Bcell", "Pericyte", "SmoothMuscle", "Fibroblast", "Endothelial", "NK", "Lymphoid", "Myeloid", "Nonmalig", "Cancer")
Dout <- Dout[row_order,]

colors <- c("#F44336", "#FF5722", "#FF9800", "#FFEB3B", "#CDDC39", "#8BC34A", "#4CAF50", "#00BCD4", "#2196F3", "#3F51B5", "#FFC107", "#283593")

pdf(paste0(figures.dir, 'EX_FIG_1H.pdf'), 6, 5.5, useDingbats=FALSE)
par(mar=c(17, 4.1, 4.1, 2.1))
barplot2(Dout, main="Composition", legend = rownames(Dout), col = colors, xlim=c(0, ncol(Dout) + 15), las=2,
        legend.text=TRUE,args.legend=list(x=ncol(Dout)+15,y=max(colSums(Dout)),bty = "n"))
dev.off()
