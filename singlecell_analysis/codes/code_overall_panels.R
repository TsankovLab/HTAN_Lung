######################@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ 
######################@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ FIGURE 1 - EXTENDED DATA 1 - OVERALL PANELS
library(ggplot2)
library(Seurat)
library(SeuratObject)
library(dplyr)
library(Matrix)
library(RColorBrewer)
library(scales)
library(gplots)
library(ComplexHeatmap)
library(data.table)
library(pheatmap)
library(rlist)
library(ggpubr)
library(circlize)
library(readxl)

######
# filtered
p53.filt <- c("1172", "1173", "1175", "1182", "BW16", "BW11", "BW14", "BW06")
WT.filt <- c("14", "1174", "1176", "1183", "BW01", "BW04", "BW05", "BW09", "BW19", "BW23") 

sampleids.filt <- c(p53.filt, WT.filt)
p53.status.filt <- c(rep("p53_mut", 8), rep("WT", 10))

other <- c("19", "1181", "1179", "1170", "17")
total.samples <- c(p53.filt, WT.filt, other)

# my outputs
data.path <- "../data/"
figures.dir <- "../../figures/"

if (!file.exists(paste0(figures.dir,'Fig1_overview/'))){dir.create(paste0(figures.dir,'Fig1_overview/'), recursive=TRUE)}

######
source(paste0(data.path, "/R_utils/plotutils.R"))
source(paste0(data.path, "/R_utils/seuratutils.R"))
source(paste0(data.path, "/R_utils/seuratutilsV3.R"))
source(paste0(data.path, "/R_utils/color.R"))



#####################@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ 
#####################@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
##################### Figure 1C

##########
# load all obj
all = get(load("../zenodo_data/all.merge.Rda"))

# celltype_2: rename celltype
all$celltype_2 = all$celltype
all$celltype_2[which(all$celltype == "Cancer")] = "Malignant"
all$celltype_2[which(all$celltype == "Nonmalig")] = "Nonmalig Epi"
all$celltype_2[which(all$subtype == "Lymphatic")] = "Lymphatic"
all$celltype_2[which(all$subtype %in% c("CAF.ADH1B", "Myofibroblast", "CAF.Complement", 
                                        "CAF.COL", "CAF.Ribo", "CAF.Adventitial", 
                                        "CAF.APOE.IGFBP", "CAF.ISG", "CAF.Cycling"))] = "CAF"
all$celltype_2[which(all$subtype %in% c("Pericyte.EMT", "SmoothMuscle.Airway", "SmoothMuscle.Vascular"))] = "SmoothMuscle"

#
genes <- c("EPCAM", "KRT19", "KRT8", "SLPI", "FXYD3", "CYB5A", "LYZ", "TYROBP", "FCER1G", "TRAC", "CD3D", 
           "IL32", "GNLY", "NKG7", "GZMB", "VWF", "CLDN5", "PECAM1", "CCL21", "FABP4", "TFF3", 
           "LUM", "COL3A1", "COL1A1", "ACTA2", "TAGLN", "MYH11", "HIGD1B", "COX4I2", "PTN", 
           "CD79A", "MS4A1", "CD79B", "IGLC2", "IGKC", "IGHA1", "TPSB2", "TPSAB1", "CPA3")

all <- NormalizeData(all)
all <- FindVariableFeatures(all)
all <- ScaleData(all, features = unique(c(VariableFeatures(all), genes))) 


Idents(all) <- "celltype_2"
Idents(all) <- factor(Idents(all), levels = rev(c("Malignant", "Nonmalig Epi", "Myeloid", "Tcell", "NK", "Endothelial", "Lymphatic",
                                                  "CAF", "SmoothMuscle", "Pericyte", "Bcell", "Plasma", "Mast")))

#
set.seed(1)
data <- all[,sample(colnames(all), 500)][["RNA"]]@data[genes,]
data <- as.matrix(data)

pal <- gradient_n_pal(brewer_pal(palette = "Spectral", direction = -1)(7))

custom_color_scale <- scale_color_gradientn(
  colours = pal(c(0, rescale(seq_along(data)), 1)), 
  limits = c(0, 4), breaks = 0:4,
  na.value = "#9E0142",
  values = c(0,rescale(seq_along(data)),1) 
)


############ figure 1C
pdf(paste0(figures.dir, "Fig1_overview/Fig1C.pdf"), useDingbats = F, width = 12, height = 4)
p = DotPlot(object = all, features = genes, scale = T) +
  theme(axis.text.x = element_text(angle = 45, hjust=1), panel.border = element_rect(colour = "black", fill=NA, size=0.5), panel.grid.major = element_line(colour = "gainsboro")) + 
  scale_colour_gradient2(low = "navy", high = "firebrick") +
  geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.5)
p
dev.off()



#####################@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ 
#####################@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
##################### Figure 1D
all = get(load("../zenodo_data/all.merge.Rda"))

# celltype_2: rename celltype
all$celltype_2 = all$celltype
all$celltype_2[which(all$celltype == "Cancer")] = "Malignant"
all$celltype_2[which(all$celltype == "Nonmalig")] = "Nonmalig Epi"
all$celltype_2[which(all$subtype == "Lymphatic")] = "Lymphatic"
all$celltype_2[which(all$subtype %in% c("CAF.ADH1B", "Myofibroblast", "CAF.Complement", 
                                        "CAF.COL", "CAF.Ribo", "CAF.Adventitial", 
                                        "CAF.APOE.IGFBP", "CAF.ISG", "CAF.Cycling"))] = "CAF"
all$celltype_2[which(all$subtype %in% c("Pericyte.EMT", "SmoothMuscle.Airway", "SmoothMuscle.Vascular"))] = "SmoothMuscle"

############ 1D right
my_color_palette<-c("#FF8D05", "#ACCF41", "#32B37C", "#70BA4C", "#283593", "#F44336",
                    "#3D4EB2", "#08B2DB", "#FFB605", "#FFDD29", "#FE5523", "#E1E239", "#268AE8")
pdf(paste0(figures.dir, "Fig1_overview/Fig1D.left.pdf"), useDingbats = F, width = 8)
p<-DimPlot(all, group.by = "celltype_2", cols = my_color_palette)
p
dev.off()

############ 1D left
all$p53_status <- "WT"
all$p53_status[all$orig.identSec %in% p53.filt] <- "p53_mut"
all$p53_status[all$orig.identSec %in% other] <- "other"

Idents(all) <- "p53_status"
Idents(all) <- factor(Idents(all), levels = c("WT", "p53_mut", "other"))
colors <- c("#0091CA", "#D8423D", "gray")
pdf(file = paste0(figures.dir, "Fig1_overview/Fig1D.right.pdf"), useDingbats = F, width = 8)
DimPlot(all, cols = colors)
dev.off()



#####################@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ 
#####################@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
##################### Figure 1E
malig <- get(load("../zenodo_data/malig.Rda"))

maligv9.allsamples = malig

maligv9.scores.df <- maligv9.allsamples@meta.data[,c("orig.identSec", colnames(maligv9.allsamples@meta.data)[52:70])]

maligv9.scores.list <- split(maligv9.scores.df, f = maligv9.scores.df$orig.identSec)
maligv9.scores.list.avg <- lapply(maligv9.scores.list, function(x) {
  y <- colMeans(x[,2:ncol(x)])
  y
})
names(maligv9.scores.list.avg) <- names(maligv9.scores.list)

###
maligv9.scores.avg.rbind <- list.rbind(maligv9.scores.list.avg)
maligv9.scores.avg.rbind <- as.data.frame(maligv9.scores.avg.rbind)
maligv9.scores.avg.rbind <- maligv9.scores.avg.rbind[c(p53.filt, WT.filt),]

maligv9.scores.avg.rbind$p53_status <- p53.status.filt
maligv9.scores.avg.rbind$p53_status <- factor(maligv9.scores.avg.rbind$p53_status, levels = c("WT", "p53_mut"))

###
comp.df = maligv9.scores.avg.rbind
colnames(comp.df) = gsub("1$", "", colnames(comp.df))

subtypes <- colnames(comp.df)[-20]
for (i in 1:(ncol(comp.df)-1))
{
  comp.df[,i] = as.numeric(comp.df[,i])
}

## update subtype names
subtypes[1] = "AT2.like"
subtypes[2] = "AT1.2.like"
subtypes[7] = "MHCII"
subtypes[10] = "StressResponse"
subtypes[11] = "Respiration.MT"

colnames(comp.df)[-20] = subtypes

### all plots
mean.plots <- list()

for (subtype in subtypes){
  mean.plot <- ggboxplot(comp.df, x = "p53_status", y = subtype, color = "p53_status",
                         palette = c("#0091CA", "#D8423D"), add = "jitter",
                         title = subtype) +
    stat_compare_means(label.x = 0.9, label.y = .9 * max(as.numeric(comp.df[,subtype]))) +
    theme(plot.title = element_text(size = 12, face = "bold")) + NoLegend()
  mean.plots[[subtype]] <- mean.plot
}

### 1E plots only
mean.plots.2 = mean.plots[[19]] ## tp53 targets

pdf(paste0(figures.dir, "Fig1_overview/Fig1E.pdf"), width = 3, height = 4, useDingbats=FALSE)
mean.plots.2
dev.off()



#####################@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ 
#####################@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
##################### Figure 1F and Extended 1D
all = get(load("../zenodo_data/all.merge.Rda"))

# celltype_2: rename celltype
all$celltype_2 = all$celltype
all$celltype_2[which(all$celltype == "Cancer")] = "Malignant"
all$celltype_2[which(all$celltype == "Nonmalig")] = "Nonmalig Epi"
all$celltype_2[which(all$subtype == "Lymphatic")] = "Lymphatic"
all$celltype_2[which(all$subtype %in% c("CAF.ADH1B", "Myofibroblast", "CAF.Complement", 
                                        "CAF.COL", "CAF.Ribo", "CAF.Adventitial", 
                                        "CAF.APOE.IGFBP", "CAF.ISG", "CAF.Cycling"))] = "CAF"
all$celltype_2[which(all$subtype %in% c("Pericyte.EMT", "SmoothMuscle.Airway", "SmoothMuscle.Vascular"))] = "SmoothMuscle"

all.v2.meta = all@meta.data

CD4T_ID = rownames(all.v2.meta)[which(all.v2.meta$subtype %in% c("CD4.Naive.CM", "CD4.TRM","Treg", "TFH"))]
CD8T_ID = rownames(all.v2.meta)[which(all.v2.meta$subtype %in% c("CD8.IFN", "CD8.GZMK", "CD8.TRM", "T.Exhausted"))]

#
all$celltype_3 = all$celltype_2
all$celltype_3[which(all$subtype %in% c("CD4.Naive.CM", "CD4.TRM","Treg", "TFH"))] = "CD4_T"
all$celltype_3[which(all$subtype %in% c("CD8.IFN", "CD8.GZMK", "CD8.TRM", "T.Exhausted"))] = "CD8_T"
all$celltype_3[which(all$celltype_3 == "Tcell")] = "T_other"


############ Figure 1F and Extended 1D left part
Idents(all) <- "orig.ident"
unsorted <- subset(all, idents = c("14MA","14ML","17M","19_Norm","MGH_1170","MGH1172a","MGH1172b",
                                      "MGH1173a","MGH1173b","MGH1174a","MGH1174b","MGH1175a","MGH1175b","MGH1176","MGH1179a","MGH1181",
                                      "MGH1182","MGH1183","NSC001BLT_FT","NSC004bBN_FT","NSC005bLT_FT","NSC006bLT_FT_V2","NSC009bLT_FT_T1",
                                      "NSC011bLT_FT_T1","NSC014bLT_FT1","NSC016bLT_FT_T1","NSC019bLTcd45n","NSC023bLT")) 

comp.df <- as.matrix(table(unsorted$orig.identSec, unsorted$celltype_2))
comp.df  <- prop.table(comp.df , 1)

comp.df <- comp.df[sampleids.filt,]

############
p53_status <- p53.status.filt
comp.df <- cbind(comp.df, p53_status)
comp.df = as.data.frame(comp.df)
comp.df$p53_status <- factor(comp.df$p53_status, levels = c("WT", "p53_mut"))


subtypes <- colnames(comp.df)[-14]
for (i in 1:(ncol(comp.df)-1))
{
  comp.df[,i] = as.numeric(comp.df[,i])
}

mean.plots <- list()

for (subtype in subtypes){
  mean.plot <- ggboxplot(comp.df, x = "p53_status", y = subtype, color = "p53_status",
                         palette = c("#0091CA", "#D8423D"), add = "jitter",
                         title = subtype) +
    stat_compare_means(label.x = 0.9, label.y = .9 * max(as.numeric(comp.df[,subtype]))) +
    theme(plot.title = element_text(size = 12, face = "bold")) + NoLegend()
  mean.plots[[subtype]] <- mean.plot
}


############ Figure 1F
mean.plots.2 = list()
mean.plots.2[[1]] = mean.plots[[3]] ## Endothelial
mean.plots.2[[2]] = mean.plots[[10]]   ## Pericyte

pdf(paste0(figures.dir, "Fig1_overview/Fig1F.pdf"), width = 6, height = 3, useDingbats=FALSE)
cowplot::plot_grid(plotlist = mean.plots.2, ncol=2)
dev.off()


############ Extended 1D left
mean.plots.2 = list()
mean.plots.2[[1]] = mean.plots[[1]] 
mean.plots.2[[2]] = mean.plots[[2]]   
mean.plots.2[[3]] = mean.plots[[5]] 
mean.plots.2[[4]] = mean.plots[[6]] 
mean.plots.2[[5]] = mean.plots[[7]] 
mean.plots.2[[6]] = mean.plots[[8]] 
mean.plots.2[[7]] = mean.plots[[12]] 
mean.plots.2[[8]] = mean.plots[[13]] 

pdf(paste0(figures.dir, "Fig1_overview/Ext_Fig1D.pdf"), width = 12, height = 6, useDingbats=FALSE)
cowplot::plot_grid(plotlist = mean.plots.2, ncol=4)
dev.off()

comp.df1 = comp.df

############ Extended 1D right (CD8 and CD4)
comp.df <- as.matrix(table(unsorted$orig.identSec, unsorted$celltype_3))  ## has t cell annotations
comp.df  <- prop.table(comp.df , 1)

comp.df <- comp.df[sampleids.filt,]

############
p53_status <- p53.status.filt
comp.df <- cbind(comp.df, p53_status)
comp.df = as.data.frame(comp.df)
comp.df$p53_status <- factor(comp.df$p53_status, levels = c("WT", "p53_mut"))


subtypes <- colnames(comp.df)[-16]
for (i in 1:(ncol(comp.df)-1))
{
  comp.df[,i] = as.numeric(comp.df[,i])
}

############
mean.plots <- list()

for (subtype in c("CD4_T", "CD8_T")){
  mean.plot <- ggboxplot(comp.df, x = "p53_status", y = subtype, color = "p53_status",
                         palette = c("#0091CA", "#D8423D"), add = "jitter",
                         title = subtype) +
    stat_compare_means(label.x = 0.9, label.y = .9 * max(as.numeric(comp.df[,subtype]))) +
    theme(plot.title = element_text(size = 12, face = "bold")) + NoLegend()
  mean.plots[[subtype]] <- mean.plot
}

pdf(paste0(figures.dir, "Fig1_overview/Ext_Fig1D_part2.pdf"), width = 3, height = 6, useDingbats=FALSE)
cowplot::plot_grid(plotlist = mean.plots, ncol=1)
dev.off()

comp.df2 = comp.df



#####################@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ 
#####################@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
##################### Extended 1A
info = read.csv(paste0(data.path, "CNV_corr.csv"))   
rownames(info) = info[,1]
info = info[,-1]

############
pdf(paste0(figures.dir,"Fig1_overview/Ext_Fig1A.pdf"), width = 8, height = 8, useDingbats=F)
p<-pheatmap(info, cluster_rows = F, cluster_cols = F,
            breaks = seq(min(info), max(info), length.out = 100))
p
dev.off()

#####################@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ 
#####################@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
##################### Extended 1B
all = get(load("../zenodo_data/all.merge.Rda"))

############
pdf(file = paste0(figures.dir, "Fig1_overview/Ext_Fig1B.pdf"), useDingbats = F, width = 8)
DimPlot(all, group.by = "orig.identSec")
dev.off()



#####################@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ 
#####################@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
##################### Extended 1C
all = get(load("../zenodo_data/all.merge.Rda"))

# celltype_2: rename celltype
all$celltype_2 = all$celltype
all$celltype_2[which(all$celltype == "Cancer")] = "Malignant"
all$celltype_2[which(all$celltype == "Nonmalig")] = "Nonmalig Epi"
all$celltype_2[which(all$subtype == "Lymphatic")] = "Lymphatic"
all$celltype_2[which(all$subtype %in% c("CAF.ADH1B", "Myofibroblast", "CAF.Complement", 
                                        "CAF.COL", "CAF.Ribo", "CAF.Adventitial", 
                                        "CAF.APOE.IGFBP", "CAF.ISG", "CAF.Cycling"))] = "CAF"
all$celltype_2[which(all$subtype %in% c("Pericyte.EMT", "SmoothMuscle.Airway", "SmoothMuscle.Vascular"))] = "SmoothMuscle"

# Stacked bar plot
obj <- all

Idents(obj) <- obj$orig.identSec

my_color_palette <- c("#F44336", "#FE5523", "#FF8D05", "#FFDD29", "#E1E239",
                      "#ACCF41", "#70BA4C", "#32B37C", "#08B2DB", "#268AE8", "#3D4EB2", "#FFB605",
                      "#283593")

#
Dout <- table(obj$celltype_2,Idents(obj))
Dout <- Dout[,c(WT.filt,p53.filt, other)]
Dout <- prop.table(Dout, 2)

names<-c("Malignant", "Nonmalig Epi", "Myeloid", "Tcell", "NK", "Endothelial", "Lymphatic",
         "CAF", "SmoothMuscle", "Pericyte", "Bcell", "Plasma", "Mast")
Dout <- Dout[rev(names),]
colnames(Dout) = c("P14", "MGH1174", "MGH1176", "MGH1183", "BWH01", "BWH04", "BWH05", "BWH09",
                   "BWH19", "BWH23", "MGH1172", "MGH1173", "MGH1175", "MGH1182", "BWH16", "BWH11",
                   "BWH14", "BWH06", "P19", "MGH1181", "MGH1179", "MGH1170", "P17")

############
pdf(paste0(figures.dir, "Fig1_overview/Ext_Fig1C.pdf"), height = 5, width = 10, useDingbats=FALSE)
barplot2(Dout, main="Cell type composition", legend = rownames(Dout),col = my_color_palette, xlim=c(0, ncol(Dout) + 15), las=2,
         legend.text=TRUE,args.legend=list(x=ncol(Dout)+15,y=max(colSums(Dout)),bty = "n"))
dev.off()
