######################@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ 
######################@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ FIGURE 3 - EXTENDED DATA 6 - ENDOTHELIAL PANELS
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

if (!file.exists(paste0(figures.dir,'Fig3_stroma/'))){dir.create(paste0(figures.dir,'Fig3_stroma/'), recursive=TRUE)}

######
source(paste0(data.path, "/R_utils/plotutils.R"))
source(paste0(data.path, "/R_utils/seuratutils.R"))
source(paste0(data.path, "/R_utils/seuratutilsV3.R"))
source(paste0(data.path, "/R_utils/color.R"))



#####################@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ 
#####################@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
##################### Figure 3A
end <- get(load("./dropbox_data/end.final.nodoublets.Rda"))

###########
pdf(paste0(figures.dir, "Fig3_stroma/Fig3A.pdf"), useDingbats = F, width = 10)
p = DimPlot(end, group.by = "subtype")
p
dev.off()



#####################@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ 
#####################@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
##################### Figure 3B
end <- get(load("./dropbox_data/end.final.nodoublets.Rda"))

comp.df <- as.matrix(table(end$orig.identSec, end$subtype))
comp.df  <- prop.table(comp.df , 1)

comp.df <- comp.df[sampleids.filt,]

###
p53_status <- p53.status.filt
comp.df <- cbind(comp.df, p53_status)
comp.df = as.data.frame(comp.df)
comp.df$p53_status <- factor(comp.df$p53_status, levels = c("WT", "p53_mut"))

###
subtypes <- colnames(comp.df)[-10]
for (i in 1:(ncol(comp.df)-1))
{
  comp.df[,i] = as.numeric(comp.df[,i])
}

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

### 3B plots only
mean.plots.2 = list()
mean.plots.2[[1]] = mean.plots[[1]] ## Aerocyte
mean.plots.2[[2]] = mean.plots[[2]]   ## Arterial

pdf(paste0(figures.dir, "Fig3_stroma/Fig3B.pdf"), width = 6, height = 3, useDingbats=FALSE)
cowplot::plot_grid(plotlist = mean.plots.2, ncol=2)
dev.off()



#####################@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ 
#####################@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
##################### 3C
df <- get(load(paste0(data.path, "lig_rec_differential_highlvl_df_not_allzeros.Rda")))

df <- df[grepl("Endothelial", df$cell.pair, fixed=TRUE),]

df2 <- df[!grepl("Nonmalig", df$cell.pair, fixed = TRUE),]

#end
int_pairs <- c("VEGFA:KDR", "VEGFA:FLT1", "TNFSF10:TNFRSF10D", "TNFSF10:RIPK1", "TNFRSF10B:TNFSF10", "TNFRSF10A:TNFSF10",
               "TGFB2:TGFBR3", "SPP1:a9b1 complex", "PlexinA2_complex1:SEMA3A", "PODXL:SELL", "NRP2:VEGFA", "NRP2:SEMA3C", "NRP2:PGF",
               "NRP1:VEGFA", "NRP1:SEMA3A", "NRP1:PGF", "NOTCH1:JAG1", "NCR3:BAG6", "NAMPT:P2RY6", "LGALS9:PTPRK", "LGALS9:HAVCR2",
               "JAG1:NOTCH4", "FN1:a8b1 complex", "FLT1:PGF", "FLT1 complex:VEGFA", "FLT1 complex:PGF", "EPHB2:EFNB2",
               "EPHB2:EFNB1", "EFNB2:EPHB3", "EFNA1:EPHA1", "COPA:P2RY6", "COL6A2:a1b1 complex", "COL5A3:a10b1 complex", "COL16A1:a2b1 complex",
               "COL16A1:a10b1 complex", "COL15A1:a10b1 complex", "COL12A1:a10b1 complex", "CD55:ADGRE5", "CD46:JAG1", "CCL5:ACKR1",
               "ADORA2B:ENTPD1", "ACKR1:CCL17")
df2 <- df2[df2$int_name %in% int_pairs,]

range <- as.numeric(max(abs(df2$p53_wt_prop_diff)))

colors <- c(brewer_pal(palette = "Spectral", direction = -1)(7))
pal <- gradient_n_pal(colors)
custom_color_scale <- scale_fill_gradientn(
  colours = pal(c(0, rescale(seq_along(df2$p53_wt_prop_diff)), 1)), 
  values = c(0,rescale(seq_along(df2$p53_wt_prop_diff)),1), 
  limits = c(-range,range)
) 

############
pdf(file = paste0(figures.dir, "Fig3_stroma/Fig3C.pdf"), useDingbats = F, height = 15, width = 8.3)
p<-ggplot(data = df2, mapping = aes(x=cell.pair, y=int_name, color=p53_wt_prop_diff, size=p53_wt_neglogpval_proptest))+geom_point(shape = 21, aes(colour = as.factor(signif), fill = p53_wt_prop_diff))+
  scale_colour_manual(values=c("00FFFFFF", "darkgray", "black")) + theme_minimal() +
  custom_color_scale +
  theme(text = element_text(size=12), strip.text = element_text(size=22,face='bold'))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
p
dev.off()



#####################@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ 
#####################@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
##################### Extended 6B
end <- get(load("./dropbox_data/end.final.nodoublets.Rda"))

genes <- c("HPGD", "EDNRB", "S100A4", "IGFBP3", "DKK2", "FBLN5", "FCN3", "CD36", "IL7R", "CCL21", "TFF3", "FABP4", "CPE", "C7", "CLU", "SPRY1", "ZNF385D",
           "OLFM1", "COL4A1", "PLVAP", "ESM1", "STMN1", "TOP2A", "MKI67", "ISG15", "CXCL10", "CXCL11")   

DefaultAssay(end) <- "RNA"
end <- NormalizeData(end)
end <- FindVariableFeatures(end)
end <- ScaleData(end, features = unique(VariableFeatures(end), genes))

names<-c("Aerocyte", "Arterial", "Capillary", "Lymphatic", "Pulmonary.Venous", "Systemic.Venous", "VEC.COL4A1", "VEC.Cycling", "VEC.IFN")

Idents(end) <- 'subtype'
Idents(end) <- factor(Idents(end), levels = rev(names)) 

###########
pdf(paste0(figures.dir, "Fig3_stroma/Ext_Fig6B.pdf"), useDingbats = F, width = 14, height = 5)
p = DotPlot(object = end, features = genes, scale = T) +
  theme(axis.text.x = element_text(angle = 45, hjust=1), panel.border = element_rect(colour = "black", fill=NA, size=0.5), panel.grid.major = element_line(colour = "gainsboro")) + 
  scale_colour_gradient2(low = "navy", high = "firebrick") +
  geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.5)
p
dev.off()



#####################@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ 
#####################@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
##################### Extended 6C
end <- get(load("./dropbox_data/end.final.nodoublets.Rda"))

DefaultAssay(end) <- "RNA"
end <- NormalizeData(end)

###
pal <- gradient_n_pal(brewer_pal(palette = "Spectral", direction = -1)(7))

genes <- c("PECAM1", "ACKR1", "GJA5", "PLVAP", "CCL21", "MKI67", "EDNRB", "CA4", "COL4A1")

set.seed(1)
data <- end[,sample(colnames(end), 500)][["RNA"]]@data[genes,]
data <- as.matrix(data)

p <- FeaturePlot(end, features = genes, combine = FALSE, pt.size = 0.01)

custom_color_scale <- scale_color_gradientn(
  colours = pal(c(0, rescale(seq_along(data)), 1)), 
  limits = c(0, 4), breaks = 0:4,
  na.value = "#9E0142",
  values = c(0,rescale(seq_along(data)),1) 
) 

p <- lapply(p, function (x) x + custom_color_scale)
for(i in 1:length(p)) {p[[i]] <- p[[i]] + NoAxes()}

pdf(paste0(figures.dir, "Fig3_stroma/Ext_Fig6C.pdf"), width = 15, height = 10, useDingbats = F)
cowplot::plot_grid(plotlist = p, ncol=3)
dev.off()



#####################@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ 
#####################@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
##################### Extended 6D
end <- get(load("./dropbox_data/end.final.nodoublets.Rda"))

# Stacked bar plot
obj <- end

Idents(obj) <- obj$orig.identSec

my_color_palette <- rev(hue_pal()(length(unique(obj$subtype))))
Dout <- table(obj$subtype,Idents(obj))

Dout <- Dout[,c(WT.filt,p53.filt, other)]
Dout <- prop.table(Dout, 2)

names<-c("Aerocyte", "Arterial", "Capillary", "Lymphatic", "Pulmonary.Venous", "Systemic.Venous", "VEC.COL4A1", "VEC.Cycling", "VEC.IFN")
Dout <- Dout[rev(names),]

############
pdf(paste0(figures.dir, "Fig3_stroma/Ext_Fig6D.pdf"), height = 5, width = 14, useDingbats=FALSE)
par(mar=c(5, 4.1, 4.1, 2.1))
barplot2(Dout, main="Cell type composition", legend = rownames(Dout),col = my_color_palette, xlim=c(0, ncol(Dout) + 15), las=2,
         legend.text=TRUE,args.legend=list(x=ncol(Dout)+15,y=max(colSums(Dout)),bty = "n"))
dev.off()



#####################@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ 
#####################@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
##################### Extended 6F_top
end <- get(load("./dropbox_data/end.final.nodoublets.Rda"))

genes <- c('PGF', 'FN1','ITGA9','TNFRSF10D','HAVCR2','COL6A2','TNFRSF10B','NAMPT','NRP2','ITGB1','BAG6','NRP1','NOTCH1','COPA','CD55','FLT1','PODXL',
           'ENTPD1','EFNA1','KDR','TGFBR3','ITGA2','PLXNA2','ITGA10','TNFRSF10A','EFNB2','EFNB1','LGALS9','CD46','NOTCH4','ACKR1')

DefaultAssay(end) <- "RNA"
end <- NormalizeData(end)
end <- FindVariableFeatures(end)
end <- ScaleData(end, features = unique(VariableFeatures(end), genes))

names<-c("VEC.Cycling", "VEC.IFN","Aerocyte", "VEC.COL4A1","Capillary", "Pulmonary.Venous", "Systemic.Venous", "Arterial", "Lymphatic")

Idents(end) <- 'subtype'
Idents(end) <- factor(Idents(end), levels = rev(names)) 

###########
pdf(paste0(figures.dir, "Fig3_stroma/Ext_Fig6F_top.pdf"), useDingbats = F, width = 14, height = 5)
p = DotPlot(object = end, features = genes, scale = T) +
  theme(axis.text.x = element_text(angle = 45, hjust=1), panel.border = element_rect(colour = "black", fill=NA, size=0.5), panel.grid.major = element_line(colour = "gainsboro")) + 
  scale_colour_gradient2(low = "navy", high = "firebrick") +
  geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.5)
p
dev.off()



#####################@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ 
#####################@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
##################### Extended 6F_bottom
end <- get(load("./dropbox_data/end.final.nodoublets.Rda"))

genes <- c('PGF', 'TNFSF10', 'FN1','ITGA9','TNFRSF10D','HAVCR2','COL6A2','TNFRSF10B','NAMPT','NRP2','ITGB1','BAG6','NRP1','NOTCH1','COPA','CD55','FLT1','PODXL',
           'ENTPD1','EFNA1','KDR','TGFBR3','ITGA2','PLXNA2','ITGA10','TNFRSF10A','EFNB2','EFNB1','LGALS9','CD46','NOTCH4','ACKR1')
obj = end

####
obj.subtype.list <- SplitObject(obj, split.by = "subtype")
avglog2fc.list <- list()
pval.list <- list()

for (i in 1:length(obj.subtype.list)) {
  subtype <- obj.subtype.list[[i]]
  
  Idents(subtype) <- "orig.identSec"
  subtype <- AverageExpression(subtype, slots = "counts", assays = "RNA")[[1]]
  subtype <- subtype[,c(p53.filt[p53.filt %in% colnames(subtype)], WT.filt[WT.filt %in% colnames(subtype)])]
  
  subtype <- CreateSeuratObject(subtype)
  
  num.p53 <- length(p53.filt[p53.filt %in% colnames(subtype)])
  num.wt <- length(WT.filt[WT.filt %in% colnames(subtype)])
  
  subtype$p53_status <- c(rep("p53_mut", num.p53), rep("WT", num.wt))
  subtype$p53_status <- factor(subtype$p53_status, levels = c("WT", "p53_mut"))
  
  Idents(subtype) <- "p53_status"
  subtype <- NormalizeData(subtype) 
  subtype[["RNA"]]@counts <- subtype[["RNA"]]@counts[genes[genes %in% rownames(subtype[["RNA"]]@counts)],]
  subtype[["RNA"]]@data <- subtype[["RNA"]]@data[genes[genes %in% rownames(subtype[["RNA"]]@data)],]
  
  markers <- FindMarkers(subtype, ident.1 = "p53_mut", ident.2 = "WT", logfc.threshold = 0, min.pct = 0, min.cells.feature = 0, min.cells.group = 0)
  
  markers <- markers[genes,]
  print(head(markers))
  
  avglog2fc.list[[i]] <- as.data.frame(t(as.data.frame(markers$avg_log2FC)))
  pval.list[[i]] <- as.data.frame(t(as.data.frame(markers$p_val)))
  
}

##########
avglog2fc.df <- as.data.frame(rbindlist(avglog2fc.list))
pvals.df <- as.data.frame(rbindlist(pval.list))
dim(avglog2fc.df)
dim(pvals.df)
rownames(avglog2fc.df) <- names(obj.subtype.list)
rownames(pvals.df) <- names(obj.subtype.list)
colnames(avglog2fc.df) <- genes
colnames(pvals.df) <- genes

pvals.df <- as.matrix(pvals.df)
pvals.df[is.nan(pvals.df)] <- NA
pvals.df <- as.data.frame(pvals.df)
pvals.df[pvals.df <= 0.05] <- "**"
pvals.df[pvals.df > 0.05 & pvals.df <= 0.1] <- "*"

pvals.df[is.na(pvals.df)] <- ""
pvals.df[pvals.df > 0.1] <- ""

range <- 1.5

avglog2fc.df[,names(which(sapply(avglog2fc.df, function(x)all(is.na(x)))))] <- NULL # remove cols with all NAs
pvals.df <- pvals.df[,colnames(avglog2fc.df)] # remove cols with all NAs

names<-c("VEC.Cycling", "VEC.IFN","Aerocyte", "VEC.COL4A1","Capillary", "Pulmonary.Venous", "Systemic.Venous", "Arterial", "Lymphatic")
avglog2fc.df = avglog2fc.df[names,]
pvals.df = pvals.df[names,]

##########
pdf(paste0(figures.dir, "Fig3_stroma/Ext_Fig6F_bottom.pdf"), useDingbats = F, width = 20, height = 8) 
pheatmap(avglog2fc.df, display_numbers = pvals.df, breaks = seq(-range, range, length.out = 100),  cellheight=13,cellwidth=13, 
         cluster_rows = F,cluster_cols = F,
         color = colorRampPalette(c("navy", "white", "firebrick"))(100))
dev.off()