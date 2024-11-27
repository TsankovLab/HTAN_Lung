######################@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ 
######################@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ FIGURE 3 - EXTENDED DATA 7 - MESENCHYMAL PANELS
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
data.path <- "./data/"
figures.dir <- "../../figures/"

if (!file.exists(paste0(figures.dir,'Fig3_stroma/'))){dir.create(paste0(figures.dir,'Fig3_stroma/'), recursive=TRUE)}

######
source(paste0(data.path, "/R_utils/plotutils.R"))
source(paste0(data.path, "/R_utils/seuratutils.R"))
source(paste0(data.path, "/R_utils/seuratutilsV3.R"))
source(paste0(data.path, "/R_utils/color.R"))



#####################@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ 
#####################@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
##################### Figure 3E
mes <- get(load("/dropbox_data/mes.final.nodoublets.Rda"))

###########
pdf(paste0(figures.dir, "Fig3_stroma/Fig3E.pdf"), useDingbats = F, width = 10)
p = DimPlot(mes, group.by = "subtype")
p
dev.off()



#####################@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ 
#####################@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
##################### Figure 3F
meta = get(load(paste0(data.path, "Fig3F_metadata_labeltrans.Rda")))

conf_mtx <- matrix(0, ncol = length(unique(meta$predicted.id)), nrow = length(unique(meta$subtype)))
colnames(conf_mtx) <- unique(meta$predicted.id)
rownames(conf_mtx) <- unique(meta$subtype)
for (actual.subtype in rownames(conf_mtx)){
  for (predicted.subtype in colnames(conf_mtx)){
    pair_freq <- sum(meta$subtype==actual.subtype & meta$predicted.id==predicted.subtype)
    conf_mtx[actual.subtype,predicted.subtype] <- pair_freq
  }
}
library(LICORS)
conf_mtx <- t(normalize(conf_mtx, byrow=TRUE))

conf_mtx2 = conf_mtx[c("Vascular Smooth Muscle","Airway Smooth Muscle","Fibromyocyte","Lipofibroblast","Alveolar Fibroblast",
                       "Adventitial Fibroblast","Pericyte","Myofibroblast"),
                     c("Pericyte", "Pericyte.EMT","SmoothMuscle.Airway","SmoothMuscle.Vascular","CAF.ADH1B","CAF.Adventitial","Myofibroblast",
                       "CAF.COLs","CAF.Complement","CAF.ISGs", "CAF.Cycling","CAF.Ribo","CAF.APOE")]

pdf(paste0(figures.dir, "Fig3_stroma/Fig3F.pdf"), height = 5, width = 7, useDingbats=FALSE)
par(mar = c(10, 10, 10, 10))
heatmap.2(conf_mtx2, tracecol=NA, margins=c(12,14), col=bluered(100), Rowv=F,Colv=F,dendrogram="none")
dev.off()



#####################@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ 
#####################@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
##################### Figure 3G
# regulon heatmap
df <- read.csv(paste0(data.path, "mes_pyscenic_subtype_zscore.csv", row.names = 1))
df2 <- dcast(df, regulon ~ subtype, value.var = "Z")
rownames(df2) <- df2$regulon
df2$regulon <- NULL
df2 <- t(df2)
df2 <- as.matrix(df2)
library(gplots)
tfs <- c("GATA5", "FOXO3", "FOXP2", "RARG", "ELK3", "ZNF433", "TWIST1", "CREB3L1", "SOX11", "SOX4", "SPIB", "E2F1", "MYBL2", "FOSL1", "IRF7", "STAT1", "ETV7", "GATA6", "NFATC4", "ZNF426", "KLF5", "TFCP2L1", "REL", "ZNF527", "YY2", "ERG", "MLX", "RBAK", "FOXA1", "NFATC2", "ZNF671", "MEF2C", "SOX15", "JUNB", "NR4A3", "FOXN2", "MYC", "TEAD4")
df2 <- df2[,tfs]

rownames(df2) = c("CAF.ADH1B", "CAF.APOE", "CAF.COLs", "CAF.Complement",
                  "CAF.Cycling", "CAF.ISGs", "CAF.Adventitial", "CAF.Ribo", "Myofibroblast",
                  "Pericyte", "Pericyte.PTN", "SmoothMuscle.Airway", "SmoothMuscle.Vascular")

############
pdf(paste0(figures.dir, "Fig3_stroma/Fig3G.pdf"), useDingbats = F, width = 10, height = 3)
p<-Heatmap(df2, cluster_rows = FALSE, cluster_columns = FALSE)
p
dev.off()



#####################@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ 
#####################@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
##################### Figure 3H
mes <- get(load("/dropbox_data/mes.final.nodoublets.Rda"))

comp.df <- as.matrix(table(mes$orig.identSec, mes$subtype))
comp.df  <- prop.table(comp.df , 1)

comp.df <- comp.df[sampleids.filt,]

###
p53_status <- p53.status.filt
comp.df <- cbind(comp.df, p53_status)
comp.df = as.data.frame(comp.df)
comp.df$p53_status <- factor(comp.df$p53_status, levels = c("WT", "p53_mut"))

###
subtypes <- colnames(comp.df)[-14]
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

### 3H plots only
mean.plots.2 = list()
mean.plots.2[[1]] = mean.plots[[10]] ## Pericyte
mean.plots.2[[2]] = mean.plots[[1]]   ## CAF.ADH1B

pdf(paste0(figures.dir, "Fig3_stroma/Fig3H.pdf"), width = 6, height = 3, useDingbats=FALSE)
cowplot::plot_grid(plotlist = mean.plots.2, ncol=2)
dev.off()




#####################@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ 
#####################@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
########## Extended 7A
mes <- get(load("/dropbox_data/mes.final.nodoublets.Rda"))

genes <- c("PTGDS", "ADH1B", "RGCC", "TIMP1", "RGS2", "APOE", "CCDC80", "C7", "C3", "COL1A2", "COL3A1", "COL1A1", "TMEM176B", "TMEM176A", "CYP1B1", "TUBA1B", "STMN1",
           "MKI67", "ISG15", "IL32", "IFI6", "RPS19", "RPS17", "RPS4Y1", "ELN", "MT1X", "ACTG2", "HIGD1B", "COX4I2", "PTN", "COL18A1", "PLXDC1", "KCNJ8", "TAGLN", "DSTN",
           "MYH11", "ACTA2", "ID4", "FILIP1L")

DefaultAssay(mes) <- "RNA"
mes <- NormalizeData(mes)
mes <- FindVariableFeatures(mes)
mes <- ScaleData(mes, features = unique(VariableFeatures(mes), genes))

names<-c("CAF.ADH1B","CAF.APOE","CAF.Adventitial","CAF.COLs","CAF.Complement","CAF.Cycling","CAF.ISGs",
         "CAF.Ribo","Myofibroblast","Pericyte","Pericyte.EMT","SmoothMuscle.Airway","SmoothMuscle.Vascular")

Idents(mes) <- 'subtype'
Idents(mes) <- factor(Idents(mes), levels = rev(names)) 

pdf(paste0(figures.dir, "Fig3_stroma/Ext_Fig7A.pdf"), useDingbats = F, width = 14, height = 5)
p = DotPlot(object = mes, features = genes, scale = T) +
  theme(axis.text.x = element_text(angle = 45, hjust=1), panel.border = element_rect(colour = "black", fill=NA, size=0.5), panel.grid.major = element_line(colour = "gainsboro")) + 
  scale_colour_gradient2(low = "navy", high = "firebrick") +
  geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.5)
p
dev.off()



#####################@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ 
#####################@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
##################### Extended 7B
mes <- get(load("/dropbox_data/mes.final.nodoublets.Rda"))

DefaultAssay(mes) <- "RNA"
mes <- NormalizeData(mes)

##
pal <- gradient_n_pal(brewer_pal(palette = "Spectral", direction = -1)(7))

genes <- c("COL1A1", "MYH11", "COX4I2", "PTN", "ISG15", "MKI67", "ADH1B", "FGFR4", "SCX", "PI16", "IGFBP6", "TNFSF13B")

set.seed(1)
data <- mes[,sample(colnames(mes), 500)][["RNA"]]@data[genes,]
data <- as.matrix(data)

p <- FeaturePlot(mes, features = genes, combine = FALSE, pt.size = 0.01)

custom_color_scale <- scale_color_gradientn(
  colours = pal(c(0, rescale(seq_along(data)), 1)), 
  limits = c(0, 4), breaks = 0:4,
  na.value = "#9E0142",
  values = c(0,rescale(seq_along(data)),1) 
) 

p <- lapply(p, function (x) x + custom_color_scale)
for(i in 1:length(p)) {p[[i]] <- p[[i]] + NoAxes()}

pdf(paste0(figures.dir, "Fig3_stroma/Ext_Fig7B.pdf"), width = 15, height = 10, useDingbats = F)
cowplot::plot_grid(plotlist = p, ncol=4)
dev.off()




#####################@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ 
#####################@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#################### Extended 7C
mes <- get(load("/dropbox_data/mes.final.nodoublets.Rda"))

# Stacked bar plot
obj <- mes

Idents(obj) <- obj$orig.identSec

my_color_palette <- rev(hue_pal()(length(unique(obj$subtype))))
Dout <- table(obj$subtype,Idents(obj))

Dout <- Dout[,c(WT.filt,p53.filt, other)]
Dout <- prop.table(Dout, 2)

names<-c("CAF.ADH1B","CAF.APOE","CAF.Adventitial","CAF.COLs","CAF.Complement","CAF.Cycling","CAF.ISGs",
         "CAF.Ribo","Myofibroblast","Pericyte","Pericyte.EMT","SmoothMuscle.Airway","SmoothMuscle.Vascular")
Dout <- Dout[rev(names),]

############
pdf(paste0(figures.dir, "Fig3_stroma/Ext_Fig7C.pdf"), height = 5, width = 14, useDingbats=FALSE)
par(mar=c(5, 4.1, 4.1, 2.1))
barplot2(Dout, main="Cell type composition", legend = rownames(Dout),col = my_color_palette, xlim=c(0, ncol(Dout) + 15), las=2,
         legend.text=TRUE,args.legend=list(x=ncol(Dout)+15,y=max(colSums(Dout)),bty = "n"))
dev.off()



#####################@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ 
#####################@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
##################### Extended Figure 7E 
meta = get(load(paste0(data.path, "meta.caf.modulescore.top20.paper.top20.ours.rna.assay.Rda"))

paper.score = meta[,c(36:43)]
our.score = meta[,c(44:56)]

res = cor(paper.score, our.score, method = "p")

##
library(ComplexHeatmap)
library(circlize)

res2 = res[c(1,2,4,6,7,8),c(1,2,3,4,5,7,8,9)]
res2 = t(res2)
rownames(res2) = gsub("Adm20.", "", rownames(res2))
res2 = res2[c("CAF.COLs", "CAF.Complement", "CAF.Adventitial", "CAF.ADH1B", "CAF.APOE", "Myofibroblast", "CAF.ISGs", "CAF.Ribo"),
            c("c2.CAF.infla", "c4.CAF.adi", "c7.CAF.PN", "c1.CAF.myo", "c8.CAF.ap", "c6.CAF.EndMT" )]

pdf(paste0(figures.dir, "Fig3_stroma/Ext_Fig7E.pdf"), useDingbats = F, width = 8, height = 8) 
Heatmap(res2, cluster_rows = F, cluster_columns = F, name = "Pearson.cor")
dev.off()



#####################@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ 
#####################@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
##################### Extended 7I
df <- get(load(paste0(data.path, "lig_rec_differential_highlvl_df_not_allzeros.Rda")))

df <- df[grepl("Mesenchymal", df$cell.pair, fixed=TRUE),]

df2 <- df[!grepl("Nonmalig", df$cell.pair, fixed = TRUE),]

#mes
int_pairs <- c("VEGFA:KDR", "VEGFA:FLT1", "VEGFA:EPHB2", "TNFRSF1A:GRN", "TNC:a4b1 complex", "TGFB3:TGFbeta receptor2",
               "TGFB2:TGFbeta receptor2", "TGFB1:aVb6 complex", "PLAUR:a4b1 complex", "NRP2:VEGFA", "NRP1:VEGFB", "NRP1:VEGFA",
               "NRP1:SEMA3A", "NRP1:PGF", "NAMPT:P2RY6", "LGALS9:PTPRK", "JAG1:NOTCH3", "JAG1:NOTCH2", "FN1:a8b1 complex",
               "FN1:a5b1 complex", "FLT1:VEGFB", "FLT1:PGF", "FLT1 complex:VEGFA", "FBN1:a5b1 complex", "COPA:P2RY6", "COL6A2:a1b1 complex",
               "COL16A1:a2b1 complex", "COL16A1:a10b1 complex", "COL15A1:a10b1 complex", "COL12A1:a10b1 complex", "CCL2:CCR10",
               "AXL:IL15RA")
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
pdf(file = paste0(figures.dir, "Fig3_stroma/Ext_Fig7I.pdf"), useDingbats = F, height = 15, width = 8.3)
p<-ggplot(data = df2, mapping = aes(x=cell.pair, y=int_name, color=p53_wt_prop_diff, size=p53_wt_neglogpval_proptest))+geom_point(shape = 21, aes(colour = as.factor(signif), fill = p53_wt_prop_diff))+
  scale_colour_manual(values=c("00FFFFFF", "darkgray", "black")) + theme_minimal() +
  custom_color_scale +
  theme(text = element_text(size=12), strip.text = element_text(size=22,face='bold'))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
p
dev.off()



#####################@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ 
#####################@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
##################### Extended 7J_top
genes <- c("CCL2", "ITGA8", "ITGB1", "TGFBR2", "PTPRK", "COPA", "TNFRSF1A", "TGFB3", "VEGFB", "ITGA1", "NRP1", "ACVR1", "COL16A1", "NOTCH2", "FBN1", "PGF",
           "NOTCH3", "AXL", "FN1", "TNC", "ITGA5", "COL15A1", "NAMPT", "VEGFA", "PLAUR", "COL12A1")

mes <- get(load("/dropbox_data/mes.final.nodoublets.Rda"))
DefaultAssay(mes) <- "RNA"
mes <- NormalizeData(mes)
mes <- FindVariableFeatures(mes)
mes <- ScaleData(mes, features = unique(VariableFeatures(mes), genes))

names<-c("CAF.ADH1B","Myofibroblast","CAF.Ribo","Pericyte","Pericyte.EMT","SmoothMuscle.Vascular",
         "CAF.ISGs","CAF.Complement","CAF.Adventitial","SmoothMuscle.Airway","CAF.APOE","CAF.COLs","CAF.Cycling")

Idents(mes) <- 'subtype'
Idents(mes) <- factor(Idents(mes), levels = rev(names)) 

pdf(paste0(figures.dir, "Fig3_stroma/Ext_Fig7J_top.pdf"), useDingbats = F, width = 14, height = 5)
p = DotPlot(object = mes, features = genes, scale = T) +
  theme(axis.text.x = element_text(angle = 45, hjust=1), panel.border = element_rect(colour = "black", fill=NA, size=0.5), panel.grid.major = element_line(colour = "gainsboro")) + 
  scale_colour_gradient2(low = "navy", high = "firebrick") +
  geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.5)
p
dev.off()



#####################@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ 
#####################@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
##################### Extended 7J_bottom
mes <- get(load("/dropbox_data/mes.final.nodoublets.Rda"))

genes <- c("CCL2", "ITGA8", "ITGB1", "TGFBR2", "PTPRK", "COPA", "TNFRSF1A", "TGFB3", "VEGFB", "ITGA1", "NRP1", "ACVR1", "COL16A1", "NOTCH2", "FBN1", "PGF",
           "NOTCH3", "AXL", "FN1", "TNC", "ITGA5", "COL15A1", "NAMPT", "VEGFA", "PLAUR", "COL12A1")

obj = mes

########## 
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

range <- 2

avglog2fc.df[,names(which(sapply(avglog2fc.df, function(x)all(is.na(x)))))] <- NULL # remove cols with all NAs
pvals.df <- pvals.df[,colnames(avglog2fc.df)] # remove cols with all NAs

names<-c("CAF.ADH1B","Myofibroblast","CAF.Ribo","Pericyte","Pericyte.EMT","SmoothMuscle.Vascular",
         "CAF.ISGs","CAF.Complement","CAF.Adventitial","SmoothMuscle.Airway","CAF.APOE","CAF.COLs","CAF.Cycling")

avglog2fc.df = avglog2fc.df[names,]
pvals.df = pvals.df[names,]

##########
pdf(paste0(figures.dir, "Fig3_stroma/Ext_Fig7J_bottom.pdf"), useDingbats = F, width = 20, height = 8) 
pheatmap(avglog2fc.df, display_numbers = pvals.df, breaks = seq(-range, range, length.out = 100),  cellheight=13,cellwidth=13, 
         cluster_rows = F,cluster_cols = F,
         color = colorRampPalette(c("navy", "white", "firebrick"))(100))
dev.off()