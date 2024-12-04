library(Seurat)
library(pals)
library(RColorBrewer)
library(ggplot2)
library(scales)
library(ggpubr)
library(dplyr)
library(Matrix)
library(SeuratObject)
library(gplots)
library(ComplexHeatmap)
library(data.table)
library(pheatmap)
library(rlist)
library(circlize)
library(readxl)
library(tidyverse)

# my outputs
data.path <- "../data/"
figures.dir <- "../../figures/"
my.dir = '/ahg/regevdata/projects/ICA_Lung/Thinh/Projects/Workspace/Processed/p53.revision.oct2024/'
#figures.dir = paste0(my.dir, 'plots/')
source.data.dir = paste0(my.dir, 'source_data/')

######
source(paste0(data.path, "/R_utils/plotutils.R"))
source(paste0(data.path, "/R_utils/seuratutils.R"))
source(paste0(data.path, "/R_utils/seuratutilsV3.R"))
source(paste0(data.path, "/R_utils/color.R"))

#Genesets
cc.genes <- readLines(con = paste0(user.path, "/genelists/regev_lab_cell_cycle_genes.txt"))
s.genes <- cc.genes[1:43]
g2m.genes <- cc.genes[44:97]

myeloid = readRDS (file.path ('./dropbox_data','myeloid_pDC_mast_integrated.rds'))#myeloid = NormalizeData (myeloid)
###

# Set sample groups
p53.filt <- c("1172", "1173", "1175", "1182", "BW16", "BW11", "BW14", "BW06")
WT.filt <- c("14", "1174", "1176", "1183", "BW01", "BW04", "BW05", "BW09", "BW19", "BW23") 

sampleids.filt <- c(p53.filt, WT.filt)
p53.status.filt <- c(rep("p53_mut", 8), rep("WT", 10))

other <- c("19", "1181", "1179", "1170", "17")
total.samples <- c(p53.filt, WT.filt, other)


#### F4A ####
DimPlot(myeloid, group.by = "subtype")
ggsave(file.path(figures.dir,'F4_A.png'))


#### F4B ####            
pal <- gradient_n_pal(brewer_pal(palette = "Spectral", direction = -1)(7))
genes <- c("RETN", "AZU1", "S100A12", "LILRB2", "CLEC9A", "FCER1A", "LAMP3", "MKI67", "APOE", "CXCL10", "CCL4", "FABP4", "FOLR2", "SPP1", "TREM2","LILRA4","TPSB2")

data <- myeloid[,sample(colnames(myeloid), 500)][["RNA"]]@data[genes,]
data <- as.matrix (data)

DefaultAssay (myeloid) = 'RNA'
p <- FeaturePlot(myeloid, features = genes, combine = FALSE, pt.size = 0.01)

custom_color_scale <- scale_color_gradientn(
    colours = pal(c(0, rescale(seq_along(data)), 1)), # <- extra 0, 1 for out-of-bounds
    limits = c(0, 4), breaks = 0:4,
    na.value = "#9E0142",
    values = c(0,rescale(seq_along(data)),1) # <- extra 0, 1 again # also try changing to 7
) ### what is different about this?

p <- lapply(p, function (x) x + custom_color_scale)
for(i in 1:length(p)) {p[[i]] <- p[[i]] + NoLegend() + NoAxes()}

cowplot::plot_grid(plotlist = p, ncol=6)
ggsave(file.path(figures.dir,'myeloid.fp.spectral.max4.v2_rep2.png'), width = 14, height = 7)

#write.csv (do.call (cbind, lapply(p, function(x) x$data)), file.path (projdir, 'F4_B.csv'))
### feature plots are not space efficient

#### F4C_S7C ####
p53 <- c("1172", "1173", "1175", "1182", "BW16", "BW11", "BW14", "BW06")
WT <- c("14", "1174", "1176", "1183", "BW01", "BW04", "BW05", "BW09", "BW19", "BW23")
p53.filt <- c("1172", "1173", "1175", "1182", "BW16", "BW11", "BW14", "BW06")
WT.filt <- c("14", "1174", "1176", "1183", "BW01", "BW04", "BW05", "BW09", "BW19", "BW23") # no 1170

comp.df <- as.matrix(table(myeloid$orig.identSec, myeloid$subtype))
p53.status.filt <- c(rep("WT", 10), rep("p53_mut", 8))# no 1170
comp.df <- comp.df[c(WT.filt, p53.filt), ]
comp.df  <- prop.table(comp.df , 1)
comp.df <- as.data.frame.matrix(comp.df)
comp.df$p53_status <- p53.status.filt
comp.df$p53_status <- factor(comp.df$p53_status, levels = c("WT", "p53_mut"))

subtypes <- colnames(comp.df)

mean.plots <- list()

for (subtype in subtypes){
  mean.plot <- ggboxplot(comp.df, x = "p53_status", y = subtype, color = "p53_status",
               palette = c("#0091CA", "#D8423D"), add = "jitter",
               title = subtype) +
               # ylim(-.5,.5) +
               stat_compare_means(label.x = 0.9, label.y = .9 * max(as.numeric(comp.df[,subtype]))) +
               # stat_compare_means(label.x = 0.9, label.y = max(module.scores[,module])*(9.0/10), method = "t.test") +
               theme(plot.title = element_text(size = 12, face = "bold")) + NoLegend()
  #ggsave(paste0(figures.dir.p53, "median.module.scores.v3.top20.png"), width = 6, height = 6, type="cairo")
  mean.plots[[subtype]] <- mean.plot
}

#figures.dir <- paste0(proj.path, "/Results/10x_All_190615/P53_paper/")
# myeloid subtype boxplots

pdf(file.path(figures.dir,'F4C_S7C.pdf'), width = 18, height = 7.5, useDingbats=FALSE)
cowplot::plot_grid(plotlist = mean.plots, ncol=8)
dev.off()

#write.csv (do.call (cbind, lapply (mean.plots, function(x) x$data)), file.path (projdir, 'F4_C_S7_C.csv'))

#### F4D SPP1 proportion ####
p53 <- c("1172", "1173", "1175", "1182", "BW16", "BW11", "BW14", "BW06")
WT <- c("14", "1174", "1176", "1183", "BW01", "BW04", "BW05", "BW09", "BW19", "BW23")
p53.filt <- c("1172", "1173", "1175", "1182", "BW16", "BW11", "BW14", "BW06")
sampleids.filt <- c(WT, p53)


genes <- c("SPP1")
for (gene in genes) {

  # if gene, uncomment line below
#   df <- data.frame(sample = factor(broad$orig.identSec, levels = sampleids.filt), subtype = broad$subtype, value = broad[["RNA"]]@counts[gene,])
  df <- data.frame(
    sample = factor(myeloid$orig.identSec, levels = sampleids.filt), 
    subtype = myeloid$subtype,
    p53_status = ifelse (myeloid$orig.identSec %in% p53.filt, 'p53_mut','WT'),
    value = myeloid[["RNA"]]@counts[gene,])
  df = aggregate (df$value, by = list(sampleID= df$sample, subtype = df$subtype, p53_status = df$p53_status), mean)
  df$x = log10(df$x +1)

subtypes = unique (df$subtype)
  mean.plots <- list()
  print(head(df))
  for (subtype in subtypes) {
    mean.plot <- ggboxplot(df[df$subtype == subtype,], x = "p53_status", y = 'x', color = "p53_status",
                 palette = c("#0091CA", "#D8423D"), add = "jitter",
                 title = subtype) +
                 # ylim(-.5,.5) +
                 stat_compare_means(label.x = 0.9, label.y = .9 * max(as.numeric(df[df$subtype == subtype,'x']), na.rm = T), na.rm = T) +
                 # stat_compare_means(label.x = 0.9, label.y = max(module.scores[,module])*(9.0/10), method = "t.test") +
                 theme(plot.title = element_text(size = 12, face = "bold")) + NoLegend()
    #ggsave(paste0(figures.dir.p53, "median.module.scores.v3.top20.png"), width = 6, height = 6, type="cairo")
    mean.plots[[subtype]] <- mean.plot
  }


  pdf(file.path(figures.dir,'F4D_right_SPP1_single_cell_p53_mut_vs_wt.pdf'), width = 15, height = 6, useDingbats=FALSE)

  print(cowplot::plot_grid(plotlist = mean.plots, ncol=6))
  dev.off()
}


#### F4D pathway enrichment ####
q_values <- c(11.60205999,
11.60205999,
10.31695296,
9.024108864,
6.514278574,
6.514278574,
5.321481621,
5.1123827,
4.247183569,
3.32330639)

names <- c('HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION',

'HALLMARK_HYPOXIA',
'HALLMARK_GLYCOLYSIS',
'HALLMARK_MTORC1_SIGNALING',
'HALLMARK_COMPLEMENT',
'HALLMARK_INFLAMMATORY_RESPONSE',
'HALLMARK_KRAS_SIGNALING_UP',
'HALLMARK_COAGULATION',
'HALLMARK_IL2_STAT5_SIGNALING',
'HALLMARK_PROTEIN_SECRETION')


mydf <- data.frame(A=names, B=q_values)

pdf (file.path(figures.dir,'F4D_enrichment.pdf'))
par(mar=c(6,16,4,4))
barplot(rev(mydf$B), main="SPP1 Correlated Genes", horiz=TRUE,
 names.arg=rev(names), las = 1, cex.names=0.6,
col = ifelse(rev(mydf$B) < 0, "#b0002c", "#2463ac"), xlab = "Log q value", cex.axis=1)
dev.off()


Idents(myeloid) <- "orig.identSec"
myeloid.p53 <- subset(myeloid, idents = p53.filt)
corr <- cor(as.matrix(t(myeloid.p53[["RNA"]]@counts["SPP1",,drop = F])), as.matrix(t(myeloid.p53[["RNA"]]@counts)))
corr <- t(corr)
corr <- corr[order(-corr[,"SPP1"]),]
#write.csv(corr, file = file.path(projdir, "F4_D_enrichment.csv"))



#### F4G ####
all = get (load (file.path('./dropbox_data','10x_All_integratedRef16U.Rda')))
all3 = get (load (file.path('./dropbox_data','10x_All_integratedRef16U_updated_metadata.Rda')))
all@meta.data = cbind (all@meta.data, all3)

genes <- c("CXCR3", "CXCL11")

# all <- NormalizeData(all)
all$celltype2 = all$celltype_2
all$celltype2[all$celltype2 %in% c('SmoothMuscle','CAF')] = 'Mesenchymal'
all$celltype2[all$celltype2 %in% c('Endothelial','Lymphatic')] = 'Endothelial'

Idents(all) <- "celltype2"
Idents(all) = factor (Idents(all), levels = rev(c('Malignant','Nonmalig Epi','Endothelial','Mesenchymal','Pericyte','Mast','Bcell','Plasma','NK','Tcell','Myeloid')))
# Idents(all) <- "subtype"
pdf(file.path(figures.dir,'F4G.pdf'), useDingbats = F, width = 4, height = 4)
  dp = DotPlot(object = all, features = genes, scale = T) +
  theme(axis.text.x = element_text(angle = 45, hjust=1), panel.border = element_rect(colour = "black", fill=NA, size=0.5), panel.grid.major = element_line(colour = "gainsboro")) +
    scale_colour_gradient2(low = "navy", high = "firebrick") +
    geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.5)
    dp
dev.off()
###


#### F4E ####
meta = get(load(paste0(source.data.dir, "Fig4_myeloid/Fig4E_modulescore.Rda")))

####
genes <- c("hypoxia.hallmark1", "emt.hallmark1", "infg.response.hallmark1")

####
obj = myeloid
obj.subtype.list <- SplitObject(obj, split.by = "subtype")
avglog2fc.list <- list()
pval.list <- list()

for (i in 1:length(obj.subtype.list)) {
  subtype <- obj.subtype.list[[i]]
  
  ### 
  signatures <- t(subtype@meta.data[,genes])
  signatures.obj <- CreateSeuratObject(signatures)
  signatures.obj$orig.identSec <- subtype$orig.identSec
  Idents(signatures.obj) <- "orig.identSec"
  subtype <- AverageExpression(signatures.obj, slot = "counts", assays = "RNA")[[1]]
  subtype <- subtype[,c(p53.filt[p53.filt %in% colnames(subtype)], WT.filt[WT.filt %in% colnames(subtype)])]
  ###
  subtype <- CreateSeuratObject(subtype)
  
  num.p53 <- length(p53.filt[p53.filt %in% colnames(subtype)])
  num.wt <- length(WT.filt[WT.filt %in% colnames(subtype)])
  
  subtype$p53_status <- c(rep("p53_mut", num.p53), rep("WT", num.wt))
  subtype$p53_status <- factor(subtype$p53_status, levels = c("WT", "p53_mut"))
  
  Idents(subtype) <- "p53_status"
  
  # for signatures
  counts <- subtype[["RNA"]]@counts
  avglog2fc_vec <- apply(counts, 1, function(x) {
    avglog2fc <- median(x[1:num.p53], na.rm = T) - median(x[(num.p53 + 1): (num.p53 + num.wt)], na.rm = T) # median diff
    avglog2fc
  })
  pvals <- apply(counts, 1, function(x) {
    pval <- wilcox.test(x[1:num.p53], x[(num.p53 + 1): (num.p53 + num.wt)], na.rm = T)$p.value # wilcox.pval
    pval
  })
  print(length(avglog2fc_vec))
  print(class(avglog2fc_vec))
  print(length(pvals))
  avglog2fc.list[[i]] <- as.data.frame(t(as.data.frame(avglog2fc_vec)))
  pval.list[[i]] <- as.data.frame(t(as.data.frame(pvals)))
}


###
avglog2fc.df <- as.data.frame(rbindlist(avglog2fc.list))
pvals.df <- as.data.frame(rbindlist(pval.list))
dim(avglog2fc.df)
dim(pvals.df)
rownames(avglog2fc.df) <- names(obj.subtype.list)
rownames(pvals.df) <- names(obj.subtype.list)
colnames(avglog2fc.df) <- genes # new
colnames(pvals.df) <- genes  # new

pvals.df <- as.matrix(pvals.df)
pvals.df[is.nan(pvals.df)] <- NA
pvals.df <- as.data.frame(pvals.df)
pvals.df[pvals.df <= 0.05] <- "**"
pvals.df[pvals.df > 0.05 & pvals.df <= 0.1] <- "*"

pvals.df[is.na(pvals.df)] <- ""
pvals.df[pvals.df > 0.1] <- ""

###
range <- max(abs(avglog2fc.df), na.rm = T)
range

###
avglog2fc.df[,names(which(sapply(avglog2fc.df, function(x)all(is.na(x)))))] <- NULL # remove cols with all NAs
pvals.df <- pvals.df[,colnames(avglog2fc.df)] # remove cols with all NAs

colnames(avglog2fc.df) = c("Hypoxia", "EMT", "IFNG.Response")
colnames(pvals.df) = c("Hypoxia", "EMT", "IFNG.Response") 

avglog2fc.df.plot = avglog2fc.df[,c("IFNG.Response", "Hypoxia", "EMT")]
pvals.df.plot = pvals.df[,c("IFNG.Response", "Hypoxia", "EMT")]

############
pdf(paste0(figures.dir, "Fig4_myeloid/Fig4E.pdf"), useDingbats = F, width = 4, height = 4)
pheatmap(avglog2fc.df.plot, display_numbers = pvals.df.plot, breaks = seq(-range, range, length.out = 100),  cellheight=13,cellwidth=13, 
         cluster_rows = T,cluster_cols = F
)
dev.off()




#### F4F ####
df <- get(load("/ahg/regevdata/projects/lungCancerBueno/Results/10x_All_190615/P53_paper/lig_rec_differential_highlvl_df_not_allzeros.Rda"))

df <- df[grepl("Myeloid", df$cell.pair, fixed=TRUE),]

df2 <- df[!grepl("Nonmalig", df$cell.pair, fixed = TRUE),]

#end
int_pairs <- c("VEGFA:FLT1", "VEGFA:EPHB2", "TNFRSF10B:TNFSF10", "TGFB2:TGFbeta receptor1", "SPP1:a9b1 complex",
               "SPP1:a4b1 complex", "SPP1:CD44", "PLXNB1:SEMA4D", "NRP2:SEMA3C", "NRP2:PGF", "NRP1:VEGFB",
               "NRP1:VEGFA", "NRP1:SEMA3A", "NRP1:PGF", "MIF:TNFRSF14", "LGALS9:PTPRK", "LGALS9:HAVCR2",
               "GPR37:PSAP", "FLT1:VEGFB", "FLT1:PGF", "CXCL11:CXCR3", "CXCL10:CXCR3", "CTLA4:CD86",
               "COPA:P2RY6", "CD52:SIGLEC10", "CD40:TNFSF13B", "CCL5:ACKR1", "ANXA1:FPR3", "ADRB2:VEGFB", "ACKR1:CCL17")
df2 <- df2[df2$int_name %in% int_pairs,]

range <- as.numeric(max(abs(df2$p53_wt_prop_diff)))

colors <- c(brewer_pal(palette = "Spectral", direction = -1)(7))
pal <- gradient_n_pal(colors)
custom_color_scale <- scale_fill_gradientn(
  colours = pal(c(0, rescale(seq_along(df2$p53_wt_prop_diff)), 1)), 
  values = c(0,rescale(seq_along(df2$p53_wt_prop_diff)),1), #
  limits = c(-range,range)
) 

############
pdf(file = paste0(figures.dir, "Fig4_myeloid/Fig4F.pdf"), useDingbats = F, height = 15, width = 8.3)
p<-ggplot(data = df2, mapping = aes(x=cell.pair, y=int_name, color=p53_wt_prop_diff, size=p53_wt_neglogpval_proptest))+geom_point(shape = 21, aes(colour = as.factor(signif), fill = p53_wt_prop_diff))+
  scale_colour_manual(values=c("00FFFFFF", "darkgray", "black")) + theme_minimal() +
  custom_color_scale +
  theme(text = element_text(size=12), strip.text = element_text(size=22,face='bold'))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
p
dev.off()




### F4H ####
myeloid = readRDS("/ahg/regevdata/projects/lungCancerBueno/Results/10x_All_190615/TP53_paper/myeloid_pDC_mast_integrated.rds")
DefaultAssay(myeloid) = "RNA"

genes <- c('SPP1', 'CCL17', 'PSAP', 'SEMA4D', 'P2RY6', 'TGFBR2', 'LGALS9', 'VEGFB', 'TGFBR1', 'NRP2', 'NRP1', 'CD86', 'CXCL11', 'SIGLEC10', 'VEGFA',
           'TNFSF13B', 'FPR3', 'TNFRSF14', 'ANXA1', 'FLT1', 'TNFSF10', 'CXCL10', 'MIF', 'CCL5')

obj = myeloid


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


###
avglog2fc.df <- as.data.frame(rbindlist(avglog2fc.list))
pvals.df <- as.data.frame(rbindlist(pval.list))
dim(avglog2fc.df)
dim(pvals.df)
rownames(avglog2fc.df) <- names(obj.subtype.list)
rownames(pvals.df) <- names(obj.subtype.list)
colnames(avglog2fc.df) <- genes # new
colnames(pvals.df) <- genes  # new

pvals.df <- as.matrix(pvals.df)
pvals.df[is.nan(pvals.df)] <- NA
pvals.df <- as.data.frame(pvals.df)
pvals.df[pvals.df <= 0.05] <- "**"
pvals.df[pvals.df > 0.05 & pvals.df <= 0.1] <- "*"

pvals.df[is.na(pvals.df)] <- ""
pvals.df[pvals.df > 0.1] <- ""

###
range <- 2

avglog2fc.df[,names(which(sapply(avglog2fc.df, function(x)all(is.na(x)))))] <- NULL # remove cols with all NAs
pvals.df <- pvals.df[,colnames(avglog2fc.df)] # remove cols with all NAs

names<-c("TAM.FABP4", "DC2", "TAM.CXCL", "Mreg.DC", "Myeloid.cycling", "DC3", "DC1", "CD16.Mono", "Mono-Mac", "TAM.SPP1", "CD14.Mono",
         "TAM.AZU1", "TAM.FOLR2","TAM.APOE", "Mast", "Plasmacytoid.DC")

avglog2fc.df = avglog2fc.df[names,]
pvals.df = pvals.df[names,]

############
pdf(paste0(figures.dir, "Fig4_myeloid/Fig4H.pdf"), useDingbats = F, width = 20, height = 8) 
pheatmap(avglog2fc.df, display_numbers = pvals.df, breaks = seq(-range, range, length.out = 100),  cellheight=13,cellwidth=13, 
         cluster_rows = F,cluster_cols = F,
         color = colorRampPalette(c("navy", "white", "firebrick"))(100))
dev.off()


### F5D ####

# T-malignant ----------------------------------------
df <- get(load(file.path(data.path, "lig_rec_differential_htan_df_not_allzeros.Rda")))
pairs1 <- c("PVR_TIGIT", "TIGIT_NECTIN3", "TIGIT_NECTIN2", "PDCD1_FAM3C", "PDCD1_CD274", "KLRB1_CLEC2D", "LGALS9_HAVCR2", "LGALS9_CD47")
print(all(pairs1 %in% df$interacting_pair))

df1 <- df |>
    filter(interacting_pair %in% pairs1 & cell.pair %in% c("Cancer.Tcell", "Tcell.Cancer")) |>
    mutate(interacting_pair = factor(interacting_pair, levels = rev(pairs1)))

# manually remove the unwanted entries
df1 <- df1 |>
    filter(!(interacting_pair == "LGALS9_CD47" & cell.b != "Cancer"))

# df1 |>   
#     select(interacting_pair, p53_wt_neglogpval_fishertest) |> 
#     write.csv(file.path(data.path, "5d_bar.t_malignant.csv"))

p.bar.t_malig <- df1 |>   
    select(interacting_pair, p53_wt_neglogpval_fishertest) |>
    ggplot(aes(x = p53_wt_neglogpval_fishertest, y = interacting_pair)) +
    geom_bar(stat = "identity", fill = "gray30") +
    geom_vline(xintercept = 1.3, linetype = "dashed") +
    theme_classic() +
    xlab(expression(-log[10]("p-value"))) +  # Correctly format the axis label
    ylab("Ligand-Receptor pair")

# T-myeloid ------------------------------------------
pairs2 <- c("CCL4L2_VSIR", "TNF_VSIR", "PDCD1_PDCD1LG2", "PDCD1_FAM3C", "PDCD1_CD274", "HLA-E_KLRC1", "KLRB1_CLEC2D", "LGALS9_HAVCR2", "CTLA4_CD86", "CTLA4_CD80", "LGALS9_CD47")
print(all(pairs2 %in% df$interacting_pair))

df2 <- df |>
    filter(interacting_pair %in% pairs2 & cell.pair %in% c("Myeloid.Tcell", "Tcell.Myeloid")) |>
    mutate(interacting_pair = factor(interacting_pair, levels = rev(pairs2)))

# manually remove the unwanted entries
df2 <- df2 |>
    filter(!(interacting_pair == "LGALS9_HAVCR2" & cell.a != "Myeloid")) |>
    filter(!(interacting_pair == "LGALS9_CD47" & cell.a != "Myeloid")) |>
    filter(!(interacting_pair == "TNF_VSIR" & cell.b != "Myeloid")) |>
    filter(!(interacting_pair == "CCL4L2_VSIR" & cell.a != "Myeloid"))

# df2 |>   
#     select(interacting_pair, p53_wt_neglogpval_fishertest) |>
#     write.csv(file.path(data.path, "5d_bar.t_myeloid.csv"))

p.bar.t_myeloid <- df2 |>   
    select(interacting_pair, p53_wt_neglogpval_fishertest) |>
    ggplot(aes(x = p53_wt_neglogpval_fishertest, y = interacting_pair)) +
    geom_bar(stat = "identity", fill = "gray30") +
    geom_vline(xintercept = 1.3, linetype = "dashed") +
    theme_classic() +
    xlab(expression(-log[10]("p-value"))) +  # Correctly format the axis label
    ylab("Ligand-Receptor pair")

# generate patchwork for bar plots
p.5d.bar <- p.bar.t_malig + p.bar.t_myeloid + patchwork::plot_layout(ncol = 1, heights = c(1, 11/8))
ggsave(filename = file.path(figures.dir, "5d_bar.pdf"), 
       plot = p.5d.bar, 
       device = "pdf", 
       width = 5, 
       height = 5 * 2)
# ----------------------------------------------------
# dot plots
# ----------------------------------------------------
comb.df <- get(load(file.path(data.path, "viz_comb.df.Rda")))
pairs1 <- c("TIGIT--PVR", "TIGIT--NECTIN3", "TIGIT--NECTIN2", "PDCD1--FAM3C", "PDCD1--CD274", "KLRB1--CLEC2D", "HAVCR2--LGALS9", "CD47--LGALS9")
pairs2 <- c("VSIR--CCL4L2", "TNF--VSIR", "PDCD1--PDCD1LG2", "PDCD1--FAM3C", "PDCD1--CD274", "KLRC1--HLA-E", "KLRB1--CLEC2D", "HAVCR2--LGALS9", "CTLA4--CD86", "CTLA4--CD80", "CD47--LGALS9")
selectedInteractions <- c(pairs1, pairs2)
selectedCelltypes <- c("Lymphoid--cancer", "Lymphoid--Myeloid")

comb.df.filtered <- comb.df |>
    mutate(cell.pair = cell.pair.x) |>
    filter(lr.pair %in% selectedInteractions & cell.pair %in% selectedCelltypes)

# Add TP53 mutation information
p53_samples <- c("1172", "1173", "1175", "1182", "BW16", "BW11", "BW14", "BW06")
wt_samples <- c("14", "1174", "1176", "1183", "BW01", "BW04", "BW05", "BW09", "BW19", "BW23")
comb.df.filtered <- comb.df.filtered |>
    mutate(p53 = case_when(sampleid %in% p53_samples ~ "TP53mut", .default = "WT"), 
           p53 = factor(p53, levels = c("WT", "TP53mut")), 
           cell.pair = factor(cell.pair, levels = c("Lymphoid--cancer", "Lymphoid--Myeloid")), 
           sampleid = factor(sampleid, levels = c(wt_samples, p53_samples)))

# T-malignant ----------------------------------------
p.t_malig <- comb.df.filtered |>
    drop_na(sampleid) |>
    filter(lr.pair %in% pairs1 & cell.pair == selectedCelltypes[1]) |>
    ggplot(aes(x = sampleid, y = lr.pair, size = neg_log10pval, color = log2mean)) +
    geom_point() +
    scale_size_continuous(range = c(1, 10), breaks = c(0.5, 1, 2, 3)) +
    scale_color_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0, limits = c(-2, 3)) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        text = element_text(size = 12),
        strip.text = element_text(size = 14, face = "bold")) +
    labs(title = "T cell-malignant",
       x = "Sample ID",
       y = "Ligand-Receptor pair",
       size = "-log10(p-value)",
       color = "Log2 Mean Expression")

# comb.df.filtered |>
#     drop_na(sampleid) |>
#     filter(lr.pair %in% pairs1 & cell.pair == selectedCelltypes[1]) |>
#     select(sampleid, lr.pair, neg_log10pval, log2mean) |>
#     write.csv(file.path(data.path, "5d_dot.t_malignant.csv"))

# T-myeloid ------------------------------------------
p.t_myeloid <- comb.df.filtered |>
    drop_na(sampleid) |>
    filter(lr.pair %in% pairs2 & cell.pair == selectedCelltypes[2]) |>
    ggplot(aes(x = sampleid, y = lr.pair, size = neg_log10pval, color = log2mean)) +
    geom_point() +
    scale_size_continuous(range = c(1, 10), breaks = c(0.5, 1, 2, 3)) +
    scale_color_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0, limits = c(-2, 3)) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        text = element_text(size = 12),
        strip.text = element_text(size = 14, face = "bold")) +
    labs(title = "T cell-myeloid",
       x = "Sample ID",
       y = "Ligand-Receptor pair",
       size = "-log10(p-value)",
       color = "Log2 Mean Expression")

# comb.df.filtered |>
#     drop_na(sampleid) |>
#     filter(lr.pair %in% pairs2 & cell.pair == selectedCelltypes[2]) |>
#     select(sampleid, lr.pair, neg_log10pval, log2mean) |>
#     write.csv(file.path(data.path, "5d_dot.t_myeloid.csv"))

# generate patchwork for dot plots
p.5d.dot <- patchwork::wrap_plots(p.t_malig, p.t_myeloid, nrow = 2) + patchwork::plot_layout(heights = c(1, 11/8))
p.5d.dot

ggsave(filename = file.path(figures.dir, "5d_dot.pdf"), 
       plot = p.5d.dot, 
       device = "pdf", 
       width = 10, 
       height = 5 * 2)
# ----------------------------------------------------




### S7A ####
meta_df = table (data.frame (sample = myeloid$orig.identSec, subtype = myeloid$subtype))
meta_df = as.data.frame (prop.table (meta_df, 1))
meta_df$sample = factor(meta_df$sample, levels = c(
  '14',
  '1174',
  '1176',
  '1183',
  'BW01',
  'BW04',
  'BW05',
  'BW09',
  'BW19',
  'BW23',
  '1172',
  '1173',
  '1175',
  '1182',
  'BW16',
  'BW11',
  'BW14',
  'BW06',
  '19',
  '1181',
  '1179',
  '1170',
  '17'))

bp = ggplot (meta_df, aes (x = sample, y = Freq, fill = subtype)) + 
  geom_bar (position = 'stack', stat = 'identity', color='black') +
  theme(
      axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
      axis.line =element_line(colour = 'black', size = .1),
        axis.ticks = element_line(colour = "black", size = .1),
      panel.background = element_blank()#,
    #panel.border = element_blank(),
    #panel.grid.major = element_blank(),
    #panel.grid.minor = element_blank()
  )

pdf (file.path (figures.dir,'S7A.pdf'), height=4)
bp
dev.off()

#write.csv (bp$data, file.path (projdir, 'S7A.csv'))


### S8H ####
Idents(myeloid) <- "subtype"
markers <- FindAllMarkers(myeloid, only.pos = T, max.cells.per.ident = 1000)

markers <- myeloid.markers
markers[,"cluster"] <- as.character(markers[,"cluster"]) # in case this is a factor

library(dplyr)

markers <- markers %>%
  arrange(cluster, desc(avg_log2FC)) %>% 
  group_by(cluster) %>%
  slice(1:10)
  

# dir.create(paste0(figures.dir, "Figures_v4/Fig5_myeloid_bmast/"))
#write.csv(markers, file = file.path(projdir, "S7_I_myeloid.markers.top10.csv"))

order <- rev(unique(markers$cluster)) # already sorted

Idents(myeloid) <- factor(Idents(myeloid), levels = order)
genes <- unique(markers$gene)
Idents(myeloid) <- "subtype"
Idents (myeloid) = factor (Idents(myeloid), levels = rev(c(
  'TAM.FABP4',
  'DC2',
  'TAM.CXCL',
  'Mreg.DC',
  'Myeloid.cycling',
  'DC3',
  'DC1',
  'CD16.Mono',
  'Mono-Mac',
  'TAM.SPP1',
  'CD14.Mono',
  'TAM.AZU1',
  'TAM.FOLR2',
  'TAM.APOE',
  'Mast',
  'Plasmacytoid.DC')))

#myeloid <- FindVariableFeatures(myeloid)

mye_markers = ('SPP1
CCL17
PSAP
SEMA4D
P2RY6
TGFBR2
LGALS9
VEGFB
TGFBR1
NRP2
NRP1
CD86
CXCL11
SIGLEC10
VEGFA
TNFSF13B
FPR3
TNFRSF14
ANXA1
FLT1
TNFSF10
CXCL10
MIF
CCL5'
)

mye_markers = unlist(strsplit(mye_markers,'\n'))
DefaultAssay (myeloid) = 'RNA'
myeloid <- ScaleData(myeloid, features = c(VariableFeatures(myeloid),mye_markers)) # add genes of interest to the heatmap

pdf(file.path(figures.dir, 'reproducing_figures','S7I.pdf'), useDingbats = F, width = 10, height = 6)
 dp= DotPlot(object = myeloid, features = mye_markers, scale = T) +
  theme(axis.text.x = element_text(angle = 45, hjust=1), panel.border = element_rect(colour = "black", fill=NA, size=0.5), panel.grid.major = element_line(colour = "gainsboro")) + 
    #scale_color_gradientn(colours = rev(brewer.pal(11,"Spectral"))) +
      scale_colour_gradient2(low = "navy", high = "firebrick") +
    geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.5)
    dp
dev.off()
#write.csv (dp$data, file.path (projdir,'S7_I.csv'))


#### S8G ####
gene = 'SPP1'
gene_cor = cor (t(as.matrix(myeloid@assays$RNA@counts[gene,, drop=F])), as.matrix(t(myeloid@assays$RNA@counts)))
gene_cor = as.data.frame (t(gene_cor))
top_cor_genes = head (gene_cor[order (-gene_cor[,1]),, drop=F],16)
top_cor_genes = top_cor_genes[rownames(top_cor_genes) != gene,,drop=F]
top_cor_genes$gene = factor (rownames(top_cor_genes), levels = rev(rownames(top_cor_genes)))
bp = ggplot (top_cor_genes, aes (x = SPP1, y = gene)) + geom_bar(stat = 'identity') + theme_classic()

pdf (file.path (figures.dir, 'reproducing_figures', 'S7H.pdf'))
bp
dev.off()
