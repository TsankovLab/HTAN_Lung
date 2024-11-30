library (Seurat)
library(pals)
library(RColorBrewer)
library(ggplot2)
library(scales)
library(ggpubr)
library(dplyr)
library(Matrix)

# my outputs
data.path <- "./data/"
figures.dir <- "../../figures/"

######
source(paste0(data.path, "/R_utils/plotutils.R"))
source(paste0(data.path, "/R_utils/seuratutils.R"))
source(paste0(data.path, "/R_utils/seuratutilsV3.R"))
source(paste0(data.path, "/R_utils/color.R"))

myeloid = readRDS (file.path ('./dropbox_data','myeloid_pDC_mast_integrated.rds'))#myeloid = NormalizeData (myeloid)
###

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


### FIGURE S7A ####
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

### FIGURE S7I ####
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


#### S7H ####
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
