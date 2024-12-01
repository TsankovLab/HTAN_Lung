######################@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ 
######################@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ FIGURE 2 - EXTENDED DATA 3 - MALIGNANT PANELS
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
data.path <- "./data/"
figures.dir <- "../../figures/"

if (!file.exists(paste0(figures.dir,'Fig2_malig/'))){dir.create(paste0(figures.dir,'Fig2_malig/'), recursive=TRUE)}

######
source(paste0(data.path, "/R_utils/plotutils.R"))
source(paste0(data.path, "/R_utils/seuratutils.R"))
source(paste0(data.path, "/R_utils/seuratutilsV3.R"))
source(paste0(data.path, "/R_utils/color.R"))



#####################@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ 
#####################@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
##################### Figure 2A
malig <- get(load("./dropbox_data/maligv9.allsamplesv4.Rda"))

Idents(malig)  <- "subtype"
malig.names<-c("AT2.like", "AT1.2.like", "CC.G2M", "CC.S", "Ciliated", 
                 "Glycolysis.Hypox", "MHCII", "Hypoxia", "IFNG", "StressResponse", 
                 "Respiration.MT", "Metallothionein", "OxPhos", "PEMT", "Ribosome", 
                 "Secretory", "Senescence", "TNFA.NFKB")
Idents(malig) <- factor(Idents(malig), levels = malig.names)

###########
pdf(paste0(figures.dir, "Fig2_malig/Fig2A.pdf"), useDingbats = F, width = 10)
p = DimPlot(malig)
p
dev.off()



#####################@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ 
#####################@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
##################### Figure 2B
malig <- get(load("./dropbox_data/maligv9.allsamplesv4.Rda"))

genes <- c("SFTPA2", "SFTPA1", "SFTPC", "AKR1C1", "PTTG1", "CDKN3", "HIST1H4C", "PCNA", "CAPS", "CETN2", "PKM", "LDHA", "HLA-DRA", "HLA-DPA1", 'NDRG1',
           "SLC2A1", "ISG15", "IFI6", "HSPA1A", "JUN", "MT-ND5", "MT-ND6", "MT2A", "MT1X", "ATP5I", "SEPP1", "LAMC2", "VIM", "RPL32", "RPL29","LCN2", "SLPI", "S100A4", "KRT19", "NFKBIA", "CXCL2")

DefaultAssay(malig) <- "RNA"
malig <- NormalizeData(malig)
malig <- FindVariableFeatures(malig)
malig <- ScaleData(malig, features = unique(VariableFeatures(malig), genes))

#
Idents(malig)  <- "subtype"

malig.names<-c("AT2.like", "AT1.2.like", "CC.G2M", "CC.S", "Ciliated", 
                 "Glycolysis.Hypox", "MHCII", "Hypoxia", "IFNG", "StressResponse", 
                 "Respiration.MT", "Metallothionein", "OxPhos", "PEMT", "Ribosome", 
                 "Secretory", "Senescence", "TNFA.NFKB")

Idents(malig) <- factor(Idents(malig), levels = rev(malig.names)) 


###########
pdf(paste0(figures.dir, "Fig2_malig/Fig2B.pdf"), useDingbats = F, width = 14, height = 5)
p = DotPlot(object = malig, features = genes, scale = T) +
  theme(axis.text.x = element_text(angle = 45, hjust=1), panel.border = element_rect(colour = "black", fill=NA, size=0.5), panel.grid.major = element_line(colour = "gainsboro")) + 
  scale_colour_gradient2(low = "navy", high = "firebrick") +
  geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.5)
p
dev.off()



#####################@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ 
#####################@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
##################### Figure 2C
res = get(load(paste0(data.path, "malig.corr.cells.Rda")))
rownames(res)=colnames(res)=c("AT2.like", "AT1.2.like", "CC.G2M", "CC.S", "Ciliated", 
                              "Glycolysis.Hypox", "MHCII", "Hypoxia", "IFNG", "StressResponse", 
                              "Respiration.MT", "Metallothionein", "OxPhos", "PEMT", "Ribosome", 
                              "Secretory", "Senescence", "TNFA.NFKB")
ann = c("Metallothionein","Hypoxia","Senescence","Glycolysis.Hypox","OxPhos","CC.G2M", "CC.S", "Ciliated","PEMT","Secretory","TNFA.NFKB","IFNG",
        "Ribosome","StressResponse","AT2.like","MHCII","Respiration.MT","AT1.2.like")
res.ann = res[ann,ann]

###########
pdf(paste0(figures.dir,"Fig2_malig/Fig2C.pdf"), width = 10, height = 9, useDingbats=F)
p<-Heatmap(res.ann, cluster_rows = F, cluster_columns = F)
p
dev.off()




#####################@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ 
#####################@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
##################### Figure 2E (left)
malig <- get(load("./dropbox_data/maligv9.allsamplesv4.Rda"))

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

### 2E (left) plots only
mean.plots.2 = list()
mean.plots.2[[1]] = mean.plots[[19]] ## tp53 targets
mean.plots.2[[2]] = mean.plots[[1]]   ## AT2.like
mean.plots.2[[3]] = mean.plots[[3]]   ## CC.G2M
mean.plots.2[[4]] = mean.plots[[6]]   ## Glyco.hypox
mean.plots.2[[5]] = mean.plots[[14]]   ## PEMT

pdf(paste0(figures.dir, "Fig2_malig/Fig2E_left.pdf"), width = 3, height = 15, useDingbats=FALSE)
cowplot::plot_grid(plotlist = mean.plots.2, ncol=1)
dev.off()



#####################@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ 
#####################@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
##################### Figure 2G
lung.sigs.avg.rbind <- get(load(paste0(data.path, "lung.sigs.avg.rbind.Rda")))

snps <- log10(c(151, 412, 43, 694, 73, 28, 125, 66, 657, 295, 468, 57, 4227, 947, 104, 66, 84, 380, 443, 182, 544, 347, 64))
p53_status <- c("WT", "WT", "WT", "WT", "WT", "WT", "WT", "WT", "WT", "WT",
                "Mutant", "Mutant", "Mutant", "Mutant", "Mutant", "Mutant", "Mutant", "Mutant", 
                "Mutant", "Mutant", "WT", "WT", "Mutant")
histology <- c(rep("LUAD", 18), rep("LUSC", 2), rep("NE", 1), rep("Mucinous/Colloid", 2))

lung.sigs.avg.rbind <- lung.sigs.avg.rbind[,c(WT.filt, p53.filt, other)]
ann = c("P14", "MGH1174", "MGH1176", "MGH1183", "BWH01", "BWH04", "BWH05", "BWH09", 
        "BWH19", "BWH23", "MGH1172", "MGH1173", "MGH1175", "MGH1182", "BWH16", "BWH11", 
        "BWH14", "BWH06", "P19", "MGH1181", "MGH1179", "MGH1170", "P17")
colnames(lung.sigs.avg.rbind) = ann

###
dark2 <- brewer.pal(n = 4, name = "Dark2")

ha = HeatmapAnnotation(
  histology = histology,
  p53_status = p53_status,
  col = list(p53_status = c("Mutant" = "#0091CA", "WT" = "#D8423D"),
             histology = c("LUAD" = dark2[1], "LUSC" = dark2[2],
                           "NE" = dark2[3], "Mucinous/Colloid" = dark2[4]))
)

col_fun <- colorRamp2(c(min(snps), max(snps)), c("white", "#4B0082"))
hb = HeatmapAnnotation(
  snps = snps,
  col = list(snps = col_fun)
)

############ 
pdf(paste0(figures.dir, "Fig2_malig/Fig2G.pdf"), useDingbats = F, width = 8, height = 3)
p<-Heatmap(lung.sigs.avg.rbind, cluster_rows = F, cluster_columns = F, top_annotation = ha, bottom_annotation = hb)
p
dev.off()



#####################@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ 
#####################@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
##################### Figure 2I (leftmost)
ccat.df <- get(load(paste0(data.path, "ccat.df.sc.nonmalig.wt.mut.subset.Rda"))) 

ccat.df$p53_status <- as.character(ccat.df$p53_status)
ccat.df$p53_status[ccat.df$p53_status %in% "nonmalig"] <- "Normal"
ccat.df$p53_status <- factor(ccat.df$p53_status, levels = c("Normal", "WT", "p53_mut"))
my_comparisons <- list(c("Normal", "WT"), c("WT", "p53_mut"), c("Normal", "p53_mut"))

############ 
pdf(paste0(figures.dir, "Fig2_malig/Fig2I_left.pdf"), width = 4, height = 4)
ggviolin(ccat.df, x = "p53_status", y = "ccat.v", color = "black", fill = "p53_status",
         title = "Single cell entropy", palette = c("grey", "#0091CA", "#D8423D"),
         add = "boxplot", add.params = list(fill = "white")) +
  NoLegend() +
  stat_compare_means(comparisons = my_comparisons)
dev.off()




#####################@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ 
#####################@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
##################### Figure 2J
ccat.df <- get(load(paste0(data.path, "ccat.df_smartseq2.Rda")))

ccat.df = ccat.df[which(ccat.df$timesimple %in% c("01_T_early_ND", "05_K_30w_ND", "08_KP_30w_ND")),]

ccat.df$K_vs_KP <- as.character(ccat.df$K_vs_KP)
ccat.df$K_vs_KP[ccat.df$K_vs_KP %in% "T"] <- "Normal"
ccat.df$K_vs_KP <- factor(ccat.df$K_vs_KP, levels = c("Normal", "K", "KP"))
my_comparisons <- list(c("Normal", "K"), c("K", "KP"), c("Normal", "KP"))

############
pdf(paste0(figures.dir, "Fig2_malig/Fig2J.pdf"), width = 4, height = 4)
ggviolin(ccat.df, x = "K_vs_KP", y = "ccat.v", color = "black", fill = "K_vs_KP",
         title = "K vs KP mouse entropy", palette = c("grey", "#0091CA", "#D8423D"),
         add = "boxplot", add.params = list(fill = "white")) +
  NoLegend() +
  stat_compare_means(comparisons = my_comparisons)
dev.off()



#####################@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ 
#####################@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
##################### Figure 2K
ccat.df <- get(load(paste0(data.path, "ccat.df.tp53.celllines.Rda")))

#
ccat.df$impact <- factor(ccat.df$impact, levels = c("WT TP53", "WT-like TP53", "Impactful I TP53", "Impactful II TP53",
                                                    "R280K TP53 I", "R158L TP53 II", "H179R TP53 II", "Y205C TP53 II", "S241F TP53 II"))

my_comparisons <- list(c("WT TP53", "WT-like TP53"), c("WT TP53", "Impactful I TP53"),
                       c("WT TP53", "Impactful II TP53"), c("WT TP53", "R280K TP53 I"),
                       c("WT TP53", "R158L TP53 II"), c("WT TP53", "H179R TP53 II"),
                       c("WT TP53", "Y205C TP53 II"), c("WT TP53", "S241F TP53 II"))

my_color_palette <- c("#0091CA", "#8DCEE7", "#F2C0BE", "#EBA09E", "#E1716D", "#C63C37", "#A2312D", "#7E2623", "#6C211E") 

############
pdf(paste0(figures.dir, "Fig2_malig/Fig2K.pdf"), width = 4.5, height = 4.5)  
ggviolin(ccat.df, x = "impact", y = "ccat.v", fill = "impact",
         palette = my_color_palette,
         add = "boxplot", add.params = list(fill = "white"))+
  rotate_x_text(angle = 45) + NoLegend() +
  stat_compare_means(comparisons = my_comparisons)
dev.off()



#####################@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ 
#####################@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
##################### Extended Figure 3A
malig <- get(load("./dropbox_data/maligv9.allsamplesv4.Rda"))

###########
pdf(paste0(figures.dir, "Fig2_malig/Ext_Fig2A.pdf"), useDingbats = F, width = 10)
p = DimPlot(malig, group.by = "orig.identSec")
p
dev.off()



#####################@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ 
#####################@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
##################### Extended Figure 3B
library(clusterProfiler)
library(msigdbr)
library(org.Hs.eg.db)

#####
degs <- read.csv(paste0(data.path, "Top100.Cancer.subtype.markers.wilcox.rankedlogfc.csv"))
deg.epi = degs
deg.epi$X <- NULL
se = which(is.na(deg.epi$cluster))
deg.epi = deg.epi[-se,]

### set up the gene set
H_t2g <- msigdbr(species = "Homo sapiens", category = "H") %>%
  dplyr::select(gs_name, gene_symbol)
head(H_t2g)

########### 
deg.epi$cluster <- gsub("\\.", "_", deg.epi$cluster)

xx <- compareCluster(gene ~ cluster, data = deg.epi, fun = enricher,
                     TERM2GENE=H_t2g)

########### 
p <- dotplot(xx, x="cluster") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
                                      axis.text.y = element_text(size = 6, vjust = 0.5, hjust=1))
pdf(paste0(figures.dir, "Fig2_malig/Ext_Fig2B.pdf"), useDingbats = F, width = 8, height = 10)
print(p)
dev.off()



#####################@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ 
#####################@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
##################### Extended Figure 3C
res = get(load(paste0(data.path, "corr.pts.Rda")))
rownames(res)=colnames(res)=c("AT2.like", "AT1.2.like", "CC.G2M", "CC.S", "Ciliated", 
                              "Glycolysis.Hypox", "MHCII", "Hypoxia", "IFNG", "StressResponse", 
                              "Respiration.MT", "Metallothionein", "OxPhos", "PEMT", "Ribosome", 
                              "Secretory", "Senescence", "TNFA.NFKB")
ann = c("Glycolysis.Hypox","Metallothionein","Hypoxia","CC.G2M", "CC.S","PEMT","Senescence","OxPhos", "Ciliated","IFNG","Ribosome","Secretory","TNFA.NFKB",
        "MHCII","AT2.like","Respiration.MT","StressResponse","AT1.2.like")
res.ann = res[ann,ann]


pdf(paste0(figures.dir,"Fig2_malig/Ext_Fig2C.pdf"), width = 10, height = 9, useDingbats=F)
p<-Heatmap(res.ann, cluster_rows = F, cluster_columns = F)
p
dev.off()



#####################@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ 
#####################@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
##################### Extended Figure 3D 
meta = get(load(paste0(data.path, "meta_30_20_top_nmf_genes.Rda")))

malig.score = meta[,c(52:69)]
colnames(malig.score) = substr(colnames(malig.score), 1, (nchar(colnames(malig.score))-1))
cnmf.score = meta[,c(71:(ncol(meta)-1))]

res = cor(malig.score, cnmf.score, method = "p")

res2 = res[c("AT.SFTPC", "MT.respiration", "HLA.D.AT2", "AT.SFTPA2.1", "JUN.FOS", "TNFA.NFKB", 
             "Ribosome", "IFNG", "Secretory", "Ciliated", "PEMT", "Metallothionein", 
             "Hypoxia", "Glycolysis.Hypox", "Senescence", "OxPhos", "CC.S",
             "CC.G2M"),
           c("CN29", "CN1", "CN3", "CN23", "CN22", "CN12", "CN26", "CN19", "CN13", "CN30", "CN21", "CN25", "CN8", "CN5", "CN9", "CN10", "CN20", "CN7", "CN27", "CN2",
             "CN14", "CN28", "CN18", "CN17", "CN16", "CN11", "CN15", "CN6", "CN4", "CN24")]
rownames(res2) = c("AT1.2.like", "Respiration.MT", "MHCII", "AT2.like", "StressResponse", "TNFA.NFKB", 
                   "Ribosome", "IFNG", "Secretory", "Ciliated", "PEMT", "Metallothionein", 
                   "Hypoxia", "Glycolysis.Hypox", "Senescence", "OxPhos", "CC.S",
                   "CC.G2M")

############
pdf(paste0(figures.dir, "Fig2_malig/Ext_Fig3D_newfigname.pdf"), width = 10, height = 8) 
Heatmap(res2, cluster_rows = F, cluster_columns = F, name = "Pearson.cor")
dev.off()



#####################@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ 
#####################@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
##################### Extended Figure 3H
info = read_excel(paste0(data.path, "TableS1.xlsx"))
info = as.data.frame(info)
print(info[,c("Patient","SNPs","CNAs")])

p53.filt.v2 <- c("MGH1172", "MGH1173", "MGH1175", "MGH1182", "BWH16", "BWH11", "BWH14", "BWH06")
WT.filt.v2 <- c("P14", "MGH1174", "MGH1176", "MGH1183", "BWH01", "BWH04", "BWH05", "BWH09", "BWH19", "BWH23") 

comp.df = info[,c("Patient","SNPs","CNAs")]
comp.df[,"SNPs"] = log10(comp.df[,"SNPs"])
colnames(comp.df)[2] = "log10(SNPs)"
comp.df$p53_status = NA
comp.df$p53_status[which(comp.df$Patient %in% p53.filt.v2)] = "TP53_mut"
comp.df$p53_status[which(comp.df$Patient %in% WT.filt.v2)] = "TP53_wt"
se = which(is.na(comp.df$p53_status))
comp.df = comp.df[-se,]
comp.df$p53_status = factor(comp.df$p53_status, levels = c("TP53_wt", "TP53_mut"))

###########
subtypes = c("log10(SNPs)", "CNAs")
mean.plots = list()

for (subtype in subtypes){
  mean.plot <- ggboxplot(comp.df, x = "p53_status", y = subtype, color = "p53_status",
                         palette = c("#0091CA", "#D8423D"), add = "jitter",
                         title = subtype) +
    stat_compare_means(label.x = 0.9, label.y = .9 * max(as.numeric(comp.df[,subtype]))) +
    theme(plot.title = element_text(size = 12, face = "bold")) + NoLegend()
  mean.plots[[subtype]] <- mean.plot
}


pdf(paste0(figures.dir, "Fig2_malig/Ext_Fig3H.pdf"), width = 5, height = 3, useDingbats=FALSE)
cowplot::plot_grid(plotlist = mean.plots, ncol=2)
dev.off()



#####################@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ 
#####################@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
##################### Extended Figure 3I
maligv9.scores.df <- get(load(paste0(data.path, "Ext_Fig3I_modulescore.Rda")))

maligv9.scores.list <- split(maligv9.scores.df, f = maligv9.scores.df$sampleID)
maligv9.scores.list.avg <- lapply(maligv9.scores.list, function(x) {
  y <- colMeans(x[,2:ncol(x)])
  y
})
names(maligv9.scores.list.avg) <- names(maligv9.scores.list)

###
p53.filt.val = c("3", "4", "6", "RU661_TUMOR", "RU1057", "RU1061", "p023t")
WT.filt.val = c("RU1128","RU1137","RU1038","SSN05","SSN24","SSN27","p034t","p031t","p032t","p024t","p019t","p018t","RU684_TUMOR","RU653_TUMOR","RU676_TUMOR")

maligv9.scores.avg.rbind <- list.rbind(maligv9.scores.list.avg)
maligv9.scores.avg.rbind <- as.data.frame(maligv9.scores.avg.rbind)
maligv9.scores.avg.rbind <- maligv9.scores.avg.rbind[c(p53.filt.val, WT.filt.val),]

maligv9.scores.avg.rbind$p53_status <- c(rep("p53_mut", length(p53.filt.val)), rep("WT", length(WT.filt.val)))
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

##
mean.plots <- list()

for (subtype in subtypes){
  mean.plot <- ggboxplot(comp.df, x = "p53_status", y = subtype, color = "p53_status",
                         palette = c("#0091CA", "#D8423D"), add = "jitter",
                         title = subtype) +
    stat_compare_means(label.x = 0.9, label.y = .9 * max(as.numeric(comp.df[,subtype]))) +
    theme(plot.title = element_text(size = 12, face = "bold")) + NoLegend()
  mean.plots[[subtype]] <- mean.plot
}

############
mean.plots.2 = list()
mean.plots.2[[1]] = mean.plots[[19]] ## tp53 targets
mean.plots.2[[2]] = mean.plots[[1]]   ## AT2.like
mean.plots.2[[3]] = mean.plots[[3]]   ## CC.G2M
mean.plots.2[[4]] = mean.plots[[6]]   ## Glyco.hypox
mean.plots.2[[5]] = mean.plots[[14]]   ## PEMT

# 
pdf(paste0(figures.dir, "Fig2_malig/Ext_Fig3I.pdf"), width = 12, height = 4, useDingbats=FALSE)
cowplot::plot_grid(plotlist = mean.plots.2, ncol=5)
dev.off()



#####################@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ 
#####################@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
##################### Extended Figure 4E
maligv9.scores.df <- get(load(paste0(data.path, "addmodulescore.maligprog.A549celllines.Rda")))  

############ tp53
my_color_palette <- c("#0091CA", "#8DCEE7", "#F2C0BE", "#EBA09E", "#E1716D", "#C63C37", "#A2312D", "#7E2623", "#6C211E") # p53 impact
my_comparisons <- list( c("WT TP53", "WT-like TP53"), c("WT TP53", "Impactful I TP53"), c("WT TP53", "Impactful II TP53"),
                        c("WT TP53", "R280K TP53 I"), c("WT TP53", "R158L TP53 II"), c("WT TP53", "H179R TP53 II"), c("WT TP53", "Y205C TP53 II"),
                        c("WT TP53", "S241F TP53 II"))

maligv9.scores.df.subset = maligv9.scores.df[grep("_tp53", rownames(maligv9.scores.df)),]
maligv9.scores.df.subset = maligv9.scores.df.subset[which(maligv9.scores.df.subset$impact != "na"),]


comp.df = maligv9.scores.df.subset
colnames(comp.df) = gsub("1$", "", colnames(comp.df))

subtypes <- colnames(comp.df)[-20]
for (i in 1:(ncol(comp.df)-1))
{
  comp.df[,i] = as.numeric(comp.df[,i])
}

comp.df$impact = factor(comp.df$impact, levels = c("WT TP53", "WT-like TP53", "Impactful I TP53", "Impactful II TP53", "R280K TP53 I", "R158L TP53 II",
                                                   "H179R TP53 II", "Y205C TP53 II", "S241F TP53 II"))

## update subtype names
c("AT.SFTPA2.1", "AT.SFTPC", "CC.G2M", "CC.S", "Ciliated", "Glycolysis.Hypox", 
  "HLA.D.AT2", "Hypoxia", "IFNG", "JUN.FOS", "MT.respiration", 
  "Metallothionein", "OxPhos", "PEMT", "Ribosome", "Secretory", 
  "Senescence", "TNFA.NFKB", "p53_targets")
subtypes[1] = "AT2.like"
subtypes[2] = "AT1.2.like"
subtypes[7] = "MHCII"
subtypes[10] = "StressResponse"
subtypes[11] = "Respiration.MT"

colnames(comp.df)[-20] = subtypes

###
mean.plots <- list()

for (subtype in subtypes){
  mean.plot <- ggviolin(comp.df, x = "impact", y = subtype, color = "black", fill = "impact",
                        title = subtype, palette = my_color_palette,
                        add = "boxplot", add.params = list(fill = "white"))+
    rotate_x_text(angle = 45) + NoLegend() +
    stat_compare_means(comparisons = my_comparisons)
  mean.plots[[subtype]] <- mean.plot
}

###
mean.plots.2 = list()
mean.plots.2[[1]] = mean.plots[[19]] ## tp53 targets
mean.plots.2[[2]] = mean.plots[[1]]   ## AT.SFTPA2.1
mean.plots.2[[3]] = mean.plots[[3]]   ## CC.G2M
mean.plots.2[[4]] = mean.plots[[6]]   ## Glyco.hypox
mean.plots.2[[5]] = mean.plots[[14]]   ## PEMT

# 
pdf(paste0(figures.dir, "Fig2_malig/Ext_Fig4E.pdf"), width = 22, height = 5, useDingbats=FALSE)
cowplot::plot_grid(plotlist = mean.plots.2, ncol=5)
dev.off()



#####################@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ 
#####################@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
##################### Extended Figure 4F
maligv9.scores.df <- get(load(paste0(data.path, "addmodulescore.maligprog.A549celllines.Rda")))  

############ kras
my_color_palette <- c("#0091CA", "#8DCEE7", "#9E9AC8", "#756BB1", "#54278F", "#3F007D") # kras impact

my_comparisons <- list( c("WT KRAS", "WT-like KRAS"), c("WT KRAS", "Impactful I KRAS"), c("WT KRAS", "Impactful II KRAS"),
                        c("WT KRAS", "Impactful III KRAS"), c("WT KRAS", "Impactful IV KRAS"))

maligv9.scores.df.subset = maligv9.scores.df[grep("_kras", rownames(maligv9.scores.df)),]
maligv9.scores.df.subset = maligv9.scores.df.subset[which(maligv9.scores.df.subset$impact %in% c("WT KRAS", "WT-like KRAS", "Impactful I KRAS", "Impactful II KRAS", "Impactful III KRAS", "Impactful IV KRAS")),]

comp.df = maligv9.scores.df.subset
colnames(comp.df) = gsub("1$", "", colnames(comp.df))

subtypes <- colnames(comp.df)[-20]
for (i in 1:(ncol(comp.df)-1))
{
  comp.df[,i] = as.numeric(comp.df[,i])
}

comp.df$impact = factor(comp.df$impact, levels = c("WT KRAS", "WT-like KRAS", "Impactful I KRAS", "Impactful II KRAS", "Impactful III KRAS", "Impactful IV KRAS"))

## update subtype names
c("AT.SFTPA2.1", "AT.SFTPC", "CC.G2M", "CC.S", "Ciliated", "Glycolysis.Hypox", 
  "HLA.D.AT2", "Hypoxia", "IFNG", "JUN.FOS", "MT.respiration", 
  "Metallothionein", "OxPhos", "PEMT", "Ribosome", "Secretory", 
  "Senescence", "TNFA.NFKB", "p53_targets")
subtypes[1] = "AT2.like"
subtypes[2] = "AT1.2.like"
subtypes[7] = "MHCII"
subtypes[10] = "StressResponse"
subtypes[11] = "Respiration.MT"

colnames(comp.df)[-20] = subtypes

###
mean.plots <- list()

for (subtype in subtypes){
  mean.plot <- ggviolin(comp.df, x = "impact", y = subtype, color = "black", fill = "impact",
                        title = subtype, palette = my_color_palette,
                        add = "boxplot", add.params = list(fill = "white"))+
    rotate_x_text(angle = 45) + NoLegend()
  stat_compare_means(comparisons = my_comparisons)
  mean.plots[[subtype]] <- mean.plot
}

###
mean.plots.2 = list()
mean.plots.2[[1]] = mean.plots[[19]] ## tp53 targets
mean.plots.2[[2]] = mean.plots[[1]]   ## AT2.like
mean.plots.2[[3]] = mean.plots[[3]]   ## CC.G2M

###
pdf(paste0(figures.dir, "Fig2_malig/Ext_Fig4F.pdf"), width = 12, height = 5, useDingbats=FALSE)



#####################@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ 
#####################@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
##################### Extended Figure 5A (left-programs)
maligv9.scores.df <- get(load(paste0(data.path, "addmodulescore.maligprog.KPmouse.smartseq.Rda")))

comp.df = maligv9.scores.df
colnames(comp.df) = gsub("1$", "", colnames(comp.df))

subtypes <- colnames(comp.df)[-20]
for (i in 1:(ncol(comp.df)-1))
{
  comp.df[,i] = as.numeric(comp.df[,i])
}

colnames(comp.df)[20] = "status"
comp.df$status = as.character(comp.df$status)
comp.df$status[which(comp.df$status == "T")] = "Normal"
comp.df$status = factor(comp.df$status, levels = c("Normal", "K", "KP"))

#######
## update subtype names
c("AT.SFTPA2.1", "AT.SFTPC", "CC.G2M", "CC.S", "Ciliated", "Glycolysis.Hypox", 
  "HLA.D.AT2", "Hypoxia", "IFNG", "JUN.FOS", "MT.respiration", 
  "Metallothionein", "OxPhos", "PEMT", "Ribosome", "Secretory", 
  "Senescence", "TNFA.NFKB", "p53_targets")
subtypes[1] = "AT2.like"
subtypes[2] = "AT1.2.like"
subtypes[7] = "MHCII"
subtypes[10] = "StressResponse"
subtypes[11] = "Respiration.MT"

colnames(comp.df)[-20] = subtypes

###
my_comparisons = list(c("Normal", "K"), c("Normal", "KP"), c("K", "KP"))
mean.plots <- list()

for (subtype in subtypes){
  mean.plot <- ggviolin(comp.df, x = "status", y = subtype, color = "black", fill = "status",
                        title = subtype, palette = c("grey", "#0091CA", "#D8423D"),
                        add = "boxplot", add.params = list(fill = "white"))+
    rotate_x_text(angle = 45) + NoLegend() +
    stat_compare_means(comparisons = my_comparisons)
  mean.plots[[subtype]] <- mean.plot
}

############ 
mean.plots.2 = list()
mean.plots.2[[1]] = mean.plots[[19]] ## tp53 targets
mean.plots.2[[2]] = mean.plots[[1]]   ## AT.SFTPA2.1
mean.plots.2[[3]] = mean.plots[[4]]   ## CC.S
mean.plots.2[[4]] = mean.plots[[6]]   ## Glyco.hypox
mean.plots.2[[5]] = mean.plots[[8]]   ## Hypox
mean.plots.2[[6]] = mean.plots[[14]]   ## PEMT

# 
pdf(paste0(figures.dir, "Fig2_malig/Ext_Fig5A_left.pdf"), width = 20, height = 5, useDingbats=FALSE) 
cowplot::plot_grid(plotlist = mean.plots.2, ncol=6)
dev.off()



#####################@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ 
#####################@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
##################### Extended Figure 5A (Right-Hif1a)
comp.df = get(load(paste0(data.path, "exp_Hif1a_kp_mouse_30wks.Rda")))  

#############
colnames(comp.df)[2] = "status"
comp.df$status = as.character(comp.df$status)
comp.df$status[which(comp.df$status == "T")] = "Normal"
comp.df$status = factor(comp.df$status, levels = c("Normal", "K", "KP"))

my_comparisons = list(c("Normal", "K"), c("Normal", "KP"), c("K", "KP"))

mean.plot <- ggviolin(comp.df, x = "status", y = "HIF1A", color = "black", fill = "status",
                      palette = c("grey", "#0091CA", "#D8423D"),
                      add = "boxplot", add.params = list(fill = "white"))+
  rotate_x_text(angle = 45) + NoLegend() +
  stat_compare_means(comparisons = my_comparisons)



#####################@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ 
#####################@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
##################### Extended Figure 5D
res = get(load(paste0(data.path, "ccat.df.malig.subset.Rda")))

#
res = data@meta.data[,c("ccat.v", "subtype")]

subtype = unique(res$subtype)
res.ann = matrix(nrow=length(subtype), ncol = 2)
rownames(res.ann) = subtype
colnames(res.ann) = c("ccat.v", "subtype")

for (i in 1:length(subtype))
{
  tmp = res[which(res$subtype == subtype[i]),]
  res.ann[i,1] = mean(tmp[,1])
  res.ann[i,2] = subtype[i]  
}

dat = as.numeric(res.ann[,1])
dat.plot = matrix(dat, ncol = 1)
rownames(dat.plot) = rownames(res.ann)

#############
pdf(paste0(figures.dir, "Fig2_malig/Ext_Fig5D.pdf"), width = 3, height = 4)
p<-Heatmap(dat.plot, cluster_rows = T, cluster_columns = F)
p
dev.off()



#####################@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ 
#####################@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
##################### Extended Figure 5E
conf_mtx = get(load(paste0(data.path, "conf_mtx.Rda")))
conf_mtx <- as.data.frame(conf_mtx, drop = F)
conf_mtx = conf_mtx[which(!conf_mtx$HPCS==0),,drop=F]

conf_mtx$types <- rownames(conf_mtx)
conf_mtx$x <- "HPCS"

############ 
pdf(file = paste0(figures.dir, "Fig2_malig/Ext_Fig5E.pdf"), useDingbats = F, width = 2.3, height = 4)
ggplot(conf_mtx, aes(fill=types, y=HPCS, x=x)) +
  geom_bar(position="stack", stat="identity") +
  scale_fill_manual(values=material.heat(7)) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                                     panel.background = element_blank(), axis.line = element_line(colour = "black"))
dev.off()



#####################@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ 
#####################@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
##################### Extended Figure 5F
ccat.df = get(load(paste0(data.path,'ccat.df.kras.celllines.Rda')))

#
ccat.df$impact <- factor(ccat.df$impact, levels = c("WT KRAS", "WT-like KRAS",
                                                    "Impactful I KRAS", "Impactful II KRAS",
                                                    "Impactful III KRAS", "Impactful IV KRAS"))

my_comparisons <- list(c("WT KRAS", "WT-like KRAS"), c("WT KRAS", "Impactful I KRAS"),
                       c("WT KRAS", "Impactful II KRAS"), c("WT KRAS", "Impactful III KRAS"),
                       c("WT KRAS", "Impactful IV KRAS"))

my_color_palette <- c("#0091CA", "#8DCEE7", "#9E9AC8", "#756BB1", "#54278F", "#3F007D") # kras impact

############
pdf(paste0(figures.dir, "Fig2_malig/Ext_Fig5F.pdf"), width = 5, height = 5)  
ggviolin(ccat.df, x = "impact", y = "ccat.v", fill = "impact",
         palette = my_color_palette,
         add = "boxplot", add.params = list(fill = "white"))+
  rotate_x_text(angle = 45) + NoLegend() +
  stat_compare_means(comparisons = my_comparisons)
dev.off()