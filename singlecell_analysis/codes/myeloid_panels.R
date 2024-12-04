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


######
source(paste0(data.path, "/R_utils/plotutils.R"))
source(paste0(data.path, "/R_utils/seuratutils.R"))
source(paste0(data.path, "/R_utils/seuratutilsV3.R"))
source(paste0(data.path, "/R_utils/color.R"))

#Genesets
cc.genes <- readLines(con = paste0(user.path, "/genelists/regev_lab_cell_cycle_genes.txt"))
s.genes <- cc.genes[1:43]
g2m.genes <- cc.genes[44:97]

myeloid = readRDS (file.path ('../dropbox_data','myeloid_pDC_mast_integrated.rds'))#myeloid = NormalizeData (myeloid)
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

#### F4C_S8C ####
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

pdf(file.path(figures.dir,'F4C_S8C.pdf'), width = 18, height = 7.5, useDingbats=FALSE)
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
all = get (load (file.path('../dropbox_data','all.merge.new.cci.nodoublets.v4.Rda')))

genes <- c("CXCR3", "CXCL11")

# all <- NormalizeData(all)

Idents(all) <- "celltype"
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
myeloid = readRDS (file.path ('../dropbox_data','myeloid_pDC_mast_integrated.rds'))
### add module score for emt and hypoxia 
hypoxia.hallmark <- c("ACKR3","ADM","ADORA2B","AK4","AKAP12","ALDOA","ALDOB","ALDOC","AMPD3","ANGPTL4","ANKZF1","ANXA2","ATF3","ATP7A","B3GALT6","B4GALNT2","BCAN","BCL2","BGN","BHLHE40","BNIP3L","BRS3","BTG1","CA12","CASP6","CAV1","CAVIN1","CAVIN3","CCN1","CCN2","CCN5","CCNG2","CDKN1A","CDKN1B","CDKN1C","CHST2","CHST3","CITED2","COL5A1","CP","CSRP2","CXCR4","DCN","DDIT3","DDIT4","DPYSL4","DTNA","DUSP1","EDN2","EFNA1","EFNA3","EGFR","ENO1","ENO2","ENO3","ERO1A","ERRFI1","ETS1","EXT1","F3","FAM162A","FBP1","FOS","FOSL2","FOXO3","GAA","GALK1","GAPDH","GAPDHS","GBE1","GCK","GCNT2","GLRX","GPC1","GPC3","GPC4","GPI","GRHPR","GYS1","HAS1","HDLBP","HEXA","HK1","HK2","HMOX1","HOXB9","HS3ST1","HSPA5","IDS","IER3","IGFBP1","IGFBP3","IL6","ILVBL","INHA","IRS2","ISG20","JMJD6","JUN","KDELR3","KDM3A","KIF5A","KLF6","KLF7","KLHL24","LALBA","LARGE1","LDHA","LDHC","LOX","LXN","MAFF","MAP3K1","MIF","MT1E","MT2A","MXI1","MYH9","NAGK","NCAN","NDRG1","NDST1","NDST2","NEDD4L","NFIL3","NOCT","NR3C1","P4HA1","P4HA2","PAM","PCK1","PDGFB","PDK1","PDK3","PFKFB3","PFKL","PFKP","PGAM2","PGF","PGK1","PGM1","PGM2","PHKG1","PIM1","PKLR","PKP1","PLAC8","PLAUR","PLIN2","PNRC1","PPARGC1A","PPFIA4","PPP1R15A","PPP1R3C","PRDX5","PRKCA","PYGM","RBPJ","RORA","RRAGD","S100A4","SAP30","SCARB1","SDC2","SDC3","SDC4","SELENBP1","SERPINE1","SIAH2","SLC25A1","SLC2A1","SLC2A3","SLC2A5","SLC37A4","SLC6A6","SRPX","STBD1","STC1","STC2","SULT2B1","TES","TGFB3","TGFBI","TGM2","TIPARP","TKTL1","TMEM45A","TNFAIP3","TPBG","TPD52","TPI1","TPST2","UGP2","VEGFA","VHL","VLDLR","WSB1","XPNPEP1","ZFP36","ZNF292")

emt.hallmark <- c("ABI3BP","ACTA2","ADAM12","ANPEP","APLP1","AREG","BASP1","BDNF","BGN","BMP1","CADM1","CALD1","CALU","CAP2","CAPG","CCN1","CCN2","CD44","CD59","CDH11","CDH2","CDH6","COL11A1","COL12A1","COL16A1","COL1A1","COL1A2","COL3A1","COL4A1","COL4A2","COL5A1","COL5A2","COL5A3","COL6A2","COL6A3","COL7A1","COL8A2","COLGALT1","COMP","COPA","CRLF1","CTHRC1","CXCL1","CXCL12","CXCL6","CXCL8","DAB2","DCN","DKK1","DPYSL3","DST","ECM1","ECM2","EDIL3","EFEMP2","ELN","EMP3","ENO2","FAP","FAS","FBLN1","FBLN2","FBLN5","FBN1","FBN2","FERMT2","FGF2","FLNA","FMOD","FN1","FOXC2","FSTL1","FSTL3","FUCA1","FZD8","GADD45A","GADD45B","GAS1","GEM","GJA1","GLIPR1","GPC1","GPX7","GREM1","HTRA1","ID2","IGFBP2","IGFBP3","IGFBP4","IL15","IL32","IL6","INHBA","ITGA2","ITGA5","ITGAV","ITGB1","ITGB3","ITGB5","JUN","LAMA1","LAMA2","LAMA3","LAMC1","LAMC2","LGALS1","LOX","LOXL1","LOXL2","LRP1","LRRC15","LUM","MAGEE1","MATN2","MATN3","MCM7","MEST","MFAP5","MGP","MMP1","MMP14","MMP2","MMP3","MSX1","MXRA5","MYL9","MYLK","NID2","NNMT","NOTCH2","NT5E","NTM","OXTR","P3H1","PCOLCE","PCOLCE2","PDGFRB","PDLIM4","PFN2","PLAUR","PLOD1","PLOD2","PLOD3","PMEPA1","PMP22","POSTN","PPIB","PRRX1","PRSS2","PTHLH","PTX3","PVR","QSOX1","RGS4","RHOB","SAT1","SCG2","SDC1","SDC4","SERPINE1","SERPINE2","SERPINH1","SFRP1","SFRP4","SGCB","SGCD","SGCG","SLC6A8","SLIT2","SLIT3","SNAI2","SNTB1","SPARC","SPOCK1","SPP1","TAGLN","TFPI2","TGFB1","TGFBI","TGFBR3","TGM2","THBS1","THBS2","THY1","TIMP1","TIMP3","TNC","TNFAIP3","TNFRSF11B","TNFRSF12A","TPM1","TPM2","TPM4","VCAM1","VCAN","VEGFA","VEGFC","VIM","WIPF1","WNT5A")

infg.response.hallmark <- c("ADAR", "APOL6", "ARID5B", "ARL4A", "AUTS2", "B2M", "BANK1", "BATF2", "BPGM", "BST2", "BTG1", "C1R", "C1S", "CASP1", "CASP3", "CASP4", "CASP7", "CASP8", "CCL2", "CCL5", "CCL7", "CD274", "CD38", "CD40", "CD69", "CD74", "CD86", "CDKN1A", "CFB", "CFH", "CIITA", "CMKLR1", "CMPK2", "CMTR1", "CSF2RB", "CXCL10", "CXCL11", "CXCL9", "DDX58", "DDX60", "DHX58", "EIF2AK2", "EIF4E3", "EPSTI1", "FAS", "FCGR1A", "FGL2", "FPR1", "GBP4", "GBP6", "GCH1", "GPR18", "GZMA", "HELZ2", "HERC6", "HIF1A", "HLA-A", "HLA-B", "HLA-DMA", "HLA-DQA1", "HLA-DRB1", "HLA-G", "ICAM1", "IDO1", "IFI27", "IFI30", "IFI35", "IFI44", "IFI44L", "IFIH1", "IFIT1", "IFIT2", "IFIT3", "IFITM2", "IFITM3", "IFNAR2", "IL10RA", "IL15", "IL15RA", "IL18BP", "IL2RB", "IL4R", "IL6", "IL7", "IRF1", "IRF2", "IRF4", "IRF5", "IRF7", "IRF8", "IRF9", "ISG15", "ISG20", "ISOC1", "ITGB7", "JAK2", "KLRK1", "LAP3", "LATS2", "LCP2", "LGALS3BP", "LY6E", "LYSMD2", "MARCHF1", "METTL7B", "MT2A", "MTHFD2", "MVP", "MX1", "MX2", "MYD88", "NAMPT", "NCOA3", "NFKB1", "NFKBIA", "NLRC5", "NMI", "NOD1", "NUP93", "OAS2", "OAS3", "OASL", "OGFR", "P2RY14", "PARP12", "PARP14", "PDE4B", "PELI1", "PFKP", "PIM1", "PLA2G4A", "PLSCR1", "PML", "PNP", "PNPT1", "PSMA2", "PSMA3", "PSMB10", "PSMB2", "PSMB8", "PSMB9", "PSME1", "PSME2", "PTGS2", "PTPN1", "PTPN2", "PTPN6", "RAPGEF6", "RBCK1", "RIPK1", "RIPK2", "RNF213", "RNF31", "RSAD2", "RTP4", "SAMD9L", "SAMHD1", "SECTM1", "SELP", "SERPING1", "SLAMF7", "SLC25A28", "SOCS1", "SOCS3", "SOD2", "SP110", "SPPL2A", "SRI", "SSPN", "ST3GAL5", "ST8SIA4", "STAT1", "STAT2", "STAT3", "STAT4", "TAP1", "TAPBP", "TDRD7", "TNFAIP2", "TNFAIP3", "TNFAIP6", "TNFSF10", "TOR1B", "TRAFD1", "TRIM14", "TRIM21", "TRIM25", "TRIM26", "TXNIP", "UBE2L6", "UPP1", "USP18", "VAMP5", "VAMP8", "VCAM1", "WARS1", "XAF1", "XCL1", "ZBP1", "ZNFX1")

DefaultAssay(myeloid) = "RNA"
myeloid <- NormalizeData(myeloid)
myeloid <- AddModuleScore(myeloid, list(hypoxia.hallmark), name = "hypoxia.hallmark")
myeloid <- AddModuleScore(myeloid, list(emt.hallmark), name = "emt.hallmark")
myeloid <- AddModuleScore(myeloid, list(infg.response.hallmark), name = "infg.response.hallmark")

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
pdf(paste0(figures.dir, "F4E.pdf"), useDingbats = F, width = 4, height = 4)
pheatmap(avglog2fc.df.plot, display_numbers = pvals.df.plot, breaks = seq(-range, range, length.out = 100),  cellheight=13,cellwidth=13, 
         cluster_rows = T,cluster_cols = F
)
dev.off()


#### F4F ####
df <- get(load(paste0(data.path, "lig_rec_differential_highlvl_df_not_allzeros.Rda")))

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
pdf(file = paste0(figures.dir, "F4F.pdf"), useDingbats = F, height = 15, width = 8.3)
p<-ggplot(data = df2, mapping = aes(x=cell.pair, y=int_name, color=p53_wt_prop_diff, size=p53_wt_neglogpval_proptest))+geom_point(shape = 21, aes(colour = as.factor(signif), fill = p53_wt_prop_diff))+
  scale_colour_manual(values=c("00FFFFFF", "darkgray", "black")) + theme_minimal() +
  custom_color_scale +
  theme(text = element_text(size=12), strip.text = element_text(size=22,face='bold'))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
p
dev.off()




### F4H ####
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
pdf(paste0(figures.dir, "F4H.pdf"), useDingbats = F, width = 20, height = 8) 
pheatmap(avglog2fc.df, display_numbers = pvals.df, breaks = seq(-range, range, length.out = 100),  cellheight=13,cellwidth=13, 
         cluster_rows = F,cluster_cols = F,
         color = colorRampPalette(c("navy", "white", "firebrick"))(100))
dev.off()



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
