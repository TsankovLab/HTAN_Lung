library (Seurat)
library(pals)
library(RColorBrewer)
library(ggplot2)
library(scales)
library(dplyr)

# # Load Lymphoid object
# figures.dir <- paste0(proj.path, "/Results/10x_All_190615/Immune.LymTnG400v11/")
# load(file = paste0(figures.dir, "tcell.final.nodoublets.Rda")) # lymphoid

### Set paths ####
data.path <- file.path('..','data')
figures.dir <- "../../figures/"

### Load functions ####
source(file.path(data.path, "R_utils','plotutils.R"))
source(file.path(data.path, "R_utils','seuratutils.R"))
source(file.path(data.path, "R_utils','seuratutilsV3.R"))
source(file.path(data.path, "R_utils','color.R"))

### Load data ####
tcell = readRDS (file.path ('./dropbox_data','NKTcells.rds'))
load(file = file.path ('./dropbox_data', 'bcell.final.nodoublets.Rda')) # bmast
tnk <- get(load(file.path('./dropbox_data','tnk.validation.Rda'))) # validation cohort



#### F5AB ####
dp = DimPlot (tcell, group.by = 'subtype')
dp
ggsave(file.path(figures.dir, "F5A.png"), width = 7, height = 7)

colors <- material.heat(11)
colors <- jet.colors(100)
colors<-c("#191970","#121285","#0C0C9A","#0707B0","#0101C5","#0014CF","#0033D3","#0053D8","#0072DD","#0092E1","#00B2E6",
            "#00D1EB","#23E8CD","#7AF17B","#D2FA29","#FFEB00","#FFC300","#FF9B00","#FF8400","#FF7800","#FF6B00","#FF5F00","#FF5300",
            "#FF4700","#F73B00","#EF2E00","#E62300","#DD1700","#D50B00","#CD0000") # Carolyn Porter's color palette


pal <- gradient_n_pal(rev(c("#D53E4F", "#FC8D59", "#FEE08B", "#3288BD", "#0060a4")))
genes <- c("CD8A", "CD4", "MKI67", "GNLY", "FOXP3", "HAVCR2", "CTLA4", "PDCD1", "TIGIT")
data <- tcell[,sample(colnames(tcell), 500)][["RNA"]]@data[genes,]
data <- as.matrix(data)

data <- as.matrix(tcell[["RNA"]]@data[genes,])

p <- FeaturePlot(tcell, features = c("CD8A", "CD4", "MKI67", "GNLY", "FOXP3", "HAVCR2", "CTLA4", "PDCD1", "TIGIT"), combine = FALSE, pt.size = 0.001)
# custom_color_scale <- scale_color_gradientn(colours = colors, limits = c(0,4))

custom_color_scale <- scale_color_gradientn(
    colours = pal(c(0, rescale(seq_along(data)), 1)), # <- extra 0, 1 for out-of-bounds
    limits = c(0, 6), breaks = 0:6,
    na.value = "#9E0142",
    values = c(0,rescale(seq_along(data)),1) # <- extra 0, 1 again # also try changing to 7
) ### what is different about this?

p <- lapply(p, function (x) x + custom_color_scale)
for(i in 1:length(p)) {p[[i]] <- p[[i]] + NoLegend() + NoAxes()}


cowplot::plot_grid(plotlist = p, ncol=3)
ggsave(file.path(figures.dir, "F5B.png"), width = 15, height = 15)

### feature plots are not space efficient
#write.csv (do.call (cbind, lapply(p,function(x) x$data)), file.path (projdir, 'F5_A_B.csv'))

#### F5C_S8C ####
p53 <- c("1172", "1173", "1175", "1182", "BW16", "BW11", "BW14", "BW06")
WT <- c("14", "1174", "1176", "1183", "BW01", "BW04", "BW05", "BW09", "BW19", "BW23")
p53.filt <- c("1172", "1173", "1175", "1182", "BW16", "BW11", "BW14", "BW06")
WT.filt <- c("14", "1174", "1176", "1183", "BW01", "BW04", "BW05", "BW09", "BW19", "BW23") # no 1170

comp.df_disc <- as.matrix(table(tcell$orig.identSec, tcell$subtype))
p53.status.filt <- c(rep("WT", 10), rep("p53_mut", 8))# no 1170
comp.df_disc <- comp.df_disc[c(WT.filt, p53.filt), ]
comp.df_disc  <- prop.table(comp.df_disc , 1)
comp.df_disc <- as.data.frame.matrix(comp.df_disc)
comp.df_disc$p53_status <- p53.status.filt
comp.df_disc$p53_status <- factor(comp.df_disc$p53_status, levels = c("WT", "p53_mut"))

subtypes <- colnames(comp.df_disc)


mean.plots <- list()

for (subtype in subtypes){
  mean.plot <- ggboxplot(comp.df_disc, x = "p53_status", y = subtype, color = "p53_status",
               palette = c("#0091CA", "#D8423D"), add = "jitter",
               title = subtype) +
               # ylim(-.5,.5) +
               stat_compare_means(label.x = 0.9, label.y = .9 * max(as.numeric(comp.df_disc[,subtype]))) +
               # stat_compare_means(label.x = 0.9, label.y = max(module.scores[,module])*(9.0/10), method = "t.test") +
               theme(plot.title = element_text(size = 12, face = "bold")) + NoLegend()
  #ggsave(paste0(figures.dir.p53, "median.module.scores.v3.top20.png"), width = 6, height = 6, type="cairo")
  mean.plots[[subtype]] <- mean.plot
}

            
pdf(file.path(figures.dir,'F5C_S8C.pdf'), width = 18, height = 7.5, useDingbats=FALSE)
cowplot::plot_grid(plotlist = mean.plots, ncol=8)
dev.off()

#write.csv (do.call (cbind, lapply (mean.plots, function(x) x$data)), file.path(projdir, 'F5_C_S8_C.csv'))



#### S8C ####    

comp.df <- as.matrix(table(tnk$sampleID, tnk$predicted.id))
comp.df  <- prop.table(comp.df , 1)
comp.df <- as.data.frame.matrix(comp.df)
comp.df$p53_status = tnk$p53_status[match(rownames(comp.df), tnk$sampleID)]
comp.df$p53_status <- factor(comp.df$p53_status, levels = c("WT", "mut"))

subtypes <- colnames(comp.df)
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
pdf(file.path(figures.dir,'S8_C_validation_cohort.pdf'), width = 18, height = 7.5, useDingbats=FALSE)
cowplot::plot_grid(plotlist = mean.plots, ncol=8)
dev.off()

#write.csv (do.call (cbind, lapply (mean.plots[c('T.Exhausted','TFH')], function(x) x$data)), file.path (projdir, 'S8_C.csv'))

comp.df = rbind (comp.df_disc[,c('p53_status','TFH')], comp.df[,c('p53_status','TFH')])
comp.df$p53_status[comp.df$p53_status == 'mut']= 'p53_mut'


subtype= 'TFH'
  mean.plot <- ggboxplot(comp.df, x = "p53_status", y = subtype, color = "p53_status",
               palette = c("#0091CA", "#D8423D"), add = "jitter",
               title = subtype) +
               # ylim(-.5,.5) +
               stat_compare_means(label.x = 0.9, label.y = .9 * max(as.numeric(comp.df[,subtype]))) +
               # stat_compare_means(label.x = 0.9, label.y = max(module.scores[,module])*(9.0/10), method = "t.test") +
               theme(plot.title = element_text(size = 12, face = "bold")) + NoLegend()
  #ggsave(paste0(figures.dir.p53, "median.module.scores.v3.top20.png"), width = 6, height = 6, type="cairo")
  #mean.plots[[subtype]] <- mean.plot

        
pdf(file.path(figures.dir,'S8C.pdf'), width = 3, height = 3.5)
mean.plot
dev.off()

#write.csv (mean.plot$data, file.path (projdir, 'S8_C_cohort_TFH_combined.csv'))



#### F5E ####   
names <- c("TIGIT+CTLA4-PDCD1-",
  "TIGIT+CTLA4+PDCD1-",
  "TIGIT+CTLA4-PDCD1+",
  "TIGIT+CTLA4+PDCD1+",
  "TIGIT-CTLA4+PDCD1-",
  "TIGIT-CTLA4-PDCD1+",
  "TIGIT-CTLA4+PDCD1+")
  # "CTLA-PDCD1-TIGIT-")
tcell$CTLA4 <- tcell[["RNA"]]@counts["CTLA4",] > 0
tcell$PDCD1 <- tcell[["RNA"]]@counts["PDCD1",] > 0
tcell$TIGIT <- tcell[["RNA"]]@counts["TIGIT",] > 0
tcell@meta.data[,"TIGIT+CTLA4+PDCD1+"] <- (tcell$TIGIT==T & tcell$CTLA4 ==T & tcell$PDCD1==T)
tcell@meta.data[,"TIGIT+CTLA4+PDCD1-"] <- (tcell$TIGIT==T & tcell$CTLA4 ==T & tcell$PDCD1==F)
tcell@meta.data[,"TIGIT+CTLA4-PDCD1+"] <- (tcell$TIGIT==T & tcell$CTLA4 ==F & tcell$PDCD1==T)
tcell@meta.data[,"TIGIT+CTLA4-PDCD1-"] <- (tcell$TIGIT==T & tcell$CTLA4 ==F & tcell$PDCD1==F)
tcell@meta.data[,"TIGIT-CTLA4+PDCD1+"] <- (tcell$TIGIT==F & tcell$CTLA4 ==T & tcell$PDCD1==T)
tcell@meta.data[,"TIGIT-CTLA4+PDCD1-"] <- (tcell$TIGIT==F & tcell$CTLA4 ==T & tcell$PDCD1==F)
tcell@meta.data[,"TIGIT-CTLA4-PDCD1+"] <- (tcell$TIGIT==F & tcell$CTLA4 ==F & tcell$PDCD1==T)
tcell@meta.data[,"TIGIT-CTLA4-PDCD1-"] <- (tcell$TIGIT==F & tcell$CTLA4 ==F & tcell$PDCD1==F)

tcell$p53_status <- "WT"
tcell$p53_status[tcell$orig.identSec %in% p53.filt] <- "p53_mut"
#tcell$p53_status[tcell$orig.identSec %in% other] <- "other"

Idents(tcell) <- "p53_status"
tcell.p53 <- subset(tcell, idents = "p53_mut")
tcell.wt <- subset(tcell, idents = "WT")
sums.p53 <- c()
sums.wt <- c()
for (name in names) {
  sums.p53  <- c(sums.p53 , (sum(tcell.p53@meta.data[,name])/(nrow(tcell.p53@meta.data))))
  sums.wt <- c(sums.wt , (sum(tcell.wt@meta.data[,name])/(nrow(tcell.wt@meta.data))))
}
sums.p53
sums.wt

sums <- cbind(sums.wt, sums.p53)
colnames(sums) <- c("WT", "p53_mut")
rownames(sums) <- names
sums <- as.data.frame(sums)
sums <- sums[order(sums$p53_mut),]
sums <- sums[rev(names),]

sequential <- brewer.pal(7, "Spectral")
         
pdf(file.path(figures.dir, "F5E.pdf"), useDingbats = F, width = 5, height = 6)
barplot(as.matrix(sums),
 # names.arg = wide$Seed, # x-axis labels
 cex.names = 1, # makes x-axis labels small enough to show all
 col = sequential, # colors
 las=2,
 # xlab = "Seed Source",
 # ylab = "Height, Feet",
 xlim = c(0,10), # these two lines allow space for the legend
 ylim = c(0,0.5),
 width = 1) # these two lines allow space for the legend
legend("bottomright",
 legend = rev(rownames(sums)), #in order from top to bottom
 fill = sequential[7:1])#, # 6:1 reorders so legend order matches graph
 #title = "Years")
dev.off()

#write.csv (sums, file.path (projdir, 'F5_E.csv'))


#### F5H ####            
all_malig = data.frame (
        row.names = c(
        'PVR+CD274+',
        'PVR+CD274-',
        'PVR-CD274+',
        'PVR-CD274-'),WT=c(0.004290886,
     0.072632191,
     0.024315023,
     0.898761901),p53_mut=c(0.03597884,
        0.15806878,
        0.05740741,
        0.74854497))
cycling_malig = data.frame (
        row.names = c(
        'PVR+CD274+',
        'PVR+CD274-',
        'PVR-CD274+',
        'PVR-CD274-'),
        WT=c(
        0.01237964,
        0.10729023,
        0.03507565,
        0.84525447),
        p53_mut=c(
        0.1402525,
        0.3029453,
        0.0743338,
        0.4824684))

sums <- all_malig

rownames(sums) <- sums$X
sums$X <- NULL
sequential <- brewer.pal(4, "Spectral")
# sequential <- rev(sequential)

#### F5H malignant cells ####            
pdf (file.path(figures.dir,'F5H_malignant_cells.pdf'))
barplot(as.matrix(sums),
        # names.arg = wide$Seed, # x-axis labels
        cex.names = 1, # makes x-axis labels small enough to show all
        col = sequential, # colors
        las=2,
        # xlab = "Seed Source",
        ylab = "Proportion of cycling cancer cells",
        xlim = c(0,10), # these two lines allow space for the legend
        ylim = c(0,1),
        width = 1, # these two lines allow space for the legend
        legend = rownames(sums), #in order from top to bottom
        fill = sequential[4:1],
        )#, # 6:1 reorders so legend order matches graph
#title = "Years")
dev.off()

#write.csv (sums, file.path (projdir, 'F5_H_malignant_cells.csv'))

#### F5H cycling cells ####            
sums <- cycling_malig
pdf (file.path(figures.dir,'F5H_cycling_cells.pdf'))
barplot(as.matrix(sums),
        # names.arg = wide$Seed, # x-axis labels
        cex.names = 1, # makes x-axis labels small enough to show all
        col = sequential, # colors
        las=2,
        # xlab = "Seed Source",
        ylab = "Proportion of cycling cancer cells",
        xlim = c(0,10), # these two lines allow space for the legend
        ylim = c(0,1),
        width = 1, # these two lines allow space for the legend
        legend = rownames(sums), #in order from top to bottom
        fill = sequential[4:1],
        )#, # 6:1 reorders so legend order matches graph
#title = "Years")
dev.off()

#write.csv (sums, file.path (projdir, 'F5_H_cycling_cells.csv'))

#### S8A cycling cells ####            
Idents(tcell) <- "subtype"
markers <- FindAllMarkers(tcell, only.pos = T, max.cells.per.ident = 1000)

markers[,"cluster"] <- as.character(markers[,"cluster"]) # in case this is a factor

markers <- markers %>%
  arrange(cluster, desc(avg_log2FC)) %>% 
  group_by(cluster) %>%
  dplyr::slice(1:10)
  

# dir.create(paste0(figures.dir, "Figures_v4/Fig5_myeloid_bmast/"))
#write.csv(markers, file = file.path(projdir, "Figures_v4/Fig6_tnk/tnk.markers.top10_rep.csv"))
#markers <- read.csv(file = file.path(projdir, "Figures_v4/Fig6_tnk/tnk.markers.top10_rep.csv"))

#order <- rev(unique(markers$cluster)) # already sorted

#Idents(tcell) <- factor(Idents(tcell), levels = order)

genes <- c(
  "IL7R", 
  "CCR7", 
  "S100A4", 
  "CD52", 
  "GZMK", 
  "EOMES",
  "ISG15", 
  "IFI6", 
  "ZNF683",
  "CD8B",
  "KLRC1",
  "AREG",
  "FGFBP2",
  "SPON2",
  "KLRC2",
  "STMN1", 
  "MKI67", 
  "LAG3", 
  "HAVCR2", 
  "HSPA1A", 
  "JUN", 
  "CXCL13", 
  "NR3C1", 
  "IL32", 
  "FOXP3")
#genes <- unique(markers$gene)

tcell <- FindVariableFeatures(tcell)
tcell <- ScaleData(tcell, features = unique(c(VariableFeatures(tcell), genes))) # add genes of interest to the heatmap

data <- tcell[,sample(colnames(tcell), 500)][["RNA"]]@data[genes,]
data <- as.matrix(data)

pdf(file.path(figures.dir,"S8A.pdf"), useDingbats = F, width = 15, height = 6)
dp = DotPlot(object = tcell, features = genes, scale = T, assay = "RNA") +
  theme(axis.text.x = element_text(angle = 45, hjust=1), panel.border = element_rect(colour = "black", fill=NA, size=0.5), panel.grid.major = element_line(colour = "gainsboro")) + 
    scale_color_gradientn(colours = rev(brewer.pal(11,"Spectral"))) +
    geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.5)
dp
dev.off()

#write.csv (dp$data, file.path (projdir, 'S8_A.csv'))

#### S8B ####            
table (bcell$subtype)
bcell = bcell [, !bcell$subtype %in% c('Mast','Plasmacytoid.DC')]

pdf (file.path (figures.dir,'S8B.pdf'))
dp = DimPlot (bcell, group.by = 'subtype')
dp
dev.off()
#write.csv (dp$data, file.path (projdir, 'S8_B'))



#### F8E ####      
genes <- c(
  "CXCL13", 
  "KLRB1", 
  "KLRC1", 
  "ITGB2", 
  "TGFB1", 
  "CTLA4",
  "IFNG", 
  "CXCR3", 
  "TNF",
  "ITGAL",
  "ADGRE5",
  "CD96",
  "CD47",
  "CLEC2D",
  "CD44",
  "SELL", 
  "LAMP1", 
  "TNFRSF14", 
  "TIGIT", 
  "HAVCR2", 
  "PDCD1")
#genes <- unique(markers$gene)

tcell <- FindVariableFeatures(tcell)
tcell <- ScaleData(tcell, features = unique(c(VariableFeatures(tcell), genes))) # add genes of interest to the heatmap

data <- tcell[,sample(colnames(tcell), 500)][["RNA"]]@data[genes,]
data <- as.matrix(data)
Idents(tcell) = tcell$subtype
Idents(tcell) = factor (Idents(tcell), levels = rev (c('NK.CD56.dim','NK.CD56.bright','T.Stress','CD4.Naive.CM','CD8.GZMK','T.Cycling','T.Exhausted','Treg','TFH','NKT','CD8.TRM','CD8.IFN','CD4.TRM')))
      
pdf(file.path(figures.dir,"S8E.pdf"), useDingbats = F, width = 10, height = 5)
 dp = DotPlot(object = tcell, features = genes, scale = T, assay = "RNA") +
  theme(axis.text.x = element_text(angle = 45, hjust=1), panel.border = element_rect(colour = "black", fill=NA, size=0.5), panel.grid.major = element_line(colour = "gainsboro")) + 
    scale_color_gradientn(colours = rev(brewer.pal(11,"Spectral"))) +
    geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.5)
    dp
dev.off()

#write.csv (dp$data, file.path (projdir, 'S8_E.csv'))



