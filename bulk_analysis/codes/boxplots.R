######################################################
#################### TSANKOV LAB #####################
######################################################
# ----------------------------------------------------
# load the libraries
# ----------------------------------------------------
library(Seurat)
library(tidyverse)
library(ggpubr)
# ----------------------------------------------------
# set directories
# ----------------------------------------------------
data.path <- "./data/"
figures.dir <- "../../figures/"
system(paste0("mkdir -p ", figures.dir))
# ----------------------------------------------------
# load the objects
# ----------------------------------------------------
luad.2018 <- readRDS(file.path(data.path, "luad.2018.processed.rds"))
cptac.luad.prot <- readRDS(file.path(data.path, "cptac.luad.prot.processed.rds"))
# ----------------------------------------------------
# load the dataframe
# ----------------------------------------------------
df.pvals.tcga <- read_csv(file.path(data.path, "tcga.regression_plots.p53_kras_egfr_stk11_CancerProp.filtered.csv")) |>
    select(-"...1")
# ----------------------------------------------------
# custom function
# ----------------------------------------------------
createBoxplots <- function(plotVector, 
                           plotDF, 
                           outfileName, 
                           xVar = "p53_status", 
                           comparisonList = list(c("WT", "mut")), 
                           df.pvals = NULL,
                           numCol = 4, 
                           w = 3.5, 
                           h = 3) {
    plotlist <- list()

    for (plotVar in plotVector) {
        plotVar.clean <- str_remove_all(string = plotVar, pattern = "`")
        if (!plotVar.clean %in% colnames(plotDF)) {
            message(paste0(plotVar.clean, " not found in plotDF!"))
            next
        }
        message(paste0("Adding plot and saving csv file for ", plotVar, "..."))
        saveDF <- plotDF |>
            drop_na(all_of(xVar)) |>
            select(all_of(c(xVar, plotVar.clean)))
        # write.csv(saveDF, file.path(paste0(data.path, outfileName, ".", make.names(plotVar.clean), ".csv")))
        
        if (length(levels(plotDF[[xVar]])) == 2) {
            
            if (is.null(df.pvals)) {
                # for cptac
                plotlist[[plotVar]] <- plotDF |>
                    drop_na(all_of(xVar)) |>
                    ggboxplot(x = xVar, 
                              y = plotVar,
                              color = xVar, 
                              palette = c("#0091CA", "#D8423D"),
                              title = plotVar) +
                    stat_compare_means(comparisons = comparisonList) +
                    scale_x_discrete(labels = function(x) paste0(x, "\n(n = ", table(plotDF[[xVar]])[x], ")")) +
                    rotate_x_text(0) +
                    xlab(NULL) +
                    theme(plot.title = element_text(size = 12, face = "bold"), legend.position = "none")
            } else {
                # for tcga
                pvals_vector <- setNames(df.pvals$pval, df.pvals$program)
                pvals_manual <- pvals_vector[[plotVar.clean]]
                plotlist[[plotVar]] <- plotDF |>
                    drop_na(all_of(xVar)) |>
                    ggboxplot(x = xVar, 
                              y = plotVar,
                              color = xVar, 
                              palette = c("#0091CA", "#D8423D"),
                              title = plotVar) +
                    # stat_compare_means(comparisons = comparisonList) +
                    scale_x_discrete(labels = function(x) paste0(x, "\n(n = ", table(plotDF[[xVar]])[x], ")")) +
                    rotate_x_text(0) +
                    xlab(NULL) +
                    theme(plot.title = element_text(size = 12, face = "bold"), legend.position = "none") +
                    ggsignif::geom_signif(comparisons = list(c("WT", "mut")),  # Replace "WT" and "mut" with actual levels if different
                                          annotations = paste0("p = ", signif(pvals_manual, digits = 3)),
                                          y_position = max(plotDF[[plotVar.clean]], na.rm = TRUE) * 1.1,  # Adjust as needed
                                          tip_length = 0.02,  # Adjust the length of the horizontal line tips
                                          textsize = 3)
            }
        } else {
            plotlist[[plotVar]] <- plotDF |>
                drop_na(all_of(xVar)) |>
                ggboxplot(x = xVar, 
                          y = plotVar,
                          color = xVar, 
                          title = plotVar) +
                stat_compare_means(comparisons = comparisonList) +
                scale_x_discrete(labels = function(x) paste0(x, "\n(n = ", table(plotDF[[xVar]])[x], ")")) +
                rotate_x_text(90) +
                xlab(NULL) +
                theme(plot.title = element_text(size = 12, face = "bold"), legend.position = "none")
        }
    }
    
    ggsave(filename = file.path(paste0(figures.dir, outfileName, ".pdf")), 
           plot = patchwork::wrap_plots(plotlist, ncol = numCol), 
           device = "pdf", 
           width = w * numCol, 
           height = h * ceiling(length(plotlist) / 4))
}

# ----------------------------------------------------
# generate figures
# ----------------------------------------------------
# FIG 1F
v.1f <- c("Endothelial", "Pericytes")

createBoxplots(plotVector = v.1f, w = 3, h = 4, df.pvals = df.pvals.tcga,
               xVar = "p53_status", 
               comparisonList = list(c("WT", "mut")),
               numCol = 2, 
               plotDF = luad.2018@meta.data, 
               outfileName = "1f")

# FIG 2E
v.2e <- c("p53_targets", "AT2", "CC.G2M", "`Glycolysis/Hypox`", "pEMT")

createBoxplots(plotVector = v.2e, df.pvals = df.pvals.tcga,
               xVar = "p53_status", 
               comparisonList = list(c("WT", "mut")),
               numCol = 5, 
               plotDF = luad.2018@meta.data, 
               outfileName = "2e")

# FIG 4D
v.4d <- c("SPP1")

createBoxplots(plotVector = v.4d, w = 3, h = 4, df.pvals = df.pvals.tcga,
               xVar = "p53_status", 
               comparisonList = list(c("WT", "mut")),
               numCol = 1, 
               plotDF = luad.2018@meta.data, 
               outfileName = "4d")

# EXTENDED FIG 1D
v.s1d <- c("B.cells", "Fibroblast", "Cancer", "Mast", "`T.CD4+`", # row 1
           "Myeloid", "NK.cells", "Smooth.muscle", "T.cells", "`T.CD8+`") # row 2

createBoxplots(plotVector = v.s1d, df.pvals = df.pvals.tcga,
               xVar = "p53_status", 
               comparisonList = list(c("WT", "mut")),
               numCol = 5, 
               plotDF = luad.2018@meta.data, 
               outfileName = "s1d")

# EXTENDED FIG 4A,C,D
v.s4acd <- c("p53_targets", "AT2", "CC.G2M", "`Glycolysis/Hypox`", "pEMT") # row 1

levels.tmp <- levels(luad.2018@meta.data$comut)
createBoxplots(plotVector = v.s4acd, 
               xVar = "comut", 
               comparisonList = lapply(setdiff(levels.tmp, "WT"), function(x) c("WT", x))[-6],
               numCol = 5, 
               plotDF = luad.2018@meta.data, 
               outfileName = "s4a")

levels.tmp <- levels(luad.2018@meta.data$DNE_LOFclass)
createBoxplots(plotVector = v.s4acd, 
               xVar = "DNE_LOFclass", 
               comparisonList = lapply(setdiff(levels.tmp, "WT"), function(x) c("WT", x)),
               numCol = 5, 
               plotDF = luad.2018@meta.data, 
               outfileName = "s4c")

levels.tmp <- levels(luad.2018@meta.data$impact)
createBoxplots(plotVector = v.s4acd, 
               xVar = "impact", 
               comparisonList = lapply(setdiff(levels.tmp, "WT"), function(x) c("WT", x))[-3],
               numCol = 5, 
               plotDF = luad.2018@meta.data, 
               outfileName = "s4d")

# EXTENDED FIG 4G
v.s4g <- c("p53_targets", "AT2", "CC.G2M", "Glycolysis.Hypox", "pEMT")
createBoxplots(plotVector = v.s4g, 
               xVar = "p53_status", 
               comparisonList = list(c("WT", "mut")),
               numCol = 5, 
               plotDF = cptac.luad.prot@meta.data, 
               outfileName = "s4g")

# EXTENDED FIG 4H
v.s4h <- c("p53_targets", "AT2", "CC.G2M")

levels.tmp <- levels(cptac.luad.prot@meta.data$DNE_LOFclass)
createBoxplots(plotVector = v.s4h, w = 3.5, h = 5,
               xVar = "DNE_LOFclass", 
               comparisonList = lapply(setdiff(levels.tmp, "WT"), function(x) c("WT", x)),
               numCol = 3, 
               plotDF = cptac.luad.prot@meta.data, 
               outfileName = "s4hA")

levels.tmp <- levels(cptac.luad.prot@meta.data$impact)
createBoxplots(plotVector = v.s4h, w = 3.5, h = 5,
               xVar = "impact", 
               comparisonList = lapply(setdiff(levels.tmp, "WT"), function(x) c("WT", x))[-3],
               numCol = 3, 
               plotDF = cptac.luad.prot@meta.data, 
               outfileName = "s4hB")

# EXTENDED FIG 6A
v.s6a <- c("Endothelial")

levels.tmp <- levels(luad.2018@meta.data$comut)
createBoxplots(plotVector = v.s6a, w = 3, h = 5,
               xVar = "comut", 
               comparisonList = lapply(setdiff(levels.tmp, "WT"), function(x) c("WT", x))[-6],
               numCol = 1, 
               plotDF = luad.2018@meta.data, 
               outfileName = "s6aA")

levels.tmp <- levels(luad.2018@meta.data$DNE_LOFclass)
createBoxplots(plotVector = v.s6a, w = 3, h = 5,
               xVar = "DNE_LOFclass", 
               comparisonList = lapply(setdiff(levels.tmp, "WT"), function(x) c("WT", x)),
               numCol = 1, 
               plotDF = luad.2018@meta.data, 
               outfileName = "s6aB")

# EXTENDED FIG 6E
v.s6e <- c("Aerocyte", "Arterial")

createBoxplots(plotVector = v.s6e, w = 3, h = 4, df.pvals = df.pvals.tcga,
               xVar = "p53_status", 
               comparisonList = list(c("WT", "mut")),
               numCol = 2, 
               plotDF = luad.2018@meta.data, 
               outfileName = "s6e")

# EXTENDED FIG 7G
v.s7g <- c("CAF.ADH1B")

createBoxplots(plotVector = v.s7g, w = 3, h = 4, df.pvals = df.pvals.tcga,
               xVar = "p53_status", 
               comparisonList = list(c("WT", "mut")),
               numCol = 1, 
               plotDF = luad.2018@meta.data, 
               outfileName = "s7g")

# EXTENDED FIG 7H
v.s7h <- c("Pericytes")

levels.tmp <- levels(luad.2018@meta.data$comut)
createBoxplots(plotVector = v.s7h, w = 3, h = 5, 
               xVar = "comut", 
               comparisonList = lapply(setdiff(levels.tmp, "WT"), function(x) c("WT", x))[-6],
               numCol = 1, 
               plotDF = luad.2018@meta.data, 
               outfileName = "s7hA")

levels.tmp <- levels(luad.2018@meta.data$DNE_LOFclass)
createBoxplots(plotVector = v.s7h, w = 3, h = 5,
               xVar = "DNE_LOFclass", 
               comparisonList = lapply(setdiff(levels.tmp, "WT"), function(x) c("WT", x)),
               numCol = 1, 
               plotDF = luad.2018@meta.data, 
               outfileName = "s7hB")

# EXTENDED FIG 8B
v.s8b <- c("CXCL9", "CXCL10", "CXCL11")

createBoxplots(plotVector = v.s8b, w = 3, h = 4, df.pvals = df.pvals.tcga,
               xVar = "p53_status", 
               comparisonList = list(c("WT", "mut")),
               numCol = 3, 
               plotDF = luad.2018@meta.data, 
               outfileName = "s8b")

# EXTENDED FIG 8C
v.s8c <- c("TAM.CXCL")

createBoxplots(plotVector = v.s8c, w = 3, h = 4, df.pvals = df.pvals.tcga,
               xVar = "p53_status", 
               comparisonList = list(c("WT", "mut")),
               numCol = 1, 
               plotDF = luad.2018@meta.data, 
               outfileName = "s8c")

# EXTENDED FIG 8D,E
v.s8de <- c("SPP1", "CXCL9", "CXCL10", "CXCL11")

levels.tmp <- levels(luad.2018@meta.data$comut)
createBoxplots(plotVector = v.s8de, w = 3, h = 5,
               xVar = "comut", 
               comparisonList = lapply(setdiff(levels.tmp, "WT"), function(x) c("WT", x))[-6],
               numCol = 4, 
               plotDF = luad.2018@meta.data, 
               outfileName = "s8d")

levels.tmp <- levels(luad.2018@meta.data$DNE_LOFclass)
createBoxplots(plotVector = v.s8de, w = 3, h = 5,
               xVar = "DNE_LOFclass", 
               comparisonList = lapply(setdiff(levels.tmp, "WT"), function(x) c("WT", x)),
               numCol = 4, 
               plotDF = luad.2018@meta.data, 
               outfileName = "s8eA")

levels.tmp <- levels(luad.2018@meta.data$impact)
createBoxplots(plotVector = v.s8de, w = 3, h = 5,
               xVar = "impact", 
               comparisonList = lapply(setdiff(levels.tmp, "WT"), function(x) c("WT", x))[-3],
               numCol = 4, 
               plotDF = luad.2018@meta.data, 
               outfileName = "s8eB")

# EXTENDED FIG 9D
v.s9d <- c("T.Exhausted", "CD8.GZMK")

createBoxplots(plotVector = v.s9d, w = 3, h = 4, df.pvals = df.pvals.tcga,
               xVar = "p53_status", 
               comparisonList = list(c("WT", "mut")),
               numCol = 2, 
               plotDF = luad.2018@meta.data, 
               outfileName = "s9d")

# EXTENDED FIG 9F
v.s9f <- c("CD274", "CD86", "PVR", "PDCD1", "CTLA4", "TIGIT")

createBoxplots(plotVector = v.s9f, df.pvals = df.pvals.tcga,
               xVar = "p53_status", 
               comparisonList = list(c("WT", "mut")),
               numCol = 6, 
               plotDF = luad.2018@meta.data, 
               outfileName = "s9f")

# EXTENDED FIG 9G
v.s9g <- c("PVR", "CD274")
createBoxplots(plotVector = v.s9g, w = 3, h = 4, 
               xVar = "p53_status", 
               comparisonList = list(c("WT", "mut")),
               numCol = 2, 
               plotDF = cptac.luad.prot@meta.data, 
               outfileName = "s9g")
# ----------------------------------------------------
# smoking pack years
# ----------------------------------------------------
smoking_comparison_p53 <- read_csv(file.path(data.path, "smoking_comparison_p53.csv"))
p53.filt <- c("MGH1172", "MGH1173", "MGH1175", "MGH1182", "BWH16", "BWH11", "BWH14", "BWH06")
WT.filt <- c("P14", "MGH1174", "MGH1176", "MGH1183", "BWH01", "BWH04", "BWH05", "BWH09", "BWH19", "BWH23") # no MGH1170


xVar <- "p53_status"
plotVar <- "pack_yrs"
comparisonList = list(c("WT", "mut"))
plotDF <- smoking_comparison_p53 |>
#     mutate(p53_status = case_when(Patient %in% c(p53.filt, "P17", "P19", "MGH1181") ~ "mut", Patient %in% c(WT.filt, "MGH1179") ~ "WT", .default = NA)) |>
    mutate(p53_status = case_when(Patient %in% p53.filt ~ "mut", Patient %in% WT.filt ~ "WT", .default = NA)) |>
    rename("pack_yrs" = "Pack yrs")

p.smoking_pack_yrs <- plotDF |>
    drop_na(all_of(xVar)) |>
    ggboxplot(x = xVar, 
              y = plotVar,
              color = xVar, 
              palette = c("#0091CA", "#D8423D"),
              add = "jitter",
              title = plotVar) +
    stat_compare_means(comparisons = comparisonList) +
    scale_x_discrete(labels = function(x) paste0(x, "\n(n = ", table(plotDF[[xVar]])[x], ")")) +
    rotate_x_text(0) +
    xlab(NULL) +
    ylab("Pack years") +
    theme(plot.title = element_text(size = 12, face = "bold"), legend.position = "none")
p.smoking_pack_yrs

ggsave(filename = file.path(paste0(figures.dir, "smoking.pack_years.pdf")), 
       plot = p.smoking_pack_yrs, 
       device = "pdf",
       width = 3, 
       height = 4)

# plotDF |>
#     drop_na(all_of(xVar)) |>
#     select(Patient, pack_yrs, p53_status) |>
#     write.csv(file.path(data.path, "smoking.pack_years.csv"), row.names = F)
# # ----------------------------------------------------