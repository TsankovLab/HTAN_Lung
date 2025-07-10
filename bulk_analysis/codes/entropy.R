######################################################
#################### TSANKOV LAB #####################
######################################################
# ----------------------------------------------------
# load the libraries
# ----------------------------------------------------
library(tidyverse)
library(Seurat)
library(scuttle)
library(SCENT)
library(ggpubr)
# ----------------------------------------------------
# set directories
# ----------------------------------------------------
data.path <- "../data/"
figures.dir <- "../../figures/"
system(paste0("mkdir -p ", figures.dir))
# ----------------------------------------------------
# load the objects
# ----------------------------------------------------
luad.2018 <- readRDS(file.path(data.path, "luad.2018.processed.rds"))
combined.luad.2018 <- readRDS(file.path(data.path, "combined.luad.2018.rds"))
cptac.luad.prot <- readRDS(file.path(data.path, "cptac.luad.prot.processed.rds"))
# ----------------------------------------------------
print(combn(levels(combined.luad.2018@meta.data$p53_status), m = 2, simplify = F))
# ----------------------------------------------------
# custom function
# ----------------------------------------------------
createEntropyViolinplots <- function(plotVar, 
                                     plotDF,
                                     titleStr = "Bulk entropy TCGA",
                                     ylabStr = "Entropy",
                                     xVar = "p53_status", 
                                     outfileName = "noFigName",
                                     comparisonList = list(c("WT", "mut"))) {
    # Check if plotVar exists in plotDF columns
    if (!str_remove_all(string = plotVar, pattern = "`") %in% colnames(plotDF)) {
        return (message(paste0(plotVar, "doesn't exist in plotting data... Try again!")))
    }
    message(paste0("Adding plot and saving csv file for ", plotVar, "..."))
    saveDF <- plotDF |>
        drop_na(all_of(xVar)) |>
        select(all_of(c(xVar, plotVar)))
    # write.csv(saveDF, file.path(paste0(data.path, outfileName, ".", make.names(plotVar), ".csv")))
    
    # Determine the color palette based on xVar levels
    if (length(unique(plotDF[[xVar]])) < 4) {
        if (length(unique(plotDF[[xVar]])) == 2) {
            tmpPalette <- c("#00AFBB", "#FC4E07")
        } else if (length(unique(plotDF[[xVar]])) == 3) {
            tmpPalette <- c("grey", "#00AFBB", "#FC4E07")
        }

        # Generate the violin plot
        plot <- plotDF |>
            drop_na(all_of(xVar)) |>
            ggviolin(
                x = xVar, 
                y = plotVar, 
                fill = xVar, 
                title = titleStr,
                palette = tmpPalette,
                scale = "width",
                add = "boxplot", 
                add.params = list(fill = "white")
            ) +
            rotate_x_text(angle = 0) +
            stat_compare_means(comparisons = comparisonList) +
            scale_x_discrete(labels = function(x) {
                paste0(x, "\n(n = ", table(plotDF[[xVar]])[x], ")")
            }) +
            theme(legend.position = "none") +
            ylab(ylabStr) +
            xlab(NULL)
    } else {
        # Generate the violin plot without specifying a color palette
        plot <- plotDF |>
            drop_na(all_of(xVar)) |>
            ggviolin(
                x = xVar, 
                y = plotVar, 
                fill = xVar, 
                title = titleStr,
                scale = "width",
                add = "boxplot", 
                add.params = list(fill = "white")
            ) +
            rotate_x_text(angle = 90) +
            stat_compare_means(comparisons = comparisonList) +
            scale_x_discrete(labels = function(x) {
                paste0(x, "\n(n = ", table(plotDF[[xVar]])[x], ")")
            }) +
            theme(legend.position = "none") +
            ylab(ylabStr) +
            xlab(NULL)
    }

    return (plot)
}
# ----------------------------------------------------
# generate figures
# ----------------------------------------------------
# FIG 3H
p.3h.tcga <- createEntropyViolinplots(plotVar = "ccat", outfileName = "3h.tcga",
                         plotDF = combined.luad.2018@meta.data, 
                         xVar = "p53_status", 
                         comparisonList = combn(levels(combined.luad.2018@meta.data[["p53_status"]]), m = 2, simplify = FALSE))
p.3h.cptac <- createEntropyViolinplots(plotVar = "ccat.v", outfileName = "3h.cptac",
                         plotDF = cptac.luad.prot@meta.data, 
                         titleStr = "Bulk entropy CPTAC", 
                         xVar = "p53_status", 
                         comparisonList = list(c("WT", "mut")))
p.3h.tcga + p.3h.cptac
ggsave(filename = file.path(paste0(figures.dir, "3h.pdf")), 
       plot = (p.3h.tcga + p.3h.cptac), 
       device = "pdf", 
       width = 3 * 2, 
       height = 4)


# # EXTENDED FIG 4C
# levels.tmp <- levels(luad.2018@meta.data$comut)
# p.s4c <- createEntropyViolinplots(plotVar = "ccat", outfileName = "s4c",
#                          titleStr = NULL,
#                          plotDF = luad.2018@meta.data, 
#                          xVar = "comut", 
#                          comparisonList = lapply(setdiff(levels.tmp, "WT"), function(x) c("WT", x)))
# ggsave(filename = file.path(paste0(figures.dir, "s4c.pdf")), 
#        plot = p.s4c, 
#        device = "pdf", 
#        width = 3.5, 
#        height = 5)


# FIG 3K
levels.tmp <- levels(combined.luad.2018@meta.data$DNE_LOFclass)
p.3kA <- createEntropyViolinplots(plotVar = "ccat", outfileName = "3kA",
                         titleStr = NULL,
                         plotDF = combined.luad.2018@meta.data, 
                         xVar = "DNE_LOFclass", 
                         comparisonList = lapply(setdiff(levels.tmp, "WT"), function(x) c("WT", x)))
levels.tmp <- levels(combined.luad.2018@meta.data$impact)
p.3kB <- createEntropyViolinplots(plotVar = "ccat", outfileName = "3kB",
                         titleStr = NULL,
                         plotDF = combined.luad.2018@meta.data, 
                         xVar = "impact", 
                         comparisonList = lapply(setdiff(levels.tmp, "WT"), function(x) c("WT", x)))
p.3kA + p.3kB + patchwork::plot_layout(axes = "collect_y")
ggsave(filename = file.path(paste0(figures.dir, "3k.pdf")), 
       plot = (p.3kA + p.3kB + patchwork::plot_layout(axes = "collect_y")), 
       device = "pdf", 
       width = 4 * 2, 
       height = 5)

# EXTENDED FIG 5C
# # ----------------------------------------------------
# # new entropy plot with comut4 column
# # ----------------------------------------------------
# df <- luad.2018@meta.data
# df.krasSubsets <- df |> 
#     filter(comut3 != "unassigned") |>
#     mutate(comut4 = comut3) |>
#     select(ccat, comut4)
# df <- df |>
#     mutate(comut4 = comut2) |>
#     select(ccat, comut4)

# df <- bind_rows(df, df.krasSubsets)
# levels.comut4 <- c("WT", "EGFR", "KRAS", "KRAS-G12C", "KRAS-G12D", "KRAS-G12V", "TP53", "TP53_EGFR", "TP53_KRAS")
# df$comut4 <- factor(df$comut4, levels = levels.comut4)

# levels.tmp <- levels(df$comut4)
# p.new_entropy.comut4 <- createEntropyViolinplots(plotVar = "ccat", outfileName = "new_entropy.comut4",
#                                                  titleStr = NULL,
#                                                  plotDF = df, 
#                                                  xVar = "comut4", 
#                                                  comparisonList = lapply(setdiff(levels.tmp, "WT"), function(x) c("WT", x)))
                                                                         
# ggsave(filename = file.path(paste0(figures.dir, "new_entropy.comut4.pdf")), 
#        plot = p.new_entropy.comut4, 
#        device = "pdf", 
#        width = 5, 
#        height = 5)   
# # ----------------------------------------------------
# # FIG 3J- previous plot but with adjacent normals
# # ----------------------------------------------------
df <- luad.2018@meta.data
df.krasSubsets <- df |> 
    filter(comut3 != "unassigned") |>
    mutate(comut4 = comut3) |>
    select(ccat, comut4)
df <- df |>
    mutate(comut4 = comut2) |>
    select(ccat, comut4)
df <- bind_rows(df, df.krasSubsets)

df.adj_normal <- combined.luad.2018@meta.data |>
    filter(pancancer_type == "LUAD_adj_normal") |>
    mutate(comut4 = "Normal") |>
    select(ccat, comut4)
df <- bind_rows(df, df.adj_normal)

levels.comut4 <- c("Normal", "WT", "EGFR", "KRAS", "KRAS-G12C", "KRAS-G12D", "KRAS-G12V", "TP53", "TP53_EGFR", "TP53_KRAS")
df$comut4 <- factor(df$comut4, levels = levels.comut4)

levels.tmp <- levels(df$comut4)
p.new_entropy.comut4.adjacentNormal <- createEntropyViolinplots(plotVar = "ccat", outfileName = "new_entropy.comut4.adjacentNormal",
                                                                titleStr = NULL,
                                                                plotDF = df, 
                                                                xVar = "comut4", 
                                                                comparisonList = lapply(setdiff(levels.tmp, "WT"), function(x) c("WT", x)))
                                                                         
ggsave(filename = file.path(paste0(figures.dir, "3j.pdf")), 
       plot = p.new_entropy.comut4.adjacentNormal, 
       device = "pdf", 
       width = 6, 
       height = 5)                                                                                                                                               
# ----------------------------------------------------
