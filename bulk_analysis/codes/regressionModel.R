######################################################
#################### TSANKOV LAB #####################
######################################################
# ----------------------------------------------------
# load the libraries
# ----------------------------------------------------
library(Seurat)
library(tidyverse)
library(scales)
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
df <- luad.2018@meta.data
# ----------------------------------------------------
# custom functions
# ----------------------------------------------------
# custom function #1: a function to convert WT to 0 and mut to 1
convertMutStatus <- function(column) {
  column <- ifelse(column == "WT", 0, ifelse(column == "mut", 1, column))
  return (column)
}
# custom function #2: a function to generate linear modeling results
generateModelResults <- function(df, programs, formula) {
    final_model_results <- list()

    for (i in seq_along(programs)) {
        message(programs[[i]])
        # Fit the specified linear model formula
        result <- lm(as.formula(paste0("df[[programs[[i]]]] ~ ", formula)), data = df)
        numVar <- nrow(coef(summary(result)))
        coefs <- coef(summary(result))
        mutations <- rownames(coefs)[2:numVar]
        coef_vals <- coefs[2:numVar, 1]
        t_vals <- coefs[2:numVar, 3]
        pvals <- coefs[2:numVar, 4]

        mutations <- toupper(str_remove_all(string = mutations, pattern = "_status1|_status"))
        final_model_results[[i]] <- data.frame(row.names = paste0(mutations, "|", programs[[i]]),
                                               mutation = mutations,
                                               coef = coef_vals,
                                               t = t_vals,
                                               pval = pvals,
                                               program = programs[[i]])

        # Wilcoxon test for uncorrected TP53
        group_wt <- df[[programs[[i]]]][df$p53_status == 0]
        group_mut <- df[[programs[[i]]]][df$p53_status == 1]

        # Check if both groups have enough non-missing values
        if (length(na.omit(group_wt)) > 0 && length(na.omit(group_mut)) > 0) {
            wilcox_result <- wilcox.test(group_mut, group_wt)
            uncorrected_pval <- wilcox_result$p.value
            uncorrected_tval <- (median(group_mut, na.rm = TRUE) - median(group_wt, na.rm = TRUE))
            
            mutationStr <- "uncorrected_p53"
            uncorrectedp53_result <- data.frame(row.names = paste0(mutationStr, "|", programs[[i]]),
                                                mutation = mutationStr,
                                                coef = NA,
                                                t = uncorrected_tval,
                                                pval = uncorrected_pval,
                                                program = programs[[i]])

            final_model_results[[i]] <- rbind(final_model_results[[i]], uncorrectedp53_result)
        } else {
            message("Skipping Wilcoxon test for ", programs[[i]], " due to insufficient data.")
        }
    }

    final_model_results_df <- do.call(rbind, final_model_results)
    final_model_results_df$neglogpval <- -log10(final_model_results_df$pval)
    return(final_model_results_df)
}
# custom function #3: a function to generate a dot plot using linear modeling results as an input
generateRegressionPlot <- function(final_model_results_df, selected_programs, selected_mutations) {
    combined <- final_model_results_df
    combined <- combined[combined$program %in% selected_programs, ]
    combined <- combined[!combined$mutation %in% "CANCER_PROP", ]
    if (!"neglogpval" %in% colnames(combined)) {
        combined$neglogpval <- -log10(combined$pval)
    }
    if (!"signif" %in% colnames(combined)) {
        # Set significance levels
        combined$signif <- 0
        combined$signif[combined$neglogpval >= 1] <- 1
        combined$signif[combined$neglogpval >= 1.3] <- 2
        combined$signif <- factor(combined$signif, levels = c(0, 1, 2))
    }
    # Adjust factor levels for mutation and program columns
    combined$mutation <- factor(combined$mutation, levels = selected_mutations)
    combined$program <- factor(combined$program, rev(selected_programs))
    
    if (!"neglogpvaldir" %in% colnames(combined)) {
        # Calculate direction of coefficients for color scale
        combined$direction <- combined$t
        combined$direction[combined$direction > 0] <- 1
        combined$direction[combined$direction < 0] <- -1
        combined$neglogpvaldir <- combined$neglogpval * combined$direction
    }

    range <- max(abs(range(combined$neglogpval)))
    colors <- c("navyblue", "white", "firebrick")
    pal <- gradient_n_pal(colors)
    if ("t" %in% colnames(combined)) {
        custom_color_scale <- scale_fill_gradientn(colours = pal(c(0, rescale(seq_along(combined$t)), 1)),
                                                   values = c(0, rescale(seq_along(combined$t)), 1),
                                                   limits = c(-range, range))
    } else {
        custom_color_scale <- scale_fill_gradientn(colours = pal(c(0, rescale(seq_along(combined$neglogpvaldir)), 1)),
                                                   values = c(0, rescale(seq_along(combined$neglogpvaldir)), 1),
                                                   limits = c(-range, range))
    }

    # Generate the plot
    regressionPlot <- ggplot(data = combined, mapping = aes(x = mutation, y = program, color = neglogpvaldir, size = neglogpval)) +
        geom_point(shape = 21, aes(colour = as.factor(signif), fill = neglogpvaldir)) +
        scale_colour_manual(values = c("00FFFFFF", "gray", "black")) +
        custom_color_scale +
        theme_minimal() +
        theme(text = element_text(size = 12),
              strip.text = element_text(size = 22, face = 'bold'),
              panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5),
              axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
        xlab(NULL) + 
        ylab(NULL)

    return(regressionPlot)
}
# ----------------------------------------------------
# data preprocessing
# ----------------------------------------------------
# List of columns to transform
colsToTransform <- c("p53_status", "kras_status", "keap1_status", "rbm10_status", "ptprd_status", "egfr_status", "stk11_status")
# Convert specified columns to character type
df[colsToTransform] <- lapply(df[colsToTransform], as.character)
df[colsToTransform] <- lapply(df[colsToTransform], convertMutStatus)
# ----------------------------------------------------
# generate figures
# ----------------------------------------------------
# FIG 3B
formula <- "p53_status + egfr_status + kras_g12c_status + kras_g12d_status + kras_g12v_status + stk11_status + keap1_status + rbm10_status + ptprd_status + Cancer_prop"
programs <- c("p53_targets", "AT.SFTPA2,1", "CC.G2M", "Glycolysis/Hypox", "pEMT", "Endothelial", "Pericytes")

finalModelResultsDF <- generateModelResults(df, programs, formula)
# finalModelResultsDF |>
#     filter((program %in% selectedPrograms) & (mutation %in% selectedMutations)) |>
#     write.csv(file.path(data.path, "3b.csv"))
# Selected programs and mutations for the plot
selectedPrograms <- c("p53_targets", "AT.SFTPA2,1", "CC.G2M", "Glycolysis/Hypox", "pEMT", "Endothelial", "Pericytes")
selectedMutations <- c("uncorrected_p53", "EGFR", "KRAS_G12C", "KRAS_G12D", "KRAS_G12V", "STK11", "KEAP1", "RBM10", "PTPRD", "P53")

p.3b <- generateRegressionPlot(finalModelResultsDF, selectedPrograms, selectedMutations)
ggsave(filename = file.path(paste0(figures.dir, "3b.pdf")), 
       plot = p.3b, 
       device = "pdf", 
       width = 5, 
       height = 6.25)


# FIG 3C
# ---- add code here ---


# EXTENDED FIG 2F: regression plot added during final revision; new_regression_plot.csv file received from Will
df.s2f <- read_csv(file.path(data.path, "new_regression_plot.csv")) |>
    column_to_rownames("...1")

selectedPrograms <- sort(unique(df.s2f$program))
selectedMutations <- c("uncorrected", "Histology corrected", "Stage corrected", "Histology + stage corrected")

p.s2f <- generateRegressionPlot(df.s2f, selectedPrograms, selectedMutations)
ggsave(filename = file.path(paste0(figures.dir, "s2f.pdf")), 
       plot = p.s2f, 
       device = "pdf", 
       width = 5, 
       height = 6.25)


# # EXTENDED FIG 4B old- taken out during revision 
# formula <- "p53_status + kras_status + egfr_status + stk11_status + Cancer_prop"
# programs <- c("p53_targets", "AT.SFTPA2,1", "CC.G2M", "Glycolysis/Hypox", "pEMT", "Endothelial", "Pericytes")

# finalModelResultsDF <- generateModelResults(df, programs, formula)
# # finalModelResultsDF |>
# #     filter((program %in% selectedPrograms) & (mutation %in% selectedMutations)) |>
# #     write.csv(file.path(data.path, "s3b.csv"))
# # Selected programs and mutations for the plot
# selectedPrograms <- c("p53_targets", "AT.SFTPA2,1", "CC.G2M", "Glycolysis/Hypox", "pEMT", "Endothelial", "Pericytes")
# selectedMutations <- c("uncorrected_p53", "EGFR", "KRAS", "STK11", "P53")

# p.s3b <- generateRegressionPlot(finalModelResultsDF, selectedPrograms, selectedMutations)
# ggsave(filename = file.path(paste0(figures.dir, "s4b_old.pdf")), 
#        plot = p.s3b, 
#        device = "pdf", 
#        width = 4, 
#        height = 6.25)
# ----------------------------------------------------
# correct for p53, kras, egfr, stk11, and tumor purity 
# to generate p-values for boxplots
# ----------------------------------------------------
v.1f <- c("Endothelial", "Pericytes")
v.2e <- c("p53_targets", "AT2", "CC.G2M", "`Glycolysis/Hypox`", "pEMT")
v.4d <- c("TAM.SPP1")
v.s1d <- c("B.cells", "Fibroblast", "Cancer", "Mast", "`T.CD4+`", # row 1
           "Myeloid", "NK.cells", "Smooth.muscle", "T.cells", "`T.CD8+`") # row 2
v.s3acd <- c("p53_targets", "AT2", "CC.G2M", "`Glycolysis/Hypox`", "pEMT") # row 1
v.s3g <- c("p53_targets", "AT2", "CC.G2M", "Glycolysis.Hypox", "pEMT")
v.s3h <- c("p53_targets", "AT2", "CC.G2M")
v.s5a <- c("Endothelial")
v.s5e <- c("Aerocyte", "Arterial")
v.s6f <- c("CAF.ADH1B")
v.s6g <- c("Pericytes")
v.s7b <- c("CXCL9", "CXCL10", "CXCL11")
v.s7c <- c("TAM.CXCL")
v.s7de <- c("SPP1", "CXCL9", "CXCL10", "CXCL11")
v.s8d <- c("T.Exhausted", "CD8.GZMK")
v.s8f <- c("CD274", "CD86", "PVR", "PDCD1", "CTLA4", "TIGIT")

v.boxplots <- c(v.1f, v.2e, v.4d, v.s1d, v.s3acd, v.s3g, v.s3h, v.s5a, v.s5e, v.s6f, v.s6g, v.s7b, v.s7c, v.s7de, v.s8d, v.s8f)
# remove ` from strings
v.boxplots <- unique(str_remove_all(string = v.boxplots, pattern = "`"))
v.boxplots

formula <- "p53_status + kras_status + egfr_status + stk11_status + Cancer_prop"
programs <- v.boxplots

# Generate the results
finalModelResultsDF <- generateModelResults(df, programs, formula)
# write.csv(finalModelResultsDF, file.path(data.path, "tcga.regression_plots.p53_kras_egfr_stk11_CancerProp.csv"))

finalModelResultsDF.filt <- finalModelResultsDF |>
    filter(mutation == "P53") |>
    select(program, pval)
# write.csv(finalModelResultsDF.filt, file.path(data.path, "tcga.regression_plots.p53_kras_egfr_stk11_CancerProp.filtered.csv"))
# ----------------------------------------------------
