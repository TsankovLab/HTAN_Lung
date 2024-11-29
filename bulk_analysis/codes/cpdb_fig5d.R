######################################################
#################### TSANKOV LAB #####################
######################################################
# ----------------------------------------------------
# load the libraries
# ----------------------------------------------------
library(Seurat)
library(tidyverse)
# ----------------------------------------------------
# set directories
# ----------------------------------------------------
data.path <- "./data/"
figures.dir <- "../../figures/"
system(paste0("mkdir -p ", figures.dir))
# ----------------------------------------------------
# bar plots
# ----------------------------------------------------
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