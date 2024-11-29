######################################################
#################### TSANKOV LAB #####################
######################################################
# ----------------------------------------------------
# load the libraries
# ----------------------------------------------------
library(tidyverse)
library(Seurat)
library(survival)
library(survminer)
library(ggpubr)
# ----------------------------------------------------
# set directories
# ----------------------------------------------------
data.path <- "./data/"
figures.dir <- "../../figures/"
system(paste0("mkdir -p ", figures.dir))
# ----------------------------------------------------
# KM-curves for EXTENDED FIG 2d, 2e, 7g ------------------
# ----------------------------------------------------
df.km <- read_csv(file.path(data.path, "df.km.csv")) |> column_to_rownames("...1")
# ----------------------------------------------------
# generate figures
# ----------------------------------------------------
# EXTENDED FIG 2d ------------------------------------
# ----------------------------------------------------
# EXTENDED FIG 2dA
fit <- survfit(Surv(Overall.Survival..Months., Overall.Survival.Status) ~ pEMT_top25_vs_mid50_vs_bottom25, data = df.km)

p.s2dA <- ggsurvplot(fit,
                     risk.table.col = "strata",
                     linetype = "solid",
                     surv.median.line = "h",
                     pval = T,
                     pval.method = T,
                     ggtheme = theme_classic(),
                     censor = F)
p.s2dA <- p.s2dA$plot + 
    labs(color = NULL) +
    scale_color_discrete(labels = c("Low", "Medium", "High")) +
    xlab("Time (months)")

# EXTENDED FIG 2dB
fit <- survfit(Surv(Overall.Survival..Months., Overall.Survival.Status) ~ Glycolysis.Hypox_top25_vs_mid50_vs_bottom25, data = df.km)

p.s2dB <- ggsurvplot(fit,
                     risk.table.col = "strata",
                     linetype = "solid",
                     surv.median.line = "h",
                     pval = T,
                     pval.method = T,
                     ggtheme = theme_classic(),
                     censor = F)
p.s2dB <- p.s2dB$plot + 
    labs(color = NULL) +
    scale_color_discrete(labels = c("Low", "Medium", "High")) +
    xlab("Time (months)")

# EXTENDED FIG 2dC
fit <- survfit(Surv(Overall.Survival..Months., Overall.Survival.Status) ~ CC.G2M_top25_vs_mid50_vs_bottom25, data = df.km)

p.s2dC <- ggsurvplot(fit,
                     risk.table.col = "strata",
                     linetype = "solid",
                     surv.median.line = "h",
                     pval = T,
                     pval.method = T,
                     ggtheme = theme_classic(),
                     censor = F)
p.s2dC <- p.s2dC$plot + 
    labs(color = NULL) +
    scale_color_discrete(labels = c("Low", "Medium", "High")) +
    xlab("Time (months)")

# patchwork for EXTENDED FIG 2d
p.s2d <- p.s2dA + p.s2dB + p.s2dC + 
    patchwork::plot_layout(axis_titles = "collect_y", guides = "collect")

ggsave(filename = file.path(figures.dir, "s2d.pdf"), 
       plot = p.s2d, device = "pdf", width = 3.5 * 3, height = 3)


# EXTENDED FIG 7g
fit <- survfit(Surv(Overall.Survival..Months., Overall.Survival.Status) ~ SPP1_top50_vs_bottom50, data = df.km)
p.s7g <- ggsurvplot(fit,
           # conf.int = TRUE,
           risk.table.col = "strata", # Change risk table color by groups
           linetype = "solid", # Change line type by groups
           surv.median.line = "h", # Specify median survival
           pval = T,
           pval.method = T,
           ggtheme = theme_classic(), # Change ggplot2 theme
           # palette = c("#E7B800", "#2E9FDF")
           censor = F)

p.s7g <- p.s7g$plot + 
    labs(color = NULL) +
    scale_color_manual(labels = c("Low SPP1", "High SPP1"),
                       values = c("#0091CA", "#D8423D")) +
    xlab("Time (months)")

ggsave(filename = file.path(figures.dir, "s7g.pdf"), 
       plot = p.s7g, device = "pdf", width = 3, height = 3.5)
# ----------------------------------------------------
# EXTENDED FIG 2e
# ----------------------------------------------------
fit <- survfit(Surv(Overall.Survival..Months., Overall.Survival.Status) ~ p53_status2, data = df.km)
print(fit)

p.s2e <- ggsurvplot(fit, 
           # conf.int = TRUE,
           risk.table.col = "strata", # Change risk table color by groups
           pval = T,
           pval.method = T,
           ggtheme = theme_classic(), # Change ggplot2 theme
           )
p.s2e <- p.s2e$plot + 
    labs(color = NULL) +
    scale_color_manual(labels = c("TP53 WT", "TP53 mut"),
                       values = c("#0091CA", "#D8423D")) +
    xlab("Time (months)")

ggsave(filename = file.path(figures.dir, "s2e.pdf"), 
       plot = p.s2e, device = "pdf", width = 3, height = 3.5)

# ----------------------------------------------------
# FIG 2d ---------------------------------------------
# ----------------------------------------------------
df <- read_csv(file.path(data.path, "survival_maligv9_pvals_tp53_status_top50_vs_bottom50.csv")) |>
    column_to_rownames("...1")
df$info <- rownames(df)
df <- df |> 
    mutate(outcome = case_when(posneg == -1 ~ "Worse", posneg == 1 ~ "Improved", .default = NA), 
           neglog10pval = -log10(pvals) * posneg) |>
    mutate(classification = case_when(str_detect(string = info, pattern = "top50_vs_bottom50") ~ "top50_vs_bottom50", 
                                      str_detect(string = info, pattern = "top25_vs_bottom25") ~ "top25_vs_bottom25", 
                                      str_detect(string = info, pattern = "top25_vs_mid50_vs_bottom25") ~ "top25_vs_mid50_vs_bottom25", 
                                      .default = NA)) |> 
    mutate(name = str_remove_all(string = info, pattern = paste0("_", classification)))  |> 
    filter(classification == "top50_vs_bottom50") |>
    select(all_of(c("pvals", "neglog10pval", "name", "outcome")))|> 
    filter(!name == "SPP1") |> 
    mutate(name = case_when(name == "HLA_D.AT2" ~ "MHCII", 
                            name == "AT.SFTPA2.1" ~ "AT2-like", 
                            name == "AT.SFTPC" ~ "AT1/2-like", 
                            name == "S100A.KRT.Trachea.ECM.KRASup" ~ "Senescence", 
                            .default = name))

p.2d <- ggdotchart(df, x = "name", y = "neglog10pval",
           color = "outcome",                                # Color by groups
           palette = c("#0091CA", "#D8423D"), # Custom color palette
           sorting = "descending",                       # Sort value in descending order
           add = "segments",                             # Add segments from y = 0 to dots
           add.params = list(color = "lightgray", size = 2), # Change segment color and size
           group = "outcome",                                # Order by groups
           dot.size = 8,                                 # Large dot size
           label = round(df$neglog10pval,1),                        # Add mpg values as dot labels
           font.label = list(color = "white", size = 8, 
                             vjust = 0.5),               # Adjust label parameters
           ggtheme = theme_pubr()                        # ggplot2 theme
  ) +
  geom_hline(yintercept = 0, linetype = 1, color = "black") +
  geom_hline(yintercept = -1.3, linetype = 2, color = "lightgray") +
  geom_hline(yintercept = 1.3, linetype = 2, color = "lightgray") +
    xlab(NULL)

ggsave(filename = file.path(figures.dir, "2d.pdf"), 
       plot = p.2d, device = "pdf", width = 7, height = 7)
# ----------------------------------------------------