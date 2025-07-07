data.path <- "./data/"
figures.dir <- "../../figures/"
system(paste0("mkdir -p ", figures.dir))

library(ggplot2)
library(ggpubr)
library(tidyr)

df <- read.csv(paste0(data.path, "Naomi_barplot_input2.csv"))

df[, -1] <- log10(df[, -1]+1)

df_long <- pivot_longer(df, cols = -Sample, names_to = "Variable", values_to = "Value")
df_long$Sample <- factor(df_long$Sample, levels = c("HTAPP_NSC001bLT_FT","HTAPP_NCS004bLT_FT", "HTAPP_NCS006bLT_FT", "HTAPP_NCS011bLT_FT_"))
df_long$Sample2 <- df_long$Sample
df_long$Sample <- "WT"
df_long$Sample[df_long$Sample2 %in% c("HTAPP_NCS006bLT_FT", "HTAPP_NCS011bLT_FT_")] <- "MUT"
df_long$Sample <- factor(df_long$Sample, levels = c("WT", "MUT"))
df_long$Variable <- factor(df_long$Variable, levels = c("CD8_PD1", "Total_CD8", "Tumor_PDL1", "Total_Tumor"))

# Create paired boxplots
pdf(paste0(figures.dir, 'FIG_7B.pdf'))
ggplot(df_long, aes(x = Variable, y = Value, col = Sample)) +
  geom_boxplot(position = position_dodge(width = 0.8), width = 0.6, alpha = 0.8) +
  scale_color_manual(values = c("MUT" = "#D8423D", "WT" = "#0091CA"), name = "Sample") +
  labs(title = "Paired Boxplots", x = "Column Names", y = "Quantity") +
  stat_compare_means(label.x = 0.9, label.y = .9 * max(as.numeric(df_long$Value)), method='wilcox.test', size = 2.5) +
  theme_bw()
dev.off()
