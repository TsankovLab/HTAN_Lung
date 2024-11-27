data.path <- "./data/"
figures.dir <- "../../figures/"
system(paste0("mkdir -p ", figures.dir))

library(ggplot2)
library(ggpubr)

rois <- read.csv(paste0(data.path, "jason_input_spp1.csv"))
rois$p53_status[rois$p53_status %in% "WT"] <- "TP53 wt"
rois$p53_status[rois$p53_status %in% "mut"] <- "TP53 mut"
rois$p53_status <- factor(rois$p53_status, levels = c("TP53 wt", "TP53 mut"))

temp <- rois
temp <- temp[temp$phenotype_label %in% "% of Total cells that are CD31+ within 50um of CD68+OPN+",]

p <- ggplot(temp, aes(x = p53_status, y = percent, color = p53_status)) + 
  geom_boxplot(position = position_dodge(width = 0.3), width = 0.6, alpha = 0.8) +
  scale_color_manual(values = c("TP53 wt" = "#0091CA", "TP53 mut" = "#D8423D")) +
  theme_minimal() +
  labs(title = "% of Total cells that are CD31+ within 50um of CD68+OPN+", x = "p53 Status", y = "Percent") +
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = "none") +
  stat_compare_means(method = "wilcox.test", label = "p.format")
pdf(paste0(figures.dir, "FIG_7G.pdf"))
print(p)
dev.off()