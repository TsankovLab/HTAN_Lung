data.path <- "./data/"
figures.dir <- "../../figures/"
system(paste0("mkdir -p ", figures.dir))

library(ggplot2)
library(ggpubr)

rois <- read.csv(paste0(data.path, "IP_plot_input_rois_Jason.csv"))

temp <- rois
temp <- temp[temp$ONCOTREE_PRIMARY_DIAGNOSIS %in% "LUAD",]
temp$p53_status <- temp$TP53_status
temp$p53_status <- factor(temp$p53_status, levels = c("TP53 wt", "TP53 mut"))
temp <- temp[temp$label %in% c("CYTOK PDL1+"),]
temp <- temp[temp$phenotype_label %in% "CD8+PD1+ % of Total",]
temp$percent <- temp$numerator/temp$denominator

# Create the plot
p <- ggplot(temp, aes(x = p53_status, y = percent, color = p53_status)) + 
  geom_boxplot(position = position_dodge(width = 0.3), width = 0.6) +
  scale_color_manual(values = c("TP53 wt" = "#0091CA", "TP53 mut" = "#D8423D")) +
  facet_grid(. ~ name) +
  theme_minimal() +
  labs(title = "Percent by p53 Status", x = "p53 Status", y = "Percent") +
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = "none") +
  stat_compare_means(method = "wilcox.test", label = "p.format", label.x = 0.9, label.y = .9 * max(as.numeric(temp[,"percent"])))

system(paste0("mkdir -p ", figures.dir))
# Show the plot
pdf(paste0(figures.dir, "FIG_6D.pdf"))
print(p)
dev.off()