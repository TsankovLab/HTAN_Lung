data.path <- "./data/"
figures.dir <- "../../figures/"
system(paste0("mkdir -p ", figures.dir))

library(ggplot2)
library(ggpubr)

temp <- read.csv(paste0(data.path, "Jason_generate_plots.csv"))
temp <- temp[temp$measured %in% "True",]
temp[,c("frame_name", "region_label", "measured", "label", "minimum_cells", "TP53")] <- NULL
temp$name <- factor(temp$name, levels = c('Direct contact', 'Adjacent', 'Proximity (<50um)'))
temp$p53_status <- "TP53WT"
temp$p53_status[temp$sample_name %in% c("HTAPP_NCS006bLT_FT", "HTAPP_NCS011bLT_FT_")] <- "TP53mut"
temp$p53_status <- factor(temp$p53_status, levels = c("TP53WT", "TP53mut"))
temp <- temp[temp$phenotype_label %in% "CD8+PD1+ % of Total",]

p <- ggplot(temp, aes(x = p53_status, y = percent, color = p53_status)) + 
  geom_point(aes(size = denominator_coverage_percent), 
             shape = 21, 
             position = position_dodge(width = 0.3)) +
  scale_color_manual(values = c("TP53WT" = "#0091CA", "TP53mut" = "#D8423D")) +
  facet_grid(. ~ name) +
  theme_minimal() +
  labs(title = "Percent by p53 Status", x = "p53 Status", y = "Percent") +
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = "bottom") +
  stat_compare_means(method = "wilcox.test", label = "p.format")
pdf(paste0(figures.dir, 'FIG_6C.pdf'))
print(p)
dev.off()