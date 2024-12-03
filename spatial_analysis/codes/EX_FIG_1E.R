data.path <- "./data/"
figures.dir <- "../../figures/"
system(paste0("mkdir -p ", figures.dir))

library(ggplot2)
library(ggpubr)
library(gplots)
library(data.table)

load(paste0(data.path, "spatial_list.Rda"))

ncounts_list <- list()
for (sample in names(spatial_list)){
    obj <- spatial_list[[sample]]
    ncounts_list[[sample]] <- obj@meta.data$nCount_Spatial
}

ncounts_list <- lapply(ncounts_list, log10)


df <- do.call(rbind, lapply(names(ncounts_list), function(name) {
  data.frame(Group = name, Value = ncounts_list[[name]])
}))

desired_order <- c('1174.C', '1174.D', 'BW09.A', 'BW09.B', 'BW09.C', 'BW09.D', 'BW19.C', 'BW19.D', 'BW23.A', 'BW23.B', 'BW11.C', 'BW11.D', 'BW14.A', 'BW14.B', 'BW16.A', 'BW16.B', '1179.A', '1179.B', '1179.C', '1179.D')

df$Group <- factor(df$Group, levels = desired_order)

pdf(paste0(figures.dir, 'EX_FIG_1E.pdf'), useDingbats = F, width = 10, height = 10)
ggplot(df, aes(x = Group, y = Value)) +
  geom_boxplot(
    outlier.shape = NA,
  ) +
  theme_minimal() +
    theme(
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
      ) +
  labs(title = "Boxplot with Variable Group Sizes",
       x = "Group",
       y = "Values")
dev.off()