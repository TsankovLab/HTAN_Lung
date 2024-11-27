data.path <- "./data/"
figures.dir <- "../../figures/"
system(paste0("mkdir -p ", figures.dir))

library(Seurat)
library(ggplot2)
library(ggpubr)
library(patchwork)

load(paste0(data.path, "spatial_list.Rda"))

spatial_list <- spatial_list[c(1:4, 7:16)]

obj <- spatial_list[[1]]

spectral_colors <- scale_fill_gradientn(colors = rev(RColorBrewer::brewer.pal(11, "Spectral")),
                                        limits = c(0, 1),  # Fix scale limits
                                        breaks = seq(0, 1, by = 0.2))

plots <- lapply(c("Pericyte", "NK", "Myeloid", "Bcell", "Mast", "Plasma"), function(feature) {
  plot <- SpatialFeaturePlot(obj, features = feature, 
                             min.cutoff = 0, 
                             max.cutoff = 1) +
    spectral_colors +  # Apply the spectral color palette
    theme(legend.key.width = unit(1, "cm"),  # Adjust legend width
          legend.key.height = unit(0.5, "cm")) # Adjust legend height
  return(plot)
})
combined_plot <- wrap_plots(plots)

pdf(paste0(figures.dir, 'EX_FIG_1F.pdf'), useDingbats = F, width = 10, height = 10)
print(combined_plot)
dev.off()

plots <- lapply(c("Pericyte", "NK", "Myeloid", "Bcell", "Mast", "Plasma"), function(feature) {
  plot <- SpatialFeaturePlot(obj, features = feature, alpha=c(0, 0), 
                             min.cutoff = 0, 
                             max.cutoff = 1) +
    spectral_colors +  # Apply the spectral color palette
    theme(legend.key.width = unit(1, "cm"),  # Adjust legend width
          legend.key.height = unit(0.5, "cm")) # Adjust legend height
  return(plot)
})
combined_plot <- wrap_plots(plots)
pdf(paste0(figures.dir, 'EX_FIG_1F_no_spots.pdf'), useDingbats = F, width = 10, height = 10)
print(combined_plot)
dev.off()

plots <- lapply(c("Pericyte", "NK", "Myeloid", "Bcell", "Mast", "Plasma"), function(feature) {
  plot <- SpatialFeaturePlot(obj, features = feature, image.alpha=0, 
                             min.cutoff = 0, 
                             max.cutoff = 1) +
    spectral_colors +  # Apply the spectral color palette
    theme(legend.key.width = unit(1, "cm"),  # Adjust legend width
          legend.key.height = unit(0.5, "cm")) # Adjust legend height
  return(plot)
})
combined_plot <- wrap_plots(plots)
pdf(paste0(figures.dir, 'EX_FIG_1F_no_hist.pdf'), useDingbats = F, width = 10, height = 10)
print(combined_plot)
dev.off()