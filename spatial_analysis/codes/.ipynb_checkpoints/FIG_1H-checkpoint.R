data.path <- "./data/"
figures.dir <- "../../figures/"
system(paste0("mkdir -p ", figures.dir))

library(ggplot2)
library(ggpubr)

load(paste0(data.path, "spatial_list.Rda"))

end.prop.list <- list()
peri.prop.list <- list()
for (i in 1:length(spatial_list)) {
  meta <- spatial_list[[i]]@meta.data[,14:24] # high lvl
  end <- meta$Endothelial
  peri <- meta$Pericyte
  end.prop.list[[i]] <- end
  peri.prop.list[[i]] <- peri
}
end.p53 <- unlist(end.prop.list[1:4])
end.wt <- unlist(end.prop.list[5:14])
peri.p53 <- unlist(peri.prop.list[1:4])
peri.wt <- unlist(peri.prop.list[5:14])
end.df <- data.frame(Endothelial = c(end.p53, end.wt), p53_status = c(rep("p53_mut", length(end.p53)), rep("WT", length(end.wt))) )
end.df$p53_status <- factor(end.df$p53_status, levels = c("WT", "p53_mut"))
peri.df <- data.frame(Pericyte = c(end.p53, end.wt), p53_status = c(rep("p53_mut", length(peri.p53)), rep("WT", length(peri.wt))) )
peri.df$p53_status <- factor(peri.df$p53_status, levels = c("WT", "p53_mut"))

df.t <- peri.df

pdf(paste0(figures.dir, "FIG_1H_peri.pdf"), width = 3, height = 3, useDingbats=FALSE)
ggboxplot(df.t, x = "p53_status", y = "Pericyte", color = "p53_status", outlier.shape = NA,
                       palette = c("#0091CA", "#D8423D"),
                       title = "Pericyte") +
stat_compare_means(label.x = 0.9, label.y = .9 * max(as.numeric(df.t[,"Pericyte"]))) +
theme(plot.title = element_text(size = 12, face = "bold")) + theme(legend.position = "none")
dev.off()

df.t <- end.df

pdf(paste0(figures.dir, "FIG_1H_endo.pdf"), width = 3, height = 3, useDingbats=FALSE)
ggboxplot(df.t, x = "p53_status", y = "Endothelial", color = "p53_status", outlier.shape = NA,
                       palette = c("#0091CA", "#D8423D"),
                       title = "Endothelial") +
stat_compare_means(label.x = 0.9, label.y = .9 * max(as.numeric(df.t[,"Endothelial"]))) +
theme(plot.title = element_text(size = 12, face = "bold")) + theme(legend.position = "none")
dev.off()