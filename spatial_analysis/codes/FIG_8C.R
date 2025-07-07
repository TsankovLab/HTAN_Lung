data.path <- "./data/"
figures.dir <- "../../figures/"
system(paste0("mkdir -p ", figures.dir))

library(ComplexHeatmap)

df <- read.csv(paste0(data.path, 'NMF_cell_type_values_15.csv'), row.names=1)

df <- t(apply(df, 1, function(x) x / sum(x)))
df <- t(as.data.frame(df))

rownames(df) <- paste0("NMF_", 0:14)

col_order <- c(
  "Lymphatic", "Pericyte", "VEC.COL4A1", "CAF.ADH1B", "VEC.Cycling", "VEC.IFN",
  "Pulmonary.Venous", "Arterial", "Capillary", "Aerocyte", "Mast cell", "Pericyte.EMT",
  "Sm.Muscle.Vascular", "T/NK.Cycling", "NKT", "Plasmacytoid.DC", "CD8.GZMK", "CD4.TRM",
  "CD8.TRM", "DC1", "Mreg.DC", "T.Stress", "T.Exhausted", "NK.CD56.bright", "CD16.Mono",
  "NK.CD56.dim", "CD14.Mono", "DC2", "TAM.FOLR2", "TAM.CXCLs", "DC3", "Mono-Mac",
  "TAM.APOE", "TAM.FABP4", "Myeloid.Cycling", "TAM.AZU1", "CAF.COLs", "TAM.SPP1",
  "Myofibroblast", "CAF.APOE", "CAF.Ribo", "CAF.Adventitial", "CAF.Complement",
  "CAF.ISGs", "CAF.Cycling", "Ciliated", "Glycolysis.hypoxia", "CC.G2/M", "Hypoxia",
  "CC.S", "pEMT", "Secretory", "OxPhos", "MHCII", "StressResponse", "TNFA.NFKB", "IFNG",
  "Metallothionein", "Senescence", "B.MZ", "Plasma.IGHA", "Plasma.IGHG", "Treg", "TFH",
  "CD8.IFN", "Nonmalignant epi", "B.Follicular", "B.Cycling", "Sm.Muscle.Airway",
  "Systemic.Venous", "AT2-like", "CD4.Naive.CM", "Respiration.MT", "AT1/2-like", "Ribosome"
)

df <- df[,col_order]

pdf(paste0(figures.dir, "FIG_8C.pdf"), width = 20, height = 5, useDingbats=FALSE)
Heatmap(
  df,
  name = "Participation in Factor",
  col = colorRampPalette(RColorBrewer::brewer.pal(9, "RdPu"))(100),
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  row_names_side = "left",   # Place row labels on the left
  column_names_rot = 90      # Rotate column labels to vertical
)
dev.off()
