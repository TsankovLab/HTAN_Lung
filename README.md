# HTAN_Lung

This repository contains custom code to reproduce main analyses in our manuscript, "A cellular and spatial atlas of TP53-associated tissue remodeling in lung adenocarcinoma."

Processed scRNA-seq data, spatial transcriptomic data, and mutation information can be downloaded at https://www.dropbox.com/scl/fo/1v6hek4v53r8njk6sqt80/h?rlkey=dp7qa0fcyx9njcxa895y6dh9z&dl=0. Please refer to the reporting summary document for the passcode to download the files.

System Requirements:
- Operating System: Tested on RedHat 7.7 x86_64.
- R Version: 4.1.0
- R Packages:
  - Seurat (v4.1.0) (https://satijalab.org/seurat/)
  - ComplexHeatmap (v2.10.0) (https://bioconductor.org/packages/release/bioc/html/ComplexHeatmap.html)
  - clusterProfiler (v4.2.2) (https://bioconductor.org/packages/release/bioc/html/clusterProfiler.html)
  - SCENT (v1.0.3) (https://github.com/aet21/SCENT)
  - ggpubr (v0.4.0) (https://rdrr.io/cran/ggpubr/)
  - ggplot2 (v3.3.6) (https://ggplot2.tidyverse.org)
  - cgdsr (v1.3.0) (https://rdrr.io/cran/cgdsr)
  - survival (v3.2.11) (https://rdrr.io/cran/survival/)
- Additional Hardware: No special hardware required.

Installation Guide:

```R
install.packages("Seurat", version = "4.1.0")
install.packages("ggpubr", version = "0.4.0")
install.packages("ggplot2", version = "3.3.6")
install.packages("cgdsr", version = "1.3.0")
install.packages("survival", version = "3.2.11")

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("ComplexHeatmap", version = "2.10.0")
BiocManager::install("clusterProfiler", version = "4.2.2")

if (!requireNamespace("devtools", quietly = TRUE))
    install.packages("devtools")
devtools::install_github("aet21/SCENT@v1.0.3")
```

For Demo and Instructions for use, please refer to the figures_main.R file.



