library(slam)  # NB.var.genes run apply function in sparse matrices for finding row variance
library(MASS)  # NB.var.genes uses fitdistr

WriteDataToPortal <- function(countData, outCluster, outDir='./sc_portal_prep', outStub='outPortal',
                              padCluster=T,writeDataTable=T) {
  #' This function takes a Seurat object and clustering information and prepares files that
  #' can be uploaded to the single cell portal at https://portals.broadinstitute.org/single_cell
  #' This code works with Seurat v2.0 and was modified from Caroline Porter seurat.utils_cbmp.R.
  #' 
  #' @param countData Seurat object with gene-cell expression matrix and tSNE data
  #' @param outCluster Vector of cluster assignments for each cell in countData. The cluster labels 
  #' will likely come from countData@ident, and can be either numeric or named subsets. 
  #' These labels are placed in both the tSNE and metadata text files. 
  #' @return Three files are generated. A zipped file containing the gene-cell 
  #' expression matrix. A txt file containing number of genes and number of UMIs/cell.
  #' A txt file containing tSNE representation and cluster assignments.
  
  require(data.table)
  require(R.utils)
  
  options(datatable.showProgress=T)
  
  ct = outCluster
  
  # Write expression matrix
  if (writeDataTable) {
    out.ExpTableName = sprintf('%s/%s_expData.txt',outDir,outStub)
    dataMat = data.frame(Matrix::as.matrix(countData@data))
    data.table::fwrite(as.list(c("GENE",colnames(countData@data))),file=out.ExpTableName,showProgress = F,col.names = F)
    data.table::fwrite(dataMat,file=out.ExpTableName,row.names = T,col.names = F,append = T,showProgress = T)
    gzip(out.ExpTableName)
  }
  
  # Write meta data (nGene, nUMI, cluster)
  out.metaTableName = sprintf('%s/%s_metaData.txt',outDir,outStub)
  #metaData = data.table(countData@meta.data[,1:2],keep.rownames = T)
  #metaData = data.table(countData@meta.data[,c("nGene", "nUMI")], keep.rownames = T)
  metaData = data.table(countData@meta.data[,c("nGene", "nUMI", "orig.ident","orig.identSec","Phase")], keep.rownames = T)
  #metaData = data.table(countData@meta.data[,c("nGene", "orig.ident", "experiment", "case", "therapy", "clonal", "annotate")], keep.rownames = T)
  
  metaData[,"typeID"] = countData@ident
  
  if (padCluster) {
    metaData[,"SeuratClustering"] = sprintf("%02d",ct)
  } else {
    metaData[,"SeuratClustering"] = as.character(ct)
  }
  data.table::fwrite(as.list(c("NAME",colnames(metaData[,-1]))),file=out.metaTableName)
  zTypes = as.list(gsub("character","group",gsub("integer|double","numeric",sapply(metaData[1,-1],typeof))));
  data.table::fwrite(as.list(c("TYPE",zTypes)),file=out.metaTableName,append = T)
  data.table::fwrite(metaData,file=out.metaTableName,append = T)
  
  # Write tSNE data (tSNE_1, tSNE_2, cluster)
  out.clustTableName = sprintf('%s/%s_tSNE.txt',outDir,outStub)
  metaData = data.table(row.names = countData@cell.names,X=countData@dr$tsne@cell.embeddings[, "tSNE_1"],
                        Y=countData@dr$tsne@cell.embeddings[, "tSNE_2"])
  
  if (padCluster) {
    metaData[,"SeuratClustering"] = sprintf("%02d",ct)
  } else {
    metaData[,"SeuratClustering"] = as.character(ct)
  }
  data.table::fwrite(as.list(c("NAME",colnames(metaData[,-1]))),file=out.clustTableName)
  zTypes = as.list(gsub("character","group",gsub("integer|double","numeric",sapply(metaData[1,-1],typeof))));
  data.table::fwrite(as.list(c("TYPE",zTypes)),file=out.clustTableName ,append = T)
  data.table::fwrite(metaData,file=out.clustTableName,append = T)
}

# test[[1]] <- NB.var.genes(test[[1]], x.high.cutoff = 100, x.low.cutoff = 0.01, do.text = TRUE, num.sd = 1.3)
NB.var.genes <- function(object, cells.use=NULL, genes.use = NULL, do.plot=TRUE,set.var.genes=TRUE,
                              x.low.cutoff=0.005, x.high.cutoff=3, diffCV.cutoff=NULL,num.sd=NULL, 
                              cex.use=0.5,cex.text.use=0.5,do.spike=FALSE,pch.use=16, col.use="black", 
                              spike.col.use="red",do.ident=FALSE, do.text=TRUE, reads.use=FALSE, 
                              cut.quantile=1, sort.results=TRUE) {
  #' This function was written by Karthik Shekhar. I slightly modified this to work with Seurat2 objects
  #' I replaced object@count.data with object@raw.data (raw UMI) and I left object@data alone (log 
  #' trasformed UMI). I replaced object@mean.var with object@hvg.info. I named the dataframe columns in
  #' mv.df mean and dispersion.
  #' 
  print("Identifying variable genes based on UMI Counts")
  cells.use=set.ifnull(cells.use,colnames(object@data))
  genes.use=set.ifnull(genes.use, rownames(object@data))
  if (!reads.use){
    count.data=object@raw.data[genes.use,cells.use]  # modified by Orr
  } else{
    count.data = object@reads.data[genes.use,cells.use]  # this option will not work with Seurat2
  }
  
  # Empirical mean, var and CV (modified by Orr Ashenberg to deal with sparse matrices)
  # mean_emp = apply(count.data, 1, mean)  # creates a dense matrix
  # var_emp = apply(count.data, 1, var)  # creates a dense matrix
  mean_emp = Matrix::rowMeans(count.data)
  var_emp = rowapply_simple_triplet_matrix(as.simple_triplet_matrix(count.data), FUN = var)
  
  genes.use=names(mean_emp)[mean_emp > 0]
  mean_emp = mean_emp[genes.use]
  var_emp = var_emp[genes.use]
  cv_emp = sqrt(var_emp) / mean_emp
  # NB sampling
  a=Matrix::colSums(count.data)
  a = a[a <= quantile(a,cut.quantile)]
  size_factor =  a/ mean(a)
  fit=fitdistr(size_factor, "Gamma")
  if (!reads.use){
    hist(size_factor, 50, probability=TRUE, xlab="N_UMI/<N_UMI>")
  } else {
    hist(size_factor, 50, probability=TRUE, xlab="N_Reads/<N_Reads>")
  }
  curve(dgamma(x, shape=fit$estimate[1], rate=fit$estimate[2]),from=0, to=quantile(size_factor, 0.95), add=TRUE, col="red",
        main="Gamma dist fit for size factor")
  text(5,0.6, paste("shape = ", round(fit$estimate[1],2)))
  text(5,0.5, paste("rate = ", round(fit$estimate[2],2)))
  
  # Gamma distributions of individual genes are just scaled versions. If X ~ Gamma(a,b)
  # then cX ~ Gamma(a, b/c)
  a_i = rep(fit$estimate[1], length(mean_emp)); names(a_i) = names(mean_emp)
  b_i = fit$estimate[2] / mean_emp; names(b_i) = names(mean_emp)
  mean_NB = a_i / b_i; var_NB = a_i*(1+b_i) / (b_i^2)
  cv_NB = sqrt(var_NB)/mean_NB
  diffCV = log(cv_emp) - log(cv_NB)
  
  hist(diffCV,500, main="Select a delta-logCV cutoff for variable gene: ", xlab="delta-logCV")
  
  if (!is.null(num.sd)){
    diffCV.cutoff = mean(diffCV) + num.sd*sd(diffCV)
  }
  if (is.null(diffCV.cutoff)){
    diffCV.cutoff = readline("Select a delta-logCV cutoff (genes with a higher value will be considered):")
    diffCV.cutoff = as.numeric(diffCV.cutoff)
  }
  
  
  print(paste0("Using diffCV = ", diffCV.cutoff, " as the cutoff"))
  abline(v=diffCV.cutoff)
  Sys.sleep(4)
  
  print(paste0("Considering only genes with mean counts less than ", x.high.cutoff, " and more than ", x.low.cutoff))
  pass.cutoff=names(diffCV)[which(diffCV > diffCV.cutoff & (mean_emp > x.low.cutoff & mean_emp < x.high.cutoff))]
  print(paste0("Found ", length(pass.cutoff), " variable genes"))
  mv.df=data.frame(gene.mean=mean_emp,gene.dispersion=cv_emp,gene.dispersion.scaled=0)  # modified by Orr
  rownames(mv.df)=names(mean_emp)
  object@hvg.info=mv.df  # modified by Orr

  if (do.spike) spike.genes=grep("^ERCC", rownames(count.data), value=TRUE)
  if (do.plot) {
    
    plot(mean_emp,cv_emp,pch=pch.use,cex=cex.use,col="black",xlab="Mean Counts",ylab="CV (counts)", log="xy")
    curve(sqrt(1/x), add=TRUE, col="red", log="xy", lty=2, lwd=2)
    or = order(mean_NB)
    lines(mean_NB[or], cv_NB[or], col="magenta", lwd=2)
    points(mean_emp[pass.cutoff], cv_emp[pass.cutoff], col="blue", pch=16, cex=cex.use)
    
    if (do.spike) points(mean_emp[spike.genes],cv_emp[spike.genes],pch=16,cex=cex.use,col=spike.col.use)
    if(do.text) text(mean_emp[pass.cutoff],cv_emp[pass.cutoff],pass.cutoff,cex=cex.text.use)
    
    if (do.ident) {
      identify(mean_emp,cv_emp,labels = names(mean_emp))
    }
  }
  if (set.var.genes) { 
    object@var.genes=pass.cutoff
    if (sort.results) {  # Orr added Seurat code
      object@hvg.info <- object@hvg.info[order(
        object@hvg.info$gene.dispersion,
        decreasing = TRUE
      ),]
    }
    return(object)
  }
  if (!set.var.genes) return(pass.cutoff)
}

set.ifnull=function(x,y) {
  if(is.null(x)) x=y
  return(x)
}

QCSS2 <- function(genome.mapping.file, transcriptome.mapping.file, qc.out.file, genic.out.file, sampleid, gcdata = NULL, 
                  exp.title = NULL)  {
  #' This function takes the genome_mapping_summary and transcriptome_mapping_summary files produced
  #' when running RNASeq-QC on SmartSeq2 data in the KCO RNAseq pipeline. It summarizes some of these 
  #' quality metrics related to exon content and number of genes per cell.
  #' If a Seurat object is provided, the function updates a few of the quality metrics. 
  #' The Seurat object contains the TPM gene expression matrix that has been filtered 
  #' for low-quality genes and cells. The function also plots the fraction of reads mapping
  #' to genic, exonic, intronic, and intergenic regions for each experimental replicate.
  #' 
  #' @param genome.mapping.file This file is generated by
  #' /seq/regev_genome_portal/SOFTWARE/KCO/RNASEQ_pipeline/generate_sample_summary_stats.pl.
  #' It contains information on how well reads aligned to the genome.
  #' @param transcriptome.mapping.file This file is generated by
  #' /seq/regev_genome_portal/SOFTWARE/KCO/RNASEQ_pipeline/util/summarize_rnaseqQC_results.pl.
  #' It contains information on how well genome-mapped reads align to known transcript features.
  #' @param qc.out.file File name for tsv file containing quality control metrics that have been 
  #' collected across cells.
  #' @param genic.out.file File name for plot of fraction of reads mapping to different genomic regions.
  #' @param sampleid Vector of sample names that should match the sample names in @param genome.mapping.file 
  #' and @param transcriptome.mapping.file. In addition if using @param gcdata, these names should match those
  #' in the gcdata@meta.data slot sample column. 
  #' @param gcdata A Seurat objet. The @meta.data slot should have a column named sample, which has 
  #' the sample name for each cell.
  #' @param exp.title Plot title for @param genic.out.file 
  #' @return A list with two elements: a dataframe with the quality control metrics, named table, and a
  #' ggplot object of the fraction of reads mapping to different genomic regions, named plot.genic .
  
  # Read genome_mapping_summary, taking specific columns
  qc.genome.vals <- read.table(genome.mapping.file, sep = "\t", stringsAsFactors = F)
  colnames(qc.genome.vals) <- strsplit(readLines(genome.mapping.file, n = 1), "\t")[[1]]
  qc.genome.vals <- qc.genome.vals[, c(1, 2, 5:8, 13)]
  colnames(qc.genome.vals)[1] <- "sample"  # used for merging qc.genome.vals and qc.transcriptome.vals
  colnames(qc.genome.vals)[7] <- "Number of Genes"
  
  # Read transcriptome_mapping_summary, taking specific columns
  qc.transcriptome.vals <- read.table(transcriptome.mapping.file, skip = 2, sep = "\t", stringsAsFactors = F)
  qc.colnames <- strsplit(readLines(transcriptome.mapping.file, n = 2)[2], "\t")[[1]][c(7:8, 10:11, 13, 28:31)]
  qc.transcriptome.vals <- qc.transcriptome.vals[, c(7:8, 10:11, 13, 28:31)]
  colnames(qc.transcriptome.vals) <- qc.colnames
  colnames(qc.transcriptome.vals)[1] <- "sample"  # used for merging qc.genome.vals and qc.transcriptome.vals
  
  # Merge genome and transcriptome mapping summaries based on the sample names, which are shared.
  qc <- merge(qc.genome.vals, qc.transcriptome.vals, by = "sample")
  qc <- qc[match(sampleid, qc$sample), ]  # sort rows using order of samples in sampleid
  
  # If Seurat object available, update a few quality metrics after the filtering.
  if (!is.null(gcdata)) {  
    qc$`Number of Genes` <- NULL  # no longer need the estimate from genome_mapping_summary
    qc <- merge(gcdata@meta.data[, c("sample", "nGene")], qc, by = "sample")
    colnames(qc)[2] <- "Number of Genes"
    qc <- qc[match(gcdata@meta.data$sample, qc$sample), ]  # sort rows using order of samples in Seurat object
    qc$sample <- as.character(qc$sample)  # do not use factor as this leads to warnings when adding the median values later
  }
  
  # Plot fraction of reads mapping to different genomic regions.
  # Store fraction of reads mapping to different genomic regions in dataframe.
  qc.colnames <- c("genic", "exonic", "intronic", "intergenic")
  qc.genic <- 100 * qc[ ,c("Intragenic Rate", "Exonic Rate", "Intronic Rate", "Intergenic Rate")]
  colnames(qc.genic) <- qc.colnames
  
  # Gather QC dataframe into genomic region : percent mapping key-value pairs for making violin plots. 
  qc.violin <- gather(qc.genic, "region", "percent", c(1:4))
  # Set order in which genic regions and samples are plotted in barplot
  qc.violin$region <- factor(qc.violin$region, levels = qc.colnames[1:4]) 
  p <- ggplot(data = qc.violin, mapping = aes(x = region, y = percent)) + geom_violin(fill = "grey", color = "grey30") + 
    geom_jitter(height = 0, size = 1) +
    scale_y_continuous(breaks=seq(0,100,10), limits = c(0, 100)) + 
    ylab("percent bases mapping") + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + ggtitle(exp.title)
  ggsave(genic.out.file, width = 6, height = 4)
  
  # Add median value of each column, other than the first column which contains the names for each sample.
  median.vals <- apply(qc[-c(1)], 2, median)
  qc <- rbind(c("median", median.vals), qc)
  
  # Write qc table to file.
  write.table(t(qc), qc.out.file, sep = "\t", quote = F)
  
  # Return qc table and genic plot objects.
  list(table = qc, plot.genic = p)
}

QCDropSeq <- function(exon.files, gene.files, qc.out.file, genic.out.file, gcdata = NULL, 
                      exp.labels = NULL, exp.title = NULL)  {
  #' This function takes two files produced by the DropSeq pipeline and 
  #' extracts quality metrics related to exon content and number of genes per cell barcode.
  #' If a Seurat object is provided, the function updates a few of the quality metrics. 
  #' The Seurat object contains the digital gene expression matrix that has been filtered 
  #' for low-quality genes and cells. The function also plots the fraction of reads mapping
  #' to genic, exonic, intronic, and intergenic regions for each experimental replicate.
  #' 
  #' @param exon.files Vector of file names with percent of bases encoding introns and exons. Names 
  #' of vector elements should be the sample names.
  #' @param gene.files Vector of file names with total number of genes and transcript counts (total UMI) 
  #' per cell barcode. Names of vector elements should be the sample names.
  #' @param qc.out.file File name for tsv file containing quality control metrics that have been 
  #' concatenated across several samples.
  #' @param genic.out.file File name for plot of fraction of reads mapping to different genomic regions.
  #' @param gcdata A list of Seurat objects, one for each sample. This should follow the same order of
  #' samples as in @param exon.file and @param gene.file. The names of the gcdata list elements 
  #' should be the sample names.
  #' @param exp.labels Alternative labels to use in @param genic.out.file for the sample names. These
  #' should have the same length as the sample names in @param exon.files
  #' @param exp.title Plot title for @param genic.out.file 
  #' @return A list with two elements: a dataframe with the quality control metrics, named table, and a
  #' ggplot object of the fraction of reads mapping to different genomic regions, named plot.genic .
  qc <- data.frame()
  sampleid <- names(exon.files)  # names of experimental replicates
  for (sample in sampleid) {
    # Read exon information
    lines <- readLines(exon.files[sample]) 
    qc.exon.vals <- strsplit(lines[8], "\t")[[1]][11:16]
    qc.exon.vals <- round(as.numeric(qc.exon.vals), digits = 2)
    names(qc.exon.vals) <- strsplit(lines[7], "\t")[[1]][11:16]
    
    # Fraction of bases passing filter that aligned = PF_ALIGNED_BASES/PF_BASES
    qc.aligned <- strsplit(lines[8], "\t")[[1]][1:2] # PF_BASES PF_ALIGNED_BASES 
    qc.aligned <- as.numeric(qc.aligned)
    qc.aligned <- round(qc.aligned[2] / qc.aligned[1], digits = 2)
    names(qc.aligned) <- "PCT_GENIC_BASES"
    
    # Read total gene and UMI per barcode summary
    qc.gene.vals <- read.table(gene.files[sample], skip = 2, header = T)
    qc.gene.vals <- round(c(nrow(qc.gene.vals), median(qc.gene.vals$NUM_GENES), 
                            median(qc.gene.vals$NUM_TRANSCRIPTS)))
    names(qc.gene.vals) <- c("Estimated Number of Cells", "Median Genes per Cell", 
                             "Median UMI Counts per Cell")
    
    qc <- rbind(qc, t(data.frame(c(qc.gene.vals, qc.aligned, qc.exon.vals))))  # add qc for this sample to dataframe
  }
  rownames(qc) <- sampleid

  # If Seurat object available, update a few quality metrics after the filtering.
  if (!is.null(gcdata)) {  
    stopifnot(all.equal(names(gcdata), names(exon.files)))  # check the samples correspond to one another
    qc$`Estimated Number of Cells` <- sapply(gcdata, function(data) length(data@ident))
    qc$`Median Genes per Cell` <- sapply(gcdata, function(data) as.integer(median(data@meta.data$nGene)))
    qc$`Median UMI Counts per Cell` <- sapply(gcdata, function(data) as.integer(median(data@meta.data$nUMI)))
    qc$`Total Genes Detected` <- sapply(gcdata, function(data) nrow(data@data))
    qc <- qc[, c(1:3, 11, 4:10)]
  }
  
  # Write qc table to file
  write.table(t(qc), qc.out.file, sep = "\t", quote = F)
  
  # Plot fraction of reads mapping to different genomic regions.
  # Store fraction of reads mapping to different genomic regions in dataframe.
  qc.colnames <- c("genic", "exonic", "intronic", "intergenic", "ribosomal", "replicate")
  qc.genic <- 100 * qc[,c("PCT_GENIC_BASES", "PCT_CODING_BASES", "PCT_INTRONIC_BASES", 
                          "PCT_INTERGENIC_BASES", "PCT_RIBOSOMAL_BASES")]
  qc.genic$replicate <- sampleid
  colnames(qc.genic) <- qc.colnames

  # Gather QC dataframe into genomic region : percent mapping key-value pairs for making bar plots. 
  qc.bar <- gather(qc.genic, "region", "percent", c(1:5))
  # Set order in which genic regions and samples are plotted in barplot
  qc.bar$region <- factor(qc.bar$region, levels = qc.colnames[1:5])  
  qc.bar$replicate <- factor(qc.bar$replicate, levels = sampleid)
  if (is.null(exp.labels)) {exp.labels <- sampleid}
  p <- ggplot(data = qc.bar) + geom_bar(stat = "identity", mapping = aes(x = region, y = percent, 
                                                                         fill = replicate),
                                        position = position_dodge(), color = "black") + 
    scale_y_continuous(breaks=seq(0,100,10), limits = c(0, 100)) + 
    ylab("percent bases mapping") + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
    scale_fill_discrete(name="Experimental\nCondition", breaks=sampleid, labels=exp.labels) + 
    ggtitle(exp.title)
  ggsave(genic.out.file, width = 6, height = 4)
   
  list(table = qc, plot.genic = p)
}

QC10x <- function(qc.files, qc.out.file, genic.out.file, gcdata = NULL, 
                  exp.labels = NULL, exp.title = NULL) {
  #' This function takes the metrics_summary.csv files made from several runs of Cell Ranger count 
  #' and concantenates them into a single file. If a Seurat object is provided, the function 
  #' updates a few of the quality metrics. The Seurat object contains the digital gene expression matrix 
  #' that has been filtered for low-quality genes and cells. The function also plots the fraction of 
  #' reads mapping to genic, exonic, intronic, and intergenic regions for each experimental replicate.
  #' 
  #' Each metrics_summary.csv file contains information on each sequenced molecule.
  #' [molecule_info.h5](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/output/molecule_info) 
  #' has information on each molecule from the 10x run. Each molecule is a unique (cell barcode, UMI, gene) 
  #' tuple and there is other associated information like number of reads and quality.
  #' 
  #' @param qc.files Vector of metrics_summary.csv file names, one for each sample, created by 
  #' Cell Ranger count. The names of the vector elements should be the sample names.
  #' @param qc.out.file File name for tsv file containing quality control metrics that have been 
  #' concatenated across several samples.
  #' @param genic.out.file File name for plot of fraction of reads mapping to different genomic regions.
  #' @param gcdata A list of Seurat objects, one for each sample. This should follow the same order of
  #' samples as in @param qc.files. The names of the gcdata list elements should be the sample names.
  #' @param exp.labels Alternative labels to use in @param genic.out.file for the sample names. These
  #' should have the same length as the sample names in @param exon.files
  #' @param exp.title Plot title for @param genic.out.file 
  #' @return A list with two elements: a dataframe with the quality control metrics, named table, and a
  #' ggplot object of the fraction of reads mapping to different genomic regions, named plot.genic .
  
  # Iterate over quality control metrics files and collect the metrics into a single dataframe
  qc <- data.frame()
  sampleid <- names(qc.files)   # names of experimental replicates
  for (sample in sampleid) {
    qc.add <- read.table(qc.files[sample], header = T, sep = ",", check.names = F)
    
    # Calculate percent of reads mapping to genic regions
    qc.genic <- unlist(qc.add[,c("Reads Mapped Confidently to Exonic Regions", 
                                      "Reads Mapped Confidently to Intronic Regions", 
                                      "Reads Mapped Confidently to Intergenic Regions",
                                      "Reads Mapped Antisense to Gene")])
    qc.genic <- sum(as.numeric(sub("%", "", qc.genic)))
    qc.genic <- paste0(qc.genic, "%")
    names(qc.genic) <- "Reads Mapped Confidently to Genic Regions"

    qc <- rbind(qc, cbind(qc.add[,c(1:5)], t(qc.genic), qc.add[,c(6:18)]))
  }
  rownames(qc) <- sampleid

  # If Seurat object available, update a few quality metrics after the filtering.
  if (!is.null(gcdata)) {  
    stopifnot(all.equal(names(gcdata), names(qc.files)))  # check the samples correspond to one another
    qc$`Estimated Number of Cells` <- sapply(gcdata, function(data) length(data@ident))
    qc$`Median Genes per Cell` <- sapply(gcdata, function(data) as.integer(median(data@meta.data$nGene)))
    qc$`Median UMI Counts per Cell` <- sapply(gcdata, function(data) as.integer(median(data@meta.data$nUMI)))
    qc$`Total Genes Detected` <- sapply(gcdata, function(data) nrow(data@data))
  }
  
  # Write to file
  write.table(t(qc), qc.out.file, sep = "\t", quote = F)
  
  # Plot fraction of reads mapping to different genomic regions.
  # Store fraction of reads mapping to different genomic regions in dataframe.
  qc.colnames <- c("genic", "exonic", "intronic", "intergenic", "replicate")
  # Remove percent symbols from fraction mapping reads.
  # If there is only one sample, the result of apply is a vector, which then needs to be transposed.
  if (length(sampleid) > 1) {
    qc.genic <- data.frame(apply(qc[,c("Reads Mapped Confidently to Genic Regions",
                                         "Reads Mapped Confidently to Exonic Regions", 
                                         "Reads Mapped Confidently to Intronic Regions", 
                                         "Reads Mapped Confidently to Intergenic Regions")], 
                               2, function(x) as.numeric(sub("%", "", x)))) 
  } else {
  qc.genic <- data.frame(t(apply(qc[,c("Reads Mapped Confidently to Genic Regions",
                                       "Reads Mapped Confidently to Exonic Regions", 
                                       "Reads Mapped Confidently to Intronic Regions", 
                                       "Reads Mapped Confidently to Intergenic Regions")], 
                                 2, function(x) as.numeric(sub("%", "", x))))) 
  }
  qc.genic$replicate <- sampleid
  colnames(qc.genic) <- qc.colnames

  # Gather QC dataframe into genomic region : percent mapping key-value pairs for making bar plots. 
  qc.bar <- gather(qc.genic, "region", "percent", c(1:4))
  # Set order in which genic regions and samples are plotted in barplot
  qc.bar$region <- factor(qc.bar$region, levels = qc.colnames[1:4])  
  qc.bar$replicate <- factor(qc.bar$replicate, levels = sampleid)
  if (is.null(exp.labels)) {exp.labels <- sampleid}
  p <- ggplot(data = qc.bar) + geom_bar(stat = "identity", mapping = aes(x = region, y = percent, 
                                                                         fill = replicate),
                                        position = position_dodge(), color = "black") + 
    scale_y_continuous(breaks=seq(0,100,10), limits = c(0, 100)) + 
    ylab("percent bases mapping") + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
    scale_fill_discrete(name="Experimental\nCondition", breaks=sampleid, labels=exp.labels) + 
    ggtitle(exp.title)
  ggsave(genic.out.file, width = 6, height = 4)
  
  list(table = qc, plot.genic = p)
}

MergeMultipleSeuratObjects <- function(gcdata, project) {
  #' This function uses the MergeSeurat function to merge two or more Seurat objects.
  #' Important note: do.normalize = F flag in Seurat::MergeSeurat() means that no normalization
  #' or scaling information (gcdata@data and gcdata@scale.data) exists in the merged Seurat object.
  #' gcdata@scale.data = NULL and gcdata@data = gcdata@raw.data. Normalization and scaling must be
  #' performed on the merged object before PCA.
  #' 
  #' @param gcdata List of Seurat objects where the object elements are named by their sample.
  #' @param project Project name to give merged Seurat object.
  #' @return Merged Seurat object.
  if (length(gcdata) == 1) {  # if list of Seurat objects only has 1 element, no merge needs to be done
    return(gcdata)
  }
  
  sampleid <- names(gcdata)  # sample names
  
  # Iteratively merge Seurat objects two at a time.
  gcdata.merge <- MergeSeurat(object1 = gcdata[[1]], object2 = gcdata[[2]], 
                              add.cell.id1 = sampleid[1], add.cell.id2 = sampleid[2], 
                              project = project, do.normalize = F)
  for (i in seq(gcdata)[-c(1, length(gcdata))]) {  # skip first and last index
    gcdata.merge <- MergeSeurat(object1 = gcdata.merge, object2 = gcdata[[i+1]], 
                                add.cell.id2 = sampleid[i+1], 
                                project = project, do.normalize = F)  
  }
  
  # Seurat has issues naming the identity class when there are underscores in the class names. 
  # This is an issue in the @ident slot, where it will only keep part of the name. This can be 
  # modified using names.field and names.delim but is annoying. For some reason, 
  # @meta.data$orig.ident does not have this issue, so we use it to name the @ident class for each cell.
  gcdata.merge <- SetIdent(gcdata.merge, ident.use = gcdata.merge@meta.data$orig.ident)
  gcdata.merge
}

MultipleFeaturePlot <- function(gcdata, features="nGene", feature.type="meta", ncols=ceiling(sqrt(length(features))), pt.size=1, same.scale=FALSE) {
  #' This function takes a feature from the @meta.data of a Seurat object, or from the un-scaled gene expression 
  #' @data of a Seurat object, and plots it on a tSNE with a blue to yellow to red gradient of colors. 
  #' If more than one feature is supplied, the graphs are sub-plotted. Metadata features and genes cannot be plotted
  #' simultaneously. This extends the Seurat::FeaturePlot() function
  #' 
  #' @param gcdata Seurat object with @data and @meta.data slots.
  #' @param features Either a vector of features to be plotted from @meta.data, or a vector of gene names to be
  #' plotted from @data.
  #' @param feature.type Either "gene" or "meta" to indicate the type of feature, @data or @meta.data respectively.
  #' @param ncols Number of columns desired for facet wrap when making subplots.
  #' @param pt.size Size of points in plot.
  #' @param same.scale TRUE if color bar should be scaled based on entire dataset rather than just 
  #' the expression for plotted features. This is only used when @param feature.type = "gene".
  #' @return ggplot object.
  #' 
  #' Written by Caroline Porter and color functions provided by Sam Riesenfield in original function plot.tsne.feature. 
  #' I modified the code to better deal with plotting more than one @meta.data feature. 
  #' 
  #' TO DO: Add code to check that the gene or feature actually exists in the data set, throw error if not.
  #' Add option to scale colar bar baed on entire dataset, rather than the plotted sub-set
  
  # Load required libraries 
  library(tidyr)
  
  # Get color gradient (THIS COULD BE UPDATED TO ALLOW ADDITIONAL INPUTS PER SAM'S SCRIPT)
  # source("/ahg/regevdata/users/cporter/code/colrs.fromSamRiesen.R")
  # colors <- get.hmap.col()
  colors<-c("#191970","#121285","#0C0C9A","#0707B0","#0101C5","#0014CF","#0033D3","#0053D8","#0072DD","#0092E1","#00B2E6",
            "#00D1EB","#23E8CD","#7AF17B","#D2FA29","#FFEB00","#FFC300","#FF9B00","#FF8400","#FF7800","#FF6B00","#FF5F00","#FF5300",
            "#FF4700","#F73B00","#EF2E00","#E62300","#DD1700","#D50B00","#CD0000")
  
  # Error out if user does not specify correct feature type 
  if (feature.type!="meta" & feature.type!="gene") {
    stop("feature type must be 'meta' or 'gene'")
  }
  
  # If plotting gene expression, verify that the desired genes to plot exist in the @data slot
  if (feature.type == "gene") {
    features <- features[which(features %in% rownames(gcdata@data))]
  }
  if (!length(features)) {
    stop("None of the genes requested for plotting are in gcdata@data.")
  }
  
  # Collect feature info either from @meta.data or @data
  if (feature.type=="meta") { 
    feature.info <- as.matrix(gcdata@meta.data[,features])   # column format
    if (length(features) == 1) {
      colnames(feature.info) <- features
    }        
  } else if (feature.type=="gene") { 
    feature.info <- as.matrix(gcdata@data[features,])  # row format
    if (length(features) > 1) {
      feature.info <- t(feature.info)  # transpose into more convenient column format.
    }        
    colnames(feature.info) <- features  # only required when length of features is 1.
  }
  
  # Build data frame of feature info and tSNE coordinates 
  tmp.df <- data.frame(feature.info, gcdata@dr$tsne@cell.embeddings)
  plot.df <- gather(tmp.df, name, val, 1:length(features), factor_key=TRUE)
  
  # Set scales for color bar 
  if (same.scale==TRUE & feature.type=="gene") {
    lower=min(gcdata@data) 
    upper=max(gcdata@data)
  } else if (same.scale==TRUE & feature.type=="meta") {
    stop("same.scale can only be set to true when plotting gene features")
  } else {
    lower=min(plot.df$val) 
    upper=max(plot.df$val)
  }
  
  # Color tSNE plot by gene expression or by meta data
  p <- ggplot(plot.df, aes(x=tSNE_1, y=tSNE_2)) + geom_point(aes(color=val), alpha=0.8, shape=16, size=pt.size) + 
    theme(aspect.ratio = 1) + scale_color_gradientn(colors=colors, limits=c(lower, upper))
  p <- p + theme(aspect.ratio=1, text = element_text(size=10), axis.text=element_text(size=6), 
                 strip.text.x = element_text(margin = ggplot2::margin(.1, 0, .1, 0, "cm")),
                 strip.text = element_text(size=12)) + 
    facet_wrap( ~ name, ncol=ncols)
  return(p)
}

RankModuleScore <- function(gcdata, name.cluster, gs, plotdir, ctrl.size=10) {
  #' The purpose of this function is to help somewhat automate assigning cell types to
  #' cell clusters. For each cell cluster, the cluster is scored by annotated sets of genes
  #' using Seurat::AddModuleScore(). Each cell is individually scored, and then the average 
  #' score for all cells within a cluster is recorded. The gene sets that score the highest 
  #' for a given cluster can potentially reveal what cell type that cluster is. This function 
  #' requires a Seurat object containing gene expression values, cells that have been clustered, and a 
  #' tSNE dimensional reduction to display the gene set scores mapped onto the cells. 
  #' The function also requires a list of gene sets to score each cell cluster with. After
  #' running this function, it is good to visually inspect the created tSNE plots to see how specific
  #' the best scoring gene sets are for each cell cluster of interest. 
  #' 
  #' This function uses dplyr functions, so make sure the plyr package is detached before running.
  #' 
  #' @param gcdata Seurat object with @data and @meta.data slots. Clustering must have been performed
  #' and the cluster identities must be stored as a column of gcdata@meta.data. In addition,
  #' there must be a tSNE dimensional reduction slot for making plots.
  #' @param name.cluster Name of column in gcdata@meta.data with cluster labels. This defines the
  #' cell clusters that are scored by each gene set.
  #' @param gs This is a GSA (R pacckage) object containing the names of each gene set and a list
  #' of genes that define that corresponding gene set. The object can be created using  
  #' gs <- GSA.read.gmt(file.gmt). gmt files are further described 
  #' [here](https://software.broadinstitute.org/cancer/software/gsea/wiki/index.php/Data_formats).
  #' @param plotdir For each gene set, a tSNE plot is colored by the scores from that gene set and placed
  #' in this plotting directory.
  #' @param ctrl.size When calculating the average expression level of each gene set, Seurat::AddModuleScore()
  #' subtracts off expression from a control set of genes. This parameters defines how many control genes to
  #' take from each expression bin.
  #' @return The function returns the dataframe cluster.scores. The first column is the cluster names from
  #' @param name.cluster, and the subsequent columns are the mean gene set scores for all cells within 
  #' that cluster.
  
  # Create directory to store tSNE plots colored by mean gene set expression.
  dir.create(plotdir, showWarnings = FALSE) 
  
  # Create data frame where first column is cluster number and
  # subsequent columns are mean scores for a given gene set for all cells in given cluster.
  cluster.scores <- data.frame("clusters" = unique(gcdata@meta.data[,name.cluster]))
  names(cluster.scores)[1] <- name.cluster

  # Iterate over each gene set and score each cell.
  for (i in seq(gs$geneset.names)) {
    name <- gs$geneset.names[i]
    geneset <- gs$genesets[[i]]
    print(c(name, geneset))
    
    # Score cells by gene sets.
    gcdata <- Seurat::AddModuleScore(gcdata, genes.list = list(geneset), ctrl.size = ctrl.size, 
                                     enrich.name = name)
    oldname <- paste0(name, "1")  # Remove 1 from end of name
    names(gcdata@meta.data)[names(gcdata@meta.data) == paste0(name, "1")] <- name
    
    # Rank clusters by average gene set score for cells in the cluster.
    # http://dplyr.tidyverse.org/articles/programming.html
    summ_name <- paste0("mean_", name)
    add.scores <- gcdata@meta.data %>% group_by(!!as.name(name.cluster)) %>% 
      summarise(!!summ_name := mean(!!as.name(name)))
    cluster.scores <- merge(cluster.scores, add.scores, by=name.cluster)  # need common name for cluster column

    # tSNE plot colored by gene set scores.
    MultipleFeaturePlot(gcdata, features = name, feature.type = "meta", pt.size = 2)
    ggsave(paste0(plotdir, "tSNE_", name, ".png"), width = 8, height = 8, dpi = 200)
  }
  
  # Order clusters by their cluster number. 
  cluster.scores <- cluster.scores %>% dplyr::arrange(as.numeric(as.character(!!as.name(name.cluster))))
  return(cluster.scores)
}
