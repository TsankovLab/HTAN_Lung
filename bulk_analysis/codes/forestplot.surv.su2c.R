#####################@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ 
#####################@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
##################### FIGURE 6E AND EXTENDED 10B 
##################### FORESTPLOT OF SURVIVAL ANALYSES OF SU2C-MARK COHORT
library(forestplot)

# my outputs
data.path <- "./data/"
figures.dir <- "../../figures/"

if (!file.exists(paste0(figures.dir,'Fig6/'))){dir.create(paste0(figures.dir,'Fig6/'), recursive=TRUE)}


##################### 6E (left)
#####################
res.surv <- get(load(paste0(data.path, "surv.results.mut.adeno.pd1.su2c.Rda")))

##
coxph.pvalPFS = res.surv$coxph.pvalPFS
lbPFS = res.surv$lbPFS
hrPFS = res.surv$hrPFS
ubPFS = res.surv$ubPFS
sub = res.surv$sub
num_high =  res.surv$num_high
num_low = res.surv$num_low

coxph.pvalPFS.ann = NULL
for (i in 1:length(coxph.pvalPFS))
{
  tmp = coxph.pvalPFS[i]
  coxph.pvalPFS.ann[i] = as.character(cut(tmp,
                                          breaks = c(-Inf, 0.001, 0.1, Inf),
                                          labels = c(paste("p =", formatC(tmp, 1, format = "e")), 
                                                     paste("p =", formatC(signif(tmp, 2), 2, flag = "#")), 
                                                     "p > 0.1")))
}

sub2 = c("EGFR", "STK11", "KEAP1", "KRAS", "TP53", "  Impactful 1", "  Impactful 2", "  DNE,LOF", "  notDNE,LOF", "  notDNE,notLOF",
         "  Missense", "  Nonsense", "  Splice Site", "  Unique Splice Site", "  Frameshift")
sub2 = paste0(sub2, " (", num_high, ")")

res1 = data.frame(lbPFS, hrPFS, ubPFS)
res2 = data.frame(sub2, coxph.pvalPFS.ann)


##############
pdf(paste0(figures.dir, "Fig6/Fig6E_left.pdf"), useDingbats = F, 9, 9)
forestplot(res2, as.matrix(res1), graph.pos = 2, new_page = F,
           xlab = "Hazard ratio",
           xticks = c(0, 1, 2, 3, 4),
           graphwidth = unit(90, "mm"),
           vertices = TRUE,
           col = fpColors(lines = "black"),
           boxsize = 0.3,
           txt_gp = fpTxtGp(xlab = gpar(cex = 2.25), 
                            ticks = gpar(cex = 1.75),
                            title = gpar(cex = 2.5, fontface = "bold"),
                            cex = 1.5),
           align = c("l", "l"),
           zero = 1)
dev.off()



##################### Exteneded 10B (left)
#####################
res.surv <- get(load(paste0(data.path, "surv.results.mut.adeno.all.su2c.Rda")))

##
coxph.pvalPFS = res.surv$coxph.pvalPFS
lbPFS = res.surv$lbPFS
hrPFS = res.surv$hrPFS
ubPFS = res.surv$ubPFS
sub = res.surv$sub
num_high =  res.surv$num_high
num_low = res.surv$num_low

coxph.pvalPFS.ann = NULL
for (i in 1:length(coxph.pvalPFS))
{
  tmp = coxph.pvalPFS[i]
  coxph.pvalPFS.ann[i] = as.character(cut(tmp,
                                          breaks = c(-Inf, 0.001, 0.1, Inf),
                                          labels = c(paste("p =", formatC(tmp, 1, format = "e")), 
                                                     paste("p =", formatC(signif(tmp, 2), 2, flag = "#")), 
                                                     "p > 0.1")))
}

sub2 = c("EGFR", "STK11", "KEAP1", "KRAS", "TP53", "  Impactful 1", "  Impactful 2", "  DNE,LOF", "  notDNE,LOF", "  notDNE,notLOF",
         "  Missense", "  Nonsense", "  Splice Site", "  Unique Splice Site", "  Frameshift")
sub2 = paste0(sub2, " (", num_high, ")")

res1 = data.frame(lbPFS, hrPFS, ubPFS)
res2 = data.frame(sub2, coxph.pvalPFS.ann)


##############
pdf(paste0(figures.dir, "Fig6/Ext_Fig10B_left.pdf"), useDingbats = F, 9, 9)
forestplot(res2, as.matrix(res1), graph.pos = 2, new_page = F,
           xlab = "Hazard ratio",
           xticks = c(0, 1, 2, 3, 4),
           graphwidth = unit(90, "mm"),
           vertices = TRUE,
           col = fpColors(lines = "black"),
           boxsize = 0.3,
           txt_gp = fpTxtGp(xlab = gpar(cex = 2.25), 
                            ticks = gpar(cex = 1.75),
                            title = gpar(cex = 2.5, fontface = "bold"),
                            cex = 1.5),
           align = c("l", "l"),
           zero = 1)
dev.off()



##################### 6E (right)
#####################
res.surv <- get(load(paste0(data.path, "surv.results.mut.sq.pd1.su2c.Rda")))

##
coxph.pvalPFS = res.surv$coxph.pvalPFS
lbPFS = res.surv$lbPFS
hrPFS = res.surv$hrPFS
ubPFS = res.surv$ubPFS
sub = res.surv$sub
num_high =  res.surv$num_high
num_low = res.surv$num_low

coxph.pvalPFS.ann = NULL
for (i in 1:length(coxph.pvalPFS))
{
  tmp = coxph.pvalPFS[i]
  coxph.pvalPFS.ann[i] = as.character(cut(tmp,
                                          breaks = c(-Inf, 0.001, 0.1, Inf),
                                          labels = c(paste("p =", formatC(tmp, 1, format = "e")), 
                                                     paste("p =", formatC(signif(tmp, 2), 2, flag = "#")), 
                                                     "p > 0.1")))
}

sub2 = c("EGFR", "STK11", "KEAP1", "KRAS", "TP53", "  Impactful 1", "  Impactful 2", "  DNE,LOF", "  notDNE,LOF", "  notDNE,notLOF",
         "  Missense", "  Nonsense", "  Splice Site", "  Unique Splice Site", "  Frameshift")
sub2 = paste0(sub2, " (", num_high, ")")

res1 = data.frame(lbPFS, hrPFS, ubPFS)
res2 = data.frame(sub2, coxph.pvalPFS.ann)

##
res1 = res1[-1,]
res2 = res2[-1,]

## set range hr >17 -> 17
res1.mod = res1
se = which(res1.mod[,"ubPFS"] > 17)
res1.mod[se,"ubPFS"] = 17


##############
pdf(paste0(figures.dir, "Fig6/Fig6E_right.pdf"), useDingbats = F, 9, 9)
forestplot(res2, as.matrix(res1.mod), graph.pos = 2, new_page = F,
           xlab = "Hazard ratio",
           xticks = c(0, 1, 5, 10, 15, 17),
           graphwidth = unit(90, "mm"),
           vertices = TRUE,
           col = fpColors(lines = "black"),
           boxsize = 0.3,
           txt_gp = fpTxtGp(xlab = gpar(cex = 2.25), 
                            ticks = gpar(cex = 1.75),
                            title = gpar(cex = 2.5, fontface = "bold"),
                            cex = 1.5),
           align = c("l", "l"),
           zero = 1)
dev.off()



##################### Exteneded 10B (right))
##################### note that egfr has only 1 mut sample, so the regression couldn't run
res.surv <- get(load(paste0(data.path, "surv.results.mut.sq.all.su2c.Rda")))

##
coxph.pvalPFS = res.surv$coxph.pvalPFS
lbPFS = res.surv$lbPFS
hrPFS = res.surv$hrPFS
ubPFS = res.surv$ubPFS
sub = res.surv$sub
num_high =  res.surv$num_high
num_low = res.surv$num_low

coxph.pvalPFS.ann = NULL
for (i in 1:length(coxph.pvalPFS))
{
  tmp = coxph.pvalPFS[i]
  coxph.pvalPFS.ann[i] = as.character(cut(tmp,
                                          breaks = c(-Inf, 0.001, 0.1, Inf),
                                          labels = c(paste("p =", formatC(tmp, 1, format = "e")), 
                                                     paste("p =", formatC(signif(tmp, 2), 2, flag = "#")), 
                                                     "p > 0.1")))
}

sub2 = c("EGFR", "STK11", "KEAP1", "KRAS", "TP53", "  Impactful 1", "  Impactful 2", "  DNE,LOF", "  notDNE,LOF", "  notDNE,notLOF",
         "  Missense", "  Nonsense", "  Splice Site", "  Unique Splice Site", "  Frameshift")
sub2 = paste0(sub2, " (", num_high, ")")

res1 = data.frame(lbPFS, hrPFS, ubPFS)
res2 = data.frame(sub2, coxph.pvalPFS.ann)

##
res1 = res1[-1,]
res2 = res2[-1,]

## set range hr >18 -> 18
res1.mod = res1
se = which(res1.mod[,"ubPFS"] > 18)
res1.mod[se,"ubPFS"] = 18


##############
pdf(paste0(figures.dir, "Fig6/Ext_Fig10B_right.pdf"), useDingbats = F, 9, 9)
forestplot(res2, as.matrix(res1.mod), graph.pos = 2, new_page = F,
           xlab = "Hazard ratio",
           xticks = c(0, 1, 5, 10, 15, 18),
           graphwidth = unit(90, "mm"),
           vertices = TRUE,
           col = fpColors(lines = "black"),
           boxsize = 0.3,
           txt_gp = fpTxtGp(xlab = gpar(cex = 2.25), 
                            ticks = gpar(cex = 1.75),
                            title = gpar(cex = 2.5, fontface = "bold"),
                            cex = 1.5),
           align = c("l", "l"),
           zero = 1)
dev.off()