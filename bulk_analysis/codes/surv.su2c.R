#####################@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ 
#####################@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
##################### FIGURE 6E AND EXTENDED 10B 
##################### SURVIVAL ANALYSES OF SU2C-MARK COHORT
library(readxl)
library(survival)
library(forestplot)

# my outputs
data.path <- "./data/"
figures.dir <- "../../figures/"

if (!file.exists(paste0(figures.dir,'Fig6/'))){dir.create(paste0(figures.dir,'Fig6/'), recursive=TRUE)}

#####################
#####################
clin = read_excel(paste0(data.path, "su2c_clinical_info.xlsx"))
clin = as.data.frame(clin)
clin = clin[-1,]
colnames(clin) = clin[1,]
clin = clin[-1,]

clin$OS_Years = as.numeric(clin$Harmonized_OS_Days) / 365
clin$PFS_Years = as.numeric(clin$Harmonized_PFS_Days) / 365
clin$OS_Status = as.numeric(clin$Harmonized_OS_Event)
clin$PFS_Status = as.numeric(clin$Harmonized_PFS_Event)
clin.bkup1 = clin
se = which(clin$WES_All == 1)
clin = clin[se,]
rownames(clin) = clin$Harmonized_SU2C_WES_Tumor_Sample_ID_v2

########
p53 = get(load(paste0(data.path, "su2c.mut.info.p53.Rda"))
mymut = get(load(paste0(data.path, "su2c.mut.info.selectgenes1.Rda"))

old.clin = clin
all(rownames(clin) == rownames(p53))
all(rownames(clin) == rownames(mymut))
clin = cbind(clin, mymut[,-c(3,4)], p53)

######## convert silent mutations to WT
se = which(clin$EGFR.all == "Silent")
clin$EGFR_status[se] = "Wildtype"
se = which(clin$KEAP1.all == "Silent")
clin$KEAP1_status[se] = "Wildtype"
se = which(clin$KRAS.all == "Silent")
clin$KRAS_status[se] = "Wildtype"
se = which(clin$STK11.all == "Silent")  ## 0
se = which(clin$TP53.all == "Silent")  ## 0



##################### 6E (left)
##################### adeno + anti-pd1
se = which(clin$Agent_PD1 %in% c("Nivolumab", "Pembrolizumab") & clin$Histology_Harmonized == "Adeno")  ## adeno+ anti-pd-1
clin.subset = clin[se,]

### mut of common genes
clin.subset1 = clin.subset[,c("EGFR_status", "STK11_status", "KEAP1_status", "KRAS_status", 
                              "TP53_status", "OS_Years", "PFS_Years", "OS_Status", "PFS_Status")]

sub = colnames(clin.subset1)[1:5]

survdiff.pvalOS = coxph.pvalOS = survdiff.pvalPFS = coxph.pvalPFS = rep(0, length(sub))
hrOS = hrPFS = lbOS = ubOS = lbPFS = ubPFS = num_high = num_low = rep(0, length(sub))

for (k in 1:length(sub))
{
  cat("\r", k)
  
  label = rep(NA, nrow(clin.subset1))
  se = which(clin.subset1[,k] == "Wildtype")
  label[se] = 0
  se = which(clin.subset1[,k] == "Mutants")
  label[se] = 1
  
  xx = cbind(clin.subset1, label)
  se = !is.na(xx$label)
  xx = xx[se,]
  
  a = table(label)
  num_high[k] = as.numeric(a[names(a)==1])
  num_low[k] = as.numeric(a[names(a)==0])
  
  mycox = coxph(Surv(OS_Years, OS_Status)~label, xx)
  mycox = summary(mycox)
  coxph.pvalOS[k] = mycox$coefficients[5]
  tmp = mycox$conf.int
  hrOS[k] = tmp[1]
  lbOS[k] = tmp[3]
  ubOS[k] = tmp[4]
  mycox = coxph(Surv(PFS_Years, PFS_Status)~label, xx)
  mycox = summary(mycox)
  coxph.pvalPFS[k] = mycox$coefficients[5]
  tmp = mycox$conf.int
  hrPFS[k] = tmp[1]
  lbPFS[k] = tmp[3]
  ubPFS[k] = tmp[4]
  
}

res.surv = data.frame(sub, num_high, num_low, coxph.pvalOS, hrOS, lbOS, ubOS, coxph.pvalPFS, hrPFS, lbPFS, ubPFS)
res.surv1 = res.surv


### mut of p53 specific subsets

se.list = list()
se.list[[1]] = which(clin.subset$Impact1 == 1)
se.list[[2]] = which(clin.subset$Impact2 == 1)
se.list[[3]] = which(clin.subset$DNE_LOFclass == "DNE_LOF")
se.list[[4]] = which(clin.subset$DNE_LOFclass == "notDNE_LOF")
se.list[[5]] = which(clin.subset$DNE_LOFclass == "notDNE_notLOF")
se.list[[6]] = grep("Missense_Mutation", clin.subset$TP53.all)
se.list[[7]] = grep("Nonsense_Mutation", clin.subset$TP53.all)
se.list[[8]] = grep("Splice_Site", clin.subset$TP53.all)
se.list[[9]] = which(clin.subset$TP53.all == "Splice_Site")
se.list[[10]] = grep("Frame_Shift_Ins|Frame_Shift_Del", clin.subset$TP53.all)


sub = c("Impact1", "Impact2", "DNE_LOF", "notDNE_LOF", "notDNE_notLOF",
        "Missense", "Nonsense", "Splice_Site_All", "Splice_Site_Unique", "Frameshift")

survdiff.pvalOS = coxph.pvalOS = survdiff.pvalPFS = coxph.pvalPFS = rep(0, length(sub))
hrOS = hrPFS = lbOS = ubOS = lbPFS = ubPFS = num_high = num_low = rep(0, length(sub))

for (k in 1:length(sub))
{
  cat("\r", k)
  
  label = rep(NA, nrow(clin.subset))
  se = which(clin.subset$TP53.all == "Wildtype")
  label[se] = 0
  se = se.list[[k]]
  label[se] = 1
  
  xx = cbind(clin.subset, label)
  se = !is.na(xx$label)
  xx = xx[se,]
  
  a = table(label)
  num_high[k] = as.numeric(a[names(a)==1])
  num_low[k] = as.numeric(a[names(a)==0])
  
  mycox = coxph(Surv(OS_Years, OS_Status)~label, xx)
  mycox = summary(mycox)
  coxph.pvalOS[k] = mycox$coefficients[5]
  tmp = mycox$conf.int
  hrOS[k] = tmp[1]
  lbOS[k] = tmp[3]
  ubOS[k] = tmp[4]
  mycox = coxph(Surv(PFS_Years, PFS_Status)~label, xx)
  mycox = summary(mycox)
  coxph.pvalPFS[k] = mycox$coefficients[5]
  tmp = mycox$conf.int
  hrPFS[k] = tmp[1]
  lbPFS[k] = tmp[3]
  ubPFS[k] = tmp[4]
  
}

res.surv = data.frame(sub, num_high, num_low, coxph.pvalOS, hrOS, lbOS, ubOS, coxph.pvalPFS, hrPFS, lbPFS, ubPFS)
res.surv2 = res.surv

###
res.surv = rbind(res.surv1, res.surv2)

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
##################### adeno + all treatments
se = which(clin$Histology_Harmonized == "Adeno") ## adeno + all treatments
clin.subset = clin[se,]

### mut of common genes
clin.subset1 = clin.subset[,c("EGFR_status", "STK11_status", "KEAP1_status", "KRAS_status", 
                              "TP53_status", "OS_Years", "PFS_Years", "OS_Status", "PFS_Status")]

sub = colnames(clin.subset1)[1:5]

survdiff.pvalOS = coxph.pvalOS = survdiff.pvalPFS = coxph.pvalPFS = rep(0, length(sub))
hrOS = hrPFS = lbOS = ubOS = lbPFS = ubPFS = num_high = num_low = rep(0, length(sub))

for (k in 1:length(sub))
{
  cat("\r", k)
  
  label = rep(NA, nrow(clin.subset1))
  se = which(clin.subset1[,k] == "Wildtype")
  label[se] = 0
  se = which(clin.subset1[,k] == "Mutants")
  label[se] = 1
  
  xx = cbind(clin.subset1, label)
  se = !is.na(xx$label)
  xx = xx[se,]
  
  a = table(label)
  num_high[k] = as.numeric(a[names(a)==1])
  num_low[k] = as.numeric(a[names(a)==0])
  
  mycox = coxph(Surv(OS_Years, OS_Status)~label, xx)
  mycox = summary(mycox)
  coxph.pvalOS[k] = mycox$coefficients[5]
  tmp = mycox$conf.int
  hrOS[k] = tmp[1]
  lbOS[k] = tmp[3]
  ubOS[k] = tmp[4]
  mycox = coxph(Surv(PFS_Years, PFS_Status)~label, xx)
  mycox = summary(mycox)
  coxph.pvalPFS[k] = mycox$coefficients[5]
  tmp = mycox$conf.int
  hrPFS[k] = tmp[1]
  lbPFS[k] = tmp[3]
  ubPFS[k] = tmp[4]
  
}

res.surv = data.frame(sub, num_high, num_low, coxph.pvalOS, hrOS, lbOS, ubOS, coxph.pvalPFS, hrPFS, lbPFS, ubPFS)
res.surv1 = res.surv


### mut of p53 specific subsets

se.list = list()
se.list[[1]] = which(clin.subset$Impact1 == 1)
se.list[[2]] = which(clin.subset$Impact2 == 1)
se.list[[3]] = which(clin.subset$DNE_LOFclass == "DNE_LOF")
se.list[[4]] = which(clin.subset$DNE_LOFclass == "notDNE_LOF")
se.list[[5]] = which(clin.subset$DNE_LOFclass == "notDNE_notLOF")
se.list[[6]] = grep("Missense_Mutation", clin.subset$TP53.all)
se.list[[7]] = grep("Nonsense_Mutation", clin.subset$TP53.all)
se.list[[8]] = grep("Splice_Site", clin.subset$TP53.all)
se.list[[9]] = which(clin.subset$TP53.all == "Splice_Site")
se.list[[10]] = grep("Frame_Shift_Ins|Frame_Shift_Del", clin.subset$TP53.all)


sub = c("Impact1", "Impact2", "DNE_LOF", "notDNE_LOF", "notDNE_notLOF",
        "Missense", "Nonsense", "Splice_Site_All", "Splice_Site_Unique", "Frameshift")

survdiff.pvalOS = coxph.pvalOS = survdiff.pvalPFS = coxph.pvalPFS = rep(0, length(sub))
hrOS = hrPFS = lbOS = ubOS = lbPFS = ubPFS = num_high = num_low = rep(0, length(sub))

for (k in 1:length(sub))
{
  cat("\r", k)
  
  label = rep(NA, nrow(clin.subset))
  se = which(clin.subset$TP53.all == "Wildtype")
  label[se] = 0
  se = se.list[[k]]
  label[se] = 1
  
  xx = cbind(clin.subset, label)
  se = !is.na(xx$label)
  xx = xx[se,]
  
  a = table(label)
  num_high[k] = as.numeric(a[names(a)==1])
  num_low[k] = as.numeric(a[names(a)==0])
  
  mycox = coxph(Surv(OS_Years, OS_Status)~label, xx)
  mycox = summary(mycox)
  coxph.pvalOS[k] = mycox$coefficients[5]
  tmp = mycox$conf.int
  hrOS[k] = tmp[1]
  lbOS[k] = tmp[3]
  ubOS[k] = tmp[4]
  mycox = coxph(Surv(PFS_Years, PFS_Status)~label, xx)
  mycox = summary(mycox)
  coxph.pvalPFS[k] = mycox$coefficients[5]
  tmp = mycox$conf.int
  hrPFS[k] = tmp[1]
  lbPFS[k] = tmp[3]
  ubPFS[k] = tmp[4]
  
}

res.surv = data.frame(sub, num_high, num_low, coxph.pvalOS, hrOS, lbOS, ubOS, coxph.pvalPFS, hrPFS, lbPFS, ubPFS)
res.surv2 = res.surv

###
res.surv = rbind(res.surv1, res.surv2)

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
##################### squamous + anti-pd1
se = which(clin$Agent_PD1 %in% c("Nivolumab", "Pembrolizumab") & clin$Histology_Harmonized == "Squamous")  ## sq+ anti-pd-1
clin.subset = clin[se,]

### mut of common genes
clin.subset1 = clin.subset[,c("EGFR_status", "STK11_status", "KEAP1_status", "KRAS_status", 
                              "TP53_status", "OS_Years", "PFS_Years", "OS_Status", "PFS_Status")]

sub = colnames(clin.subset1)[1:5]

survdiff.pvalOS = coxph.pvalOS = survdiff.pvalPFS = coxph.pvalPFS = rep(0, length(sub))
hrOS = hrPFS = lbOS = ubOS = lbPFS = ubPFS = num_high = num_low = rep(0, length(sub))

for (k in 2:length(sub))  ## egfr has no muts (pdl1), only 1 mut (all treatments)
{
  cat("\r", k)
  
  label = rep(NA, nrow(clin.subset1))
  se = which(clin.subset1[,k] == "Wildtype")
  label[se] = 0
  se = which(clin.subset1[,k] == "Mutants")
  label[se] = 1
  
  xx = cbind(clin.subset1, label)
  se = !is.na(xx$label)
  xx = xx[se,]
  
  a = table(label)
  num_high[k] = as.numeric(a[names(a)==1])
  num_low[k] = as.numeric(a[names(a)==0])
  
  mycox = coxph(Surv(OS_Years, OS_Status)~label, xx)
  mycox = summary(mycox)
  coxph.pvalOS[k] = mycox$coefficients[5]
  tmp = mycox$conf.int
  hrOS[k] = tmp[1]
  lbOS[k] = tmp[3]
  ubOS[k] = tmp[4]
  mycox = coxph(Surv(PFS_Years, PFS_Status)~label, xx)
  mycox = summary(mycox)
  coxph.pvalPFS[k] = mycox$coefficients[5]
  tmp = mycox$conf.int
  hrPFS[k] = tmp[1]
  lbPFS[k] = tmp[3]
  ubPFS[k] = tmp[4]
  
}

res.surv = data.frame(sub, num_high, num_low, coxph.pvalOS, hrOS, lbOS, ubOS, coxph.pvalPFS, hrPFS, lbPFS, ubPFS)
res.surv1 = res.surv

### mut of p53 specific subsets

se.list = list()
se.list[[1]] = which(clin.subset$Impact1 == 1)
se.list[[2]] = which(clin.subset$Impact2 == 1)
se.list[[3]] = which(clin.subset$DNE_LOFclass == "DNE_LOF")
se.list[[4]] = which(clin.subset$DNE_LOFclass == "notDNE_LOF")
se.list[[5]] = which(clin.subset$DNE_LOFclass == "notDNE_notLOF")
se.list[[6]] = grep("Missense_Mutation", clin.subset$TP53.all)
se.list[[7]] = grep("Nonsense_Mutation", clin.subset$TP53.all)
se.list[[8]] = grep("Splice_Site", clin.subset$TP53.all)
se.list[[9]] = which(clin.subset$TP53.all == "Splice_Site")
se.list[[10]] = grep("Frame_Shift_Ins|Frame_Shift_Del", clin.subset$TP53.all)


sub = c("Impact1", "Impact2", "DNE_LOF", "notDNE_LOF", "notDNE_notLOF",
        "Missense", "Nonsense", "Splice_Site_All", "Splice_Site_Unique", "Frameshift")

survdiff.pvalOS = coxph.pvalOS = survdiff.pvalPFS = coxph.pvalPFS = rep(0, length(sub))
hrOS = hrPFS = lbOS = ubOS = lbPFS = ubPFS = num_high = num_low = rep(0, length(sub))

for (k in 1:length(sub))
{
  cat("\r", k)
  
  label = rep(NA, nrow(clin.subset))
  se = which(clin.subset$TP53.all == "Wildtype")
  label[se] = 0
  se = se.list[[k]]
  label[se] = 1
  
  xx = cbind(clin.subset, label)
  se = !is.na(xx$label)
  xx = xx[se,]
  
  a = table(label)
  num_high[k] = as.numeric(a[names(a)==1])
  num_low[k] = as.numeric(a[names(a)==0])
  
  mycox = coxph(Surv(OS_Years, OS_Status)~label, xx)
  mycox = summary(mycox)
  coxph.pvalOS[k] = mycox$coefficients[5]
  tmp = mycox$conf.int
  hrOS[k] = tmp[1]
  lbOS[k] = tmp[3]
  ubOS[k] = tmp[4]
  mycox = coxph(Surv(PFS_Years, PFS_Status)~label, xx)
  mycox = summary(mycox)
  coxph.pvalPFS[k] = mycox$coefficients[5]
  tmp = mycox$conf.int
  hrPFS[k] = tmp[1]
  lbPFS[k] = tmp[3]
  ubPFS[k] = tmp[4]
  
}

res.surv = data.frame(sub, num_high, num_low, coxph.pvalOS, hrOS, lbOS, ubOS, coxph.pvalPFS, hrPFS, lbPFS, ubPFS)
res.surv2 = res.surv

###
res.surv = rbind(res.surv1, res.surv2)

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
##################### 
se = which(clin$Histology_Harmonized == "Squamous") ## sq + all treatments
clin.subset = clin[se,]

### mut of common genes
clin.subset1 = clin.subset[,c("EGFR_status", "STK11_status", "KEAP1_status", "KRAS_status", 
                              "TP53_status", "OS_Years", "PFS_Years", "OS_Status", "PFS_Status")]

sub = colnames(clin.subset1)[1:5]

survdiff.pvalOS = coxph.pvalOS = survdiff.pvalPFS = coxph.pvalPFS = rep(0, length(sub))
hrOS = hrPFS = lbOS = ubOS = lbPFS = ubPFS = num_high = num_low = rep(0, length(sub))

for (k in 2:length(sub))  ## egfr has no muts (pdl1), only 1 mut (all treatments)
{
  cat("\r", k)
  
  label = rep(NA, nrow(clin.subset1))
  se = which(clin.subset1[,k] == "Wildtype")
  label[se] = 0
  se = which(clin.subset1[,k] == "Mutants")
  label[se] = 1
  
  xx = cbind(clin.subset1, label)
  se = !is.na(xx$label)
  xx = xx[se,]
  
  a = table(label)
  num_high[k] = as.numeric(a[names(a)==1])
  num_low[k] = as.numeric(a[names(a)==0])
  
  mycox = coxph(Surv(OS_Years, OS_Status)~label, xx)
  mycox = summary(mycox)
  coxph.pvalOS[k] = mycox$coefficients[5]
  tmp = mycox$conf.int
  hrOS[k] = tmp[1]
  lbOS[k] = tmp[3]
  ubOS[k] = tmp[4]
  mycox = coxph(Surv(PFS_Years, PFS_Status)~label, xx)
  mycox = summary(mycox)
  coxph.pvalPFS[k] = mycox$coefficients[5]
  tmp = mycox$conf.int
  hrPFS[k] = tmp[1]
  lbPFS[k] = tmp[3]
  ubPFS[k] = tmp[4]
  
}

res.surv = data.frame(sub, num_high, num_low, coxph.pvalOS, hrOS, lbOS, ubOS, coxph.pvalPFS, hrPFS, lbPFS, ubPFS)
res.surv1 = res.surv

### mut of p53 specific subsets

se.list = list()
se.list[[1]] = which(clin.subset$Impact1 == 1)
se.list[[2]] = which(clin.subset$Impact2 == 1)
se.list[[3]] = which(clin.subset$DNE_LOFclass == "DNE_LOF")
se.list[[4]] = which(clin.subset$DNE_LOFclass == "notDNE_LOF")
se.list[[5]] = which(clin.subset$DNE_LOFclass == "notDNE_notLOF")
se.list[[6]] = grep("Missense_Mutation", clin.subset$TP53.all)
se.list[[7]] = grep("Nonsense_Mutation", clin.subset$TP53.all)
se.list[[8]] = grep("Splice_Site", clin.subset$TP53.all)
se.list[[9]] = which(clin.subset$TP53.all == "Splice_Site")
se.list[[10]] = grep("Frame_Shift_Ins|Frame_Shift_Del", clin.subset$TP53.all)


sub = c("Impact1", "Impact2", "DNE_LOF", "notDNE_LOF", "notDNE_notLOF",
        "Missense", "Nonsense", "Splice_Site_All", "Splice_Site_Unique", "Frameshift")

survdiff.pvalOS = coxph.pvalOS = survdiff.pvalPFS = coxph.pvalPFS = rep(0, length(sub))
hrOS = hrPFS = lbOS = ubOS = lbPFS = ubPFS = num_high = num_low = rep(0, length(sub))

for (k in 1:length(sub))
{
  cat("\r", k)
  
  label = rep(NA, nrow(clin.subset))
  se = which(clin.subset$TP53.all == "Wildtype")
  label[se] = 0
  se = se.list[[k]]
  label[se] = 1
  
  xx = cbind(clin.subset, label)
  se = !is.na(xx$label)
  xx = xx[se,]
  
  a = table(label)
  num_high[k] = as.numeric(a[names(a)==1])
  num_low[k] = as.numeric(a[names(a)==0])
  
  mycox = coxph(Surv(OS_Years, OS_Status)~label, xx)
  mycox = summary(mycox)
  coxph.pvalOS[k] = mycox$coefficients[5]
  tmp = mycox$conf.int
  hrOS[k] = tmp[1]
  lbOS[k] = tmp[3]
  ubOS[k] = tmp[4]
  mycox = coxph(Surv(PFS_Years, PFS_Status)~label, xx)
  mycox = summary(mycox)
  coxph.pvalPFS[k] = mycox$coefficients[5]
  tmp = mycox$conf.int
  hrPFS[k] = tmp[1]
  lbPFS[k] = tmp[3]
  ubPFS[k] = tmp[4]
  
}

res.surv = data.frame(sub, num_high, num_low, coxph.pvalOS, hrOS, lbOS, ubOS, coxph.pvalPFS, hrPFS, lbPFS, ubPFS)
res.surv2 = res.surv

###
res.surv = rbind(res.surv1, res.surv2)

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