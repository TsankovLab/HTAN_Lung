######################################################
#################### TSANKOV LAB #####################
######################################################
# ----------------------------------------------------
# load the libraries
# ----------------------------------------------------
library(tidyverse)
# ----------------------------------------------------
# set directories
# ----------------------------------------------------
data.path <- "../data/"
figures.dir <- "../../figures/"
system(paste0("mkdir -p ", figures.dir))
# ----------------------------------------------------
# load the data
# ----------------------------------------------------
p53.filt <- c("MGH1172", "MGH1173", "MGH1175", "MGH1182", "BWH16", "BWH11", "BWH14", "BWH06")
WT.filt <- c("P14", "MGH1174", "MGH1176", "MGH1183", "BWH01", "BWH04", "BWH05", "BWH09", "BWH19", "BWH23") # no MGH1170
comparisonList = list(c("WT", "mut"))

df <- read_csv(file.path(data.path, "pathology_reports_metadata.csv")) |>
    mutate(p53_status = case_when(sampleIDs %in% p53.filt ~ "mut", sampleIDs %in% WT.filt ~ "WT", .default = NA)) |>
    drop_na(all_of("p53_status")) |>
    mutate(stageV2 = case_when(stage %in% c("IAI", "IA2", "IA3") ~ "IA", .default = stage))
# ----------------------------------------------------
# Draw contingency tables
# ----------------------------------------------------
# Stage
stage_table <- table(df$p53_status, df$stage)
print(stage_table)
# Stage V2
stagev2_table <- table(df$p53_status, df$stageV2)
print(stagev2_table)
# Simplified Stage
stageSimplified_table <- table(df$p53_status, df$stageSimplified)
print(stageSimplified_table)
# Histologic Grade
histologicGrade_table <- table(df$p53_status, df$histologicGrade)
print(histologicGrade_table)
# ----------------------------------------------------
# Run statistical tests for each variable
# ----------------------------------------------------
## Chi-squared test ----------------------------------

stage_chi <- chisq.test(stage_table)
cat("Chi-Square Test for Stage:\n")
print(stage_chi)

stagev2_chi <- chisq.test(stagev2_table)
cat("Chi-Square Test for StageV2:\n")
print(stagev2_chi)

stageSimplified_chi <- chisq.test(stageSimplified_table)
cat("Chi-Square Test for Simplified Stage:\n")
print(stageSimplified_chi)

histologicGrade_chi <- chisq.test(histologicGrade_table)
cat("Chi-Square Test for Histologic Grade:\n")
print(histologicGrade_chi)

## Fisher exact test ---------------------------------

stage_fisher <- fisher.test(stage_table)
cat("Fisher's Exact Test for Stage:\n")
print(stage_fisher)

stagev2_fisher <- fisher.test(stagev2_table)
cat("Fisher's Exact Test for Stage V2:\n")
stagev2_fisher

stageSimplified_fisher <- fisher.test(stageSimplified_table)
cat("Fisher's Exact Test for Simplified Stage:\n")
print(stageSimplified_fisher)

histologicGrade_fisher <- fisher.test(histologicGrade_table)
cat("Fisher's Exact Test for Histologic Grade:\n")
print(histologicGrade_fisher)

## Chi-squared test with Monte carlo simulation ------

# For Stage
stage_chi_mc <- chisq.test(stage_table, simulate.p.value = TRUE, B = 10000)
cat("Chi-Square Test with Monte Carlo Simulation for Stage:\n")
print(stage_chi_mc)

stagev2_chi_mc <- chisq.test(stagev2_table, simulate.p.value = TRUE, B = 10000)
cat("Chi-Square Test with Monte Carlo Simulation for Stage V2:\n")
stagev2_chi_mc

# For Simplified Stage
stageSimplified_chi_mc <- chisq.test(stageSimplified_table, simulate.p.value = TRUE, B = 10000)
cat("Chi-Square Test with Monte Carlo Simulation for Simplified Stage:\n")
print(stageSimplified_chi_mc)

# For Histologic Grade
histologicGrade_chi_mc <- chisq.test(histologicGrade_table, simulate.p.value = TRUE, B = 10000)
cat("Chi-Square Test with Monte Carlo Simulation for Histologic Grade:\n")
print(histologicGrade_chi_mc)
# ----------------------------------------------------
