#!/usr/bin/env Rscript

library(SKAT)
require(methods)
library(SPAtest)
library(RSpectra)

args <- commandArgs(trailingOnly = TRUE)

if (length(args) == 0) {
  stop("No input prefix provided. Exiting.")
}

wd <- "~/runs/senkkon/2025/AMP_PD_PRS/SKAT-O_run_ver2/"
setwd(wd)

cohort <- args[1]  # Input prefix

# Define file paths
File.Bed   <- paste0("filtered_snp/", cohort, ".bed")
File.Bim   <- paste0("filtered_snp/", cohort, ".bim")
File.Fam   <- paste0("filtered_snp/", cohort, ".fam")
File.SetID <- paste0("filtered_snp/", cohort, ".SetID")
File.SSD   <- paste0("filtered_snp/", cohort, ".SSD")
File.Info  <- paste0("filtered_snp/", cohort, ".info")
File.Cov   <- paste0("filtered_snp/covar_", cohort, ".txt")
File.Results.SKATO  <- paste0("SKAT_AMP_PD/", cohort, ".results.skato")
File.Results.BURDEN <- paste0("SKAT_AMP_PD/", cohort, ".results.burden")

# Generate SSD and SetID
Generate_SSD_SetID(File.Bed, File.Bim, File.Fam, File.SetID, File.SSD, File.Info)

# Open SSD file
SSD.INFO <- Open_SSD(File.SSD, File.Info)

# Print summary
print(paste("Samples:", SSD.INFO$nSample, "Sets:", SSD.INFO$nSets))

# Check for covariates
if (file.exists(File.Cov)) {
  message("Running SKAT-O with covariates")
  FAM <- Read_Plink_FAM_Cov(File.Fam, File.Cov, Is.binary = TRUE, cov_header = TRUE)
  y <- FAM$Phenotype
  Age <- FAM$Age
  Sex <- FAM$Sex.y
  obj <- SKAT_Null_Model(y ~ Sex + Age, out_type = "D")
} else {
  message("Running SKAT-O without covariates")
  FAM <- Read_Plink_FAM(File.Fam, Is.binary = TRUE)
  y <- FAM$Phenotype
  obj <- SKAT_Null_Model(y ~ 1, out_type = "D")
}

# Run SKAT-O
out.skato <- SKATBinary.SSD.All(SSD.INFO, obj, method = "optimal.adj")
write.table(out.skato$results, file = File.Results.SKATO, col.names = TRUE, row.names = FALSE)

print(paste("SKAT-O results saved to:", File.Results.SKATO))

