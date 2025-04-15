library(dplyr)

args <- commandArgs(trailingOnly = TRUE)

if (length(args) == 0) {
  stop("No .fam file provided for covariance filtering. Exiting.")
}

# Set working directory
wd <- '~/runs/senkkon/2025/AMP_PD_PRS/SKAT-O_run_ver2/'
setwd(wd)

# Path to covariates file
covariates_file <- "~/runs/senkkon/2025/AMP_PD_PRS/covar_AMP_PD.txt"

# Read covariates file
covariates <- read.table(covariates_file, header = TRUE, stringsAsFactors = FALSE)

# Process each provided .fam file
for (fam_file in args) {
  print(paste("Processing:", fam_file))
  
  input_file <- read.table(fam_file, header = FALSE, stringsAsFactors = FALSE)
  
  out_fam <- input_file %>%
    select(FID = V1, IID = V2) %>%
    mutate(FID = trimws(FID), IID = trimws(IID))  # Remove extra spaces
  
  out_covariates <- covariates %>%
    mutate(FID = trimws(FID), IID = trimws(IID)) %>%
    inner_join(out_fam, by = c("FID", "IID"))
  
  out_file <- paste0("./filtered_snp/covar_", sub(pattern = "(.*)\\..*$", replacement = "\\1", basename(fam_file)), ".txt")
  write.table(out_covariates, file = out_file, quote = FALSE, row.names = FALSE, sep = "\t")
  
  print(paste("Filtered covariates saved to:", out_file))
}


