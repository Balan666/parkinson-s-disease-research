#based on PRS scores filter 25% top and bottom genomes

library(data.table)

prs_data <- fread("~/runs/senkkon/2025/AMP_PD_PRS/SKAT-O_run_ver1/AMP_PD_plink_ver2_res.0.4.profile")

# Sort by PRS scores
prs_data <- prs_data[order(SCORE)]

# Split into top 25% and bottom 25%
n <- nrow(prs_data)
top_25 <- prs_data[(0.75 * n):n, ]
bottom_25 <- prs_data[1:(0.25 * n), ]

write.table(top_25[, .(FID, IID)], "~/runs/senkkon/2025/AMP_PD_PRS/SKAT-O_run_ver1/inter_files/PRS_top_25_fids.txt", col.names = FALSE, row.names = FALSE, sep = "\t")
write.table(bottom_25[, .(FID, IID)], "~/runs/senkkon/2025/AMP_PD_PRS/SKAT-O_run_ver1/inter_files/PRS_bottom_25_fids.txt", col.names = FALSE, row.names = FALSE, sep = "\t")