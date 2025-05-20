library(dplyr)
library(tidyr)
library(purrr)
library(ggplot2)

results_dir <- "../PRS_MAF_01"
all_chr_aggregated <- data.frame()
PREFIX <- 'PRS_bottom_25'
CHR <- 2
CHR <- 3
for (PREFIX in c('PRS_bottom_25', 'PRS_top_25', 'PRS_b_case_top_control_chr')) { 
  
  for (CHR in 1:22) {
    merged <- data.frame()
    
    pattern <- paste0(PREFIX, "_chr", CHR, ".(CADD|LOF|nonsyn|ALL|ENCODE)\\.results\\.skato")
    files <- list.files(results_dir, pattern = pattern, full.names = TRUE)
    
    if (length(files) >= 4) {
      group_names <- sub(".*\\.(.*?)\\.results\\.skato", "\\1", files)
      
      tables <- map2(files, group_names, function(file, group) {
        df <- read.table(file, header = TRUE, sep = " ")
        if (!is.null(df) && "SetID" %in% names(df) && "P.value" %in% names(df)) {
          df <- df[, c("SetID", "P.value")]
          colnames(df)[2] <- group
          return(df)
        } else {
          return(NULL)
        }
      })
      
      tables <- tables[!sapply(tables, is.null)]
      
      if (length(tables) >= 2) {
        merged <- reduce(tables, full_join, by = "SetID")
        
        for (group in c('ALL', 'CADD', 'ENCODE', 'LOF', 'nonsyn')) {
          if (group %in% names(merged)) {
            merged[[paste0(group, "_FDR")]] <- p.adjust(merged[[group]], method = "fdr")
          }
        }
        
        merged$CHR <- CHR
        merged$PREFIX <- PREFIX
        all_chr_aggregated <- rbind(all_chr_aggregated, merged)
        
      } else {
        cat("Not enough valid tables to merge for", PREFIX, "chr", CHR, "\n")
      }
      
    } else {
      cat("Not enough files found for", PREFIX, "chr", CHR, "\n")
    }
  }
  filtered_results <- all_chr_aggregated[
    apply(all_chr_aggregated[ , grepl("_FDR$", names(all_chr_aggregated))], 1, function(row) any(row < 0.05, na.rm = TRUE)),
  ]
  write.table(filtered_results, paste0("aggregated/", PREFIX, "_filtered_table.txt"), sep = "\t", row.names = FALSE, quote = FALSE)
  
  output_file <- paste0("aggregated/", PREFIX, "_multi_group_fdr_table.txt")
  write.table(all_chr_aggregated, output_file, sep = "\t", row.names = FALSE, quote = FALSE)
}
