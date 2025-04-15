library(dplyr)
library(tidyr)
library(purrr)

results_dir <- "/home/balan/Documents/PRS/results/SKAT-O_AMP_PD/results/RESULTS"

for (PREFIX in c('PRS_bottom_25', 'PRS_top_25')) {  
  all_chr_aggregated <- data.frame()
  
  for (CHR in 1:22) {
    
    pattern <- paste0(PREFIX, "_chr", CHR, ".*\\.results\\.skato")
    files <- list.files(results_dir, pattern = pattern, full.names = TRUE)
    
    if (length(files) >= 4) {
      group_names <- sub(".*\\.(.*?)\\.results\\.skato", "\\1", files)
      
      tables <- map2(files, group_names, function(file, group) {
        df <- tryCatch(read.table(file, header = TRUE, sep = " "), error = function(e) NULL)
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
        
        for (group in c('CADD', 'ENCODE', 'LOF', 'nonsyn')) {
          if (group %in% names(merged)) {
            merged[[paste0(group, "_FDR")]] <- p.adjust(merged[[group]], method = "fdr")
          }
        }
        
        merged$CHR <- CHR
        merged$PREFIX <- PREFIX
        all_chr_aggregated <- rbind(all_chr_aggregated, merged)
        all_chr_aggregated <- all_chr_aggregated %>%
          arrange(all_chr_aggregated[2:5])
        all_chr_aggregated <- na.omit(all_chr_aggregated) #%>%
          #filter((CADD.x < 0.05) & (ENCODE.x < 0.05)  & (nonsyn.x < 0.05))
        #all_chr_aggregated <- all_chr_aggregated[rowSums(is.na(all_chr_aggregated)) != ncol(all_chr_aggregated),]
        
      } else {
        cat("Not enough valid tables to merge for", PREFIX, "chr", CHR, "\n")
      }
      
    } else {
      cat("Not enough files found for", PREFIX, "chr", CHR, "\n")
    }
  }
  
  output_file <- paste0("/home/balan/Documents/PRS/results/SKAT-O_AMP_PD/", PREFIX, "_multi_group_fdr_table.txt")
  #write.table(all_chr_aggregated, output_file, sep = "\t", row.names = FALSE, quote = FALSE)
}

