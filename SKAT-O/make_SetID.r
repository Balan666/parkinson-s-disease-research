args <- commandArgs(trailingOnly = TRUE)

if (length(args) == 0) {
  stop("No input file provided. Exiting.")
}

input_file_path <- args[1]

print(paste("Processing SetID for:", input_file_path))

input_file <- read.table(input_file_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
input_file <- input_file[!duplicated(input_file$Otherinfo6), ]

output_file <- paste0("./filtered_snp/", sub(pattern = "(.*)\\..*$", replacement = "\\1", basename(input_file_path)), ".SetID")

write.table(input_file[, c("Gene.refGene", "Otherinfo6")], output_file, sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)

print(paste("SetID file saved to:", output_file))
