rm(list=ls())


library(data.table)


recode_genotype <- function(genotype, ref, alt) {
  if (is.na(genotype) || genotype == "MISSING") {
    return("./.")  
  }
  alleles <- unlist(strsplit(genotype, "/")) 
  recoded_alleles <- ifelse(alleles == ref, 0, ifelse(alleles == alt, 1, "."))  
  return(paste(recoded_alleles, collapse = "/")) 
}

process_genotype_data <- function(genotype_data, snp_info, individual_col) {
  for (i in 2:ncol(genotype_data)) {  
    snp_id <- names(genotype_data)[i]
    ref_allele <- snp_info[match(snp_id, snp_info$ID), REF]
    alt_allele <- snp_info[match(snp_id, snp_info$ID), ALT]
    genotype_data[[i]] <- sapply(genotype_data[[i]], recode_genotype, ref_allele, alt_allele)
  }
  transposed_genotype_data <- t(as.matrix(genotype_data[, -1, with = FALSE]))  
  transposed_genotype_data_df <- as.data.frame(transposed_genotype_data)
  colnames(transposed_genotype_data_df) <- genotype_data[[individual_col]]
  transposed_genotype_data_df$ID <- snp_info$ID
  final_genotype_matrix <- cbind(
    snp_info[, .(CHROM, POS, ID, REF, ALT)],  
    QUAL = rep(".", nrow(snp_info)), 
    FILTER = rep("PASS", nrow(snp_info)), 
    INFO = rep(".", nrow(snp_info)), 
    FORMAT = rep(".", nrow(snp_info)), 
    transposed_genotype_data_df[, -ncol(transposed_genotype_data_df)]  
  )
  return(final_genotype_matrix)
}

LITDLGS47031 <- fread("LITDLGS47031.csv")
snp_info <- fread("BONNE_LISTE_SNP_CHIP_27877_SNPS_beye_header2.txt")
genotype_data47031 <- LITDLGS47031[, c("INDIVIDUAL", snp_info$ID), with = FALSE]
final_genotype_matrix47031 <- process_genotype_data(genotype_data47031, snp_info, "INDIVIDUAL")
cat("First dataset processed.\n")
head(final_genotype_matrix47031[, 1:10])

LITDLGS43763 <- fread("LITDLGS43763.csv")
snps <- as.character(LITDLGS43763[9:.N, INDIVIDUAL])  
cat("Nombre de SNPs extraits :", length(snps), "\n")
genotypes <- transpose(LITDLGS43763[9:.N, -1, with = FALSE])  
individuals <- as.character(LITDLGS43763[3, -1])  
genotypes <- data.table(INTERNALCODE = individuals, genotypes)
if (length(names(genotypes)) - 1 == length(snps)) {
  setnames(genotypes, old = names(genotypes)[-1], new = snps)
} else {
  stop("Mismatch: Nombre de colonnes génotypiques et SNPs n'est pas le même.")
}
snp_cols <- grep("^PAV_", names(genotypes), value = TRUE) 
LITDLGS43763 <- genotypes[, !..snp_cols, with = FALSE]  
genotype_data43763 <- LITDLGS43763[, c("INTERNALCODE", snp_info$ID), with = FALSE]
final_genotype_matrix43763 <- process_genotype_data(genotype_data43763, snp_info, "INTERNALCODE")
cat("Second dataset processed.\n")
head(final_genotype_matrix43763[, 1:10])

fixed_cols <- c("CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT") 
sample_cols1 <- setdiff(names(final_genotype_matrix47031), fixed_cols) 
sample_cols2 <- setdiff(names(final_genotype_matrix43763), fixed_cols) 

duplicate_samples <- intersect(sample_cols1, sample_cols2)
if (length(duplicate_samples) > 0) {
  cat("Les échantillons suivants sont en double :\n", paste(duplicate_samples, collapse = ", "), "\n")
  final_genotype_matrix43763 <- final_genotype_matrix43763[, !duplicate_samples, with = FALSE]
}


merged_data <- merge(final_genotype_matrix47031, final_genotype_matrix43763, by = fixed_cols, all = TRUE)

fwrite(merged_data, "merged_LITDLGS47031_LITDLGS43763.csv")
cat("Merged dataset exported.\n")