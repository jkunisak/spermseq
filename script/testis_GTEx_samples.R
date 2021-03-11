#############################################################################
##### Function to get testis samples from the GTEx transcript tpm file #####
#############################################################################
## Load in packages
library(data.table)
library(dplyr)

## testis_GTEx_samples function
testis_GTEx_samples <- function(GTEx_file=GTEx_file, GTEx_samples=GTEx_samples) {
  
  ## Read in the GTEx expression data
  df <- fread(GTEx_file, header = TRUE, sep = "\t")
  
  ## Read in the GTEx samples of interest 
  cols <- c("transcript_id", "gene_id", fread(GTEx_samples, header=FALSE)$V1)
  
  GTEx_testis_data <- df[, which(colnames(df) %in% cols), with=FALSE]
  
  return(GTEx_testis_data)
}


