################################################################
##### Obtain APPRIS transcript data from genes of interest #####
################################################################

## Load in packages
library(data.table)
library(dplyr)

## Function to read and manipulate appris data
appris <- function(appris_file, genes) {
  appris_colnames <- c("gene_name", "gene_id", "transcript_id", "CCDS", "status")
  appris <- fread(appris_file, header=FALSE) %>% `colnames<-` (appris_colnames)
  
  ## Get the transcript ids from the genes of interest 
  appris_genes <- appris[gene_name %in% genes]

  return(appris_genes)
}