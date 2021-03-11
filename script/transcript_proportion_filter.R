########################################################################################################
##### Function to remove exons that are below threshold proportion of transcripts for a given gene #####
########################################################################################################
## Load in necessary packages
library(data.table)
library(dplyr)


transcript_proportion_filter <- function(merged_exons) {
  ## Get the genes
  genes <- unique(merged_exons$gene_id)
  gene <- genes[7]
  
  ## Go through each gene
  final <- rbindlist(lapply(genes, function(gene, merged_exons) {
    ## Store merged_exon data as temp dataset
    temp <- merged_exons[gene_id == gene]
    
    ## Get the number of transcripts for the given gene
    total_transcript_num <- (temp[gene_id == gene]$transcript_id %>% strsplit(x = ., split = " | ", fixed = TRUE)) %>% unlist() %>% unique() %>% length()
    temp$total_transcript_num <- total_transcript_num
    
    ## Get the proportion of transcripts the exon is found in 
    temp$transcript_proportion <- temp$transcript_num / total_transcript_num
    
    return(temp)
  }, merged_exons=merged_exons)) 
  
  return(final)
}