#####################################################################
##### Obtain Ensembl GTF data from APPRIS principal transcripts #####
#####################################################################

## Load in packages
library(data.table)
library(dplyr)
library(pbmcapply)

## Function to read and manipulate appris data
gtf_manipulation <- function(gtf, df) {
  #browser()
  
  ## Convert df to list 
  if (class(df)[1] != "list") {
    df <- df %>% transpose(.) %>% as.list(.) %>% unname(.)
  }
  r <- df[[1]]
  
  final <- rbindlist(pbmclapply(mc.cores = 8, df, function(r, gtf) {
    ## Print the row
    #print(r)
    
    ## Get gene_name, gene_id, and transcript_id information
    gene_name <- r[1]
    gene_id <- r[grep("ENSG", r)]
    transcript_id <- r[grep("ENST", r)]
    
    ## Get the gtf information on the transcript of interest 
    gtf_r <- gtf[grep(pattern = transcript_id, x = gtf$info)]
    
    ## Get the exon lengths and total exon lengths
    if (nrow(gtf_r) > 0) {
      ## Get the protein_id (use for mapping protein coordinates to genomic coordinates)
      protein_id <- ((gsub(pattern = "\"; ", replacement = "|", x = gtf_r$info[1]) %>% 
                        gsub(pattern = " \"", replacement = ": ", x = .) %>% strsplit(., split = "|", fixed = TRUE))[[1]] %>%
                       .[grep(pattern = "protein_id", x = .)] %>% strsplit(., split = "protein_id: "))[[1]][2]
      
      ## Get the necessary information for the transcript
      transcript_df <- gtf_r[, c("chr", "start", "stop")] %>% 
        mutate(transcript_id = transcript_id, gene_name = gene_name, gene_id = gene_id, protein_id=protein_id)
      
      ## Get the exon lengths 
      transcript_df$exon_length <- transcript_df$stop - transcript_df$start + 1
      
      ## Get the total exon length
      transcript_df$total_exon_length <- sum(transcript_df$exon_length)
      
      ## Get the total number of exons for the transcript
      transcript_df$exon_num <- nrow(transcript_df)
    }
    
    ## Fill in data if transcript_df is empty (perhaps non-coding?); not in gtf file
    if (nrow(gtf_r) == 0) {
      transcript_df <- data.table("chr"=0, "start"=0, "stop"=0,
                                  "transcript_id"=transcript_id,
                                  "gene_name"=gene_name,
                                  "gene_id"=gene_id,
                                  "protein_id"="nonCDS",
                                  "exon_length"=0,
                                  "total_exon_length"=0,
                                  "exon_num"=0)
    }
    
    ## Keep the APPRIs and GTEx data
    temp <- data.table("CCDS"=r[4],
                       "status"=r[5],
                       "avg_GTEx_TPM"=as.numeric(r[6]),
                       "appris"=as.numeric(r[7]),
                       "GTEx"=as.numeric(r[8]))
    
    #print(cbind(transcript_df, temp))
    return(cbind(transcript_df, temp))

  }, gtf=gtf))
  
  return(final)
}