###############################################################
##### Obtain GTEx testes transcript-level expression data #####
###############################################################

## Load in packages
library(data.table)
library(dplyr)
library(beeswarm)
library(biomaRt)

## Function to read and manipulate appris and GTEx data
gtex_transcript_tpm_manipulation <- function(appris_df_tpm, testis_GTEx_tpm, moderate_gene_ids) {
  
  ## Get all of the gene_ids in the entire testis GTEx tpm file
  gene_id_GTEx <- unique(testis_GTEx_tpm$gene_id) %>% gsub("\\..*", "", .)
  
  ## Remove moderate gene ids 
  if (length(moderate_gene_ids) > 0) {
    gene_id_GTEx <- gene_id_GTEx[! gene_id_GTEx %in% moderate_gene_ids]
  }
  
  ## Go through each gene_id and get the TPM values for GTEx only transcripts that are highly expressed
  #g_id="ENSG00000163041"
  appris_gtex_transcript_tpm <- rbindlist(lapply(gene_id_GTEx, function(g_id, testis_GTEx_tpm, appris_df_tpm) {
    print(g_id)
    #print(g_id)
    #browser()
    ## Get all of the testis transcripts in the GTEx TPM file 
    GTEx_gene_df <- testis_GTEx_tpm[grep(g_id, testis_GTEx_tpm$gene_id)]
    GTEx_transcripts <- GTEx_gene_df$transcript_id %>% gsub("\\..*", "", .)
    
    ## Get the appris gene data  
    appris_df_tpm_gene <- appris_df_tpm[grep(g_id, appris_df_tpm$gene_id)]
    appris_transcripts <- appris_df_tpm_gene$transcript_id
    
    ## Get the appris only transcripts
    appris_only_transcripts <- appris_transcripts[!(appris_transcripts %in% GTEx_transcripts)]
    
    
    ####################################################################################
    ##### SANITY: Check if the appris-only transcript has a GTEx binary value of 0 #####
    ####################################################################################
    if (length(appris_only_transcripts) > 0) {
      GTEx_val <- appris_df_tpm[transcript_id %in% appris_only_transcripts]$GTEx
      if (any(GTEx_val != 0)) {
        warning_message <- paste0("One of the appris only transcripts: ", paste0(c(appris_only_transcripts), collapse = ", "), 
                                  "appears to also be found in GTEx, but a mean TPM score could not be calculated...")
        print(warning_message)
      } 
    }
    
    ## Get the GTEx only transcripts 
    GTEx_only_transcripts <- GTEx_transcripts[!(GTEx_transcripts %in% appris_transcripts)]
    
    ## Get the GTEx TPM values for each GTEx-only transcript --> if a transcript is 
    ## greater than the minimum APPRIS-GTEx TPM value for a given gene:  
    #t_id <- "ENST00000336199"
    GTEx_only <- rbindlist(lapply(GTEx_only_transcripts, function(t_id, GTEx_gene_df, appris_df_tpm_gene){
      #print(t_id)
      #browser()
      ## Get the mean TPM for the GTEx-only transcript
      GTEx_only_transcript_tpm <- GTEx_gene_df[grep(t_id, GTEx_gene_df$transcript_id), c(3:ncol(GTEx_gene_df)), with=FALSE] %>% 
        as.numeric() %>% mean()
      
      ## Get the minimum APPRIS-GTEx transcript TPM mean value 
      ## Some transcripts are not expressed/quantified in the GTEx dataset and will have NA values 
      min_appris_GTEx <- appris_df_tpm_gene$avg_GTEx_TPM %>% na.omit() %>% as.numeric() %>% min()
      if (length(min_appris_GTEx) == 0 ) {
        min_appris_GTEx <- 1
      }
      
      ##### SANITY: Check if the GTEx-only transcripts DO NOT have a TPM value greater than the thresholds
      if (!(GTEx_only_transcript_tpm >= min_appris_GTEx & GTEx_only_transcript_tpm >= 1)) {
        warning_message <- paste0("The GTEx-only-transcript, ", t_id, 
                                  " did not meet the threshold avg TPM cutoffs...")
        print(warning_message)
      }
      
      ## Check to see if the mean GTEx-only transcript TPM is greater than the minimum of a APPRIS-GTEx TPM 
      if (GTEx_only_transcript_tpm >= min_appris_GTEx & GTEx_only_transcript_tpm >= 1) {
        
        return(data.table("gene_name"=appris_df_tpm_gene$gene_name[1],
                          "gene_id"=appris_df_tpm_gene$gene_id[1],
                          "transcript_id"=t_id,
                          "CCDS"="not_in_appris",
                          "status"="minor",
                          "avg_GTEx_TPM"=round(GTEx_only_transcript_tpm, digits = 2),
                          "appris"=0,
                          "GTEx"=1))
      }
      
    }, 
    GTEx_gene_df=GTEx_gene_df, appris_df_tpm_gene=appris_df_tpm_gene))
    
    ## Combine the appris-GTEx and GTEx-only data
    appris_df_tpm_gene$appris <- 1
    if (nrow(GTEx_gene_df) > 0) {
      final <- rbind(appris_df_tpm_gene, GTEx_only)
      return(final)
    } 
    if (nrow(GTEx_gene_df) == 0) {
      return(appris_df_tpm_gene)
    }
    
  }, 
  testis_GTEx_tpm=testis_GTEx_tpm, appris_df_tpm=appris_df_tpm))
  
  return(appris_gtex_transcript_tpm)

}