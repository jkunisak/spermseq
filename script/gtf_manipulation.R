#####################################################################
##### Obtain Ensembl GTF data from APPRIS principal transcripts #####
#####################################################################

## Load in packages
library(data.table)
library(dplyr)
library(pbmcapply)
library(biomaRt)

## Function to read and manipulate appris data
gtf_manipulation <- function(gtf, df) {
  #browser()

  ## Define the mart to use 
  #ensembl = useMart("ensembl", dataset = "hsapiens_gene_ensembl") ## Ensembl Genes 103
  ensembl <- useMart("ENSEMBL_MART_ENSEMBL", host = "http://nov2020.archive.ensembl.org", dataset = "hsapiens_gene_ensembl") ## Ensembl Genes 102
  
  # df = appris_gtex_transcript_tpm
  # df <- df[c(225:230),]
  
  ## Convert df to list 
  if (class(df)[1] != "list") {
    df_list <- df %>% transpose(.) %>% as.list(.) %>% unname(.)
  }
  r <- df_list[[169]]
  
  final <- rbindlist(pbmclapply(mc.cores = 8, df_list, function(r, gtf, ensembl) {
    ## Print the row
    print(r)
    
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
      
      transcript_biotype <- ((gsub(pattern = "\"; ", replacement = "|", x = gtf_r$info[1]) %>% 
                                gsub(pattern = " \"", replacement = ": ", x = .) %>% strsplit(., split = "|", fixed = TRUE))[[1]] %>%
                               .[grep(pattern = "transcript_biotype", x = .)] %>% strsplit(., split = "transcript_biotype: "))[[1]][2]
      
      ####################################################################
      ##### SANITY: Check to see if the transcript is protein coding #####
      ####################################################################
      
      ## Define the attributes (output of the query)
      attributes <- c("transcript_biotype", "ensembl_peptide_id")
      
      ## Define the filters (input to the query)
      filters <- "ensembl_transcript_id"
      
      ## Run the biomaRt query
      biotype_info <- getBM(attributes = attributes, 
                            filters = filters,
                            values = transcript_id,
                            mart = ensembl, useCache = FALSE)
      print(biotype_info)
      
      ## Check if there are no values in biomaRt
      if (nrow(biotype_info) == 0) {
        warning_message <- paste0("WARNING: ", gene_name, "'s transcript: ", transcript_id, 
                                  " was labelled as protein coding by the GTF file but is not in biomart... Changing the ensembl_peptide_id label to the gtf protein_id: ", protein_id,  
                                  " and changing the transcript_type value to what is in the GTF file: ", transcript_biotype)
        print(warning_message)
        
        biotype_info <- data.table()
        biotype_info$transcript_biotype <- transcript_biotype
        biotype_info$ensembl_peptide_id <- protein_id
      }
      
      ## Look at cases where the transcript biotype is NOT protein_coding
      if (biotype_info$transcript_biotype != "protein_coding") {
        warning_message <- paste0("WARNING: ", gene_name, "'s transcript: ", transcript_id, 
                                  " is labelled as functionally important in APPRIS and/or has high testis-specific GTEx average TPM values. ", 
                                  "However, the transcript biotype is not identified as protein_coding using biomaRt...")
        print(warning_message)
      }
      
      ## Look at cases where there is no value in the transcript_biotype field from biomart
      if (length(biotype_info$transcript_biotype) == 0 | is.na(biotype_info$transcript_biotype)) {
        warning_message <- paste0("WARNING: ", gene_name, "'s transcript: ", transcript_id, 
                                  " is labelled as functionally important in APPRIS and/or has high testis-specific GTEx average TPM values. ", 
                                  "However, the transcript is not found in biomaRt...")
        print(warning_message)
      }
      
      ## Look at cases where there is no value in the ensembl_protein_id field from biomart
      if (length(biotype_info$ensembl_peptide_id) == 0 | is.na(biotype_info$ensembl_peptide_id)) {
        warning_message <- paste0("WARNING: ", gene_name, "'s transcript: ", transcript_id, 
                                  " is labelled as functionally important in APPRIS and/or has high testis-specific GTEx average TPM values. ", 
                                  "However, the transcript does not have a protein_id in biomaRt...", 
                                  " Will use protein id from gtf file.")
        print(warning_message)
        biotype_info$ensembl_peptide_id <- protein_id
        
      }
      
      ###########################################################################
      ##### SANITY: Check if the protein_ids from the gtf and biomaRt match #####
      ###########################################################################
      if (protein_id != biotype_info$ensembl_peptide_id) {
        warning_message <- paste0("WARNING: ", gene_name, "'s transcript: ", transcript_id, 
                                  " has discrepancies in the protein_id according to the GTF file and biomaRt... BiomaRt = ", 
                                  biotype_info$ensembl_peptide_id, " while GTF = ", protein_id, 
                                  "... Will use the GTF's protein_id value.")
        print(warning_message)
      }
      
      
      
      ## Get the necessary information for the transcript
      transcript_df <- gtf_r[, c("chr", "start", "stop")] %>% 
        mutate(transcript_id = transcript_id, gene_name = gene_name, gene_id = gene_id, 
               protein_id=protein_id, transcript_biotype=transcript_biotype,
               gtf=1)
      
      ## Get the exon lengths 
      transcript_df$exon_length <- transcript_df$stop - transcript_df$start + 1
      
      ##############################################################
      ##### SANITY: Check if the exon lengths are not negative #####
      ##############################################################
      
      if (any(transcript_df$exon_length <= 0)) {
        warning_message <- paste0("WARNING: ", gene_name, "'s transcript: ", transcript_id, 
                                  " has a negative (or = 0) exon length value of ", transcript_df$exon_length[which(transcript_df$exon_length <= 0)], 
                                  "... Please check if this is intended or not...")
        print(warning_message)
      }
      
      ## Get the total exon length
      transcript_df$total_exon_length <- sum(transcript_df$exon_length)
      
      ## Get the total number of exons for the transcript
      transcript_df$exon_num <- nrow(transcript_df)
    }
    
    ## Fill in data if transcript_df is empty (perhaps non-coding?); not in gtf file
    if (nrow(gtf_r) == 0) {

      ######################################################################################################################      
      ##### SANITY: Use biomaRt to get the transcript type (should NOT be protein coding if it has reached this point) #####
      ######################################################################################################################
      
      ## Define the attributes (output of the query)
      attributes <- c("transcript_biotype", "ensembl_peptide_id")
      
      ## Define the filters (input to the query)
      filters <- "ensembl_transcript_id"
      
      ## Run the biomaRt query
      biotype_info <- getBM(attributes = attributes, 
                            filters = filters,
                            values = transcript_id,
                            mart = ensembl, useCache = FALSE)
      print(biotype_info)
      
      ## Check if there are no values in biomaRt
      if (nrow(biotype_info) == 0) {
        warning_message <- paste0("WARNING: ", gene_name, "'s transcript: ", transcript_id, 
                                  " not in biomart or the GTF file... Changing the transcript_biotype and ensembl_peptide_id label to not_in_gtf_and_biomart")
        print(warning_message)
        
        biotype_info <- data.table()
        biotype_info$transcript_biotype <- "not_in_gtf_and_biomart"
        biotype_info$ensembl_peptide_id <- "not_in_gtf_and_biomart"
      }
      
      ## Change NA values
      if (length(biotype_info$transcript_biotype) == 0 | is.na(biotype_info$ensembl_peptide_id)) {
        print(paste0("WARNING: ", gene_name, "'s transcript: ", transcript_id, 
                     " did not have an ensembl_peptide_id in biomaRt (returned NA) and the GTF file. Changing the ensembl_peptide_id label to not_in_gtf_and_biomart"))
        biotype_info$ensembl_peptide_id <- "not_in_gtf_and_biomart"
      }
      
      ## See if the biotype_info$ensembl_peptide_id is not in biomaRt (e.g. ENST00000442234, a transcript of LIG4) 
      if (length(biotype_info$transcript_biotype) == 0 | is.na(biotype_info$transcript_biotype)) {
        print(paste0("WARNING: ", gene_name, "'s transcript: ", transcript_id, 
                     " was not found in the GTF file or biomaRt (returned empty value or NA when querying for transcript type). ",
                     "Changing the biotype label (i.e. transcript_biotype) to not_in_gtf_and_biomart"))
        biotype_info$transcript_biotype <- "not_in_gtf_and_biomart"
        #biotype_info$ensembl_peptide_id <- "not_in_gtf_and_biomart"
      }
      
      ## See if the biotype_info$ensembl_peptide_id is protein_coding 
      if (biotype_info$transcript_biotype == "protein_coding") {
        print(paste0("WARNING: ", gene_name, "'s transcript: ", transcript_id, 
                     " was labelled as protein coding by biomaRt but is not in the GTF file. Changing the ensembl_peptide_id label to protein_coding_not_in_gtf"))
        #biotype_info$transcript_biotype <- "protein_coding_not_in_gtf"
        biotype_info$ensembl_peptide_id <- "protein_coding_not_in_gtf"
      }
      
      
      transcript_df <- data.table("chr"=0, "start"=0, "stop"=0,
                                  "transcript_id"=transcript_id,
                                  "gene_name"=gene_name,
                                  "gene_id"=gene_id,
                                  "protein_id"=biotype_info$ensembl_peptide_id,
                                  "transcript_biotype"=biotype_info$transcript_biotype,
                                  "gtf"=0,
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

  }, gtf=gtf, ensembl=ensembl))
  
  return(final)
}