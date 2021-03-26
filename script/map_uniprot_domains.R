##################################################################################################
##### Function to map protein coordinates of uniprot function domains to genomic coordinates #####
##################################################################################################
## Load in necessary packages
library(data.table)
library(dplyr)
library(ensembldb)
library(EnsDb.Hsapiens.v86)
library(AnnotationHub)
library(biomaRt)
options(warn=1)

# ## Read in uniprot domian data for transcripts that did not have functional domains from biomaRt
# uniprot <- fread("~/git/spermseq/twinstrand_target_gene_set/output/add_uniprot_domains.txt", header = TRUE)

## Test data
# row <- uniprot[16]
# gene <- row$gene_name
# start <- row$p_start
# end <- row$p_end
# transcript <- row$transcript_id
# domain <- paste0(uniprot$uniprot_id, ":", uniprot$domain_name)

## Write the function
map_uniprot_domains <- function(uniprot, gtf_df, ensembl) {
  
  ah <- AnnotationHub()
  query(ah, "EnsDb.Hsapiens.v101")
  edb <- ah[["AH83216"]]
  
  transcripts <- uniprot$transcript_id %>% unique
  
  df <- gtf_df
  
  ## Iterate through each transcript
  #transcript <- transcripts[1]
  uniprot_mapping_final <- rbindlist(lapply(transcripts, function(transcript, uniprot, df, ensembl, edb){
   # browser()
    ## Get the uniprot transcript information
    uniprot_transcript_df <- uniprot[transcript_id == transcript]
    
    ############################################################################################
    ##### SANITY: Check if the exon positions used in the edb are the same as the GTF file #####
    ############################################################################################
    ## Get the coding sequence exons for the gene of interest (x[2])
    cds <- cdsBy(edb, by = "tx", filter = ~ tx_id == transcript)
    
    #############################################################
    #### Check 1: see if there are the same number of exons #####
    #############################################################
    cds_exons <- cds[[1]]
    gtf_exons <- df[transcript_id == transcript]
    cds_exon_num <- length(cds_exons@elementMetadata$exon_rank)
    gtf_exon_num <- nrow(gtf_exons)
    if (cds_exon_num != gtf_exon_num) {
      warning_message <- paste0("The number of exons in gene's", gene, " transcript: ", transcript, " does not match in the edb used (v101): ", cds_exon_num, 
                                " and the gtf file (v102): ", gtf_exon_num, "...PLEASE LOOK INTO THIS...")
      print(warning_message)
      stop()
    }
    
    ###################################################
    ## Check 2: see if the positions of the exons match 
    ###################################################
    exon_num <- 1
    for (exon_num in 1:cds_exon_num) {
      ## Get the start/stop coordinates in the edb and gtf file
      cds_start_stop <- paste0(cds_exons@ranges@start[exon_num], "-", (cds_exons@ranges@start[exon_num] + cds_exons@ranges@width[exon_num] -1))
      gtf_start_stop <- paste0(gtf_exons$start[exon_num], "-", gtf_exons$stop[exon_num])
      
      ## Perform the comparison
      if (cds_start_stop != gtf_start_stop) {
        warning_message <- paste0("The start/stop coordinates of exon ", exon_num, " of gene's ", unique(gtf_exons$gene_name), " transcript: ", transcript, " are not equivalent in the ",
                                  "edb dataset (v101):", cds_start_stop, " and the gtf file (v102):", gtf_start_stop, "... PLEASE LOOK INTO THIS")
        print(warning_message)
      }
    }
    
    ## Merge intervals of uniprot genomic coordinates and gtf exon coordinates
    gtf_exons$pfam_exons <- 0
    gtf_exons$pfam_id <- NA
    
    ## Iterate through each uniprot domain
    domain_num <- 1
    for (domain_num in 1:length(uniprot_transcript_df$uniprot_id)) {
      ## Get the start/end positions of the uniprot domains 
      start <- uniprot_transcript_df$p_start[domain_num]
      end <- uniprot_transcript_df$p_end[domain_num]
      
      ## Get the gene of interest 
      gene <- uniprot_transcript_df$gene_name[domain_num]
      
      ## Get the domain 
      domain <- paste0(uniprot_transcript_df$uniprot_id[domain_num], ":", uniprot_transcript_df$domain_name[domain_num])
      
      #############################################################################
      ##### Use biomaRt to identify protein id for genomic coordinate mapping #####
      #############################################################################
      ## Get the protein_id for the transcript of interest 
      ## Define the attributes (output)
      attributes <- c("ensembl_peptide_id")
      
      ## Define the filters (input to the query)
      filters <- "ensembl_transcript_id"
      
      ## Use biomaRt to derive the protein_id 
      info <- getBM(attributes = attributes, filters = filters, values = transcript, 
                    mart = ensembl, useCache = FALSE)
      
      ## Get the genomic coordinates of the 
      uniprot_p_coord <- IRanges(start = start, end = end, names = info$ensembl_peptide_id)
      uniprot_g_coord <- proteinToGenome(uniprot_p_coord, edb)
      
      ## If the cds could not be found --> return the entire transcript
      if (any(uniprot_g_coord[[1]]@elementMetadata$cds_ok) == FALSE) {
        gtf_exons$pfam_id <- as.character(gtf_exons$pfam_id)
        gtf_exons$pfam_id <- as.character("including entire cds for this transcript")
        return(gtf_exons)
      }
      
      ####################################################################################
      ##### SANITY: Check if the transcript actually translates to a protein product #####
      ####################################################################################
      if (length(uniprot_g_coord[[1]]@ranges) == 0) {
        warning_message <- paste0("The transcript: ", transcript, " expressed from ", gene, " does not have a coding sequence in the ensembldb... skipping...")
        print(warning_message)
      }
      
      #####################################################################################
      ##### SANITY: Check to see if multiple protein_ids are mapped to the transcript #####
      #####################################################################################
      if (length(uniprot_g_coord) > 1) {
        warning_message <- paste0("Transcript information (below) mapped to ", length(uniprot_g_coord), " genes...", 
                                  paste0(info$ensembl_peptide_id, collapse = " | "), "... FIX THIS...")
        print(row)
        print(warning_message)
      }
      
      #################################################################################
      ##### Determine which gtf exons overlap with the uniprot domain of interest #####
      #################################################################################
      ## Get uniprot genomic coordinates 
      uniprot_start <- uniprot_g_coord[[1]]@ranges@start
      uniprot_end <- uniprot_g_coord[[1]]@ranges@start + uniprot_g_coord[[1]]@ranges@width - 1
      
      exon_num <- 1
      if (length(uniprot_g_coord[[1]]@ranges) > 0) {
        for (exon_num in 1:length(uniprot_g_coord[[1]]@ranges)) {
          start <- uniprot_start[exon_num]
          end <- uniprot_end[exon_num]
          print(domain)
          
          ## See if either the pfam genomic start/stop coordinate is within an exon 
          gstart_exon <- which((start >= gtf_exons$start) & (start <= gtf_exons$stop))
          gstop_exon <- which((end >= gtf_exons$start) & (end <= gtf_exons$stop)) 
          
          ## Check if the pfam domain is not found in any exon 
          if (length(gstart_exon) == 0 & length(gstop_exon) == 0) {
            message(paste0(domain, " is not found in any exons of ", transcript, " (gene: ", gene, ")"))
          }
          
          ## Check if multiple exons correspond to a single pfam domain --> shouldn't because the pfam output describes the genomic coordinate of each exon 
          if (length(gstart_exon) > 1 | length(gstop_exon) > 1) {
            message(paste0("Multiple exons correspond to ", domain, " in ", transcript, " (gene: ", gene, ")"))
            stop()
          }
          
          ## Where the uniprot domain start and stop coordinates are within an exon 
          if (length(gstart_exon) == 1 & length(gstop_exon) == 1) {
            if (gstart_exon == gstop_exon) {
              ## See if the pfam domain column is na --> if so, put in pfam domain, if not, add "|" 
              ## Also chceck to make sure that the pfam domain is not already annotated in the transcript at a given exon
              
              ## See above: the pfam_id column in the gtf_exons column is already initiated to NAs
              if (is.na(gtf_exons$pfam_id[gstart_exon])) {
                gtf_exons$pfam_id[gstart_exon] <- domain
              }
              if (!is.na(gtf_exons$pfam_id[gstart_exon]) & length(grep(domain, x = gtf_exons$pfam_id[gstart_exon])) == 0) {
                gtf_exons$pfam_id[gstart_exon] <- paste0(gtf_exons$pfam_id[gstart_exon], "|", domain)
              }
              gtf_exons$pfam_exons[gstart_exon] <- gtf_exons$pfam_exons[gstart_exon] + 1
            }
          }
          
          ## Where only uniprot domain start is within an exon
          if ((length(gstart_exon) == 1 & length(gstop_exon) != 1)) {
            message(paste0("Only start coordinate of uniprot domain: ", domain, " is within an exon"))
            ## See if the pfam domain column is na --> if so, put in pfam domain, if not, add "|" 
            if (is.na(gtf_exons$pfam_id[gstart_exon])) {
              gtf_exons$pfam_id[gstart_exon] <- domain
            }
            ## Also check to make sure that the uniprot domain is not already annotated in the transcript at a given exon
            if (!is.na(gtf_exons$pfam_id[gstart_exon]) & length(grep(domain, x = gtf_exons$pfam_id[gstart_exon])) == 0) {
              gtf_exons$pfam_id[gstart_exon] <- paste0(gtf_exons$pfam_id[gstart_exon], "|", domain)
            }
            gtf_exons$pfam_exons[gstart_exon] <- gtf_exons$pfam_exons[gstart_exon] + 1
          }
        }
        
        ## Where only the uniprot domain stop coordinate is within an exon 
        if (length(gstart_exon) != 1 & length(gstop_exon) == 1) {
          message(paste0("Only stop coordinate of uniprot domain: ", domain, " is within an exon"))
          gtf_exons$pfam_exons[gstop_exon] <- gtf_exons$pfam_exons[gstop_exon] + 1
          
          ## See if the pfam domain column is na --> if so, put in pfam domain, if not, add "|" 
          if (is.na(gtf_exons$pfam_id[gstop_exon])) {
            gtf_exons$pfam_id[gstop_exon] <- domain
          }
          if (!is.na(gtf_exons$pfam_id[gstop_exon]) & length(grep(domain, x = gtf_exons$pfam_id[gstop_exon])) == 0) {
            gtf_exons$pfam_id[gstop_exon] <- paste0(gtf_exons$pfam_id[gstop_exon], "|", domain)
          }
          gtf_exons$pfam_exons[gstop_exon] <- gtf_exons$pfam_exons[gstop_exon] + 1
        } 
      }
    }
    return(gtf_exons)
  }, 
  uniprot=uniprot, df=df, ensembl=ensembl, edb=edb))
  
  return(uniprot_mapping_final)
}