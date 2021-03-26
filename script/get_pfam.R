########################################################################
##### Function to merge exons that belong to different transcripts #####
########################################################################
library(data.table)
library(dplyr)
library(ensembldb)
library(EnsDb.Hsapiens.v86)
library(AnnotationHub)
library(biomaRt)
options(warn = 1)

#df <- gtf_df
get_pfam <- function(df, ensembl, edb) {
  ###########################################################################
  ##### SANITY: Initiate file to store transcripts without pfam domains #####
  ###########################################################################
  cat("", file = "~/git/spermseq/twinstrand_target_gene_set/output/gene_transcripts_no_pfam.txt", append = FALSE)
  cat("", file = "~/git/spermseq/twinstrand_target_gene_set/output/gene_transcripts_no_CDS_of_expected_length.txt", append = FALSE)
    
  ## Select for genes/transcripts that are found in the gtf, appris, and GTEx files
  df <- df[gtf == 1 & transcript_biotype == "protein_coding" & (GTEx == 1 | appris == 1)]
  ## df[!(gtf == 1 & transcript_biotype == "protein_coding" & (GTEx == 1 | appris == 1))] ## Sanity check to see what was left out
   
  
  ## Get the transcript_ids
  transcript_ids <- df[, c("transcript_id", "gene_name"), with=FALSE] %>% unique() %>% transpose(.) %>% as.list(.) %>% unname(.)
  
  ## Iterate through each transcript and get the pfam domain information for each exon
  x <- transcript_ids[[99]] ## FGFR2 
  x <- transcript_ids[[201]] ## USP9Y
  x <- transcript_ids[[2]]
  exon_pfam_df <- rbindlist(lapply(transcript_ids, function(x, df, ensembl, edb) {
    print(x)
    transcript_df <- df[transcript_id == x[1]]
    
    ## Initiate necessary pfam columns 
    transcript_df <- cbind(transcript_df, data.table("pfam_exons"=0,
                                         "pfam_id"=NA))
    
    ## Iterate through each transcript to get: 
    # 1) start/stop of each exon  
    # 2) pfam domains of the transcript 
    
    #####################################################
    ##### Get the functional domains with ensembldb #####
    #####################################################
    # pfam <- proteins(edb, filter = ~ tx_id %in% x[1] & protein_domain_source == "pfam",
    #                  columns = c("protein_domain_id", "prot_dom_start", "prot_dom_end"))
    
    ############################################################################
    ##### Use biomaRt to get the Pfam domain IDs/start/stop and protein id #####
    ############################################################################
    
    ## Define the attributes (output of the query)
    attributes <- c("pfam", "pfam_start", "pfam_end", "ensembl_peptide_id")

    ## Define the filters (input to the query)
    filters <- "ensembl_transcript_id"

    ## Run the biomaRt query
    pfam <- getBM(attributes = attributes,
                          filters = filters,
                          values = x[1],
                          mart = ensembl, useCache = FALSE)
    
    ###############################################################################
    ##### SANITY: Check if there are no pfam domains for the given transcript #####
    ###############################################################################
    if (nrow(pfam) == 0 | any(is.na(pfam))) {
      print(paste0("gene: ", x[2], " with transcript: ", x[1], " do not have any pfam information"))

      dummy <- paste0(x, collapse = "\t")
      write(dummy, file = "~/git/spermseq/twinstrand_target_gene_set/output/gene_transcripts_no_pfam.txt", append=TRUE)
      return(transcript_df)
    }
    
    ## Convert pfam domain start/end to genomic coordinates (biomaRt)
    pfam_pstart <- pfam$pfam_start
    pfam_pend <- pfam$pfam_end
    pfam_pid <- pfam$ensembl_peptide_id
    pfam_id <- pfam$pfam
    
    ## Convert pfam domain start/end to genomic coordinates (ensembldb)
    # pfam_pstart <- pfam@listData$prot_dom_start
    # pfam_pend <- pfam@listData$prot_dom_end
    # pfam_pid <- pfam@listData$protein_id
    # pfam_id <- pfam@listData$protein_domain_id
    pfam_genomic <- IRanges(start = pfam_pstart, end = pfam_pend, names = pfam_pid)
    pfam2g <- proteinToGenome(pfam_genomic, edb)
    
    ## Get the pfam genomic coordinates and pfam_id --> save as a data.table()
    ## Iterate through each domain (p) and each exon contributing to a given domain (w)
    pfam_df <- data.table()
    for (p in 1:length(pfam_id)) {
      for (w in 1:nrow(pfam2g[[p]]@elementMetadata)) {
        ## Get the boolean to see if the correct CDS was found
        cds_ok_bool <- pfam2g[[p]]@elementMetadata$cds_ok[w]
        
        ## If the correct CDS was found: make the dataset
        if (cds_ok_bool == TRUE) {
          pfamID <- pfam_id[p]
          pfam_gstart <- pfam2g[[p]]@ranges@start[w]
          pfam_gstop <- pfam2g[[p]]@ranges@start[w] + pfam2g[[p]]@ranges@width[w] - 1
          pfam_df <- rbind(pfam_df, data.table("pfam_id"=pfamID,
                                               "pfam_genomic_start"=pfam_gstart,
                                               "pfam_genomic_stop"=pfam_gstop)) 
        }
        
        ############################################################################################################
        ##### SANITY: Check if there was a 'Could not find a CDS whith the expected length for protein:' error #####
        ############################################################################################################
        if (cds_ok_bool == FALSE) {
          warning_message <- paste0("Could not find a CDS whith the expected length for protein: ", pfam_pid,
                                    " and its transcript: ", x[1], " at row: ", w, "... including entire cds for this transcript")
          print(warning_message)
          print(pfam_df)
          print(transcript_df)
          dummy <- paste0(x, collapse = "\t")
          write(dummy, "~/git/spermseq/twinstrand_target_gene_set/output/gene_transcripts_no_CDS_of_expected_length.txt", append = TRUE)
          
          transcript_df$pfam_id <- "including entire cds for this transcript"
          return(transcript_df)
        }
      }
    }
    
    ############################################################################################
    ##### SANITY: Check if the exon positions used in the edb are the same as the GTF file #####
    ############################################################################################
    ## Get the coding sequence exons for the gene of interest (x[2])
    cds <- cdsBy(edb, by = "tx", filter = ~ tx_id == x[1])
    
    #############################################################
    #### Check 1: see if there are the same number of exons #####
    cds_exons <- cds[[1]]
    gtf_exons <- df[transcript_id == x[1]]
    cds_exon_num <- length(cds_exons@elementMetadata$exon_rank)
    gtf_exon_num <- nrow(gtf_exons)
    if (cds_exon_num != gtf_exon_num) {
      warning_message <- paste0("The number of exons in gene's", x[2], " transcript: ", x[1], " does not match in the edb used (v101): ", cds_exon_num, 
                                " and the gtf file (v102): ", gtf_exon_num, "...PLEASE LOOK INTO THIS...")
      print(warning_message)
      stop()
    }
    ###################################################
    ## Check 2: see if the positions of the exons match 
    exon_num <- 1
    for (exon_num in 1:cds_exon_num) {
      ## Get the start/stop coordinates in the edb and gtf file
      cds_start_stop <- paste0(cds_exons@ranges@start[exon_num], "-", (cds_exons@ranges@start[exon_num] + cds_exons@ranges@width[exon_num] -1))
      gtf_start_stop <- paste0(gtf_exons$start[exon_num], "-", gtf_exons$stop[exon_num])
      
      ## Perform the comparison
      if (cds_start_stop != gtf_start_stop) {
        warning_message <- paste0("The start/stop coordinates of exon ", exon_num, " of gene's ", x[2], " transcript: ", x[1], " are not equivalent in the ",
                                  "edb dataset (v101):", cds_start_stop, " and the gtf file (v102):", gtf_start_stop, "... PLEASE LOOK INTO THIS")
        print(warning_message)
      }
    }
    
    
    ## Go through each pfam domain and identify which exons reside in those regions (at least partially)
    i=6
    for (i in 1:nrow(pfam_df)) {
      #print(i)
      ## Get information of individual pfam domain
      pfam_domain <- pfam_df$pfam_id[i]
      pfam_gstart <- pfam_df$pfam_genomic_start[i]
      pfam_gstop <- pfam_df$pfam_genomic_stop[i]
      
      ## See if either the pfam genomic start/stop coordinate is within an exon 
      pfam_gstart_exon <- which((pfam_gstart >= transcript_df$start) & (pfam_gstart <= transcript_df$stop))
      pfam_gstop_exon <- which((pfam_gstop >= transcript_df$start) & (pfam_gstop <= transcript_df$stop)) 
      
      ## Check if the pfam domain is not found in any exon 
      if (length(pfam_gstart_exon) == 0 & length(pfam_gstop_exon) == 0) {
        message(paste0(pfam_domain, " is not found in any exons of ", x[1], " (gene: ", x[2], ")"))
      }
      
      ## Check if multiple exons correspond to a single pfam domain --> shouldn't because the pfam output describes the genomic coordinate of each exon 
      if (length(pfam_gstart_exon) > 1 | length(pfam_gstop_exon) > 1) {
        message(paste0("Multiple exons correspond to ", pfam_domain, " in ", x[1], " (gene: ", x[2], ")"))
        stop()
      }
      ## Where the pfam domain start and stop coordinates are within an exon 
      if (length(pfam_gstart_exon) == 1 & length(pfam_gstop_exon) == 1) {
        if (pfam_gstart_exon == pfam_gstop_exon)
        ## See if the pfam domain column is na --> if so, put in pfam domain, if not, add "|" 
        ## Also chceck to make sure that the pfam domain is not already annotated in the transcript at a given exon
          
        ## See above: the pfam_id column in the transcript_df column is already initiated to NAs
        if (is.na(transcript_df$pfam_id[pfam_gstart_exon])) {
          transcript_df$pfam_id[pfam_gstart_exon] <- pfam_domain
        }
        if (!is.na(transcript_df$pfam_id[pfam_gstart_exon]) & length(grep(pfam_domain, x = transcript_df$pfam_id[pfam_gstart_exon])) == 0) {
          transcript_df$pfam_id[pfam_gstart_exon] <- paste0(transcript_df$pfam_id[pfam_gstart_exon], "|", pfam_domain)
        }
        transcript_df$pfam_exons[pfam_gstart_exon] <- transcript_df$pfam_exons[pfam_gstart_exon] + 1
      }
      
      if ((length(pfam_gstart_exon) == 1 & length(pfam_gstop_exon) != 1)) {
        message(paste0("Only start coordinate of pfam domain: ", pfam_domain, " is within an exon"))
        ## See if the pfam domain column is na --> if so, put in pfam domain, if not, add "|" 
        ## Also chceck to make sure that the pfam domain is not already annotated in the transcript at a given exon
        if (is.na(transcript_df$pfam_id[pfam_gstart_exon])) {
          transcript_df$pfam_id[pfam_gstart_exon] <- pfam_domain
        }
        if (!is.na(transcript_df$pfam_id[pfam_gstart_exon]) & length(grep(pfam_domain, x = transcript_df$pfam_id[pfam_gstart_exon])) == 0) {
          transcript_df$pfam_id[pfam_gstart_exon] <- paste0(transcript_df$pfam_id[pfam_gstart_exon], "|", pfam_domain)
        }
        transcript_df$pfam_exons[pfam_gstart_exon] <- transcript_df$pfam_exons[pfam_gstart_exon] + 1
      }
      
      ## Where only the pfam domain stop coordinate is within an exon 
      if (length(pfam_gstart_exon) != 1 & length(pfam_gstop_exon) == 1) {
        message(paste0("Only stop coordinate of pfam domain: ", pfam_domain, " is within an exon"))
        transcript_df$pfam_exons[pfam_gstop_exon] <- transcript_df$pfam_exons[pfam_gstop_exon] + 1
        
        ## See if the pfam domain column is na --> if so, put in pfam domain, if not, add "|" 
        if (is.na(transcript_df$pfam_id[pfam_gstop_exon])) {
          transcript_df$pfam_id[pfam_gstop_exon] <- pfam_domain
        }
        if (!is.na(transcript_df$pfam_id[pfam_gstop_exon]) & length(grep(pfam_domain, x = transcript_df$pfam_id[pfam_gstop_exon])) == 0) {
          transcript_df$pfam_id[pfam_gstop_exon] <- paste0(transcript_df$pfam_id[pfam_gstop_exon], "|", pfam_domain)
        }
        transcript_df$pfam_exons[pfam_gstop_exon] <- transcript_df$pfam_exons[pfam_gstop_exon] + 1
        #transcript_df$pfam_id[pfam_gstop_exon] <- paste0(transcript_df$pfam_id[pfam_gstop_exon], "|", pfam_domain)
      }
      # print(i)
      # print(pfam_gstart_exon)
      # print(pfam_gstop_exon)
      # print(transcript_df$pfam_exons)
    }
    
    return(transcript_df)
  }, df=df, ensembl=ensembl, edb=edb))
  return(exon_pfam_df)
}