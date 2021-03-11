########################################################################
##### Function to merge exons that belong to different transcripts #####
########################################################################
library(data.table)
library(dplyr)
library(ensembldb)
library(EnsDb.Hsapiens.v86)

#df <- gtf_df
get_pfam <- function(df) {
  ## Get rid of rows that do not have GTEx_TPM/gtf CDS values 
  df <- df[complete.cases(df),]
  df <- df[chr != 0 & chr != start]
  
  ## Get the transcript_ids
  transcript_ids <- df[, c("transcript_id", "gene_name"), with=FALSE] %>% unique() %>% transpose(.) %>% as.list(.) %>% unname(.)
  
  ## Iterate through each transcript and get the pfam domain information for each exon
  x <- transcript_ids[[30]]
  exon_pfam_df <- rbindlist(lapply(transcript_ids, function(x, df) {
    print(x)
    transcript_df <- df[transcript_id == x[1]]
    
    ## Initiate necessary pfam columns 
    transcript_df <- cbind(transcript_df, data.table("pfam_exons"=0,
                                         "pfam_id"=NA))
    
    ## Get the ensembl human GRCh38 annotation dataaset
    edb <- EnsDb.Hsapiens.v86
    
    ## Iterate through each transcript to get: 
    # 1) start/stop of each exon  
    # 2) pfam domains of the transcript 
    
    ## Get the functional domains with ensembldb
    pfam <- proteins(edb, filter = ~ tx_id %in% x[1] & protein_domain_source == "pfam",
                     columns = c("protein_domain_id", "prot_dom_start", "prot_dom_end"))
    if (nrow(pfam) == 0) {
      print(paste0("gene: ", x[2], " with transcript: ", x[1], " do not have any pfam information"))
      return(transcript_df)
    }
    
    ## Convert pfam domain start/end to genomic coordinates 
    pfam_pstart <- pfam@listData$prot_dom_start
    pfam_pend <- pfam@listData$prot_dom_end
    pfam_pid <- pfam@listData$protein_id
    pfam_id <- pfam@listData$protein_domain_id
    pfam_genomic <- IRanges(start = pfam_pstart, end = pfam_pend, names = pfam_pid)
    pfam2g <- proteinToGenome(pfam_genomic, edb)
    
    ## Get the pfam genomic coordinates and pfam_id --> save as a data.table()
    pfam_df <- data.table()
    for (p in 1:length(pfam_id)) {
      pfamID <- pfam_id[p]
      pfam_gstart <- pfam2g[[p]]@ranges@start
      pfam_gstop <- pfam2g[[p]]@ranges@start + pfam2g[[p]]@ranges@width - 1
      
      pfam_df <- rbind(pfam_df, data.table("pfam_id"=pfamID,
                                           "pfam_genomic_start"=pfam_gstart,
                                           "pfam_genomic_stop"=pfam_gstop))
      
    }
    
    ## Go through each pfam domain and identify which exons reside in those regions (at least partially)
    for (i in 1:nrow(pfam_df)) {
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
        if (is.na(transcript_df$pfam_id[pfam_gstart_exon])) {
          transcript_df$pfam_id[pfam_gstart_exon] <- pfam_domain
          transcript_df$pfam_exons[pfam_gstart_exon] <- transcript_df$pfam_exons[pfam_gstart_exon] + 1
        }
        if (!is.na(transcript_df$pfam_id[pfam_gstart_exon]) & length(grep(pfam_domain, x = transcript_df$pfam_id[pfam_gstart_exon])) == 0) {
          transcript_df$pfam_id[pfam_gstart_exon] <- paste0(transcript_df$pfam_id[pfam_gstart_exon], "|", pfam_domain)
          transcript_df$pfam_exons[pfam_gstart_exon] <- transcript_df$pfam_exons[pfam_gstart_exon] + 1
          
        }
      }
      
      if ((length(pfam_gstart_exon) == 1 & length(pfam_gstop_exon) != 1)) {
        message(paste0("Only start coordinate of pfam domain: ", pfam_domain, " is within an exon"))
        ## See if the pfam domain column is na --> if so, put in pfam domain, if not, add "|" 
        ## Also chceck to make sure that the pfam domain is not already annotated in the transcript at a given exon
        if (is.na(transcript_df$pfam_id[pfam_gstart_exon])) {
          transcript_df$pfam_id[pfam_gstart_exon] <- pfam_domain
          transcript_df$pfam_exons[pfam_gstart_exon] <- transcript_df$pfam_exons[pfam_gstart_exon] + 1
        }
        if (!is.na(transcript_df$pfam_id[pfam_gstart_exon]) & length(grep(pfam_domain, x = transcript_df$pfam_id[pfam_gstart_exon])) == 0) {
          transcript_df$pfam_id[pfam_gstart_exon] <- paste0(transcript_df$pfam_id[pfam_gstart_exon], "|", pfam_domain)
          transcript_df$pfam_exons[pfam_gstart_exon] <- transcript_df$pfam_exons[pfam_gstart_exon] + 1
          
        }
      }
      
      ## Where only the pfam domain stop coordinate is within an exon 
      if (length(pfam_gstart_exon) != 1 & length(pfam_gstop_exon) == 1) {
        message(paste0("Only stop coordinate of pfam domain: ", pfam_domain, " is within an exon"))
        transcript_df$pfam_exons[pfam_gstop_exon] <- transcript_df$pfam_exons[pfam_gstop_exon] + 1
        
        ## See if the pfam domain column is na --> if so, put in pfam domain, if not, add "|" 
        if (is.na(transcript_df$pfam_id[pfam_gstop_exon])) {
          transcript_df$pfam_id[pfam_gstop_exon] <- pfam_domain
        }
        if (!is.na(transcript_df$pfam_id[pfam_gstart_exon]) & length(grep(pfam_domain, x = transcript_df$pfam_id[pfam_gstop_exon])) == 0) {
          transcript_df$pfam_id[pfam_gstop_exon] <- paste0(transcript_df$pfam_id[pfam_gstop_exon], "|", pfam_domain)
        }
        transcript_df$pfam_id[pfam_gstop_exon] <- paste0(transcript_df$pfam_id[pfam_gstop_exon], "|", pfam_domain)
      }
      # print(i)
      # print(pfam_gstart_exon)
      # print(pfam_gstop_exon)
      # print(transcript_df$pfam_exons)
    }
    
    return(transcript_df)
  }, df=df))
  return(exon_pfam_df)
}