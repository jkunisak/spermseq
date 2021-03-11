########################################################################
##### Function to merge exons that belong to different transcripts #####
########################################################################
library(data.table)
library(dplyr)

## Function to change the pfam_ids in overlapping exons across transcripts of the same gene
merge_pfam <- function(tmp_pfam, exon_library_df, row) {
  ## Get rid of NA values
  tmp_pfam <- na.omit(tmp_pfam)
  if (length(tmp_pfam) > 0) {
    ## Iterate through all pfam_ids in the non-max transcript 
    for (i in 1:length(tmp_pfam)) {
      x <- tmp_pfam[i]
      
      ## Get the pfam_ids in the exon library
      exon_library_pfam <- (exon_library_df[row]$pfam_id %>% strsplit(x = ., split = "|", fixed = TRUE))[[1]]
      
      ## See if the pfam_id is NOT in the exon library
      if (!(x %in% exon_library_pfam)) {
        merged_pfam_ids <- paste0(c(exon_library_pfam, x), collapse = "|")
        exon_library_df[row]$pfam_id <- merged_pfam_ids
        exon_library_df[row]$pfam_exons <- exon_library_df[row]$pfam_exons + 1
      }
    }
    return(exon_library_df) 
  }
  if (length(tmp_pfam) == 0) {
    return(exon_library_df)
  }
}


#df <- gtf_df
merge_exons <- function(df) {
 
  ## Get the gene_ids
  gene_ids <- df[, c("gene_id", "gene_name"), with=FALSE] %>% unique() %>% transpose(.) %>% as.list(.) %>% unname(.)
  
  ## Merge overlapping exons 
  x <- gene_ids[[7]]
  merged_exons_df <- rbindlist(lapply(gene_ids, function(x, df) {
    gene_df <- df[gene_id == x[1]]
    #print(x)
    
    ## Get the transcript (of the gene_id) with the most number of exons 
    max_exon_transcript <- (max(gene_df$exon_num) %>% grep(., gene_df$exon_num) %>% 
      gene_df[., "transcript_id", with=FALSE] %>% unique())[[1]][1]
    
    ## Get the library with exons of the transcript with the most exons 
    cols_to_keep <- c("chr", "start", "stop", "gene_name", "gene_id", "transcript_id", "pfam_exons", "pfam_id")
    exon_library <- gene_df[which(gene_df$transcript_id == max_exon_transcript), colnames(gene_df) %in% cols_to_keep, with=FALSE]
    exon_library$info <- paste0(exon_library$chr, ":", exon_library$start, "-", exon_library$stop, " - o")
    exon_library$transcript_num <- 1

    ## Get the transcript(s) without the most number of exons
    non_max_exon_transcript <- gene_df$transcript_id[!(gene_df$transcript_id %in% max_exon_transcript)] %>% unique()
    
    ## If there is only 1 transcript for the gene_id --> return the exon_library as it is: 
    if (length(non_max_exon_transcript) == 0) {
      transcript_info <- data.table("total_exons"=nrow(exon_library),
                                    "exon_size"=(exon_library$stop - exon_library$start + 1))
      transcript_info$avg_exon_size <- round(mean(transcript_info$exon_size), digits = 2)
      transcript_info$total_nucleotides <- round(sum(transcript_info$exon_size), digits = 2)
      
      return(cbind(exon_library, transcript_info))
    }
    
    ## Go through each transcript that does not have the most number of exonns
    for (t in 1:length(non_max_exon_transcript)) {
      
      ## Get the gtf data for the non-max exon transcripts
      transcript <- non_max_exon_transcript[t]
      tmp <- gene_df[transcript_id==transcript, colnames(gene_df) %in% cols_to_keep, with=FALSE]
      
      ## Go through each row and see if the exon start/stop is already in the library
      for (r in 1:nrow(tmp)) {
        ## Get the start/stop/other info of the exon
        tmp_start <- tmp$start[r]
        tmp_stop <- tmp$stop[r]
        tmp_pfam <- (tmp$pfam_id[r] %>% strsplit(x = ., split = "|", fixed = TRUE) %>% unique())[[1]]
        tmp_info <- paste0(tmp$chr[r], ":", tmp_start, "-", tmp_stop)

        ## Check if the exon is in the exon library 
        in_library <- grep(tmp_info, exon_library$info)
        
        ## If the exon in the tmp dataset is in the exon library: 
        if (length(in_library) > 0) {
          exon_library <- merge_pfam(tmp_pfam = tmp_pfam, exon_library_df = exon_library, row = in_library) ## Change the pfam_id
          exon_library[in_library]$transcript_id <- paste0(exon_library[in_library]$transcript_id, " | ", transcript) ## Change the transcript_id
          exon_library[in_library]$transcript_num <- exon_library[in_library]$transcript_num + 1 ## Change the transcript_num
          exon_library[in_library]$info <- paste0(exon_library[in_library]$info, " | ", tmp_info, " - o") ## Change the transcript_info
        }
        
        ## If the exon in the tmp dataset is not in the exon library: 
        if (length(in_library) == 0) {
          
          ## Possibilities: 
          ## 0) exon spans multiple exons existing in the exon_library
          ## 1) start and stop encompass exon(s) in the existing exon_library
          ## 2) start is before exon(s) in the existing exon_library --> stop is the same/before
          ## 3) stop is after exons(s) in the existing exon_library --> start is the same/after
          ## 4) start and stop and wtihin an existing exon in the exon_library
          ## 5) start and/or stop are not within the bounds of the exons in the exon_library --> add it as a new exon in the library
          
          ## Function to change the exon positions if necessary 
          source("~/Google Drive/Quinlan Lab - PhD/Projects/spermseq/script/alter_exons.R")
          
          ## Situation (0): 
          sit0 <- exon_library[tmp_start <= start & tmp_stop >= stop]
          #sit0 <- exon_library[c(5,15,26)]
          if (nrow(sit0) > 1) {
            sit0_chr <- sit0$chr[1] ## Get the chr of the encompassing exon
            sit0_start <- sit0$start[1] ## Get the start of the encompassing exon
            sit0_stop <- sit0$stop[nrow(sit0)] ## Get the stop of the encompassing exon 
            sit0_transcript_temp <- unique(sit0$transcript_id) %>% strsplit(., split =" | ", fixed = TRUE) %>% unlist() ## Get the unique transcript values 
            sit0_transcript_id <- paste0(sit0_transcript_temp, " | ", transcript)
            sit0_transcript_num <- length(sit0_transcript_temp) + 1
            ## Get the unique pfam_ids
            if (length(na.omit(tmp_pfam)) > 0) {
              sit0_pfam_id_temp <- unique(sit0$pfam_id) %>%  strsplit(., split = "|", fixed = "TRUE") %>% unlist() %>% unique() %>% na.omit
              sit0_pfam_id <- unique(c(sit0_pfam_id_temp, tmp_pfam)) %>% na.omit() %>% paste(., collapse = "|")
              sit0_pfam_exons <- sum(sit0$pfam_exons) + 1
            }
            if (length(na.omit(tmp_pfam)) == 0) {
              sit0_pfam_id <- unique(sit0$pfam_id) %>%  strsplit(., split = "|", fixed = "TRUE") %>% unlist() %>% unique() %>% na.omit %>% paste(., collapse = "|")
              sit0_pfam_exons <- sum(sit0$pfam_exons)
            }

            
            ## Remove the rows that are associated with situation 0
            sit0_info <- sit0$info
            exon_library <- exon_library[!which(exon_library$info %in% sit0_info)]
            
            ## Generate exon-flattened data in situtation 0
            sit0_final <- data.table("chr"=sit0_chr, "start"=sit0_start, "stop"=sit0_stop,
                                     "transcript_id"=sit0_transcript_id,
                                     "gene_name"=x[2],
                                     "gene_id"=x[1],
                                     "pfam_exons"= sit0_pfam_exons,
                                     "pfam_id"=sit0_pfam_id,
                                     "info"=paste0(sit0_chr, ":", sit0_start, "-", sit0_stop, " - s0"),
                                     "transcript_num"=sit0_transcript_num)
            exon_library <- rbind(exon_library, sit0_final)
          }
        
          
          ## Situation (1): 
          sit1 <- exon_library[tmp_start < start & tmp_stop > stop]$info
          #sit1 <- exon_library$info[15]
          if (length(sit1) > 0 & nrow(sit0) == 0) {
            #print(paste0("length: ", length(sit1)))
            exon_library <- alter_exon(exon_library = exon_library, situation = "situation 1", tmp_info = tmp_info, 
                                       tmp_start = tmp_start, tmp_stop = tmp_stop, tmp_pfam = tmp_pfam, 
                                       sit_row = grep(sit1, exon_library$info), transcript = transcript)}
      
          ## Situation (2):
          sit2 <- exon_library[tmp_start < start & tmp_stop >= start & tmp_stop <= stop]$info
          if (length(sit2) > 0 & nrow(sit0) == 0) {
            #print(paste0("length: ", length(sit2)))
            exon_library <- alter_exon(exon_library = exon_library, situation = "situation 2", tmp_info = tmp_info, 
                                       tmp_start = tmp_start, tmp_stop = tmp_stop, tmp_pfam = tmp_pfam,
                                       sit_row = grep(sit2, exon_library$info), transcript = transcript)}
          
          ## Situation (3): 
          sit3 <- exon_library[tmp_start >= start & tmp_start <= stop & tmp_stop > stop]$info
          if (length(sit3) > 0 & nrow(sit0) == 0) {
            #print(paste0("length: ", length(sit3)))
            exon_library <- alter_exon(exon_library = exon_library, situation = "situation 3", tmp_info = tmp_info, 
                                       tmp_start = tmp_start, tmp_stop = tmp_stop, tmp_pfam = tmp_pfam,
                                       sit_row = grep(sit3, exon_library$info),
                                       transcript = transcript)}

          ## Situation (4): 
          sit4 <- exon_library[tmp_start >= start & tmp_stop <= stop]$info
          if (length(sit4) > 0 & nrow(sit0) == 0) {
            #print(paste0("t: ", t))
            #print(paste0("r: ", r))
            #print(paste0("length: ", length(sit4)))
            exon_library <- alter_exon(exon_library = exon_library, situation = "situation 4", tmp_info = tmp_info,
                                       tmp_start = tmp_start, tmp_stop = tmp_stop, tmp_pfam = tmp_pfam, 
                                       sit_row = grep(sit4, exon_library$info),
                                       transcript = transcript)}

          ## Situation (5): 
          if (length(c(sit1, sit2, sit3, sit4)) == 0 & nrow(sit0) == 0) {
            sit5_df <- data.table("chr"=tmp$chr[r],
                                  "start"=tmp_start,
                                  "stop"=tmp_stop, 
                                  "transcript_id"=transcript, 
                                  "gene_name"=x[2],
                                  "gene_id"=x[1],
                                  "pfam_exons"=tmp$pfam_exons[r],
                                  "pfam_id"=tmp_pfam,
                                  "info"=tmp_info,
                                  "transcript_num"=1)
            exon_library <- rbind(exon_library, sit5_df)
          }
          
          
        }
        
      }
    }
    
    ## Store important information for the transcripts
    transcript_info <- data.table("total_exons"=nrow(exon_library),
                                  "exon_size"=exon_library$stop - exon_library$start + 1)
    transcript_info$avg_exon_size <- round(mean(transcript_info$exon_size), digits = 2)
    transcript_info$total_nucleotides <- round(sum(transcript_info$exon_size), digits = 2)
    
    final <- cbind(exon_library, transcript_info)
    #print(paste0(x[2], ": ", x[1], " has ncol: ", ncol(final)))
    return(final)
    
  }, df=df), fill = TRUE)
  
  return(merged_exons_df)
}