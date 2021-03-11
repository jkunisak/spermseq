#####################################################
##### Change the positions of overlapping exons #####
#####################################################
library(data.table)
library(dplyr)

alter_exon <- function(exon_library, situation, tmp_info, tmp_start, tmp_stop, tmp_pfam, sit_row, transcript) {
  ## Change the start/stop position if needed
  if (situation == "situation 1") {
    exon_library[sit_row]$start <- tmp_start
    exon_library[sit_row]$stop <- tmp_stop
    tmp_info <- paste0(tmp_info, " - s1")
  }
  if (situation == "situation 2") {
    exon_library[sit_row]$start <- tmp_start
    tmp_info <- paste0(tmp_info, " - s2")
  }
  if (situation == "situation 3") {
    exon_library[sit_row]$stop <- tmp_stop
    tmp_info <- paste0(tmp_info, " - s3")
  }
  if (situation == "situation 4") {
    tmp_info <- paste0(tmp_info, " - s4")
  }
  
  ## Add the transcript_ids that have each exon
  exon_library[sit_row]$transcript_id <- paste0(exon_library[sit_row]$transcript_id, " | ", transcript)
  
  ## Add how the situation of the transcript's exon
  exon_library[sit_row]$info <- paste0(exon_library[sit_row]$info, " | ", tmp_info)
  
  ## Add 1 to the number of transcripts with the exon
  exon_library[sit_row]$transcript_num <- exon_library[sit_row]$transcript_num + 1
  
  ## Merge pfam domains across transcripts
  if (length(na.omit(tmp_pfam)) > 0) {
    exon_library_pfam <- unique(exon_library[sit_row]$pfam_id) %>% strsplit(., split = "|", fixed = TRUE) %>% unlist() %>% na.omit
    exon_library[sit_row]$pfam_id <- unique(c(exon_library_pfam, tmp_pfam)) %>% na.omit() %>% paste(., collapse = "|")
    
    ## Add 1 to the number of times an exon contributes to a different pfam domain 
    exon_library[sit_row]$pfam_exons <- exon_library[sit_row]$pfam_exons + 1
  }
  if (length(na.omit(tmp_pfam)) == 0) {
    return(exon_library)
  }
  
  return(exon_library)

}