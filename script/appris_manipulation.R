################################################################
##### Obtain APPRIS transcript data from genes of interest #####
################################################################

## Load in packages
library(data.table)
library(dplyr)
library(biomaRt)

## Function to read and manipulate appris data
appris_manipulation <- function(appris_file, genes) {
  appris_colnames <- c("gene_name", "gene_id", "transcript_id", "CCDS", "status")
  appris <- fread(appris_file, header=FALSE) %>% `colnames<-` (appris_colnames)
  
  ## Get the appris principal isoform data for the genes of interest 
  appris_genes <- appris[gene_name %in% genes]
  
  #################################################################################
  ##### SANITY: Check See if there are any genes that do not have APPRIS data #####
  #################################################################################
  
  not_in_appris <- genes[! genes %in% unique(appris_genes$gene_name)]
  message(paste0("The following genes do not have any APPRIS isoform information: ", paste(c(not_in_appris), collapse = ", ")))
  
  ## IF there are genes not in APPRIS --> get the Gene ID
  if (length(not_in_appris) > 0){
    
    ## Get the ensembl human GRCh38 annotation dataaset
    ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")
    
    for (i in 1:length(not_in_appris)) {
      ## Get the gene id 
      gene_id <- getBM(filters = c("external_gene_name"), 
                       values = list(not_in_appris[i]), 
                       attributes = c("ensembl_gene_id"), 
                       mart = ensembl)[[1]]
      
      ## Make dummy data and rbind with original dataset
      appris_genes <- rbind(appris_genes, data.table("gene_name"=not_in_appris[i],
                                     "gene_id"=gene_id,
                                     "transcript_id"=NA,
                                     "CCDS"=NA,
                                     "status"="minor"))
    }
  }
  return(appris_genes)
}

## Function to get GTEx TPM data for appris isoforms 
get_appris_tpm <- function(appris_df, testis_GTEx_tpm) {
  
  ## Convert each row to its own list 
  if (class(appris_df)[1] == "data.table") {
    appris_df <- appris_df %>% transpose(.) %>% as.list(.) %>% unname(.)
  }
  
  ## Go through each appris isoform and get the average GTEx expression (TPM) in testis samples 
  x <- appris_df[[1]]
  appris_df_tpm <- lapply(appris_df, function(x, testis_GTEx_tpm) {
    print(paste0("Getting average GTEx TPM data for ", x[1], "'s transcript: ", x[3]))
    
    ##########################################################################################
    ##### SANITY: Check to see if the gene has an APPRIS designation (i.e. is not minor) #####
    ##########################################################################################
    
    if (x[5] == "minor") {
      warning_message <- paste0("Could not get transcript daata for ", x[1], "'s transcript: ", x[3], " as it has a \"minor\" isoform tag in APPRIS")
      print(warning_message)
    }
    
    if (x[5] != "minor") {
      ## Get the transcript's mean TPM from GTEx testis samples
      x[length(x) + 1] <- grep(x[3], testis_GTEx_tpm$transcript_id) %>% 
        testis_GTEx_tpm[.,(3:ncol(testis_GTEx_tpm)), with=FALSE]  %>% as.numeric(.) %>% mean(.) %>% round(., digits = 2)
      
      ## Assign appris binary
      x[7] <- 1
      
      ## Determine if the appris transcript is in the gtex dataset
      if (is.na(x[6])) {
        x[8] <- 0
      }
      if (!is.na(x[6])) {
        x[8] <- 1
      }
      
    }

    return(x)  
  }, testis_GTEx_tpm=testis_GTEx_tpm)
  
  return(appris_df_tpm)
}  
  
  
  
