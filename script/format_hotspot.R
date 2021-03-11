########################################################
##### Function to format hotspot mutation datasets #####
########################################################
## Load in packages
library(data.table)
library(dplyr)
library(biomaRt)
library(ensembldb)
library(pbmcapply)
library(rtracklayer)

## Define the function 
format_3D_hotspot <- function(hotspot_file, genes, merged_exons) {
  
  ## Read in the hotspot mutation file 
  hotspot_df <- fread(hotspot_file) 
  
  ## Get the genes of interest 
  hotspot_df_genes <- hotspot_df[Gene %in% genes, c("Gene", "Amino_Acid_Position"), with = FALSE]
  
  ## Get the ensembl human GRCh38 annotation dataaset
  edb <- EnsDb.Hsapiens.v86
  
  ## Iterate through the hotspot mutated genes to identify exons with mutation 
  hotspot_genes <- unique(hotspot_df_genes$Gene) 
  
  ## Get merged_exon data from non-hotspot genes
  merged_exons_non3d_hotspot <- merged_exons[!(gene_name %in% hotspot_genes)]
  merged_exons_non3d_hotspot$hotspot_3d_exon_num <- 0
  
  ## Define database to use (Ensembl Genes 102) and the dataset to use (Human genes GRCh38.p13)
  ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")
  
  #x <- hotspot_genes[4]
  merged_exons_3d_hotspot <- rbindlist(pbmclapply(mc.cores = 8, hotspot_genes, function(x, hotspot_df_genes, merged_exons, ensembl) {
    #print(x)
    ## Get the hotspot AA positions for the gene of interest 
    df <- unique(hotspot_df_genes[Gene == x])
    aa_pos <- df$Amino_Acid_Position

    ## Get the gene_id for the gene of interest 
    merged_exons_gene <- merged_exons[gene_name == x]
    merged_exons_gene$hotspot_3d_exon_num <- 0
  
    ## Get the transcript_ids of the gene_id of interest 
    t_ids <- merged_exons[gene_name == x, "transcript_id"][[1]]  %>% strsplit(x = ., split = " | ", fixed = TRUE) %>% unlist() %>% unique()
    
    ###################################################################      
    ##### Map hotspot mutation AA position to genomic coordinates #####
    ###################################################################
    ## Iterate through each aa position 
    for (w in 1:length(aa_pos)) {
      pos <- aa_pos[w]
      #print(w)
      ## Iterate through each transcript --> get protein_id and genomic coordinates of hotspot mutation
      for (i in 1:length(t_ids)) {
        #print(i)
        transcript <- t_ids[i]  
        
        ###################################################################      
        ##### Map hotspot mutation AA position to genomic coordinates #####
        ###################################################################
        ## Get the protein_id with biomaRt
        protein_id <- getBM(filters = c("ensembl_transcript_id"), 
                            values = list(transcript), 
                            attributes = c("ensembl_peptide_id"), 
                            mart = ensembl, useCache = FALSE)[[1]]
        
        ## Use the protein_id and hotspot mutation position to get genomic coordinates
        tmp <- IRanges(start = pos, end = pos, names = protein_id)
        p2g <- proteinToGenome(tmp, edb)
        
        ## Get the genomic coordinates 
        hotspot_genomic_start <- p2g[[protein_id]]@ranges@start
        hotspot_genomic_stop <- p2g[[protein_id]]@ranges@width + hotspot_genomic_start - 1
        
        hotspot_genomic_coord <- unique(c(hotspot_genomic_start, hotspot_genomic_stop))
        
        ## Iterate through each genomic coordinate and term which exon the hotspot start/stop is found in 
        for (y in 1:length(hotspot_genomic_coord)) {
          ## Get the coordinate 
          coord <- hotspot_genomic_coord[y]
          
          ## Get the row (i.e. exon) the hotspot coordinate is found in
          row <- which((coord >= merged_exons_gene$start) & (coord <= merged_exons_gene$stop))
          
          ## Change length(start/stop_in_exon) == 0 --> to 0
          if (length(row) == 0) {
            ## Do nothing
            
          }
          
          ## See if the start/stops are within multiple exons (i.e. start in 2 exons, stop in 2 exons)
          if (length(row) > 1) {
            print(paste0("AA: ", pos, "; transcript: ", transcript, "; gene_id: ", g_id, " - start/stop of hotspot mutation found in 1+ exon"))
            stop()
          } 
          
          # ## See if hotspot is not within an exon
          # if (row == 0) {
          #   ## Do nothing 
          # }
          
          ## See if both start/stop are within an exon 
          if (length(row) == 1) { ## Start/stop correspond to the same exon
            merged_exons_gene[row]$hotspot_3d_exon_num <- merged_exons_gene[row]$hotspot_3d_exon_num + 1 
          }
        }
      }
    }
    print(merged_exons_gene)
    return(merged_exons_gene)
  }, 
  hotspot_df_genes=hotspot_df_genes, merged_exons=merged_exons, ensembl=ensembl))
  
  return(rbind(merged_exons_non3d_hotspot, merged_exons_3d_hotspot))
}

format_hotspot <- function(hotspot_file, genes, merged_exons_3d_hotspot) {
  ## Read in the data 
  hotspot_df <- fread(hotspot_file, header = TRUE) %>% `colnames<-`(c("Hugo_Symbol", "Entrez_Gene_Id", "Chromosome", "Start_Position", "End_Position", "Transcript_ID"))

  ## Initiate the hotspot_transcript_count column
  merged_exons_3d_hotspot$hotspot_exon_num <- 0
  
  ## Define database to use (Ensembl Genes 102) and the dataset to use (Human genes GRCh38.p13)
  ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")
  
  ## Change H3F3A and H3F3B to H3-3A and H3-3B
  hotspot_df$Hugo_Symbol <- gsub(pattern = "H3F3A", replacement = "H3-3A", x = hotspot_df$Hugo_Symbol)
  hotspot_df$Hugo_Symbol <- gsub(pattern = "H3F3B", replacement = "H3-3B", x = hotspot_df$Hugo_Symbol)
  
  ## Get the hotspot genes
  hotspot_genes <- hotspot_df[Hugo_Symbol %in% genes]$Hugo_Symbol %>% unique() 
  
  
  ## Get the genes not in the hotspot dataset
  non_hotspot_genes <- genes[-which(genes %in% hotspot_genes)]
  non_hotspot_genes_df <- merged_exons_3d_hotspot[gene_name %in% non_hotspot_genes]
  
  ## Iterate through each hotspot gene
  ## Define the gene name 
  gene <- hotspot_genes[1]
  final <- rbindlist(lapply(hotspot_genes, function(gene, hotspot_df, ensembl, merged_exons_3d_hotspot){

    ## Get the MAF data for the entrez id
    hotspot_df_gene <- hotspot_df[Hugo_Symbol == gene]
    
    if (length(unique(hotspot_df_gene$Transcript_ID)) > 1) {
      print(gene)
      print("multiple_transcripts")
    }
    
    ## Get the genomic coordinates of the gene-MAF dataset
    ranges <- IRanges(start = hotspot_df_gene$Start_Position, 
                      end = hotspot_df_gene$End_Position, 
                      names = hotspot_df_gene$Transcript_ID)
    hg19_coord <- GRanges(seqnames=c(paste0("chr", hotspot_df_gene$Chromosome[1])), ranges)
    
    ## Get the liftover chain file 
    chain <- import.chain("~/Google Drive/Quinlan Lab - PhD/Projects/spermseq/twinstrand_target_gene_set/data/hg19ToHg38.over.chain")
    
    ## Get the liftover results 
    results <- data.table(liftOver(hg19_coord, chain))
    
    ## Store hg38 results  
    hg38_transcript <- ranges@NAMES
    hg38_coord_start <- results$V1@unlistData@ranges@start
    hg38_coord_end <- hg38_coord_start + results$V1@unlistData@ranges@width - 1
    
    ## Make into data.table
    hg38_coord_start <- data.table("hg38"=hg38_coord_start, "transcript_id"=hg38_transcript)
    hg38_coord_end <- data.table("hg38"=hg38_coord_end, "transcript_id"=hg38_transcript)
    hg38_coord <- unique(rbind(hg38_coord_start, hg38_coord_end)) ## Combine the coordinates 
    
    ## Go through each hg38_coord and see if it can be found in an exon of the gene of interest
    merged_exons_3d_hotspot_gene <- merged_exons_3d_hotspot[gene_name == gene]
    
    ## Skip if there the hotspot gene is not in the merge_exon dataset
    if (nrow(merged_exons_3d_hotspot_gene) == 0 ){
      return(data.table())
    }
    
    ## Iterate through each hg38 hotspot mutation
    for (w in 1:nrow(hg38_coord)) {
      ## Get the coordinate
      coord <- hg38_coord$hg38[w]
      transcript <- hg38_coord$transcript_id[w]
      
      ## Get the row (i.e. exon) in the gene dataset that the coord is found in 
      row <- which(which((coord >= merged_exons_3d_hotspot_gene$start & coord <= merged_exons_3d_hotspot_gene$stop)) == 
                     grep(transcript, merged_exons_3d_hotspot_gene$transcript_id)) 
      
      ## Change length(start/stop_in_exon) == 0 --> to 0
      if (length(row) == 0) {
        ## Do nothing
        
      }
      
      ## See if the start/stops are within multiple exons (i.e. start in 2 exons, stop in 2 exons)
      if (length(row) > 1) {
        print(paste0("AA: ", pos, "; transcript: ", transcript, "; gene_id: ", g_id, " - start/stop of hotspot mutation found in 1+ exon"))
        merged_exons_3d_hotspot_gene[row]$hotspot_exon_num <- merged_exons_3d_hotspot_gene[row]$hotspot_exon_num + 1 
        #stop()
      } 
      
      ## See if both start/stop are within an exon 
      if (length(row) == 1) { ## Start/stop correspond to the same exon
        merged_exons_3d_hotspot_gene[row]$hotspot_exon_num <- merged_exons_3d_hotspot_gene[row]$hotspot_exon_num + 1 
      }
    }
    
    ## Return the data
    return(merged_exons_3d_hotspot_gene)
    
  }, 
  hotspot_df=hotspot_df, ensembl=ensembl, merged_exons_3d_hotspot=merged_exons_3d_hotspot))
  
  final <- rbind(final, non_hotspot_genes_df)

  return(final)
}
