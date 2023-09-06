#' A function for converting gtf to psl 
#'
#' Converts a gtf to psl, depending on the output filename extension;
#' gtf exons need to be grouped by transcript and sorted by coordinate w/in a transcript
#' this function is adapted from Flair https://github.com/BrooksLabUCSC/flair
#'
#' @param gtf_file_path User provided gtf file 
#' @param psl_output_file Name for the PSL file to be outputted. Default is "converted_gtf.psl"
#' @param chrom_sizes Chrom sizes file for psl conversion, recommended
#' @param filter_CDS if "yes" then will only include CDS exons in the PSL output file
#' @return psl formatted file
#' @examples
#' # Example with gtf file 
#'
#' @name isoviz_gtf_to_psl
#' @export

isoviz_gtf_to_psl = function(gtf_file_path, 
                             psl_output_file="converted_gtf.psl", 
                             specie="human", 
                             filter_CDS="no",
                             chrom_sizes=NULL){
  
  # Read the GTF file 
  # Use import.gff to read the GFF file into a GRanges object
  gr <- import.gff(gtf_file_path)
  gtf_data <- as.data.table(gr)
  
  # Extract relevant columns
  print("Keeping only protein-coding and lncRNA genes")
  gtf_data = dplyr::filter(gtf_data, gene_type %in% c("protein_coding", "lncRNA"))
  
  # keep only the following columns 
  exons_df = gtf_data %>% dplyr::select(seqnames, type, start, end, 
                                        strand, transcript_id, gene_id)
  
  colnames(exons_df) = c("chrom", "ty", "start", "end", "strand", "transcript_id", "gene_id")
  exons_df$start = exons_df$start-1

  # Filter for ty == "CDS"
  if (filter_CDS == "yes"){
    z <- which(exons_df$ty == "CDS")
    cds_df <- exons_df[z,]
    print("Filtered to just include CDS entries")
    
  } else{
    z <- which(exons_df$ty == "exon")
    cds_df = exons_df[z,] # keep all exons
    print("Filtered to include all exons")
  }
  
  # Group and summarize exon information
  cols_to_convert <- c("chrom", "ty", "strand", "transcript_id", "gene_id")
  cds_df = as.data.frame(cds_df)
  # Use lapply to convert specified columns
  cds_df[cols_to_convert] <- lapply(cds_df[cols_to_convert], as.character)
  
  exons_grouped <- cds_df %>%
    dplyr::group_by(transcript_id)
  
  # Individual transcript groups with corresponding exons 
  individual_groups <- exons_grouped %>%
    group_split()
  
  # Open the PSL output file for writing
  con <- file(psl_output_file, "w")
  
  print("Starting to sort through transcripts and write entries to PSL file!")
  
  # Initialize the progress bar
  pb <- progress_bar$new(
    format = "[:bar] :percent :current/:total (:eta)",
    total = length(individual_groups)
  )
  
  print("Initializing progress bar!")
  
  for (i in seq_along(individual_groups)) {
    
    pb$tick()  # Advance the progress bar
    test_ex <- as.data.frame(individual_groups[[i]])
    blockstarts <- list()
    blocksizes <- list()
    
    # Append start coords of CDS in transcript
    blockstarts <- c(blockstarts, test_ex$start)
    # Append sizes of CDS in transcript 
    blocksizes <- c(blocksizes, test_ex$end - test_ex$start)
    blockcount <- length(blockstarts)
    
    # Check if the start of the first CDS block is greater than the start of the next one 
    if (blockcount > 1 && blockstarts[[1]] > blockstarts[[2]]) {
      # If so, reverse the sizes and the starts
      blocksizes <- rev(blocksizes)
      blockstarts <- rev(blockstarts)
    }
    
    # get main start site of transcript (first CDS start)
    tstart <- blockstarts[[1]]
    # get last position in transcript (final blockstart + final blocksize)
    tend <- blockstarts[[length(blockstarts)]] + blocksizes[[length(blocksizes)]]
    # total sizes of CDS blocks in transcripts 
    qsize <- sum(unlist(blocksizes))
    # combine transcript and gene ID
    qname <- paste(test_ex$transcript_id[1], '_', test_ex$gene_id[1], sep = "")
    pos <- 0
    qstarts <- list(pos)
    
    if (length(blocksizes) > 1) {
      for (i in 1:(length(blocksizes) - 1)) {
        b <- blocksizes[[i]]
        pos <- pos + b
        qstarts <- c(qstarts, pos)
      }
    }
    
    if (length(blocksizes) == 1) {
      pos <- pos + blocksizes[[1]]
      qstarts <- c(qstarts, pos)
    }
    
    qstarts_str <- paste(qstarts, collapse = ",")
    qstarts_str <- paste(qstarts_str, ",", sep = "")
    
    blocksizes_str <- paste(blocksizes, collapse = ",")
    blocksizes_str <- paste(blocksizes_str, ",", sep = "")
    
    blockstarts_str <- paste(blockstarts, collapse = ",")
    blockstarts_str <- paste(blockstarts_str, ",", sep = "")
    
    # Construct the PSL formatted output line
    psl_line <- c(
      0, 0, 0, 0, 0, 0, 0, 0, test_ex$strand[1], qname, qsize, 0, qsize,
      test_ex$chrom[1], 0, tstart, tend, blockcount, blocksizes_str, qstarts_str, blockstarts_str)
    
    # Convert psl_line to a tab-separated string
    psl_line_str <- paste(psl_line, collapse = "\t")
    
    if (nzchar(psl_line_str)) {
      # Write the line to the file
      writeLines(psl_line_str, con)
    }
  
  }
  
  close(con)
  print("Done generating PSL file")
}
