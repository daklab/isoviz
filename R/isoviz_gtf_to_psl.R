#' A function for converting gtf to psl 
#'
#' Converts a gtf to psl, depending on the output filename extension;
#' gtf exons need to be grouped by transcript and sorted by coordinate w/in a transcript
#' this function is adapted from Flair https://github.com/BrooksLabUCSC/flair
#'
#' @param gtf_file_path User provided gtf file 
#' @param psl_output_file Name for the PSL file to be outputted. Default is "converted_gtf.psl"
#' @param chrom_sizes Chrom sizes file for psl conversion, recommended
#' @return psl formatted file
#' @examples
#' # Example with gtf file 
#'
#' @name isoviz_gtf_to_psl
#' @export

#  library(GenomicFeatures)
#  library(progress)
#  library(data.table)
#  library(R.utils)      
#  library(dplyr)

# gtf_file_path = "/Users/kisaev/Downloads/gencode.v41.basic.annotation.gtf.gz"

isoviz_gtf_to_psl = function(gtf_file_path, psl_output_file="converted_gtf.psl", specie="human", chrom_sizes=NULL){
  
  # Read the GTF file using fread
  gtf_data <- fread(gtf_file_path, header = FALSE)
  
  # Filter out lines that start with '#'
  gtf_data <- gtf_data[!grepl("^#", V1)]  
  
  # Extract relevant columns
  print("Extracting info from gtf file to compile all exon information...")
  
  if (specie=="human") {
    attributes <- gtf_data$V9
    trans_id <- gsub('.*"(ENST\\d+\\.\\d+)".*', '\\1', attributes)
    this_transcript <- sub(".*\"(ENST\\d+\\.\\d+)\".*", "\\1", trans_id)
    gene_ids <- gsub('.*"(ENSG\\d+\\.\\d+)".*', '\\1', attributes)}
  
  if (specie=="mouse") {
    attributes <- gtf_data$V9
    trans_id <- gsub('.*"(ENSMUST\\d+\\.\\d+)".*', '\\1', attributes)
    this_transcript <- sub(".*\"(ENSMUST\\d+\\.\\d+)\".*", "\\1", trans_id)
    gene_ids <- gsub('.*"(ENSMUSG\\d+\\.\\d+)".*', '\\1', attributes)}
  
  # Create an empty data frame to store exon information
  exons_df <- data.frame(
    chrom = gtf_data$V1,
    ty = gtf_data$V3,
    start = as.integer(gtf_data$V4) - 1,
    end = as.integer(gtf_data$V5),
    strand = gtf_data$V7,
    trans_id = this_transcript,
    transcript_id = this_transcript, 
    gene_id = gene_ids
  )
  
  # Filter for ty == "CDS"
  z <- which(exons_df$ty == "CDS")
  cds_df <- exons_df[z,]

  print("Filtered to just include CDS entries")
  
  # Group and summarize exon information
  exons_grouped <- cds_df %>%
    dplyr::group_by(transcript_id)
  
  # Individual transcript groups with corresponding CDS 
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
    qname <- paste(test_ex$trans_id[1], '_', test_ex$gene_id[1], sep = "")
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
      #cat("\n", file = con)
    }
  
  }
  
  close(con)
  print("Done generating PSL file")
}
