#' A function for obtaining all necessary coordinates for leafcutter integration
#' and plotting isoforms given a gtf file or psl file from long read sequencing
#'
#' Inputs gencode gtf or psl file, outputs exon coords (for plotting) and
#' intron coords (for leafcutter integration) with ENSG/ENST names plus
#' meaningful names (ex RBFOX2)
#' @param file_path Gencode PSL file! (can be custom PSL file though)
#' @param gene_trans Dataframe with ensembl gene and transcript ids and names
#' @param input_type Either 'psl' or 'gtf' for input file type. Default is psl.
#' @param count_file User can provide a read count file per transcript from a long read experiment
#' @param min_count Require a specific number of reads per transcript. Only works if providing count file. Default is 5. 
#' @return list of exon and intron coordinates in that order.
#' @examples
#' # Example with psl file and pre-loaded gene to transcript conversions
#' psl_file <- system.file("data", "gencode.v41.basic.annotation.psl", package="isoviz")
#' gene_trans <- system.file("data", "gencode_v41_gene_transcript_convert.txt", package="isoviz")
#'
#' results <- isoviz_coords(psl_file, gene_trans, input_type = "psl")
#' print(head(results)[1]) #exons coords
#' print(head(results)[2]) #introns coords
#'
#' @name isoviz_coords
#' @export


isoviz_coords = function(file_path, gene_trans,
                         count_file = "", min_count = 5){
  
  print("Please ensure your input file is a PSL file~")
  genome_data <- fread(file_path)
  
  # Set options to avoid scientific notation
  options(scipen = 999)
  
  genome_data$gene_id = sapply(genome_data$V10, function(x){strsplit(x, "_")[[1]][2]})
  genome_data$trans_id = sapply(genome_data$V10, function(x){strsplit(x, "_")[[1]][1]})
  
  # Get coordinates for each "block"
  genome_data$blocksizes = sapply(genome_data$V19, function(x){paste(strsplit(x, ",")[[1]], collapse=",")})
  genome_data$blockstarts = sapply(genome_data$V21, function(x){paste(strsplit(x, ",")[[1]], collapse=",")})
  genome_data$strand = genome_data$V9
  genome_data$chr = genome_data$V14
  genome_data$start = genome_data$V16
  genome_data$end = genome_data$V17
  genome_data$transcript_length = genome_data$end - genome_data$start

  # Clean up and seperate blocks into rows
  # Need to do different adjustments here to get intron coords that will line up with leafcutter
  genome_data = genome_data %>% dplyr::select(chr, start, end, trans_id, gene_id,strand, blocksizes, blockstarts,
                                              transcript_length) %>% separate_rows(blocksizes, blockstarts)
  
  genome_data$blockstarts=as.numeric(genome_data$blockstarts)
  genome_data$blocksizes=as.numeric(genome_data$blocksizes)
  genome_data$blockends = genome_data$blockstarts + genome_data$blocksizes
  genome_data = as.data.table(genome_data)
  genome_data = unique(genome_data)
  
  transcript_ids = unique(genome_data$trans_id)
  print(paste(length(unique(genome_data$trans_id)), "transcripts with at least one intron!"))
  
  # Split the data frame by trans_id
  list_of_data <- split(genome_data, genome_data$trans_id)
  
  # Process each data frame in the list
  result <- lapply(list_of_data, function(data){
    new_ends <- data$blockstarts[-1] + 1
    len <- nrow(data)
    new_starts <- data$blockends[-len]
    n <- len - 1
    
    list(
      trans_id = rep(data$trans_id[1], n),
      intron_starts = new_starts,
      intron_ends = new_ends
    )
  })
  
  # Convert the list of lists back to three vectors
  trans_id <- unlist(lapply(result, `[[`, "trans_id"))
  intron_starts <- unlist(lapply(result, `[[`, "intron_starts"))
  intron_ends <- unlist(lapply(result, `[[`, "intron_ends"))
  intron_data = data.frame(trans_id, intron_starts, intron_ends) # 723,230
  intron_data = unique(intron_data)
  print(paste(length(unique(intron_data$trans_id)), "transcripts with at least one intron!"))
  
  # Convert to gene and transcript names (should this be placed somewhere else?)
  convert = as.data.table(read_tsv(gene_trans, col_names = TRUE, show_col_types = FALSE))
  convert = unique(convert)
  
  trans_info = as.data.table(genome_data %>% dplyr::select(chr, trans_id, gene_id, strand) %>% 
    distinct() %>% filter(chr != "chrY", chr != "chrM"))
  trans_info = unique(trans_info)
  
  intron_data %<>% left_join(trans_info,  by = "trans_id") %>%
    left_join(dplyr::select(convert, "gene_id", "gene_name", "gene_type") %>% 
                distinct, by = "gene_id") %>%
    left_join(convert, by = c("gene_id", "trans_id", "gene_name", "gene_type")) %>% 
    dplyr::select(chr, intron_starts, intron_ends, gene_id, trans_id, strand, gene_name, transcript_name, gene_type, transcript_type)
  
  id = intron_data %>% 
      dplyr::select(gene_id, trans_id) %>% 
      distinct() %>% 
      dplyr::arrange(gene_id, desc(trans_id)) %>%
      group_by(gene_id) %>% 
      dplyr::mutate(id = 1:n()) %>% 
      ungroup()
  
  intron_data$transcript_name[intron_data$transcript_name == ""] <- NA
  intron_data %<>% left_join(id, by = c("gene_id", "trans_id"))  %>% 
    mutate(transcript_name = ifelse(is.na(transcript_name), paste0(gene_name, "-novel", id), transcript_name)) %>% 
    ungroup() %>% dplyr::select(-id)
  
  # if a count file is provided, incorporate that here and filter
  if(count_file != ""){
    counts = read_tsv(count_file, col_names = FALSE)
    counts$gene_id = sapply(counts$X1, function(x){strsplit(x, "_")[[1]][2]})
    counts$trans_id = sapply(counts$X1, function(x){strsplit(x, "_")[[1]][1]})
    counts %<>% dplyr::select(gene_id, trans_id, read_counts = X2)
    intron_data %<>% left_join(counts, by = c("gene_id", "trans_id"))
    
    intron_data %<>% filter(read_counts >= min_count)
    
    # tpm
    per_million = sum(counts$read_counts)/1000000
    intron_data %<>% mutate(tpm = round(read_counts/per_million, 2))
    
  }
  
  # get file of all transcripts and their transcript names, mostly care about the novel id's here
  novel_ids = intron_data %>% dplyr::select(trans_id, novel_name = transcript_name) %>% distinct()
  
  print("Printing dimensions of intron data file")
  print(dim(intron_data))
  
  # Final intron dataset
  iso_intron_data = as.data.table(intron_data)
  
  # Final exon dataset
  genome_data %<>% 
    left_join(dplyr::select(convert, "gene_id", "gene_name", "gene_type") %>% distinct, by = "gene_id") %>%
    left_join(convert, by = c("gene_id", "trans_id", "gene_name", "gene_type")) %>% 
    filter(trans_id %in% intron_data$trans_id) %>%
    left_join(novel_ids, by = "trans_id") %>%
    mutate(transcript_name = ifelse(is.na(transcript_name), novel_name, transcript_name)) %>% dplyr::select(-novel_name)
  iso_exon_data = as.data.table(genome_data)
  
  return(list(iso_exon_data, iso_intron_data))
}
