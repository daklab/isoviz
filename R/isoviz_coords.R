#' A function for obtaining all necessary coordinates for leafcutter integration
#' and plotting isoforms given a gtf file or psl file from long read sequencing
#'
#' Inputs gencode gtf or psl file, outputs exon coords (for plotting) and
#' intron coords (for leafcutter integration) with ENSG/ENST names plus
#' meaningful names (ex RBFOX2)
#' @param file_path Gencode gtf or psl file
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
                         input_type="psl", count_file = "", min_count = 5){
  if(input_type=="gtf"){
    isoviz_gtf_to_psl(file_path, psl_output_file = "converted_gtf.psl", chrom_sizes=NULL)
    genome_data <- fread(psl_output_file)
  }
  
  if(input_type=="psl"){
    #print(paste("The input file being used is"), file_path)
    genome_data <- fread(file_path)
  }
    
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

  print(paste(length(unique(genome_data$trans_id)), "transcripts with at least one intron!"))

  # for each transcript
  # exon end = intron start (except for last exon)
  # exon start = intron end (except for the first exon)

  intron_starts = genome_data %>%
          group_by(trans_id) %>%
          dplyr::slice(1:(n()-1)) %>%
          dplyr::select(blockends)

  intron_ends = genome_data %>%
          group_by(trans_id) %>%
          dplyr::slice(2:n()) %>%
          dplyr::select(blockstarts)

  intron_ends$trans_id = NULL
  intron_data = cbind(intron_starts, intron_ends)
  intron_data = as.data.frame(intron_data)
  colnames(intron_data) = c("trans_id", "intron_starts", "intron_ends")

  # Important note: this no longer has single exon transcripts
  intron_data$intron_ends = intron_data$intron_ends+1
  print(paste(length(unique(intron_data$trans_id)), "transcripts with at least one intron!"))

  # Convert to gene and transcript names (should this be placed somewhere else?)
  convert = read_tsv(gene_trans, col_names = TRUE)

  # Should avoid selecting columns this way
  #trans_info = genome_data %>% dplyr::select(1, 4, 5, 6) %>% distinct() %>% filter(chr != "chrY", chr != "chrM")
  #intron_data = intron_data %>% left_join(trans_info,  by = "trans_id") %>%
  #  left_join(convert, by = c("gene_id", "trans_id")) %>% select(4, 2, 3, 5, 1, 6, 7, 8, 9, 10)
  trans_info = genome_data %>% dplyr::select(chr, trans_id, gene_id, strand) %>% distinct() %>% filter(chr != "chrY", chr != "chrM")
  intron_data %<>% left_join(trans_info,  by = "trans_id") %>%
    left_join(dplyr::select(convert, "gene_id", "gene_name", "gene_type") %>% distinct, by = "gene_id") %>%
    left_join(convert, by = c("gene_id", "trans_id", "gene_name", "gene_type")) %>% 
    dplyr::select(chr, intron_starts, intron_ends, gene_id, trans_id, strand, gene_name, transcript_name, gene_type, transcript_type)
    
  id = intron_data %>% dplyr::select(gene_id, trans_id) %>% distinct() %>% arrange(gene_id, desc(trans_id)) %>%
    group_by(gene_id) %>% mutate(id = 1:n()) %>% ungroup()
    
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
