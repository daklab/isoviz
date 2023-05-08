#' A function for obtaining all necessary coordinates for leafcutter integration
#' and plotting isoforms given a gtf file or psl file from long read sequencing
#'
#' Inputs gencode gtf or psl file, outputs exon coords (for plotting) and
#' intron coords (for leafcutter integration) with ENSG/ENST names plus
#' meaningful names (ex RBFOX2)
#' @param file_path Gencode gtf or psl file
#' @param gene_trans Dataframe with ensembl gene and transcript ids and names
#' @return list of exon and intron coordinates in that order.
#' @examples
#' isoviz_coords()
#' @name isoviz_coords
#' @import data.table
#' @export

isoviz_coords = function(file_path, gene_trans,
                         input_type="psl"){

  if(input_type=="psl"){

    #print(paste("The input file being used is"), file_path)
    genome_data <- fread(file_path)
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
    # exon start = intron end (execept for the first exon)

    intron_starts = genome_data %>%
            group_by(trans_id) %>%
            slice(1:(n()-1)) %>%
            select(blockends)

    intron_ends = genome_data %>%
            group_by(trans_id) %>%
            slice(2:n()) %>%
            select(blockstarts)

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
    trans_info = genome_data %>% select(1, 4, 5, 6) %>% distinct() %>% filter(chr != "chrY", chr != "chrM")
    intron_data = intron_data %>% left_join(trans_info,  by = "trans_id") %>%
      left_join(convert, by = c("gene_id", "trans_id")) %>% select(4, 2, 3, 5, 1, 6, 7, 8, 9, 10)

    # Final intron dataset
    iso_intron_data = as.data.table(intron_data)

    # Final exon dataset
    iso_exon_data = as.data.table(genome_data %>% left_join(convert, by = c("gene_id", "trans_id")))

    } #if statement for psl file

  if(input_type=="gtf"){
  }

  return(list(iso_exon_data, iso_intron_data))
}
