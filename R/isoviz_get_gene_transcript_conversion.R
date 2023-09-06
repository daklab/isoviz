#' A function for converting gtf to psl 
#'
#' Converts a gtf to psl, depending on the output filename extension;
#' gtf exons need to be grouped by transcript and sorted by coordinate w/in a transcript
#' this function is adapted from Flair https://github.com/BrooksLabUCSC/flair
#'
#' @param gtf_file_path User provided gtf file 
#' @param output_file Name for the txt file to be outputted. Default is "gene_transcript_convert.txt"
#' @return table with gene ids, trancsript ids, names... 
#' @examples
#' # Example with gtf file 
#'
#' @name isoviz_get_gene_transcript_conversion
#' @export

isoviz_get_gene_transcript_conversion = function(gtf_file, output_file="gene_transcript_convert.txt"){
  
  gr <- import.gff(gtf_file)
  gtf_data <- as.data.table(gr)
  print("Keeping only protein-coding and lncRNA genes")
  gtf_data = dplyr::filter(gtf_data, gene_type %in% c("protein_coding", "lncRNA"))
  
  # keep only the following columns 
  exons_df = gtf_data %>% dplyr::select(gene_id, transcript_id, gene_name,
                                        transcript_name, gene_type, transcript_type)
  exons_df = unique(exons_df)
  
  # filter rows where transcript_id is NA 
  exons_df = dplyr::filter(exons_df, !(is.na(transcript_id)))
  colnames(exons_df)[which(colnames(exons_df) == "transcript_id")] = "trans_id"

  # Print the result
  print("The number of transcripts in each category is:")
  print(table(exons_df$gene_type))
  write.table(exons_df, file=output_file, quote=F, sep="\t", row.names=F)
  print(paste("Saving the conversion file to: ", output_file, sep=" "))
}
