#' A function for converting gtf to psl 
#'
#' Converts a gtf to psl, depending on the output filename extension;
#' gtf exons need to be grouped by transcript and sorted by coordinate w/in a transcript
#' this function is adapted from Flair https://github.com/BrooksLabUCSC/flair
#'
#' @param gtf_file_path User provided gtf file 
#' @param output_file Name for the PSL file to be outputted. Default is "converted_gtf.psl"
#' @return table with gene ids, trancsript ids, names... 
#' @examples
#' # Example with gtf file 
#'
#' @name isoviz_get_gene_transcript_conversion
#' @export

isoviz_get_gene_transcript_conversion = function(gtf_file, output_file="gene_transcript_convert.txt"){
  
  gtf = read_tsv(gtf_file, col_names = FALSE, skip = 5) %>% dplyr::filter(X3 == "transcript") 

  attributes <- gtf$X9
  print(length(attributes))
  
  # Function to extract specific attributes
  extract_attributes <- function(attr_string) {
    attr_list <- strsplit(attr_string, "; ")[[1]]
    attr_pairs <- sapply(attr_list, function(attr) {
      key_value <- strsplit(attr, " ", fixed = TRUE)
      list(key = key_value[[1]][1], value = key_value[[1]][2])
    })
    
    attr_df <- as.data.frame(matrix(unlist(attr_pairs), ncol = 2, byrow = TRUE))
    colnames(attr_df) <- c("Attribute", "Value")
    specific_attributes <- c("gene_id", "gene_name", "transcript_id", "transcript_name", "gene_type", "transcript_type")
    specific_data <- attr_df[attr_df$Attribute %in% specific_attributes, ]
    specific_data <- spread(specific_data, Attribute, Value)
    #specific_data <- specific_data %>%
    #  mutate_all(~ str_replace(., "^\"|\"$", "")) %>% 
    #  mutate_all(~ str_replace(., "^\"|\"$", ""))
    # clean up values 
    return(specific_data)
  }
  
  # Apply the function to all attributes
  all_extracted <- llply(attributes, extract_attributes, .progress="text")
  
  # Combine all extracted attributes into a single data frame
  final_df <- as.data.table(ldply(all_extracted))
  final_df <- final_df %>%
    final_df(~ str_replace(., "^\"|\"$", "")) %>% 
    final_df(~ str_replace(., "^\"|\"$", ""))
  
  # need to rename "transcript_id" column to "trans_id"
  # Print the result
  print(final_df)
  write.table(final_df, file=output_file, quote=F, sep="\t", row.names=F)
  print(paste("Saving the conversion file to: ", output_file, sep=" "))
}
