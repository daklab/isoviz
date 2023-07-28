#' A function for obtaining junction-level guide RNA efficiency predictions for a specific gene or junction
#' @param gene ensembl ID for now
#' @param cell_type cell type of interest, default=hESC
#' @param guides_per_junction the number of highest ranking guides to print out, minimum of 1 and maximum of 8, default is 8
#' @param output_format Either return the table as a dataframe or flextable, default=dataframe
#' @param include_specific_junctions provide list of junctions for guide RNA prediction
#' @param leafcutter_input output from running isoviz_minicutter
#' @param output_format set to dataframe to obtain a regular dataframe object, to obtain a nice table, choose "nice_table"
#' @return table
#' @examples
#' isoviz_get_guide_predictions()
#' @name isoviz_get_guide_predictions
#' @export

isoviz_get_guide_predictions = function(gene = "ENSG00000100320", cell_type = "hESC",
                                 guides_per_junction = 8,
                                 include_specific_junctions = c(),
                                 leafcutter_input, 
                                 output_format = "dataframe"){

  print("Loading gencode sgRNA prediction file")
  # 1. load most up-to-date gencode_predictions file (this will be updated)
  gencode_predictions <- system.file("data", "230728_gencodev41_final_TIGER_predictions_to_compare.rda.gz", package="isoviz")
  gencode_predictions <- gzfile(gencode_predictions, "rb")
  load(gencode_predictions)
  gencode_predictions = tiger_pred_out

  # ranking of prediction score -> values 0 to 1
  # min_score = min(gencode_predictions$predicted_lfc)
  # max_score = max(gencode_predictions$predicted_lfc)
  # gencode_predictions %<>% mutate(score = 1 - (predicted_lfc - min_score) / (max_score - min_score)) %>% arrange(score)
  
  # filter here if the user wants to filter for specific junctions
  if(length(include_specific_junctions) != 0){
    filtered_predictions = gencode_predictions %>% filter(junc_id %in% include_specific_junctions)
  }else{
    filtered_predictions = gencode_predictions %>% filter(grepl(gene, gene_id))
  }

  gene_name = unique(filtered_predictions$gene_name)

  # leafcutter junction counts minicutter input load
  juncs_recluster = leafcutter_input
  juncs_recluster %<>% select(everything(), junc.counts = readcount, Usage = usage_ratio) %>% dplyr::group_by(cluster_idx) %>%
    dplyr::mutate(cluster.counts = sum(junc.counts), Usage = round(Usage, 2)) %>% ungroup()

  all_guide_info = filtered_predictions %>%
    left_join(juncs_recluster, by = c("chr" = "chrom", "start", "end", "strand")) %>%
    arrange(junc_id, TIGER_score)

  # filter based on parameter for number of guides to report and change NA in counts to 0
  filtered_guide_info = all_guide_info %>% dplyr::group_by(junc_id) %>% slice_max(TIGER_score, n = guides_per_junction) %>% ungroup()
  filtered_guide_info$Transcripts <- gsub(paste0(gene_name, "-"), "", filtered_guide_info$transcript_isoforms)
  filtered_guide_info %<>% dplyr::select(Cluster = cluster_idx, JuncID = junc_id, Category = junction_category, Transcripts, Counts = junc.counts, Usage,
                                  GuideSeq = guide_sequence, TIGER_score)
  filtered_guide_info %<>% dplyr::mutate(Counts = ifelse(is.na(Counts), 0, Counts))
  #print(head(filtered_guide_info))
  print("Closed gencode sgRNA prediction file")
  #print(gene_name)
  #print(guides_per_junction)
  # don't add sig digits
  #usetext <- function(x) {
  #  formatC(x, format = "f", digits = 0)
  #}

  #table_return = nice_table(filtered_guide_info, col.format.custom = 1:6,
  #                          format.custom = usetext,
  #           title = paste0(gene_name, ": Top ", guides_per_junction, " gRNAs per Junction"))
  if(output_format == "dataframe"){
    table_return = filtered_guide_info
  }
  if(output_format == "nice_table"){
    table_return = nice_table(filtered_guide_info,
                              title = paste0(gene_name, ": Top ", guides_per_junction, " gRNAs per Junction"))
  }

  return(table_return)
}

#junctions = c("junc188584", "junc188587", "junc188588", "junc188585")
#guide_table = isoviz_get_guide_predictions(gene = "ENSG00000100320", include_all_juncs = FALSE, include_specific_junctions = c("junc188584", "junc188587", "junc188588", "junc188585"))
#flextable::save_as_docx(guide_table, path = "~/isoviz/plots/230420_rbfox2_guides_table.docx")
