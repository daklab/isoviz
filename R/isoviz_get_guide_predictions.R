#' A function for obtaining junction-level guide RNA efficiency predictions for a specific gene or junction
#' Importantly, this function only works for annotated junctions in whatever gencode version you are using
#' @param gene ensembl ID for now
#' @param cell_type cell type of interest, default=hESC
#' @param guides_per_junction the number of highest ranking guides to print out, minimum of 1 and maximum of 8, default is 8
#' @param include_all_juncs print guides for all junctions of the gene, default is TRUE
#' @param include_specific_junctions provide list of junctions for guide RNA prediction, can only specify if include_all_juncs is FALSE
#' @param leafcutter_input output from running isoviz_minicutter
#' @return table
#' @examples
#' isoviz_get_guide_predictions()
#' @name isoviz_get_guide_predictions
#' @export

# Updated: Function to get guide RNA predictions for specific junctions
# merge junction info and predictions
# generate this in advance and have as a table
#v41_info = read_tsv("~/knowles_lab/cas13_share/gencode_v41_all_junction_guides.txt", col_names = TRUE) %>%
#  mutate(guide_sequence = as.character(guide_sequence)) # 2,255,040
#predictions = read_csv("~/knowles_lab/cas13_share/230428_predictions/combined-corrected-predictions-gencode.csv") %>%
#  select(guide_id, predicted_lfc) # 2,254,874

#gencode_predictions = v41_info %>% left_join(predictions, by = "guide_id")
#save(gencode_predictions, file = "~/for_Karin/gencode_predictions.rda")

# #' @import gencode_predictions, juncs_recluster
#require(rempsyc)

junctions = c("junc188584", "junc188587", "junc188588", "junc188585")

isoviz_get_guide_predictions = function(gene = "ENSG00000100320", cell_type = "hESC",
                                 guides_per_junction = 8,
                                 include_all_juncs = TRUE,
                                 include_specific_junctions = c(),
                                 leafcutter_input){

  # 1. load most up-to-date gencode_predictions file (this will be updated)
  gencode_predictions <- system.file("data", "gencode_predictions.rda.gz", package="isoviz")
  load(gencode_predictions)

  # ranking of prediction score -> values 0 to 1
  min_score = min(gencode_predictions$predicted_lfc)
  max_score = max(gencode_predictions$predicted_lfc)
  gencode_predictions %<>% mutate(score = 1 - (predicted_lfc - min_score) / (max_score - min_score)) %>% arrange(score)

  # filter here if the user wants to filter for specific junctions
  if(include_all_juncs == FALSE){
    if(length(include_specific_junctions) == 0){
      stop("Must include list of junctions if include_all_juncs == FALSE")
    } else{
      filtered_predictions = gencode_predictions %>% filter(junc_id %in% include_specific_junctions)
    }
  } else{
    filtered_predictions = gencode_predictions %>% filter(grepl(gene, gene_id))
  }

  gene_name = unique(filtered_predictions$gene_name)

  # leafcutter junction counts minicutter input load
  juncs_recluster = leafcutter_input
  juncs_recluster %<>% select(everything(), junc.counts = readcount, Usage = usage_ratio) %>% dplyr::group_by(cluster_idx) %>%
    dplyr::mutate(cluster.counts = sum(junc.counts), Usage = round(Usage, 2)) %>% ungroup()

  all_guide_info = filtered_predictions %>%
    left_join(juncs_recluster, by = c("chr" = "chrom", "junc_start" = "start", "junc_end" = "end", "strand")) %>%
    arrange(junc_id, score)

  # filter based on parameter for number of guides to report and change NA in counts to 0
  filtered_guide_info = all_guide_info %>% dplyr::group_by(junc_id) %>% slice_max(score, n = guides_per_junction) %>% ungroup()
  filtered_guide_info$Transcripts <- gsub(paste0(gene_name, "-"), "", filtered_guide_info$transcript_isoforms)
  filtered_guide_info %<>% dplyr::select(Cluster = cluster_idx, JuncID = junc_id, Category = junction_category, Transcripts, Counts = junc.counts, Usage,
                                  GuideSeq = guide_sequence, Score = score)
  filtered_guide_info %<>% dplyr::mutate(Counts = ifelse(is.na(Counts), 0, Counts))

  # don't add sig digits
  fun <- function(x) {formatC(x, format = "f", digits = 0)}

  return(nice_table(filtered_guide_info, col.format.custom = 1:6, format.custom = "fun",
                    title = paste0(gene_name, ": Top ", guides_per_junction, " gRNAs per Junction")))
}

#guide_table = get_guide_predictions(gene = "ENSG00000100320", include_all_juncs = FALSE, include_specific_junctions = c("junc188584", "junc188587", "junc188588", "junc188585"))

#flextable::save_as_docx(guide_table, path = "~/isoviz/plots/230420_rbfox2_guides_table.docx")
