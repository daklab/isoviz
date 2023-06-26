#' A function for mapping junctions to transcript isoforms
#' @param gene_intron_coords obtained intron coordinates for gene of interest
#' @param leafcutter_input output from running isoviz_minicutter
#' @param cell_type cell type of interest, default=hESC
#' @param gencode_intron_all_data expanded annotations of introns/transcripts
#' @return mapping
#' @examples
#' isoviz_map_junctions()
#' @name isoviz_map_junctions
#' @export

#library(data.table)
#library(plyr)
#library(GenomicFeatures)
#library(tidyr)

# map junctions function
isoviz_map_junctions = function(cell_type = "hESC", gene_intron_coords, leafcutter_input, gencode_intron_all_data){

  # note make sure to use the junction IDs from the gencode_intron_all_data file
  # join intron_coords with leafcutter data and mark any unannotated junctions
  juncs_recluster = leafcutter_input
  juncs_recluster$name = NULL
  print("removing junction IDs obtained from clustering introns, using default names from gencode file")

  # intron_coords rename columns
  gene_iso_data = gene_intron_coords

  # get total counts for clusters across their junctions
  juncs_recluster %<>% dplyr::select(everything(), junc.counts = readcount) %>% dplyr::group_by(cluster_idx) %>% dplyr::mutate(cluster.counts = sum(junc.counts)) %>% ungroup()

  # since leafcutter junctions dont have gene annotations, need to combine here first
  gene_cluster = gene_iso_data %>%
    left_join(juncs_recluster, by = c("chr" = "chrom", "intron_starts" = "start", "intron_ends" = "end", "strand")) %>%
    arrange(cluster_idx)

  # get additional columns from the gencode_intron_all_data file
  gene_cluster = gene_cluster %>%
    left_join(gencode_intron_all_data, by=c("chr" = "chr", "intron_starts" = "junc_start", "intron_ends"="junc_end", "strand", "gene_id", "gene_name", "gene_type"))

  # the column 'transcript_isoforms' has pre-mapped transcript isoforms that map those junctions/introns
  name = unique(gene_cluster$gene_name)
  gene_id = unique(gene_cluster$gene_id)
  gene_strand = unique(gene_cluster$strand)

  # to get information on leafcutter junctions that don't match annotation, still need to plot them
  # but they just won't get transcripts assigned to them
  all_gene_clusters = juncs_recluster %>% dplyr::filter(cluster_idx %in% gene_cluster$cluster_idx) %>%
    left_join(gene_iso_data, by = c("chrom" = "chr", "start" = "intron_starts", "end" = "intron_ends", "strand")) %>%
    filter(is.na(gene_id)) %>%
    dplyr::mutate(gene_name = unique(gene_cluster$gene_name), strand = unique(gene_cluster$strand),
           gene_id = unique(gene_cluster$gene_id), junction_category = "unknown",
           trans_id = "unknown", junc_id = paste0("unk.", row_number())) %>%
    dplyr::select(chr = chrom, intron_starts = start, intron_ends = end, everything()) %>%
    bind_rows(gene_cluster) %>%
    dplyr::mutate(junc.counts = ifelse(is.na(cluster_idx), 0, junc.counts), cluster.counts = ifelse(is.na(cluster_idx), 0, cluster.counts)) %>%
    arrange(cluster_idx)

  n = all_gene_clusters %>% distinct(cluster_idx, junc_id) %>% dplyr::group_by(cluster_idx) %>%
    dplyr::mutate(junc.per.cluster = ifelse(is.na(cluster_idx), 1, n())) %>% dplyr::ungroup()

  all_gene_clusters %<>% left_join(n, by = c("cluster_idx", "junc_id")) %>%
    mutate(junc.usage = ifelse(junc.counts == 0, 0, round((junc.counts/cluster.counts)*100, digits = 0)), cell_line = paste0(cell_type))

  # create list of expressed isoforms
  # the logic is that if an isoform has a fully unique junction with no counts,
  # it is unlikely to be expressed
  # I want to be careful here since some junctions won't be represented because
  # of technical issues so only eliminate isoforms where none of their unique
  # junctions have counts, avoid filtering on partially unique for now

  non_expressed_genes = all_gene_clusters %>%
    separate_longer_delim(transcript_isoforms, delim = ",") %>% dplyr::group_by(transcript_isoforms, junction_category) %>%
    mutate(isoform_counts = ifelse(junction_category == "fully_unique", sum(junc.counts), junc.counts)) %>%
    filter(junction_category == "fully_unique" & isoform_counts == 0) %>% distinct(transcript_isoforms) %>% ungroup()

  uniquely_targetable = all_gene_clusters %>%
    separate_longer_delim(transcript_isoforms, delim = ",") %>% filter(junction_category == "fully_unique") %>%
    distinct(transcript_isoforms)

  partially_targetable = all_gene_clusters %>%
    separate_longer_delim(transcript_isoforms, delim = ",") %>% filter(junction_category == "partial_unique") %>%
    distinct(transcript_isoforms)

  output = all_gene_clusters %>%
    separate_longer_delim(transcript_isoforms, delim = ",") %>%
    mutate(transcript_targetable = ifelse(transcript_isoforms %in% uniquely_targetable$transcript_isoforms, "Uniquely",
                                          ifelse(transcript_isoforms %in% partially_targetable$transcript_isoforms, "Partially", "None"))) %>%
    mutate(isoform_expressed = ifelse(transcript_isoforms %in% non_expressed_genes$transcript_isoforms, "Unlikely", "Likely"))

  output$output_id = NULL
  output = unique(output)

  if(nrow(output) == 0){
    print(paste0("No junction reads detected for ", name, " in ", cell_type))
  }

  return(output)
}
