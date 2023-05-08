#' A function for mapping junctions to transcript isoforms
#' @param gene ensembl ID for now
#' @param cell_type cell type of interest, default=hESC
#' @return plot
#' @examples
#' isoviz_plot_juncs_to_iso()
#' @name isoviz_plot_juncs_to_iso
#' @import data.table
#' @export

# map junctions function
map_all_junctions = function(gene = "ENSG00000100320", cell_type = "hESC", intron_coords, leafcutter_input){

  # join intron_coords with leafcutter data and mark any unannotated junctions

  # since leafcutter junctions dont have gene annotations, need to combine here first
  juncs_recluster %<>% select(everything(), junc.counts = readcount) %>%
    group_by(cluster_idx) %>% mutate(cluster.counts = sum(junc.counts)) %>%
    ungroup()

  gene_cluster = gene_iso_data %>%
    left_join(juncs_recluster, by = c("chr" = "chrom", "junc_start" = "start", "junc_end" = "end", "strand")) %>%
    arrange(cluster_idx) # 24

  # check minicutter output
  juncs_recluster %>% filter(cluster_idx == "134447")

  name = unique(gene_cluster$gene_name)
  gene_strand = unique(gene_cluster$strand)

  # to get information on leafcutter junctions that don't match annotation
  all_gene_clusters = juncs_recluster %>% filter(cluster_idx %in% gene_cluster$cluster_idx) %>%
    left_join(gene_iso_data, by = c("chrom" = "chr", "start" = "junc_start", "end" = "junc_end", "strand")) %>%
    filter(is.na(gene_id)) %>%
    mutate(gene_name = unique(gene_cluster$gene_name), strand = unique(gene_cluster$strand),
           gene_id = unique(gene_cluster$gene_id), junction_category = "unknown",
           transcript_isoforms = "unknown", junc_id = paste0("unk.", row_number())) %>%
    select(chr = chrom, junc_start = start, junc_end = end, everything()) %>%
    bind_rows(gene_cluster) %>%
    mutate(junc.counts = ifelse(is.na(cluster_idx), 0, junc.counts), cluster.counts = ifelse(is.na(cluster_idx), 0, cluster.counts)) %>%
    arrange(cluster_idx)

  n = all_gene_clusters %>% distinct(cluster_idx, junc_id) %>% group_by(cluster_idx) %>%
    mutate(junc.per.cluster = ifelse(is.na(cluster_idx), 1, n())) %>% ungroup()

  all_gene_clusters %<>% left_join(n, by = c("cluster_idx", "junc_id")) %>%
    mutate(junc.usage = ifelse(junc.counts == 0, 0, round((junc.counts/cluster.counts)*100, digits = 0)), cell_line = paste0(cell_type))

  #test
  all_gene_clusters %>% filter(grepl("209", transcript_isoforms)) %>% arrange(junc_start)

  # create list of expressed isoforms
  # the logic is that if an isoform has a fully unique junction with no counts, it is unlikely to be expressed
  # I want to be careful here since some junctions won't be represented because of technical issues so only eliminate isoforms where none of their unique junctions have counts, avoid filtering on partially unique for now
  non_expressed_genes = all_gene_clusters %>%
    separate_longer_delim(transcript_isoforms, delim = ",") %>% group_by(transcript_isoforms, junction_category) %>%
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
    mutate(isoform_expressed = ifelse(transcript_isoforms %in% non_expressed_genes$transcript_isoforms, "Unlikely", "Likely")) # 135

  if(nrow(output) == 0){
    print(paste0("No junction reads detected for ", name, " in ", cell_type))
  }

  return(output)
}
