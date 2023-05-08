#' A function for obtaining a comprehensive visualization of gene transcript
#' isoforms with junction coverage across intron clusters information from RNA-seq
#' and leafcutter clustering
#' @param gene ensembl ID for now
#' @param cell_type cell type of interest, default=hESC
#' @param junction_usage minimum junction usage, default is 5
#' @param include_common_juncs show junctions common to all isoforms?
#' @param filter_expressed_isoforms include only isoforms that are detected
#' @param filter_all_100percent filter junctions that are always included (is this right?)
#' @return plot
#' @examples
#' isoviz_plot_juncs_to_iso()
#' @name isoviz_plot_juncs_to_iso
#' @import data.table
#' @export

isoviz_plot_juncs_to_iso = function(gene = "ENSG00000099875", cell_type = "hESC", junction_usage = 5, include_common_juncs = FALSE, filter_expressed_isoforms = TRUE, filter_all_100percent = FALSE){

  # call map_junction function here
  df = map_all_junctions(gene, cell_type)

  # filter here based on input parameters- also unknown junctions that belong to common junctions
  if(include_common_juncs == FALSE){
    remove_cluster = df %>% filter(junction_category == "common") %>% select(cluster_idx) %>% distinct()
    df %<>% filter(!cluster_idx %in% remove_cluster$cluster_idx)
  }

  if(filter_expressed_isoforms == TRUE){
    df %<>% filter(isoform_expressed == "Likely")

    # create list of expressed transcripts
    expressed_isoforms = df %>% filter(isoform_expressed == "Likely", transcript_isoforms != "unknown") %>%
      separate(col = transcript_isoforms, into = c("name", "transcript")) %>% select(transcript) %>% distinct() %>%
      arrange(transcript) %>%
      mutate(isoforms = paste0(transcript, collapse = ",")) %>% select(isoforms) %>% distinct()
  }

  to_remove = df %>% filter(junc.per.cluster == 2 & (junc.usage < junction_usage | junc.usage > (100 - junction_usage)))
  df %<>% filter(!junc_id %in% to_remove$junc_id, junc.counts > 0)

  # get info after filtering
  name = unique(df$gene_name)
  transcripts = unique(df$transcript_isoforms)

  # filter transcripts for plotting too
  gene_data = iso_data_complete %>% filter(transcript_name %in% df$transcript_isoforms)

  # setting parameters for the plot based on length of the gene and number of isoforms
  length = max(gene_data$end) - min(gene_data$start)
  min_start = min(gene_data$start)
  n = gene_data %>% group_by(trans_id) %>% tally() %>% nrow()
  y_max = n + 1
  max_end = max(gene_data$end)

  # order the isoforms by length for plotting
  order_plot = gene_data %>% group_by(transcript_name, transcript_length) %>% tally() %>%
    dplyr::arrange(desc(transcript_length), desc(n))
  order_plot$trans_order = 1:nrow(order_plot)
  order_plot = order_plot %>% dplyr::select(transcript_name, trans_order)

  to_plot_ordered = gene_data %>% full_join(order_plot, by = c("transcript_name")) %>%
    arrange(trans_order, blockstarts)
  to_plot_ordered = to_plot_ordered %>% mutate(seg_end = end)

  # make isoform plot
  p1 = ggplot() +
    geom_segment(aes(x = to_plot_ordered$start, y = to_plot_ordered$trans_order,
                     xend = to_plot_ordered$seg_end, yend = to_plot_ordered$trans_order)) +
    geom_rect(data = to_plot_ordered, mapping = aes(xmin = blockstarts, xmax = blockends, ymin = trans_order - 0.3, ymax = trans_order + 0.3)) +
    xlim(min_start, max_end) + ylim(0,y_max) +
    geom_text(aes(x = Inf, y = to_plot_ordered$trans_order, hjust = -0.1, label = to_plot_ordered$transcript_name), size = 3, check_overlap = TRUE) +
    theme_bw() + ylab("") + xlab("") +
    theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(),
          axis.text.y=element_blank(), axis.ticks.y=element_blank(), legend.position = "top",
          plot.title = element_text(hjust = 0.5), plot.margin = margin(0.1,1,0,0.1, "in")) +
    ggtitle(paste0(cell_type, ": ", name, " Junction to Isoform Map (" , to_plot_ordered$strand[1], " strand )")) +
    coord_cartesian(xlim = c(min_start, max_end), clip = 'off')

  # new code to plot leafcutter cluster info underneath isoform level data
  # aggregate isoform information
  text = df %>% separate(col = transcript_isoforms, into = c("name", "transcript")) %>%
    mutate(transcript = ifelse(is.na(transcript), "unknown", transcript)) %>%
    group_by(cluster_idx, junc_id) %>%
    arrange(cluster_idx, junc_id, transcript) %>%
    mutate(isoforms = paste0(transcript, collapse = ",")) %>%
    select(1:5, junc_id, junction_category, junc.usage, isoforms, junc.per.cluster) %>% distinct() %>% ungroup() %>%
    mutate(text_plot = paste0(junc.counts, " reads (", junc.usage, "%) : ", isoforms), alpha = junc.usage/100)

  clusters = unique(df$cluster_idx)
  introns = df %>% select(-gene_id, -gene_name, -transcript_isoforms) %>%
    distinct() %>% arrange(cluster_idx, junc.usage)
  introns %<>% left_join(text)

  # gets rid of gencode partially unique junctions that are common in this cell type (100% usage)
  if(filter_expressed_isoforms == TRUE & filter_all_100percent == FALSE){
    remove_junc = introns %>% filter(junc.usage == "100" & isoforms %in% expressed_isoforms$isoforms)
    introns %<>% filter(!junc_id %in% remove_junc$junc_id)
  }else{
    introns %<>% filter(junc.usage != "100")
  }

  introns$junc.order = 1:nrow(introns)

  # for guides - remove this or print out?
  for_guides = introns %>% filter(grepl("junc", junc_id))

  p2 = ggplot() +
    geom_rect(data = introns, mapping = aes(xmin = junc_start, xmax = junc_end, ymin = junc.order - 0.2, ymax = junc.order + 0.2, fill = as.character(cluster_idx), alpha = alpha)) +
    geom_text(data = introns, aes(x = junc_end, y = junc.order, label = text_plot), hjust = "inward", size = 3) +
    scale_fill_npg() +
    xlim(min_start, max_end) + theme_bw() +
    xlab("Hg38 Genomic Position (bp)") + ylab("") +
    geom_text(data = introns, aes(x = Inf, y = junc.order, hjust = -0.1, label = junc_id), size = 3) +
    theme(axis.text.y=element_blank(), axis.ticks.y=element_blank(), legend.position = "none", plot.margin = margin(0,1,0,0.1, "in")) +
    coord_cartesian(xlim = c(min_start, max_end), clip = 'off')

  # plot dimenstions based on number of isoforms and junctions
  t = (n + nrow(introns))
  h = t/3

  # p1
  if(n <= 4){
    iso_h = n*1.3
  }else{
    iso_h = n*0.8
  }

  # p2
  if(nrow(introns) > 6){
    junc_h = nrow(introns) * 0.8
  }else{
    junc_h = nrow(introns) * 1.2
  }
  p1_p2 = plot_grid(p1, p2, ncol = 1, align = "v", rel_heights = c(iso_h, junc_h))

  return(p1_p2)
}
