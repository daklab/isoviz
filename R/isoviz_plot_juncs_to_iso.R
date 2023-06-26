#' A function for obtaining a comprehensive visualization of gene transcript
#' isoforms with junction coverage across intron clusters information from RNA-seq
#' and leafcutter clustering
#' @param mapped_junctions output from running isoviz_map_junctions
#' @param cell_type cell type of interest, default=hESC
#' @param junction_usage minimum junction usage as a percentage, default is 5
#' @param include_all_juncs show all junctions on plot, will still filter out based on junction usage parameter, default is TRUE
#' @param include_specific_junctions provide list of junctions to plot in second panel, can only specify if include_all_juncs is FALSE
#' @return plot
#' @examples
#' isoviz_plot_juncs_to_iso()
#' @name isoviz_plot_juncs_to_iso
#' @import data.table
#' @import magrittr
#' @import ggsci
#' @import cowplot
#' @export

library(data.table)
library(dplyr)
library(tidyr)
library(magrittr)
library(ggsci)
library(cowplot)

#gene_data=rbfox2_exons

isoviz_plot_juncs_to_iso = function(mapped_junctions, gene_data,
                                    gene_introns,
                                    cell_type = "hESC",
                                    junc_usage = 5,
                                    intron_scale = "no",
                                    intron_scale_width = 10,
                                    include_all_juncs = TRUE,
                                    include_specific_junctions = c()){

  # call isoviz_map_junctions first and use output from that
  df = mapped_junctions

  # filter here if the user wants to filter for specific junctions
  if(include_all_juncs == FALSE){
    if(length(include_specific_junctions) == 0){
      stop("Must include list of junctions if include_all_juncs == FALSE")
    } else{
      df %<>% filter(junc_id %in% include_specific_junctions)
    }
  }

  # filter based on junction usage
  df = dplyr::filter(df, junc.usage >= junc_usage)

  # setting parameters for the plot based on length of the gene and number of isoforms
  length = max(gene_data$end) - min(gene_data$start)
  min_start = min(gene_data$start)
  max_end = max(gene_data$end)
  n_transcripts = gene_data %>% group_by(trans_id) %>% tally() %>% nrow()
  y_max = n_transcripts + 1
  strand = unique(gene_data$strand)

  if(intron_scale == "yes"){
    # also remove NA isoform from mapped junctions if we want to scale
    df = dplyr::filter(df, !(is.na(transcript_isoforms)))
    old_exons = gene_data
    # rescaled gene_data
    scaled_res = isoviz_rescale_introns(gene_introns, gene_data, intron_scale_width)
    gene_data = scaled_res[[1]]
    int_g = scaled_res[[2]]
    gene_data$blockstarts = gene_data$new_e_start
    gene_data$blockends = gene_data$new_e_end

    # save old exon starts and ends
    gene_data$old_e_start = gene_data$start
    gene_data$old_e_end = gene_data$end

    gene_data$start = gene_data$min_start  # transcript start
    gene_data$end =  gene_data$max_end  # transcript end

    min_start = min(gene_data$start)
    max_end = max(gene_data$end)
  }

  # order the isoforms by length for plotting
  order_plot = gene_data %>% group_by(transcript_name, transcript_length) %>% tally() %>%
    dplyr::arrange(desc(transcript_length), desc(n))
  order_plot$trans_order = 1:nrow(order_plot)
  order_plot = order_plot %>% dplyr::select(transcript_name, trans_order)

  to_plot_ordered = gene_data %>% full_join(order_plot, by = c("transcript_name")) %>%
    arrange(trans_order, blockstarts) %>% mutate(seg_end = end)

  # make isoform - level plot
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
    ggtitle(paste0(cell_type, ": ", df$gene_name[1], " Junction to Isoform Map (" , to_plot_ordered$strand[1], " strand )")) +
    coord_cartesian(xlim = c(min_start, max_end), clip = 'off')

  # code to plot leafcutter cluster info underneath isoform level data
  # aggregate isoform information
  df = unique(df)

  text = df %>% tidyr::separate(col = transcript_name, into = c("name", "transcript")) %>%
    mutate(transcript = ifelse(is.na(transcript), "unknown", transcript)) %>%
    group_by(cluster_idx, junc_id) %>%
    arrange(cluster_idx, junc_id, transcript) %>%
    dplyr::mutate(isoforms = toString(unique(transcript))) %>%
    dplyr::select(1:5, junc_id, junction_category, junc.usage, isoforms, junc.per.cluster) %>% dplyr::distinct() %>% ungroup() %>%
    dplyr::mutate(text_plot = paste0(junc.counts, " reads (", junc.usage, "%) : ", isoforms), alpha = junc.usage/100)

  text = unique(text)
  clusters = unique(df$cluster_idx)

  introns = df %>% dplyr::select(-gene_id, -gene_name, -transcript_isoforms, -transcript_targetable, -isoform_expressed) %>%
    distinct() %>% arrange(cluster_idx, junc.usage)
  introns %<>% left_join(text)

  # get colors
  col_n = nrow(introns %>% distinct(cluster_idx))
  mycols <- rep(pal_npg("nrc", alpha = 1)(8), length.out = col_n)

  # get updated coordinates based on plot above
  if(intron_scale == "yes"){
    #fix intron coordinates
    #if(strand == "-"){
    #  colnames(int_g)[2:3] = c("intron_ends", "intron_starts")} #is this always the case?
    #if(strand == "+"){
    #  colnames(int_g)[2:3] = c("intron_starts", "intron_ends")}
    #}
    #introns = merge(introns, int_g)
    if(strand == "-"){
      introns$new_intron_start = introns$intron_ends
      introns$new_intron_end = introns$intron_starts
      introns$intron_ends = introns$new_intron_end
      introns$intron_starts = introns$new_intron_start
      introns$new_intron_start = NULL
      introns$new_intron_end = NULL
    }

    introns$length = abs(introns$intron_ends - introns$intron_starts)
    introns$rescaled_length = introns$length / intron_scale_width

    # merge with exon data to map new start and end points
    if(strand == "-"){
      get_coords = unique(to_plot_ordered[,c("new_e_start", "old_e_start", "transcript_name")])
      colnames(get_coords)[2] = "new_intron_end"
      introns$new_intron_end = introns$intron_ends
      get_coords = unique(merge(introns, get_coords, by=c("transcript_name", "new_intron_end")))

      # for each junction just want one set of coordinates for start and end
      # so new_intron_end is old_exon_start
      get_coords$new_intron_end = get_coords$new_e_start
      get_coords$new_intron_start = get_coords$new_intron_end - get_coords$rescaled_length}

    if(strand == "+"){
        get_coords = unique(to_plot_ordered[,c("new_e_end", "old_e_end", "transcript_name")])
        colnames(get_coords)[2] = "new_intron_start"
        introns$new_intron_start = introns$intron_starts
        get_coords = unique(merge(introns, get_coords, by=c("transcript_name", "new_intron_start")))

        # for each junction just want one set of coordinates for start and end
        get_coords$new_intron_start = get_coords$new_e_end
        get_coords$new_intron_end = get_coords$new_intron_start + get_coords$rescaled_length}

    introns = get_coords
    introns$intron_starts = introns$new_intron_start
    introns$intron_ends = introns$new_intron_end
    }

  # get unique entries only
  introns = unique(introns %>% dplyr::select("chr", "intron_starts", "intron_ends",
                            "cluster_idx", "text_plot", "junc_id"))

  # reorder
  introns = as.data.table(introns)
  introns = introns[order(cluster_idx)] #, junc.usage)]
  introns$junc.order = 1:nrow(introns)

  if(intron_scale == "yes"){
    p2 = ggplot() +
      geom_rect(data = introns, mapping = aes(xmin = intron_starts, xmax = intron_ends, ymin = junc.order - 0.2, ymax = junc.order + 0.2, fill = as.character(cluster_idx))) +
      scale_fill_manual(values = mycols) +
      xlim(min_start, max_end) + theme_bw() +
      xlab("Re-scaled Genomic Position (bp)") + ylab("") +
      geom_text(data = introns, aes(x = Inf, y = junc.order, hjust = -0.1, label = junc_id), size = 3) +
      theme(axis.text.y=element_blank(), axis.ticks.y=element_blank(), legend.position = "none", plot.margin = margin(0,1,0,0.1, "in")) +
      coord_cartesian(xlim = c(min_start, max_end), clip = 'off')

  }
  if(intron_scale == "no"){
  p2 = ggplot() +
    geom_rect(data = introns, mapping = aes(xmin = intron_starts, xmax = intron_ends, ymin = junc.order - 0.2, ymax = junc.order + 0.2, fill = as.character(cluster_idx))) +
    geom_text(data = introns, aes(x = intron_ends, y = junc.order, label = text_plot), hjust = "inward", size = 3) +
    scale_fill_manual(values = mycols) +
    xlim(min_start, max_end) + theme_bw() +
    xlab("Hg38 Genomic Position (bp)") + ylab("") +
    geom_text(data = introns, aes(x = Inf, y = junc.order, hjust = -0.1, label = junc_id), size = 3) +
    theme(axis.text.y=element_blank(), axis.ticks.y=element_blank(), legend.position = "none", plot.margin = margin(0,1,0,0.1, "in")) +
    coord_cartesian(xlim = c(min_start, max_end), clip = 'off')}

  # plot dimenstions based on number of isoforms and junctions
  total_height = (n_transcripts + nrow(introns))

  # p1
  if(n_transcripts/total_height <= 0.3){
    iso_h = n_transcripts*1.5
  }else{
    iso_h = n_transcripts
  }

  # p2
  if(n_transcripts/total_height <= 0.3){
    junc_h = nrow(introns) * 0.8
  }else{
    junc_h = nrow(introns) * 1.2
  }
  p1_p2 = plot_grid(p1, p2, ncol = 1, align = "v", rel_heights = c(n_transcripts/total_height, nrow(introns)/total_height))

  return(p1_p2)
}
