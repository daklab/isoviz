#' A function for plotting transcript isoforms for a specific gene
#'
#' This function allows you to visualize all annotated isoforms for a given gene.
#' @param gene Input your gene of interest using ENSG for now. Defaults to "ENSG00000100320" (RBFOX2)
#' @return A plot.
#' @examples
#' plot_isoforms("ENSG00000100320")
#' @name plot_isoforms

load(file='data/iso_exon_data.rda')

plot_isoforms = function(gene="ENSG00000100320"){

  gene_trans =which(str_detect(genomedata$gene_id, gene))
  gene_data = genomedata[gene_trans,]

  # setting parameters for the plot based on length of the gene and number of isoforms
  length = max(gene_data$end) - min(gene_data$start)
  min_start = min(gene_data$start)
  n = gene_data %>% group_by(trans_id) %>% tally() %>% nrow()
  y_max = n + 1
  max_end = max(gene_data$end)

  # order the isoforms by length for plotting
  order_plot = gene_data %>% group_by(trans_id, transcript_length) %>% tally() %>%
    dplyr::arrange(desc(transcript_length), desc(n))
  order_plot$trans_order = 1:nrow(order_plot)
  order_plot = order_plot %>% dplyr::select(trans_id, trans_order)

  to_plot_ordered = gene_data %>% full_join(order_plot, by = c("trans_id")) %>%
    arrange(trans_order, blockstarts) %>%
    mutate(y1 = trans_order - 0.2, y2 = trans_order + 0.2)
  to_plot_ordered = to_plot_ordered %>% mutate(seg_end = end)

  # make plot
  p1 = ggplot() +
    geom_segment(aes(x = to_plot_ordered$start, y = to_plot_ordered$trans_order,
                     xend = to_plot_ordered$seg_end, yend = to_plot_ordered$trans_order)) +
    geom_rect(data = to_plot_ordered, mapping = aes(xmin = blockstarts, xmax = blockends, ymin = y1, ymax = y2)) +
    xlim(min_start, max_end) + ylim(0,y_max) +
    geom_text(aes(x = Inf, y = to_plot_ordered$trans_order, hjust = 0, label = to_plot_ordered$trans_id)) +
    theme_bw() + ylab("") +
    theme(axis.text.y=element_blank(), axis.ticks.y=element_blank(), legend.position = "top") +
    xlab("Hg38 Genomic Position (bp)") +
    ggtitle(paste0(gene, " (", to_plot_ordered$strand[1], ")")) +
    theme(axis.text = element_text(size = 10), plot.title = element_text(hjust = 0.5)) +
    coord_cartesian(clip = 'off') +
    theme(plot.margin = margin(0.5,1.5,0.5,0.5, "in"))
  print(p1)
}
