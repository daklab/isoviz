gene_data = rescaled_coords

# setting parameters for the plot based on length of the gene and number of isoforms
length = max(gene_data$transcript_length) 
min_start = min(gene_data$new_e_start)
#min_start=0
n = gene_data %>% group_by(transcript_name) %>% tally() %>% nrow()
y_max = n + 1
max_end = max(gene_data$new_e_end)

# order the isoforms by length for plotting
order_plot = gene_data %>% group_by(transcript_name, transcript_length) %>% tally() %>%
  dplyr::arrange(desc(transcript_length), desc(n))
order_plot$trans_order = 1:nrow(order_plot)
order_plot = order_plot %>% dplyr::select(transcript_name, trans_order)

to_plot_ordered = gene_data %>% full_join(order_plot, by = c("transcript_name")) %>%
  arrange(trans_order, new_e_start) %>%
  mutate(y1 = trans_order - 0.2, y2 = trans_order + 0.2)
to_plot_ordered = to_plot_ordered %>% mutate(seg_end = new_e_end)
  
# make plot
p1 = ggplot() +
    geom_segment(aes(x = to_plot_ordered$new_e_start, y = to_plot_ordered$trans_order,
                     xend = to_plot_ordered$new_e_end, yend = to_plot_ordered$trans_order)) +
    geom_rect(data = to_plot_ordered, mapping = aes(xmin = new_e_start, xmax = new_e_end, ymin = y1, ymax = y2)) +
    xlim(min_start, max_end) + ylim(0,y_max) +
    geom_text(aes(x = Inf, y = to_plot_ordered$trans_order, hjust = 0, label = to_plot_ordered$transcript_name)) +
    theme_bw() + ylab("") +
    theme(axis.text.y=element_blank(), axis.ticks.y=element_blank(), legend.position = "top") +
    xlab("Re-scaled coordinates") +
    ggtitle(paste0(" (", to_plot_ordered$strand[1], ")")) +
    theme(axis.text = element_text(size = 10), plot.title = element_text(hjust = 0.5)) +
    coord_cartesian(clip = 'off') +
    theme(plot.margin = margin(0.5,1.5,0.5,0.5, "in"))

print(p1)
