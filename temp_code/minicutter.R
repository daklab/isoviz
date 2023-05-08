# pulled from David's github on 4/20/23
# edited regtools input on 4/21/23

require(tidyverse)
require(Matrix)
require(igraph)

leafcutter_one_step = function(juncs) {

  splice_sites = bind_rows(
    juncs %>% select(chrom, strand, position = start),
    juncs %>% select(chrom, strand, position = end)) %>%
    distinct() %>%
    arrange(chrom, strand, position) %>%
    mutate(idx = 1:n())

  juncs = juncs %>%
    left_join(splice_sites, by = c(chrom = "chrom", strand = "strand", start = "position")) %>%
    left_join(splice_sites, by = c(chrom = "chrom", strand = "strand", end = "position"), suffix = c("_start","_end"))

  nss = nrow(splice_sites)

  intron_connectivity <- sparseMatrix(
    i = juncs$idx_start,
    j = juncs$idx_end,
    dims = c(nss, nss),
    x = 1,
    symmetric = TRUE
  )

  g = graph_from_adjacency_matrix(intron_connectivity, "undirected")
  components(g)$membership[juncs$idx_start]
}

leafcutter_refine = function(juncs, min_usage_ratio = 0.001) {
  juncs_filtered = juncs %>% filter(usage_ratio >= min_usage_ratio)
  juncs_filtered$cluster_idx = leafcutter_one_step(juncs_filtered)
  juncs_filtered
}

calculate_usage_ratios = function(juncs) {
  juncs %>% group_by(cluster_idx) %>% mutate(usage_ratio = readcount / sum(readcount)) %>% ungroup() # calculate usage ratios
}

leafcutter_two_step = function(juncs, min_usage_ratio = 0.001) {
  juncs$cluster_idx = leafcutter_one_step(juncs) # get inital clustering
  juncs = juncs %>% calculate_usage_ratios()
  leafcutter_refine(juncs, min_usage_ratio = min_usage_ratio)
}

plot_cluster_sizes = function(juncs) {
  cluster_sizes = juncs %>% group_by(cluster_idx) %>% summarize(n = n()) %>% ungroup()
  ta = table(cluster_sizes$n)
  tibble(
    cluster_size = as.numeric(names(ta)),
    num_clusters = as.numeric(ta)
  ) %>%
    ggplot(aes(cluster_size, num_clusters)) + geom_point() + xlab("cluster size") + ylab("number clusters") + scale_y_log10()
}

min_junction_reads = 1

# code to get correct coords from regtools output
juncs <- read_tsv(
  "~/knowles_lab/Megan/230320_cas13_validation/junctions/HEK-NFYA-rep2-empty_v41_basic.junc",
  col_names = FALSE, col_types = "cddcdcddndcc") %>%
  filter(X5 >= min_junction_reads) %>%
  separate(X11, into = c("five.p", "three.p"), sep = ",") %>%
  mutate(five.p = as.integer(five.p), three.p = as.integer(three.p)) %>%
  mutate(start = X2 + five.p, end = X3 - three.p + 1) %>%
  select(chrom = X1, start, end, name = X4, readcount = X5, strand = X6)

juncs$cluster_idx = leafcutter_one_step(juncs)
plot_cluster_sizes(juncs)

juncs = juncs %>% calculate_usage_ratios()
juncs_recluster = leafcutter_refine(juncs, min_usage_ratio = 0.01) # multiple rounds of refinement don't seem to make much difference
plot_cluster_sizes(juncs_recluster)

save(juncs_recluster, file = "~/isoviz/data/HEK293_all_lfc_junctions.rda")
