#' Helper function to plot cluster sizes
#' plot_cluster_sizes
#' This function takes a list of cluster sizes and plots a bar chart of the sizes.
#'
#' @param juncs A list of integers representing the sizes of each cluster
#'
#' @return A plot
#'
#' @examples
#' plot_cluster_sizes(juncs)
#'

plot_cluster_sizes = function(juncs) {
  cluster_sizes = juncs %>% group_by(cluster_idx) %>% dplyr::summarize(n = n()) %>% ungroup()
  ta = table(cluster_sizes$n)
  tibble(
    cluster_size = as.numeric(names(ta)),
    num_clusters = as.numeric(ta)) %>%
    ggplot(aes(cluster_size, num_clusters)) + geom_point() + xlab("cluster size") + ylab("number clusters") + scale_y_log10()
}

#' Helper function to get junction usage ratios
#' calculate_usage_ratios
#' This function takes
#'
#' @param juncs A list of integers representing the sizes of each cluster
#'
#' @return A plot
#'
#' @examples
#' calculate_usage_ratios(juncs)
#'

calculate_usage_ratios = function(juncs) {
    juncs %>% group_by(cluster_idx) %>%
    mutate(usage_ratio = readcount / sum(readcount)) %>%
    ungroup() # calculate usage ratios
}

#' leafcutter_one_step function
#' A function for obtaining leafcutter intron clusters
#'
#' Inputs junction coordinates from regtools
#' @param juncs_file Gencode gtf or psl file
#' @return intron clusters
#' @examples
#' leafcutter_one_step()
#' @name leafcutter_one_step
#' @export

leafcutter_one_step = function(juncs) {

  juncs = juncs %>% dplyr::select("chrom", "strand", "start", "end", "name", "readcount")
  splice_sites = bind_rows(
    juncs %>% dplyr::select(chrom, strand, position = start),
    juncs %>% dplyr::select(chrom, strand, position = end)) %>%
    distinct() %>%
    arrange(chrom, strand, position) %>%
    dplyr::mutate(idx = 1:n())

  juncs = juncs %>%
    left_join(splice_sites, by = c(chrom = "chrom", strand = "strand", start = "position")) %>%
    left_join(splice_sites, by = c(chrom = "chrom", strand = "strand", end = "position"),
              suffix = c("_start","_end"))

  nss = nrow(splice_sites)

  intron_connectivity <- sparseMatrix(
    i = juncs$idx_start,
    j = juncs$idx_end,
    dims = c(nss, nss),
    x = 1,
    symmetric = TRUE
  )

  g = graph_from_adjacency_matrix(intron_connectivity, "undirected")
  juncs$cluster_idx = components(g)$membership[juncs$idx_start]
  return(juncs)
}

#' isoviz_minicutter function
#' A function for obtaining leafcutter intron clusters
#'
#' Inputs junction coordinates from regtools
#' @param juncs_file Gencode gtf or psl file
#' @return intron clusters
#' @examples
#' isoviz_minicutter()
#' @name isoviz_minicutter
#' @export

isoviz_minicutter = function(juncs_file, min_junction_reads=1, plot_summary=TRUE, min_usage_ratio=0.001) {

  juncs = fread(juncs_file)

  # add column names via https://regtools.readthedocs.io/en/latest/commands/junctions-extract/
  colnames(juncs) = c("chrom", "start", "end", "junc.name", "score", "strand", "thickStart", "thickEnd",
                      "itemRgb", "blockCount", "blockSizes", "blockStarts")

  juncs = juncs %>%
    filter(score >= min_junction_reads) %>%
    separate(blockSizes, into = c("five.p", "three.p"), sep = ",") %>%
    mutate(five.p = as.integer(five.p), three.p = as.integer(three.p)) %>%
    mutate(start = start + five.p, end = end - three.p + 1) %>%
    dplyr::select(chrom = chrom, start, end, name = junc.name, readcount = score, strand = strand)

  juncs = leafcutter_one_step(juncs)

  if(plot_summary){
    print("Printing summary of intron-cluster sizes = how many junctions across in clusters")
    p = plot_cluster_sizes(juncs)
    print(p)
  }

  juncs = juncs %>% calculate_usage_ratios()

  # refine clusters (multiple rounds of refinement don't seem to make much difference)
  juncs_filtered = as.data.table(juncs %>% filter(usage_ratio >= min_usage_ratio))
  juncs_recluster = leafcutter_one_step(juncs_filtered) # multiple rounds of refinement don't seem to make much difference

  if(plot_summary){
    print("Printing summary of intron-cluster sizes = how many junctions across in clusters")
    p=plot_cluster_sizes(juncs_filtered) + ggtitle("Post cluster refinemenet and removing lowly used junctions")
    print(p)
  }

  return(juncs)
}
