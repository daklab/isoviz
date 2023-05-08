#' A function for obtaining leafcutter intron clusters
#'
#' Inputs junction coordinates from regtools
#' @param juncs Gencode gtf or psl file
#' @return intron clusters
#' @examples
#' minicutter()
#' @name minicutter
#' @import tidyverse
#' @import Matrix
#' @import igraph
#' @export

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
  return(components)
}
