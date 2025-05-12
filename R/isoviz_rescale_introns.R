#' A function for reducing the length of introns when visualizing isoforms
#' through re-scaling all coordinates in the gene.
#'
#' Re-scaling intron coordinates while perserving exon alignment by only
#' impacting the intron width and not the exons. These re-scaled coordinates are
#' meant only for plotting!
#'
#' @param introns #description
#' @param exons #description
#' @param width_rescale #description
#' @param transcripts_group Column name containing transcript name
#' @return 'data.table()' containing the re-scaled coordinates of 'introns'
#' and 'exons' for all transcripts in input
#' @examples
#' isoviz_rescale_introns()
#' @name isoviz_rescale_introns
#' @export

#library(data.table)
#library(plyr)
#library(GenomicFeatures)

#introns = rbfox2_introns
#exons = rbfox2_exons
#exons = gene_data # has filtered isoforms
#introns = gene_introns # has all isoforms

isoviz_rescale_introns = function(introns, exons,
                                  width_rescale=50) {
  
  # filter intron isoform data based on exon data
  introns %<>% filter(transcript_name %in% exons$transcript_name)
  
  # very important to consider strand!
  strand = introns$strand[1]

  # order transcripts to figure out which is the reference/starting transcript
  start_pos = unique(exons[,c("start", "transcript_name")])[order(start)][1]
  starts = unique(exons[,c("start", "end", "transcript_name", "transcript_length")])[order(start)]

  # print which transcript gets chosen as the reference 

  if(strand == "-"){
    new_ends = introns$intron_starts
    new_starts = introns$intron_ends
    introns$intron_starts = new_starts
    introns$intron_ends = new_ends
    new_ends_ex = exons$blockstarts
    new_starts_ex = exons$blockends
    exons$blockstarts = new_starts_ex
    exons$blockends = new_ends_ex
    new_ends = exons$start
    new_starts = exons$end
    exons$start = new_starts
    exons$end = new_ends
    starts = unique(exons[,c("start", "end", "transcript_name", "transcript_length")])[order(-start)]
    start_pos = unique(exons[,c("start", "end", "transcript_name")])[order(-start)][1]
  }

  introns$length = abs(introns$intron_starts-introns$intron_ends)
  introns$rescaled_length = introns$length / width_rescale
  transcripts_names = unique(introns$transcript_name)

  if(strand == "-"){
    starts = starts[order(-start)]}

  # Add number of exons per transcript
  exon_counts = exons %>%
    dplyr::group_by(transcript_name) %>%
    dplyr::summarize(num_exons = n())

  # Merge exon counts to start table
  starts = left_join(starts, exon_counts, by = "transcript_name")

  # Count how many transcripts share each start position
start_counts = starts %>%
  dplyr::group_by(start) %>%
  dplyr::summarize(n_transcripts = dplyr::n())

# Join with starts metadata
starts = starts %>%
  dplyr::left_join(start_counts, by = "start")

    z = which(starts$start == starts$start[1])
  starts$new_starts = ""

  # All transcripts that start with start_pos$start will get a new_blockstart set to 0
  z = which(exons$blockstarts == start_pos$start)
  exons$new_blockstart = ""
  exons$new_blockstart[z] = 0 #0 only if the start is the same as gene start

  # convert introns and exons to granges objects
  int_g = introns %>% dplyr::select(chr, intron_starts, intron_ends, strand, transcript_name, rescaled_length)
  ex_g = exons %>% dplyr::select(chr, blockstarts, blockends, strand, blocksizes, transcript_name)
  colnames(int_g)[1:3] = c("chr", "start", "end")
  colnames(ex_g)[1:3] = c("chr", "start", "end")
  z = which(ex_g$start == start_pos$start)
  ex_g$new_start = ""

# Strand-aware choice of most frequent upstream start
if (strand == "+") {
  most_shared_start = starts %>%
    dplyr::arrange(desc(n_transcripts), start) %>%
    dplyr::slice(1) %>%
    dplyr::pull(start)
}

if (strand == "-") {
  most_shared_start = starts %>%
    dplyr::arrange(desc(n_transcripts), desc(start)) %>%
    dplyr::slice(1) %>%
    dplyr::pull(start)
}

# Now filter candidates and choose best scoring reference transcript
ref_transcript = starts %>%
  dplyr::filter(start == most_shared_start) %>%
  dplyr::arrange(desc(num_exons), desc(transcript_length)) %>%
  dplyr::pull(transcript_name) %>%
  .[1]

cat("Ref transcript:", ref_transcript, "\n")
  # Find reference transcript --> longest transcript with blockstart=0
  # Print refernece trancsript that's chosen s

  # Internal function for changing exon coordinates within each transcript
  # according to corresponding introns reduced lengths (update!!!)
  # incorrect to use intron reduced lengths because they are just relative
  # while actual exonic coordinates need to be conserved 

  .get_rescaled_txs = function(trans_name, ref="ref", ref_trans="none") {
    introns_shortened = filter(int_g, transcript_name == trans_name)
    exons_remake = filter(ex_g, transcript_name == trans_name)
    # ensure that everything is ordered correctly
    if (strand == "-"){
      introns_shortened = introns_shortened[order(-start)]
      exons_remake = exons_remake[order(-start)]}
    if (strand == "+"){
      introns_shortened = introns_shortened[order(start)]
      exons_remake = exons_remake[order(start)]}

    # give new coordinate to first exon for each transcript relative to the start
    # of the most upstream transcript

    if(ref=="ref"){ # for reference transcript only
      exons_remake$new_e_start = starts[starts$transcript_name == trans_name]$start - starts$start[1]
      exons_remake$new_e_start = abs(as.numeric(exons_remake$new_e_start))
      exons_remake$new_e_end =  exons_remake$new_e_start + exons_remake$blocksizes[1]
      print("getting ref coordinates")
      for (i in 2:nrow(exons_remake)){
        # update the other exons
        exons_remake$new_e_start[i] = exons_remake$new_e_end[i-1] + introns_shortened$rescaled_length[i-1]
        exons_remake$new_e_end[i] = exons_remake$new_e_start[i] + exons_remake$blocksizes[i]}

      exons_remake$new_transcript_length = max(exons_remake$new_e_end) - min(exons_remake$new_e_start)
      exons_remake$new_start = min(exons_remake$new_e_start)
      exons_remake$new_end = max(exons_remake$new_e_end)
      exons_remake = exons_remake %>% dplyr::select(chr, new_e_start, new_e_end, strand, new_start, new_end, start, end, transcript_name, blocksizes)
    }

    # for non reference transcripts, figure out which positions in transcript
    # are the same as reference

    else{
      exons_remake$new_e_start <- NA
      exons_remake$new_e_end <- NA

      for (i in 1:nrow(exons_remake)) {
        # Find matching row(s) in ref based on start or end position
        matching_row_start = ref_trans$start == exons_remake$start[i]
        match_row_end = ref_trans$end == exons_remake$end[i]

        # Check if any matching row(s) found
        if (any(matching_row_start)) {
          # Assign corresponding new_e_start and new_e_end values to exons_remake
          exons_remake$new_e_start[i] <- ref_trans$new_e_start[matching_row_start]
        }
        if (any(match_row_end)) {
          # Assign corresponding new_e_start and new_e_end values to exons_remake
          exons_remake$new_e_end[i] <- ref_trans$new_e_end[match_row_end]
        }
      }
    }

    exons_remake$exon_num = 1:nrow(exons_remake)
    return(exons_remake)
  }

  # first get new coords for ref_trans
  ref = .get_rescaled_txs(ref_transcript, ref="ref")
  unref_trans=transcripts_names[-which(transcripts_names == ref_transcript)]
  ref_attach_to = ref
  ref_attach_to$new_end = NULL
  #colnames(ref_attach_to) = c("chr", "start", "end", "strand", "blocksizes",
   #                           "transcript_name" , "new_start","new_e_start",
    #                          "new_e_end", "exon_num")
  ref_attach_to = ref_attach_to[1,]
  ref_attach_to[1,] = 0

  for (i in 1:length(unref_trans)){
    trans = .get_rescaled_txs(unref_trans[i], ref="not_ref", ref)
    ref_attach_to = rbind(ref_attach_to, trans)
  }
  ref_attach_to = ref_attach_to[-1,]
  ref_attach_to$new_e_start = as.numeric(ref_attach_to$new_e_start)
  ref_attach_to$blocksizes = as.numeric(ref_attach_to$blocksizes)
  ref_attach_to$fill_ends = ref_attach_to$new_e_start + ref_attach_to$blocksizes
  ref_attach_to$fill_starts = ref_attach_to$new_e_end - ref_attach_to$blocksizes

  ref_attach_to$new_e_start[is.na(ref_attach_to$new_e_start)] = ref_attach_to$fill_starts[is.na(ref_attach_to$new_e_start)]
  ref_attach_to$new_e_end[is.na(ref_attach_to$new_e_end)] = ref_attach_to$fill_ends[is.na(ref_attach_to$new_e_end)]

  ref_attach_to$fill_starts = NULL
  ref_attach_to$fill_ends = NULL

  # for rows that missing new start and end need to find their nearest exons that do have it
  # and use scaled introns coords to fix it
nas = filter(ref_attach_to, is.na(new_e_start))
if(!(dim(nas)[1] == 0)){

  z = which(is.na(ref_attach_to$new_e_start))
  ref_attach_to = ref_attach_to[-z,]
  nas = nas[order(-exon_num)]

  for(i in 1:nrow(nas)) {
    trans = nas$transcript_name[i]
    exon_start = nas$start[i]
    exon_end = nas$end[i]

    # Look for same start in other transcripts
    matched_start = ref_attach_to %>% filter(start == exon_start & !is.na(new_e_start)) %>% pull(new_e_start)
    matched_end = ref_attach_to %>% filter(end == exon_end & !is.na(new_e_end)) %>% pull(new_e_end)

    if (length(matched_start) > 0) {
      nas$new_e_start[i] = matched_start[1]
      nas$new_e_end[i] = nas$new_e_start[i] + nas$blocksizes[i]
      ref_attach_to = rbind(ref_attach_to, nas[i])
    } else if (length(matched_end) > 0) {
      nas$new_e_end[i] = matched_end[1]
      nas$new_e_start[i] = nas$new_e_end[i] - nas$blocksizes[i]
      ref_attach_to = rbind(ref_attach_to, nas[i])
    } else {
      # Use intron-based position estimation
      introns_shortened = filter(int_g, transcript_name == trans)
      introns_shortened = introns_shortened[order(ifelse(strand == "-", -start, start))]
      introns_shortened$intron_num = 1:nrow(introns_shortened)

      prev_exons = ref_attach_to %>%
        filter(transcript_name == trans) %>%
        arrange(exon_num)

      exon_n = nas$exon_num[i]

      # Case 1: Between exons (use previous exon end + intron)
      if (exon_n > 1 && (exon_n - 1) <= nrow(prev_exons)) {
        prev_end = prev_exons$new_e_end[exon_n - 1]
        intron_len = introns_shortened$rescaled_length[exon_n - 1]
        nas$new_e_start[i] = prev_end + intron_len
        nas$new_e_end[i] = nas$new_e_start[i] + nas$blocksizes[i]
        ref_attach_to = rbind(ref_attach_to, nas[i])

      # Case 2: Before known exon (use next exon start - intron)
      } else if (exon_n <= nrow(introns_shortened)) {
        next_exon = prev_exons$new_e_start[exon_n + 1]
        intron_len = introns_shortened$rescaled_length[exon_n]
        nas$new_e_end[i] = next_exon - intron_len
        nas$new_e_start[i] = nas$new_e_end[i] - nas$blocksizes[i]
        ref_attach_to = rbind(ref_attach_to, nas[i])

         # 🔥 New Case 3: Last exon beyond known range — extend from last exon
        }  else if (exon_n > nrow(prev_exons)) {
    last_exon_end = max(prev_exons$new_e_end, na.rm = TRUE)
    last_intron_len = introns_shortened$rescaled_length[exon_n - 1]
    nas$new_e_start[i] = last_exon_end + last_intron_len
    nas$new_e_end[i] = nas$new_e_start[i] + nas$blocksizes[i]
  }
    ref_attach_to = rbind(ref_attach_to, nas[i])
    }
  }
}
 # (!(dim(nas)[1] == 0))
  ref_attach_to$new_start = NULL
  ref = ref %>% dplyr::select(chr, start, end, strand, blocksizes, transcript_name, new_e_start, new_e_end, exon_num)
  ref_attach_to = rbind(ref_attach_to, ref)

  min_starts = ref_attach_to %>% dplyr::group_by(transcript_name) %>% dplyr::summarize(min_start = min(new_e_start))
  max_exons = ref_attach_to %>% dplyr::group_by(transcript_name) %>% dplyr::summarize(max_end = max(new_e_end))

  trans_lengths = merge(min_starts, max_exons)
  ref_attach_to = merge(ref_attach_to, trans_lengths)
  ref_attach_to$transcript_length = ref_attach_to$max_end - ref_attach_to$min_start   # max exon end and min exon start for each transcript
  return(list(ref_attach_to, int_g))

  } #end of isoviz_rescale_introns