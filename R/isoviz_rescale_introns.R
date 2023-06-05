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
#' @import data.table
#' @import tidyr
#' @import dplyr
#' @import readr
#' @export

library(data.table)
library(plyr)
library(GenomicFeatures)

#introns = rbfox2_introns
#exons = rbfox2_exons

isoviz_rescale_introns = function(introns, exons,
                          width_rescale=50) {

  # very important to consider strand!
  strand = introns$strand[1]
  print(strand)

  # order transcripts to figure out which is the reference/starting transcript
  start_pos = unique(exons[,c("start", "transcript_name")])[order(start)][1]
  starts = unique(exons[,c("start", "end", "transcript_name", "transcript_length")])[order(start)]

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

  z = which(starts$start == starts$start[1])
  starts$new_starts = ""
  starts$new_starts[z] = 0 #0 only if the start is the same as gene start

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
  ex_g$new_start[z] = 0

  # Find reference transcript --> longest transcript with blockstart=0
  if(strand == "-"){
    ref_transcript = filter(starts, new_starts == 0)[order(end)]$transcript_name[1] #if strand is - need to change if strand is +
  }
  if(strand == "+"){
    ref_transcript = filter(starts, new_starts == 0)[order(-end)]$transcript_name[1]
  }

  # Internal function for changing exon coordinates within each transcript
  # according to corresponding introns reduced lengths

  .get_rescaled_txs = function(trans_name, ref="ref", ref_trans="none") {
    introns_shortened = filter(int_g, transcript_name == trans_name)
    exons_remake = filter(ex_g, transcript_name == trans_name)
    print(trans_name)
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
  colnames(ref_attach_to) = c("chr", "start", "end", "strand", "blocksizes",
                              "transcript_name" , "new_start","new_e_start",
                              "new_e_end", "exon_num")
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

  for(i in 1:nrow(nas)){
    print(i)
    trans = nas$transcript_name[i]
    exons = filter(ref_attach_to, transcript_name == trans) #all the non NA exons
    exon_n = nas$exon_num[i]
    introns_shortened = filter(int_g, transcript_name == trans)
    if (strand == "-"){
      introns_shortened = introns_shortened[order(-start)]
      exons = exons[order(-start)]}
    if (strand == "+"){
      introns_shortened = introns_shortened[order(start)]
      exons = exons[order(start)]}
    introns_shortened$intron_num = 1:nrow(introns_shortened)

    # get rescaled length of NA exon and next exon
    val_intron = filter(introns_shortened, intron_num == exon_n)$rescaled_length
    next_start = filter(exons, exon_num == exon_n + 1)$new_e_start

    if(length(val_intron)==0){ #most likely last exon
      val_intron = filter(introns_shortened, intron_num == exon_n-1)$rescaled_length
      prev_end = filter(exons, exon_num == exon_n -1)$new_e_end #actually previous end
      if(!(is.na(prev_end))){
        new_e_start = prev_end + val_intron
        nas$new_e_start[i]= new_e_start
        nas$new_e_end[i] = nas$new_e_start[i] + nas$blocksizes[i]
        ref_attach_to = rbind(ref_attach_to, nas[i])
      }
    }

    if(!(length(next_start)==0)){
        new_e_ending = next_start - val_intron
        nas$new_e_end[i]= new_e_ending
        nas$new_e_start[i] = nas$new_e_end[i] - nas$blocksizes[i]
        ref_attach_to = rbind(ref_attach_to, nas[i])}
    }

  } # (!(dim(nas)[1] == 0))

  ref_attach_to$new_start = NULL
  ref = ref %>% dplyr::select(chr, start, end, strand, blocksizes, transcript_name, new_e_start, new_e_end, exon_num)
  ref_attach_to = rbind(ref_attach_to, ref)
  min_starts = ref_attach_to %>% dplyr::group_by(transcript_name) %>% dplyr::summarize(min_start = min(new_e_start))
  max_exons = ref_attach_to %>% dplyr::group_by(transcript_name) %>% dplyr::summarize(max_end = max(new_e_end))
  trans_lengths = merge(min_starts, max_exons)
  ref_attach_to = merge(ref_attach_to, trans_lengths)
  ref_attach_to$transcript_length = ref_attach_to$max_end - ref_attach_to$min_start   # max exon end and min exon start for each transcript
  return(ref_attach_to)
  } #end of isoviz_rescale_introns

