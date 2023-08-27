#' A function for summarizing junction and transcript level information from GENCODE annotated transcripts or long read sequencing transcripts.
#' This is especially useful to identify unique junctions that are specific to a single transcript.
#'
#' Inputs gencode intron rda file or long read intron rda file (output from isoviz_coords), outputs summary table
#' and summary plot
#' 
#' @param iso_intron_data This is output from isoviz_coords.
#' @param type Either 'gencode' or 'longread' depending on the initial source of transcript information. Default is gencode.
#' @param plot_path This is optional, but if you would like a summary plot, give the path and name of the output. Include proper extension such as '.pdf' or .'png'. Default is no plot.
#' @return table summarizing junction and transcript level information.
#' 
#' @name isoviz_junction_to_transcript_summary
#' @export

isoviz_junction_to_transcript_summary = function(iso_intron_data, type = "gencode", plot_path = ""){
  
  # check input file
  if(ncol(iso_intron_data) != 10){
    stop("Input file has a different number of columns than expected. Columns should be: chr, intron_starts, intron_ends, gene_id, trans_id, strand, gene_name, transcript_name, gene_type, transcript_type")
  }
  
  # filter for protein coding and lncRNA genes
  iso_intron_pc = iso_intron_data %>% dplyr::filter(gene_type == "protein_coding" | gene_type == "lncRNA")
  
  # filter for protein-coding and lncRNA transcripts- we can only do this for gencode annotations since we wont have this info for novel transcripts in long read data
  if(type == "gencode"){
    iso_intron_pc %<>% dplyr::filter(transcript_type == "protein_coding" | transcript_type == "lncRNA")
  }
  
  # remove transcript ids and collapse on junction coords- this gives the total number of junctions in your annotations
  # assign universal ids for junction
  collapsed_junctions = iso_intron_pc %>% dplyr::select(-trans_id, -transcript_name, -gene_type, -transcript_type) %>% 
    distinct() %>% dplyr::select(chr, start = intron_starts, end = intron_ends, gene_id, gene_name, strand) %>% 
    dplyr::group_by(chr, start, end, gene_id, gene_name, strand) %>%
    dplyr::mutate(junc_id = paste0("junc", cur_group_id())) %>% ungroup() # 281,880 junctions- how many show any uniqueness?

  # to address how many annotated isoforms have a truly unique junction
  # first filter the genes with more than one transcript
  # 85,586 transcripts for 33,417 genes represented, 17,035 genes have more than 1 isoform
  n_df = iso_intron_pc %>% dplyr::select(gene_id, trans_id) %>% distinct() %>% dplyr::group_by(gene_id) %>% dplyr::mutate(n_transcripts_per_gene = n()) %>% ungroup()
  n_genes = nrow(n_df %>% dplyr::group_by(gene_id) %>% tally())
  print(paste0("There are ", nrow(collapsed_junctions), " junctions across ", nrow(n_df), " transcripts for ", n_genes, " protein-coding and lncRNA genes."))

  # define junction type
  # the logic here has to be in the right order for the ifelse statement
  nrow(iso_intron_pc) # 696,609
  iso_pc_df = iso_intron_pc %>% full_join(n_df, by = c("gene_id", "trans_id")) %>%
    dplyr::group_by(chr, intron_starts, intron_ends, gene_id, strand, gene_name, n_transcripts_per_gene) %>% 
    dplyr::mutate(n_trans_per_junc = n(), commonality = round((n_trans_per_junc/n_transcripts_per_gene)*100, 0)) %>%
    dplyr::mutate(junction_category = ifelse(n_transcripts_per_gene == 1 & n_trans_per_junc == 1, "single_isoform", 
                                      ifelse(n_trans_per_junc == 1, "fully_unique", ifelse(commonality == 100, "common", "partially_unique")))) %>% ungroup()
  
  # total genes, total transcripts
  g = iso_pc_df %>% dplyr::select(gene_id, gene_type) %>% distinct() %>% group_by(gene_type) %>% tally()
  t = iso_pc_df %>% dplyr::select(transcript_name, gene_type) %>% distinct() %>% group_by(gene_type) %>% tally() 
  j = iso_pc_df %>% dplyr::select(chr, intron_starts, intron_ends, strand, gene_id, gene_name, gene_type, junction_category) %>% distinct() %>% group_by(gene_type) %>% tally()
  
  print(paste0("LncRNA summary: There are ", j$n[j$gene_type == "lncRNA"], " junctions across ", t$n[t$gene_type == "lncRNA"], " transcripts for ", g$n[g$gene_type == "lncRNA"], " lncRNA genes."))
  print(paste0("Protein-coding gene summary: There are ", j$n[j$gene_type == "protein_coding"], " junctions across ", t$n[t$gene_type == "protein_coding"], " transcripts for ", g$n[g$gene_type == "protein_coding"], " protein-coding genes."))
  
  junction_df = iso_pc_df %>% dplyr::select(chr, junc_start = intron_starts, junc_end = intron_ends, gene_id, 
                                  gene_name, strand, transcript_name, gene_type, junction_category) %>% 
    dplyr::group_by(chr, junc_start, junc_end, gene_id, gene_name, strand, gene_type, junction_category) %>%
    dplyr::mutate(transcript_isoforms = paste0(transcript_name, collapse = ",")) %>% ungroup() %>% dplyr::select(-transcript_name) %>% distinct()
  
  # get lists of transcripts for each junction category
  unique_list = iso_pc_df %>% dplyr::filter(junction_category == "fully_unique") %>% dplyr::select(gene_id, gene_name, trans_id, transcript_name, gene_type) %>% 
    distinct() %>% dplyr::select(trans_id) %>% mutate(Isoform_targetable = "fully_unique")
  
  partial_list = iso_pc_df %>% dplyr::filter(junction_category == "partially_unique") %>% 
    dplyr::select(gene_id, gene_name, trans_id, transcript_name, gene_type) %>% distinct() %>% dplyr::select(trans_id) %>% dplyr::mutate(Isoform_targetable = "partially_unique")
  
  common_list = iso_pc_df %>% dplyr::filter(junction_category == "common") %>% 
    dplyr::select(gene_id, gene_name, trans_id, transcript_name, gene_type) %>% distinct() %>% dplyr::select(trans_id) %>% dplyr::mutate(Isoform_targetable = "common")
  
  single_isoform_list = iso_pc_df %>% dplyr::filter(junction_category == "single_isoform") %>% 
    dplyr::select(gene_id, gene_name, trans_id, transcript_name, gene_type) %>% distinct() %>% dplyr::select(trans_id) %>% dplyr::mutate(Isoform_targetable = "single_isoform")
  
  # order matters here for ifelse. my transcripts will have a combination of junctions. But if a transcript has even 1 unique junction, then it is uniquely targetable, so we start with that
  transcript_df = iso_pc_df %>% dplyr::mutate(Isoform_targetable = ifelse(trans_id %in% single_isoform_list$trans_id, "single_isoform", 
                                                    ifelse(trans_id %in% unique_list$trans_id, "fully_unique", 
                                                           ifelse(trans_id %in% partial_list$trans_id, "partially_unique",
                                                                  ifelse(trans_id %in% common_list$trans_id, "common", "No_category")))))
  
  summary = transcript_df %>% left_join(junction_df, by = c("chr", "intron_starts" = "junc_start", "intron_ends" = "junc_end", "gene_id", "strand", "gene_name", "gene_type", "junction_category")) %>%
    dplyr::select(-n_transcripts_per_gene, -n_trans_per_junc, -commonality) # 696,609
  
  if(plot_path != ""){
    
    # junction-level
    j_summary = transcript_df %>% dplyr::select(chr, intron_starts, intron_ends, strand, gene_id, gene_name, gene_type, category = junction_category) %>% distinct() %>% 
      group_by(gene_type, category) %>% tally() %>% ungroup() %>%
      dplyr::mutate(level = "Junction-level")
    
    t_summary = transcript_df %>% dplyr::select(trans_id, gene_type, category = Isoform_targetable) %>% distinct() %>% group_by(gene_type, category) %>% tally() %>% 
      ungroup() %>% dplyr::mutate(level = "Transcript-level")
    
    # Make summary plot
    
    df = rbind(j_summary, t_summary)
    df$gene_type = factor(df$gene_type, levels = c("lncRNA", "protein_coding"))
    df$category = factor(df$category, levels = c("fully_unique", "partially_unique", "single_isoform", "common"))
    
    max_n = max(df$n)
    
    ggplot(df) + geom_bar(aes(x = gene_type, y = n, fill = category, group = category), position="dodge",stat="identity") + coord_flip() + facet_wrap(~level) + theme_bw() + 
      scale_fill_manual(values = c("#F39B7FFF", "#3C5488FF", "#00A087FF", "#4DBBD5FF"), name = "Category", labels = c("Fully Unique", "Partially Unique", "Single Isoform", "Common")) + 
      labs(x = "", y = "", title = "Summary of Transcript Annotations") + 
      geom_text(aes(x = gene_type, y = ifelse(n < max_n, n, max_n*0.75), label = n, group = category), hjust = 0, size = 4, position = position_dodge(width = 0.9), inherit.aes = TRUE) +
      theme(axis.text = element_text(size = 12), text = element_text(size=14), plot.title = element_text(hjust = 0.5), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
            legend.position = "top") +
      guides(fill = guide_legend(title = "", nrow = 2, byrow = TRUE))
    
    ggsave(plot_path, width = 6, height = 4)
    
  }
  
  return(summary)

}
