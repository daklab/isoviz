## code to prepare `genomedata` dataset goes here

# data-raw/genomedata/R
# code to prepare `genomedata` dataset goes here

require(readr)
require(tidyverse)
require(stringr)

mydataset <- fread("data-raw/gencode.v41.basic.annotation.psl")

# Data cleaning code here...

# loop creates one line for each annotated exon,
#combining the blocksizes column witht the start position column in the psl file

mydataset$gene_id = sapply(mydataset$V10, function(x){strsplit(x, "_")[[1]][2]})
mydataset$trans_id = sapply(mydataset$V10, function(x){strsplit(x, "_")[[1]][1]})

# Get coordinates for each exon "block"
mydataset$blocksizes = sapply(mydataset$V19, function(x){paste(strsplit(x, ",")[[1]], collapse=",")})
mydataset$blockstarts = sapply(mydataset$V21, function(x){paste(strsplit(x, ",")[[1]], collapse=",")})
mydataset$strand = mydataset$V9
mydataset$chr = mydataset$V14
mydataset$start = mydataset$V16
mydataset$end = mydataset$V17
mydataset$transcript_length = mydataset$end - mydataset$start

# Clean up and seperate blocks into rows
mydataset = mydataset %>% dplyr::select(chr, start, end, trans_id, gene_id,strand, blocksizes, blockstarts,
                            transcript_length) %>% separate_rows(blocksizes, blockstarts)
mydataset$blockstarts=as.numeric(mydataset$blockstarts)
mydataset$blocksizes=as.numeric(mydataset$blocksizes)
mydataset$blockends = mydataset$blockstarts + mydataset$blocksizes

iso_exon_data = as.data.table(mydataset)

# Save to package
usethis::use_data(iso_exon_data, overwrite = TRUE)
