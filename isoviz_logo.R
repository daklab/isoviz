library(hexSticker)

# get image
file_path <- "/Users/kisaev/Documents/github/isoviz/inst/figures/Readme_RBFOX2.png"

library(ggplot2)
library(isoviz)

## basic example code
file_path <- system.file("data", "gencode.v41.basic.annotation.psl", package="isoviz")
gene_trans <- system.file("data", "gencode_v41_gene_transcript_convert.txt", package="isoviz")

all_coordinates <- isoviz_coords(file_path, gene_trans, input_type="psl") #use default genome .psl file

exon_coords <- all_coordinates[[1]]
intron_coords <- all_coordinates[[2]]

rbfox2_exons <- filter(exon_coords, gene_name == "RBFOX2")
rbfox2_introns <- filter(intron_coords, gene_name == "RBFOX2")

rescaled_coords = isoviz_rescale_introns(rbfox2_introns, rbfox2_exons, width_rescale=10)

# load junctions for cell type of interest or input your own
junctions <- system.file("data", "hESC-MKNK2-G1_v41_basic.junc", package="isoviz")

# run minicutter to get clusters
intron_clusts <- isoviz_minicutter(juncs_file = junctions)
print(head(intron_clusts))

intron_annotations <- system.file("data", "gencode_intron_all_data.rda", package="isoviz")
load(intron_annotations)

mapped_junctions = isoviz_map_junctions(cell_type = "hESC", rbfox2_introns, intron_clusts, gencode_intron_all_data)

p = isoviz_plot_juncs_to_iso(mapped_junctions, rbfox2_exons, rbfox2_introns,
                         cell_type = "hESC",
                         junc_usage = 5, #min junc usage to be included
                         include_all_juncs = FALSE,
                         include_specific_junctions = c("junc178149",
                                                        "junc178148",
                                                        "junc178145"))

sticker(p,
        package="isoviz",
        p_color = "orange",
        p_size=30, s_x=1.2, s_y=0.5, s_width=2, s_height=1.5,
        h_fill="white", h_color="#f39c12",
        filename="inst/figures/isoviz.png")
