library(hexSticker)

# get image
#file_path <- system.file("figures", "rbfox2_hesc_filtered.png", package="isoviz")
file_path <- "/Users/kisaev/Documents/github/isoviz/inst/figures/rbfox2_hesc_filtered.png"

sticker(file_path,
             package="isoviz",
             p_color = "orange",
             p_size = 30,
             h_fill="white", h_color="#f39c12",
             filename="inst/figures/logo.png")
