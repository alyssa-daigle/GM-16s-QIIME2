pond_name_mapping <- c(
  "MP-1" = "Mill Pond",
  "ODR-2" = "Dairy Farm\nPond 1",
  "ODR-3" = "Dairy Farm\nPond 2",
  "TF-1" = "Thompson Farm\nPond 1",
  "TF-2" = "Thompson Farm\nPond 2",
  "UM-1" = "Upper \nMill Pond"
)

pond_levels <- c(
  "Mill Pond",
  "Upper \nMill Pond",
  "Dairy Farm\nPond 1",
  "Dairy Farm\nPond 2",
  "Thompson Farm\nPond 1",
  "Thompson Farm\nPond 2"
)

expt2_custom_colors <- c(
  "#44AA99",
  "#CC6677",
  "#88CCEE",
  "#117733",
  "#DDCC77",
  "#332288"
)

pca_theme <- function(base_size = 7, base_family = "") {
  theme(
    legend.position = "right",
    legend.title = element_markdown(size = 10),
    legend.text = element_text(size = 8),
    strip.text = element_text(size = 5.7)
  )
}

gravityFalls_colors <- c(
  "#474747FF", # gray
  "#8B4513FF", # brown
  "#D2B48CFF", # tan (brownish)
  "#000000FF", # black
  "#417BA1FF", # blue
  "hotpink", # pink
  "#FFFF2EFF", # yellow
  "#345634FF", # dark green
  "#8B0000FF", # dark red
  "#E2725B", # orange
  "#93C0D5FF", # light blue
  "#9248A7FF", # purple
  "#1C8859FF", # green
  "pink2", # pink
  "#8FBC8FFF" # light green
)
