pond_name_mapping <- c(
  "MP-1" = "Mill Pond",
  "ODR-2" = "Dairy Farm\nPond 1",
  "ODR-3" = "Dairy Farm\nPond 2",
  "TF-1" = "Thompson Farm\nPond 1",
  "TF-2" = "Thompson Farm\nPond 2",
  "UM-1" = "Upper Mill Pond"
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
