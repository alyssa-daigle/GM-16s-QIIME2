invisible(
  c(
    "MCMCglmm", #
    "phyloseq", #
    "ggplot2", #
    "tidyr",
    "purrr", #
    "dplyr",
    "cowplot",
    "vegan",
    "tibble",
    "ggtext",
    "dotenv",
    "compositions",
    "corncob"
  ) |>
    lapply(function(x) {
      if (suppressMessages(!require(x, character.only = TRUE))) {
        install.packages(x)
        library(x, character.only = TRUE)
      }
    })
)
