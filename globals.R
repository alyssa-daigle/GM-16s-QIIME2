packages <- c(
  "MCMCglmm",
  "phyloseq",
  "ggplot2",
  "tidyr",
  "purrr",
  "dplyr",
  "cowplot",
  "vegan",
  "tibble",
  "ggtext",
  "dotenv",
  "compositions",
  "corncob"
)

invisible(
  lapply(packages, function(pkg) {
    if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
      install.packages(pkg, quiet = TRUE)
      suppressPackageStartupMessages(
        library(pkg, character.only = TRUE, quietly = TRUE)
      )
    } else {
      suppressPackageStartupMessages(
        library(pkg, character.only = TRUE, quietly = TRUE)
      )
    }
  })
)
