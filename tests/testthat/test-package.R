test_that("the public exports are visible", {
  exports <- c(
    "as_da", "as_enrichment", "clean_label", "dedup_terms", "load_gene_sets", "plot_bias_ring", "plot_scatter",
    "plot_theme", "plot_volcano_ring", "read_example", "read_study", "run_enrichment", "write_plot", "write_table"
  )
  expect_setequal(getNamespaceExports("enrichVolcano"), exports)
})

test_that("R code is plain ASCII, so figures print on any device", {
  r_dir <- test_path("..", "..", "R")
  skip_if_not(dir.exists(r_dir), "R sources not available (installed package)")
  files <- list.files(r_dir, pattern = "\\.R$", full.names = TRUE)
  non_ascii <- function(f) any(grepl("[^\\x01-\\x7F]", readLines(f), perl = TRUE, useBytes = TRUE))
  offending <- Filter(non_ascii, files)
  expect_identical(basename(offending), character(0))
})
