test_that("the public exports are visible", {
  exports <- c(
    "as_enrichment", "dedup", "enrichment", "ev_clean_label", "nes_scatter",
    "volcano_ring", "volcano_ring_theme"
  )
  expect_true(all(exports %in% getNamespaceExports("enrichVolcano")))
})

test_that("R code is plain ASCII, so figures print on any device", {
  r_dir <- test_path("..", "..", "R")
  skip_if_not(dir.exists(r_dir), "R sources not available (installed package)")
  files <- list.files(r_dir, pattern = "\\.R$", full.names = TRUE)
  non_ascii <- function(f) any(grepl("[^\\x01-\\x7F]", readLines(f), perl = TRUE, useBytes = TRUE))
  offending <- Filter(non_ascii, files)
  expect_identical(basename(offending), character(0))
})
