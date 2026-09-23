test_that("the public exports are visible", {
  exports <- c(
    "enrichment", "ev_clean_label", "volcano_ring", "volcano_ring_grid",
    "volcano_ring_theme"
  )
  expect_true(all(exports %in% getNamespaceExports("enrichVolcano")))
  expect_true(!is.null(getS3method("print", "volcano_ring_grid")))
})
