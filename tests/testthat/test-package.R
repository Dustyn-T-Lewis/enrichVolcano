test_that("the four public exports are visible", {
  exports <- c("volcano_ring", "volcano_ring_grid", "volcano_ring_theme")
  for (fn in exports) {
    expect_true(exists(fn, envir = asNamespace("enrichVolcano")))
  }
  expect_true(!is.null(getS3method("print", "volcano_ring_grid")))
})
