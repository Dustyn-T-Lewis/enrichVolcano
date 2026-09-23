source(test_path("fixtures/make_toy.R"))

vcols_default <- function() {
  list(gene = "gene", logfc = "logFC", pval = "P.Value", padj = "padj", has_padj = TRUE)
}

test_that("validate_volc_df accepts a tidy frame", {
  expect_silent(validate_volc_df(make_toy_volc(), vcols_default()))
})

test_that("validate_volc_df rejects non-data.frame input", {
  expect_error(
    validate_volc_df(list(), vcols_default()),
    class = "enrichVolcano_input_error"
  )
})

test_that("validate_volc_df rejects an empty frame", {
  v <- make_toy_volc()[0, ]
  expect_error(
    validate_volc_df(v, vcols_default()),
    class = "enrichVolcano_input_error"
  )
})

test_that("validate_volc_df rejects non-numeric logFC", {
  v <- make_toy_volc()
  v$logFC <- as.character(v$logFC)
  expect_error(
    validate_volc_df(v, vcols_default()),
    class = "enrichVolcano_column_error"
  )
})

test_that("validate_volc_df rejects a non-numeric P-value column", {
  v <- make_toy_volc()
  v$P.Value <- as.character(v$P.Value)
  expect_error(
    validate_volc_df(v, vcols_default()),
    class = "enrichVolcano_column_error"
  )
})

test_that("validate_volc_df rejects P-values outside [0, 1]", {
  v <- make_toy_volc()
  v$P.Value[1] <- 1.5
  expect_error(
    validate_volc_df(v, vcols_default()),
    class = "enrichVolcano_data_error"
  )
})

test_that("validate_ring_geometry rejects a bad label_headroom", {
  for (bad in list(-1, c(0.5, 1), "0.5", NA_real_)) {
    expect_error(
      validate_ring_geometry(4.8, 4.0, c(0.4, 1.6), label_headroom = bad),
      class = "enrichVolcano_param_error"
    )
  }
})

test_that("validate_grid_spacing rejects a negative or non-scalar knob", {
  expect_error(
    validate_grid_spacing(-1, 2, 26),
    class = "enrichVolcano_param_error"
  )
  expect_error(
    validate_grid_spacing(1.5, c(2, 3), 26),
    class = "enrichVolcano_param_error"
  )
})
