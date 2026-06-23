source(test_path("fixtures/make_toy.R"))

test_that("input error class fires for non-data.frame volc_df", {
  expect_error(
    volcano_ring(list(), make_toy_enrich()),
    class = "enrichVolcano_input_error"
  )
})

test_that("input error class fires for non-data.frame enrich_df", {
  expect_error(
    volcano_ring(make_toy_volc(), list()),
    class = "enrichVolcano_input_error"
  )
})

test_that("input error class fires when grid input is malformed", {
  expect_error(
    volcano_ring_grid(42, list()),
    class = "enrichVolcano_input_error"
  )
})

test_that("column error class fires for a missing volc column", {
  v <- make_toy_volc()
  v$logFC <- NULL
  expect_error(
    volcano_ring(v, make_toy_enrich()),
    class = "enrichVolcano_column_error"
  )
})

test_that("column error class fires for a missing enrich column", {
  e <- make_toy_enrich()
  e$NES <- NULL
  expect_error(
    volcano_ring(make_toy_volc(), e),
    class = "enrichVolcano_column_error"
  )
})

test_that("data error class fires for impossible padj", {
  v <- make_toy_volc()
  v$P.Value[1] <- 2
  expect_error(
    volcano_ring(v, make_toy_enrich()),
    class = "enrichVolcano_data_error"
  )
})

test_that("data error class fires when every NES is NA", {
  e <- make_toy_enrich()
  e$NES <- NA_real_
  expect_error(
    volcano_ring(make_toy_volc(), e),
    class = "enrichVolcano_data_error"
  )
})

test_that("param error class fires for a bad nes_stops length", {
  expect_error(
    volcano_ring_theme(nes_stops = c(-1, 1)),
    class = "enrichVolcano_param_error"
  )
})

test_that("invalid colour fires its own classed error", {
  expect_error(
    volcano_ring(make_toy_volc(), make_toy_enrich(), disc_color = "not_a_real_colour"),
    class = "ev_invalid_colour"
  )
})

test_that("every classed error inherits enrichVolcano_error", {
  err <- tryCatch(
    volcano_ring(list(), make_toy_enrich()),
    error = function(e) e
  )
  expect_s3_class(err, "enrichVolcano_error")
})
