source(test_path("fixtures/make_toy.R"))

test_that("input error class fires for non-data.frame da", {
  expect_error(
    plot_volcano_ring(list(), make_toy_ring_enrichment()),
    class = "enrichVolcano_input_error"
  )
})

test_that("ring_radius must be a single positive number", {
  for (bad in list(0, -3, c(1, 2), "4", NA_real_)) {
    expect_error(
      plot_volcano_ring(make_toy_da(), make_toy_ring_enrichment(), ring_radius = bad),
      class = "enrichVolcano_param_error"
    )
  }
})

test_that("arc_height_range must be c(min, max) with 0 <= min <= max", {
  for (bad in list(1.6, c(1.6, 0.05), c(-1, 1), c(0.1, NA))) {
    expect_error(
      plot_volcano_ring(make_toy_da(), make_toy_ring_enrichment(), arc_height_range = bad),
      class = "enrichVolcano_param_error"
    )
  }
})

test_that("ring_radius below volcano_radius warns about overflow", {
  expect_warning(
    suppressMessages(plot_volcano_ring(make_toy_da(), make_toy_ring_enrichment(),
      ring_radius = 2, volcano_radius = 3.5
    )),
    class = "enrichVolcano_param_warning"
  )
})

test_that("column error class fires for a missing volc column", {
  v <- make_toy_da()
  v$logFC <- NULL
  expect_error(
    plot_volcano_ring(v, make_toy_ring_enrichment()),
    class = "enrichVolcano_column_error"
  )
})

test_that("data error class fires for impossible padj", {
  v <- make_toy_da()
  v$p[1] <- 2
  expect_error(
    plot_volcano_ring(v, make_toy_ring_enrichment()),
    class = "enrichVolcano_data_error"
  )
})

test_that("param error class fires for a bad score_stops length", {
  expect_error(
    plot_theme(score_stops = c(-1, 1)),
    class = "enrichVolcano_param_error"
  )
})

test_that("invalid colour fires its own classed error", {
  expect_error(
    plot_volcano_ring(make_toy_da(), make_toy_ring_enrichment(), disc_colour = "not_a_real_colour"),
    class = "ev_invalid_colour"
  )
})

test_that("every classed error inherits enrichVolcano_error", {
  err <- tryCatch(
    plot_volcano_ring(list(), make_toy_ring_enrichment()),
    error = function(e) e
  )
  expect_s3_class(err, "enrichVolcano_error")
})
