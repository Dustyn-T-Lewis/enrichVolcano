source(test_path("fixtures/make_toy.R"))

# Cross-platform font metrics make exact-pixel diffs flaky on CI. Pin these
# to the local development environment via skip_on_ci(); reviewers can run
# them locally with `testthat::snapshot_review()`.

test_that("default palette toy snapshot is stable", {
  skip_on_ci()
  skip_if_not_installed("vdiffr")
  p <- suppressMessages(volcano_ring(make_toy_da(), make_toy_ring_enrichment(),
    databases = NULL, term_threshold = 1, n_terms = Inf,
    title = "toy default"
  ))
  vdiffr::expect_doppelganger("toy-default", p)
})

test_that("viridis palette toy snapshot is stable", {
  skip_on_ci()
  skip_if_not_installed("vdiffr")
  p <- suppressMessages(volcano_ring(make_toy_da(), make_toy_ring_enrichment(),
    databases = NULL, term_threshold = 1, n_terms = Inf,
    title = "toy viridis",
    theme = volcano_ring_theme(palette = "viridis")
  ))
  vdiffr::expect_doppelganger("toy-viridis", p)
})

test_that("magnitude = 'size' toy snapshot is stable", {
  skip_on_ci()
  skip_if_not_installed("vdiffr")
  p <- suppressMessages(volcano_ring(make_toy_da(), make_toy_ring_enrichment(),
    databases = NULL, term_threshold = 1, n_terms = Inf,
    title = "toy size", magnitude = "size"
  ))
  vdiffr::expect_doppelganger("toy-magnitude-size", p)
})
