source(test_path("fixtures/make_toy.R"))

# Cross-platform font metrics make exact-pixel diffs flaky on CI. Pin these
# to the local development environment via skip_on_ci(); reviewers can run
# them locally with `testthat::snapshot_review()`.

test_that("default palette toy snapshot is stable", {
  skip_on_ci()
  skip_if_not_installed("vdiffr")
  p <- suppressMessages(volcano_ring(make_toy_volc(), make_toy_enrich(),
    title = "toy default"
  ))
  vdiffr::expect_doppelganger("toy-default", p)
})

test_that("viridis palette toy snapshot is stable", {
  skip_on_ci()
  skip_if_not_installed("vdiffr")
  p <- suppressMessages(volcano_ring(make_toy_volc(), make_toy_enrich(),
    title = "toy viridis",
    theme = volcano_ring_theme(palette = "viridis")
  ))
  vdiffr::expect_doppelganger("toy-viridis", p)
})

test_that("magnitude = 'size' toy snapshot is stable", {
  skip_on_ci()
  skip_if_not_installed("vdiffr")
  p <- suppressMessages(volcano_ring(make_toy_volc(), make_toy_enrich(),
    title = "toy size", magnitude = "size"
  ))
  vdiffr::expect_doppelganger("toy-magnitude-size", p)
})

test_that("2-panel grid toy snapshot is stable", {
  skip_on_ci()
  skip_if_not_installed("vdiffr")
  g <- suppressMessages(volcano_ring_grid(
    list(A = make_toy_volc(seed = 1L), B = make_toy_volc(seed = 2L)),
    list(A = make_toy_enrich(seed = 1L), B = make_toy_enrich(seed = 2L))
  ))
  vdiffr::expect_doppelganger("toy-grid-2x", g$plot)
})
