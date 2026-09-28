# Cross-platform font metrics make exact-pixel diffs flaky on CI. Pin these
# to the local development environment via skip_on_ci(); reviewers can run
# them locally with `testthat::snapshot_review()`.

toy_ring <- function(...) {
  suppressMessages(plot_volcano_ring(make_toy_da(), make_toy_ring_enrichment(),
    databases = NULL, term_threshold = 1, n_terms = Inf, ...
  ))
}

snapshots <- list(
  "toy-default" = list(title = "toy default"),
  "toy-viridis" = list(title = "toy viridis", theme = plot_theme(palette = "viridis")),
  "toy-okabe" = list(title = "toy okabe", theme = plot_theme(palette = "okabe")),
  "toy-magnitude-size" = list(title = "toy size", magnitude = "size"),
  "toy-gene-labels" = list(title = "toy gene labels", label_mode = "top_per_direction", label_n = 3)
)

for (name in names(snapshots)) {
  test_that(paste(name, "snapshot is stable"), {
    skip_on_ci()
    skip_if_not_installed("vdiffr")
    vdiffr::expect_doppelganger(name, do.call(toy_ring, snapshots[[name]]))
  })
}
