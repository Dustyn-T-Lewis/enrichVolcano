test_that("ring_plot returns a ggplot using geom_arc_bar", {
  enrich <- tibble::tibble(
    contrast = "ctr", database = "h", mode = "fgsea",
    pathway = paste0("HALLMARK_TEST_", 1:5),
    pval = c(0.001, 0.002, 0.005, 0.01, 0.02),
    padj = c(0.005, 0.01, 0.02, 0.05, 0.1),
    NES = c(2, 1.5, -1, -2, 1),
    size = c(50, 40, 30, 60, 20),
    direction = c("up", "up", "down", "down", "up"),
    leading_edge = rep("A;B;C", 5),
    dedup_kept = TRUE
  )
  p <- ring_plot(enrich, contrast = "ctr", max_terms = 5)
  expect_s3_class(p, "ggplot")
  layer_geoms <- vapply(p$layers, function(l) class(l$geom)[1], character(1))
  expect_true(any(grepl("ArcBar", layer_geoms)))
})

test_that("ring_plot warns above 12 terms, aborts above 25", {
  enrich_big <- tibble::tibble(
    contrast = "ctr", database = "h", mode = "fgsea",
    pathway = paste0("PATH_", 1:30),
    pval = runif(30), padj = runif(30),
    NES = rnorm(30), size = rep(10, 30),
    direction = rep("up", 30),
    leading_edge = rep("A;B", 30),
    dedup_kept = TRUE
  )
  expect_warning(ring_plot(enrich_big, contrast = "ctr", max_terms = 15),
                 class = "ev_max_terms_warn")
  expect_error(ring_plot(enrich_big, contrast = "ctr", max_terms = 30),
               class = "ev_max_terms_exceeded")
})

test_that("ring_plot renders placeholder when zero significant pathways", {
  enrich_empty <- tibble::tibble(
    contrast = character(0), database = character(0), mode = character(0),
    pathway = character(0), pval = numeric(0), padj = numeric(0),
    NES = numeric(0), size = integer(0), direction = character(0),
    leading_edge = character(0), dedup_kept = logical(0)
  )
  p <- ring_plot(enrich_empty, contrast = "ctr", max_terms = 10)
  expect_s3_class(p, "ggplot")
})
