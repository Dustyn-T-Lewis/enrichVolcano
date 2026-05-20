test_that("ev_compose returns enrichVolcano-classed patchwork", {
  data <- make_ranked_input()
  data <- pi_score(data, variant = "eq2")
  paths <- make_mini_pathways()
  enrich <- suppressWarnings(ev_enrich(
    data, contrast = "ctr",
    databases = list(test = paths),
    enrich_mode = "fgsea", rank_by = "pi_eq2",
    nperm = 100, min_size = 5, max_size = 100
  ))
  enrich <- ev_collapse(enrich, method = "jaccard", cutoff = 0.5)
  v <- list(ctr = ev_volcano(data, "ctr"))
  r <- list(ctr = ring_plot(enrich, "ctr"))
  p <- ev_compose(v, r)
  expect_s3_class(p, "enrichVolcano")
  expect_s3_class(p, "patchwork")
})

test_that("ev_compose attaches ev_data attribute when supplied", {
  v <- list(ctr = ggplot2::ggplot())
  r <- list(ctr = ggplot2::ggplot())
  p <- ev_compose(v, r, data = list(validated_input = tibble::tibble(x = 1)))
  expect_named(attr(p, "ev_data"), "validated_input")
})

test_that("ev_compose aborts on mismatched names", {
  v <- list(a = ggplot2::ggplot())
  r <- list(b = ggplot2::ggplot())
  expect_error(ev_compose(v, r), class = "ev_compose_name_mismatch")
})
