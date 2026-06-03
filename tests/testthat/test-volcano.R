test_that("ev_volcano returns a ggplot for a single contrast", {
  data <- make_ranked_input()
  data <- pi_score(data, variant = "eq2")
  p <- ev_volcano(data, contrast = "ctr")
  expect_s3_class(p, "ggplot")
})

test_that("ev_volcano applies significance classes via ev_class column", {
  data <- make_ranked_input()
  data <- pi_score(data, variant = "eq2")
  p <- ev_volcano(data, contrast = "ctr",
                  p_threshold = 0.05, logfc_threshold = log2(1.5))
  layer_data <- p$data
  expect_true("ev_class" %in% colnames(layer_data))
  expect_true(all(layer_data$ev_class %in% c("up", "down", "ns")))
})

test_that("ev_volcano labels top N proteins via ggrepel", {
  data <- make_ranked_input()
  data <- pi_score(data, variant = "eq2")
  p <- ev_volcano(data, contrast = "ctr",
                  label_mode = "top_total", label_n = 3)
  has_repel <- any(vapply(p$layers, function(l) {
    inherits(l$geom, "GeomTextRepel") || grepl("repel", class(l$geom)[1], ignore.case = TRUE)
  }, logical(1)))
  expect_true(has_repel)
})

test_that("ev_volcano aborts on unknown contrast", {
  data <- make_ranked_input()
  data <- pi_score(data, variant = "eq2")
  expect_error(ev_volcano(data, contrast = "nonexistent"),
               class = "ev_unknown_contrast")
})
