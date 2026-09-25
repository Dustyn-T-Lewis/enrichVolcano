# Three up, one down, all significant; B has half of A's score.
uneven_enrichment <- function() {
  tbl <- data.frame(
    contrast = "toy", term = paste0("HALLMARK_TOY_", c("A", "B", "C", "D", "E")),
    score = c(2.4, 1.2, 0.6, -1.8, 0.3), padj = c(0.01, 0.01, 0.01, 0.01, 0.5), size = 20
  )
  suppressMessages(as_enrichment(tbl, enrichment_test = "custom", score_type = "NES"))
}

bias <- function(x = uneven_enrichment(), ...) {
  suppressMessages(plot_bias_ring(make_toy_da(), x, ...))
}

arc_layer <- function(p) {
  geoms <- vapply(p$layers, function(l) class(l$geom)[1], character(1))
  p$layers[[utils::tail(which(geoms == "GeomArcBar"), 1)]]
}

test_that("plot_bias_ring draws every significant term and no term names", {
  p <- bias()
  expect_s3_class(p, "ggplot")
  expect_setequal(arc_layer(p)$data$term, paste0("HALLMARK_TOY_", c("A", "B", "C", "D")))
  drawn <- unlist(lapply(ggplot2::ggplot_build(p)$data, function(d) if ("label" %in% names(d)) as.character(d$label)))
  expect_false(any(grepl("Toy", drawn)))
})

test_that("fill is the score over the strongest score, on a -1 to 1 scale", {
  p <- bias()
  fill <- arc_layer(p)$data$score
  expect_equal(max(abs(fill)), 1)
  expect_equal(sort(fill), sort(c(2.4, 1.2, 0.6, -1.8) / 2.4))
  expect_identical(p$scales$get_scales("fill")$limits, c(-1, 1))
  expect_identical(p$scales$get_scales("fill")$name, "NES / max")
})

test_that("arc height is proportional to the normalised score", {
  ring <- arc_layer(bias(arc_height_range = c(0.4, 1.6)))$data
  base <- 4.8 + 0.55
  top <- stats::setNames(ring$arc_r1_var - base, ring$term)
  expect_equal(top[["HALLMARK_TOY_A"]], 1.6)
  expect_equal(top[["HALLMARK_TOY_B"]], 0.8)
})

test_that("the subtitle counts up and down terms", {
  expect_identical(bias()$labels$subtitle, "3 up, 1 down of 4 significant terms")
  expect_identical(bias(subtitle = "mine")$labels$subtitle, "mine")
})

test_that("a contrast with no significant terms still draws and says why", {
  expect_message(
    p <- plot_bias_ring(make_toy_da(), uneven_enrichment(), term_threshold = 1e-9),
    class = "enrichVolcano_empty_ring"
  )
  expect_s3_class(p, "ggplot")
})

test_that("layout arguments pass through and anything else is refused", {
  expect_s3_class(bias(ring_radius = 5.2, point_size = 0.5), "ggplot")
  expect_error(bias(colour = "red"), "colour", class = "enrichVolcano_param_error")
  expect_error(
    plot_bias_ring(make_toy_da(), uneven_enrichment(), "toy", NULL, TRUE, 0.05, NULL, NULL, plot_theme(), 5),
    "unnamed",
    class = "enrichVolcano_param_error"
  )
})

test_that("toy bias ring snapshot is stable", {
  skip_on_ci()
  skip_if_not_installed("vdiffr")
  vdiffr::expect_doppelganger("toy-bias-ring", bias(title = "toy bias"))
})
