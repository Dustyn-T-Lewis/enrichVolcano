source(test_path("fixtures/make_toy.R"))

scatter <- function(x = make_toy_scatter_enrichment(), ...) {
  suppressMessages(plot_scatter(x, "A", "B", databases = NULL, ...))
}

layer_labels <- function(p) {
  unlist(lapply(ggplot2::ggplot_build(p)$data, function(d) {
    if ("label" %in% names(d)) as.character(d$label) else character(0)
  }))
}

test_that("significance classes are named after the contrasts", {
  s <- scatter_significance(c(0.01, 0.01, 0.2, 0.2), c(0.01, 0.2, 0.01, 0.2), 0.05, "Young", "Old")
  expect_identical(as.character(s), c("Both", "Young only", "Old only", "NS"))
  expect_identical(levels(s), c("Young only", "Old only", "Both", "NS"))
})

test_that("a missing padj counts as not significant", {
  s <- scatter_significance(c(NA, 0.01), c(0.01, NA), 0.05, "A", "B")
  expect_identical(as.character(s), c("B only", "A only"))
})

test_that("quadrant counts and concordance match the hand count", {
  x <- make_toy_scatter_enrichment()@results
  a <- x[x$contrast == "A", ]
  b <- x[x$contrast == "B", ]
  sig <- scatter_significance(a$padj, b$padj, 0.05, "A", "B") != "NS"
  counts <- quadrant_counts(a$score, b$score, sig)
  expect_identical(counts, c(top_right = 4L, top_left = 1L, bottom_left = 2L, bottom_right = 2L))
  expect_equal(concordance_fraction(counts), 6 / 9)
})

test_that("points on an axis belong to no quadrant, and no points give NA", {
  counts <- quadrant_counts(c(0, 1), c(1, 0), c(TRUE, TRUE))
  expect_identical(sum(counts), 0L)
  expect_true(is.na(concordance_fraction(counts)))
})

test_that("the correlation carries a CI when the correlation package is there", {
  skip_if_not_installed("correlation")
  x <- c(1, 2, 3, 4, 5, 6, 7, 8)
  y <- c(2, 1, 4, 3, 6, 5, 8, 7)
  out <- score_correlation(x, y)
  ref <- correlation::cor_test(data.frame(x, y), "x", "y", method = "spearman")
  expect_equal(out$rho, ref$rho)
  expect_equal(out$ci, c(ref$CI_low, ref$CI_high))
  expect_equal(out$p, ref$p)
})

test_that("without the correlation package it reports rho and p only", {
  x <- c(1, 2, 3, 4, 5, 6, 7, 8)
  y <- c(2, 1, 4, 3, 6, 5, 8, 7)
  expect_message(out <- score_correlation(x, y, with_ci = FALSE), class = "enrichVolcano_no_ci")
  ref <- suppressWarnings(stats::cor.test(x, y, method = "spearman", exact = FALSE))
  expect_equal(out$rho, unname(ref$estimate))
  expect_equal(out$p, ref$p.value)
  expect_null(out$ci)
})

test_that("fewer than three points give no correlation and a shorter subtitle", {
  expect_null(score_correlation(c(1, 2), c(2, 1)))
  expect_identical(scatter_subtitle(NULL, NA, "concordant", n = 2, n_sig = 0), "2 terms, 0 significant")
})

test_that("the subtitle states rho, CI, p, concordance and counts", {
  s <- scatter_subtitle(list(rho = 0.62, ci = c(0.41, 0.77), p = 2e-5), 6 / 9, "concordant", n = 12, n_sig = 9)
  expect_match(s, "^rho = 0.62 \\[0.41, 0.77\\]")
  expect_match(s, "p < 0.001")
  expect_match(s, "67% concordant")
  expect_match(s, "12 terms, 9 significant")
  s <- scatter_subtitle(list(rho = 0.1, ci = NULL, p = 0.42), NA, "concordant", n = 5, n_sig = 0)
  expect_match(s, "^rho = 0.10, p = 0.42")
  expect_false(grepl("concordant", s))
})

test_that("plot_scatter returns a ggplot labelled with score type and contrasts", {
  p <- scatter()
  expect_s3_class(p, "ggplot")
  expect_identical(p$labels$x, "NES (A)")
  expect_identical(p$labels$y, "NES (B)")
})

test_that("corner labels name each quadrant and carry its count", {
  labs <- layer_labels(scatter())
  expect_true(all(c("Concordant up\nn = 4", "Discordant\nn = 1", "Concordant down\nn = 2") %in% labs))
})

test_that("reversal framing flips the reference line, names and headline share", {
  p <- scatter(comparison = "reversal")
  labs <- layer_labels(p)
  expect_true(all(c("Exacerbated\nn = 4", "Reversed\nn = 1") %in% labs))
  abline <- Filter(function(d) "slope" %in% names(d), ggplot2::ggplot_build(p)$data)[[1]]
  expect_identical(abline$slope, -1)
  expect_match(paste(deparse(p$labels$subtitle), collapse = ""), "33% reversed")
  expect_error(scatter(comparison = "sideways"), class = "rlang_error")
})

test_that("only significant terms of at least label_min_size are labelled, up to label_n", {
  labs <- layer_labels(scatter())
  expect_true(any(grepl("T01", labs)))
  expect_false(any(grepl("T12", labs)))
  expect_false(any(grepl("T06", labs)))
  labs_few <- layer_labels(scatter(label_n = 2))
  expect_identical(sum(grepl("^T[0-9]+$", labs_few)), 2L)
})

test_that("terms missing from one contrast are dropped with a note", {
  x <- make_toy_scatter_enrichment()
  r <- x@results
  x@results <- r[!(r$contrast == "B" & r$term == "T01"), ]
  expect_message(plot_scatter(x, "A", "B", databases = NULL), "1 term")
})

test_that("collapse hides terms redundant somewhere and representative nowhere", {
  x <- make_toy_scatter_enrichment()
  r <- x@results
  r$dedup_status <- NA_character_
  r$dedup_status[r$term == "T02" & r$contrast == "A"] <- "redundant"
  r$dedup_status[r$term == "T03" & r$contrast == "A"] <- "redundant"
  r$dedup_status[r$term == "T03" & r$contrast == "B"] <- "kept"
  x@results <- r
  pts <- function(p) {
    d <- Filter(function(d) "shape" %in% names(d), ggplot2::ggplot_build(p)$data)
    sum(vapply(d, nrow, integer(1)))
  }
  expect_identical(pts(scatter(x)), 11L)
  expect_identical(pts(scatter(x, collapse = FALSE)), 12L)
})

test_that("contrasts must be two different ones that exist", {
  x <- make_toy_scatter_enrichment()
  expect_error(plot_scatter(x, "A", "A"), class = "enrichVolcano_input_error")
  expect_error(plot_scatter(x, "A", "Z"), "Z", class = "enrichVolcano_input_error")
  expect_error(plot_scatter(make_toy_enrich(), "A", "B"), class = "enrichVolcano_input_error")
})

test_that("colour_by and shape_by take any column, or NULL", {
  expect_s3_class(scatter(colour_by = "database", shape_by = NULL), "ggplot")
  expect_s3_class(scatter(colour_by = NULL), "ggplot")
  expect_s3_class(scatter(colour_by = "size"), "ggplot")
  expect_error(scatter(colour_by = "nope"), class = "enrichVolcano_column_error")
  expect_error(scatter(shape_by = "nope"), class = "enrichVolcano_column_error")
})

test_that("significance falls back to nominal p when padj is absent", {
  x <- make_toy_scatter_enrichment()
  r <- x@results
  r$p <- r$padj
  r$padj <- NA_real_
  x@results <- r
  expect_message(p <- plot_scatter(x, "A", "B", databases = NULL), "nominal")
  expect_true(any(grepl("n = 4", layer_labels(p))))
})

test_that("the toy scatter looks the same", {
  skip_on_ci()
  skip_if_not_installed("vdiffr")
  skip_if_not_installed("correlation")
  vdiffr::expect_doppelganger("toy-scatter", scatter())
})

test_that("the subtitle draws rho as a plotmath symbol", {
  expect_true(is.call(subtitle_math("rho = 0.5, p = 0.1 | 3 terms, 1 significant")))
  expect_identical(subtitle_math("2 terms, 0 significant"), "2 terms, 0 significant")
})

test_that("the scatter prints on a plain pdf device", {
  pdf(NULL)
  on.exit(dev.off())
  expect_no_error(print(scatter()))
})

test_that("terms without a set size are drawn at one size", {
  x <- make_toy_scatter_enrichment()
  r <- x@results
  r$size <- NA_real_
  x@results <- r
  p <- scatter(x)
  expect_null(p$scales$get_scales("size"))
  expect_s3_class(p, "ggplot")
})

test_that("more than five shapes is refused", {
  x <- make_toy_scatter_enrichment()
  r <- x@results
  r$group <- rep(letters[1:6], length.out = nrow(r) / 2)
  x@results <- r
  expect_error(scatter(x, shape_by = "group"), class = "enrichVolcano_param_error")
})

test_that("more than eight colour groups fall back to the default palette", {
  x <- make_toy_scatter_enrichment()
  r <- x@results
  r$group <- rep(sprintf("g%02d", 1:12), 2)
  x@results <- r
  expect_s3_class(scatter(x, colour_by = "group")$scales$get_scales("fill"), "ScaleDiscrete")
})

test_that("an unlabelled database column draws one shape and no shape legend", {
  x <- make_toy_scatter_enrichment()
  r <- x@results
  r$database <- NA_character_
  x@results <- r
  expect_null(scatter(x)$scales$get_scales("shape"))
})

test_that("tied scores on one axis draw without rho instead of failing", {
  x <- make_toy_scatter_enrichment()
  r <- x@results
  r$score[r$contrast == "A"] <- 1.5
  r$direction[r$contrast == "A"] <- "up"
  x@results <- r
  p <- suppressWarnings(scatter(x))
  expect_false(grepl("rho", paste(deparse(p$labels$subtitle), collapse = "")))
})

test_that("contrasts with no shared terms give a clear error", {
  x <- make_toy_scatter_enrichment()
  r <- x@results
  r$term[r$contrast == "B"] <- paste0(r$term[r$contrast == "B"], "_b")
  x@results <- r
  expect_error(suppressMessages(plot_scatter(x, "A", "B")), "no terms", class = "enrichVolcano_input_error")
})

test_that("the scatter defaults to every database", {
  expect_null(formals(plot_scatter)$databases)
})
