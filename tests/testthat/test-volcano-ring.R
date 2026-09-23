source(test_path("fixtures/make_toy.R"))

# Every toy term on the ring, no database filter: the pre-1.0 behaviour.
ring <- function(v = make_toy_volc(), x = make_toy_ring_enrichment(), ...) {
  suppressMessages(volcano_ring(v, x, databases = NULL, term_threshold = 1, n_terms = Inf, ...))
}

fill_title <- function(p) p$scales$get_scales("fill")$name

test_that("volcano_ring returns a ggplot from an enrichment", {
  expect_s3_class(ring(), "ggplot")
})

test_that("volcano_ring builds without leading edges (no ticks)", {
  x <- make_toy_ring_enrichment()
  r <- x@results
  r$leading_edge <- rep(list(character(0)), nrow(r))
  x@results <- r
  expect_s3_class(ring(x = x), "ggplot")
})

test_that("volcano_ring accepts user-supplied volcano column names", {
  v <- make_toy_volc()
  names(v) <- sub("^P\\.Value$", "pval", names(v))
  names(v) <- sub("^logFC$", "lfc", names(v))
  names(v) <- sub("^gene$", "id", names(v))
  expect_s3_class(ring(v, gene_col = "id", logfc_col = "lfc", pval_col = "pval"), "ggplot")
})

test_that("volcano_ring respects the magnitude = 'size' switch", {
  expect_s3_class(ring(magnitude = "size"), "ggplot")
})

test_that("volcano_ring works with the okabe palette", {
  expect_s3_class(ring(theme = volcano_ring_theme(palette = "okabe")), "ggplot")
})

test_that("volcano_ring runs on the bundled YvO example with default selection", {
  da <- read.csv(system.file("extdata", "examples", "yvo_da.csv.gz", package = "enrichVolcano"))
  ex <- suppressMessages(as_enrichment(read.csv(
    system.file("extdata", "examples", "yvo_fgsea.csv.gz", package = "enrichVolcano")
  )))
  da <- da[da$contrast == "Training_Young", ]
  names(da)[names(da) == "adj.P.Val"] <- "padj"
  p <- suppressWarnings(suppressMessages(volcano_ring(da, ex, contrast = "Training_Young")))
  expect_s3_class(p, "ggplot")
  expect_identical(fill_title(p), "NES")
})

test_that("the legend is titled with the score type", {
  v <- make_toy_volc()
  fry <- utils::read.csv(test_path("fixtures", "limma_fry.csv"), row.names = 1)
  x <- suppressMessages(as_enrichment(list(toy = fry)))
  p <- ring(v, x)
  expect_identical(fill_title(p), "signed -log10(FDR)")
})

test_that("a multi-contrast object needs a contrast, and a real one", {
  x <- make_toy_ring_enrichment(c("A", "B"))
  expect_error(volcano_ring(make_toy_volc(), x), class = "enrichVolcano_input_error")
  expect_error(volcano_ring(make_toy_volc(), x, contrast = "C"), "C", class = "enrichVolcano_input_error")
  expect_s3_class(ring(x = x, contrast = "B"), "ggplot")
})

test_that("volcano_ring refuses a plain enrichment table", {
  expect_error(
    volcano_ring(make_toy_volc(), make_toy_enrich()),
    "as_enrichment",
    class = "enrichVolcano_input_error"
  )
})

test_that("default selection drops non-significant terms from the ring", {
  labels_in <- function(p) {
    unlist(lapply(ggplot2::ggplot_build(p)$data, function(d) {
      if ("label" %in% names(d)) as.character(d$label) else character(0)
    }))
  }
  all_terms <- labels_in(ring())
  selected <- labels_in(suppressMessages(volcano_ring(make_toy_volc(), make_toy_ring_enrichment(), databases = NULL)))
  expect_true(any(grepl("Toy E", all_terms)))
  expect_false(any(grepl("Toy E", selected)))
})

test_that("a ring with no terms left still draws and says why", {
  expect_message(
    p <- volcano_ring(make_toy_volc(), make_toy_ring_enrichment(), databases = NULL, term_threshold = 1e-9),
    class = "enrichVolcano_empty_ring"
  )
  expect_s3_class(p, "ggplot")
})

test_that("disc_color draws a tinted central disc", {
  expect_s3_class(ring(disc_color = "grey70"), "ggplot")
})

test_that("volc_sig_col drives point significance independently of padj", {
  expect_s3_class(ring(volc_sig_col = "pi_eq2"), "ggplot")
})

test_that("x_scale and y_scale compress the volcano point cloud", {
  point_span <- function(p, axis) {
    pts <- Filter(function(d) "shape" %in% names(d), ggplot2::ggplot_build(p)$data)
    diff(range(unlist(lapply(pts, `[[`, axis)), na.rm = TRUE))
  }
  full <- ring()
  narrow <- ring(x_scale = 0.5)
  short <- ring(y_scale = 0.5)
  expect_lt(point_span(narrow, "x"), point_span(full, "x"))
  expect_lt(point_span(short, "y"), point_span(full, "y"))
})

test_that("arc_order = 'nes' reorders arcs by absolute score within a half", {
  # P: most significant, weakest score; Q: least significant, strongest score.
  e <- data.frame(
    term = c("P", "Q"), score = c(1.2, 3.0),
    padj = c(0.001, 0.04), size = c(20L, 20L),
    stringsAsFactors = FALSE
  )
  gl <- replicate(nrow(e), character(0), simplify = FALSE)
  mag <- -log10(e$padj)
  by_padj <- ev_ring_geometry(e, "term", "padj", "score", mag, gl, order_by = "padj")
  by_nes <- ev_ring_geometry(e, "term", "padj", "score", mag, gl, order_by = "nes")
  expect_equal(by_padj$term, c("P", "Q"))
  expect_equal(by_nes$term, c("Q", "P"))
})

test_that("arc_height_range sets the shortest and tallest arc extent", {
  e <- make_toy_ring_enrichment()@results
  mag <- -log10(pmax(e$padj, 1e-300))
  wide <- ev_ring_geometry(e, "term", "padj", "score", mag, e$leading_edge,
    min_height = 0.1, max_height = 3
  )
  heights <- wide$arc_r1_var - ev_ring_r_outer
  expect_lte(max(heights), 3)
  expect_gte(min(heights), 0.1 - 1e-8)
})

test_that("show_counts = FALSE drops the up/down count badges", {
  v <- make_toy_volc()
  labels_of <- function(p) {
    unlist(lapply(ggplot2::ggplot_build(p)$data, function(d) {
      if ("label" %in% names(d)) as.character(d$label) else character(0)
    }))
  }
  with_counts <- labels_of(ring(v, show_counts = TRUE))
  without <- labels_of(ring(v, show_counts = FALSE))
  n_up <- as.character(sum(v$logFC > 0))
  expect_true(n_up %in% with_counts)
  expect_false(n_up %in% without)
})

grid <- function(volc, x = make_toy_ring_enrichment(c("A", "B")), ...) {
  suppressMessages(volcano_ring_grid(volc, x, databases = NULL, term_threshold = 1, n_terms = Inf, ...))
}

test_that("volcano_ring_grid composes a 2-panel layout from named lists", {
  g <- grid(list(A = make_toy_volc(seed = 1L), B = make_toy_volc(seed = 2L)))
  expect_s3_class(g, "volcano_ring_grid")
  expect_named(g, c("plot", "data"))
  expect_s3_class(g$plot, "patchwork")
  expect_named(g$data, c("A", "B"))
})

test_that("volcano_ring_grid splits a single data.frame by contrast", {
  d <- rbind(
    transform(make_toy_volc(seed = 1L), contrast = "A"),
    transform(make_toy_volc(seed = 2L), contrast = "B")
  )
  g <- grid(d)
  expect_named(g$data, c("A", "B"))
})

test_that("volcano_ring_grid aborts when a contrast is missing on either side", {
  expect_error(
    grid(list(A = make_toy_volc()), contrasts = c("A", "B")),
    class = "enrichVolcano_input_error"
  )
  expect_error(
    grid(list(A = make_toy_volc(), C = make_toy_volc())),
    "C",
    class = "enrichVolcano_input_error"
  )
})

test_that("volcano_ring_grid aborts when contrasts cannot be inferred", {
  expect_error(grid(list(make_toy_volc())), class = "enrichVolcano_input_error")
})

test_that("volcano_ring_grid aborts when split-by-contrast has no column", {
  v <- make_toy_volc()
  v$contrast <- NULL
  expect_error(grid(v), class = "enrichVolcano_input_error")
})

test_that("volcano_ring_grid refuses a plain enrichment table", {
  expect_error(
    volcano_ring_grid(list(A = make_toy_volc()), make_toy_enrich()),
    class = "enrichVolcano_input_error"
  )
})

test_that("print.volcano_ring_grid returns the object invisibly", {
  g <- grid(list(A = make_toy_volc(seed = 1L), B = make_toy_volc(seed = 2L)))
  pdf(NULL)
  on.exit(dev.off())
  expect_invisible(print(g))
})

test_that("a fry ring's fill scale spans its scores instead of squishing at 3", {
  fry <- utils::read.csv(test_path("fixtures", "limma_fry.csv"), row.names = 1)
  x <- suppressMessages(as_enrichment(list(toy = fry)))
  p <- ring(x = x)
  limits <- p$scales$get_scales("fill")$limits
  expect_equal(max(limits), max(abs(x@results$score[x@results$padj < 0.05])))
})

test_that("the ring holds at most n_terms terms in total", {
  arcs <- function(p) {
    geoms <- vapply(p$layers, function(l) class(l$geom)[1], character(1))
    length(unique(ggplot2::layer_data(p, which(geoms == "GeomArcBar")[1])$group))
  }
  x <- make_toy_ring_enrichment()
  expect_identical(arcs(suppressMessages(volcano_ring(make_toy_volc(), x, databases = NULL, n_terms = 3))), 3L)
  expect_identical(formals(volcano_ring)$n_terms, 12)
})

test_that("results without set sizes draw arcs by padj instead of failing", {
  ora <- data.frame(term = c("S1", "S2", "S3"), padj = c(0.001, 0.01, 0.02), direction = c("up", "down", "up"))
  x <- as_enrichment(list(toy = ora), enrichment_test = "ora")
  expect_s3_class(suppressMessages(volcano_ring(make_toy_volc(), x)), "ggplot")
  expect_error(
    suppressMessages(volcano_ring(make_toy_volc(), x, magnitude = "size")),
    class = "enrichVolcano_param_error"
  )
})

test_that("the default draws from every database", {
  expect_null(formals(volcano_ring)$databases)
  x <- make_toy_ring_enrichment()
  r <- x@results
  r$database <- "KEGG"
  x@results <- r
  expect_s3_class(suppressMessages(volcano_ring(make_toy_volc(), x)), "ggplot")
})

test_that("a hand-picked term found in two databases is drawn once", {
  x <- make_toy_ring_enrichment()
  r <- x@results
  twin <- r[1, ]
  twin$database <- "Other"
  twin$padj <- 0.3
  x@results <- rbind(r, twin)
  picked <- ring_terms(x@results, 0.05, 12, terms = r$term[1])
  expect_identical(nrow(picked), 1L)
  expect_identical(picked$padj, r$padj[1])
})

test_that("grid subtitles can be named by contrast or given in order", {
  volc <- list(A = make_toy_volc(seed = 1L), B = make_toy_volc(seed = 2L))
  subtitles_of <- function(g) c(g$plot[[1]]$labels$subtitle, g$plot[[2]]$labels$subtitle)
  expect_identical(subtitles_of(grid(volc, subtitles = c(B = "second", A = "first"))), c("first", "second"))
  expect_identical(subtitles_of(grid(volc, subtitles = c("first", "second"))), c("first", "second"))
})
