source(test_path("fixtures/make_toy.R"))

test_that("volcano_ring returns a ggplot from toy frames", {
  p <- suppressMessages(volcano_ring(make_toy_volc(), make_toy_enrich()))
  expect_s3_class(p, "ggplot")
})

test_that("volcano_ring builds without a tick column (silent ticks)", {
  e <- make_toy_enrich()
  e$leading_edge <- NULL
  p <- suppressMessages(volcano_ring(make_toy_volc(), e))
  expect_s3_class(p, "ggplot")
})

test_that("volcano_ring accepts user-supplied column names", {
  v <- make_toy_volc()
  names(v) <- sub("^P\\.Value$", "pval", names(v))
  names(v) <- sub("^logFC$", "lfc", names(v))
  names(v) <- sub("^gene$", "id", names(v))
  p <- suppressMessages(volcano_ring(
    v, make_toy_enrich(),
    gene_col = "id", logfc_col = "lfc", pval_col = "pval"
  ))
  expect_s3_class(p, "ggplot")
})

test_that("volcano_ring respects the magnitude = 'size' switch", {
  p <- suppressMessages(volcano_ring(
    make_toy_volc(), make_toy_enrich(),
    magnitude = "size"
  ))
  expect_s3_class(p, "ggplot")
})

test_that("volcano_ring works with the okabe palette", {
  p <- suppressMessages(volcano_ring(
    make_toy_volc(), make_toy_enrich(),
    theme = volcano_ring_theme(palette = "okabe")
  ))
  expect_s3_class(p, "ggplot")
})

test_that("volcano_ring runs against the bundled YvO fixture", {
  da <- read.csv(system.file("extdata", "examples", "yvo_da.csv",
    package = "enrichVolcano"
  ))
  en <- read.csv(system.file("extdata", "examples", "yvo_enrichment.csv",
    package = "enrichVolcano"
  ))
  da <- da[da$contrast == da$contrast[1], ]
  en <- en[en$contrast == en$contrast[1], ]
  names(da)[names(da) == "adj.P.Val"] <- "padj"
  p <- suppressWarnings(suppressMessages(volcano_ring(da, en)))
  expect_s3_class(p, "ggplot")
})

test_that("arc_order = 'nes' reorders arcs by absolute NES within a half", {
  # P: most significant, weakest NES; Q: least significant, strongest NES.
  e <- data.frame(
    pathway = c("P", "Q"), NES = c(1.2, 3.0),
    padj = c(0.001, 0.04), size = c(20L, 20L),
    stringsAsFactors = FALSE
  )
  gl <- replicate(nrow(e), character(0), simplify = FALSE)
  mag <- -log10(e$padj)
  by_padj <- ev_ring_geometry(e, "pathway", "padj", "NES", mag, gl,
    order_by = "padj"
  )
  by_nes <- ev_ring_geometry(e, "pathway", "padj", "NES", mag, gl,
    order_by = "nes"
  )
  expect_equal(by_padj$pathway, c("P", "Q"))
  expect_equal(by_nes$pathway, c("Q", "P"))
})

test_that("arc_height_range sets the shortest and tallest arc extent", {
  e <- make_toy_enrich()
  gl <- replicate(nrow(e), character(0), simplify = FALSE)
  mag <- -log10(pmax(e$padj, 1e-300))
  wide <- ev_ring_geometry(e, "pathway", "padj", "NES", mag, gl,
    min_height = 0.1, max_height = 3
  )
  heights <- wide$arc_r1_var - ev_ring_r_outer
  expect_lte(max(heights), 3)
  expect_gte(min(heights), 0.1 - 1e-8)
})

test_that("show_counts = FALSE drops the up/down count badges", {
  v <- make_toy_volc()
  e <- make_toy_enrich()
  labels_of <- function(p) {
    unlist(lapply(ggplot2::ggplot_build(p)$data, function(d) {
      if ("label" %in% names(d)) as.character(d$label) else character(0)
    }))
  }
  with_counts <- labels_of(suppressMessages(volcano_ring(v, e, show_counts = TRUE)))
  without <- labels_of(suppressMessages(volcano_ring(v, e, show_counts = FALSE)))
  n_up <- as.character(sum(v$logFC > 0))
  expect_true(n_up %in% with_counts)
  expect_false(n_up %in% without)
})

test_that("volcano_ring_grid composes a 2-panel layout from named lists", {
  v1 <- make_toy_volc(seed = 1L)
  v2 <- make_toy_volc(seed = 2L)
  e1 <- make_toy_enrich(seed = 1L)
  e2 <- make_toy_enrich(seed = 2L)
  g <- suppressMessages(volcano_ring_grid(
    list(A = v1, B = v2),
    list(A = e1, B = e2)
  ))
  expect_s3_class(g, "volcano_ring_grid")
  expect_named(g, c("plot", "data"))
  expect_equal(names(g$data), c("A", "B"))
})

test_that("volcano_ring_grid splits a single data.frame by contrast", {
  d <- rbind(
    transform(make_toy_volc(seed = 1L), contrast = "A"),
    transform(make_toy_volc(seed = 2L), contrast = "B")
  )
  e <- rbind(
    transform(make_toy_enrich(seed = 1L), contrast = "A"),
    transform(make_toy_enrich(seed = 2L), contrast = "B")
  )
  g <- suppressMessages(volcano_ring_grid(d, e))
  expect_equal(names(g$data), c("A", "B"))
})

test_that("volcano_ring_grid aborts when a contrast is missing in volc_dfs", {
  expect_error(
    volcano_ring_grid(
      list(A = make_toy_volc(seed = 1L)),
      list(A = make_toy_enrich(seed = 1L), B = make_toy_enrich(seed = 2L)),
      contrasts = c("A", "B")
    ),
    class = "enrichVolcano_input_error"
  )
})

test_that("volcano_ring_grid aborts when contrasts cannot be inferred", {
  expect_error(
    volcano_ring_grid(list(make_toy_volc()), list(make_toy_enrich())),
    class = "enrichVolcano_input_error"
  )
})

test_that("volcano_ring_grid aborts when split-by-contrast has no column", {
  v <- make_toy_volc()
  v$contrast <- NULL
  expect_error(
    volcano_ring_grid(v, make_toy_enrich()),
    class = "enrichVolcano_input_error"
  )
})

test_that("print.volcano_ring_grid returns the object invisibly", {
  d <- rbind(
    transform(make_toy_volc(seed = 1L), contrast = "A"),
    transform(make_toy_volc(seed = 2L), contrast = "B")
  )
  e <- rbind(
    transform(make_toy_enrich(seed = 1L), contrast = "A"),
    transform(make_toy_enrich(seed = 2L), contrast = "B")
  )
  g <- suppressMessages(volcano_ring_grid(d, e))
  pdf(NULL)
  on.exit(dev.off(), add = TRUE)
  expect_invisible(print(g))
})
