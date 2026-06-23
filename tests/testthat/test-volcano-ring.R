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
