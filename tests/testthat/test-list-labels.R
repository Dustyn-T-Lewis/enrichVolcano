two_contrasts <- function() enrichment(results = make_toy_results(), metadata = make_toy_metadata())

drawn_labels <- function(p) {
  unlist(lapply(ggplot2::ggplot_build(p)$data, function(d) if ("label" %in% names(d)) as.character(d$label)))
}

test_that("list_labels lists each drawn term once, with the contrasts that draw it", {
  out <- list_labels(two_contrasts())
  expect_named(out, c("database", "term", "contrasts", "clean", "short", "label"))
  expect_identical(out$term, c("HALLMARK_TOY_A", "HALLMARK_TOY_B"))
  expect_identical(out$contrasts, c("A;B", "A;B"))
  expect_identical(out$label, out$short)
})

test_that("list_labels lists the terms the ring draws", {
  x <- make_toy_ring_enrichment()
  out <- list_labels(x, databases = NULL, n_terms = 3)
  p <- suppressMessages(plot_volcano_ring(make_toy_da(), x, databases = NULL, n_terms = 3))
  drawn <- intersect(drawn_labels(p), clean_label(x@results$term))
  expect_setequal(drawn, clean_label(out$term))
})

test_that("contrast, term_threshold and n_terms narrow or widen the list", {
  out <- list_labels(two_contrasts(), contrast = "A", term_threshold = 1, n_terms = Inf)
  expect_identical(out$term, c("HALLMARK_TOY_A", "HALLMARK_TOY_B", "SET_C"))
  expect_identical(unique(out$contrasts), "A")
})

test_that("shipped terms show both styles, with line breaks written as \\n", {
  tbl <- data.frame(
    contrast = "K", score = c(2, -1), padj = c(0.01, 0.02), size = 20,
    term = c("HALLMARK_OXIDATIVE_PHOSPHORYLATION", "GOSLIM_GENERATION_OF_PRECURSOR_METABOLITES_AND_ENERGY")
  )
  x <- suppressMessages(as_enrichment(tbl, enrichment_test = "custom", score_type = "NES"))
  out <- list_labels(x)
  expect_identical(out$clean[1], "Oxidative Phosphorylation")
  expect_identical(out$short[1], "OXPHOS")
  expect_identical(out$short[2], "Precursor Metabolites\\n& Energy")
})

test_that("an edited write_labels file relabels the ring", {
  x <- make_toy_ring_enrichment()
  file <- withr::local_tempfile(fileext = ".csv")
  expect_invisible(out <- write_labels(x, file, databases = NULL))
  expect_identical(out, file)
  edited <- utils::read.csv(file)
  edited$label[1] <- "Edited\\nname"
  p <- suppressMessages(plot_volcano_ring(make_toy_da(), x, databases = NULL, labels = edited))
  expect_true("Edited\nname" %in% drawn_labels(p))
})

test_that("list_labels needs an enrichment object and a contrast it holds", {
  expect_error(list_labels(data.frame()), class = "enrichVolcano_input_error")
  expect_error(list_labels(two_contrasts(), contrast = "Z"), class = "enrichVolcano_input_error")
})
