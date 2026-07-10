source(test_path("fixtures/make_toy.R"))

test_that("resolve_volc_cols returns the four resolved names + has_padj flag", {
  v <- make_toy_volc()
  cols <- resolve_volc_cols(v, "gene", "logFC", "P.Value", "padj")
  expect_equal(cols$gene, "gene")
  expect_equal(cols$logfc, "logFC")
  expect_equal(cols$pval, "P.Value")
  expect_equal(cols$padj, "padj")
  expect_true(cols$has_padj)
})

test_that("resolve_volc_cols sets has_padj = FALSE when padj is absent", {
  v <- make_toy_volc()
  v$padj <- NULL
  cols <- resolve_volc_cols(v, "gene", "logFC", "P.Value", "padj")
  expect_false(cols$has_padj)
})

test_that("resolve_volc_cols aborts on a missing required column", {
  v <- make_toy_volc()
  expect_error(
    resolve_volc_cols(v, "gene", "no_such_lfc", "P.Value", "padj"),
    class = "enrichVolcano_column_error"
  )
})

test_that("resolve_enrich_cols accepts the toy enrichment", {
  e <- make_toy_enrich()
  cols <- resolve_enrich_cols(e, "pathway", "NES", "padj", "size", "neg_log_padj")
  expect_equal(cols$term, "pathway")
  expect_equal(cols$nes, "NES")
})

test_that("resolve_enrich_cols only requires size_col when magnitude = 'size'", {
  e <- make_toy_enrich()
  e$size <- NULL
  expect_silent(
    resolve_enrich_cols(e, "pathway", "NES", "padj", "size", "neg_log_padj")
  )
  expect_error(
    resolve_enrich_cols(e, "pathway", "NES", "padj", "size", "size"),
    class = "enrichVolcano_column_error"
  )
})

test_that("resolve_genes_col auto-detects leading_edge by default", {
  e <- make_toy_enrich()
  out <- suppressMessages(resolve_genes_col(e, NULL, NULL))
  expect_type(out, "list")
  expect_equal(out[[1]], c("G1", "G2", "G3", "G4"))
})

test_that("resolve_genes_col auto-detects leadingEdge list-column", {
  e <- make_toy_enrich_listcol()
  out <- suppressMessages(resolve_genes_col(e, NULL, NULL))
  expect_equal(out[[1]], c("G1", "G2", "G3", "G4"))
})

test_that("resolve_genes_col auto-detects core_enrichment with '/' sep", {
  e <- make_toy_enrich_cp()
  out <- suppressMessages(resolve_genes_col(e, NULL, NULL))
  expect_equal(out[[2]], c("G5", "G6"))
})

test_that("resolve_genes_col auto-detects Genes column", {
  e <- make_toy_enrich_enrichr()
  out <- suppressMessages(resolve_genes_col(e, NULL, NULL))
  expect_equal(out[[3]], c("G11", "G12"))
})

test_that("resolve_genes_col returns NULL when no tick column is present", {
  e <- make_toy_enrich()
  e$leading_edge <- NULL
  out <- suppressMessages(resolve_genes_col(e, NULL, NULL))
  expect_null(out)
})

test_that("resolve_genes_col respects a user-supplied column + sep override", {
  e <- make_toy_enrich()
  names(e)[names(e) == "leading_edge"] <- "my_genes"
  out <- suppressMessages(resolve_genes_col(e, "my_genes", ";"))
  expect_equal(out[[1]], c("G1", "G2", "G3", "G4"))
})

test_that("resolve_genes_col returns a user-supplied list-column as-is", {
  e <- make_toy_enrich_listcol()
  names(e)[names(e) == "leadingEdge"] <- "my_genes"
  out <- resolve_genes_col(e, "my_genes", NULL)
  expect_equal(out[[1]], c("G1", "G2", "G3", "G4"))
})

test_that("resolve_genes_col aborts when the override column is missing", {
  e <- make_toy_enrich()
  expect_error(
    resolve_genes_col(e, "no_such_col", ";"),
    class = "enrichVolcano_column_error"
  )
})
