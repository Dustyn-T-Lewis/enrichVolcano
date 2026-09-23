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
