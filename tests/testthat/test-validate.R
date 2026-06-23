source(test_path("fixtures/make_toy.R"))

vcols_default <- function() {
  list(gene = "gene", logfc = "logFC", pval = "P.Value", padj = "padj", has_padj = TRUE)
}
ecols_default <- function() {
  list(term = "pathway", nes = "NES", padj = "padj", size = "size")
}

test_that("validate_volc_df accepts a tidy frame", {
  expect_silent(validate_volc_df(make_toy_volc(), vcols_default()))
})

test_that("validate_volc_df rejects non-data.frame input", {
  expect_error(
    validate_volc_df(list(), vcols_default()),
    class = "enrichVolcano_input_error"
  )
})

test_that("validate_volc_df rejects an empty frame", {
  v <- make_toy_volc()[0, ]
  expect_error(
    validate_volc_df(v, vcols_default()),
    class = "enrichVolcano_input_error"
  )
})

test_that("validate_volc_df rejects non-numeric logFC", {
  v <- make_toy_volc()
  v$logFC <- as.character(v$logFC)
  expect_error(
    validate_volc_df(v, vcols_default()),
    class = "enrichVolcano_column_error"
  )
})

test_that("validate_volc_df rejects P-values outside [0, 1]", {
  v <- make_toy_volc()
  v$P.Value[1] <- 1.5
  expect_error(
    validate_volc_df(v, vcols_default()),
    class = "enrichVolcano_data_error"
  )
})

test_that("validate_enrich_df accepts a tidy frame", {
  expect_silent(validate_enrich_df(make_toy_enrich(), ecols_default()))
})

test_that("validate_enrich_df rejects all-NA NES", {
  e <- make_toy_enrich()
  e$NES <- NA_real_
  expect_error(
    validate_enrich_df(e, ecols_default()),
    class = "enrichVolcano_data_error"
  )
})

test_that("validate_enrich_df warns on duplicate terms", {
  e <- make_toy_enrich()
  e$pathway[2] <- e$pathway[1]
  expect_warning(
    validate_enrich_df(e, ecols_default()),
    class = "enrichVolcano_input_error"
  )
})

test_that("validate_enrich_df rejects non-data.frame input", {
  expect_error(
    validate_enrich_df(list(), ecols_default()),
    class = "enrichVolcano_input_error"
  )
})

test_that("validate_enrich_df rejects an empty frame", {
  e <- make_toy_enrich()[0, ]
  expect_error(
    validate_enrich_df(e, ecols_default()),
    class = "enrichVolcano_input_error"
  )
})

test_that("validate_enrich_df rejects non-numeric NES", {
  e <- make_toy_enrich()
  e$NES <- as.character(e$NES)
  expect_error(
    validate_enrich_df(e, ecols_default()),
    class = "enrichVolcano_column_error"
  )
})

test_that("validate_enrich_df rejects non-numeric padj", {
  e <- make_toy_enrich()
  e$padj <- as.character(e$padj)
  expect_error(
    validate_enrich_df(e, ecols_default()),
    class = "enrichVolcano_column_error"
  )
})

test_that("validate_enrich_df rejects padj outside [0, 1]", {
  e <- make_toy_enrich()
  e$padj[1] <- 1.5
  expect_error(
    validate_enrich_df(e, ecols_default()),
    class = "enrichVolcano_data_error"
  )
})

test_that("dedup_by_term passes through when no duplicates exist", {
  e <- make_toy_enrich()
  expect_identical(dedup_by_term(e, ecols_default()), e)
})

test_that("dedup_by_term keeps the lowest-padj row per term", {
  e <- make_toy_enrich()
  e2 <- rbind(e, transform(e[1, ], padj = 0.5))
  out <- dedup_by_term(e2, ecols_default())
  expect_equal(nrow(out), 5L)
  expect_equal(out$padj[out$pathway == e$pathway[1]], e$padj[1])
})
