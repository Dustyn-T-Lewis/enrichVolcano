test_that("ev_enrich fgsea mode runs and returns expected schema", {
  data <- make_ranked_input()
  pathways <- make_mini_pathways()
  result <- suppressWarnings(ev_enrich(
    data, contrast = "ctr",
    databases = list(hallmark = pathways),
    enrich_mode = "fgsea",
    rank_by = "pi_eq2",
    nperm = 1000, seed = 42, min_size = 5, max_size = 100
  ))
  expect_s3_class(result, "tbl_df")
  expected_cols <- c("contrast", "database", "pathway", "pval", "padj",
                     "NES", "size", "leading_edge", "mode", "direction")
  expect_true(all(expected_cols %in% colnames(result)))
  expect_equal(unique(result$mode), "fgsea")
})

test_that("ev_enrich returns empty tibble (not error) when no pathways match", {
  data <- make_ranked_input()
  empty_paths <- list(NEVER_MATCH = c("ZZ1", "ZZ2"))
  result <- suppressWarnings(ev_enrich(
    data, contrast = "ctr",
    databases = list(test = empty_paths),
    enrich_mode = "fgsea",
    rank_by = "pi_eq2",
    nperm = 100, min_size = 5
  ))
  expect_s3_class(result, "tbl_df")
  expect_equal(nrow(result), 0)
})
