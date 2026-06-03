test_that("ev_enrich fgsea mode runs and returns expected schema", {
  data <- make_ranked_input()
  pathways <- make_mini_pathways()
  result <- suppressWarnings(ev_enrich(
    data, contrast = "ctr",
    databases = list(hallmark = pathways),
    enrich_mode = "fgsea",
    rank_by = "signed_p",
    nperm = 1000, seed = 42, min_size = 5, max_size = 100
  ))
  expect_s3_class(result, "tbl_df")
  expected_cols <- c("contrast", "database", "pathway", "pval", "padj",
                     "NES", "size", "leading_edge", "mode", "direction")
  expect_true(all(expected_cols %in% colnames(result)))
  expect_equal(unique(result$mode), "fgsea")
})

test_that("ev_enrich ORA mode carries fold_enrichment on its rows", {
  data <- make_ranked_input()
  pathways <- make_mini_pathways()
  result <- suppressWarnings(ev_enrich(
    data, contrast = "ctr",
    databases = list(hallmark = pathways),
    enrich_mode = "ora",
    min_size = 5, max_size = 100
  ))
  expect_s3_class(result, "tbl_df")
  expect_equal(unique(result$mode), "ora")
  expect_true("fold_enrichment" %in% colnames(result))
  expect_true(is.numeric(result$fold_enrichment))
  expect_true(all(result$fold_enrichment >= 0, na.rm = TRUE))
  # ORA has no GSEA enrichment score
  expect_true(all(is.na(result$NES)))
})

test_that("ev_enrich both modes union columns: fold_enrichment for ORA, ES for fgsea", {
  data <- make_ranked_input()
  pathways <- make_mini_pathways()
  result <- suppressWarnings(ev_enrich(
    data, contrast = "ctr",
    databases = list(hallmark = pathways),
    enrich_mode = c("fgsea", "ora"),
    rank_by = "signed_p", nperm = 1000, seed = 42,
    min_size = 5, max_size = 100
  ))
  expect_true(all(c("fold_enrichment", "ES") %in% colnames(result)))
  expect_true(all(is.na(result$fold_enrichment[result$mode == "fgsea"])))
  expect_false(any(is.na(result$fold_enrichment[result$mode == "ora"])))
})

test_that("ev_enrich returns empty tibble (not error) when no pathways match", {
  data <- make_ranked_input()
  empty_paths <- list(NEVER_MATCH = c("ZZ1", "ZZ2"))
  result <- suppressWarnings(ev_enrich(
    data, contrast = "ctr",
    databases = list(test = empty_paths),
    enrich_mode = "fgsea",
    rank_by = "signed_p",
    nperm = 100, min_size = 5
  ))
  expect_s3_class(result, "tbl_df")
  expect_equal(nrow(result), 0)
})

test_that("include_terms restricts pathway universe before testing", {
  data <- make_ranked_input()
  paths <- list(
    HALLMARK_GLYCOLYSIS = paste0("G", 1:20),
    HALLMARK_HYPOXIA    = paste0("G", 21:40),
    MITO_TRANSPORT      = paste0("G", 41:60)
  )
  out_all <- ev_enrich(data, contrast = "ctr",
                       databases = list(test = paths),
                       enrich_mode = "fgsea",
                       min_size = 5, max_size = 100,
                       include_terms = NULL, exclude_terms = NULL)
  out_mito <- ev_enrich(data, contrast = "ctr",
                        databases = list(test = paths),
                        enrich_mode = "fgsea",
                        min_size = 5, max_size = 100,
                        include_terms = "^MITO",
                        exclude_terms = NULL)
  expect_setequal(unique(out_mito$pathway), "MITO_TRANSPORT")
  expect_gt(nrow(out_all), nrow(out_mito))
})
