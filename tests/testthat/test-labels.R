make_df <- function() {
  tibble::tibble(
    gene = paste0("G", 1:6),
    symbol = c("AAA", "BBB", NA, "DDD", "EEE", "FFF"),
    uniprot = paste0("P0000", 1:6),
    logFC = c(3, -2, 1.5, -3, 0.1, 2.5),
    P.Value = c(1e-5, 1e-4, 1e-3, 1e-6, 0.2, 1e-2),
    adj.P.Val = c(1e-4, 1e-3, 1e-2, 1e-5, 0.5, 0.05),
    contrast = "C1"
  )
}

test_that("ev_label_text prefers symbol, falls back to accession then gene", {
  df <- make_df()
  expect_identical(ev_label_text(df),
                   c("AAA", "BBB", "P00003", "DDD", "EEE", "FFF"))
})

test_that("mode none labels nothing", {
  df <- make_df()
  out <- ev_select_labels(df, mode = "none", n = 10, rank_by = "significance",
                          genes = NULL, p_col = "adj.P.Val",
                          p_threshold = 0.05, logfc_threshold = log2(1.5))
  expect_equal(nrow(out), 0)
})

test_that("top_per_direction returns up to n per direction, significant only", {
  df <- make_df()
  out <- ev_select_labels(df, mode = "top_per_direction", n = 1,
                          rank_by = "significance", genes = NULL,
                          p_col = "adj.P.Val", p_threshold = 0.05,
                          logfc_threshold = log2(1.5))
  # most-significant up = G1 (adj 1e-4), most-significant down = G4 (adj 1e-5)
  expect_setequal(out$gene, c("G1", "G4"))
})

test_that("top_total ranks across both directions", {
  df <- make_df()
  out <- ev_select_labels(df, mode = "top_total", n = 2,
                          rank_by = "significance", genes = NULL,
                          p_col = "adj.P.Val", p_threshold = 0.05,
                          logfc_threshold = log2(1.5))
  expect_setequal(out$gene, c("G4", "G1"))   # smallest adj.P.Val overall
})

test_that("rank_by logfc selects biggest movers", {
  df <- make_df()
  out <- ev_select_labels(df, mode = "top_total", n = 2,
                          rank_by = "logfc", genes = NULL,
                          p_col = "adj.P.Val", p_threshold = 0.05,
                          logfc_threshold = log2(1.5))
  expect_setequal(out$gene, c("G1", "G4"))   # |logFC| 3 and 3
})

test_that("all_significant returns every threshold-passing point", {
  df <- make_df()
  out <- ev_select_labels(df, mode = "all_significant", n = 10,
                          rank_by = "significance", genes = NULL,
                          p_col = "adj.P.Val", p_threshold = 0.05,
                          logfc_threshold = log2(1.5))
  # G5 fails (ns), others pass logFC + p
  expect_setequal(out$gene, c("G1", "G2", "G3", "G4", "G6"))
})

test_that("explicit matches symbol, accession, or gene", {
  df <- make_df()
  out <- ev_select_labels(df, mode = "explicit", n = 10,
                          rank_by = "significance",
                          genes = c("AAA", "P00004"), p_col = "adj.P.Val",
                          p_threshold = 0.05, logfc_threshold = log2(1.5))
  expect_setequal(out$gene, c("G1", "G4"))
})
