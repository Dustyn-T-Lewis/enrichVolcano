test_that("ev_collapse jaccard drops near-duplicate pathway", {
  enrich <- tibble::tibble(
    contrast = "c", database = "db", mode = "fgsea",
    pathway = c("P1", "P2"),
    pval = c(0.001, 0.01), padj = c(0.01, 0.05),
    NES = c(2, 1.5), size = c(10, 10), direction = c("up", "up"),
    leading_edge = c("A;B;C;D;E;F;G;H;I;J", "A;B;C;D;E;F;G;H;I;K")
  )
  result <- ev_collapse(enrich, method = "jaccard", cutoff = 0.5,
                       scope = "within_db", keep_by = "padj")
  expect_equal(sum(result$dedup_kept), 1)
  expect_true(result$dedup_kept[result$pathway == "P1"])
})

test_that("ev_collapse keeps both when Jaccard < cutoff", {
  enrich <- tibble::tibble(
    contrast = "c", database = "db", mode = "fgsea",
    pathway = c("P1", "P2"),
    pval = c(0.001, 0.01), padj = c(0.01, 0.05),
    NES = c(2, 1.5), size = c(10, 10), direction = c("up", "up"),
    leading_edge = c("A;B;C;D;E", "F;G;H;I;J")
  )
  result <- ev_collapse(enrich, method = "jaccard", cutoff = 0.5,
                       scope = "within_db", keep_by = "padj")
  expect_equal(sum(result$dedup_kept), 2)
})

test_that("ev_collapse cross_db scope merges across DBs", {
  enrich <- tibble::tibble(
    contrast = "c", database = c("db1", "db2"), mode = "fgsea",
    pathway = c("P1", "P2"),
    pval = c(0.001, 0.01), padj = c(0.01, 0.05),
    NES = c(2, 1.5), size = c(10, 10), direction = c("up", "up"),
    leading_edge = c("A;B;C;D;E", "A;B;C;D;F")
  )
  result <- ev_collapse(enrich, method = "jaccard", cutoff = 0.5,
                       scope = "cross_db", keep_by = "padj")
  expect_equal(sum(result$dedup_kept), 1)
})

test_that("ev_collapse method='both' applies jaccard then collapsePathways", {
  data <- make_ranked_input()
  paths <- make_mini_pathways()
  enrich <- suppressWarnings(ev_enrich(
    data, contrast = "ctr",
    databases = list(test = paths),
    enrich_mode = "fgsea", rank_by = "signed_p",
    nperm = 100, min_size = 5, max_size = 100
  ))
  result <- ev_collapse(enrich, method = "both",
                       collapse_pval_threshold = 0.05)
  expect_true("dedup_kept" %in% colnames(result))
})
