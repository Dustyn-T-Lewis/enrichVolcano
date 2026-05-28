test_that("ev_collapse jaccard drops near-duplicate pathway", {
  enrich <- tibble::tibble(
    contrast = "c", database = "db", mode = "fgsea",
    pathway = c("P1", "P2"),
    pval = c(0.001, 0.01), padj = c(0.01, 0.05),
    NES = c(2, 1.5), size = c(10, 10), direction = c("up", "up"),
    leading_edge = c("A;B;C;D;E;F;G;H;I;J", "A;B;C;D;E;F;G;H;I;K")
  )
  # sig_threshold = NA preserves original "dedup all rows" intent:
  # fixture has P2 padj = 0.05 which is NOT < 0.05 (gated out by default).
  result <- ev_collapse(enrich, method = "jaccard", cutoff = 0.5,
                       scope = "within_db", keep_by = "padj",
                       sig_threshold = NA)
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
  # sig_threshold = NA preserves original "dedup all rows" intent:
  # fixture has P2 padj = 0.05 which is NOT < 0.05 (gated out by default).
  result <- ev_collapse(enrich, method = "jaccard", cutoff = 0.5,
                       scope = "within_db", keep_by = "padj",
                       sig_threshold = NA)
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
  # sig_threshold = NA preserves original "dedup all rows" intent:
  # fixture has P2 padj = 0.05 which is NOT < 0.05 (gated out by default).
  result <- ev_collapse(enrich, method = "jaccard", cutoff = 0.5,
                       scope = "cross_db", keep_by = "padj",
                       sig_threshold = NA)
  expect_equal(sum(result$dedup_kept), 1)
})

test_that("ev_collapse method='jaccard_then_collapse' applies jaccard then collapsePathways", {
  data <- make_ranked_input()
  paths <- make_mini_pathways()
  enrich <- suppressWarnings(ev_enrich(
    data, contrast = "ctr",
    databases = list(test = paths),
    enrich_mode = "fgsea", rank_by = "signed_p",
    nperm = 100, min_size = 5, max_size = 100
  ))
  result <- ev_collapse(enrich, method = "jaccard_then_collapse",
                       collapse_pval_threshold = 0.05)
  expect_true("dedup_kept" %in% colnames(result))
})

test_that('method = "both" is deprecated alias for jaccard_then_collapse', {
  enrich <- tibble::tibble(
    contrast = "c", database = "db", mode = "fgsea",
    pathway = c("P1", "P2"),
    pval = c(0.001, 0.01), padj = c(0.01, 0.05),
    NES = c(2, 1.5), size = c(10, 10), direction = c("up", "up"),
    leading_edge = c("A;B;C;D;E;F;G;H;I;J", "A;B;C;D;E;F;G;H;I;K")
  )
  rlang::local_options(lifecycle_verbosity = "warning")
  expect_warning(
    out_both <- ev_collapse(enrich, method = "both", cutoff = 0.5,
                            scope = "within_db", keep_by = "padj",
                            sig_threshold = NA),
    class = "lifecycle_warning_deprecated"
  )
  out_new <- ev_collapse(enrich, method = "jaccard_then_collapse",
                         cutoff = 0.5, scope = "within_db",
                         keep_by = "padj", sig_threshold = NA)
  expect_identical(out_both$dedup_kept, out_new$dedup_kept)
})

test_that('sig_threshold gates dedup to significant pathways only', {
  enrich <- tibble::tibble(
    contrast = "c", database = "db", mode = "fgsea",
    pathway = c("P1", "P2", "P3", "P4"),
    pval = c(0.001, 0.001, 0.2, 0.3),
    padj = c(0.01, 0.02, 0.20, 0.40),
    NES = c(2, 1.9, 1.5, 1.3),
    size = c(10, 10, 10, 10),
    direction = "up",
    leading_edge = c(
      "A;B;C;D;E;F;G;H;I;J",
      "A;B;C;D;E;F;G;H;I;K",
      "A;B;C;D;E;F;G;H;I;L",
      "A;B;C;D;E;F;G;H;I;M"
    )
  )
  out <- ev_collapse(enrich, method = "jaccard", cutoff = 0.5,
                     scope = "within_db", sig_threshold = 0.05,
                     keep_by = "padj")
  expect_true(out$dedup_kept[out$pathway == "P3"])
  expect_true(out$dedup_kept[out$pathway == "P4"])
  expect_true(out$dedup_kept[out$pathway == "P1"])
  expect_false(out$dedup_kept[out$pathway == "P2"])
})

test_that('sig_threshold = NA disables gating (dedup all rows)', {
  enrich <- tibble::tibble(
    contrast = "c", database = "db", mode = "fgsea",
    pathway = c("P1", "P2"),
    pval = c(0.001, 0.001), padj = c(0.20, 0.30),
    NES = c(2, 1.5), size = c(10, 10), direction = "up",
    leading_edge = c("A;B;C;D;E;F;G;H;I;J", "A;B;C;D;E;F;G;H;I;K")
  )
  out <- ev_collapse(enrich, method = "jaccard", cutoff = 0.5,
                     scope = "within_db", sig_threshold = NA,
                     keep_by = "padj")
  expect_equal(sum(out$dedup_kept), 1)
})

test_that('NA padj is treated as non-significant (kept)', {
  enrich <- tibble::tibble(
    contrast = "c", database = "db", mode = "fgsea",
    pathway = c("P1", "P2"),
    pval = c(0.001, 0.001),
    padj = c(0.01, NA_real_),
    NES = c(2, 1.5), size = c(10, 10), direction = "up",
    leading_edge = c("A;B;C;D;E;F;G;H;I;J", "A;B;C;D;E;F;G;H;I;K")
  )
  out <- ev_collapse(enrich, method = "jaccard", cutoff = 0.5,
                     scope = "within_db", sig_threshold = 0.05,
                     keep_by = "padj")
  expect_true(all(out$dedup_kept))
})

test_that('collapse_then_jaccard runs collapse first, then Jaccard on survivors', {
  # Two glycolysis-like pathways that share leading edge (collapse drops one),
  # plus a third with high Jaccard to the survivor (Jaccard drops it).
  enrich <- tibble::tibble(
    contrast = "c", database = "db", mode = "fgsea",
    pathway = c("HALLMARK_GLYCOLYSIS", "REACTOME_GLYCOLYSIS",
                "WP_GLYCOLYSIS_REDUNDANT"),
    pval = c(0.001, 0.002, 0.003),
    padj = c(0.01, 0.02, 0.03),
    NES = c(2.5, 2.0, 1.8),
    size = c(40, 35, 30),
    direction = "up",
    leading_edge = c(
      "ENO1;PKM;GAPDH;HK2;ALDOA;PFKM;LDHA;TPI1;PGAM1;PGK1",
      "ENO1;PKM;GAPDH;HK2;ALDOA;PFKM;LDHA;TPI1;PGAM1;ENO2",
      "ENO1;PKM;GAPDH;HK2;ALDOA;PFKM;LDHA;TPI1;PGAM1;ENO3"
    )
  )
  pw_list <- list(
    HALLMARK_GLYCOLYSIS      = strsplit(enrich$leading_edge[1], ";")[[1]],
    REACTOME_GLYCOLYSIS      = strsplit(enrich$leading_edge[2], ";")[[1]],
    WP_GLYCOLYSIS_REDUNDANT  = strsplit(enrich$leading_edge[3], ";")[[1]]
  )
  ranks <- stats::setNames(
    seq(3, by = -0.01, length.out = length(unique(unlist(pw_list)))),
    unique(unlist(pw_list))
  )
  attr(enrich, "ev_pathways") <- pw_list
  attr(enrich, "ev_stats") <- list(c = ranks)

  out <- ev_collapse(enrich, method = "collapse_then_jaccard",
                     cutoff = 0.5, scope = "within_db",
                     sig_threshold = 0.05, keep_by = "padj")
  expect_true(out$dedup_kept[1])
  expect_equal(sum(out$dedup_kept), 1)
})
