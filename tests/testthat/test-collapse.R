# Set-similarity metrics (Merico 2010 PMID:21085593, EnrichmentMap): a small set
# fully contained in a large one has low Jaccard but Overlap = 1. The combined
# coefficient (0.5*Jaccard + 0.5*Overlap) is the Cytoscape EnrichmentMap default.
make_containment_pair <- function() {
  tibble::tibble(
    contrast = "c", database = "db", mode = "fgsea",
    pathway = c("BIG", "SUBSET"),
    pval = c(0.001, 0.01), padj = c(0.01, 0.02),
    NES = c(2, 1.8), size = c(10, 3), direction = c("up", "up"),
    # SUBSET (A;B;C) ⊂ BIG: Jaccard = 3/10 = 0.30; Overlap = 3/3 = 1.00
    leading_edge = c("A;B;C;D;E;F;G;H;I;J", "A;B;C")
  )
}

test_that("similarity defaults to jaccard and is exposed as a parameter", {
  expect_identical(eval(formals(ev_collapse)$similarity), c("jaccard", "overlap", "combined"))
  enr <- make_containment_pair()
  expect_identical(
    ev_collapse(enr, method = "jaccard", cutoff = 0.5, sig_threshold = NA),
    ev_collapse(enr, method = "jaccard", similarity = "jaccard",
                cutoff = 0.5, sig_threshold = NA)
  )
})

test_that("overlap coefficient drops a contained set that Jaccard keeps", {
  enr <- make_containment_pair()
  jac <- ev_collapse(enr, method = "jaccard", similarity = "jaccard",
                     cutoff = 0.5, sig_threshold = NA)
  ovl <- ev_collapse(enr, method = "jaccard", similarity = "overlap",
                     cutoff = 0.5, sig_threshold = NA)
  expect_equal(sum(jac$dedup_kept), 2)          # Jaccard 0.30 < 0.5 -> keep both
  expect_equal(sum(ovl$dedup_kept), 1)          # Overlap 1.00 > 0.5 -> drop SUBSET
  expect_true(ovl$dedup_kept[ovl$pathway == "BIG"])
})

test_that("combined coefficient = weighted Jaccard + overlap (Merico 2010)", {
  enr <- make_containment_pair()
  # combined = 0.5*0.30 + 0.5*1.00 = 0.65 > 0.5 -> drop SUBSET
  cmb <- ev_collapse(enr, method = "jaccard", similarity = "combined",
                     cutoff = 0.5, sig_threshold = NA)
  expect_equal(sum(cmb$dedup_kept), 1)
  # weight entirely on Jaccard reproduces the Jaccard-only verdict (0.30 < 0.5)
  cmb_j <- ev_collapse(enr, method = "jaccard", similarity = "combined",
                       combined_weight = 1, cutoff = 0.5, sig_threshold = NA)
  expect_equal(sum(cmb_j$dedup_kept), 2)
})

test_that("staged scope runs within_db then cross_db on survivors", {
  # P1/P2 are redundant WITHIN db1; P1/P3 are redundant ACROSS db1/db2.
  enr <- tibble::tibble(
    contrast = "c", database = c("db1", "db1", "db2"), mode = "fgsea",
    pathway = c("P1", "P2", "P3"),
    pval = c(0.001, 0.01, 0.02), padj = c(0.005, 0.01, 0.02),
    NES = c(2, 1.8, 1.6), size = c(10, 10, 10), direction = "up",
    leading_edge = c("A;B;C;D;E", "A;B;C;D;F", "A;B;C;D;G")
  )
  staged <- ev_collapse(enr, method = "jaccard", scope = c("within_db", "cross_db"),
                        cutoff = 0.5, sig_threshold = NA)
  # within_db drops P2 (vs P1 in db1); cross_db then drops P3 (vs P1) -> only P1.
  expect_equal(sum(staged$dedup_kept), 1)
  expect_true(staged$dedup_kept[staged$pathway == "P1"])

  # within_db alone keeps P1 and the cross-db P3 (different DBs, never compared).
  within_only <- ev_collapse(enr, method = "jaccard", scope = "within_db",
                             cutoff = 0.5, sig_threshold = NA)
  expect_equal(sum(within_only$dedup_kept), 2)
})

test_that("per-stage cutoff vector is accepted and validated", {
  enr <- make_containment_pair()
  expect_silent(ev_collapse(enr, method = "jaccard", scope = c("within_db", "cross_db"),
                            cutoff = c(0.9, 0.375), sig_threshold = NA))
  expect_error(
    ev_collapse(enr, method = "jaccard", scope = c("within_db", "cross_db"),
                cutoff = c(0.5, 0.5, 0.5), sig_threshold = NA),
    class = "ev_bad_cutoff")
})

test_that("keep_by='NES' keeps the strongest term by magnitude, not signed NES", {
  # Two redundant DOWN pathways; the strongly-down one (NES -3) must be the
  # representative, not the weakly-down one (NES -0.5).
  enr <- tibble::tibble(
    contrast = "c", database = "db", mode = "fgsea",
    pathway = c("STRONG_DOWN", "WEAK_DOWN"),
    pval = c(0.001, 0.04), padj = c(0.01, 0.04),
    NES = c(-3, -0.5), size = c(10, 10), direction = c("down", "down"),
    leading_edge = c("A;B;C;D;E;F;G;H;I;J", "A;B;C;D;E;F;G;H;I;K")
  )
  res <- ev_collapse(enr, method = "jaccard", keep_by = "NES",
                     cutoff = 0.5, sig_threshold = NA)
  expect_equal(sum(res$dedup_kept), 1)
  expect_true(res$dedup_kept[res$pathway == "STRONG_DOWN"])
})

test_that("dedup never collapses an up- and a down-regulated pathway together", {
  # Identical leading edge but opposite direction = distinct biology; keep both.
  enr <- tibble::tibble(
    contrast = "c", database = "db", mode = "fgsea",
    pathway = c("UP", "DOWN"),
    pval = c(0.001, 0.01), padj = c(0.01, 0.02),
    NES = c(2, -2), size = c(10, 10), direction = c("up", "down"),
    leading_edge = c("A;B;C;D;E", "A;B;C;D;E")
  )
  res <- ev_collapse(enr, method = "jaccard", cutoff = 0.5, sig_threshold = NA)
  expect_equal(sum(res$dedup_kept), 2)
})

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

test_that("sig_threshold gates dedup to significant pathways only", {
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

test_that("sig_threshold = NA disables gating (dedup all rows)", {
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

test_that("NA padj is treated as non-significant (kept)", {
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

test_that("collapse_then_jaccard runs collapse first, then Jaccard on survivors", {
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

test_that("ev_collapse default method is collapse_then_jaccard", {
  fn <- ev_collapse
  default <- eval(formals(fn)$method)[1]
  expect_identical(default, "collapse_then_jaccard")
})

test_that("collapse_then_jaccard: Jaccard prunes survivors when >=2 survive", {
  # Exercises R/collapse.R lines 110-116: collapse keeps >=2 pathways, then
  # Jaccard further dedups them.
  #
  # Trick: ev_collapse_fgsea uses the supplied pathway *gene sets* to re-run
  # fgsea internally — so we make those structurally independent (3 separated
  # signal blocks in 500-gene background). collapsePathways then keeps all 3.
  # Separately, the *leading_edge* strings (which feed Jaccard) are crafted so
  # PW_A and PW_B share 19/20 genes (J ~= 0.90 > cutoff 0.5) and PW_C is
  # disjoint. After collapse-then-Jaccard, we expect {PW_A, PW_C} kept.
  set.seed(42)
  big_n <- 500
  ranks <- stats::setNames(stats::rnorm(big_n, 0, 1), paste0("G", seq_len(big_n)))
  ranks[paste0("G", 1:20)]    <- seq(5, 4, length.out = 20)
  ranks[paste0("G", 150:169)] <- seq(5, 4, length.out = 20)
  ranks[paste0("G", 350:369)] <- seq(5, 4, length.out = 20)

  pw_list <- list(
    PW_A = paste0("G", 1:20),
    PW_B = paste0("G", 150:169),
    PW_C = paste0("G", 350:369)
  )
  enrich <- tibble::tibble(
    contrast = "c", database = "db", mode = "fgsea",
    pathway = c("PW_A", "PW_B", "PW_C"),
    pval = c(0.001, 0.002, 0.003),
    padj = c(0.01, 0.02, 0.03),
    NES = c(2.5, 2.4, 2.3),
    size = c(20, 20, 20),
    direction = "up",
    leading_edge = c(
      # A and B leading edges share 19/20 (Jaccard ~ 0.90)
      paste(paste0("G", 1:20), collapse = ";"),
      paste(c(paste0("G", 1:19), "G170"), collapse = ";"),
      # C disjoint
      paste(paste0("G", 350:369), collapse = ";")
    )
  )
  attr(enrich, "ev_pathways") <- pw_list
  attr(enrich, "ev_stats") <- list(c = ranks)

  out <- suppressWarnings(
    ev_collapse(enrich, method = "collapse_then_jaccard",
                cutoff = 0.5, scope = "within_db",
                sig_threshold = 0.05, keep_by = "padj")
  )
  expect_true(out$dedup_kept[out$pathway == "PW_A"])
  expect_false(out$dedup_kept[out$pathway == "PW_B"])
  expect_true(out$dedup_kept[out$pathway == "PW_C"])
})
