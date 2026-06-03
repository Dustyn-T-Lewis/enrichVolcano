test_that("collapse_fgsea removes redundant pathways sharing leading-edge genes", {
  skip_if_not_installed("fgsea")
  set.seed(99)
  # Strong, consistent signal in genes G1..G25 so fgsea calls the overlapping sets
  data <- tibble::tibble(
    gene     = paste0("G", 1:200),
    contrast = "ctr",
    logFC    = c(rnorm(25, 4, 0.3), rnorm(175, 0, 1)),
    P.Value  = c(runif(25, 0, 1e-4), runif(175, 0.2, 1))
  )
  data <- pi_score(data, variant = "eq2")
  # Two near-identical gene sets (should collapse) + two distinct sets
  paths <- list(
    SET_A     = paste0("G", 1:25),
    SET_A_DUP = paste0("G", 1:24),   # ~96% overlap with SET_A
    SET_B     = paste0("G", 60:90),
    SET_C     = paste0("G", 120:160)
  )
  enrich <- ev_enrich(
    data, contrast = "ctr",
    databases    = list(test = paths),
    enrich_mode  = "fgsea", rank_by = "signed_p",
    min_size = 1, max_size = 1000, nperm = 1000
  )
  expect_false(is.null(attr(enrich, "ev_pathways")))
  expect_false(is.null(attr(enrich, "ev_stats")))
  collapsed <- ev_collapse(enrich, method = "collapse_fgsea")
  expect_true("dedup_kept" %in% colnames(collapsed))
  # Mechanism check: collapse must keep <= all, and the kept set is a valid subset
  expect_lte(sum(collapsed$dedup_kept), nrow(enrich))
  expect_true(all(collapsed$pathway[collapsed$dedup_kept] %in% enrich$pathway))
})

test_that("collapse method='jaccard_then_collapse' applies jaccard then fgsea", {
  skip_if_not_installed("fgsea")
  skip_if_not_installed("msigdbr")

  data <- make_ranked_input()
  data <- pi_score(data, variant = "eq2")
  m <- msigdbr::msigdbr(species = "Homo sapiens", collection = "H")
  paths <- split(m$gene_symbol, m$gs_name)[1:5]

  enrich <- ev_enrich(
    data, contrast = "ctr",
    databases = list(hallmark_sub = paths),
    enrich_mode = "fgsea", rank_by = "signed_p",
    min_size = 1, max_size = 1000, nperm = 1000
  )

  both <- ev_collapse(enrich, method = "jaccard_then_collapse", cutoff = 0.3)
  jac  <- ev_collapse(enrich, method = "jaccard", cutoff = 0.3)
  expect_lte(sum(both$dedup_kept), sum(jac$dedup_kept),
             label = "'jaccard_then_collapse' should keep <= jaccard alone")
})

test_that("collapse_fgsea warns and falls back when ev_pathways attribute missing", {
  skip_if_not_installed("fgsea")
  enrich <- tibble::tibble(
    contrast = "ctr", database = "test",
    pathway = paste0("P", 1:5),
    pval = runif(5), padj = runif(5), NES = rnorm(5), size = 10,
    leading_edge = paste0("G", 1:5), mode = "fgsea",
    direction = "up"
  )
  # sig_threshold = NA bypasses the sig gate so collapse_fgsea is reached
  # and can emit its missing-attribute warning regardless of fixture padj.
  expect_warning(
    out <- ev_collapse(enrich, method = "collapse_fgsea", sig_threshold = NA),
    class = "ev_collapse_no_pathways"
  )
  expect_true(all(out$dedup_kept))
})

test_that("ev_collapse_fgsea flattens db-keyed ev_pathways from ev_enrich", {
  # Mimic the structure ev_enrich attaches: list keyed by db name, each
  # value a named list of pathway -> gene character vectors. Build two
  # genuinely independent pathways: HALLMARK_GLYCOLYSIS enriches at the
  # top of the ranked list (G1-G40, ranks +3 -> +1) and HALLMARK_HYPOXIA
  # at the bottom (G61-G100, ranks -1 -> -3), separated by buffer genes
  # G41-G60 with ~0 ranks. With disjoint leading edges and opposite
  # directions, collapsePathways should retain both.
  set.seed(11)
  enrich <- tibble::tibble(
    contrast = "c", database = c("hallmark", "hallmark"),
    mode = "fgsea",
    pathway = c("HALLMARK_GLYCOLYSIS", "HALLMARK_HYPOXIA"),
    pval = c(0.001, 0.002), padj = c(0.01, 0.02),
    NES = c(2.5, -2.0), size = c(40, 40), direction = c("up", "down"),
    leading_edge = c(
      paste0("G", 1:10, collapse = ";"),
      paste0("G", 91:100, collapse = ";")
    )
  )
  db_keyed <- list(
    hallmark = list(
      HALLMARK_GLYCOLYSIS = paste0("G", 1:40),
      HALLMARK_HYPOXIA    = paste0("G", 61:100)
    )
  )
  ranks <- c(
    stats::setNames(seq(3, 1, length.out = 40), paste0("G", 1:40)),
    stats::setNames(rnorm(20, 0, 0.1), paste0("G", 41:60)),
    stats::setNames(seq(-1, -3, length.out = 40), paste0("G", 61:100))
  )
  attr(enrich, "ev_pathways") <- db_keyed
  attr(enrich, "ev_stats") <- list(c = ranks)

  out <- ev_collapse(enrich, method = "collapse_fgsea", cutoff = 0.5,
                     scope = "within_db", sig_threshold = 0.05,
                     keep_by = "padj")
  # Both pathways have disjoint leading edges → collapsePathways should
  # keep both. Regression check: with the unfixed indexing, this would
  # drop everything because pathway_list[fg_subset$pathway] returns NULLs.
  expect_true(all(out$dedup_kept))
})

test_that("ev_collapse_fgsea: missing stats for contrast warns and keeps all", {
  # Triggers lines 168-172: ev_stats present but with no entry for this contrast.
  enrich <- tibble::tibble(
    contrast = "missing_ctr", database = "db", mode = "fgsea",
    pathway = c("P1", "P2"),
    pval = c(0.001, 0.002), padj = c(0.01, 0.02),
    NES = c(2, 1.5), size = c(10, 10), direction = "up",
    leading_edge = c("A;B;C;D;E", "F;G;H;I;J")
  )
  attr(enrich, "ev_pathways") <- list(P1 = LETTERS[1:5], P2 = LETTERS[6:10])
  # Stats keyed by a different contrast — lookup returns NULL.
  attr(enrich, "ev_stats") <- list(other_ctr = stats::setNames(rnorm(10), LETTERS[1:10]))
  expect_warning(
    out <- ev_collapse(enrich, method = "collapse_fgsea", sig_threshold = NA),
    class = "ev_collapse_no_pathways"
  )
  expect_true(all(out$dedup_kept))
})

test_that("ev_collapse_fgsea: returns TRUE for groups with no fgsea rows", {
  # Triggers line 164: nrow(fg_subset) == 0 early-return when mode != fgsea.
  enrich <- tibble::tibble(
    contrast = "c", database = "db", mode = "ora",
    pathway = c("P1", "P2"),
    pval = c(0.001, 0.002), padj = c(0.01, 0.02),
    NES = c(NA_real_, NA_real_), size = c(10, 10), direction = "up",
    leading_edge = c("A;B;C", "D;E;F")
  )
  attr(enrich, "ev_pathways") <- list(P1 = LETTERS[1:5], P2 = LETTERS[6:10])
  attr(enrich, "ev_stats") <- list(c = stats::setNames(rnorm(10), LETTERS[1:10]))
  out <- ev_collapse(enrich, method = "collapse_fgsea", sig_threshold = NA)
  expect_true(all(out$dedup_kept))
})

test_that("jaccard_then_collapse never produces empty redundancy clusters", {
  # Regression for the AND-merge defect: the old code passed the full sig_idx
  # (not just Jaccard survivors) to ev_collapse_fgsea and ANDed the result
  # against the existing dedup_kept flags. When collapsePathways and Jaccard's
  # keep_by tie-breaker disagreed on which member represents a redundant
  # cluster, the AND zeroed out both, producing clusters with no kept rows.
  # After the fix, collapse runs only on Jaccard survivors and writes back
  # positionally, so every redundancy cluster keeps >= 1 representative.
  skip_if_not_installed("fgsea")
  enrich <- tibble::tibble(
    contrast = "c", database = "db", mode = "fgsea",
    pathway = c("PW_A_small_pval", "PW_B_large_NES",
                "PW_C_small_pval", "PW_D_large_NES"),
    # Pval ordering disagrees with NES ordering on which row represents each
    # redundant cluster (A-B share a leading edge; C-D share another).
    pval = c(0.0001, 0.002, 0.0002, 0.003),
    padj = c(0.001, 0.02, 0.002, 0.03),
    NES = c(2.0, 3.0, 2.1, 3.2),
    size = c(20, 20, 20, 20),
    direction = "up",
    leading_edge = c(
      # A and B share 9/10 leading-edge genes (Jaccard ~ 0.82 > 0.5).
      paste(paste0("G", 1:10), collapse = ";"),
      paste(c(paste0("G", 1:9), "G11"), collapse = ";"),
      # C and D share 9/10 leading-edge genes (Jaccard ~ 0.82 > 0.5).
      paste(paste0("G", 50:59), collapse = ";"),
      paste(c(paste0("G", 50:58), "G60"), collapse = ";")
    )
  )
  pw_list <- list(
    db = list(
      PW_A_small_pval = paste0("G", 1:20),
      PW_B_large_NES  = paste0("G", 1:20),
      PW_C_small_pval = paste0("G", 50:70),
      PW_D_large_NES  = paste0("G", 50:70)
    )
  )
  ranks <- c(
    stats::setNames(seq(3, 1, length.out = 20), paste0("G", 1:20)),
    stats::setNames(rnorm(29, 0, 0.1), paste0("G", 21:49)),
    stats::setNames(seq(3, 1, length.out = 21), paste0("G", 50:70))
  )
  attr(enrich, "ev_pathways") <- pw_list
  attr(enrich, "ev_stats") <- list(c = ranks)

  out <- suppressWarnings(
    ev_collapse(enrich, method = "jaccard_then_collapse",
                cutoff = 0.5, scope = "within_db",
                sig_threshold = 0.05, keep_by = "NES")
  )
  # Both redundancy clusters must keep at least one representative — the
  # previous AND-merge could leave both clusters empty when keep_by ranking
  # and pval ranking disagreed.
  expect_gte(sum(out$dedup_kept), 2)
})

test_that("ev_collapse_fgsea resolves pathway-name collisions across databases", {
  # Regression for the flatten defect: two databases each define a pathway
  # called "GLYCOLYSIS" but with disjoint gene members. The old flatten
  # (do.call(c, unname(pathway_list))) preserved duplicate names, so name-based
  # lookup returned db1's genes for both rows. db-aware lookup must return the
  # correct member of each pair.
  skip_if_not_installed("fgsea")
  set.seed(7)
  enrich <- tibble::tibble(
    contrast = "c", database = c("db1", "db2"), mode = "fgsea",
    pathway = c("GLYCOLYSIS", "GLYCOLYSIS"),
    pval = c(0.001, 0.002), padj = c(0.01, 0.02),
    NES = c(2.5, -2.0), size = c(20, 20),
    direction = c("up", "down"),
    leading_edge = c(
      paste0("G", 1:10, collapse = ";"),
      paste0("G", 91:100, collapse = ";")
    )
  )
  # Disjoint gene sets — collapsePathways must NOT collapse them when each
  # row is paired with its real gene set, regardless of name collision.
  db_keyed <- list(
    db1 = list(GLYCOLYSIS = paste0("G", 1:20)),
    db2 = list(GLYCOLYSIS = paste0("G", 81:100))
  )
  ranks <- c(
    stats::setNames(seq(3, 1, length.out = 40), paste0("G", 1:40)),
    stats::setNames(rnorm(20, 0, 0.1), paste0("G", 41:60)),
    stats::setNames(seq(-1, -3, length.out = 40), paste0("G", 61:100))
  )
  attr(enrich, "ev_pathways") <- db_keyed
  attr(enrich, "ev_stats") <- list(c = ranks)

  out <- ev_collapse(enrich, method = "collapse_fgsea", cutoff = 0.5,
                     scope = "within_db", sig_threshold = 0.05,
                     keep_by = "padj")
  # Both pathways should survive because their real gene sets are disjoint.
  # Under the pre-fix flatten the second row's row was scored against db1's
  # gene set and the test would have dropped one or both rows.
  expect_true(all(out$dedup_kept))
})

test_that("ev_collapse_fgsea scores each contrast against its own stats vector", {
  # Regression for scope='global' + multi-contrast: ev_collapse_fgsea must
  # partition fg_subset by contrast and re-run collapsePathways per contrast
  # with the matching stats. Previously it pulled stats for contrast[1] only,
  # so contrast 2's pathway rows were scored against contrast 1's ranks.
  skip_if_not_installed("fgsea")
  enrich <- tibble::tibble(
    contrast = c("ctr_up", "ctr_dn"),
    database = c("db", "db"),
    mode = "fgsea",
    pathway = c("PW_X", "PW_Y"),
    pval = c(0.001, 0.001), padj = c(0.01, 0.01),
    NES = c(2.5, -2.5), size = c(20, 20),
    direction = c("up", "down"),
    leading_edge = c(
      paste0("G", 1:10, collapse = ";"),
      paste0("G", 81:90, collapse = ";")
    )
  )
  attr(enrich, "ev_pathways") <- list(
    db = list(PW_X = paste0("G", 1:20),
              PW_Y = paste0("G", 81:100))
  )
  # ctr_up ranks PW_X at the top, ctr_dn ranks PW_Y at the top. If the loop
  # erroneously scored both contrasts against ctr_up's stats, PW_Y would land
  # at the bottom (negative ES) and collapsePathways behavior would differ.
  ranks_up <- c(
    stats::setNames(seq(3, 1, length.out = 20), paste0("G", 1:20)),
    stats::setNames(rnorm(60, 0, 0.1), paste0("G", 21:80)),
    stats::setNames(seq(-1, -3, length.out = 20), paste0("G", 81:100))
  )
  ranks_dn <- c(
    stats::setNames(seq(3, 1, length.out = 20), paste0("G", 81:100)),
    stats::setNames(rnorm(60, 0, 0.1), paste0("G", 21:80)),
    stats::setNames(seq(-1, -3, length.out = 20), paste0("G", 1:20))
  )
  attr(enrich, "ev_stats") <- list(ctr_up = ranks_up, ctr_dn = ranks_dn)

  out <- ev_collapse(enrich, method = "collapse_fgsea", cutoff = 0.5,
                     scope = "global", sig_threshold = 0.05,
                     keep_by = "padj")
  # Each contrast has a single pathway in its partition; both must survive.
  expect_true(all(out$dedup_kept))
})

test_that("ev_collapse_fgsea: catches errors from fgsea::collapsePathways", {
  # Triggers lines 196-203 (tryCatch error path). Force an error by handing
  # collapsePathways a stats vector containing NAs, which it rejects.
  skip_if_not_installed("fgsea")
  enrich <- tibble::tibble(
    contrast = "c", database = "db", mode = "fgsea",
    pathway = c("P1", "P2"),
    pval = c(0.001, 0.002), padj = c(0.01, 0.02),
    NES = c(2, 1.5), size = c(5, 5), direction = "up",
    leading_edge = c("A;B;C;D;E", "F;G;H;I;J")
  )
  attr(enrich, "ev_pathways") <- list(
    P1 = c("A", "B", "C", "D", "E"),
    P2 = c("F", "G", "H", "I", "J")
  )
  # Stats with NA values causes fgsea to throw.
  attr(enrich, "ev_stats") <- list(c = stats::setNames(
    c(NA_real_, NA_real_, NA_real_, NA_real_, NA_real_,
      NA_real_, NA_real_, NA_real_, NA_real_, NA_real_),
    LETTERS[1:10]
  ))
  expect_warning(
    out <- ev_collapse(enrich, method = "collapse_fgsea",
                       sig_threshold = NA),
    class = "ev_collapse_fgsea_error"
  )
  expect_true(all(out$dedup_kept))
})
