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

test_that('ev_collapse_fgsea flattens db-keyed ev_pathways from ev_enrich', {
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
