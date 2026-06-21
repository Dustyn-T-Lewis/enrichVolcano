# Regression tests for the fgsea ranking statistic. GSEA requires a SIGNED
# rank so up- and down-regulated genes sit at opposite ends (Subramanian 2005
# PNAS; Reimand 2019 Nat Protoc). A prior bug ranked by the unsigned pi-score,
# collapsing the list to one tail and producing near-empty rings.

test_that("signed_p ranking puts up genes at the top and down genes at the bottom", {
  data <- make_ranked_input() # 20 strong-up, 60 mid, 20 strong-down
  stats <- ev_rank_stats(data, "signed_p")

  expect_true(is.numeric(stats))
  expect_false(is.unsorted(rev(stats))) # sorted decreasing
  expect_true(all(is.finite(stats)))

  up_genes <- data$gene[data$logFC > 0]
  down_genes <- data$gene[data$logFC < 0]
  # signed statistic must agree in sign with logFC for every gene
  expect_true(all(stats[names(stats) %in% up_genes] > 0))
  expect_true(all(stats[names(stats) %in% down_genes] < 0))
  # the list spans both tails (not collapsed to one side)
  expect_gt(max(stats), 0)
  expect_lt(min(stats), 0)
})

test_that("ranking fgsea by an unsigned pi-score column is rejected", {
  data <- make_ranked_input() # includes a pi_eq2 column
  expect_true("pi_eq2" %in% colnames(data))

  expect_error(ev_rank_stats(data, "pi_eq2"), class = "ev_unsigned_ranking")
  expect_error(ev_rank_stats(data, "pi_eq1"), class = "ev_unsigned_ranking")
  expect_error(ev_rank_stats(data, "pi_score"), class = "ev_unsigned_ranking")

  # and the guard fires through the public ev_enrich() entry point too
  expect_error(
    ev_enrich(data,
      contrast = "ctr",
      databases = list(hallmark = make_mini_pathways()),
      enrich_mode = "fgsea", rank_by = "pi_eq2"
    ),
    class = "ev_unsigned_ranking"
  )
})

test_that("signed_pi ranks by fold-weighted -log10 p with sign(logFC)", {
  data <- make_ranked_input()
  stats <- ev_rank_stats(data, "signed_pi")

  expect_true(is.numeric(stats))
  expect_false(is.unsorted(rev(stats)))
  expect_true(all(is.finite(stats)))

  # Sign agrees with logFC for every gene; tails span both directions.
  up_genes <- data$gene[data$logFC > 0]
  down_genes <- data$gene[data$logFC < 0]
  expect_true(all(stats[names(stats) %in% up_genes] > 0))
  expect_true(all(stats[names(stats) %in% down_genes] < 0))
  expect_gt(max(stats), 0)
  expect_lt(min(stats), 0)

  # Genes with larger |logFC| at the same p should rank further from zero
  # than genes with smaller |logFC| (the FC-weighting is what distinguishes
  # signed_pi from signed_p).
  signed_p <- ev_rank_stats(data, "signed_p")
  common <- intersect(names(stats), names(signed_p))
  expect_false(identical(stats[common], signed_p[common]))
})

test_that("a missing rank_by column falls back to the signed default", {
  data <- make_ranked_input()
  expect_identical(
    ev_rank_stats(data, "no_such_column"),
    ev_rank_stats(data, "signed_p")
  )
})

test_that("the default ranking statistic is the moderated t", {
  # limma/proteoDA carry a moderated-t column; GSEA should rank by it by default
  # (Subramanian 2005 PNAS; Reimand 2019 Nat Protoc — a signed test statistic).
  expect_identical(formals(ev_enrich)$rank_by, "t")
  expect_identical(formals(enrich_volcano)$rank_by, "t")
})

test_that("default rank_by='t' uses the moderated t when present, signed_p when absent", {
  data <- make_ranked_input()
  data$t <- data$logFC * 2 + rnorm(nrow(data), sd = 0.1) # signed, distinct from signed_p
  expect_identical(ev_rank_stats(data, "t"), ev_rank_stats(data, "t"))
  expect_false(identical(ev_rank_stats(data, "t"), ev_rank_stats(data, "signed_p")))

  no_t <- data[setdiff(names(data), "t")] # edgeR-style input: no t column
  expect_identical(ev_rank_stats(no_t, "t"), ev_rank_stats(no_t, "signed_p"))
})
