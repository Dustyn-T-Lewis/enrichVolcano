test_that("YvO reproduction config runs with rank_by=t and balanced ring", {
  skip_on_cran()
  skip_if_not_installed("msigdbr")
  fixture <- testthat::test_path("fixtures", "yvo_tidy.rds")
  yvo <- readRDS(fixture)
  yvo <- yvo[!is.na(yvo$gene) & !is.na(yvo$contrast), ]
  rk <- if ("t" %in% colnames(yvo)) "t" else "signed_p"
  p <- suppressWarnings(enrich_volcano(
    yvo, contrast = "Aging", rank_by = rk,
    databases = "hallmark",
    ring = list(max_terms = 12, direction_balance = TRUE),
    disc_color = "#1B9E77", enrich_padj = 0.25
  ))
  expect_s3_class(p, "enrichVolcano")
  expect_true(is.list(attr(p, "ev_data")))
})

test_that("YvO end-to-end: enrich_volcano runs and returns patchwork", {
  skip_on_cran()
  fx <- file.path(test_path("fixtures"), "yvo_tidy.rds")
  skip_if_not(file.exists(fx), "YvO fixture not built")
  data <- readRDS(fx)
  ctrs <- intersect(unique(data$contrast), c("Aging", "Training_Young",
                                              "Training_Old", "Interaction"))
  skip_if(length(ctrs) == 0, "No expected YvO contrasts present")
  # Use mini pathway list to avoid msigdbr dependency in tests
  paths <- list(
    HALLMARK_TEST_UP = head(data$gene[data$logFC > 0], 20),
    HALLMARK_TEST_DN = head(data$gene[data$logFC < 0], 20)
  )
  p <- suppressWarnings(enrich_volcano(
    data, contrast = ctrs[1],
    databases = list(test = paths),
    p_method = "pi_eq2", p_adjust = "BH",
    enrich_mode = "fgsea", enrich_padj = 1.0
  ))
  expect_s3_class(p, "enrichVolcano")
  expect_true(!is.null(attr(p, "ev_call")))
  expect_named(attr(p, "ev_data"),
               c("validated_input", "pi_scores", "enrichment", "dedup_result"))
})

test_that('default pipeline reproduces YvO F03 panel A kept-pathway set', {
  fixture <- system.file("extdata", "examples", "yvo_tidy.rds",
                         package = "enrichVolcano")
  if (!nzchar(fixture)) fixture <- testthat::test_path("fixtures", "yvo_tidy.rds")
  skip_if_not(file.exists(fixture), "yvo_tidy fixture missing")
  res <- readRDS(fixture)
  res <- res[!is.na(res$gene) & !is.na(res$contrast), ]
  res <- pi_score(res, variant = "eq2")
  res <- adjust_p(res, method = "BH")

  # Defaults exercise collapse_then_jaccard, sig_threshold = 0.05,
  # filter_mode = "before", exclude_terms = DISEASE|CANCER|TUMOR.
  e <- ev_enrich(res, contrast = "Aging",
                 databases = c("hallmark", "go_bp"),
                 enrich_mode = "fgsea",
                 min_size = 10, max_size = 500,
                 seed = 42)
  d <- ev_collapse(e)

  kept <- d[d$dedup_kept & d$padj < 0.05, c("database", "pathway")]
  testthat::expect_snapshot_value(
    sort(paste(kept$database, kept$pathway, sep = "::")),
    style = "json2"
  )
})
