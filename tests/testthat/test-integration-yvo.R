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
