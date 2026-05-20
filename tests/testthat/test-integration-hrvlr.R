test_that("HRvLR end-to-end runs on real subsample", {
  skip_on_cran()
  fx <- file.path(test_path("fixtures"), "hrvlr_tidy.rds")
  skip_if_not(file.exists(fx), "HRvLR fixture not built")
  data <- readRDS(fx)
  skip_if(nrow(data) < 50, "HRvLR fixture too small")
  ctrs <- unique(data$contrast)
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
})
