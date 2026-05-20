test_that("ev_validate passes through tidy long input unchanged", {
  data <- make_tidy_minimal()
  result <- ev_validate(data)
  expect_s3_class(result, "tbl_df")
  expect_true(all(c("gene", "contrast", "logFC", "P.Value", "adj.P.Val") %in% colnames(result)))
  expect_equal(attr(result, "ev_source"), "tidy_long")
})

test_that("ev_validate aborts on missing required columns", {
  bad <- data.frame(gene = "x", logFC = 1)
  expect_error(ev_validate(bad), class = "ev_input_missing_column")
})

test_that("ev_validate handles wide-suffix YvO format", {
  wide <- tibble::tibble(
    gene = paste0("G", 1:5),
    logFC_Cancer_vs_Healthy = rnorm(5),
    P.Value_Cancer_vs_Healthy = runif(5),
    adj.P.Val_Cancer_vs_Healthy = runif(5),
    logFC_Training_CR = rnorm(5),
    P.Value_Training_CR = runif(5),
    adj.P.Val_Training_CR = runif(5)
  )
  result <- ev_validate(wide)
  expect_equal(attr(result, "ev_source"), "wide_suffix")
  expect_setequal(unique(result$contrast), c("Cancer_vs_Healthy", "Training_CR"))
  expect_equal(nrow(result), 10)
})

test_that("ev_validate handles limma topTable single-contrast", {
  d <- tibble::tibble(
    gene = paste0("G", 1:5),
    logFC = rnorm(5),
    AveExpr = rnorm(5),
    t = rnorm(5),
    P.Value = runif(5),
    adj.P.Val = runif(5),
    B = rnorm(5),
    contrast = "A_vs_B"
  )
  result <- ev_validate(d)
  expect_equal(attr(result, "ev_source"), "limma")
})

test_that("ev_validate handles DESeq2 results", {
  d <- tibble::tibble(
    gene = paste0("G", 1:5),
    baseMean = runif(5),
    log2FoldChange = rnorm(5),
    lfcSE = runif(5),
    stat = rnorm(5),
    pvalue = runif(5),
    padj = runif(5),
    contrast = "ctr"
  )
  result <- ev_validate(d)
  expect_equal(attr(result, "ev_source"), "DESeq2")
  expect_true("logFC" %in% colnames(result))
  expect_true("P.Value" %in% colnames(result))
  expect_true("adj.P.Val" %in% colnames(result))
})

test_that("ev_validate handles edgeR topTags", {
  d <- tibble::tibble(
    gene = paste0("G", 1:5),
    logFC = rnorm(5),
    logCPM = runif(5),
    F = rnorm(5),
    PValue = runif(5),
    FDR = runif(5),
    contrast = "ctr"
  )
  result <- ev_validate(d)
  expect_equal(attr(result, "ev_source"), "edgeR")
})

test_that("ev_validate handles DEP wide with underscore-containing contrast", {
  dep_wide <- tibble::tibble(
    gene = paste0("P", 1:5),
    Cancer_vs_Healthy_T0_diff = rnorm(5),
    Cancer_vs_Healthy_T0_p.val = runif(5),
    Cancer_vs_Healthy_T0_p.adj = runif(5),
    Cancer_vs_Healthy_T0_significant = c(TRUE, FALSE, TRUE, FALSE, TRUE)
  )
  result <- ev_validate(dep_wide)
  expect_equal(attr(result, "ev_source"), "DEP_wide")
  expect_setequal(unique(result$contrast), "Cancer_vs_Healthy_T0")
  expect_true(all(c("logFC", "P.Value", "adj.P.Val") %in% colnames(result)))
})

test_that("ev_validate handles Perseus -log10 p-values", {
  perseus <- tibble::tibble(
    `Gene names` = paste0("G", 1:5),
    `Student's T-test Difference A_B` = rnorm(5),
    `-Log Student's T-test p-value A_B` = c(2, 1, 0.5, 3, 1.5),
    `Student's T-test q-value A_B` = runif(5)
  )
  expect_message(result <- ev_validate(perseus), regexp = "back-transforming")
  expect_equal(attr(result, "ev_source"), "Perseus_wide")
  expect_equal(result$P.Value[result$gene == "G1"], 0.01, tolerance = 1e-6)
})

test_that("ev_validate handles MaxQuant ratios with BH synthesis", {
  mq <- tibble::tibble(
    gene = paste0("G", 1:5),
    `Ratio H/L normalized expr1` = c(1, 2, 0.5, 4, 0.25),
    `Ratio H/L significance(B) expr1` = c(0.01, 0.5, 0.001, 0.3, 0.1)
  )
  expect_message(result <- ev_validate(mq), regexp = "BH")
  expect_equal(attr(result, "ev_source"), "MaxQuant_wide")
  expect_true("adj.P.Val" %in% colnames(result))
})

test_that("ev_validate derives P.Value from t + df when missing", {
  d <- tibble::tibble(
    gene = paste0("G", 1:5),
    logFC = rnorm(5),
    t = c(2, -2, 0, 3, -3),
    df = rep(10, 5),
    contrast = "ctr"
  )
  expect_message(result <- ev_validate(d), regexp = "Deriving|deriving")
  expect_true("P.Value" %in% colnames(result))
  expect_equal(result$P.Value[result$gene == "G1"], 2 * pt(-2, 10), tolerance = 1e-6)
})
