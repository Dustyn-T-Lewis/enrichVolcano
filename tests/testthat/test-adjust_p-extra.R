test_that("adjust_p qvalue path runs when qvalue is installed", {
  skip_if_not_installed("qvalue")
  # Use a large, well-behaved p-value distribution so qvalue's pi0 smoother
  # (smooth.spline) has enough spread to estimate without erroring.
  set.seed(7)
  d <- tibble::tibble(
    gene = paste0("g", seq_len(1000)),
    contrast = "ctr",
    P.Value = c(runif(150, 0, 0.01), runif(850, 0, 1))
  )
  result <- adjust_p(d, method = "qvalue", by_contrast = FALSE)
  expect_true("adj.P.Val" %in% colnames(result))
  expect_true(all(result$adj.P.Val >= 0 & result$adj.P.Val <= 1))
})

test_that("adjust_p IHW path runs when IHW is installed and covariate is provided", {
  skip_if_not_installed("IHW")
  set.seed(4)
  d <- tibble::tibble(gene = paste0("g", 1:200),
                       contrast = "ctr",
                       P.Value = runif(200),
                       AveExpr = rnorm(200))
  result <- adjust_p(d, method = "IHW", covariate = "AveExpr",
                     by_contrast = FALSE)
  expect_true("adj.P.Val" %in% colnames(result))
})

test_that("adjust_p errors clearly when qvalue Suggests package is missing", {
  skip_if(requireNamespace("qvalue", quietly = TRUE), "qvalue installed; skipping miss-handler test")
  d <- tibble::tibble(gene = "g1", contrast = "ctr", P.Value = 0.5)
  expect_error(adjust_p(d, method = "qvalue"), class = "ev_missing_qvalue")
})
