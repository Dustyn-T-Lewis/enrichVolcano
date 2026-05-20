test_that("adjust_p BH matches stats::p.adjust", {
  d <- tibble::tibble(
    gene = paste0("G", 1:10),
    contrast = "c",
    logFC = rnorm(10),
    P.Value = c(0.001, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.1, 0.5, 0.9)
  )
  result <- adjust_p(d, method = "BH", by_contrast = FALSE)
  expect_equal(result$adj.P.Val, stats::p.adjust(d$P.Value, "BH"))
})

test_that("adjust_p Bonferroni matches stats::p.adjust", {
  d <- tibble::tibble(gene = c("x", "y"), contrast = "c", logFC = c(1, 1),
                       P.Value = c(0.01, 0.5))
  result <- adjust_p(d, method = "bonferroni", by_contrast = FALSE)
  expect_equal(result$adj.P.Val, c(0.02, 1))
})

test_that("adjust_p by_contrast adjusts within each group", {
  d <- tibble::tibble(
    gene = rep(paste0("G", 1:5), 2),
    contrast = rep(c("c1", "c2"), each = 5),
    logFC = rnorm(10),
    P.Value = c(rep(0.01, 5), rep(0.1, 5))
  )
  result <- adjust_p(d, method = "BH", by_contrast = TRUE)
  bh_c1 <- stats::p.adjust(rep(0.01, 5), "BH")
  bh_c2 <- stats::p.adjust(rep(0.1, 5), "BH")
  expect_equal(sort(result$adj.P.Val[result$contrast == "c1"]), sort(bh_c1))
  expect_equal(sort(result$adj.P.Val[result$contrast == "c2"]), sort(bh_c2))
})

test_that("adjust_p IHW requires covariate", {
  d <- tibble::tibble(gene = "x", contrast = "c", logFC = 1, P.Value = 0.01)
  expect_error(adjust_p(d, method = "IHW"), class = "ev_ihw_missing_covariate")
})
