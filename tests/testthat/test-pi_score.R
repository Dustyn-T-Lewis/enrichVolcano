test_that("pi_score Eq.2 formula: pi = p^|logFC|", {
  d <- tibble::tibble(gene = "x", contrast = "c", logFC = 1, P.Value = 0.05)
  result <- pi_score(d, variant = "eq2")
  expect_equal(result$pi_eq2, 0.05)
})

test_that("pi_score Eq.1 formula: pi = |log2FC| * -log10(p)", {
  d <- tibble::tibble(gene = "x", contrast = "c", logFC = 1, P.Value = 0.1)
  result <- pi_score(d, variant = "eq1")
  expect_equal(result$pi_eq1, 1 * -log10(0.1))
})

test_that("pi_score is idempotent", {
  d <- tibble::tibble(gene = "x", contrast = "c", logFC = 1, P.Value = 0.05)
  once <- pi_score(d, variant = "eq2")
  twice <- pi_score(once, variant = "eq2")
  expect_identical(once, twice)
})

test_that("pi_score handles NA p-values", {
  d <- tibble::tibble(gene = "x", contrast = "c", logFC = 1, P.Value = NA_real_)
  result <- pi_score(d, variant = "eq2")
  expect_true(is.na(result$pi_eq2))
})

test_that("pi_score aborts on unknown variant via match.arg", {
  d <- tibble::tibble(gene = "x", contrast = "c", logFC = 1, P.Value = 0.05)
  expect_error(pi_score(d, variant = "eq3"))
})
