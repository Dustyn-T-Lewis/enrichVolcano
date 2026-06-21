test_that("apply_gate on adj.P.Val applies both p and logFC thresholds", {
  d <- tibble::tibble(
    gene = c("A", "B", "C", "D"),
    logFC = c(2, 0.1, 2, -2),
    adj.P.Val = c(0.01, 0.001, 0.5, 0.001)
  )
  out <- apply_gate(d, value = "adj.P.Val", p_cut = 0.05, fc_cut = log2(1.5))
  # A: p<0.05 and |logFC|>fc_cut -> up
  # B: p<0.05 but |logFC| too small -> ns
  # C: |logFC| ok but p too large -> ns
  # D: p<0.05 and |logFC|>fc_cut, negative -> down
  expect_identical(out$sig, c(TRUE, FALSE, FALSE, TRUE))
  expect_identical(out$direction, c("up", "ns", "ns", "down"))
})

test_that("apply_gate on pi_eq2 thresholds Pi alone (no second logFC gate)", {
  # Per Xiao 2014, Pi already folds in fold change; the convention is to
  # threshold Pi directly. A small-|logFC| gene whose Pi crosses the cutoff
  # must still be flagged sig.
  d <- tibble::tibble(
    gene = c("A", "B"),
    logFC = c(0.1, 2),         # one tiny, one large
    pi_eq2 = c(0.01, 0.5)      # only A crosses 0.05
  )
  out <- apply_gate(d, value = "pi_eq2", p_cut = 0.05)
  expect_identical(out$sig, c(TRUE, FALSE))
  expect_identical(out$direction, c("up", "ns"))
})

test_that("apply_gate treats NA stat as not significant", {
  d <- tibble::tibble(
    gene = c("A", "B"),
    logFC = c(2, 2),
    adj.P.Val = c(NA_real_, 0.01)
  )
  out <- apply_gate(d, value = "adj.P.Val")
  expect_identical(out$sig, c(FALSE, TRUE))
  expect_identical(out$direction, c("ns", "up"))
})

test_that("apply_gate aborts when the chosen column is missing", {
  d <- tibble::tibble(gene = "A", logFC = 1, P.Value = 0.01)
  expect_error(apply_gate(d, value = "pi_eq2"),
               class = "ev_missing_gate_column")
})

test_that("apply_gate aborts when logFC is missing", {
  d <- tibble::tibble(gene = "A", adj.P.Val = 0.01)
  expect_error(apply_gate(d, value = "adj.P.Val"),
               class = "ev_missing_logfc")
})
