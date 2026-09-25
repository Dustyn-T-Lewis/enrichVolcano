limma_table <- function() {
  data.frame(
    logFC = c(1.2, -0.8, 0.1), AveExpr = c(20, 18, 22), t = c(4.1, -3.2, 0.3),
    P.Value = c(1e-4, 2e-3, 0.7), adj.P.Val = c(1e-3, 0.01, 0.8), B = c(1, 0, -5),
    row.names = c("P31040", "Q9UBK2", "P12345-2")
  )
}

test_that("a limma topTable reads IDs from row names", {
  d <- as_da(limma_table(), contrast = "Aging")
  expect_s3_class(d, "enrichVolcano_da")
  expect_identical(d$protein, c("P31040", "Q9UBK2", "P12345-2"))
  expect_identical(unique(d$contrast), "Aging")
  expect_equal(d$t, c(4.1, -3.2, 0.3))
  expect_equal(d$abundance, c(20, 18, 22))
  expect_identical(unique(d$rank_stat), "t")
  expect_equal(d$rank, d$t)
})

test_that("a topTable written to CSV reads IDs from its unnamed first column", {
  tbl <- limma_table()
  csv <- cbind(X = rownames(tbl), tbl)
  rownames(csv) <- NULL
  d <- as_da(csv, contrast = "Aging")
  expect_identical(d$protein, rownames(tbl))
})

test_that("proteoDA per-contrast CSVs without t rank by signed -log10(p)", {
  tbl <- data.frame(
    uniprot_id = c("P31040", "Q9UBK2"), gene_symbol = c("=\"SDHA\"", "=\"PPARGC1A\""),
    logFC = c(1, -2), P.Value = c(0.01, 0.001), adj.P.Val = c(0.05, 0.01),
    average_intensity = c(20, 19), sig.FDR = c(1, -1)
  )
  d <- as_da(tbl, contrast = "Training")
  expect_identical(d$gene, c("SDHA", "PPARGC1A"))
  expect_true(all(is.na(d$t)))
  expect_identical(unique(d$rank_stat), "signed -log10(p)")
  expect_equal(d$rank, c(2, -3))
})

test_that("MSstats results keep their Label contrasts", {
  tbl <- data.frame(
    Protein = c("P1", "P2", "P1"), Label = c("B-A", "B-A", "C-A"),
    log2FC = c(1, -1, 0.5), SE = 0.2, Tvalue = c(5, -5, 2.5), DF = 10,
    pvalue = c(1e-4, 1e-4, 0.03), adj.pvalue = c(1e-3, 1e-3, 0.06), issue = NA
  )
  d <- as_da(tbl)
  expect_identical(d$contrast, c("B-A", "B-A", "C-A"))
  expect_equal(d$t, c(5, -5, 2.5))
  expect_true(all(is.na(d$abundance)))
})

test_that("proDA, msqrob2 and prolfQua columns are recognised", {
  prod <- data.frame(
    name = c("P1", "P2"), pval = c(0.01, 0.2), adj_pval = c(0.02, 0.3), diff = c(1, -0.5),
    t_statistic = c(3, -1), se = 0.3, df = 8, avg_abundance = c(21, 19), n_approx = 6, n_obs = 6
  )
  msq <- data.frame(
    logFC = c(1, -0.5), se = 0.3, df = 8, t = c(3, -1), pval = c(0.01, 0.2),
    adjPval = c(0.02, 0.3), row.names = c("P1", "P2")
  )
  prolf <- data.frame(
    protein_Id = c("P1", "P2"), contrast = "B_vs_A", modelName = "m", avgAbd = c(21, 19),
    diff = c(1, -0.5), FDR = c(0.02, 0.3), statistic = c(3, -1), p.value = c(0.01, 0.2)
  )
  for (d in list(as_da(prod, contrast = "B_vs_A"), as_da(msq, contrast = "B_vs_A"), as_da(prolf))) {
    expect_identical(d$protein, c("P1", "P2"))
    expect_equal(d$logFC, c(1, -0.5))
    expect_equal(d$t, c(3, -1))
    expect_equal(d$padj, c(0.02, 0.3))
  }
})

test_that("ProtRank tables without p rank by signed -log10(padj)", {
  tbl <- data.frame(
    row = c("P1", "P2"), LFC = c(1, -1), `rank score` = c(9, 8), FDR = c(0.01, 0.1),
    check.names = FALSE
  )
  d <- as_da(tbl, contrast = "B_vs_A")
  expect_identical(unique(d$rank_stat), "signed -log10(padj)")
  expect_equal(d$rank, c(2, -1))
})

test_that("explicit column names override detection", {
  tbl <- data.frame(id = c("P1", "P2"), fc = c(1, -1), stat = c(3, -3), q = c(0.01, 0.02))
  d <- as_da(tbl, contrast = "A", protein = "id", logfc = "fc", t = "stat", padj = "q")
  expect_identical(d$protein, c("P1", "P2"))
  expect_equal(d$rank, c(3, -3))
})

test_that("a named list of tables takes its names as contrasts", {
  d <- as_da(list(Young = limma_table(), Old = limma_table()))
  expect_identical(unique(d$contrast), c("Young", "Old"))
  expect_identical(nrow(d), 6L)
})

test_that("a p-value of zero gives a finite rank", {
  tbl <- data.frame(protein = "P1", logFC = 1, p = 0)
  expect_true(is.finite(as_da(tbl, contrast = "A")$rank))
})

test_that("unusable tables are refused with a reason", {
  expect_error(as_da(limma_table()), "contrast", class = "enrichVolcano_input_error")
  ftest <- data.frame(BvA = 1, CvA = 2, AveExpr = 20, F = 9, P.Value = 0.01, adj.P.Val = 0.02, row.names = "P1")
  expect_error(as_da(ftest, contrast = "A"), "F-test", class = "enrichVolcano_input_error")
  wide <- data.frame(uniprot_id = "P1", logFC_A = 1, t_A = 3, P.Value_A = 0.01)
  expect_error(as_da(wide, contrast = "A"), "wide", class = "enrichVolcano_input_error")
  expect_error(as_da(wide), "wide", class = "enrichVolcano_input_error")
  no_stat <- data.frame(protein = "P1", logFC = 1)
  expect_error(as_da(no_stat, contrast = "A"), "t, p or padj", class = "enrichVolcano_column_error")
  no_fc <- data.frame(protein = "P1", t = 3)
  expect_error(as_da(no_fc, contrast = "A"), "logFC", class = "enrichVolcano_column_error")
  no_id <- data.frame(logFC = 1, t = 3)
  expect_error(as_da(no_id, contrast = "A"), "protein", class = "enrichVolcano_column_error")
})

test_that("impossible values and duplicate rows are refused", {
  tbl <- data.frame(protein = c("P1", "P2"), logFC = c(1, 2), p = c(0.1, 1.4))
  expect_error(as_da(tbl, contrast = "A"), "\\[0, 1\\]", class = "enrichVolcano_data_error")
  dup <- data.frame(protein = c("P1", "P1"), logFC = c(1, 2), p = c(0.1, 0.2))
  expect_error(as_da(dup, contrast = "A"), "duplicate", class = "enrichVolcano_data_error")
})

test_that("columns as_da does not use are kept for later", {
  tbl <- data.frame(protein = c("P1", "P2"), logFC = c(1, -1), P.Value = c(0.01, 0.02), pi_eq2 = c(0.1, 0.3))
  d <- as_da(tbl, contrast = "A", species = NULL)
  expect_equal(d$pi_eq2, c(0.1, 0.3))
  expect_false("P.Value" %in% names(d))
})
