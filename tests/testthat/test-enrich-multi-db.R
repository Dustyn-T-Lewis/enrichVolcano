test_that("enrich_volcano accepts multiple databases and binds results", {
  skip_if_not_installed("msigdbr")

  # 80 real HGNC symbols drawn from the top of the Hallmark frequency table so
  # that fgsea has enough unique genes to exceed the default min_size = 15 for
  # at least one pathway in each of hallmark and reactome.  The original
  # 10-gene pool produced only 9 unique ranked entries after deduplication,
  # below the size filter, causing fgsea to return zero rows.
  gene_pool <- c(
    "CDKN1A", "IL6", "MYC", "SERPINE1", "BTG2", "CA2", "CD44", "IDH1",
    "IRF1", "PLAUR", "RETSAT", "ATF3", "BHLHE40", "BMP2", "CCND1", "CCND2",
    "CXCL10", "ELOVL5", "FAS", "SOD1", "TAP1", "TGFB1", "TIMP1", "VEGFA",
    "ALDOA", "AURKA", "B4GALT1", "CAT", "CCND3", "CD36", "CDK1", "CDKN1B",
    "CSF1", "ECH1", "EGFR", "ENO2", "F3", "FOS", "GCH1", "GPX4",
    "HMOX1", "ICAM1", "ID2", "IDI1", "IL1B", "IL4R", "ISG20", "JUN",
    "LDHA", "LIF", "PPP1R15A", "SOD2", "TNFAIP3", "TOP2A", "ABCA1", "ALDH9A1",
    "APP", "BTG1", "CASP1", "CASP3", "CCL2", "CCL5", "CD9", "CDK4",
    "CFB", "CLU", "CRAT", "CXCL11", "CXCR4", "DDIT3", "DDIT4", "DHCR24",
    "DHCR7", "ERO1A", "GADD45A", "GADD45B", "GCLC", "GLRX", "GPX3", "GSR"
  )

  set.seed(2)
  data <- tibble::tibble(
    gene = sample(gene_pool, 200, replace = TRUE),
    contrast = "ctr",
    logFC = rnorm(200, 0, 2),
    P.Value = runif(200, 0, 0.1)
  )
  p <- suppressWarnings(enrich_volcano(
    data, contrast = "ctr",
    databases = c("hallmark", "reactome"),
    enrich_mode = "fgsea", enrich_padj = 1.0
  ))
  enrich <- attr(p, "ev_data")$enrichment
  expect_true(all(c("hallmark", "reactome") %in% unique(enrich$database)))
})

test_that("ev_collapse cutoff sweep produces monotonically decreasing kept count", {
  set.seed(3)
  enrich <- tibble::tibble(
    contrast = "ctr", database = "test",
    pathway = paste0("P", 1:10),
    pval = runif(10), padj = runif(10), NES = rnorm(10), size = 10,
    leading_edge = c(
      "G1;G2;G3", "G1;G2;G3", "G1;G2;G4",
      "G5;G6;G7", "G5;G6;G7", "G8;G9;G10",
      "G1;G2;G3;G4", "G11;G12;G13", "G11;G12;G14", "G15;G16;G17"
    ),
    mode = "fgsea", direction = "up"
  )
  kept_strict  <- sum(ev_collapse(enrich, method = "jaccard", cutoff = 0.2)$dedup_kept)
  kept_default <- sum(ev_collapse(enrich, method = "jaccard", cutoff = 0.5)$dedup_kept)
  kept_loose   <- sum(ev_collapse(enrich, method = "jaccard", cutoff = 0.9)$dedup_kept)
  expect_lte(kept_strict, kept_default)
  expect_lte(kept_default, kept_loose)
})

test_that("ev_collapse respects scope parameter", {
  enrich <- tibble::tibble(
    contrast = rep("ctr", 6), database = rep(c("db1", "db2"), each = 3),
    pathway = paste0("P", 1:6),
    pval = runif(6), padj = runif(6), NES = rnorm(6), size = 10,
    leading_edge = rep("G1;G2;G3", 6),  # all identical -> maximally redundant
    mode = "fgsea", direction = "up"
  )
  # sig_threshold = NA preserves original intent: this test isolates the
  # scope behavior on a maximally redundant fixture; random padj would
  # otherwise be filtered out by the default sig gate.
  kept_within <- sum(ev_collapse(enrich, method = "jaccard",
                                  cutoff = 0.2, scope = "within_db",
                                  sig_threshold = NA)$dedup_kept)
  kept_global <- sum(ev_collapse(enrich, method = "jaccard",
                                  cutoff = 0.2, scope = "global",
                                  sig_threshold = NA)$dedup_kept)
  # within_db: 1 per db = 2 kept; global: 1 total
  expect_equal(kept_within, 2)
  expect_equal(kept_global, 1)
})
