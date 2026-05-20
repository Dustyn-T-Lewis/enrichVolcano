test_that("enrich_volcano runs end-to-end for human / mouse / rat on hallmark", {
  skip_if_not_installed("msigdbr")

  for (species in c("human", "mouse", "rat")) {
    set.seed(1)
    data <- tibble::tibble(
      gene = sample(c("TP53", "BRCA1", "MYC", "AKT1", "EGFR", "KRAS",
                      "PTEN", "PIK3CA", "RB1", "CDKN2A"), 30, replace = TRUE),
      contrast = "ctr",
      logFC = rnorm(30, 0, 2),
      P.Value = runif(30, 0, 0.1)
    )
    p <- enrich_volcano(
      data, contrast = "ctr", species = species,
      databases = "hallmark",
      enrich_mode = "fgsea",
      enrich_padj = 1.0
    )
    expect_true(inherits(p, "enrichVolcano"),
                info = paste("species =", species))
  }
})
