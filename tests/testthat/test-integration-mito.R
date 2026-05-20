test_that("Mito end-to-end: rat-shaped data works", {
  skip_on_cran()
  fx <- file.path(test_path("fixtures"), "mito_tidy.rds")
  skip_if_not(file.exists(fx), "Mito fixture not built")
  data <- readRDS(fx)
  skip_if(nrow(data) < 50, "Mito fixture too small")
  ctrs <- unique(data$contrast)
  # MitoCarta-shaped mock paths (top up/down genes, simulating MitoPathway labels)
  paths <- list(
    OXPHOS = head(data$gene[data$logFC > 0], 30),
    TCA_Cycle = tail(data$gene[data$logFC > 0], 30),
    Apoptosis = head(data$gene[data$logFC < 0], 20)
  )
  p <- suppressWarnings(enrich_volcano(
    data, contrast = ctrs[1],
    species = "rat",
    databases = list(mitocarta3 = paths),
    p_method = "pi_eq2", p_adjust = "BH",
    enrich_mode = "fgsea", enrich_padj = 1.0
  ))
  expect_s3_class(p, "enrichVolcano")
})
