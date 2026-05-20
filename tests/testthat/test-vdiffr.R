# Visual regression for the composite. The committed SVG baseline is rendered on
# the maintainer's machine; font metrics differ across platforms, so the
# comparison runs locally only (skip_on_ci) and guards against accidental
# changes to the composite's geometry or layers.

test_that("the volcano-in-ring composite is visually stable", {
  skip_if_not_installed("vdiffr")
  skip_on_ci()  # baseline is platform-specific; compare only on the reference machine

  genes <- paste0("G", seq_len(40))
  volc <- tibble::tibble(
    gene = genes,
    contrast = "ctr",
    logFC = seq(-3, 3, length.out = 40),
    P.Value = rep(c(1e-4, 1e-3, 1e-2, 0.2, 0.6), times = 8),
    pi_eq2 = rep(c(0.001, 0.01, 0.03, 0.4, 0.8), times = 8)
  )
  enr <- tibble::tibble(
    contrast = "ctr", database = "test",
    pathway = c("HALLMARK_GLYCOLYSIS", "HALLMARK_HYPOXIA",
                "HALLMARK_MTORC1_SIGNALING", "HALLMARK_APOPTOSIS"),
    padj = c(0.001, 0.01, 0.02, 0.04),
    NES = c(2.4, 1.7, -1.5, -2.2),
    size = c(40L, 35L, 30L, 25L),
    leading_edge = c(paste0("G", 1:8, collapse = ";"),
                     paste0("G", 5:12, collapse = ";"),
                     paste0("G", 28:35, collapse = ";"),
                     paste0("G", 33:40, collapse = ";")),
    mode = "fgsea",
    direction = c("up", "up", "down", "down")
  )

  p <- ev_volcano_ring(volc, enr, title = "ctr", logfc_threshold = log2(1.5))
  vdiffr::expect_doppelganger("volcano-in-ring composite", p)
})
