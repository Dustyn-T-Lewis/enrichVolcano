test_that("ev_volcano_ring returns an enrichVolcano ggplot with arcs and labels", {
  set.seed(11)
  volc <- tibble::tibble(
    gene = paste0("G", 1:200),
    contrast = "ctr",
    logFC = c(rnorm(30, 3, 0.5), rnorm(140), rnorm(30, -3, 0.5)),
    P.Value = c(runif(30, 0, 1e-3), runif(140, 0.1, 1), runif(30, 0, 1e-3)),
    pi_eq2 = c(runif(30, 0, 0.01), runif(140, 0.2, 1), runif(30, 0, 0.01))
  )
  enr <- tibble::tibble(
    contrast = "ctr", database = "test",
    pathway = c("HALLMARK_MTORC1_SIGNALING", "HALLMARK_GLYCOLYSIS",
                "HALLMARK_HYPOXIA"),
    padj = c(0.001, 0.01, 0.04),
    NES = c(2.1, 1.6, -1.8),
    size = c(40, 35, 30),
    leading_edge = c(paste0("G", 1:15, collapse = ";"),
                     paste0("G", 5:20, collapse = ";"),
                     paste0("G", 171:185, collapse = ";"))
  )
  p <- ev_volcano_ring(volc, enr, title = "ctr")
  expect_s3_class(p, "enrichVolcano")
  expect_s3_class(p, "ggplot")
})

test_that("ev_volcano_ring degrades to volcano only when no pathways", {
  volc <- tibble::tibble(
    gene = paste0("G", 1:50), contrast = "ctr",
    logFC = rnorm(50), P.Value = runif(50), pi_eq2 = runif(50)
  )
  empty <- enr <- tibble::tibble(
    contrast = character(0), database = character(0), pathway = character(0),
    padj = double(0), NES = double(0), size = integer(0),
    leading_edge = character(0)
  )
  p <- ev_volcano_ring(volc, empty, title = "ctr")
  expect_s3_class(p, "enrichVolcano")
})

test_that("ev_clean_label strips prefixes, expands acronyms, handles MitoCarta", {
  expect_equal(enrichVolcano:::ev_clean_label("HALLMARK_MTORC1_SIGNALING"),
               "MTORC1\nSignaling")
  expect_match(enrichVolcano:::ev_clean_label("HALLMARK_E2F_TARGETS"), "E2F")
  mc <- enrichVolcano:::ev_clean_label("MITOCARTA_OXPHOS__CV_SUBUNITS")
  expect_match(mc, "Complex V")
})

test_that("ev_ring_geometry tolerates NA padj/NES from fgsea", {
  enr <- tibble::tibble(
    contrast = "ctr", database = "test",
    pathway = c("P1", "P2"), padj = c(NA, 0.01), NES = c(NA, 1.5),
    size = c(10, 10), leading_edge = c("G1;G2", "G3;G4")
  )
  g <- enrichVolcano:::ev_ring_geometry(enr)
  expect_true(all(is.finite(g$arc_r1_var)))
})
