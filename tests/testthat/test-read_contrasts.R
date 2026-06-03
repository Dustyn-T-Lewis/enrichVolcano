fx <- function() testthat::test_path("fixtures", "contrasts")

test_that("ev_read_contrasts reads a folder, filename = contrast", {
  out <- suppressMessages(ev_read_contrasts(fx(), species = "human"))
  expect_setequal(unique(out$contrast),
                  c("Cancer_vs_Healthy", "Treated_vs_Control"))
  expect_equal(nrow(out), 5)
  expect_true(all(c("uniprot", "symbol", "gene", "logFC",
                    "P.Value", "adj.P.Val") %in% colnames(out)))
  expect_identical(out$symbol[out$uniprot == "P04637"][1], "TP53")
  rep <- ev_idmap_report(out)
  expect_equal(rep$n_unmapped, 1)   # Q0FAKE9
})

test_that("ev_read_contrasts computes adjustment when none supplied", {
  out <- suppressMessages(ev_read_contrasts(fx(), species = "human"))
  expect_false(any(is.na(out$adj.P.Val)))
})

test_that("ev_read_contrasts accepts an explicit file vector", {
  files <- list.files(fx(), pattern = "\\.csv$", full.names = TRUE)
  out <- suppressMessages(ev_read_contrasts(files, species = "human"))
  expect_setequal(unique(out$contrast),
                  c("Cancer_vs_Healthy", "Treated_vs_Control"))
})

test_that("ev_read_contrasts aborts on an empty folder", {
  d <- withr::local_tempdir()
  expect_error(ev_read_contrasts(d, species = "human"),
               class = "ev_no_contrast_files")
})

test_that("ev_read_contrasts aborts on duplicate contrast names", {
  d <- withr::local_tempdir()
  sub1 <- file.path(d, "a")
  sub2 <- file.path(d, "b")
  dir.create(sub1)
  dir.create(sub2)
  src <- list.files(fx(), pattern = "Cancer", full.names = TRUE)
  file.copy(src, file.path(sub1, "X.csv"))
  file.copy(src, file.path(sub2, "X.csv"))
  expect_error(
    suppressMessages(
      ev_read_contrasts(c(file.path(sub1, "X.csv"), file.path(sub2, "X.csv")),
                        species = "human")),
    class = "ev_dup_contrast"
  )
})

test_that("ev_read_contrasts output flows through enrich_volcano end-to-end", {
  res <- suppressMessages(
    ev_read_contrasts(testthat::test_path("fixtures", "contrasts"),
                      species = "human"))
  paths <- list(S1 = c("TP53", "EGFR"), S2 = c("BRCA1", "TP53"))
  p <- suppressWarnings(suppressMessages(enrich_volcano(
    res, contrast = "Cancer_vs_Healthy",
    databases = list(test = paths), enrich_padj = 1,
    label_mode = "top_per_direction", label_n = 2
  )))
  expect_s3_class(p, "enrichVolcano")
})

test_that("raw UniProt data flows through enrich_volcano with species set", {
  df <- data.frame(
    UniProt = c("P04637", "P00533", "P38398"),
    logFC = c(2, -1.5, 1.2), P.Value = c(1e-4, 1e-3, 1e-2),
    contrast = "C1", stringsAsFactors = FALSE)
  paths <- list(S1 = c("TP53", "EGFR"), S2 = "BRCA1")
  p <- suppressWarnings(suppressMessages(enrich_volcano(
    df, contrast = "C1", species = "human",
    databases = list(test = paths), enrich_padj = 1)))
  expect_s3_class(p, "enrichVolcano")
})
