toy_study <- function() suppressMessages(read_study(test_path("fixtures", "toy_study")))

toy_sets <- function() {
  list(
    Hallmark = list(
      OXPHOS = c("SDHA", "SDHB", "COXFA4", "COX4I1", "ATP5F1A", "CS", "IDH3A", "MDH2", "CYCS", "UQCRC1"),
      GLYCOLYSIS = c("HK2", "PKM", "LDHA", "ENO1", "GAPDH", "PFKM", "ALDOA", "TPI1", "PGK1", "PGAM1")
    ),
    Mixed = list(
      HALF_A = c("SDHA", "SDHB", "COXFA4", "COX4I1", "ATP5F1A", "HK2", "PKM", "LDHA", "ENO1", "GAPDH"),
      HALF_B = c("CS", "IDH3A", "MDH2", "CYCS", "UQCRC1", "PFKM", "ALDOA", "TPI1", "PGK1", "PGAM1")
    )
  )
}

run_fgsea <- function(...) {
  skip_if_not_installed("fgsea")
  withr::local_seed(1)
  suppressMessages(run_enrichment(..., tests = "fgsea", min_size = 5))
}

test_that("fgsea through run_enrichment equals a direct fgsea call", {
  skip_if_not_installed("org.Hs.eg.db")
  s <- toy_study()
  x <- run_fgsea(s, toy_sets())
  expect_named(x, "fgsea")
  r <- x$fgsea@results
  stats <- stats::setNames(s$da$rank, s$da$gene)
  expected <- withr::with_seed(1, fgsea::fgsea(toy_sets()$Hallmark, stats, minSize = 5, maxSize = 500))
  got <- r[r$database == "Hallmark", ]
  expect_equal(got$score[match(expected$pathway, got$term)], expected$NES)
  expect_equal(got$p[match(expected$pathway, got$term)], expected$pval)
  expect_identical(x$fgsea@metadata$enrichment_test, "fgsea")
  expect_identical(x$fgsea@metadata$ranking, "t")
})

test_that("each collection is corrected on its own", {
  skip_if_not_installed("org.Hs.eg.db")
  s <- toy_study()
  r <- run_fgsea(s, toy_sets())$fgsea@results
  for (db in c("Hallmark", "Mixed")) {
    part <- r[r$database == db, ]
    expect_equal(part$padj, stats::p.adjust(part$p, "BH"))
  }
})

test_that("gene-set versions are copied into the metadata", {
  skip_if_not_installed("org.Hs.eg.db")
  sets <- toy_sets()
  attr(sets, "versions") <- list(species = "Homo sapiens", msigdbr = "26.1.0")
  x <- run_fgsea(toy_study(), sets)
  expect_identical(x$fgsea@metadata$gene_sets$msigdbr, "26.1.0")
})

test_that("a flat list of sets is one collection", {
  skip_if_not_installed("org.Hs.eg.db")
  x <- run_fgsea(toy_study(), toy_sets()$Hallmark)
  expect_identical(unique(x$fgsea@results$database), "gene sets")
})

test_that("a DA table alone is enough for fgsea", {
  skip_if_not_installed("org.Hs.eg.db")
  s <- toy_study()
  x <- run_fgsea(s$da, toy_sets())
  expect_identical(nrow(x$fgsea@results), 4L)
})

test_that("one protein per gene: the most abundant, else the first accession", {
  da <- data.frame(
    protein = c("A2", "A1", "B1", "B2"), gene = c("G1", "G1", "G2", "G2"), contrast = "C",
    logFC = c(2, -2, 1, -1), t = c(5, -5, 2, -2), p = 0.01, padj = 0.02,
    abundance = c(10, 20, NA, NA)
  )
  kept <- one_per_gene(da)
  expect_identical(sort(kept$protein), c("A1", "B1"))
})

test_that("matrix row means fill in a missing abundance", {
  skip_if_not_installed("org.Hs.eg.db")
  s <- toy_study()
  s$da$abundance <- NA_real_
  filled <- fill_abundance(s$da, s$matrix)
  expect_equal(filled$abundance, unname(rowMeans(s$matrix)[s$da$protein]))
})

test_that("bad inputs are refused", {
  skip_if_not_installed("org.Hs.eg.db")
  s <- toy_study()
  expect_error(run_enrichment(data.frame(x = 1), toy_sets()), class = "enrichVolcano_input_error")
  expect_error(run_enrichment(s, list(1, 2)), class = "enrichVolcano_input_error")
  expect_error(run_enrichment(s, toy_sets(), tests = "gsva"), class = "rlang_error")
  expect_error(
    run_enrichment(s, list(Tiny = list(A = c("SDHA", "CS"))), tests = "fgsea"),
    "no gene set",
    class = "enrichVolcano_input_error"
  )
})
