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

# The gene-level matrix and design the package should build for the toy study.
toy_expected <- function(fixed = FALSE) {
  s <- toy_study()
  m <- s$matrix
  rownames(m) <- s$da$gene[match(rownames(m), s$da$protein)]
  w <- s$weights
  rownames(w) <- rownames(m)
  f <- if (fixed) ~ 0 + group + subject else ~ 0 + group
  design <- model.matrix(f, data.frame(group = factor(s$samples$group), subject = factor(s$samples$subject)))
  colnames(design)[1:2] <- levels(factor(s$samples$group))
  cm <- limma::makeContrasts(contrasts = "Post - Pre", levels = design)
  index <- limma::ids2indices(toy_sets()$Hallmark, rownames(m))
  list(m = m, w = w, design = design, contrast = cm[, 1], index = index, subject = s$samples$subject)
}

run_limma <- function(tests, ...) {
  skip_if_not_installed("limma")
  skip_if_not_installed("org.Hs.eg.db")
  suppressMessages(run_enrichment(toy_study(), toy_sets(), tests = tests, min_size = 5, ...))
}

test_that("blocked fry equals a direct limma::fry call", {
  e <- toy_expected()
  corr <- limma::duplicateCorrelation(e$m, e$design, block = e$subject, weights = e$w)$consensus.correlation
  expected <- limma::fry(e$m, e$index, e$design, e$contrast,
    block = e$subject, correlation = corr, weights = e$w, sort = "none"
  )
  x <- run_limma("fry")$fry
  r <- x@results[x@results$database == "Hallmark", ]
  expect_equal(r$p[match(rownames(expected), r$term)], expected$PValue)
  expect_identical(x@metadata$enrichment_test, "fry")
})

test_that("blocked camera runs cameraPR on the blocked fit's t", {
  e <- toy_expected()
  corr <- limma::duplicateCorrelation(e$m, e$design, block = e$subject, weights = e$w)$consensus.correlation
  fit <- limma::lmFit(e$m, e$design, block = e$subject, correlation = corr, weights = e$w)
  fit <- limma::eBayes(limma::contrasts.fit(fit, e$contrast))
  expected <- limma::cameraPR(fit$t[, 1], e$index, sort = FALSE)
  x <- run_limma("camera")$camera
  r <- x@results[x@results$database == "Hallmark", ]
  expect_equal(r$p[match(rownames(expected), r$term)], expected$PValue)
  expect_identical(x@metadata$enrichment_test, "cameraPR")
})

test_that("with a fixed subject effect, camera runs on the full design", {
  e <- toy_expected(fixed = TRUE)
  expected <- limma::camera(e$m, e$index, e$design, e$contrast, weights = e$w, sort = FALSE)
  x <- run_limma("camera", subject_effect = "fixed")$camera
  r <- x@results[x@results$database == "Hallmark", ]
  expect_equal(r$p[match(rownames(expected), r$term)], expected$PValue)
  expect_identical(x@metadata$enrichment_test, "camera")
})

test_that("the default runs all three tests when the study has a matrix", {
  skip_if_not_installed("fgsea")
  x <- run_limma(c("fgsea", "camera", "fry"))
  expect_named(x, c("fgsea", "camera", "fry"))
})

test_that("a DA-only study falls back to fgsea unless camera or fry is asked for", {
  skip_if_not_installed("fgsea")
  skip_if_not_installed("org.Hs.eg.db")
  da <- toy_study()$da
  expect_message(x <- run_enrichment(da, toy_sets(), min_size = 5), class = "enrichVolcano_fgsea_only")
  expect_named(x, "fgsea")
  expect_error(run_enrichment(da, toy_sets(), tests = "fry"), "matrix", class = "enrichVolcano_input_error")
})

test_that("design_matrix builds groups, covariates and fixed subjects", {
  samples <- data.frame(
    sample = paste0("S", 1:6), group = rep(c("A", "B"), 3), subject = rep(c("P1", "P2", "P3"), each = 2),
    sex = rep(c("F", "M", "F"), each = 2)
  )
  d <- design_matrix(samples, "block", covariates = "sex")
  expect_identical(colnames(d), c("A", "B", "sexM"))
  d <- design_matrix(samples, "fixed", covariates = NULL)
  expect_identical(colnames(d), c("A", "B", "subjectP2", "subjectP3"))
  expect_error(design_matrix(samples, "block", covariates = "age"), "age", class = "enrichVolcano_column_error")
  samples$subject <- NULL
  expect_error(design_matrix(samples, "fixed", NULL), "subject", class = "enrichVolcano_column_error")
})

test_that("blocking without a subject column says the samples are independent", {
  skip_if_not_installed("limma")
  skip_if_not_installed("org.Hs.eg.db")
  s <- toy_study()
  s$samples$subject <- NULL
  expect_message(
    run_enrichment(s, toy_sets(), tests = "fry", min_size = 5),
    class = "enrichVolcano_no_blocking"
  )
})

test_that("a matrix with missing values is refused for fry and camera", {
  skip_if_not_installed("limma")
  skip_if_not_installed("org.Hs.eg.db")
  s <- toy_study()
  s$matrix[1, 1] <- NA
  expect_error(run_enrichment(s, toy_sets(), tests = "fry"), "1 missing", class = "enrichVolcano_data_error")
})

test_that("camera and fry say so when no set fits the size window", {
  skip_if_not_installed("limma")
  skip_if_not_installed("org.Hs.eg.db")
  expect_error(
    suppressMessages(run_enrichment(toy_study(), toy_sets(), tests = "fry", min_size = 50)),
    "no gene set",
    class = "enrichVolcano_input_error"
  )
})

test_that("camera can estimate the inter-gene correlation instead of fixing it", {
  e <- toy_expected(fixed = TRUE)
  expected <- limma::camera(e$m, e$index, e$design, e$contrast, weights = e$w, inter.gene.cor = NA, sort = FALSE)
  x <- run_limma("camera", subject_effect = "fixed", inter_gene_cor = NA)$camera
  r <- x@results[x@results$database == "Hallmark", ]
  expect_equal(r$p[match(rownames(expected), r$term)], expected$PValue)
  expect_true(is.na(x@metadata$inter_gene_cor))
})

test_that("a blocked design cannot estimate the correlation, and says so", {
  skip_if_not_installed("limma")
  skip_if_not_installed("org.Hs.eg.db")
  expect_error(
    suppressMessages(run_enrichment(toy_study(), toy_sets(), tests = "camera", min_size = 5, inter_gene_cor = NA)),
    "cameraPR",
    class = "enrichVolcano_param_error"
  )
})

test_that("tied ranks are reported once as a message, not as fgsea warnings", {
  skip_if_not_installed("fgsea")
  skip_if_not_installed("org.Hs.eg.db")
  da <- toy_study()$da
  da$rank[1:2] <- da$rank[3]
  # Older fgsea also warns that tiny toy sets have overestimated p-values.
  suppressWarnings(expect_no_warning(
    expect_message(
      run_enrichment(da, toy_sets(), tests = "fgsea", min_size = 5),
      "3 of 20",
      class = "enrichVolcano_rank_ties"
    ),
    message = "ties"
  ))
})

test_that("a missing matrix value does not hide a protein's abundance", {
  skip_if_not_installed("org.Hs.eg.db")
  s <- toy_study()
  s$da$abundance <- NA_real_
  s$matrix[1, 1] <- NA
  filled <- fill_abundance(s$da, s$matrix)
  expect_false(anyNA(filled$abundance))
})
