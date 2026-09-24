source(test_path("fixtures/make_toy.R"))

flags <- function(x) {
  r <- x@results
  stats::setNames(r$dedup_status, r$term)
}

test_that("set similarities match hand-computed values", {
  a <- paste0("G", 1:10)
  b <- paste0("G", c(1:8, 11:12))
  c4 <- paste0("G", 1:4)
  expect_equal(jaccard(a, b), 8 / 12)
  expect_equal(overlap_coefficient(a, b), 8 / 10)
  expect_equal(combined_similarity(a, b), 0.5 * 8 / 12 + 0.5 * 8 / 10)
  expect_equal(jaccard(a, c4), 4 / 10)
  expect_equal(overlap_coefficient(a, c4), 1)
  expect_identical(jaccard(a, paste0("X", 1:3)), 0)
  expect_identical(overlap_coefficient(a, character(0)), 0)
  expect_equal(jaccard(c(a, a), b), 8 / 12)
})

test_that("the default combined rule at 0.375 flags B and C under A", {
  x <- dedup_terms(make_toy_dedup_enrichment(), make_toy_gene_sets())
  f <- flags(x)
  expect_identical(
    unname(f[c("SET_A", "SET_B", "SET_C", "SET_D")]),
    c("kept", "redundant", "redundant", "kept")
  )
  r <- x@results
  expect_identical(r$merged_into[r$term == "SET_C"], "SET_A")
  expect_equal(r$similarity[r$term == "SET_B"], 0.5 * 8 / 12 + 0.5 * 8 / 10)
  expect_equal(r$similarity[r$term == "SET_C"], 0.5 * 4 / 10 + 0.5 * 1)
  expect_identical(
    x@metadata$dedup,
    list(method = "enrichmentmap", similarity = "combined", cutoff = 0.375, term_threshold = 0.05)
  )
})

test_that("Jaccard at 0.5 keeps C, whose Jaccard with A is only 0.4", {
  x <- dedup_terms(make_toy_dedup_enrichment(), make_toy_gene_sets(), similarity = "jaccard")
  f <- flags(x)
  expect_identical(unname(f[c("SET_A", "SET_B", "SET_C")]), c("kept", "redundant", "kept"))
  expect_equal(x@results$similarity[x@results$term == "SET_B"], 8 / 12)
  expect_identical(x@metadata$dedup$cutoff, 0.5)
})

test_that("an explicit cutoff overrides the rule's default", {
  x <- dedup_terms(make_toy_dedup_enrichment(), make_toy_gene_sets(), similarity = "jaccard", cutoff = 0.9)
  expect_identical(unname(flags(x)["SET_B"]), "kept")
})

test_that("identical sets in different databases are never merged", {
  x <- dedup_terms(make_toy_dedup_enrichment(), make_toy_gene_sets())
  expect_identical(unname(flags(x)["SET_E"]), "kept")
})

test_that("non-significant rows are left unflagged", {
  x <- dedup_terms(make_toy_dedup_enrichment(), make_toy_gene_sets())
  r <- x@results
  expect_true(is.na(r$dedup_status[r$term == "SET_F"]))
  expect_true(is.na(r$merged_into[r$term == "SET_F"]))
})

test_that("each contrast is walked in its own padj order", {
  x <- make_toy_dedup_enrichment()
  r1 <- x@results
  r2 <- r1
  r2$contrast <- "L"
  r2$padj[r2$term == "SET_B"] <- 0.0001
  x@results <- rbind(r1, r2)
  r <- dedup_terms(x, make_toy_gene_sets())@results
  in_l <- r[r$contrast == "L", ]
  expect_identical(in_l$dedup_status[in_l$term == "SET_B"], "kept")
  expect_identical(in_l$merged_into[in_l$term == "SET_A"], "SET_B")
})

test_that("p-values and row count are never changed", {
  x <- make_toy_dedup_enrichment()
  y <- dedup_terms(x, make_toy_gene_sets())
  expect_identical(y@results$padj, x@results$padj)
  expect_identical(nrow(y@results), nrow(x@results))
})

test_that("terms without a gene set are kept and counted", {
  sets <- make_toy_gene_sets()
  sets$SET_B <- NULL
  expect_message(x <- dedup_terms(make_toy_dedup_enrichment(), sets), "1 term")
  expect_identical(unname(flags(x)["SET_B"]), "kept")
})

test_that("earlier dedup columns are replaced, not duplicated", {
  x <- make_toy_dedup_enrichment()
  r <- x@results
  r$dedup_status <- "kept"
  r$overlap_jaccard <- 0.9
  x@results <- r
  y <- suppressMessages(dedup_terms(x, make_toy_gene_sets()))
  expect_false("overlap_jaccard" %in% names(y@results))
  expect_identical(sum(names(y@results) == "dedup_status"), 1L)
})

test_that("bad inputs are refused with classed errors", {
  x <- make_toy_dedup_enrichment()
  sets <- make_toy_gene_sets()
  expect_error(dedup_terms(x@results, sets), class = "enrichVolcano_input_error")
  expect_error(dedup_terms(x, unname(sets)), class = "enrichVolcano_input_error")
  expect_error(dedup_terms(x, "SET_A"), class = "enrichVolcano_input_error")
  expect_error(dedup_terms(x, sets, cutoff = 1.5), class = "enrichVolcano_param_error")
  expect_error(dedup_terms(x, sets, cutoff = c(0.3, 0.4)), class = "enrichVolcano_param_error")
  expect_error(dedup_terms(x, sets, similarity = "cosine"), class = "rlang_error")
})

fgsea_example <- function() {
  skip_if_not_installed("fgsea")
  skip_if_not_installed("data.table")
  pathways <- NULL
  ranks <- NULL
  utils::data("examplePathways", package = "fgsea", envir = environment())
  utils::data("exampleRanks", package = "fgsea", envir = environment())
  pathways <- get("examplePathways")[1:40]
  ranks <- get("exampleRanks")
  withr::local_seed(1)
  res <- fgsea::fgsea(pathways, ranks, minSize = 15, maxSize = 500)
  list(res = res, pathways = pathways, ranks = ranks)
}

test_that("collapse_pathways reproduces fgsea::collapsePathways", {
  fx <- fgsea_example()
  x <- suppressMessages(as_enrichment(list(A = fx$res)))

  sig <- fx$res[fx$res$padj < 0.05, ]
  sig <- sig[order(sig$pval), ]
  expected <- withr::with_seed(2, fgsea::collapsePathways(sig, fx$pathways, fx$ranks))
  y <- withr::with_seed(2, dedup_terms(x, fx$pathways, method = "collapse_pathways", stats = list(A = fx$ranks)))

  r <- y@results
  in_sig <- r$term %in% sig$pathway
  expect_setequal(r$term[in_sig & r$dedup_status == "kept"], expected$mainPathways)
  parents <- expected$parentPathways[!is.na(expected$parentPathways)]
  expect_identical(
    r$merged_into[match(names(parents), r$term)],
    unname(parents)
  )
  expect_true(all(is.na(r$dedup_status[!in_sig])))
  expect_true(all(is.na(r$similarity)))
  expect_identical(
    y@metadata$dedup,
    list(method = "collapse_pathways", similarity = NULL, cutoff = NULL, term_threshold = 0.05)
  )
})

test_that("collapse_pathways needs a ranking for every contrast", {
  fx <- fgsea_example()
  x <- suppressMessages(as_enrichment(list(A = fx$res, B = fx$res)))
  expect_error(
    dedup_terms(x, fx$pathways, method = "collapse_pathways"),
    class = "enrichVolcano_param_error"
  )
  expect_error(
    dedup_terms(x, fx$pathways, method = "collapse_pathways", stats = list(A = fx$ranks)),
    "B",
    class = "enrichVolcano_input_error"
  )
})

test_that("collapse_pathways refuses results that were not ranked GSEA", {
  fry <- utils::read.csv(test_path("fixtures", "limma_fry.csv"), row.names = 1)
  x <- suppressMessages(as_enrichment(list(A = fry)))
  sets <- list(SET_UP = "G1", SET_DOWN = "G11", SET_NULL = "G41")
  expect_error(
    dedup_terms(x, sets, method = "collapse_pathways", stats = list(A = c(G1 = 1))),
    class = "enrichVolcano_input_error"
  )
})

test_that("collapse_pathways keeps terms that have no gene set", {
  fx <- fgsea_example()
  x <- suppressMessages(as_enrichment(list(A = fx$res)))
  top <- x@results$term[which.min(x@results$padj)]
  sets <- fx$pathways[setdiff(names(fx$pathways), top)]
  withr::local_seed(2)
  expect_message(
    y <- dedup_terms(x, sets, method = "collapse_pathways", stats = list(A = fx$ranks)),
    "1 term"
  )
  expect_identical(y@results$dedup_status[y@results$term == top], "kept")
})
