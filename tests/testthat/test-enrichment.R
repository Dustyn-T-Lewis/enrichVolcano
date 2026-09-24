toy_enrichment <- function(results = make_toy_results(),
                           metadata = make_toy_metadata()) {
  enrichment(results = results, metadata = metadata)
}

test_that("a valid results table and metadata construct an enrichment", {
  x <- toy_enrichment()
  expect_true(S7::S7_inherits(x, enrichment))
  expect_identical(nrow(x@results), 6L)
  expect_identical(x@metadata$enrichment_test, "fgsea")
})

test_that("a missing required column is named in the error", {
  r <- make_toy_results()
  r$padj <- NULL
  expect_error(toy_enrichment(results = r), "padj")
})

test_that("a wrongly typed column is rejected", {
  r <- make_toy_results()
  r$score <- as.character(r$score)
  expect_error(toy_enrichment(results = r), "score.*numeric")
})

test_that("leading_edge must be a list column", {
  r <- make_toy_results()
  r$leading_edge <- "G1;G2"
  expect_error(toy_enrichment(results = r), "leading_edge.*list")
})

test_that("padj and p outside [0, 1] are rejected, NA is allowed", {
  r <- make_toy_results()
  r$padj[1] <- 1.2
  expect_error(toy_enrichment(results = r), "padj.*\\[0, 1\\]")
  r <- make_toy_results()
  r$p[2] <- -0.1
  expect_error(toy_enrichment(results = r), "`p`.*\\[0, 1\\]")
  r <- make_toy_results()
  r$padj[3] <- NA
  expect_silent(toy_enrichment(results = r))
})

test_that("direction must be up or down and agree with the score sign", {
  r <- make_toy_results()
  r$direction[1] <- "Up"
  expect_error(toy_enrichment(results = r), "direction.*up.*down")
  r <- make_toy_results()
  r$direction[1] <- "down"
  expect_error(toy_enrichment(results = r), "disagree.*score")
})

test_that("a zero or missing score does not trip the sign check", {
  r <- make_toy_results()
  r$score[1] <- 0
  r$score[2] <- NA
  expect_silent(toy_enrichment(results = r))
})

test_that("duplicate contrast / database / term rows are rejected", {
  r <- make_toy_results()
  r <- rbind(r, r[1, ])
  expect_error(toy_enrichment(results = r), "duplicate")
})

test_that("rows with a missing database still count as distinct keys", {
  r <- make_toy_results()
  r <- rbind(r, r[3, ])
  expect_error(toy_enrichment(results = r), "duplicate")
})

test_that("metadata must name a known enrichment test and a score type", {
  expect_error(
    toy_enrichment(metadata = make_toy_metadata(enrichment_test = "gsva")),
    "enrichment_test"
  )
  expect_error(
    toy_enrichment(metadata = make_toy_metadata(score_type = "")),
    "score_type"
  )
})

test_that("editing an object is validated too", {
  x <- toy_enrichment()
  bad <- x@results
  bad$direction <- "sideways"
  expect_error(x@results <- bad, "direction")
})

test_that("print summarises test, score, contrasts, databases and dedup", {
  x <- toy_enrichment()
  out <- cli::cli_fmt(print(x))
  expect_match(out, "fgsea", all = FALSE)
  expect_match(out, "NES", all = FALSE)
  expect_match(out, "A.*B", all = FALSE)
  expect_match(out, "Hallmark", all = FALSE)
  expect_match(out, "unlabelled", all = FALSE)
  expect_match(out, "not deduplicated", all = FALSE)
  x@metadata$dedup <- list(method = "enrichmentmap", similarity = "jaccard", cutoff = 0.5)
  out <- cli::cli_fmt(print(x))
  expect_match(out, "enrichmentmap.*jaccard.*0.5", all = FALSE)
})

test_that("print returns the object invisibly", {
  x <- toy_enrichment()
  expect_invisible(suppressMessages(print(x)))
})
