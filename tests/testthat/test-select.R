source(test_path("fixtures/make_toy.R"))

toy <- function() enrichment(results = make_toy_results(), metadata = make_toy_metadata())

test_that("databases keeps only the named collections", {
  r <- filter_view(toy()@results, databases = "Hallmark", collapse = TRUE)
  expect_true(all(r$database == "Hallmark"))
  expect_identical(nrow(r), 4L)
})

test_that("databases = NULL keeps every collection", {
  expect_identical(nrow(filter_view(toy()@results, databases = NULL, collapse = TRUE)), 6L)
})

test_that("databases absent from the object are named with what is there", {
  expect_error(
    filter_view(toy()@results, databases = c("KEGG", "Reactome"), collapse = TRUE),
    "Hallmark",
    class = "enrichVolcano_input_error"
  )
})

test_that("an object with no database labels skips the filter and says so", {
  r <- toy()@results
  r$database <- NA_character_
  expect_message(
    out <- filter_view(r, databases = "Hallmark", collapse = TRUE),
    class = "enrichVolcano_databases_skipped"
  )
  expect_identical(nrow(out), 6L)
})

test_that("collapse hides redundant terms and keeps unflagged ones", {
  r <- toy()@results
  r$dedup_status <- c("kept", "redundant", NA, "kept", "kept", NA)
  expect_identical(nrow(filter_view(r, databases = NULL, collapse = TRUE)), 5L)
  expect_identical(nrow(filter_view(r, databases = NULL, collapse = FALSE)), 6L)
})

test_that("ring_terms keeps significant terms, top n per direction by padj", {
  r <- data.frame(
    term = paste0("T", 1:6),
    score = c(2, 1.5, 1, -2, -1.5, -1),
    padj = c(0.03, 0.001, 0.2, 0.01, 0.04, 0.002),
    direction = c("up", "up", "up", "down", "down", "down")
  )
  expect_setequal(
    ring_terms(r, term_threshold = 0.05, n_terms = Inf, terms = NULL)$term,
    c("T1", "T2", "T4", "T5", "T6")
  )
  expect_setequal(
    ring_terms(r, term_threshold = 0.05, n_terms = 1, terms = NULL)$term,
    c("T2", "T6")
  )
})

test_that("ring_terms takes a hand-picked set as given", {
  r <- data.frame(
    term = paste0("T", 1:3), score = c(2, -1, 1),
    padj = c(0.9, 0.001, 0.5), direction = c("up", "down", "up")
  )
  expect_setequal(ring_terms(r, 0.05, 8, terms = c("T1", "T3"))$term, c("T1", "T3"))
  expect_error(ring_terms(r, 0.05, 8, terms = c("T1", "T9")), "T9", class = "enrichVolcano_input_error")
})

test_that("arc magnitude follows the score type unless set", {
  expect_identical(default_magnitude(NULL, "NES"), "neg_log_padj")
  expect_identical(default_magnitude(NULL, "signed -log10(FDR)"), "size")
  expect_identical(default_magnitude("size", "NES"), "size")
  expect_error(default_magnitude("height", "NES"), class = "rlang_error")
})

test_that("NES keeps fixed fill limits; other scores scale to the object", {
  r <- data.frame(score = c(5.2, -8.1, 0.4), padj = c(0.001, 1e-8, 0.6))
  expect_identical(default_score_limits("NES", r, term_threshold = 0.05), c(-3, 3))
  expect_identical(default_score_limits("signed -log10(FDR)", r, term_threshold = 0.05), c(-8.1, 8.1))
})

test_that("score limits fall back to every row, then to 1, when nothing is significant", {
  r <- data.frame(score = c(0.4, -0.2), padj = c(0.6, 0.9))
  expect_identical(default_score_limits("z", r, term_threshold = 0.05), c(-0.4, 0.4))
  r$score <- c(0, 0)
  expect_identical(default_score_limits("z", r, term_threshold = 0.05), c(-1, 1))
})
