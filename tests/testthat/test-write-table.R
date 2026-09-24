toy_fry <- function() {
  fry <- utils::read.csv(test_path("fixtures", "limma_fry.csv"), row.names = 1)
  suppressMessages(as_enrichment(list(toy = fry)))
}

test_that("write_table writes one row per term with leading edges joined by ;", {
  x <- make_toy_ring_enrichment()
  file <- withr::local_tempfile(fileext = ".csv")
  expect_invisible(out <- write_table(x, file))
  expect_identical(out, file)
  back <- utils::read.csv(file)
  expect_identical(nrow(back), nrow(x@results))
  expect_identical(back$leading_edge[1], paste(x@results$leading_edge[[1]], collapse = ";"))
  expect_true(all(back$enrichment_test == "custom"))
})

test_that("a written table reads back into the same results", {
  x <- make_toy_ring_enrichment()
  file <- withr::local_tempfile(fileext = ".csv")
  write_table(x, file)
  y <- suppressMessages(as_enrichment(utils::read.csv(file),
    enrichment_test = "custom", leading_edge = "leading_edge", score_type = "NES"
  ))
  core <- names(results_columns)
  expect_equal(y@results[core], x@results[core])
})

test_that("the list from run_enrichment is stacked with each test named", {
  x <- make_toy_ring_enrichment()
  file <- withr::local_tempfile(fileext = ".csv")
  write_table(list(fgsea = x, fry = toy_fry()), file)
  back <- utils::read.csv(file)
  expect_identical(unique(back$enrichment_test), c("custom", "fry"))
  expect_identical(nrow(back), nrow(x@results) + nrow(toy_fry()@results))
})

test_that("write_table refuses anything but enrichment objects", {
  file <- withr::local_tempfile(fileext = ".csv")
  expect_error(write_table(data.frame(), file), class = "enrichVolcano_input_error")
  expect_error(write_table(list(), file), class = "enrichVolcano_input_error")
  expect_error(write_table(list(a = make_toy_ring_enrichment(), b = 1), file), class = "enrichVolcano_input_error")
})
