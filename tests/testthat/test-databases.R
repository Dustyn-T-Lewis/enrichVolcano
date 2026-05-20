test_that("list_databases returns tibble with required schema", {
  db <- list_databases()
  expect_s3_class(db, "tbl_df")
  expect_true(all(c("name", "type", "species", "license") %in% colnames(db)))
  expect_true(nrow(db) >= 17)
  expect_true(all(nchar(db$license) > 0))
})

test_that("database_info prints details for a registered DB", {
  expect_message(database_info("hallmark"), regexp = "MSigDB Hallmark")
  expect_error(database_info("nope"), class = "ev_unknown_database")
})

test_that("bundled GMT files are valid when present", {
  for (db in c("mitocarta3")) {
    for (sp in c("human", "mouse", "rat")) {
      fn <- paste0(db, "_", sp, ".gmt")
      path <- system.file("extdata", "gmt", fn, package = "enrichVolcano")
      if (!nzchar(path)) {
        skip(paste("GMT not bundled yet:", fn))
      }
      paths <- fgsea::gmtPathways(path)
      expect_true(length(paths) > 0)
    }
  }
})
