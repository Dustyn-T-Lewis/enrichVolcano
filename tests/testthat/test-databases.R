test_that("build_registry.R reproduces the shipped db_registry.csv", {
  script <- testthat::test_path("..", "..", "data-raw", "build_registry.R")
  skip_if(!file.exists(script),
          "data-raw/build_registry.R not available (installed package)")
  shipped <- system.file("extdata", "db_registry.csv",
                         package = "enrichVolcano")
  script <- normalizePath(script)
  tmp <- withr::local_tempdir()
  withr::with_dir(tmp, sys.source(script, envir = new.env()))
  built <- file.path(tmp, "inst", "extdata", "db_registry.csv")
  expect_true(file.exists(built))
  expect_identical(readLines(built), readLines(shipped))
})

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
  for (db in c("mitocarta3", "corum")) {
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

test_that("go_slim is registered and its bundled GMTs are valid", {
  reg <- list_databases()
  expect_true("go_slim" %in% reg$name)
  for (sp in c("human", "mouse", "rat")) {
    f <- system.file("extdata", "gmt", paste0("go_slim_", sp, ".gmt"),
                     package = "enrichVolcano")
    skip_if(!nzchar(f), paste("go_slim GMT not bundled:", sp))
    p <- fgsea::gmtPathways(f)
    expect_gt(length(p), 30)
    expect_true(all(grepl("^GOSLIM_", names(p))))
  }
})
