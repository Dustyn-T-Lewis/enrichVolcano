test_that("ev_app launcher exists and the bundled app is present", {
  expect_true(is.function(ev_app))
  expect_true(nzchar(system.file("shiny", "app.R", package = "enrichVolcano")))
})

test_that("ev_app errors clearly when shiny is unavailable", {
  skip_if(requireNamespace("shiny", quietly = TRUE),
          "shiny installed; cannot test missing-shiny path")
  expect_error(ev_app(), class = "ev_shiny_missing")
})

test_that("bundled shiny app is constructible", {
  skip_if_not_installed("shiny")
  app <- shiny::shinyAppDir(system.file("shiny", package = "enrichVolcano"))
  expect_s3_class(app, "shiny.appobj")
})

test_that("ev_app server validates the bundled example and lists its contrasts", {
  skip_if_not_installed("shiny")
  app_dir <- system.file("shiny", package = "enrichVolcano")
  skip_if(!nzchar(app_dir), "bundled app not installed")
  shiny::testServer(app_dir, {
    session$setInputs(source = "example")
    v <- validated()
    expect_s3_class(v, "tbl_df")
    expect_true(!is.null(attr(v, "ev_source")))
    expect_gt(length(unique(v$contrast)), 0)
    expect_match(output$summary, "Detected format")
  })
})
