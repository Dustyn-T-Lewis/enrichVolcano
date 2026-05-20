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
