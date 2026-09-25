toy_panels <- function() {
  list(
    Ring = suppressMessages(plot_volcano_ring(make_toy_da(), make_toy_ring_enrichment(),
      databases = NULL, term_threshold = 1, n_terms = Inf
    )),
    Scatter = suppressMessages(plot_scatter(make_toy_scatter_enrichment(), "A", "B"))
  )
}

cairo_works <- function() {
  before <- grDevices::dev.cur()
  suppressWarnings(grDevices::cairo_pdf(withr::local_tempfile(fileext = ".pdf")))
  opened <- !identical(grDevices::dev.cur(), before)
  if (opened) grDevices::dev.off()
  opened
}

pdf_pages <- function(file) pdftools::pdf_info(file)$pages

test_that("write_plot writes the composite, then one page per panel", {
  skip_if_not(cairo_works(), "cairo cannot start here")
  skip_if_not_installed("pdftools")
  file <- withr::local_tempfile(fileext = ".pdf")
  expect_invisible(out <- write_plot(toy_panels(), file))
  expect_identical(out, file)
  expect_identical(pdf_pages(file), 3L)
  write_plot(toy_panels()["Ring"], file)
  expect_identical(pdf_pages(file), 2L)
})

test_that("without cairo, write_plot stops and says how to get it", {
  skip_if(cairo_works(), "cairo works here")
  file <- withr::local_tempfile(fileext = ".pdf")
  expect_error(write_plot(toy_panels(), file), "XQuartz", class = "enrichVolcano_device_error")
  expect_false(file.exists(file))
})

test_that("the composite letters panels in list order, with an optional caption", {
  composite <- compose_panels(toy_panels(), ncol = 2)
  expect_s3_class(composite, "patchwork")
  expect_identical(composite$patches$annotation$tag_levels, "A")
  expect_null(composite$patches$annotation$caption)
  expect_identical(compose_panels(toy_panels(), caption = "n = 20")$patches$annotation$caption, "n = 20")
  expect_s3_class(compose_panels(toy_panels(), design = "AB"), "patchwork")
})

test_that("the page gives each composite cell room for a full-size panel, in landscape", {
  expect_identical(page_size(3, ncol = 3), c(21, 7.5))
  expect_equal(page_size(3, design = "AB\nCC"), c(15 * 11 / 8.5, 15))
  expect_equal(page_size(4), c(15 * 11 / 8.5, 15))
  expect_identical(page_size(2, nrow = 1), c(14, 7.5))
  size <- page_size(1)
  expect_equal(size[1] / size[2], 11 / 8.5)
})

test_that("write_plot sizes the pages from the layout unless told otherwise", {
  skip_if_not(cairo_works(), "cairo cannot start here")
  skip_if_not_installed("pdftools")
  file <- withr::local_tempfile(fileext = ".pdf")
  write_plot(toy_panels(), file, ncol = 2)
  expect_equal(unname(unlist(pdftools::pdf_pagesize(file)[1, c("width", "height")])), c(14, 7.5) * 72, tolerance = 0.01)
  write_plot(toy_panels(), file, width = 11, height = 8.5)
  expect_equal(unname(unlist(pdftools::pdf_pagesize(file)[1, c("width", "height")])), c(11, 8.5) * 72, tolerance = 0.01)
})

test_that("panels must be a named list of ggplots", {
  p <- toy_panels()
  for (bad in list(list(), unname(p), c(p, list(Table = data.frame(a = 1))), p$Ring)) {
    expect_error(compose_panels(bad), class = "enrichVolcano_param_error")
  }
})

test_that("two-panel composite snapshot is stable", {
  skip_on_ci()
  skip_if_not_installed("vdiffr")
  vdiffr::expect_doppelganger("composite-two-panels", compose_panels(toy_panels(), ncol = 2))
})
