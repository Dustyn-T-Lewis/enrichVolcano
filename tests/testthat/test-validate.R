source(test_path("fixtures/make_toy.R"))

test_that("validate_ring_geometry rejects a bad label_headroom", {
  for (bad in list(-1, c(0.5, 1), "0.5", NA_real_)) {
    expect_error(
      validate_ring_geometry(4.8, 4.0, c(0.4, 1.6), label_headroom = bad),
      class = "enrichVolcano_param_error"
    )
  }
})
