test_that("ev_theme rejects unknown palette names", {
  expect_error(ev_theme(palette = "neon"), class = "ev_unknown_palette")
})

test_that("ev_theme default palette has expected hex codes", {
  th <- ev_theme()
  expect_equal(th$palette$up, "#D6604D")
  expect_equal(th$palette$down, "#4393C3")
  expect_length(th$palette$nes_scale, 5)
  expect_length(th$palette$nes_values, 5)
})
