test_that("ev_theme returns a list with palette + theme", {
  th <- ev_theme()
  expect_named(th, c("base_size", "theme", "palette"))
  expect_s3_class(th$theme, "theme")
  expect_named(th$palette, c("up", "down", "ns", "nes_scale", "nes_values"))
})

test_that("ev_theme palette colors are locked hex codes", {
  th <- ev_theme()
  expect_equal(th$palette$up, "#D6604D")
  expect_equal(th$palette$down, "#4393C3")
  expect_equal(th$palette$ns, "grey70")
})
