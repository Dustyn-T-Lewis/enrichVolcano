test_that("volcano_ring_theme overrides point colours without touching the rest", {
  th <- volcano_ring_theme(up = "#B2182B", down = "#2166AC", ns = "grey80")
  expect_equal(th$palette$up, "#B2182B")
  expect_equal(th$palette$down, "#2166AC")
  expect_equal(th$palette$ns, "grey80")
  # untouched palette pieces fall back to the preset
  expect_equal(th$palette$nes_scale, volcano_ring_theme()$palette$nes_scale)
})

test_that("nes_colors replaces the ramp and spreads stops across the limits", {
  th <- volcano_ring_theme(
    nes_colors = c("#053061", "white", "#67001F"),
    nes_limits = c(-2, 2)
  )
  expect_equal(th$palette$nes_scale, c("#053061", "white", "#67001F"))
  expect_equal(th$palette$nes_values, c(-2, 0, 2))
})

test_that("nes_stops must match the (possibly overridden) ramp length", {
  expect_error(
    volcano_ring_theme(
      nes_colors = c("#053061", "#67001F"),
      nes_stops = c(-3, 0, 3)
    ),
    class = "enrichVolcano_param_error"
  )
})

test_that("an invalid colour override is caught at the boundary", {
  expect_error(
    volcano_ring_theme(up = "not_a_real_colour"),
    class = "ev_invalid_colour"
  )
})
