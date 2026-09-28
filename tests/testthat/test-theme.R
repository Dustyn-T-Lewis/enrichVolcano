test_that("plot_theme overrides point colours without touching the rest", {
  th <- plot_theme(up = "#B2182B", down = "#2166AC", ns = "grey80", both = "#1B7837")
  expect_equal(th$palette$up, "#B2182B")
  expect_equal(th$palette$both, "#1B7837")
  expect_equal(th$palette$down, "#2166AC")
  expect_equal(th$palette$ns, "grey80")
  # untouched palette pieces fall back to the preset
  expect_equal(th$palette$nes_scale, plot_theme()$palette$nes_scale)
})

test_that("score_colours replaces the ramp and spreads stops across the limits", {
  th <- plot_theme(
    score_colours = c("#053061", "white", "#67001F"),
    score_limits = c(-2, 2)
  )
  expect_equal(th$palette$nes_scale, c("#053061", "white", "#67001F"))
  expect_equal(th$palette$nes_values, c(-2, 0, 2))
})

test_that("score_stops of matching length overrides the ramp positions", {
  th <- plot_theme(score_stops = c(-3, -1, 0, 1, 3))
  expect_equal(th$palette$nes_values, c(-3, -1, 0, 1, 3))
})

test_that("score_stops must match the (possibly overridden) ramp length", {
  expect_error(
    plot_theme(
      score_colours = c("#053061", "#67001F"),
      score_stops = c(-3, 0, 3)
    ),
    class = "enrichVolcano_param_error"
  )
})

test_that("an invalid colour override is caught at the boundary", {
  expect_error(
    plot_theme(up = "not_a_real_colour"),
    class = "enrichVolcano_param_error"
  )
})

test_that("white label text stays readable on every preset's up and down colours", {
  contrast <- function(col) {
    v <- grDevices::col2rgb(col)[, 1] / 255
    v <- ifelse(v <= 0.04045, v / 12.92, ((v + 0.055) / 1.055)^2.4)
    1.05 / (sum(c(0.2126, 0.7152, 0.0722) * v) + 0.05)
  }
  for (palette in c("default", "viridis", "okabe")) {
    pal <- plot_theme(palette = palette)$palette
    expect_gte(contrast(pal$up), 3, label = paste(palette, "up"))
    expect_gte(contrast(pal$down), 3, label = paste(palette, "down"))
  }
})
