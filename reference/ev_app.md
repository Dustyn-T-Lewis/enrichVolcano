# Launch the interactive enrichVolcano GUI

Opens a Shiny app to load differential abundance results, choose
databases and thresholds, and render the volcano-in-ring composite
interactively. Requires the `shiny` package (Suggests).

## Usage

``` r
ev_app(...)
```

## Arguments

- ...:

  Passed to
  [`shiny::runApp()`](https://rdrr.io/pkg/shiny/man/runApp.html) (e.g.
  `launch.browser`, `port`).

## Value

Invisibly `NULL`; called for the side effect of running the app.
