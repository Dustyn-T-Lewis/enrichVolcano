# Write panels to one PDF

Page 1 lays every panel out as a composite lettered A, B, C in list
order. Each later page holds one panel at full size, headed by its
letter and name.
[`grDevices::cairo_pdf()`](https://rdrr.io/r/grDevices/cairo.html)
embeds the fonts, so the file looks the same on any machine. On macOS,
R's cairo needs XQuartz (<https://www.xquartz.org>).

## Usage

``` r
write_plot(
  panels,
  file,
  ncol = NULL,
  nrow = NULL,
  design = NULL,
  caption = NULL,
  width = 11,
  height = 8.5
)
```

## Arguments

- panels:

  A named list of ggplots, such as the output of
  [`plot_volcano_ring()`](https://Dustyn-T-Lewis.github.io/enrichVolcano/reference/plot_volcano_ring.md)
  and
  [`plot_scatter()`](https://Dustyn-T-Lewis.github.io/enrichVolcano/reference/plot_scatter.md).
  Names head the single-panel pages.

- file:

  Path of the PDF to write.

- ncol, nrow:

  Grid of the composite, passed to
  [`patchwork::wrap_plots()`](https://patchwork.data-imaginist.com/reference/wrap_plots.html).

- design:

  A patchwork layout such as `"AB\nCC"`, used instead of `ncol` and
  `nrow`.

- caption:

  Text under the composite. `NULL` draws none.

- width, height:

  Page size in inches. The default is US Letter landscape.

## Value

`file`, invisibly.

## Examples

``` r
da <- read_example("yvo")$da
#> 2103 of 2106 accessions mapped to Homo sapiens symbols.
#> 2106 of 2106 matrix proteins have DA results.
ex <- as_enrichment(read.csv(system.file("extdata", "examples", "yvo_fgsea.csv.gz",
  package = "enrichVolcano"
)))
#> Reading "fgsea" results.
panels <- list(
  "Young, trained" = plot_volcano_ring(da, ex, contrast = "Training_Young"),
  "Old, trained" = plot_volcano_ring(da, ex, contrast = "Training_Old")
)
file <- tempfile(fileext = ".pdf")
try(write_plot(panels, file, ncol = 2))
```
