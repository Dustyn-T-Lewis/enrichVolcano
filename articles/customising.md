# Customising enrichVolcano plots

``` r

library(enrichVolcano)
```

## Themes and palettes

[`ev_theme()`](https://Dustyn-T-Lewis.github.io/enrichVolcano/reference/ev_theme.md)
returns a list with three elements: `base_size`, `theme` (a ggplot2
theme object), and `palette` (named hex codes).

``` r

th <- ev_theme()
names(th)
#> [1] "base_size" "theme"     "palette"
th$palette
#> $up
#> [1] "#D6604D"
#> 
#> $down
#> [1] "#4393C3"
#> 
#> $ns
#> [1] "grey70"
#> 
#> $nes_scale
#> [1] "#08306B" "white"   "#67000D"
```

The default palette is locked. `up = "#D6604D"` is the rust-red used
across the YvO and CvH figures. `down = "#4393C3"` is the complementary
blue. `ns = "grey70"` is the non-significant point colour.
`nes_scale = c("#08306B", "white", "#67000D")` is the diverging colour
ramp used for NES on the ring panel.

For poster fonts, pass a larger base size:

``` r

th <- ev_theme(base_size = 14)
```

The `theme` element drops into any ggplot via `+ th$theme`. To change
colours, edit the returned palette element and pass the whole list to
the volcano or ring functions via their `theme` argument:

``` r

th <- ev_theme()
th$palette$up   <- "#B2182B"
th$palette$down <- "#2166AC"

p_volcano <- ev_volcano(yvo, contrast = "Aging", theme = th)
```

`ev_theme()$palette` is a plain named list, so any standard list edit
works. Both
[`ev_volcano()`](https://Dustyn-T-Lewis.github.io/enrichVolcano/reference/ev_volcano.md)
and
[`ring_plot()`](https://Dustyn-T-Lewis.github.io/enrichVolcano/reference/ring_plot.md)
read the palette out of `theme$palette` at draw time.

## Volcano label control

The volcano panel honours two label arguments passed through the
`volcano` list to
[`enrich_volcano()`](https://Dustyn-T-Lewis.github.io/enrichVolcano/reference/enrich_volcano.md).

``` r

p <- enrich_volcano(
  data     = yvo,
  contrast = "Aging",
  volcano  = list(
    label_n  = 10,
    label_by = "pi_eq2"
  )
)
```

`label_n` is the number of top genes to label. `label_by` is the column
used to rank them: `"pi_eq2"`, `"pi_eq1"`, `"adj.P.Val"`, or `"logFC"`.
Smaller is better for pi-score and p-value columns; absolute value is
used for `logFC`.

To force-label a specific set of genes regardless of rank, call
[`ev_volcano()`](https://Dustyn-T-Lewis.github.io/enrichVolcano/reference/ev_volcano.md)
directly with `label_genes`:

``` r

p_volcano <- ev_volcano(
  yvo,
  contrast    = "Aging",
  label_n     = 10,
  label_by    = "pi_eq2",
  label_genes = c("MYH7", "ATP5F1A", "PPARGC1A")
)
```

[`enrich_volcano()`](https://Dustyn-T-Lewis.github.io/enrichVolcano/reference/enrich_volcano.md)
in v0.1 forwards `label_n` and `label_by` to
[`ev_volcano()`](https://Dustyn-T-Lewis.github.io/enrichVolcano/reference/ev_volcano.md)
but does not pass `label_genes` through the `volcano` list. Forwarding
is planned for v0.2. Until then, build the volcano panel with
[`ev_volcano()`](https://Dustyn-T-Lewis.github.io/enrichVolcano/reference/ev_volcano.md),
the ring with
[`ring_plot()`](https://Dustyn-T-Lewis.github.io/enrichVolcano/reference/ring_plot.md),
and stitch with
[`ev_compose()`](https://Dustyn-T-Lewis.github.io/enrichVolcano/reference/ev_compose.md).

`ggrepel` parameters are set inside
[`ev_volcano()`](https://Dustyn-T-Lewis.github.io/enrichVolcano/reference/ev_volcano.md):
`max.overlaps = 50`, `force = 6`, `box.padding = 0.3`. These match the
YvO 2025 figure series. They are not user-facing in v0.1; if you need
different repel behaviour, call
[`ev_volcano()`](https://Dustyn-T-Lewis.github.io/enrichVolcano/reference/ev_volcano.md)
directly and add your own `geom_text_repel()` layer:

``` r

p <- ev_volcano(yvo, contrast = "Aging")
p <- p + ggrepel::geom_text_repel(
  data = subset(yvo, gene %in% c("MYH7", "ATP5F1A")),
  aes(label = gene),
  max.overlaps = Inf
)
```

The volcano returned by
[`ev_volcano()`](https://Dustyn-T-Lewis.github.io/enrichVolcano/reference/ev_volcano.md)
is a plain ggplot, so any ggplot2 modifier works.

## Ring encodings

The ring panel maps three pathway-level columns to three visual
channels. All three are exposed through the `ring` list.

``` r

p <- enrich_volcano(
  data     = yvo,
  contrast = "Aging",
  ring     = list(
    max_terms = 8,
    order_by  = "padj",
    order_dir = "asc",
    magnitude = "neg_log_padj",
    color     = "nes"
  )
)
```

`magnitude` sets arc thickness. Options are `"neg_log_padj"` (default;
thicker = more significant), `"nes"` (thicker = stronger enrichment
regardless of direction), and `"gene_count"` (thicker = larger pathway).

`color` sets arc fill. Options are `"nes"` (default; diverging blue to
white to red across the NES range), `"neg_log_padj"` (sequential ramp),
and `"direction"` (categorical: up red, down blue).

`order_by` controls the arc layout around the ring. Options are `"padj"`
(default), `"size"`, `"NES"`, and `"alphabetic"`. Combine with
`order_dir = "asc"` or `"desc"` to flip direction.

A practical example: map arc thickness to gene count and arc fill to
NES, ordered by absolute NES:

``` r

p <- enrich_volcano(
  data     = yvo,
  contrast = "Aging",
  ring     = list(
    max_terms = 8,
    order_by  = "NES",
    order_dir = "desc",
    magnitude = "gene_count",
    color     = "nes"
  )
)
```

This layout shows the biggest pathways first and uses fill for
direction, which works well when reviewers want to see effect size at a
glance.

## Multi-contrast layouts

Pass a vector to `contrast` to build one volcano-ring pair per contrast.
[`ev_compose()`](https://Dustyn-T-Lewis.github.io/enrichVolcano/reference/ev_compose.md)
arranges the pairs into a single patchwork.

``` r

p <- enrich_volcano(
  data     = yvo,
  contrast = c("Aging", "Training_Young", "Training_Old",
               "Interaction"),
  databases = list(mock = mock_msigdb),
  facet    = list(nrow = 2, ncol = 2)
)
```

`facet$nrow` and `facet$ncol` go to
[`patchwork::plot_layout()`](https://patchwork.data-imaginist.com/reference/plot_layout.html).
The result is a 2x2 grid of volcano-ring pairs. Each pair keeps its own
contrast title; the outer layout collects guides into a single legend
strip via `plot_layout(guides = "collect")`.

The 4-contrast YvO example produces the same composite that appears in
F02 of the YvO 2025 pipeline. The fixture has the four contrasts needed;
substitute a real database for `mock_msigdb` to reproduce the published
figure.

## Term-label cleaning

Ring labels can be long. `HALLMARK_OXIDATIVE_PHOSPHORYLATION` is 35
characters and crowds the ring. `ev_default_clean_label()` strips the
database prefix (`HALLMARK_`, `REACTOME_`, `KEGG_LEGACY_`,
`KEGG_MEDICUS_`) and replaces underscores with spaces, leaving
`Oxidative Phosphorylation`.

The default fires automatically inside
[`ring_plot()`](https://Dustyn-T-Lewis.github.io/enrichVolcano/reference/ring_plot.md).
For v0.1, the clean function is fixed. A custom override
(`ring$label_clean_fn = function(x) gsub(...)`) is planned for v0.2.
Until then, post-process the ggplot returned by
[`ring_plot()`](https://Dustyn-T-Lewis.github.io/enrichVolcano/reference/ring_plot.md):

``` r

p <- ring_plot(enr, contrast = "Aging")
p$data$label <- gsub("Oxidative", "OxPhos", p$data$label)
```

This works because the label data lives on the ggplot, not in a separate
text annotation.

## Exporting for journals

[`enrich_volcano()`](https://Dustyn-T-Lewis.github.io/enrichVolcano/reference/enrich_volcano.md)
returns a patchwork. `ggsave()` works directly:

``` r

ggplot2::ggsave("fig.pdf", p, width = 180, height = 120, units = "mm")
```

A sizing gotcha: nested patchwork annotations need
`plot_layout(guides = "collect")` at the outer level to share legends.
[`ev_compose()`](https://Dustyn-T-Lewis.github.io/enrichVolcano/reference/ev_compose.md)
already sets this, so multi-contrast figures collect guides
automatically. If you wrap an
[`enrich_volcano()`](https://Dustyn-T-Lewis.github.io/enrichVolcano/reference/enrich_volcano.md)
result inside a larger patchwork, set `guides = "collect"` again at the
new outer level or the legends will multiply.

Common sizes for journal submission:

- Single column: `width = 89, height = 80, units = "mm"`.
- One and a half column: `width = 120, height = 90, units = "mm"`.
- Double column: `width = 180, height = 120, units = "mm"`.

For PNG at 300 dpi:

``` r

ggplot2::ggsave("fig.png", p, width = 180, height = 120,
                units = "mm", dpi = 300)
```

For TIFF with LZW compression (some journals require it):

``` r

ggplot2::ggsave("fig.tiff", p, width = 180, height = 120,
                units = "mm", dpi = 300, compression = "lzw")
```

If you hit issues, check
[`vignette("enrichVolcano-faq")`](https://Dustyn-T-Lewis.github.io/enrichVolcano/articles/enrichVolcano-faq.md).
