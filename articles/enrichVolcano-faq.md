# FAQ

``` r

library(enrichVolcano)
```

## Why is my ring plot empty?

Four likely causes.

The classic one is the **wrong ranking statistic**. GSEA needs a
*signed* score so up- and down-regulated genes sit at opposite ends of
the ranked list (Subramanian et al. 2005; Reimand et al. 2019). Ranking
by an unsigned pi-score collapses everything to one tail and yields
near-empty rings. The package defaults to `rank_by = "signed_p"`
(`sign(logFC) * -log10(P)`) and now refuses an unsigned pi-score with an
`ev_unsigned_ranking` error rather than silently mis-ranking. If you see
that error, drop the `rank_by` override and use the default. See
[`vignette("scoring")`](https://Dustyn-T-Lewis.github.io/enrichVolcano/articles/scoring.md)
for the full explanation.

Second, no pathways passed the `enrich_padj` filter. The default is
0.05, which is strict on a small contrast or a small fixture. Inspect
the padj distribution on the figure object:

``` r

p <- enrich_volcano(data, contrast = "ctr")
summary(attr(p, "ev_data")$enrichment$padj)
```

If the smallest padj is above 0.05, raise the cutoff:

``` r

p <- enrich_volcano(data, contrast = "ctr", enrich_padj = 0.25)
```

Third, your gene symbols do not match the database species. msigdbr
collections key on human symbols for `species = "human"`, mouse symbols
for `"mouse"`, and so on. If your input has rat symbols but you passed
`species = "human"`, the overlap is empty by construction.

Fourth, all pathways were dropped by
[`ev_collapse()`](https://Dustyn-T-Lewis.github.io/enrichVolcano/reference/ev_collapse.md).
Loosen the dedup cutoff by passing a higher value through the `dedup`
list:

``` r

p <- enrich_volcano(
  data, contrast = "ctr",
  dedup = list(method = "jaccard", cutoff = 0.9, scope = "within_db")
)
```

A cutoff close to 1 keeps almost everything. If the ring still fills,
the original cutoff was eating real biology; if it stays empty, the
issue is upstream of dedup.

## Can I use DESeq2 / edgeR / Perseus output directly?

Yes.
[`ev_validate()`](https://Dustyn-T-Lewis.github.io/enrichVolcano/reference/ev_validate.md)
auto-detects input formats from column patterns. Supported sources: tidy
long, wide, DEP `data.frame`, proteoDA result tables, limma `topTable`,
DESeq2 `results`, edgeR `topTags`, MSstats `comparisonResult`, proDA
`test_diff`, DEqMS, MaxQuant ratio outputs, and Perseus exported tables.
Pass the data frame directly to
[`enrich_volcano()`](https://Dustyn-T-Lewis.github.io/enrichVolcano/reference/enrich_volcano.md).

Perseus reports its p-value as `-Log10 p-value`.
[`ev_validate()`](https://Dustyn-T-Lewis.github.io/enrichVolcano/reference/ev_validate.md)
back-transforms with a `cli_inform` notice so you see what happened.
MaxQuant ratio outputs lack an FDR column; BH is synthesized at
validation time, also with a notice.

DEP `SummarizedExperiment` objects and proteoDA `DAList` objects are not
yet directly supported in v0.1. Extract `rowData()` (DEP) or `$results`
(proteoDA) to a data frame first:

``` r

df <- as.data.frame(SummarizedExperiment::rowData(dep_se))
p  <- enrich_volcano(df, contrast = "ctr")
```

## How do I force-label specific genes on the volcano?

In v0.1,
[`enrich_volcano()`](https://Dustyn-T-Lewis.github.io/enrichVolcano/reference/enrich_volcano.md)
does not forward `label_genes` through the `volcano` list. Call
[`ev_volcano()`](https://Dustyn-T-Lewis.github.io/enrichVolcano/reference/ev_volcano.md)
directly:

``` r

p_volcano <- ev_volcano(
  data, contrast = "ctr",
  label_n     = 10,
  label_by    = "pi_eq2",
  label_genes = c("GENE1", "GENE2")
)
```

Then build the ring and compose:

``` r

enr <- ev_enrich(data, contrast = "ctr", databases = "hallmark")
p_ring <- ring_plot(enr, contrast = "ctr")
p <- ev_compose(list(ctr = p_volcano), list(ctr = p_ring))
```

`label_genes` forwarding through
[`enrich_volcano()`](https://Dustyn-T-Lewis.github.io/enrichVolcano/reference/enrich_volcano.md)
is on the v0.2 list.

## Why don’t I see KEGG pathways even though I requested them?

Four checks.

- `msigdbr` is installed and reasonably current. KEGG_MEDICUS was added
  to MSigDB in 2024; older `msigdbr` releases ship `kegg_legacy` only.
- The species matches. KEGG via msigdbr supports the seven species in
  the registry. Pass `species = "human"` (or the matching value)
  alongside `databases = "kegg_medicus"`.
- `database_info("kegg_medicus")` confirms the database is registered
  and shows the source. If the source reads `msigdbr`, the lookup goes
  through
  [`msigdbr::msigdbr()`](https://igordot.github.io/msigdbr/reference/msigdbr.html).
- The default `exclude_terms = "DISEASE|CANCER|TUMOR"` regex inside
  [`ev_enrich()`](https://Dustyn-T-Lewis.github.io/enrichVolcano/reference/ev_enrich.md)
  strips disease-named KEGG entries out of pathway databases. KEGG
  includes hundreds of disease-named entries; the default filters them
  so they do not double-count against `disgenet` results. To keep them,
  call
  [`ev_enrich()`](https://Dustyn-T-Lewis.github.io/enrichVolcano/reference/ev_enrich.md)
  directly:

``` r

enr <- ev_enrich(
  data, contrast = "ctr",
  databases     = "kegg_medicus",
  exclude_terms = NULL
)
p_ring <- ring_plot(enr, contrast = "ctr")
```

## How do I add my own gene set file?

For v0.1, pass a named list of named lists. The outer name is the
database name; the inner names are pathway names; the values are
character vectors of gene symbols.

``` r

my_paths <- list(
  CUSTOM_SET_1 = c("GENE1", "GENE2", "GENE3"),
  CUSTOM_SET_2 = c("GENE4", "GENE5")
)
p <- enrich_volcano(
  data      = data,
  contrast  = "ctr",
  databases = list(custom = my_paths)
)
```

To load from a GMT file, read it with
[`fgsea::gmtPathways()`](https://rdrr.io/pkg/fgsea/man/gmtPathways.html)
and pass the result through the same shape:

``` r

my_paths <- fgsea::gmtPathways("my_sets.gmt")
p <- enrich_volcano(
  data      = data,
  contrast  = "ctr",
  databases = list(custom = my_paths)
)
```

A direct `custom_gmt = "path/to/sets.gmt"` argument is planned for v0.2.

## My data is rat (or mouse / fly / zebrafish / yeast / pig). What do I pass?

Set `species` and the registered databases follow:

``` r

p <- enrich_volcano(
  data      = data,
  contrast  = "ctr",
  species   = "rat",
  databases = "mitocarta3"
)
```

Supported species strings: `"human"`, `"mouse"`, `"rat"`, `"zebrafish"`,
`"fly"`, `"yeast"`, `"pig"`. msigdbr handles ortholog mapping for
fetched collections; bundled GMTs (MitoCarta, HPA, UniProt) ship
per-species files where the species is supported.

Input gene symbols must match the species namespace. Rat symbols for a
rat study, mouse symbols for a mouse study. Mixing namespaces silently
returns empty enrichment.

## Why is the rendered PDF different from the screen preview?

Three sizing rules to remember.

`ggsave()` defaults to inches. Pass `units = "mm"` (or `"cm"`)
explicitly when sizing for a journal column width specified in
millimetres.

Patchwork inherits panel sizes from its components. If the volcano and
ring panels carry very different intrinsic aspect ratios, the on-screen
preview at one device size will not match the saved file at a different
size. Set fixed sizes with `&` and `theme(aspect.ratio = 1)` if the
layout matters.

The ring uses `coord_fixed()`. Its 1:1 aspect ratio is locked at draw
time and ignores save dimensions. The ring will always be a circle; the
patchwork resizes the volcano to fill the remaining space.

For predictable output, save to PDF first at the exact width and height
you need, then check the PDF in a viewer rather than the RStudio plots
pane.

## How do I cite this package?

``` r

citation("enrichVolcano")
```

A full citation needs four references: this package, Xiao et al. 2014
for the pi-score, Korotkevich et al. 2021 for fgsea, and Benjamini and
Hochberg 1995 for BH adjustment. All four are in `inst/REFERENCES.bib`
and ship with the package install.

If you used a specific database, cite its primary source as well.
`database_info(name)$citation_key` resolves to the BibTeX key in
`REFERENCES.bib`.
