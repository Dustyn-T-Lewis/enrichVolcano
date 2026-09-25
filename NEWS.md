# enrichVolcano 2.2.0

## New features

* Hallmark and GO slim terms take their labels from a table written for the
  package, in two styles: `labels = "short"` (the default) fits two ring
  lines, and `labels = "clean"` spells every word out.
* `list_labels()` lists the terms a figure will draw, with both labels, and
  `write_labels()` saves the list to CSV. Edit its `label` column and pass the
  table back through `labels`.
* `labels` also takes a function that names every term.

## Changes

* Terms outside the table lose their prefix and underscores and are
  title-cased, with common acronyms restored. The phrase shortening and
  hand-set line breaks of 2.1.0 are gone, so some Reactome, KEGG and GO:BP
  labels are longer. Shorten them with `write_labels()` or a function.
* MitoCarta names keep only their last level.

# enrichVolcano 2.1.0

## Breaking changes

* `load_gene_sets()` no longer takes `min_size` or `max_size` and keeps every
  set. Sets were filtered twice: on genome size here, then on the genes in
  the data in `run_enrichment()`. The genome filter removed about half the GO
  slim terms, such as signalling and programmed cell death, that proteomics
  data measure well below 500 genes. GO slim results now test more terms.

## New features

* `plot_volcano_ring()` and `plot_scatter()` take `labels`, a character
  vector named by term, to replace the automatic names of chosen terms.

## Bug fixes

* fgsea seeds each contrast and collection separately from your
  `set.seed()`. A collection's results no longer shift when other
  collections run alongside it. fgsea p-values differ from 2.0.0 by
  permutation noise.
* A ring whose terms all go one way keeps them in its own half: down on the
  left, up on the right.
* The ring panel widens to fit labels at its sides, so they are no longer
  clipped or drawn under the legend.
* An empty size window in `run_enrichment()` fails before the tied-rank note.
* A scatter correlation that rounds to zero prints as 0.00, not -0.00.
* `plot_scatter()` draws contrasts with no significant terms instead of
  failing, and `clean_label(character(0))` returns `character(0)`.
* `clean_label()` keeps every word of a MitoCarta name that has no hierarchy,
  and writes complex numerals upper case: `Complex IV`, not `Complex Iv`.
* `plot_bias_ring()` refuses an unnamed extra argument instead of ignoring it.
* `write_plot()` sizes its pages from the layout, about 7 by 7.5 inches per
  composite cell, so a composite of full-size panels stays readable. `width`
  and `height` still set the page.
* A protein group such as `P31040;Q9UBK2` is mapped to its gene symbol
  through its first accession.
* A blank gene symbol counts as missing.
* The annotation packages no longer print an empty line when loaded.
* A wide table without a contrast column gets the wide-table message.
* The `mito` example study ships the Mito pipeline's gene symbols. The rat
  annotation package lacked 447 of its 4,806 proteins, which fgsea dropped.

## Documentation

* The vignette cites the paper behind each method and ends with a reference
  list. `citation("enrichVolcano")` gives the package citation. The README
  shows the workflow as a figure.

# enrichVolcano 2.0.0

2.0.0 renames the exports to one rule, verb then noun, and gives arguments
that mean the same thing one name. Old names are gone, with no aliases.
To stay on the old API, install 1.1.0:
`remotes::install_github("Dustyn-T-Lewis/enrichVolcano@v1.1.0")`.

## Migration

| 1.1.0 | 2.0.0 |
|---|---|
| `example_study(name)` | `read_example(name)` |
| `dedup(x, gene_sets)` | `dedup_terms(enrichment, gene_sets)` |
| `dedup(p_threshold = )` | `dedup_terms(term_threshold = )` |
| `volcano_ring(volc_df, enrichment)` | `plot_volcano_ring(da, enrichment)` |
| `volcano_ring(gene_col, logfc_col, pval_col, padj_col)` | removed; pass `as_da()` output |
| `volcano_ring(volc_sig_col = "pi")` | `as_da(x, padj = "pi")`, then `plot_volcano_ring()` |
| `volcano_ring(disc_color = )` | `plot_volcano_ring(disc_colour = )` |
| `volcano_ring(nes_limits = )` | `plot_volcano_ring(score_limits = )` |
| `volcano_ring(arc_order = "nes")` | `plot_volcano_ring(arc_order = "score")` |
| `volcano_ring_grid(volc_dfs, enrichment, contrasts)` | `lapply()` over contrasts, then `write_plot()` or `patchwork::wrap_plots()` |
| `nes_scatter(enrichment, x, y)` | `plot_scatter(enrichment, x, y)` |
| `nes_scatter(p_threshold = )` | `plot_scatter(term_threshold = )` |
| `nes_scatter(color_by = )` | `plot_scatter(colour_by = )` |
| `nes_scatter(max_labels = )` | `plot_scatter(label_n = )` |
| `volcano_ring_theme()` | `plot_theme()` |
| `volcano_ring_theme(nes_colors = )` | `plot_theme(score_colours = )` |
| `volcano_ring_theme(nes_limits = )` | `plot_theme(score_limits = )` |
| `volcano_ring_theme(nes_stops = )` | `plot_theme(score_stops = )` |
| `ev_clean_label(name)` | `clean_label(name)` |
| `load_gene_sets(collections = )` | `load_gene_sets(databases = )` |
| `enrichment(results, metadata)` | `as_enrichment()`; the constructor is internal |

## New features

* `plot_bias_ring()` draws every significant term of one contrast, without
  names, with fill and arc height scaled to the strongest term. It shows at
  a glance whether a contrast leans up or down.
* The ring draws the twelve most significant unique terms by default. A term
  found in two collections counts once.
* `write_plot()` writes a named list of plots to one PDF: a lettered
  composite on page 1, then one panel per page. Fonts are embedded.
* `write_table()` writes one enrichment object, or the list from
  `run_enrichment()`, to one CSV that `as_enrichment()` reads back.

## Other changes

* Plot functions refuse plain volcano tables. Pass `as_da()` output.
* Points without a gene symbol are labelled with their accession.
* An invalid colour raises `enrichVolcano_param_error`.
* `dedup_terms()` records its cutoff as `term_threshold` in
  `metadata$dedup`.
* `CITATION` and `REFERENCES.bib` are gone. `citation("enrichVolcano")`
  reads DESCRIPTION, and each method's paper is on the help page of the
  function that runs it.

# enrichVolcano 1.1.0

## New features

* `read_study()` reads a study from one Excel workbook or a folder of CSV
  files: `da_results`, plus optionally `matrix`, `samples`, `contrasts` and
  `weights`. Every link between sheets is checked by name.
* `as_da()` gives DA results from limma, limpa, proteoDA, MSstats, proDA,
  msqrob2, prolfQua and ProtRank one set of column names, ranks proteins by
  the moderated t (else the signed -log10 p-value) and looks up current gene
  symbols for UniProt accessions in `org.Hs.eg.db`, `org.Mm.eg.db` or
  `org.Rn.eg.db`.
* `load_gene_sets()` loads Hallmark, Reactome, KEGG MEDICUS and GO:BP from
  msigdbr, and the GO Consortium generic slim (pinned, release 2026-07-26),
  for human, mouse or rat, and records the versions used.
* `run_enrichment()` runs fgsea, camera and fry on a study and returns one
  `enrichment` object per test. camera and fry refit limma from the sample
  sheet and contrasts, with subjects blocked or fixed, optional covariates
  and precision weights; a blocked camera runs as `cameraPR()`.
  `inter_gene_cor` sets camera's inter-gene correlation or estimates it.
* `example_study()` loads seven example studies; three include the full
  sample-level data.

## Deprecated

* `volcano_ring()` and `volcano_ring_grid()` take `as_da()` output. Plain
  tables and the `gene_col`, `logfc_col`, `pval_col` and `padj_col`
  arguments still work with a warning and will be removed in 2.0.0.

## Other changes

* `yvo_da.csv.gz` is replaced by `example_study("yvo")`, a full study
  (matrix, weights, samples, contrasts and DA results) quantified with limpa.
* GO slim sets are biological-process terms only and are named like
  `GOSLIM_PROTEIN_FOLDING`.
* The example data are licensed CC BY 4.0
  (`inst/extdata/studies/LICENSE.md`); the code stays MIT.

# enrichVolcano 1.0.0

## Breaking changes

* Plots read an `enrichment` object built by `as_enrichment()` instead of a
  data frame described by column-name arguments. `volcano_ring()` loses
  `term_col`, `nes_col`, `size_col`, `genes_col` and `genes_sep`, and gains
  `contrast`; `volcano_ring_grid()` takes one `enrichment` for all contrasts.
* `volcano_ring()` now chooses its terms: the `n_terms` (12) most significant
  below `term_threshold`, in either direction, from `databases` (all by
  default), with redundant terms hidden (`collapse`). `terms` hand-picks
  instead.
* `magnitude` defaults to `"neg_log_padj"` for NES and `"size"` for other
  scores; the legend is titled with the score type.
* The three vignettes are merged into one, `vignette("enrichVolcano")`.
* Example data: `yvo_fgsea.csv.gz` (fgsea output for four contrasts across
  MSigDB Hallmark, KEGG MEDICUS, Reactome and GO:BP from msigdbr 26.1.0, and
  the GO Consortium generic GO slim, release 2026-07-26) and `yvo_da.csv.gz`
  replace `yvo_enrichment.csv` and `yvo_da.csv`. The unused mito and cvh
  examples are removed.

## New features

* `as_enrichment()` converts fgsea, clusterProfiler `gseaResult`, limma
  `fry`, `mroast`, `camera` and `cameraPR`, over-representation tables with a
  direction column, and custom tables into one validated `enrichment` (an S7
  class). The producing test is recognised from the columns and recorded;
  camera and cameraPR, whose columns are identical, must be named with
  `enrichment_test`. Tests without an effect size are scored as signed
  -log10(FDR). clusterProfiler results are also read when exported as a data
  frame; a long table of limma results needs a `term` column, because row
  names do not survive stacking; `dedup_status` flags carried in the input
  are kept and recorded as `precomputed`.
* `dedup()` flags redundant terms for display without dropping rows or
  changing p-values: `method = "enrichmentmap"` (EnrichmentMap combined
  coefficient at 0.375 by default, or Jaccard at 0.5) or
  `method = "collapse_pathways"` (`fgsea::collapsePathways()`).
* `nes_scatter()` plots two contrasts term by term, as concordance or, with
  `comparison = "reversal"`, as reversal, with quadrant counts, Spearman's rho
  and the share of concordant or reversed terms. `color_by` and `shape_by`
  map any per-term column.
* `ev_clean_label()` gains `width`.
* For scores other than NES, the ring's fill scale spans the object's largest
  significant score rather than squishing at 3.

## Dependencies

* Imports S7. Suggests correlation, data.table, fgsea and withr.

## Also in this release

* `volcano_ring_theme()` takes direct colour overrides: `up`, `down`, `ns`
  for the points and `nes_colors` for the arc ramp, so custom palettes no
  longer need editing the returned list by hand.
* `volcano_ring()` gains three layout controls: `arc_order` (`"padj"` or
  `"nes"`) sets the angular order of arcs within each half,
  `arc_height_range` sets the shortest and tallest arc, and `show_counts`
  toggles the up/down count badges.
* `volcano_ring_grid()` exposes the composite-layout knobs: `panel_spacing`,
  `panel_margin`, `label_headroom`, `legend_position` (`"bottom"`, `"right"`,
  `"none"`), and `legend_width`. The shared NES legend now collects along the
  bottom by default.
* `volcano_ring()` gains `ring_thickness` and `tick_width` for the
  leading-edge tick band, `label_gap` for the gap between each arc and its
  label, and `count_size` / `axis_size` for the badge and axis-annotation
  text. Pathway labels now sit a fixed distance above their own arc, so every
  leader line is the same short length.

## Internal

* Pathway-name cleaning (`ev_clean_label()` and helpers) moved from
  `volcano_ring.R` to `labels.R`; `volcano_ring_theme()` no longer returns an
  unused theme element and now honours `base_family`. No behaviour change.

# enrichVolcano 0.3.0

## Breaking changes

* enrichVolcano is now a plotting-only package. Enrichment computation has
  been removed. Compute your enrichment with `fgsea::fgseaMultilevel()`,
  `clusterProfiler::gseGO()`, `enrichR::enrichr()`, or any tool you like,
  and pass the resulting tidy table to `volcano_ring()`.

## Removed (use v0.2.0 to recover)

* Enrichment engine: `ev_enrich()`, `ev_collapse()`, `list_databases()`,
  `database_info()`, `load_go_slim()`, `apply_gate()`.
* Scoring + adjustment: `pi_score()`, `adjust_p()`.
* Input layer: `ev_read_contrasts()`, `ev_validate()`, `ev_idmap_report()`.
* Hero wrapper: `enrich_volcano()`.
* Shiny app: `ev_app()` (relocated to the sibling `enrichVolcanoApp` repo).

## Renamed

* `ev_volcano_ring()` -> `volcano_ring()`
* `ev_compose()` -> `volcano_ring_grid()` (signature changed; now takes
  paired lists of tidy DA + enrichment frames keyed by contrast, and returns
  an S3 `volcano_ring_grid` object carrying `$plot` and `$data`).
* `ev_theme()` -> `volcano_ring_theme()` (new args: `base_family`,
  `nes_limits`, `nes_stops`; palettes: `"default"`, `"viridis"`, `"okabe"`).
* `print.enrichVolcano()` -> `print.volcano_ring_grid()`.
* `ev_volcano()` and `ring_plot()` removed. The composite is the
  plot.

## New

* Column-naming arguments on `volcano_ring()`: `gene_col`, `logfc_col`,
  `pval_col`, `padj_col`, `term_col`, `nes_col`, `size_col`, `genes_col`,
  `genes_sep`. Defaults assume limma + fgsea conventional names.
* Tick-line column auto-detect: `leading_edge` (`;`-string),
  `leadingEdge` (list-col), `core_enrichment` (`/`-string), and `Genes`
  (`;`-string) all work out of the box. Override via `genes_col`.
* Structured error classes: `enrichVolcano_input_error`,
  `enrichVolcano_column_error`, `enrichVolcano_data_error`,
  `enrichVolcano_param_error`, all under the parent
  `enrichVolcano_error`.
* `magnitude = c("neg_log_padj", "size")` controls the arc-thickness
  encoding on the ring.

## Dependencies

* Imports trimmed from 18 to 8 (kept: cli, ggforce, ggplot2, ggrepel,
  patchwork, rlang, scales, stringr).
* Suggests trimmed from 15 to 6 (kept: covr, knitr, pkgdown, rmarkdown,
  testthat, vdiffr).

# enrichVolcano 0.2.0

## Breaking changes

* The default fgsea ranking statistic (`rank_by`) is now `"t"`, the limma /
  proteoDA moderated t-statistic, instead of `"signed_p"`. This is the
  recommended signed GSEA statistic for moderated linear models (Subramanian
  2005 PNAS; Reimand 2019 Nat Protoc) and matches the hand-written source
  pipelines. Inputs without a `t` column (e.g. edgeR) fall back to `"signed_p"`
  automatically, so only limma-family results change. Pass `rank_by = "signed_p"`
  to restore the previous behaviour.

## Bug fixes

* `ev_enrich(nperm = ...)` now takes effect. It is passed to
  `fgsea::fgseaMultilevel()` as `nPermSimple` (previously a no-op).
* `ev_collapse(keep_by = "NES")` now ranks representatives by |NES|, so a
  strongly down-regulated pathway is no longer dropped in favour of a weakly
  up-regulated one.
* `ev_collapse()` no longer collapses an up- and a down-regulated pathway that
  share leading-edge genes. Dedup is now within-direction.
* `adjust_p(method = "qvalue")` falls back to BH (with a warning) instead of
  erroring when a contrast has too few p-values to estimate pi0.
* `ev_validate()` guards MaxQuant ratios `<= 0` (set to `NA` with a warning)
  instead of producing `-Inf`/`NaN` log2 values.
* `ev_read_contrasts()` no longer creates a duplicate `adj.P.Val` column when
  `padj` collides with an existing one.
* `ev_compose()` matches volcano/ring lists by name (order-insensitive) rather
  than requiring identical name order.

## New features

* `ev_collapse(scope = ...)` accepts a vector of stages run in order, e.g.
  `c("within_db", "cross_db")`: each stage collapses only the survivors of the
  previous, so redundancy can be cleaned within each database first and then
  merged across databases. `cutoff` may be one value per stage. (Within/across
  database redundancy control follows Vivar 2013 ReCiPa, PMID 23758478.)
* `ev_volcano_ring()` and `enrich_volcano()` gain `nes_limits` (default
  `c(-3, 3)`) to widen the NES colour scale when enrichment exceeds |NES| = 3.
* `ev_collapse()` gains a `similarity` argument selecting the gene-overlap metric
  for the Jaccard-family steps: `"jaccard"` (default), `"overlap"`
  (\eqn{|A\cap B|/\min(|A|,|B|)}, catches a small set contained in a larger one),
  or `"combined"` (the Cytoscape EnrichmentMap coefficient
  \eqn{w\cdot Jaccard + (1-w)\cdot Overlap}; Merico 2010, PMID 21085593), with
  `combined_weight` (default 0.5). Forwarded through `enrich_volcano(dedup = ...)`.
  Default behaviour is unchanged (`similarity = "jaccard"`).
* `ev_collapse()` default method is now `"collapse_then_jaccard"`, matching the
  source pipeline that generated the YvO 2025 figures. The previous default
  `"jaccard"` remains available; the alias `"both"` is deprecated and will be
  removed in a future release.
* `ev_collapse()` gains `sig_threshold` (default 0.05) so non-significant
  pathways are no longer affected by dedup. Pass `NA` to restore the previous
  dedup-all-rows behavior.
* `ev_enrich()` gains `include_terms` and `filter_mode = c("before","display")`
  for pre-test and post-hoc pathway-name filtering. See
  `vignette("pathway-dedup")` for when to use which.
* `ev_enrich()` output now carries an `ev_filter` attribute recording the
  filter inputs and per-database pathway counts before and after filtering.

## Bug fixes

* `enrich_volcano()` no longer pins `method = "jaccard"` in its `dedup`
  default; the hero function now inherits whatever default `ev_collapse()`
  ships, so users calling the wrapper with defaults get the documented
  `"collapse_then_jaccard"` behavior.
* `ev_collapse_fgsea()` now looks up pathway gene sets by `(database,
  pathway)` rather than by name alone. The previous flatten preserved
  duplicate pathway names across databases, so a user combining (for example)
  a custom GMT and an MSigDB collection that shared a pathway name would get
  the first database's gene set for both rows.
* `ev_collapse_fgsea()` now partitions its input by contrast and runs
  `fgsea::collapsePathways` once per contrast against the matching gene-level
  rank vector. Previously, under `scope = "global"` with multiple contrasts,
  every pathway was silently scored against the first contrast's ranks.
* `ev_collapse(method = "jaccard_then_collapse")` (and the deprecated
  `"both"` alias) now passes only the Jaccard survivors to the collapse step
  and writes back positionally, mirroring `"collapse_then_jaccard"`. The
  previous AND-merge could drop both members of a redundant cluster when
  `collapsePathways` and Jaccard's `keep_by` tie-breaker disagreed on which
  row represented the cluster.
* `ev_enrich()` `ev_filter$n_pathways_after` is now structurally parallel to
  `n_pathways_before` even when every contrast is empty (previously
  `integer(0)` instead of a named integer keyed by database).
