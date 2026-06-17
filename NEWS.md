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

* `ev_enrich(nperm = ...)` now actually takes effect — it is passed to
  `fgsea::fgseaMultilevel()` as `nPermSimple` (previously a no-op).
* `ev_collapse(keep_by = "NES")` now ranks representatives by **|NES|**, so a
  strongly down-regulated pathway is no longer dropped in favour of a weakly
  up-regulated one.
* `ev_collapse()` no longer collapses an up- and a down-regulated pathway that
  share leading-edge genes — dedup is now within-direction.
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
