# enrichVolcano 0.2.0

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
