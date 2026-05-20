# Run enrichment (fgsea + ORA) against one or many databases

Run enrichment (fgsea + ORA) against one or many databases

## Usage

``` r
ev_enrich(
  data,
  contrast,
  databases,
  species = "human",
  enrich_mode = c("fgsea", "ora"),
  rank_by = "pi_eq2",
  min_size = 10,
  max_size = 500,
  background = NULL,
  exclude_terms = "DISEASE|CANCER|TUMOR",
  nperm = 10000,
  seed = 42
)
```

## Arguments

- data:

  Tidy long tibble (post
  [`pi_score()`](https://Dustyn-T-Lewis.github.io/enrichVolcano/reference/pi_score.md)
  /
  [`adjust_p()`](https://Dustyn-T-Lewis.github.io/enrichVolcano/reference/adjust_p.md)).

- contrast:

  Character; one or many contrast names.

- databases:

  Either character vector of registered DB names (resolved via the
  registry), a named list of pathway lists, or path(s) to GMT files.

- species:

  Character (for registered DBs).

- enrich_mode:

  `c("fgsea", "ora")` — either or both.

- rank_by:

  Column to rank genes by for fgsea.

- min_size, max_size:

  Gene-set size filters.

- background:

  Optional character vector of background genes for ORA; default is all
  genes in `data` for the given contrast.

- exclude_terms:

  Regex of pathway names to exclude (default removes
  DISEASE/CANCER/TUMOR terms).

- nperm:

  Permutations for fgsea (default 10000).

- seed:

  Integer seed.

## Value

Tidy tibble: contrast \| database \| pathway \| pval \| padj \| NES \|
size \| leading_edge \| mode \| direction. fgsea-mode rows also carry
`ES` and `log2err` (NA for ORA rows). The pathway list used for
enrichment is attached as `attr(., "ev_pathways")` and the per-contrast
gene-level ranking vectors as `attr(., "ev_stats")`, both required by
`ev_collapse(method = "collapse_fgsea")`.

## Details

Each entry in `databases` (or each registered name) produces fgsea
and/or ORA results. fgsea ranks genes by `rank_by` (default `pi_eq2`).
ORA uses significant genes (`pi_eq2 < 0.05` if present, else
`adj.P.Val < 0.05`) against the supplied or inferred background.
