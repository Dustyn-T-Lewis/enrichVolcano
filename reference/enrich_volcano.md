# Build a composite volcano + enrichment ring plot

Hero function for `enrichVolcano`. Runs the full pipeline: validate
input -\> compute pi-score -\> adjust p-values -\> run enrichment -\>
deduplicate -\> build volcano + ring panels -\> compose patchwork.

## Usage

``` r
enrich_volcano(
  data,
  contrast,
  species = "human",
  databases = c("hallmark", "reactome", "go_bp"),
  x = NULL,
  y = NULL,
  lab = NULL,
  custom_gmt = NULL,
  p_method = c("pi_eq2", "pi_eq1", "raw_p", "adj_p"),
  p_adjust = "BH",
  p_threshold = 0.05,
  logfc_threshold = log2(1.5),
  enrich_mode = c("fgsea", "ora"),
  enrich_padj = 0.05,
  dedup = list(method = "jaccard", cutoff = 0.5, scope = "within_db"),
  ring = list(max_terms = 10, order_by = "padj", magnitude = "neg_log_padj", color =
    "nes"),
  volcano = list(label_n = 10, label_by = NULL),
  facet = list(nrow = NULL, ncol = NULL),
  theme = ev_theme()
)
```

## Arguments

- data:

  Differential abundance results. Accepts tidy long, wide-suffix, DEP,
  proteoDA, limma, DESeq2, edgeR, MSstats, proDA, DEqMS, MaxQuant, or
  Perseus formats -
  [`ev_validate()`](https://Dustyn-T-Lewis.github.io/enrichVolcano/reference/ev_validate.md)
  auto-detects.

- contrast:

  Character; one or many contrast names.

- species:

  Character; "human" (default), "mouse", "rat", "zebrafish", "fly",
  "yeast", "pig".

- databases:

  Character vector of registered DB names (see
  [`list_databases()`](https://Dustyn-T-Lewis.github.io/enrichVolcano/reference/list_databases.md)),
  or a named list of pathway lists.

- x, y, lab:

  Column names for logFC, p-value, gene label (auto-detected if NULL).

- custom_gmt:

  Reserved for v0.2 (GMT path input).

- p_method:

  Y-axis transform: `"pi_eq2"` (default), `"pi_eq1"`, `"raw_p"`,
  `"adj_p"`.

- p_adjust:

  Method: `"BH"` (default), `"bonferroni"`, `"qvalue"`, `"IHW"`.

- p_threshold, logfc_threshold:

  Volcano significance cutoffs.

- enrich_mode:

  `c("fgsea", "ora")` - either or both.

- enrich_padj:

  Pathway significance cutoff (default 0.05).

- dedup:

  List with `method` (`"jaccard"`, `"collapse_fgsea"`, or `"both"`),
  `cutoff`, and `scope`. See
  [`ev_collapse()`](https://Dustyn-T-Lewis.github.io/enrichVolcano/reference/ev_collapse.md).

- ring:

  List with `max_terms`, `order_by`, `magnitude`, `color`.

- volcano:

  List with `label_n`, `label_by`.

- facet:

  List with `nrow`, `ncol` for patchwork outer layout.

- theme:

  Output of
  [`ev_theme()`](https://Dustyn-T-Lewis.github.io/enrichVolcano/reference/ev_theme.md).

## Value

Object of class `c("enrichVolcano", "patchwork")` with
`attr(., "ev_data")` and `attr(., "ev_call")`.
