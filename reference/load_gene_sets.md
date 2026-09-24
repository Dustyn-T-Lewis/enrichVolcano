# Load gene-set collections for enrichment

Fetches collections at the versions installed on this machine and
records those versions, so an analysis can state exactly what it tested.

## Usage

``` r
load_gene_sets(
  collections = c("Hallmark", "GO Slim"),
  species = "Homo sapiens",
  min_size = 15,
  max_size = 500
)
```

## Arguments

- collections:

  Any of `"Hallmark"`, `"GO Slim"`, `"Reactome"`, `"KEGG"`, `"GO:BP"`.

- species:

  `"Homo sapiens"`, `"Mus musculus"` or `"Rattus norvegicus"`.

- min_size, max_size:

  Keep sets with this many genes.

## Value

A named list of collections, each a named list of gene symbols, with a
`versions` attribute (msigdbr, GO slim release, annotation package,
species).

## Details

- `"Hallmark"`, `"Reactome"`, `"KEGG"` (KEGG MEDICUS) and `"GO:BP"` come
  from msigdbr; mouse and rat sets are human sets mapped through
  orthologs.

- `"GO Slim"` is the biological-process part of the GO Consortium's
  generic slim, pinned in this package (release 2026-07-26), with each
  term's genes taken from the species' annotation package, counting
  genes annotated to the term or any of its descendants. Sets are named
  like `GOSLIM_PROTEIN_FOLDING`.

Each collection is kept separate so
[`run_enrichment()`](https://Dustyn-T-Lewis.github.io/enrichVolcano/reference/run_enrichment.md)
corrects p-values within it.
