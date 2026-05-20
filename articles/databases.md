# Choosing an enrichment database

``` r

library(enrichVolcano)
```

## The 21 registered databases

[`list_databases()`](https://Dustyn-T-Lewis.github.io/enrichVolcano/reference/list_databases.md)
returns a tibble with one row per registered database. The full registry
ships with the package in `inst/extdata/db_registry.csv` and is the
single source of truth for the rest of the API.

``` r

db <- list_databases()
db[, c("name", "display_name", "type", "species", "license")]
#> # A tibble: 18 × 5
#>    name               display_name              type       species       license
#>    <chr>              <chr>                     <chr>      <chr>         <chr>  
#>  1 hallmark           MSigDB Hallmark           pathway    human;mouse;… MSigDB…
#>  2 reactome           Reactome                  pathway    human;mouse;… CC BY …
#>  3 kegg_legacy        KEGG Legacy               pathway    human;mouse;… KEGG A…
#>  4 kegg_medicus       KEGG Medicus              pathway    human;mouse;… KEGG A…
#>  5 wikipathways       WikiPathways              pathway    human;mouse;… CC BY-…
#>  6 biocarta           BioCarta                  pathway    human;mouse;… Academ…
#>  7 pid                NCI Pathway Interaction   pathway    human;mouse;… Academ…
#>  8 go_bp              GO Biological Process     ontology   human;mouse;… CC BY …
#>  9 go_cc              GO Cellular Component     ontology   human;mouse;… CC BY …
#> 10 go_mf              GO Molecular Function     ontology   human;mouse;… CC BY …
#> 11 hpo                Human Phenotype Ontology  ontology   human         CC BY …
#> 12 hallmark_oncogenic MSigDB Oncogenic (C6)     pathway    human         Academ…
#> 13 immunesigdb        ImmuneSigDB (C7)          cell_type  human;mouse   Academ…
#> 14 cell_type_sig      Cell Type Signatures (C8) cell_type  human         Academ…
#> 15 mirdb              miRDB targets             regulatory human;mouse;… CC BY-…
#> 16 tft_gtrd           TFT GTRD                  regulatory human;mouse;… CC BY …
#> 17 mitocarta3         MitoCarta 3.0             organelle  human;mouse   CC BY-…
#> 18 corum              CORUM 5.0                 complex    human;mouse;… CC BY …
```

Twenty-one rows cover seven database types. The next section walks
through each type and gives a one-line use-case.

## Database types

### Pathway

Curated biological pathways with explicit gene memberships. Examples:
`hallmark`, `reactome`, `kegg_legacy`, `kegg_medicus`, `wikipathways`,
`biocarta`, `pid`. Use these to answer the question “what biological
process is shifting?”. Hallmark is the smallest collection (50 sets) and
gives a fast first pass. Reactome and KEGG_MEDICUS are larger and more
granular; expect more results and more redundancy.

Sub-choices matter. Hallmark sets are deliberately broad (e.g.
HALLMARK_OXIDATIVE_PHOSPHORYLATION has 200 genes). Reactome breaks the
same biology into a hierarchy from “Metabolism” down through
“Respiratory electron transport” to “Complex I biogenesis”, with
gene-set sizes from a few thousand to under twenty. KEGG_LEGACY is the
pre-2024 KEGG release; KEGG_MEDICUS is the redesigned 2024 collection
with cleaner network logic. WikiPathways skews toward community-curated
metabolism and signalling. BioCarta sets are small, old, and rarely the
right primary source but useful as cross-reference. PID covers
cancer-relevant signalling specifically.

### Ontology

Controlled vocabulary trees from the Gene Ontology and HPO projects.
Examples: `go_bp`, `go_cc`, `go_mf`, `hpo`. Use for fine-grained term
enrichment when you want to drill into a specific functional axis.
Expect heavy parent-child redundancy: GO_BP can return 30 enriched terms
that all describe the same process at different depths. Pair with
[`ev_collapse()`](https://Dustyn-T-Lewis.github.io/enrichVolcano/reference/ev_collapse.md)
(see
[`vignette("pathway-dedup")`](https://Dustyn-T-Lewis.github.io/enrichVolcano/articles/pathway-dedup.md))
before plotting.

The three GO branches answer different questions. GO_BP describes
processes (“response to oxidative stress”). GO_CC describes compartments
(“mitochondrial inner membrane”). GO_MF describes molecular activities
(“ATP-binding”). For a proteomics aging study, GO_BP usually carries the
headline biology; GO_CC is the natural follow-up when localization
changes are suspected; GO_MF is best reserved for biochemistry-heavy
studies (enzyme classes, transporter families). HPO is human-only and
links genes to clinical phenotypes; use it when an aging signature
overlaps with a known disease phenotype set.

### Organelle

Subcellular localization sets. Examples: `mitocarta3`,
`hpa_subcellular`, `uniprot_subcellular`. Use when the question is
“which cellular compartment is affected?”. MitoCarta3 has 149 curated
mitochondrial sub-pathways on top of the master mitochondrial gene list.
HPA and UniProt cover the rest of the cell.

MitoCarta3 is the strongest organelle source for mitochondrial work in
human, mouse, and rat. The sub-pathways resolve OXPHOS complexes
individually, separate matrix from intermembrane space, and split
beta-oxidation by chain length. HPA_subcellular comes from the Human
Protein Atlas immunofluorescence atlas and gives compartment calls
backed by image evidence. UniProt_subcellular is broader in species
coverage but coarser; it pulls keyword annotations from SwissProt and
TrEMBL.

### Cell type

Signature sets from sorted cell populations. Examples: `immunesigdb`
(MSigDB C7), `cell_type_sig` (MSigDB C8). Use for tissue-mixing
questions: a bulk muscle biopsy contains satellite cells, immune
infiltrate, and endothelium, and a strong enrichment in one of these
signatures often points at composition shift rather than per-cell
regulation.

### Regulatory

TF and miRNA target sets. Examples: `tft_gtrd`, `mirdb`. Use to ask “is
this driven by a specific transcription factor?”. GTRD is built from
ChIP-seq data and gives a binding-based view; miRDB gives the miRNA
seed-match prediction. Treat results as hypothesis-generating.

### Disease

Disease-gene associations. Examples: `disgenet`, `opentargets_disease`.
Use for disease overlap questions: does the aging muscle proteome share
more genes with sarcopenia than with type-2 diabetes? Read with caution.
Pathway databases such as KEGG also include disease-named entries, so
the same gene can show up in both a “pathway” hit and a “disease” hit.
The default `exclude_terms = "DISEASE|CANCER|TUMOR"` regex in
[`enrich_volcano()`](https://Dustyn-T-Lewis.github.io/enrichVolcano/reference/enrich_volcano.md)
strips disease-named entries out of pathway databases so they do not
double-count.

## A decision flowchart

A short tour for picking a starting point:

- First-pass discovery: `hallmark`. Fifty sets cover most coarse
  biological states and return a readable ring panel.
- Curated canonical pathways: `reactome` or `kegg_medicus`.
- Deep ontology with dedup: `go_bp` plus
  `ev_collapse(scope = "within_db")`.
- Mitochondrial focus: `mitocarta3` plus `go_cc` for cross-check.
- Subcellular localization: `hpa_subcellular` or `uniprot_subcellular`.
- TF activity: `tft_gtrd`.
- Immune phenotype: `immunesigdb`.
- Cell-type composition: `cell_type_sig`.
- Disease overlap: `disgenet` or `opentargets_disease`.

A typical paper uses two or three sources. Hallmark for the headline
ring, Reactome or GO_BP for the deep dive, and one organelle or
regulatory source if the biology calls for it.

Beware of “throw every database at it” mode. Twenty-one databases run
against one contrast returns hundreds of overlapping hits and forces you
to triage them post-hoc. Pick two or three sources up front; the ring
panel can only show eight to twelve terms anyway. For supplementary
tables, run everything separately, save the unfiltered enrichment
object, and pick top hits per source.

## Per-database `database_info()` output

[`database_info()`](https://Dustyn-T-Lewis.github.io/enrichVolcano/reference/database_info.md)
prints the registry row for one database in a readable format. It is the
right function to call when a reviewer asks where a gene set came from.

``` r

database_info("hallmark")
#> 
#> ── MSigDB Hallmark (hallmark) ──
#> 
#> Type: pathway
#> Species: human;mouse;rat;zebrafish;fly;yeast;pig
#> Source: msigdbr
#> License: MSigDB Academic
#> Description: 50 hallmark gene sets representing well-defined biological states
```

Hallmark is the entry point for almost every analysis in this package.
The sets are coarse enough that a single contrast rarely hits more than
five with strong NES, which keeps the ring panel readable. Pair with
Reactome or GO_BP for the supplementary deep dive.

``` r

database_info("reactome")
#> 
#> ── Reactome (reactome) ──
#> 
#> Type: pathway
#> Species: human;mouse;rat;zebrafish;fly;yeast;pig
#> Source: msigdbr
#> License: CC BY 4.0
#> Description: Canonical biological pathways (Reactome)
```

Reactome ships through `msigdbr` here. It has more than 1600 sets in the
2024 release. Expect overlap with KEGG, WikiPathways, and BioCarta on
canonical processes;
[`ev_collapse()`](https://Dustyn-T-Lewis.github.io/enrichVolcano/reference/ev_collapse.md)
handles that automatically.

``` r

database_info("go_bp")
#> 
#> ── GO Biological Process (go_bp) ──
#> 
#> Type: ontology
#> Species: human;mouse;rat;zebrafish;fly;yeast;pig
#> Source: msigdbr
#> License: CC BY 4.0
#> Description: Biological process ontology
```

GO_BP is the gold standard for fine-grained biological process
enrichment. It is also the worst offender for parent-child redundancy.
Always pair with Jaccard dedup before plotting.

``` r

database_info("mitocarta3")
#> 
#> ── MitoCarta 3.0 (mitocarta3) ──
#> 
#> Type: organelle
#> Species: human;mouse
#> Source: bundled_gmt
#> License: CC BY-NC 4.0
#> Description: Mitochondrial proteins + 149 MitoPathways
```

MitoCarta3 ships as a bundled GMT file inside the package, not as a
fetched msigdbr collection. It covers human, mouse, and rat. The 149
MitoPathways subdivide the master mitochondrial gene list into oxidative
phosphorylation complexes, import machinery, beta-oxidation, and so on.
License is CC BY-NC: cite the Rath 2021 paper and stay non-commercial.

## Multi-species support

The package supports seven species: human, mouse, rat, zebrafish, fly,
yeast, and pig. The species argument flows through both fetched and
bundled databases.

For msigdbr-sourced collections (Hallmark, Reactome, KEGG, GO, etc.),
species names are resolved to scientific names via `msigdbr`’s ortholog
table. Passing `species = "rat"` resolves to `Rattus norvegicus`. Human
gene symbols are the canonical key; non-human collections are derived by
orthology mapping from human symbols using HCOP and msigdbr’s curated
tables.

Since msigdbr 10, the underlying database is selected with `db_species`
(`"HS"` or `"MM"`) and the output orthologs with `species`. This package
uses the human database (`db_species = "HS"`) and maps to your
`species`, so for a mouse study you will see a one-time note that mouse
results come from human-derived orthologs rather than the mouse-native
MM collections. That is the intended default; human-keyed sets are the
better-curated and more widely cited reference.

The mapping is one-to-many in principle. One human gene can map to two
rat paralogs; one rat gene can map to none if the orthology is
unresolved. msigdbr handles this by emitting every mapped pair into the
species-specific collection. The practical consequence is that a rat
collection can be slightly larger or smaller than its human counterpart,
and a few rare paralog-only genes will appear in rat but not human. This
is normal and rarely changes top-line enrichment.

For bundled GMT files (MitoCarta, HPA, UniProt), the package selects the
species-specific file from `inst/extdata/gmt/`. The naming convention is
`<db>_<species>.gmt`. If a species file is missing,
[`enrich_volcano()`](https://Dustyn-T-Lewis.github.io/enrichVolcano/reference/enrich_volcano.md)
aborts with a clear message rather than silently falling back to human.

A practical implication: if your gene symbols are not in the canonical
species namespace (rat symbols for a mouse study, for instance),
enrichment will be empty or noise. Re-key the input before calling
[`enrich_volcano()`](https://Dustyn-T-Lewis.github.io/enrichVolcano/reference/enrich_volcano.md).

## Custom GMT files

For v0.1, pass a named list of pathway lists directly through the
`databases` argument:

``` r

my_paths <- list(
  CUSTOM_SET_1 = c("ENO1", "PKM", "GAPDH"),
  CUSTOM_SET_2 = c("ATP5F1A", "ATP5F1B", "ATP5F1C")
)
p <- enrich_volcano(
  data       = yvo,
  contrast   = "Aging",
  databases  = list(custom = my_paths)
)
```

The outer list name (here `custom`) becomes the `database` column in the
enrichment output. The inner names become pathway labels in the ring
panel.

A full GMT path argument (`custom_gmt = "path/to/sets.gmt"`) is reserved
for v0.2. For now, read the file with
[`fgsea::gmtPathways()`](https://rdrr.io/pkg/fgsea/man/gmtPathways.html)
and pass the resulting list directly.

## Licensing

License terms vary by source and apply to redistribution, not to running
the analysis. The registry summarises each:

- **MSigDB Academic**: Hallmark and most C2-C8 collections require
  academic use. Commercial use needs a Broad Institute license.
- **MitoCarta CC BY-NC 4.0**: non-commercial with attribution to Rath et
  al. 2021.
- **GO, Reactome, UniProt**: open under CC BY 4.0. Cite the consortium.
- **KEGG Academic**: free for academic use; commercial use needs KEGG
  Inc. license.
- **DisGeNET CC BY-NC-SA 4.0**: non-commercial, share-alike, attribute
  Piñero et al. 2020.

The full license per database is in `list_databases()$license` and
`inst/extdata/db_registry.csv`. The package itself is MIT-licensed; the
licenses on derived gene-set collections are separate.

## Recipes by study type

A few starting points tied to study patterns.

Aging muscle proteomics (YvO 2025 pipeline): Hallmark for the headline
ring, Reactome for the deep dive, GO_CC for compartment shifts,
MitoCarta3 for OXPHOS detail. Four sources, run separately, combined
with `scope = "within_db"` so each source contributes a distinct slice
to the ring.

Cancer vs healthy serum proteomics (CvH 2026 pipeline): Hallmark and
Reactome as primary; HPA_subcellular to separate secreted from cell-leak
signals; ImmuneSigDB to flag immune-driven changes. Disease databases
(`disgenet`, `opentargets_disease`) are useful here for the
comparison-to-known-disease overlap.

High vs low responder training studies (HRvLR 2026 pipeline): Hallmark
and Reactome with `scope = "cross_db"` since responders and
non-responders share most biology and the goal is non-redundant top-line
hits. Add `tft_gtrd` to ask whether a TF program separates the two
groups.

H9c2 mitochondrial transplant (Mito 2026 conceptual pipeline):
MitoCarta3 as the primary source; GO_CC for membrane vs matrix splits;
Reactome for OXPHOS and TCA cycle detail. Hallmark is too coarse for the
mitochondrial focus and best dropped.

These recipes are starting points, not prescriptions. The right database
mix depends on the biology and the audience. For a methods section, name
every source and cite the registry version.

## Comparing across runs

Re-running
[`enrich_volcano()`](https://Dustyn-T-Lewis.github.io/enrichVolcano/reference/enrich_volcano.md)
with the same data and the same database set returns the same enrichment
table. `msigdbr` ships versioned releases, so locking the package
version locks the gene sets. For long-running projects, record the
`msigdbr` version alongside the analysis script:

``` r

sessionInfo()$otherPkgs$msigdbr$Version
```

Bundled GMT files are versioned through the package version. The 2024
release of MitoCarta3 ships with v0.1 of enrichVolcano. A future
Mitocarta3.1 will bump the enrichVolcano minor version.

Next: how to style your plots — see
[`vignette("customising")`](https://Dustyn-T-Lewis.github.io/enrichVolcano/articles/customising.md).
