# extdata/goslim_generic.obo — Provenance

This is the canonical **generic GO slim** from the Gene Ontology Consortium,
shipped frozen in this bundle for reproducibility.

| Field           | Value                                                          |
|-----------------|----------------------------------------------------------------|
| Source URL      | http://current.geneontology.org/ontology/subsets/goslim_generic.obo |
| Retrieved       | 2026-06-15                                                     |
| Last-Modified   | 2026-05-28                                                     |
| GO release      | go/releases/2026-05-19                                         |

This file is the authoritative, version-pinned artifact used by `load_go_slim()`
to derive biological_process slim term IDs at runtime. It is NOT hand-curated
or hand-typed — it is the GO Consortium's published generic GO slim, downloaded
directly from the official GO Consortium release endpoint and committed frozen
so that any run of the pipeline produces the same term set.
