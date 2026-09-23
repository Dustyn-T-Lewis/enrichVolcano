# Builds inst/extdata/examples/yvo_fgsea.csv.gz, the example enrichment table.
#
# Source: fgsea on limma moderated-t rankings from the YvO study, all four
# contrasts, five collections tested as separate families (min size 15,
# max 500). Built by YvO_2026/04_Figures/shared/build_fgsea_cache.R with
#   msigdbr 26.1.0 (Hallmark, KEGG MEDICUS, Reactome, GO:BP)
#   GO slim generic, go/releases/2026-07-26, genes assigned via GO.db
# The dedup_status / merged_into / overlap_jaccard columns are that script's
# EnrichmentMap-style flags (Jaccard >= 0.5, within database).
#
# Rerun when the upstream cache changes, and record the new versions in NEWS.

source_csv <- file.path(
  Sys.getenv("YVO_DIR", "~/Desktop/A_Proteomics_Analysis/2026/YvO/YvO_2026"),
  "04_Figures", "shared", "fgsea_tstat_all_v2.csv"
)
fgsea_cache <- utils::read.csv(source_csv)

out <- gzfile(here::here("inst", "extdata", "examples", "yvo_fgsea.csv.gz"), "w")
utils::write.csv(fgsea_cache, out, row.names = FALSE)
close(out)
