# Build inst/extdata/db_registry.csv from scratch.
# Not shipped to users (in .Rbuildignore).

library(tibble)
library(readr)

db_registry <- tibble::tribble(
  ~name,                ~display_name,             ~type,        ~species,                                                ~source,             ~license,            ~citation_key,        ~description_short,
  "hallmark",           "MSigDB Hallmark",         "pathway",    "human;mouse;rat;zebrafish;fly;yeast;pig",               "msigdbr",           "MSigDB Academic",   "msigdb_hallmark",    "50 hallmark gene sets representing well-defined biological states",
  "reactome",           "Reactome",                "pathway",    "human;mouse;rat;zebrafish;fly;yeast;pig",               "msigdbr",           "CC BY 4.0",         "reactome_2024",      "Canonical biological pathways (Reactome)",
  "kegg_legacy",        "KEGG Legacy",             "pathway",    "human;mouse;rat;zebrafish;fly;yeast;pig",               "msigdbr",           "KEGG Academic",     "kegg_kanehisa",      "Pre-2024 KEGG pathway maps",
  "kegg_medicus",       "KEGG Medicus",            "pathway",    "human;mouse;rat;zebrafish;fly;yeast;pig",               "msigdbr",           "KEGG Academic",     "kegg_kanehisa",      "Current KEGG Medicus pathway maps",
  "wikipathways",       "WikiPathways",            "pathway",    "human;mouse;rat;zebrafish;fly;yeast;pig",               "msigdbr",           "CC BY-SA 3.0",      "wikipathways_2024",  "Community-curated biological pathways",
  "biocarta",           "BioCarta",                "pathway",    "human;mouse;rat;zebrafish;fly;yeast;pig",               "msigdbr",           "Academic",          "biocarta_2002",      "Classic small canonical pathways",
  "pid",                "NCI Pathway Interaction", "pathway",    "human;mouse;rat;zebrafish;fly;yeast;pig",               "msigdbr",           "Academic",          "pid_schaefer_2009",  "NCI Pathway Interaction Database",
  "go_bp",              "GO Biological Process",   "ontology",   "human;mouse;rat;zebrafish;fly;yeast;pig",               "msigdbr",           "CC BY 4.0",         "go_consortium",      "Biological process ontology",
  "go_cc",              "GO Cellular Component",   "ontology",   "human;mouse;rat;zebrafish;fly;yeast;pig",               "msigdbr",           "CC BY 4.0",         "go_consortium",      "Cellular component (organelle-relevant)",
  "go_mf",              "GO Molecular Function",   "ontology",   "human;mouse;rat;zebrafish;fly;yeast;pig",               "msigdbr",           "CC BY 4.0",         "go_consortium",      "Molecular function ontology",
  "hpo",                "Human Phenotype Ontology","ontology",   "human",                                                 "msigdbr",           "CC BY 4.0",         "hpo_kohler_2021",    "Human Phenotype Ontology",
  "hallmark_oncogenic", "MSigDB Oncogenic (C6)",   "pathway",    "human",                                                 "msigdbr",           "Academic",          "msigdb_c6",          "Oncogenic signatures",
  "immunesigdb",        "ImmuneSigDB (C7)",        "cell_type",  "human;mouse",                                           "msigdbr",           "Academic",          "msigdb_c7",          "Immune cell signatures",
  "cell_type_sig",      "Cell Type Signatures (C8)","cell_type", "human",                                                 "msigdbr",           "Academic",          "msigdb_c8",          "Curated cell-type signatures",
  "mirdb",              "miRDB targets",           "regulatory", "human;mouse;rat;zebrafish;fly;yeast;pig",               "msigdbr",           "CC BY-NC 4.0",      "mirdb_2020",         "miRNA target gene sets",
  "tft_gtrd",           "TFT GTRD",                "regulatory", "human;mouse;rat;zebrafish;fly;yeast;pig",               "msigdbr",           "CC BY 4.0",         "tft_gtrd_2019",      "Transcription factor target genes",
  "mitocarta3",         "MitoCarta 3.0",           "organelle",  "human;mouse",                                           "bundled_gmt",       "CC BY-NC 4.0",      "mitocarta_rath_2021","Mitochondrial proteins + 149 MitoPathways",
  "corum",              "CORUM 5.0",               "complex",    "human;mouse;rat",                                       "bundled_gmt",       "CC BY 4.0",         "corum_tschen_2023",  "Mammalian protein complexes (proteomics-native gene sets)",
  "go_slim",            "GO Slim (generic BP)",    "ontology",   "human;mouse;rat",                                       "bundled_gmt",       "CC BY 4.0",         "go_consortium",      "62 generic GO Slim biological-process categories"
)

dir.create("inst/extdata", recursive = TRUE, showWarnings = FALSE)
readr::write_csv(db_registry, "inst/extdata/db_registry.csv")
