source(test_path("fixtures/make_toy.R"))

make_lbl_input <- function() {
  data.frame(
    gene = paste0("G", 1:10),
    logFC = c(3, 2, 1, 0.5, 0.1, -0.1, -0.5, -1, -2, -3),
    P.Value = c(1e-5, 1e-4, 1e-3, 0.01, 0.5, 0.5, 0.01, 1e-3, 1e-4, 1e-5),
    label_text = paste0("G", 1:10),
    stringsAsFactors = FALSE
  )
}

test_that("ev_select_labels mode = 'none' returns zero rows", {
  out <- ev_select_labels(
    make_lbl_input(),
    mode = "none", n = 5, rank_by = "significance",
    genes = NULL, p_col = "P.Value",
    p_threshold = 0.05, logfc_threshold = 0
  )
  expect_equal(nrow(out), 0L)
})

test_that("ev_select_labels mode = 'top_per_direction' picks n per side", {
  out <- ev_select_labels(
    make_lbl_input(),
    mode = "top_per_direction", n = 2, rank_by = "significance",
    genes = NULL, p_col = "P.Value",
    p_threshold = 0.05, logfc_threshold = 0
  )
  expect_equal(nrow(out), 4L)
  expect_true(any(out$logFC > 0))
  expect_true(any(out$logFC < 0))
})

test_that("ev_select_labels mode = 'top_total' returns the top n by significance", {
  out <- ev_select_labels(
    make_lbl_input(),
    mode = "top_total", n = 3, rank_by = "significance",
    genes = NULL, p_col = "P.Value",
    p_threshold = 0.05, logfc_threshold = 0
  )
  expect_equal(nrow(out), 3L)
})

test_that("ev_select_labels rank_by = 'logfc' sorts by |logFC|", {
  out <- ev_select_labels(
    make_lbl_input(),
    mode = "top_total", n = 2, rank_by = "logfc",
    genes = NULL, p_col = "P.Value",
    p_threshold = 0.05, logfc_threshold = 0
  )
  expect_equal(sort(out$gene), c("G1", "G10"))
})

test_that("ev_select_labels mode = 'all_significant' keeps every sig row", {
  out <- ev_select_labels(
    make_lbl_input(),
    mode = "all_significant", n = 0, rank_by = "significance",
    genes = NULL, p_col = "P.Value",
    p_threshold = 0.05, logfc_threshold = 0
  )
  expect_equal(nrow(out), 8L)
})

test_that("ev_select_labels mode = 'explicit' matches the genes vector", {
  out <- ev_select_labels(
    make_lbl_input(),
    mode = "explicit", n = 0, rank_by = "significance",
    genes = c("G2", "G9"), p_col = "P.Value",
    p_threshold = 0.05, logfc_threshold = 0
  )
  expect_equal(sort(out$gene), c("G2", "G9"))
})

test_that("ev_label_text prefers symbol, then uniprot, then gene", {
  d <- data.frame(
    gene = c("g1", "g2", "g3"),
    symbol = c("S1", NA, ""),
    uniprot = c("U1", "U2", NA),
    stringsAsFactors = FALSE
  )
  expect_equal(ev_label_text(d), c("S1", "U2", "g3"))
})

test_that("ev_clean_label strips canonical database prefixes", {
  expect_equal(
    ev_clean_label("HALLMARK_TCA_CYCLE"),
    stringr::str_wrap("TCA Cycle", width = 15)
  )
  expect_equal(
    ev_clean_label("GOBP_AUTOPHAGY"),
    stringr::str_wrap("Autophagy", width = 15)
  )
})

test_that("ev_clean_label routes MITOCARTA_ names through the leaf shortener", {
  out <- ev_clean_label("MITOCARTA_OXPHOS__CI_SUBUNITS")
  expect_true(grepl("Complex I", out))
})

# A ring arc has room for two lines. A third pushes the label box into its
# neighbours, so verbose MSigDB names need a phrase entry, not a wider wrap.
test_that("ev_clean_label keeps verbose MSigDB names within two lines", {
  verbose <- c(
    "GOBP_STRIATED_MUSCLE_CELL_DIFFERENTIATION",
    "GOBP_MEMBRANELESS_ORGANELLE_ASSEMBLY",
    "GOBP_CELLULAR_COMPONENT_ASSEMBLY_INVOLVED_IN_MORPHOGENESIS",
    "GOBP_PROTON_TRANSMEMBRANE_TRANSPORT",
    "GOBP_NUCLEOSIDE_TRIPHOSPHATE_BIOSYNTHETIC_PROCESS",
    "GOBP_RIBOSOMAL_SMALL_SUBUNIT_BIOGENESIS",
    "GOSLIM_PROTEIN_LOCALIZATION_TO_PLASMA_MEMBRANE",
    "KEGG_MEDICUS_REFERENCE_RAB7_REGULATED_MICROTUBULE_MINUS_END_DIRECTED_TRANSPORT",
    "REACTOME_SEPARATION_OF_SISTER_CHROMATIDS",
    "REACTOME_NON_INTEGRIN_MEMBRANE_ECM_INTERACTIONS",
    "REACTOME_REGULATION_OF_PD_L1_CD274_POST_TRANSLATIONAL_MODIFICATION",
    "REACTOME_ASPARAGINE_N_LINKED_GLYCOSYLATION"
  )
  lines <- lengths(strsplit(ev_clean_label(verbose), "\n", fixed = TRUE))
  expect_equal(lines, rep(2L, length(verbose)))
})

test_that("ev_clean_label capitalises gene and complex acronyms", {
  expect_equal(ev_clean_label("REACTOME_RRNA_PROCESSING"), "rRNA Processing")
  expect_equal(ev_clean_label("REACTOME_UCH_PROTEINASES"), "UCH Proteinases")
  expect_match(ev_clean_label("REACTOME_ECM_PROTEOGLYCANS"), "^ECM")
  expect_match(ev_clean_label("REACTOME_CYTOPROTECTION_BY_HMOX1"), "HMOX1")
  expect_match(
    ev_clean_label("REACTOME_REGULATION_OF_PD_L1_CD274_POST_TRANSLATIONAL_MODIFICATION"),
    "PD-L1"
  )
})

test_that("ev_clean_label does not repeat a word the acronym already carries", {
  expect_equal(ev_clean_label("GOBP_ELECTRON_TRANSPORT_CHAIN"), "ETC")
  expect_match(
    ev_clean_label("REACTOME_MITOTIC_G2_G2_M_PHASES"), "G2/M",
    fixed = TRUE
  )
})

test_that("volcano_ring with label_mode = 'top_per_direction' runs", {
  p <- suppressMessages(volcano_ring(
    make_toy_volc(), make_toy_enrich(),
    label_mode = "top_per_direction", label_n = 2
  ))
  expect_s3_class(p, "ggplot")
})

test_that("volcano_ring with label_mode = 'by_significance' runs", {
  p <- suppressMessages(volcano_ring(
    make_toy_volc(), make_toy_enrich(),
    label_mode = "by_significance", label_n = 3
  ))
  expect_s3_class(p, "ggplot")
})

test_that("volcano_ring with label_mode = 'by_genes' runs", {
  p <- suppressMessages(volcano_ring(
    make_toy_volc(), make_toy_enrich(),
    label_mode = "by_genes",
    label_genes = c("G1", "G15")
  ))
  expect_s3_class(p, "ggplot")
})
