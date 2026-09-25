make_lbl_input <- function() {
  data.frame(
    protein = paste0("P", 1:10),
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

test_that("ev_select_labels mode = 'by_significance' returns the top n by significance", {
  out <- ev_select_labels(
    make_lbl_input(),
    mode = "by_significance", n = 3, rank_by = "significance",
    genes = NULL, p_col = "P.Value",
    p_threshold = 0.05, logfc_threshold = 0
  )
  expect_equal(nrow(out), 3L)
})

test_that("ev_select_labels rank_by = 'logfc' sorts by |logFC|", {
  out <- ev_select_labels(
    make_lbl_input(),
    mode = "by_significance", n = 2, rank_by = "logfc",
    genes = NULL, p_col = "P.Value",
    p_threshold = 0.05, logfc_threshold = 0
  )
  expect_equal(sort(out$gene), c("G1", "G10"))
})

test_that("ev_select_labels mode = 'by_genes' matches symbols and accessions", {
  out <- ev_select_labels(
    make_lbl_input(),
    mode = "by_genes", n = 0, rank_by = "significance",
    genes = c("G2", "P9"), p_col = "P.Value",
    p_threshold = 0.05, logfc_threshold = 0
  )
  expect_equal(sort(out$gene), c("G2", "G9"))
})

test_that("a point without a gene symbol is labelled with its accession", {
  d <- make_lbl_input()
  d$gene[1:2] <- c(NA, "")
  out <- ev_select_labels(
    d,
    mode = "by_significance", n = 2, rank_by = "significance",
    genes = NULL, p_col = "P.Value", p_threshold = 0.05, logfc_threshold = 0
  )
  expect_setequal(out$label_text, c("P1", "G10"))
})

test_that("clean_label strips canonical database prefixes", {
  expect_equal(
    clean_label("HALLMARK_TCA_CYCLE"),
    stringr::str_wrap("TCA Cycle", width = 15)
  )
  expect_equal(
    clean_label("GOBP_AUTOPHAGY"),
    stringr::str_wrap("Autophagy", width = 15)
  )
})

test_that("clean_label routes MITOCARTA_ names through the leaf shortener", {
  out <- clean_label("MITOCARTA_OXPHOS__CI_SUBUNITS")
  expect_true(grepl("Complex I", out))
})

# A ring arc has room for two lines. A third pushes the label box into its
# neighbours, so verbose MSigDB names need a phrase entry, not a wider wrap.
test_that("clean_label keeps verbose MSigDB names within two lines", {
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
  lines <- lengths(strsplit(clean_label(verbose), "\n", fixed = TRUE))
  expect_equal(lines, rep(2L, length(verbose)))
})

test_that("clean_label capitalises gene and complex acronyms", {
  expect_equal(clean_label("REACTOME_RRNA_PROCESSING"), "rRNA Processing")
  expect_equal(clean_label("REACTOME_UCH_PROTEINASES"), "UCH Proteinases")
  expect_match(clean_label("REACTOME_ECM_PROTEOGLYCANS"), "^ECM")
  expect_match(clean_label("REACTOME_CYTOPROTECTION_BY_HMOX1"), "HMOX1")
  expect_match(
    clean_label("REACTOME_REGULATION_OF_PD_L1_CD274_POST_TRANSLATIONAL_MODIFICATION"),
    "PD-L1"
  )
})

test_that("clean_label does not repeat a word the acronym already carries", {
  expect_equal(clean_label("GOBP_ELECTRON_TRANSPORT_CHAIN"), "ETC")
  expect_match(
    clean_label("REACTOME_MITOTIC_G2_G2_M_PHASES"), "G2/M",
    fixed = TRUE
  )
})

test_that("plot_volcano_ring with label_mode = 'top_per_direction' runs", {
  p <- suppressMessages(plot_volcano_ring(
    make_toy_da(), make_toy_ring_enrichment(),
    label_mode = "top_per_direction", label_n = 2
  ))
  expect_s3_class(p, "ggplot")
})

test_that("plot_volcano_ring with label_mode = 'by_significance' runs", {
  p <- suppressMessages(plot_volcano_ring(
    make_toy_da(), make_toy_ring_enrichment(),
    label_mode = "by_significance", label_n = 3
  ))
  expect_s3_class(p, "ggplot")
})

test_that("plot_volcano_ring with label_mode = 'by_genes' runs", {
  p <- suppressMessages(plot_volcano_ring(
    make_toy_da(), make_toy_ring_enrichment(),
    label_mode = "by_genes",
    label_genes = c("G1", "G15")
  ))
  expect_s3_class(p, "ggplot")
})

test_that("the default width is unchanged, hand-placed breaks included", {
  expect_identical(clean_label("HALLMARK_HEME_METABOLISM"), "Heme\nMetabolism")
})

test_that("a wider width wraps less and skips the ring's hand-placed breaks", {
  expect_identical(clean_label("HALLMARK_HEME_METABOLISM", width = 40), "Heme Metabolism")
  long <- "GOBP_REGULATION_OF_CYTOPLASMIC_TRANSLATION_IN_RESPONSE_TO_STRESS"
  lines <- strsplit(clean_label(long, width = 40), "\n")[[1]]
  expect_true(all(nchar(lines) <= 40))
  expect_lt(length(lines), length(strsplit(clean_label(long), "\n")[[1]]))
})

test_that("width reaches vectorised and MitoCarta names", {
  out <- clean_label(c("HALLMARK_HEME_METABOLISM", "HALLMARK_MITOTIC_SPINDLE"), width = 40)
  expect_identical(out, c("Heme Metabolism", "Mitotic Spindle"))
  mito <- "MITOCARTA_OXPHOS__OXPHOS_ASSEMBLY_FACTORS"
  expect_false(grepl("\n", clean_label(mito, width = 40)))
  expect_true(grepl("\n", clean_label(mito)))
})

test_that("labels replace the cleaned names of the terms they name", {
  p <- suppressMessages(plot_volcano_ring(make_toy_da(), make_toy_ring_enrichment(),
    databases = NULL, term_threshold = 1, n_terms = Inf,
    labels = c(HALLMARK_TOY_A = "My own name for A", HALLMARK_TOY_C = "Line one\nline two")
  ))
  drawn <- unlist(lapply(ggplot2::ggplot_build(p)$data, function(d) if ("label" %in% names(d)) as.character(d$label)))
  expect_true(stringr::str_wrap("My own name for A", 15) %in% drawn)
  expect_true("Line one\nline two" %in% drawn)
  expect_true(clean_label("HALLMARK_TOY_B") %in% drawn)
})

test_that("the scatter uses the same labels", {
  p <- suppressMessages(plot_scatter(make_toy_scatter_enrichment(), "A", "B", labels = c(T01 = "First term")))
  drawn <- unlist(lapply(ggplot2::ggplot_build(p)$data, function(d) if ("label" %in% names(d)) as.character(d$label)))
  expect_true("First term" %in% drawn)
})

test_that("labels must be a named character vector", {
  for (bad in list("no names", c(a = 1), stats::setNames("x", ""))) {
    expect_error(
      plot_volcano_ring(make_toy_da(), make_toy_ring_enrichment(), labels = bad),
      class = "enrichVolcano_param_error"
    )
    expect_error(
      plot_scatter(make_toy_scatter_enrichment(), "A", "B", labels = bad),
      class = "enrichVolcano_param_error"
    )
  }
})

test_that("clean_label returns an empty vector for no names", {
  expect_identical(clean_label(character(0)), character(0))
})
