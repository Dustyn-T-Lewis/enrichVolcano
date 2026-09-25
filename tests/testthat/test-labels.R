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

test_that("clean_label reads the shipped table, short by default", {
  expect_identical(clean_label("HALLMARK_OXIDATIVE_PHOSPHORYLATION"), "OXPHOS")
  oxphos <- "HALLMARK_OXIDATIVE_PHOSPHORYLATION"
  expect_identical(clean_label(oxphos, style = "clean", width = 40), "Oxidative Phosphorylation")
  expect_identical(clean_label("HALLMARK_OXIDATIVE_PHOSPHORYLATION", style = "clean"), "Oxidative\nPhosphorylation")
  expect_identical(clean_label("GOSLIM_LIPID_METABOLIC_PROCESS", width = 40), "Lipid Metabolism")
})

test_that("a line break written into the table is kept at any width", {
  term <- "GOSLIM_GENERATION_OF_PRECURSOR_METABOLITES_AND_ENERGY"
  expect_identical(clean_label(term), "Precursor Metabolites\n& Energy")
  expect_identical(clean_label(term, width = 40), "Precursor Metabolites\n& Energy")
})

test_that("every shipped label is filled and fits two ring lines", {
  tbl <- utils::read.csv(system.file("extdata", "term_labels.csv", package = "enrichVolcano"))
  expect_false(anyDuplicated(tbl$term) > 0)
  expect_true(all(nzchar(tbl$clean) & nzchar(tbl$short)))
  lines <- lengths(strsplit(clean_label(tbl$term), "\n", fixed = TRUE))
  expect_true(all(lines <= 2))
})

test_that("a term outside the table gets the plain clean", {
  expect_identical(clean_label("REACTOME_RESPIRATORY_ELECTRON_TRANSPORT", width = 40), "Respiratory Electron Transport")
  expect_identical(clean_label("GOBP_AUTOPHAGY"), "Autophagy")
  expect_identical(clean_label("KEGG_MEDICUS_REFERENCE_GLYCOLYSIS", width = 40), "Reference Glycolysis")
  expect_identical(clean_label("MY_OWN_SET", width = 40), "My Own Set")
  expect_identical(clean_label("REACTOME_RRNA_PROCESSING"), "rRNA Processing")
  expect_identical(clean_label("REACTOME_CYTOPROTECTION_BY_HMOX1", width = 40), "Cytoprotection By HMOX1")
  expect_identical(clean_label("GOBP_AUTOPHAGY", style = "clean"), clean_label("GOBP_AUTOPHAGY"))
})

test_that("the plain clean keeps the last level of a hierarchy", {
  expect_identical(clean_label("MITOCARTA_OXPHOS>Complex_IV", width = 40), "Complex IV")
  expect_identical(clean_label("MITOCARTA_OXPHOS__OXPHOS_SUBUNITS", width = 40), "OXPHOS Subunits")
  expect_identical(clean_label("MITOCARTA_Mitochondrial_central_dogma", width = 40), "Mitochondrial Central Dogma")
})

test_that("clean_label passes missing and empty names through", {
  expect_identical(clean_label(c(NA, "", "HALLMARK_APOPTOSIS")), c(NA, "", "Apoptosis"))
})

test_that("clean_label refuses an unknown style", {
  expect_error(clean_label("HALLMARK_APOPTOSIS", style = "tiny"), class = "rlang_error")
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

test_that("labels picks the shipped style", {
  terms <- c("HALLMARK_OXIDATIVE_PHOSPHORYLATION", "HALLMARK_TOY_A")
  expect_identical(display_labels(terms, "short", 40), c("OXPHOS", "Toy A"))
  expect_identical(display_labels(terms, "clean", 40), c("Oxidative Phosphorylation", "Toy A"))
  expect_identical(display_labels(terms, NULL, 40), display_labels(terms, "short", 40))
})

test_that("an edited table renames the terms it lists and leaves the rest short", {
  terms <- c("HALLMARK_OXIDATIVE_PHOSPHORYLATION", "HALLMARK_TOY_A", "HALLMARK_TOY_B")
  edited <- data.frame(
    term = c("HALLMARK_OXIDATIVE_PHOSPHORYLATION", "HALLMARK_TOY_B"),
    label = c("Mito\\nrespiration", NA)
  )
  expect_identical(display_labels(terms, edited, 40), c("Mito\nrespiration", "Toy A", "Toy B"))
})

test_that("a function labels every term its own way", {
  terms <- c("HALLMARK_OXIDATIVE_PHOSPHORYLATION", "GOBP_AUTOPHAGY")
  lower <- function(x) tolower(sub("^[A-Z]+_", "", x))
  expect_identical(display_labels(terms, lower, 40), c("oxidative_phosphorylation", "autophagy"))
  expect_error(display_labels(terms, function(x) "one", 40), class = "enrichVolcano_param_error")
})

test_that("the plots draw a table or function label", {
  p <- suppressMessages(plot_volcano_ring(make_toy_da(), make_toy_ring_enrichment(),
    databases = NULL, term_threshold = 1, n_terms = Inf,
    labels = data.frame(term = "HALLMARK_TOY_A", label = "Edited A")
  ))
  drawn <- unlist(lapply(ggplot2::ggplot_build(p)$data, function(d) if ("label" %in% names(d)) as.character(d$label)))
  expect_true("Edited A" %in% drawn)
  p <- suppressMessages(plot_scatter(make_toy_scatter_enrichment(), "A", "B", labels = function(x) paste("Set", x)))
  drawn <- unlist(lapply(ggplot2::ggplot_build(p)$data, function(d) if ("label" %in% names(d)) as.character(d$label)))
  expect_true("Set T01" %in% drawn)
})

test_that("labels must be a style, a named vector, a term and label table, or a function", {
  bad <- list("tiny", c(a = 1), stats::setNames("x", ""), data.frame(x = 1), list(a = "b"))
  for (b in bad) {
    expect_error(
      plot_volcano_ring(make_toy_da(), make_toy_ring_enrichment(), labels = b),
      class = "enrichVolcano_param_error"
    )
    expect_error(
      plot_scatter(make_toy_scatter_enrichment(), "A", "B", labels = b),
      class = "enrichVolcano_param_error"
    )
  }
})

test_that("clean_label returns an empty vector for no names", {
  expect_identical(clean_label(character(0)), character(0))
})
