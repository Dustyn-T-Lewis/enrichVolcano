test_that("UniProt accessions map to current symbols, isoforms included", {
  skip_if_not_installed("org.Hs.eg.db")
  ids <- c("P31040", "O00483", "P31040-2", "NOT_AN_ID")
  expect_identical(unname(map_symbols(ids, "Homo sapiens")), c("SDHA", "COXFA4", "SDHA", NA))
})

test_that("mouse and rat use their own annotation packages", {
  skip_if_not_installed("org.Mm.eg.db")
  skip_if_not_installed("org.Rn.eg.db")
  expect_identical(unname(map_symbols("Q8K2B3", "Mus musculus")), "Sdha")
  expect_identical(unname(map_symbols("Q920L2", "Rattus norvegicus")), "Sdha")
})

test_that("an unsupported species is refused", {
  expect_error(map_symbols("P31040", "Danio rerio"), class = "enrichVolcano_param_error")
})

test_that("as_da looks up symbols for UniProt accessions", {
  skip_if_not_installed("org.Hs.eg.db")
  tbl <- data.frame(protein = c("O00483", "P31040"), gene = c("NDUFA4", "SDHA"), logFC = c(1, -1), p = c(0.01, 0.02))
  expect_message(d <- as_da(tbl, contrast = "A"), "2 of 2", class = "enrichVolcano_symbol_lookup")
  expect_identical(d$gene, c("COXFA4", "SDHA"))
})

test_that("species = NULL keeps the table's own symbols", {
  tbl <- data.frame(protein = c("O00483", "P31040"), gene = c("NDUFA4", "SDHA"), logFC = c(1, -1), p = c(0.01, 0.02))
  expect_identical(as_da(tbl, contrast = "A", species = NULL)$gene, c("NDUFA4", "SDHA"))
})

test_that("tables keyed by symbols rather than accessions keep them", {
  tbl <- data.frame(gene = c("SDHA", "CS"), logFC = c(1, -1), p = c(0.01, 0.02))
  d <- as_da(tbl, contrast = "A")
  expect_identical(d$gene, c("SDHA", "CS"))
  expect_identical(d$protein, c("SDHA", "CS"))
})

test_that("Hallmark loads 50 sets and records its version", {
  skip_if_not_installed("msigdbr")
  sets <- load_gene_sets("Hallmark", min_size = 1)
  expect_named(sets, "Hallmark")
  expect_length(sets$Hallmark, 50L)
  expect_true("SDHA" %in% sets$Hallmark$HALLMARK_OXIDATIVE_PHOSPHORYLATION)
  v <- attr(sets, "versions")
  expect_identical(v$species, "Homo sapiens")
  expect_identical(v$msigdbr, as.character(utils::packageVersion("msigdbr")))
})

test_that("the pinned GO slim maps to genes and records its release", {
  skip_if_not_installed("org.Hs.eg.db")
  sets <- load_gene_sets("GO Slim", min_size = 1, max_size = Inf)
  slim <- sets$`GO Slim`
  expect_gt(length(slim), 50)
  expect_true(all(grepl("^GOSLIM_[A-Z0-9_]+$", names(slim))))
  expect_true("GOSLIM_PROTEIN_FOLDING" %in% names(slim))
  expect_false("GOSLIM_MITOCHONDRION" %in% names(slim))
  expect_identical(attr(sets, "versions")$go_slim, "go/releases/2026-07-26/subsets/goslim_generic.owl")
})

test_that("set sizes are filtered to the requested window", {
  skip_if_not_installed("msigdbr")
  sets <- load_gene_sets("Hallmark", min_size = 150, max_size = 200)
  sizes <- lengths(sets$Hallmark)
  expect_true(all(sizes >= 150 & sizes <= 200))
})

test_that("mouse gene sets come back in mouse symbols", {
  skip_if_not_installed("msigdbr")
  sets <- load_gene_sets("Hallmark", species = "Mus musculus", min_size = 1)
  expect_true("Sdha" %in% sets$Hallmark$HALLMARK_OXIDATIVE_PHOSPHORYLATION)
})

test_that("unknown collections and species are refused", {
  expect_error(load_gene_sets("Pfam"), class = "enrichVolcano_param_error")
  expect_error(load_gene_sets("Hallmark", species = "Danio rerio"), class = "enrichVolcano_param_error")
})

test_that("a table with no known accession maps to nothing rather than failing", {
  skip_if_not_installed("org.Hs.eg.db")
  expect_identical(unname(map_symbols(c("A0A000", "B1B111"), "Homo sapiens")), c(NA_character_, NA_character_))
})

test_that("an accession the org package lacks keeps the table's own symbol", {
  skip_if_not_installed("org.Hs.eg.db")
  tbl <- data.frame(protein = c("A0A999Z999", "P31040"), gene = c("MYGENE", "SDHA"), logFC = c(1, -1), p = 0.01)
  expect_message(d <- as_da(tbl, contrast = "A"), "1 of 2", class = "enrichVolcano_symbol_lookup")
  expect_identical(d$gene, c("MYGENE", "SDHA"))
})

test_that("a protein group is looked up by its first accession", {
  skip_if_not_installed("org.Hs.eg.db")
  tbl <- data.frame(protein = c("P31040;Q9UBK2", "Q9UBK2"), logFC = c(1, -1), p = 0.01)
  d <- suppressMessages(as_da(tbl, contrast = "A"))
  expect_identical(d$protein, c("P31040;Q9UBK2", "Q9UBK2"))
  expect_identical(d$gene, c("SDHA", "PPARGC1A"))
})

test_that("a blank symbol counts as missing", {
  skip_if_not_installed("org.Hs.eg.db")
  tbl <- data.frame(protein = c("A0A999Z999", "P31040"), gene = c("", "SDHA"), logFC = c(1, -1), p = 0.01)
  d <- suppressMessages(as_da(tbl, contrast = "A"))
  expect_identical(d$gene, c(NA_character_, "SDHA"))
})

test_that("the symbol lookup prints no annotation startup message", {
  skip_if_not_installed("org.Rn.eg.db")
  if (isNamespaceLoaded("org.Rn.eg.db")) try(unloadNamespace("org.Rn.eg.db"), silent = TRUE)
  skip_if(isNamespaceLoaded("org.Rn.eg.db"), "org.Rn.eg.db is already loaded")
  expect_no_condition(map_symbols("P04797", "Rattus norvegicus"), class = "packageStartupMessage")
})
