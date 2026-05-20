# Parametric test: every msigdbr-backed DB in the registry must resolve
# to a non-empty pathway list for at least one supported species.
test_that("every msigdbr-backed registered DB resolves to non-empty gene set list", {
  skip_if_not_installed("msigdbr")
  reg <- list_databases()
  msigdbr_rows <- reg[reg$source == "msigdbr", ]
  expect_gt(nrow(msigdbr_rows), 0)
  for (i in seq_len(nrow(msigdbr_rows))) {
    db_name <- msigdbr_rows$name[i]
    species_list <- strsplit(msigdbr_rows$species[i], ";", fixed = TRUE)[[1]]
    species_pick <- species_list[1]  # use first supported species
    paths <- enrichVolcano:::ev_resolve_databases(db_name, species = species_pick)
    expect_true(length(paths) > 0,
                info = paste(db_name, "/", species_pick, "returned 0 pathway lists"))
    expect_true(length(paths[[1]]) > 0,
                info = paste(db_name, "/", species_pick, "has 0 gene sets"))
  }
})

test_that("ev_msigdbr_species maps all 7 supported species", {
  expect_silent(enrichVolcano:::ev_msigdbr_species("human"))
  expect_silent(enrichVolcano:::ev_msigdbr_species("mouse"))
  expect_silent(enrichVolcano:::ev_msigdbr_species("rat"))
  expect_silent(enrichVolcano:::ev_msigdbr_species("zebrafish"))
  expect_silent(enrichVolcano:::ev_msigdbr_species("fly"))
  expect_silent(enrichVolcano:::ev_msigdbr_species("yeast"))
  expect_silent(enrichVolcano:::ev_msigdbr_species("pig"))
  expect_error(enrichVolcano:::ev_msigdbr_species("unicorn"),
               class = "ev_unknown_species")
})

test_that("ev_msigdbr_category errors on unmapped name", {
  expect_error(enrichVolcano:::ev_msigdbr_category("not_a_db"),
               class = "ev_msigdbr_unmapped")
})

test_that("ev_resolve_databases warns on unsupported species and aborts on unknown DB", {
  expect_error(enrichVolcano:::ev_resolve_databases("not_a_db", species = "human"),
               class = "ev_unknown_database")
  # hpo is human-only — asking for mouse should warn + skip
  expect_warning(
    out <- enrichVolcano:::ev_resolve_databases("hpo", species = "mouse"),
    class = "ev_species_unsupported"
  )
  expect_length(out, 0)
})

test_that("ev_resolve_databases passes through a named list unchanged", {
  paths <- list(GROUP_A = paste0("G", 1:5))
  out <- enrichVolcano:::ev_resolve_databases(paths, species = "human")
  expect_identical(out, paths)
})

test_that("database keyword 'compartments' expands to compartment DBs", {
  expanded <- enrichVolcano:::ev_expand_db_keywords("compartments")
  expect_setequal(expanded, c("go_cc", "corum", "mitocarta3"))
})

test_that("database keyword 'pathways' expands and mixes with explicit names", {
  expanded <- enrichVolcano:::ev_expand_db_keywords(c("pathways", "go_cc"))
  expect_true(all(c("hallmark", "reactome", "go_bp", "go_cc") %in% expanded))
  expect_false("compartments" %in% expanded)
})

test_that("non-keyword names pass through unchanged", {
  expect_equal(enrichVolcano:::ev_expand_db_keywords(c("hallmark", "reactome")),
               c("hallmark", "reactome"))
})
