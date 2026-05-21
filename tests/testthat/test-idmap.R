test_that("ev_load_idmap returns the bundled schema and caches", {
  tab <- ev_load_idmap("human")
  expect_identical(colnames(tab), c("accession", "symbol", "is_secondary"))
  expect_gt(nrow(tab), 20000)
  expect_identical(ev_load_idmap("human"), tab)   # cache returns same object
})

test_that("ev_load_idmap aborts on unsupported species", {
  expect_error(ev_load_idmap("frog"), class = "ev_no_idmap")
})

test_that("ev_map_uniprot maps a known human accession", {
  res <- ev_map_uniprot(c("P04637", "P00533"), "human")
  expect_identical(res$symbol, c("TP53", "EGFR"))
  expect_true(all(res$mapped))
  expect_false(any(res$isoform_stripped))
})

test_that("ev_map_uniprot strips isoform suffix but preserves identity", {
  res <- ev_map_uniprot("P04637-2", "human")
  expect_identical(res$symbol, "TP53")
  expect_true(res$isoform_stripped)
  expect_identical(res$uniprot, "P04637-2")
})

test_that("ev_map_uniprot resolves a secondary accession and flags it", {
  tab <- ev_load_idmap("human")
  sec <- tab$accession[tab$is_secondary][1]
  res <- ev_map_uniprot(sec, "human")
  expect_true(res$mapped)
  expect_true(res$via_secondary)
})

test_that("ev_map_uniprot returns NA for unmapped accessions", {
  res <- ev_map_uniprot("Q0FAKE9", "human")
  expect_false(res$mapped)
  expect_true(is.na(res$symbol))
})

test_that("ev_make_idmap_report counts mapping outcomes incl. secondary", {
  tab <- ev_load_idmap("human")
  sec <- tab$accession[tab$is_secondary][1]
  res <- ev_map_uniprot(c("P04637", "P04637-2", "Q0FAKE9", sec), "human")
  rep <- ev_make_idmap_report(res)
  expect_equal(rep$n_input, 4)
  expect_equal(rep$n_mapped, 3)             # P04637, P04637-2, sec
  expect_equal(rep$n_isoform_stripped, 1)   # only P04637-2
  expect_equal(rep$n_unmapped, 1)
  expect_identical(rep$unmapped, "Q0FAKE9")
  expect_equal(rep$n_via_secondary, 1)      # sec
})

test_that("ev_is_uniprot_accession matches accessions, not gene symbols", {
  valid <- c("P04637", "O15111", "Q9Y6K9", "A0A0C5B5G6",
             "P04637-2", "P04637-10")
  expect_true(all(ev_is_uniprot_accession(valid)))
  invalid <- c("TP53", "EGFR", "GAPDH", "BRCA1", "P04637-0", "", "P0")
  expect_false(any(ev_is_uniprot_accession(invalid)))
})

test_that("ev_idmap_report returns the attached report or NULL", {
  x <- tibble::tibble(a = 1)
  expect_null(ev_idmap_report(x))
  attr(x, "ev_idmap_report") <- list(n_input = 5L)
  expect_equal(ev_idmap_report(x)$n_input, 5L)
})
