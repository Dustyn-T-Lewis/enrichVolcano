test_that("UniProt input maps to symbol and sets identity columns", {
  df <- data.frame(
    UniProt = c("P04637", "P00533-2", "Q0FAKE9"),
    logFC = c(2, -1, 0.5),
    P.Value = c(1e-4, 1e-3, 0.2),
    contrast = "C1",
    stringsAsFactors = FALSE
  )
  out <- ev_validate(df, species = "human")
  expect_identical(out$uniprot, c("P04637", "P00533-2", "Q0FAKE9"))
  expect_identical(out$symbol, c("TP53", "EGFR", NA_character_))
  # unmapped row keeps the accession as its enrichment key
  expect_identical(out$gene, c("TP53", "EGFR", "Q0FAKE9"))
  rep <- ev_idmap_report(out)
  expect_equal(rep$n_mapped, 2)
  expect_equal(rep$n_unmapped, 1)
})

test_that("UniProt input without species aborts", {
  df <- data.frame(UniProt = "P04637", logFC = 1, P.Value = 0.01)
  expect_error(ev_validate(df), class = "ev_uniprot_no_species")
})

test_that("symbol-only input is unaffected (backward compatible)", {
  df <- data.frame(gene = c("TP53", "EGFR"), logFC = c(1, -1),
                   P.Value = c(0.01, 0.02), contrast = "C1",
                   stringsAsFactors = FALSE)
  out <- ev_validate(df)
  expect_identical(out$gene, c("TP53", "EGFR"))
  expect_true(all(is.na(out$uniprot)))
  expect_identical(out$symbol, c("TP53", "EGFR"))
  expect_null(ev_idmap_report(out))
})
