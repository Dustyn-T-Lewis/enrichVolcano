test_that("enrich_volcano returns enrichVolcano patchwork end-to-end", {
  data <- make_ranked_input()
  paths <- make_mini_pathways()
  # Use direct named-list path so we don't depend on msigdbr in tests
  p <- suppressWarnings(enrich_volcano(
    data, contrast = "ctr",
    databases = list(test = paths),
    p_method = "pi_eq2", p_adjust = "BH",
    enrich_mode = "fgsea", enrich_padj = 1.0
  ))
  expect_s3_class(p, "enrichVolcano")
  expect_true(!is.null(attr(p, "ev_call")))
  expect_named(attr(p, "ev_data"),
               c("validated_input", "pi_scores", "enrichment", "dedup_result"))
})

test_that("enrich_volcano handles multiple contrasts", {
  data <- make_ranked_input()
  data2 <- data
  data2$contrast <- "ctr2"
  data_both <- dplyr::bind_rows(data, data2)
  paths <- make_mini_pathways()
  p <- suppressWarnings(enrich_volcano(
    data_both, contrast = c("ctr", "ctr2"),
    databases = list(test = paths),
    p_method = "pi_eq2", p_adjust = "BH",
    enrich_mode = "fgsea", enrich_padj = 1.0
  ))
  expect_s3_class(p, "enrichVolcano")
  expect_equal(length(unique(attr(p, "ev_data")$validated_input$contrast)), 2)
})

test_that("enrich_volcano forwards rank_by to ev_enrich", {
  data <- make_ranked_input()  # has gene/contrast/logFC/P.Value/pi_eq2, contrast "ctr"
  paths <- make_mini_pathways()
  # default signed_p works
  p1 <- suppressWarnings(enrich_volcano(
    data, contrast = "ctr", databases = list(test = paths),
    enrich_mode = "fgsea", enrich_padj = 1
  ))
  expect_s3_class(p1, "enrichVolcano")
  # explicit t ranking is accepted (add a t column)
  data$t <- data$logFC / 0.3
  p2 <- suppressWarnings(enrich_volcano(
    data, contrast = "ctr", databases = list(test = paths),
    enrich_mode = "fgsea", enrich_padj = 1, rank_by = "t"
  ))
  expect_s3_class(p2, "enrichVolcano")
  # unsigned pi-score is still rejected through the hero
  expect_error(
    enrich_volcano(data, contrast = "ctr", databases = list(test = paths),
                   enrich_mode = "fgsea", rank_by = "pi_eq2"),
    class = "ev_unsigned_ranking"
  )
})
