test_that("print.enrichVolcano emits call and summary lines", {
  data <- make_ranked_input()
  paths <- list(UP = head(data$gene[data$logFC > 0], 10),
                DN = head(data$gene[data$logFC < 0], 10))
  p <- suppressWarnings(enrich_volcano(
    data, contrast = "ctr",
    databases = list(test = paths),
    enrich_mode = "fgsea", enrich_padj = 1.0
  ))
  expect_message(print(p), "enrichVolcano")
  expect_message(print(p), "Contrasts")
  expect_message(print(p), "Proteins")
})

test_that("ev_compose errors on name mismatch", {
  v <- list(a = ggplot2::ggplot())
  r <- list(b = ggplot2::ggplot())
  expect_error(ev_compose(v, r), class = "ev_compose_name_mismatch")
})
