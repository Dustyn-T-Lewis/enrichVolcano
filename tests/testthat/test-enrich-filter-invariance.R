test_that('filter_mode = "before" filters pathway list pre-enrichment', {
  data <- make_ranked_input()
  paths <- list(
    HALLMARK_GLYCOLYSIS = paste0("G", 1:20),
    HALLMARK_HYPOXIA    = paste0("G", 21:40),
    MITO_TRANSPORT      = paste0("G", 41:60)
  )
  out_before <- ev_enrich(data, contrast = "ctr",
                          databases = list(test = paths),
                          enrich_mode = "fgsea",
                          min_size = 5, max_size = 100,
                          include_terms = "^MITO",
                          exclude_terms = NULL,
                          filter_mode = "before")
  expect_setequal(unique(out_before$pathway), "MITO_TRANSPORT")
  # padj reflects 1 test in the filtered set, so it equals raw pval.
  expect_equal(out_before$padj, out_before$pval)
})

test_that('filter_mode = "display" filters after enrichment, padj over full set', {
  data <- make_ranked_input()
  paths <- list(
    HALLMARK_GLYCOLYSIS = paste0("G", 1:20),
    HALLMARK_HYPOXIA    = paste0("G", 21:40),
    MITO_TRANSPORT      = paste0("G", 41:60)
  )
  out_full <- ev_enrich(data, contrast = "ctr",
                        databases = list(test = paths),
                        enrich_mode = "fgsea",
                        min_size = 5, max_size = 100,
                        include_terms = NULL, exclude_terms = NULL)
  out_disp <- ev_enrich(data, contrast = "ctr",
                        databases = list(test = paths),
                        enrich_mode = "fgsea",
                        min_size = 5, max_size = 100,
                        include_terms = "^MITO",
                        exclude_terms = NULL,
                        filter_mode = "display")
  mito_row_full <- out_full[out_full$pathway == "MITO_TRANSPORT", ]
  mito_row_disp <- out_disp[out_disp$pathway == "MITO_TRANSPORT", ]
  expect_equal(mito_row_disp$padj, mito_row_full$padj)
  expect_equal(nrow(out_disp), 1)
})

test_that('GSEA ranking is independent of include/exclude filters', {
  data <- make_ranked_input()
  sub <- data[data$contrast == "ctr", ]
  r1 <- ev_rank_stats(sub, "signed_p")
  r2 <- ev_rank_stats(sub, "signed_p")
  expect_identical(r1, r2)
})

test_that('ORA universe is independent of include/exclude filters', {
  skip_if_not_installed("mockery")
  data <- make_ranked_input()
  paths <- list(
    A = paste0("G", 1:20),
    B = paste0("G", 21:40)
  )
  # mockery::stub only intercepts calls *originating in* `where`; stub
  # ev_ora_one inside ev_enrich so the universe argument is observable
  # regardless of whether ev_enrich was given an include_terms filter.
  captured <- list()
  capture_fn <- function(sub, paths, db_name, ctr, background,
                         min_size, max_size) {
    bg <- if (is.null(background)) sub$gene else background
    captured[[length(captured) + 1]] <<- bg
    NULL
  }
  mockery::stub(ev_enrich, "ev_ora_one", capture_fn)
  invisible(ev_enrich(data, contrast = "ctr",
                      databases = list(test = paths),
                      enrich_mode = "ora",
                      min_size = 5, max_size = 100,
                      include_terms = NULL, exclude_terms = NULL))
  invisible(ev_enrich(data, contrast = "ctr",
                      databases = list(test = paths),
                      enrich_mode = "ora",
                      min_size = 5, max_size = 100,
                      include_terms = "^A",
                      exclude_terms = NULL))
  expect_identical(captured[[1]], captured[[2]])
})
