vdata <- function() {
  tibble::tibble(
    gene = paste0("G", 1:4),
    symbol = c("AAA", "BBB", "CCC", "DDD"),
    uniprot = paste0("P0000", 1:4),
    logFC = c(3, -2.5, 0.2, 2),
    P.Value = c(1e-5, 1e-4, 0.3, 1e-3),
    adj.P.Val = c(1e-4, 1e-3, 0.5, 1e-2),
    pi_eq2 = c(1e-4, 1e-3, 0.6, 1e-2),
    contrast = "C1"
  )
}

# The label layer is the only ggrepel GeomTextRepel layer. ev_volcano always
# adds it (empty data when nothing is labelled).
repel_layer_data <- function(p) {
  for (ly in p$layers) if (inherits(ly$geom, "GeomTextRepel")) return(ly$data)
  NULL
}

test_that("default mode none labels nothing", {
  p <- ev_volcano(vdata(), "C1")
  expect_equal(nrow(repel_layer_data(p)), 0)
})

test_that("top_per_direction labels n up and n down using label_text", {
  p <- ev_volcano(vdata(), "C1", label_mode = "top_per_direction",
                  label_n = 1)
  d <- repel_layer_data(p)
  expect_setequal(d$label_text, c("AAA", "BBB"))
})

test_that("explicit accession label resolves to symbol text", {
  p <- ev_volcano(vdata(), "C1", label_mode = "explicit",
                  label_genes = "P00004")
  d <- repel_layer_data(p)
  expect_identical(d$label_text, "DDD")
})

test_that("ev_volcano_ring labels nothing by default", {
  v <- vdata()
  enr <- tibble::tibble(contrast = "C1", database = "x", pathway = "P",
                        pval = 0.01, padj = 0.02, NES = 1.5, size = 20,
                        leading_edge = "AAA;DDD", mode = "fgsea",
                        direction = "up")
  p <- ev_volcano_ring(v, enr, title = "C1")
  expect_null(repel_layer_data(p))
})

test_that("ev_volcano_ring labels per the rule when asked", {
  v <- vdata()
  enr <- tibble::tibble(contrast = "C1", database = "x", pathway = "P",
                        pval = 0.01, padj = 0.02, NES = 1.5, size = 20,
                        leading_edge = "AAA;DDD", mode = "fgsea",
                        direction = "up")
  p <- ev_volcano_ring(v, enr, title = "C1",
                       label_mode = "top_per_direction", label_n = 1)
  d <- repel_layer_data(p)
  expect_setequal(d$label_text, c("AAA", "BBB"))
})

test_that("enrich_volcano accepts and threads label params without error", {
  data <- make_ranked_input()
  paths <- make_mini_pathways()
  p <- suppressWarnings(enrich_volcano(
    data, contrast = "ctr",
    databases = list(test = paths),
    enrich_padj = 1,
    label_mode = "top_per_direction", label_n = 5
  ))
  expect_s3_class(p, "enrichVolcano")
})
