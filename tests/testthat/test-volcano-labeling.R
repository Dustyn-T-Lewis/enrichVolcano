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
